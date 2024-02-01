import json
import glob
import os

# TODO incorporate dorado
# https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.4.1-linux-x64.tar.gz
#TODO add option to copy fastq files from online basecalling,
#TODO add option to download cpu vs gpu version of guppy,

def get_exp_info(library_name):
    pattern = os.path.join(library_name, 'report_*.json')
    files = glob.glob(pattern)
    assert len(files) == 1, 'Number of configs found != 1'

    f = open(files[0])
    run_config = json.load(f)
    #TODO which one is right?
    # flowcell = run_config['protocol_run_info']['user_info']['user_specified_flow_cell_id']
    flowcell = run_config['protocol_run_info']['meta_info']['tags']['flow cell']['string_value']

    protocol_group_id = run_config['protocol_run_info']['user_info']['protocol_group_id']
    sample_id = run_config['protocol_run_info']['user_info']['sample_id']
    experiment_name = protocol_group_id+'_'+sample_id
    kit = run_config['protocol_run_info']['meta_info']['tags']['kit']['string_value']

    return {
        'flowcell':flowcell,
        'experiment_name':experiment_name,
        'library_name':library_name,
        'kit':kit,
        'save_path':'outputs/basecalling/'+library_name+'/guppy/',
    }

# TODO extrahovat do params 
rule download_basecaller_tar:
    output: '/tmp/ont-guppy-cpu_6.4.6_linux64.tar.gz'
    shell:
        f"""
        wget -P {GLOBAL_TMPD_PATH} https://cdn.oxfordnanoportal.com/software/analysis/ont-guppy-cpu_6.4.6_linux64.tar.gz
        """
        
rule extract_basecaller:
    input:'/tmp/ont-guppy-cpu_6.4.6_linux64.tar.gz'
    output: basecaller_location
    shell:
        """
        cd /tmp ;
        tar -xf {input}
        """

rule basecalling:
    input: 
        #TODO take only pass?
        library_path = config["run_dir"] + "/fast5_pass",
        #TODO generalize to cpu or gpu
        basecaller_location = basecaller_location,
    output:
        '{library_name}/basecalling/guppy/sequencing_summary.txt'
    params:
        kit = lambda wildcards: get_exp_info(wildcards.library_name)['kit'],
        flowcell = lambda wildcards: get_exp_info(wildcards.library_name)['flowcell'],
        save_path = lambda wildcards: get_exp_info(wildcards.library_name)['save_path'],
    threads: workflow.cores * 0.75
    resources: gpus=1
    shell:
        """
        {input.basecaller_location} \
            --quiet \
            --flowcell {params.flowcell} \
            --kit {params.kit} \
            --records_per_fastq 0 \
            --trim_strategy none \
            --save_path {params.save_path} \
            --recursive \
            --num_callers {threads} \
            --chunks_per_runner 512 \
            --calib_detect \
            --input_path {input.library_path} \
            -q 0 \
            2>&1; \
        """

rule merge_fastq_files:
    input:
        '{library_name}/basecalling/guppy/sequencing_summary.txt'
    output:
        "{library_name}/basecalling/guppy/reads.fastq"
    conda:
        "../envs/merge_fastq.yaml"
    shell:
        """
        if [ -d {wildcards.library_name}/basecalling/guppy/pass ]; then cat {wildcards.library_name}/basecalling/guppy/pass/fastq_runid*.fastq > {output}; \
        else cat {wildcards.library_name}/basecalling/guppy/fastq_runid*.fastq > {output}; fi
        """

rule align_to_genome:
    input:
        reads="{library_name}/basecalling/guppy/reads.fastq"
    params: 
        reference_path = reference_path
    output:
        bam = '{library_name}/alignment/minimap2/reads-align.genome.sorted.bam',
        bai = '{library_name}/alignment/minimap2/reads-align.genome.sorted.bam.bai',
    conda:
        "../envs/alignment.yaml"
    threads: 32
    shell:
        """
		minimap2 \
			-x splice \
			-a \
			-t {threads} \
			-u b \
			-p 1 \
			--secondary=no \
			{params.reference_path} \
			{input.reads} \
			| samtools view -bh -\
			| samtools sort --threads {threads} \
			> {output.bam}  
		samtools index {output.bam}
		"""

rule SV_calling:
    input: 
        bam = '{library_name}/alignment/minimap2/reads-align.genome.sorted.bam'
    output:
        vcf = '{library_name}/sv_calling/variants.vcf'
    params: reference_path = reference_path,
    conda: 
        "../envs/svim_environment.yaml"
    shell:
        """
        svim alignment {wildcards.library_name}/sv_calling/ {input.bam} {params.reference_path} 
        """