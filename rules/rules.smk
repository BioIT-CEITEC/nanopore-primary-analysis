import json
import glob
import os

rule download_dorado_tar:
    output: expand('{GLOBAL_TMPD_PATH}/dorado-0.5.3-linux-x64.tar.gz', GLOBAL_TMPD_PATH = GLOBAL_TMPD_PATH)
    shell:
        f"""
        wget -P {GLOBAL_TMPD_PATH} https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.5.3-linux-x64.tar.gz
        """
rule extract_dorado:
    input: expand('{GLOBAL_TMPD_PATH}/dorado-0.5.3-linux-x64.tar.gz', GLOBAL_TMPD_PATH = GLOBAL_TMPD_PATH)
    output: basecaller_location
    shell:
        f"""
        cd {GLOBAL_TMPD_PATH} ;
        tar -xf {input}
        """

# Detect methylation changes based on parameter in config
METHYLATION = ""

if config["6mA_methylation"] and config["5mC_methylation"]:
    METHYLATION = ",5mCG_5hmCG,6mA"
elif config["5mC_methylation"]:
    METHYLATION = ",5mCG_5hmCG "
elif config["6mA_methylation"]:
    METHYLATION = ",6mA"

rule supafixed_basecalling_dorado:
    input: 
        pod5_path = "raw_reads/{sample_name}/{sample_name}.pod5",
        basecaller_location = basecaller_location,
    output:
        'aligned/{sample_name}/{sample_name}.bam'
    params:
        dorado_model=config["dorado_model"],
        non_empty_input=lambda wildcards, input: int(os.path.getsize(input.pod5_path) != 0), # converting to 0/1 to bash if
        reference_path = reference_path,
        sample_names = sample_names,
        dirname = 'aligned/{sample_name}',
        methylation = METHYLATION
    threads: workflow.cores * 0.75
    resources: gpus=1
    shell:
        """
        mkdir -p {params.dirname}
        if [ {params.non_empty_input} -eq 0 ]; then
            touch {output}
        elif 
            {input.basecaller_location} basecaller {params.dorado_model}{params.methylation} {input.pod5_path} --reference {params.reference_path} > {output}
        fi
        """

rule supafixed_sorting_indexing:
    input:
        'aligned/{sample_name}/{sample_name}.bam'
    output:
        sorted='aligned/{sample_name}/{sample_name}_sorted.bam'
    conda:
        "../envs/alignment.yaml" #TODO make yaml just for indexing (samtools only)
    shell:
        """
        samtools sort {input} -o {output.sorted}
        samtools index {output.sorted}
        """

# # TODO quality control, we need to convert to fastq files 
# rule convert_to_fastq:
#     input:
#         'aligned/{sample_name}/{sample_name}.bam'
#     output:
#         sorted_bam = 'aligned/{sample_name}/{sample_name}_sorted.bam',   
#         fastq = 'raw_reads/{sample_name}/{sample_name}.fastq'
#     conda:
#         "../envs/alignment.yaml" 
#     shell:
#         """
#         samtools sort -n {input} -o {output.sorted_bam}
#         samtools fastq -T "*" {output.sorted_bam} > {output.fastq}
#         """

rule create_stats_from_bam:
    input:
        'aligned/{sample_name}/{sample_name}.bam'
    output:
        stats='aligned/{sample_name}/{sample_name}_stats.txt'
    conda:
        "../envs/alignment.yaml" 
    shell:
        """
        samtools stats {input} > {output.stats}
        """

# # TODO quality control after alignment
# rule alignment_multiqc:
#     input: expand('aligned/{sample_name}/{sample_name}_stats.txt', sample_name = sample_names)
#     output: html= "qc_reports/all_samples/multiqc_report.html"
#     params: outdir = "qc_reports/all_samples"
#     conda:
#         "../envs/alignment_multiqc.yaml"
#     shell:
#         """
#         mkdir -p {params.outdir}
#         multiqc --force -o {params.outdir} {input}
#         """

rule sequencing_summary:
    input: 
        basecaller_location = basecaller_location,
        bam = 'aligned/{sample_name}/{sample_name}.bam'
    output: 'summary/{sample_name}/{sample_name}_summary.tsv'
    shell:
        """
        {input.basecaller_location} summary {input.bam} > {output}
        """