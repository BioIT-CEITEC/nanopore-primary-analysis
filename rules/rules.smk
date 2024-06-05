import json
import glob
import os

rule download_dorado_tar:
    output: '/tmp/dorado-0.5.3-linux-x64.tar.gz'
    shell:
        f"""
        wget -P {GLOBAL_TMPD_PATH} https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.5.3-linux-x64.tar.gz
        """
rule extract_dorado:
    input:'dorado-0.5.3-linux-x64.tar.gz'
    output: basecaller_location
    shell:
        f"""
        cd {GLOBAL_TMPD_PATH} ;
        tar -xf {input}
        """
        
rule supafixed_basecalling_dorado:
    input: 
        pod5_path = "{library_name}/outputs/{sample_name}/reads_merged.pod5",
        basecaller_location = basecaller_location,
    output:
        '{library_name}/outputs/{sample_name}/reads_merged.bam'
    params:
        reference_path = reference_path
    threads: workflow.cores * 0.75
    resources: gpus=1
    shell:
        """
        {input.basecaller_location} basecaller hac {input.pod5_path} --reference {params.reference_path} > {output}
        """"

rule supafixed_indexing:
    input:
        '{library_name}/outputs/{sample_name}/reads_merged.bam'
    output:
        '{library_name}/outputs/{sample_name}/reads_merged.bam.bai'
    conda:
        "../envs/alignment.yaml" #TODO make yaml just for indexing (samtools only)
    shell:
        "samtools index {output.bam}"

# TODO quality control after alignment
# rule alignment_multiqc:
#     input: bam = 'outputs/{wildcards.library_name}/{wildcards.sample_name}/basecalling/reads_merged.bam',
#     output: html= "outputs/qc_reports/alignment_multiqc/multiqc.html"
#     conda:
#         "../envs/alignment_multiqc.yaml"