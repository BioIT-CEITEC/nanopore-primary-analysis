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
        pod5_path = "{library_name}/raw_reads/{sample_name}/{sample_name}.pod5",
        basecaller_location = basecaller_location,
    output:
        '{library_name}/aligned/{sample_name}/{sample_name}.bam'
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
        '{library_name}/aligned/{sample_name}/{sample_name}.bam'
    output:
        '{library_name}/aligned/{sample_name}/{sample_name}.bam.bai'
    conda:
        "../envs/alignment.yaml" #TODO make yaml just for indexing (samtools only)
    shell:
        "samtools index {output.bam}"

# TODO quality control, we need to convert to fastq files 
rule convert_to_fastq:
    input:
        '{library_name}/aligned/{sample_name}/{sample_name}.bam'
    output:
        sorted_bam = '{library_name}/aligned/{sample_name}/{sample_name}_sorted.bam',   
        fastq = '{library_name}/raw_reads/{sample_name}/{sample_name}.fastq'
    conda:
        "../envs/alignment.yaml" 
    shell:
        """
        samtools sort -n {input} -o {output.sorted_bam}
        samtools fastq -T "*" {output.sorted_bam} > {output.fastq}
        """

# TODO quality control after alignment
rule alignment_multiqc:
    input: '{library_name}/raw_reads/{sample_name}/{sample_name}.fastq'
    output: html= "qc_reports/all_samples/multiqc.html"
    conda:
        "../envs/alignment_multiqc.yaml"
    shell:
        """
        multiqc -f -n {output.html} {input}
        """

rule sequencing_summary:
    input: '{library_name}/aligned/{sample_name}/{sample_name}.bam'
    output: '{library_name}/summary/{sample_name}/{sample_name}_summary.tsv'
    shell:
        """
        {input.basecaller_location} summary {input} > {output}
        """"