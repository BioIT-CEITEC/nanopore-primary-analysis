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
        #reference_path = reference_path,
        genome = config["organism_fasta"],
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
        else
            {input.basecaller_location} basecaller {params.dorado_model}{params.methylation} {input.pod5_path} --reference {params.genome} > {output}
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
        rm {input} # remove unsorted bam to save space
        """

rule create_stats_from_bam:
    input:
        'aligned/{sample_name}/{sample_name}_sorted.bam'
    output:
        stats='aligned/{sample_name}/{sample_name}_stats.txt'
    conda:
        "../envs/alignment.yaml" 
    shell:
        """
        samtools stats {input} > {output.stats}
        """

rule sequencing_summary:
    input: 
        basecaller_location = basecaller_location,
        bam = 'aligned/{sample_name}/{sample_name}_sorted.bam'
    output: 'summary/{sample_name}/{sample_name}_summary.tsv'
    shell:
        """
        {input.basecaller_location} summary {input.bam} > {output}
        """

rule sequencing_QC:
    input: 
        'summary/{sample_name}/{sample_name}_summary.tsv'
    output: 
        html = 'qc_reports/{sample_name}/{sample_name}_pycoQC.html',
        json = 'qc_reports/{sample_name}/{sample_name}_pycoQC.json'
    conda:
        "../envs/alignment.yaml" 
    shell:
        """
        pycoQC -f {input} -o {output.html} --json_out {output.json}
        """

rule filter_primary_chromosomes:
    input: 
        bam = 'aligned/{sample_name}/{sample_name}_sorted.bam'
    output: 
        primary_bam = 'aligned/{sample_name}/{sample_name}_primary.bam',
        primary_bam_index = 'aligned/{sample_name}/{sample_name}_primary.bam.bai'
    conda:
        "../envs/alignment.yaml" 
    shell:
        """
        samtools view -b {input.bam} \
        $(samtools idxstats {input.bam} | cut -f1 | grep -v '_' | grep -v 'chrUn' | grep -v 'alt' | grep -v 'fix' | grep -v 'random') \
        > {output.primary_bam}
        samtools index {output.primary_bam}
        """

rule sequencing_qc_qualimap:
    input: 
        bam = 'aligned/{sample_name}/{sample_name}_primary.bam'
    output: 
        html = 'qc_reports/{sample_name}/qualimap/{sample_name}_primary_qualimap.html',
        json = 'qc_reports/{sample_name}/qualimap/{sample_name}_primary_qualimap.json'
    conda:
        "../envs/qualimap.yaml"
    shell:
        """
        qualimap bamqc -bam {input.bam} -outdir qc_reports/{wildcards.sample_name}/qualimap -outformat HTML --java-mem-size=8G
        mv qc_reports/{wildcards.sample_name}/qualimapReport.json {output.json}
        mv qc_reports/{wildcards.sample_name}/qualimapReport.html {output.html}
        """
    
rule alignment_multiqc:
    input: expand('qc_reports/{sample_name}/{sample_name}_pycoQC.json', sample_name = sample_names)
    output: html= "qc_reports/all_samples/multiqc_report.html"
    params: outdir = "qc_reports/all_samples"
    conda:
        "../envs/alignment_multiqc.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        multiqc --force -o {params.outdir} {input}
        """

# samtools view aligned/test22/test22_sorted.bam | awk '{for(i=12;i<=NF;i++) if($i ~ /^MM:Z:/ || $i ~ /^ML:B:C:/) print $1, $i}' > aligned/test22/methylation_output.txt
rule create_modified_table:
    input: bam = 'aligned/{sample_name}/{sample_name}_sorted.bam'
    output: 
    #    tsv = "methylation/{sample_name}/{sample_name}_modkit.tsv",
        bed = "methylation/{sample_name}/{sample_name}_modkit.bed"
    params: outdir = "methylation/{sample_name}"
    threads: workflow.cores * 0.75
    conda:
        "../envs/methylation.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        modkit pileup {input.bam} {output.bed} --log-filepath {params.outdir}/pileup.log --filter-threshold C:0.8 --filter-threshold A:0.9 -t {threads}
        """
#modkit extract full {input.bam} {output.tsv}
