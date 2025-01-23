from pathlib import Path
import pandas as pd
import json

configfile: "config.json"
GLOBAL_REF_PATH = config["globalResources"] 
GLOBAL_TMPD_PATH = config["globalTmpdPath"]

##### BioRoot utilities - reference #####
module BR:
    snakefile: github("BioIT-CEITEC/bioroots_utilities", path="bioroots_utilities.smk",branch="master")
    config: config

use rule * from BR as other_*
config = BR.load_organism()
sample_tab = BR.load_sample()

sample_hashes = list(config["samples"].keys())
basecaller_location = os.path.join(GLOBAL_TMPD_PATH, "dorado-0.5.3-linux-x64/bin/dorado")

sample_names = []

for sample in sample_hashes:
    sample_name = config["samples"][sample]["sample_name"]
    sample_names.append(sample_name)
    #hash_to_path[sample]=os.path.join("raw_reads", sample_name, sample_name + ".pod5") #TODO add {library_name} when copy to copy_raw_data

def input_methylation(wildcard):
    if config["6mA_methylation"] or config["5mC_methylation"]:
        return expand("methylation/{sample_name}/{sample_name}_modkit.bed", sample_name = sample_names)

##### Target rules #####
rule all:
    input:
        expand('aligned/{sample_name}/{sample_name}.bam', sample_name = sample_names),
        expand('aligned/{sample_name}/{sample_name}_sorted.bam', sample_name = sample_names),
        expand('summary/{sample_name}/{sample_name}_summary.tsv', sample_name = sample_names),
        #expand("methylation/{sample_name}/{sample_name}_modkit.tsv", sample_name = sample_names),
        # expand('aligned/{sample_name}/{sample_name}.bam', sample_name = "test23"),
        # expand('aligned/{sample_name}/{sample_name}_sorted.bam', sample_name = "test23"),
        # expand('summary/{sample_name}/{sample_name}_summary.tsv', sample_name = "test23"),
        # expand('summary/{sample_name}/{sample_name}_pycoQC.html', sample_name = "test23"),
        "qc_reports/all_samples/multiqc_report.html",
        input_methylation
        
##### Modules #####   
include: "rules/rules.smk"

##### BioRoot utilities - prepare reference #####
# module PR:
#     snakefile: gitlab("bioroots/bioroots_utilities", path="prepare_reference.smk", branch="master")
#     config: config

# use rule * from PR as other_*
