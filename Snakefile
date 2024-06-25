from pathlib import Path
import pandas as pd

configfile: "config.json"
GLOBAL_REF_PATH = config["globalResources"] 
GLOBAL_TMPD_PATH = config["globalTmpdPath"]

##### BioRoot utilities - reference #####
module BR:
    snakefile: gitlab("bioroots/bioroots_utilities", path="bioroots_utilities.smk",branch="master")
    config: config

use rule * from BR as other_*
config = BR.load_organism()

# setting organism from reference
f = open(os.path.join(GLOBAL_REF_PATH,"reference_info","reference2.json"),)
reference_dict = json.load(f)
f.close()

config["species_name"] = [organism_name for organism_name in reference_dict.keys() if isinstance(reference_dict[organism_name],dict) and config["reference"] in reference_dict[organism_name].keys()][0]
config["organism"] = config["species_name"].split(" (")[0].lower().replace(" ","_")
if len(config["species_name"].split(" (")) > 1:
    config["species"] = config["species_name"].split(" (")[1].replace(")","")

# Folders
#
reference_path = os.path.join(GLOBAL_REF_PATH,config["organism"], config["reference"], "seq", config["reference"] + ".fa")
sample_hashes = list(config["samples"].keys())
basecaller_location = os.path.join(GLOBAL_TMPD_PATH, "dorado-0.5.3-linux-x64/bin/dorado")

#hash_to_path = {}

# sample_names - from BR - not really sure how that works
# sample_tab = BR.load_sample()

# wildcard_constraints:
#     sample_name = "|".join(sample_tab.sample_name)

sample_names = []
for sample in sample_hashes:
    sample_name = config["samples"][sample]["sample_name"]
    sample_names.append(sample_name)
    #hash_to_path[sample]=os.path.join("raw_reads", sample_name, sample_name + ".pod5") #TODO add {library_name} when copy to copy_raw_data

##### Target rules #####
rule all:
    input:
        expand('aligned/{sample_name}/{sample_name}.bam', sample_name = sample_names),
        "qc_reports/all_samples/multiqc_report.html"

##### Modules #####   
include: "rules/rules.smk"

##### BioRoot utilities - prepare reference #####
# module PR:
#     snakefile: gitlab("bioroots/bioroots_utilities", path="prepare_reference.smk", branch="master")
#     config: config

# use rule * from PR as other_*
