from pathlib import Path
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

##### Config processing #####
# Folders
#
reference_path = os.path.join(GLOBAL_REF_PATH,config["organism"], config["reference"], "seq", config["reference"] + ".fa")

#ref_type = list(config["libraries"].values())[0]["reference"]
library_name = list(config["libraries"].keys())[0]
sample_hashes = list(config["libraries"][library_name]["samples"].keys())

#basecaller_location = os.path.join(GLOBAL_TMPD_PATH, "ont-guppy-cpu/bin/guppy_basecaller")
basecaller_location = os.path.join(GLOBAL_TMPD_PATH, "/tmp/dorado-0.5.3-linux-x64/bin/dorado")

sample_names = []
for sample in sample_hashes:
    sample_name = config["libraries"][library_name]["samples"][sample]["sample_name"]
    sample_names.append(sample_name)

##### Target rules #####
rule all:
    input:
        expand('{library_name}/aligned/{sample_name}/{sample_name}.bam', library_name = library_name, sample_name = sample_names),
        "qc_reports/all_samples/multiqc.html"

##### Modules #####   
include: "rules/rules.smk"

##### BioRoot utilities - prepare reference #####
module PR:
    snakefile: gitlab("bioroots/bioroots_utilities", path="prepare_reference.smk",branch="master")
    config: config

use rule * from PR as other_*
