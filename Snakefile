from pathlib import Path
configfile: "config.json"
GLOBAL_REF_PATH = config["globalResources"] 
GLOBAL_TMPD_PATH = config["globalTmpdPath"]

ref_type = list(config["libraries"].values())[0]["reference"]
library_name = list(config["libraries"].keys())[0]
sample_hashes = list(config["libraries"][library_name]["samples"].keys())

#TODO how to get organism
reference_path = os.path.join(GLOBAL_REF_PATH, "homo_sapiens", ref_type, "seq", ref_type + ".fa")
#basecaller_location = os.path.join(GLOBAL_TMPD_PATH, "ont-guppy-cpu/bin/guppy_basecaller")
basecaller_location = os.path.join(GLOBAL_TMPD_PATH, "/tmp/dorado-0.5.3-linux-x64/bin/dorado")

sample_names = []
for sample in sample_hashes:
    sample_name = config["libraries"][library_name]["samples"][sample]["sample_name"]
    sample_names.append(sample_name)

include: "rules/rules.smk"

rule all:
    input:
        expand("{library_name}/outputs/{sample_name}/reads_merged.bam", library_name = library_name, sample_name = sample_names),
        "qc_reports/alignment_multiqc/multiqc.html"

        
