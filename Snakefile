from pathlib import Path
configfile: "config.json"
GLOBAL_REF_PATH = config["globalResources"] 
GLOBAL_TMPD_PATH = config["globalTmpdPath"]

#reference_directory = os.path.join(GLOBAL_REF_PATH,config["organism"],config["reference"])
ref_type = list(config["libraries"].values())[0]["reference"]
library = list(config["libraries"].keys())[0]
#TODO how to get organism
reference_path = os.path.join(GLOBAL_REF_PATH, "homo_sapiens", ref_type, "seq", ref_type + ".fa")
basecaller_location = os.path.join(GLOBAL_TMPD_PATH, "ont-guppy-cpu/bin/guppy_basecaller")

include: "rules/rules.smk"

rule all:
    input:
        #expand("{folder}/conversion_ok.txt", folder = config["run_dir"])
        expand("{library_name}/sv_calling/variants.vcf", library_name = library)
        #expand("outputs/alignment/{library_name}/minimap2/reads-align.genome.sorted.bam", library_name = config["run_dir"])
        #expand("outputs/basecalling/{library_name}/guppy/sequencing_summary.txt", library_name = config["run_dir"])
