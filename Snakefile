from pathlib import Path
configfile: "config.json"
GLOBAL_REF_PATH = config["globalResources"] 
GLOBAL_TMPD_PATH = config["globalTmpdPath"]

#reference_directory = os.path.join(GLOBAL_REF_PATH,config["organism"],config["reference"])
ref_type = list(config["libraries"].values())[0]["reference"]
library = list(config["libraries"].keys())[0]
#TODO how to get organism
reference_path = os.path.join(GLOBAL_REF_PATH, "homo_sapiens", ref_type, "seq", ref_type + ".fa")
#basecaller_location = os.path.join(GLOBAL_TMPD_PATH, "ont-guppy-cpu/bin/guppy_basecaller")
basecaller_location = os.path.join(GLOBAL_TMPD_PATH, "/tmp/dorado-0.5.3-linux-x64/bin/dorado")

library_name = list(config["libraries"].keys())[0]
sample_names = []
for sample in library_name["samples"]:
    sample_name = sample["sample_name"]
    sample_names.append(sample_name)

barcode_flag = config["libraries"][library_name]["samples"][sample_names[0]]["i7_name"]
if(barcode_flag == "NON_BARCODED"):
    is_barcoded = False
else:
    raise Exception("Not implemented yet")

include: "rules/rules.smk"

rule all:
    input:
        expand("{library_name}/outputs/{sample_name}/basecalling/reads_merged.bam", library_name = library_name, sample_name = sample_names)
        #expand("{folder}/conversion_ok.txt", folder = config["run_dir"])
        #expand("{library_name}/sv_calling/variants.vcf", library_name = library)
        #expand("outputs/alignment/{library_name}/minimap2/reads-align.genome.sorted.bam", library_name = config["run_dir"])
        #expand("outputs/basecalling/{library_name}/guppy/sequencing_summary.txt", library_name = config["run_dir"])
