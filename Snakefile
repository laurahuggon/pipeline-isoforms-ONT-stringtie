import os
import Bio
import shutil
from os import path
from pathlib import Path

from snakemake.utils import validate
from snakemake.utils import min_version
min_version("7.18")

# ----------------------------------------------------------------

configfile: "config.yml"
validate(config, schema="schema/config_schema.yaml")
workdir: path.join(config["workdir"], config["sample_name"])

WORKDIR = config["workdir"]
SNAKEDIR = path.dirname(workflow.snakefile)

sample = config["sample_name"]
sample_dir = path.join(WORKDIR, sample)


shutil.copy2(path.join(SNAKEDIR, "config.yml"), sample_dir)


in_fastq = config["reads_fastq"]

# Make sure the path is absolute
if not path.isabs(in_fastq):
    in_fastq = path.join(SNAKEDIR, in_fastq)
    assert os.path.exists(in_fastq)

# If it is a directory, gather all FASTQ files inside, otherwise treat it as a file
if os.path.isdir(in_fastq):
    in_fastq = [str(f) for f in Path(in_fastq).rglob('*') if f.suffix in ['.fastq', '.fq', '.gz']]
else:
    in_fastq = [in_fastq]  # Wrap the single file into a list

in_genome = config["genome"]
if not path.isabs(in_genome):
    in_genome = path.join(SNAKEDIR, in_genome)
    assert os.path.exists(in_genome)

in_annotation = config["annotation"]
if not path.isabs(in_annotation) and in_annotation != "":
    in_annotation = path.join(SNAKEDIR, in_annotation)

# ----------------------------------------------------------------

target_list = [
    path.join(sample_dir, "StringTie", f"{sample}_stringtie.gff"),
    path.join(sample_dir, "Mapping", f"{sample}_alignment_stats.tsv"),
]

rule all:
    input:
        target_list

# ----------------------------------------------------------------

rule concatenate_reads:
    input:
        fq = in_fastq

    output:
        fq_concat = temp(path.join(sample_dir, "processed_reads", f"{sample}_reads_concat.fq"))

    params:
        concat = "True"

    threads: config["threads"]

    shell:
        """
        pwd
        
        mkdir -p processed_reads temp_concat
        
        ls .
        
        echo "Starting unzipping..."

        if [ "{params.concat}" = "True" ]; then
            # Handle gzipped files
            find $(dirname {input.fq[0]}) -type f -name "*.gz" | \
                xargs -P {threads} -I{{}} sh -c 'zcat "{{}}" > temp_concat/$(basename "{{}}" .gz)'
                
            ls temp_concat/*.fastq | wc -l

            echo "Starting concatenation..."
            
            # Concatenate all decompressed temp files
            cat temp_concat/* > {output.fq_concat}
            
            echo "Finished concatenation"

            # Remove temp files
            rm -r temp_concat
            
            ls .
            sleep 5m

            # Concatenate non-gz files directly
            find $(dirname {input.fq[0]}) -type f -not -name "*.gz" -exec cat {{}} + >> {output.fq_concat}
        else
            # Handle non-concatenation case
            if [[ "{input.fq}" == *.gz ]]; then
                # Decompress gzipped file
                zcat {input.fq} > {output.fq_concat}
            else
                # Create a symbolic link
                ln -s $(realpath {input.fq}) {output.fq_concat}
            fi
        fi
        """

# ----------------------------------------------------------------

# # Rule for NanoStat, using concatenated FASTQ files
# rule nanostat:
#     input:
#         fq = rules.concatenate_reads.output.fq_concat if config.get("concatenate_reads", False) else in_fastq

#     output:
#         ns = "Nanostat/stat_out.txt"

#     threads: config["threads"]

#     shell:
#         """
#         NanoStat -n {output.ns} -t {threads} --tsv --fastq {input.fq}
#         """

# ----------------------------------------------------------------

rule pychopper:
    input:
        fq = rules.concatenate_reads.output.fq_concat if config.get("concatenate_reads", False) else in_fastq
    output:
        pyfq = path.join(sample_dir, "{}/Pychopper".format(sample_dir), f"{sample}_full_length_reads.fq")
    params:
        outpath = "{}/Pychopper".format(sample_dir),  # Use this for output path directly
        prefix = f"{sample}_full_length_reads.fq",
        pc = "True" if config["run_pychopper"] else "False",
        pc_opts = config["pychopper_opts"],
        kit = config["kit"]
    log: path.join(sample_dir, "{}/Pychopper".format(sample_dir), f"{sample}_pychopped.log")
    threads: config["threads"]

    run:

        shell("mkdir -p {params.outpath}")

        if params.pc == "True":
            shell(
                """
                pychopper {params.pc_opts} -t {threads} -k {params.kit} \
                -r {params.outpath}/report.pdf \
                -u {params.outpath}/unclassified.fq \
                -w {params.outpath}/rescued.fq \
                {input.fq} {output.pyfq} 2> {log}
                """
            )
        else:
            shell(
                """
                ln -s `realpath {input.fq}` {output.pyfq}
                """
            )


# ----------------------------------------------------------------

rule minimap_mapping:
    input:
        genome = config["genome"],
        fq = rules.pychopper.output.pyfq

    output:
        sam = path.join(sample_dir, "Mapping", f"{sample}_minimap.sam")

    params:
        opts = config["minimap2_opts"]

    log: path.join(sample_dir, "Mapping", f"{sample}_minimap.log")

    threads: config["threads"]

    shell:
        """
        minimap2 -ax splice {params.opts} --secondary=no -t {threads} {input.genome} {input.fq} \
        | samtools sort -@ {threads} -o {output.sam}
        """

# ----------------------------------------------------------------

rule aln_stats:
    input:
        sam=rules.minimap_mapping.output.sam

    output:
        tsv=path.join(sample_dir, "Mapping", f"{sample}_alignment_stats.tsv")

    threads: config["threads"]

    shell:
        """
        samtools stats -@ {threads} {input.sam} > {output.tsv}
        """

# ----------------------------------------------------------------

str_threads = 10
if config["threads"] < str_threads:
    str_threads = config["threads"]

rule run_stringtie:
    input:
        sam=rules.minimap_mapping.output.sam,
        fa=config["genome"]

    output:
        gff = path.join(sample_dir, "StringTie", f"{sample}_stringtie.gff")

    params:
        opts = config["stringtie_opts"],
        ann = in_annotation if config["use_guide_annotation"] else ""

    log: path.join(sample_dir, "StringTie", f"{sample}_StringTie.log")

    threads: config["threads"]

    run:
        g_flag = f"-G {params.ann}" if params.ann else ""
        shell(
            f"""
            stringtie --rf {g_flag} -L -v -p {threads} {params.opts} -o {output.gff} {input.sam} 2> {log}
            """
        )