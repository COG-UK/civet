from Bio import Phylo
from Bio import SeqIO
import csv
import collections

config["patient_stems"] = config["patient_str"].split(",")

rule all:
    input:
        expand(os.path.join(config["outdir"],"report","figures","genome_graph_{patient}.png"), patient=config["patient_stems"]),
        os.path.join(config["tempdir"],"gather_prompt.txt")


rule get_names_and_seqs:
    input:
        patient_file = os.path.join(config["tempdir"],"patient_data","{patient}_taxa.txt"),
        aligned_query_seqs = config["aligned_query_seqs"],
        background_seqs = config["background_seqs"],
        outgroup_fasta = config["outgroup_fasta"]
    output:
        aln = os.path.join(config["tempdir"], "seqs_for_snps","{patient}.fasta"),
        new_names = os.path.join(config["tempdir"], "seqs_for_snps","{patient}_names.txt")
    run:
        taxa = []
        taxa.append("Reference\n")
        with open(input.patient_file) as f:
            next(f)
            for l in f:
                taxa.append(l.strip("\n"))

        with open(output.new_names, 'w') as fw:
            fw.write("name,label\n")
            fw.write("Reference,Reference\n")
            # for name in name_maps:
                # if name in taxa:
            for name in taxa:
                fw.write(f"{name},{name}\n") #add in display name stuff

        with open(output.aln, "w") as fw:
            for record in SeqIO.parse(input.outgroup_fasta, "fasta"):
                fw.write(f">Reference\n{record.seq}\n")

            for record in SeqIO.parse(input.aligned_query_seqs, "fasta"):
                if record.id in taxa:
                    fw.write(f">{record.description}\n{record.seq}\n")

            for record in SeqIO.parse(input.background_seqs,"fasta"):
                if record.id in taxa:
                    fw.write(f">{record.description}\n{record.seq}\n")


rule make_snp_figure:
    input:
        aln = rules.get_names_and_seqs.output.aln,
        names = rules.get_names_and_seqs.output.new_names
    params:
        out_stem = os.path.join(config["outdir"],"report","figures","genome_graph_{patient}")
    output:
        os.path.join(config["outdir"],"report","figures","genome_graph_{patient}.png")
    shell:
        """
        snipit {input.aln:q} -r "Reference" -o {params.out_stem} -l {input.names}
        """

rule gather_graphs:
    input:
        expand(os.path.join(config["outdir"],"report","figures","genome_graph_{patient}.png"), patient=config["patient_stems"])
    output:
        os.path.join(config["tempdir"],"gather_prompt.txt")
    shell:
        "touch {output}"
    