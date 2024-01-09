import os
import Bio
from Bio import Entrez,SeqIO

Entrez.email = "saiponnapalli@utexas.edu"

with open("/work2/07475/vagheesh/stampede2/forOthers/forBIO321G/humanSeq.txt", "r") as f:
    with open("/work2/09269/saip/stampede2/projects/project3/humanSeq.fasta", "w") as out:
        for line in f:
            parts = line.strip().split("\t")
            accession = parts[0]
            population = parts[1]
            handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
            record = SeqIO.read(handle, "fasta")
            out.write(f">{record.id}\n{record.seq}\n")
