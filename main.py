from Bio import SeqIO
from Bio.SeqUtils import GC

filename = "C:/Users/48512/Desktop/sequence.fasta"
seq = SeqIO.read(filename, "fasta")
print(len(seq))
