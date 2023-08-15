import numpy as np
from Bio import SeqIO
from Bio.Align._aligners import PairwiseAligner
from Bio.Seq import Seq
from Levenshtein import distance
from six import unichr
import matplotlib.pyplot as plt
import numpy as np
#Sequences

filename_sequence_human = 'C:/Users/48512/Desktop/sequence.fasta'
filename_sequence_dog = 'C:/Users/48512/Desktop/sequence_dog.fasta'

sequence_human = SeqIO.read(filename_sequence_human, 'fasta')
sequence_dog = SeqIO.read(filename_sequence_dog, 'fasta')

seq_human = sequence_human.seq
seq_dog = sequence_dog.seq

seq1 = Seq('ACTTAG')
seq2 = Seq('ACTACT')

#print(distance(seq_human,seq_dog))

def delta(x,y):
    return 0 if x == y else 1
def M(seq1,seq2,i,j,k):
    return sum(delta(x,y) for x,y in zip(seq1[i:i+k],seq2[j:j+k]))

def makeMatrix(seq1,seq2,k):
    n = len(seq1)
    m = len(seq2)
    return [[M(seq1,seq2,i,j,k) for j in range(m-k+1)] for i in range(n-k+1)]

def plotMatrix(M,t, seq1, seq2, nonblank = unichr(0x25A0), blank = ' '):
    print(' |' + seq2)
    print('-'*(2 + len(seq2)))
    for label,row in zip(seq1,M):
        line = ''.join(nonblank if s < t else blank for s in row)
        print(label + '|' + line)
def dotplot(seq1,seq2,k = 1,t = 1):
    M = makeMatrix(seq1,seq2,k)
    plotMatrix(M, t, seq1,seq2)

plt.imshow(np.array(makeMatrix(seq_human,seq_dog,1)))
xt=plt.xticks(np.arange(len(list(seq_human))), list(seq_human))
yt = plt.yticks(np.arange(len(list(seq_human))), list(seq_human))
plt.show()
