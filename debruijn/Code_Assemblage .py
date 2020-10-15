#!/usr/bin/env python
# coding: utf-8

# In[85]:


import networkx as nx
import pytest as pt
import pylint as pl

def read_fastq(fichier : str):
    with open(fichier) as filin:
        for line in enumerate(filin):
            yield next(filin)[:-1]
            next(filin)
            next(filin)

def cut_kmer(seq : str,k : int):
    #Cette fonction va permettre d'obtenir des k-mers de taille désirée à partir d'un read renseigné.
    for i in range(len(seq)-k+1):
        yield seq[i:i+k]
        
if __name__=="__main__" :
    for i in read_fastq("data\eva71_plus_perfect.fq"):
        print(i)
        for j in cut_kmer(i,k=21):
            print(j,end='')
        break
        
def build_kmer_dict(fichier_fastq, taille_kmer):
    """Cette fonction va permettre de calculer les occurrences de
    chaque Kmers contenus au sein des reads issus du fastq.
    """
    liste_reads = []
    for sequence in read_fastq("data\eva71_plus_perfect.fq"):
        liste_reads.append(sequence)
    occurrence_kmers = {}
    for read in liste_reads:
        for kmer in cut_kmer(read, taille_kmer):
            if kmer in occurrence_kmers.keys():
                occurrence_kmers[kmer] += 1
            else:
                occurrence_kmers[kmer] = 1
    return occurrence_kmers


def build_graph(dico_kmers):
    """Cette fonction va permettre de créer un digraph qui permettra,
    à terme, d'aligner les reads.
    """
    dico_kmers = build_kmer_dict(read_fastq("data\eva71_plus_perfect.fq"), 21)
    graph = nx.DiGraph()
    for kmer, poids in dico_kmers.items():
        graph.add_edge(kmer[:-1], kmer[1:], weight=poids)
        return graph
build_kmer_dict(read_fastq("data\eva71_plus_perfect.fq"), 21)    


# In[ ]:




