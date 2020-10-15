#!/usr/bin/env python
# coding: utf-8

# In[185]:


import networkx as nx
import pytest as pt
import pylint as pl
import matplotlib as plt
import scipy
taille_kmer =3
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
    for i in read_fastq("data\eva71_hundred_reads.fq."):
        print(i)
        for j in cut_kmer(i,taille_kmer):
            print(j,end='')
        break
        
def build_kmer_dict(fichier_fastq, jk):
    """Cette fonction va permettre de calculer les occurrences de
    chaque Kmers contenus au sein des reads issus du fastq.
    """
    liste_reads = []
    for sequence in read_fastq("data\eva71_hundred_reads.fq."):
        liste_reads.append(sequence)
    occurrence_kmers = {}
    for read in liste_reads:
        for kmer in cut_kmer(read, taille_kmer):
            if kmer in occurrence_kmers.keys():
                occurrence_kmers[kmer] += 1
            else:
                occurrence_kmers[kmer] = 1
    return occurrence_kmers



def build_graph(occurrence_kmers):
    """Cette fonction va permettre de créer un digraph qui permettra,
    à terme, d'aligner les reads.
    """
    graph = nx.DiGraph()
    for kmer, poids in dico_kmers.items():
        graph.add_edge(kmer[:-1], kmer[1:], weight=poids)
    #nx.draw(graph)
    return graph

build_kmer_dict(read_fastq("data\eva71_two_reads.fq."), 3)
build_graph(build_kmer_dict(read_fastq("data\eva71_two_reads.fq."), 3))

def get_starting_nodes(graph):
    """Fonction qui permet de relever les noeuds d'entrée."""
    noeuds_entree = []
    for noeud in graph.nodes:
        if len(list(graph.predecessors(noeud))) == 0:
            noeuds_entree.append(noeud)
    return noeuds_entree

def get_sink_nodes(graph):
    """Fonction qui permet de relever les noeuds de sortie."""
    noeuds_sortie = []
    for noeud in graph.nodes:
        if len(list(graph.successors(noeud))) == 0:
            noeuds_sortie.append(noeud)
    return noeuds_sortie

def get_contigs(graph, noeuds_entree, noeuds_sortie):
    """Fonction permettant de générer une liste de tulpes
    contenant les contigs associés à leur taille.
    """
    contigs = []
    for noeud_depart in noeuds_entree:
        for noeud_fin in noeuds_sortie:
            for path in nx.all_simple_paths(graph,            source=noeud_depart, target=noeud_fin):
                prep_contig = path
                contig_ecrit = prep_contig[0]
                for i in range(1, len(prep_contig)):
                    contig_ecrit += prep_contig[i][-1:]
                contigs.append((contig_ecrit, len(contig_ecrit)))
    return contigs

get_contigs(graph, get_starting_nodes(graph), get_sink_nodes(graph))



# In[ ]:




