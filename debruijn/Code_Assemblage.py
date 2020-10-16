#!/usr/bin/env python
# coding: utf-8

# In[33]:


import networkx as nx
import pytest as pt
import pylint as pl
import matplotlib as plt
import scipy
import os

taille_kmer =21

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
        
def build_kmer_dict(fichier_fastq,taille_kmer):
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
    for kmer, poids in occurrence_kmers.items():
        graph.add_edge(kmer[:-1], kmer[1:], weight=poids)
    nx.draw(graph)
    return graph

#build_kmer_dict(read_fastq("data\eva71_two_reads.fq."), taille_kmer)
#graph=build_graph(build_kmer_dict(read_fastq("data\eva71_two_reads.fq."), taille_kmer))

def get_starting_nodes(graph):
    """Fonction qui permet de relever les noeuds d'entrée."""
    noeuds_entree = []
    for noeud in graph.nodes:
        if len(list(graph.predecessors(noeud))) == 0:
            noeuds_entree.append(noeud)
    return noeuds_entree

noeuds_entree = get_starting_nodes(graph)
                                   
def get_sink_nodes(graph):
    """Fonction qui permet de relever les noeuds de sortie."""
    noeuds_sortie = []
    for noeud in graph.nodes:
        if len(list(graph.successors(noeud))) == 0:
            noeuds_sortie.append(noeud)
    return noeuds_sortie
                                   
noeuds_sortie = get_sink_nodes(graph)
                                   
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


get_contigs(graph,noeuds_entree,noeuds_sortie)

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def save_contigs(liste_contigs, nom_fichier):
    """Cette fonction permet d'exporter les contigs dans un fichier
    au format FASTA
    """
    with open(nom_fichier, "w") as fichier_sortie:
        numero = 0
        for contigs in liste_contigs:
            fichier_sortie.write(">contig_{0} len={1}\n".format(numero, contigs[1]))
            fichier_sortie.write("{0}\n".format(fill(contigs[0])))
            numero += 1

save_contigs(get_contigs(graph,noeuds_entree,noeuds_sortie),"Export_contigs.fna")


# In[ ]:




