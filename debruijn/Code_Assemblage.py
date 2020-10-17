#!/usr/bin/env python
# coding: utf-8

# In[13]:


import networkx as nx
import pytest as pt
import pylint as pl
import matplotlib as plt
import scipy
import os
import statistics
import argparse




def read_fastq(fichier_fastq):
    with open(fichier_fastq,"r") as fastq:
        for line in enumerate(fastq):
            yield next(fastq)[:-1]
            next(fastq)
            next(fastq)
            


def cut_kmer(seq, taille_k):
    #Cette fonction va permettre d'obtenir des k-mers de taille désirée à partir d'un read renseigné.
    for i in range(len(seq)-taille_k+1):
        yield seq[i:i+taille_k]
        
if __name__=="__main__" :
    for i in read_fastq("data\eva71_hundred_reads.fq."):
        print(i)
        for j in cut_kmer(i,taille_k):
            print(j,end='')
        break
        
def build_kmer_dict(fichier_fastq,taille_k):
    """Cette fonction va permettre de calculer les occurrences de
    chaque Kmers contenus au sein des reads issus du fastq.
    """
    liste_reads = []
    for sequence in fichier_fastq:
        liste_reads.append(sequence)
    occurrence_kmers = {}
    for read in liste_reads:
        for kmer in cut_kmer(read, taille_k):
            if kmer in occurrence_kmers.keys():
                occurrence_kmers[kmer] += 1
            else:
                occurrence_kmers[kmer] = 1
    return occurrence_kmers



def build_graph(dict_kmers):
    """Cette fonction va permettre de créer un digraph qui permettra,
    à terme, d'aligner les reads.
    """
    graph = nx.DiGraph()
    for kmer, poids in dict_kmers.items():
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

#noeuds_entree = get_starting_nodes(graph)
                                   
def get_sink_nodes(graph):
    """Fonction qui permet de relever les noeuds de sortie."""
    noeuds_sortie = []
    for noeud in graph.nodes:
        if len(list(graph.successors(noeud))) == 0:
            noeuds_sortie.append(noeud)
    return noeuds_sortie
                                   
#noeuds_sortie = get_sink_nodes(graph)
                                   
def get_contigs(graph, start, end):
    """Fonction permettant de générer une liste de tulpes
    contenant les contigs associés à leur taille.
    """
    contigs = []
    for noeud_depart in start:
        for noeud_fin in end:
            for path in nx.all_simple_paths(graph,source=noeud_depart, target=noeud_fin):
                prep_contig = path
                contig_ecrit = prep_contig[0]
                for i in range(1, len(prep_contig)):
                    contig_ecrit += prep_contig[i][-1:]
                contigs.append((contig_ecrit, len(contig_ecrit)))
    return contigs


#get_contigs(graph,noeuds_entree,noeuds_sortie)

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

#save_contigs(get_contigs(graph,noeuds_entree,noeuds_sortie),"Contigs.fna")

def std(liste_valeurs):
    """Calcul l'écart-type de la liste de valeurs"""
    return statistics.stdev(liste_valeurs)


def path_average_weight(graph, chemin):
    """Cette fonction permet de retourner le poids moyen d'un
    chemin.
    """
    poids = 0
    nbre_edges = 0
    edges = graph.subgraph(chemin).edges(data=True)
    #On va récupérer les poids de chaque lien.
    for u_value, v_value, e_value in edges:
        poids += e_value['weight']
        nbre_edges += 1
    poids = poids/nbre_edges
    return poids



    
    
fichier_fastq = read_fastq("data\eva71_hundred_reads.fq.")
taille_k=21
   

liste_reads = []
for sequence in fichier_fastq:
    liste_reads.append(sequence)
#print (liste_reads)

liste_kmers = []
for read in liste_reads:
    for kmers in cut_kmer(read, taille_k):
        liste_kmers.append(kmers)
#print(liste_kmers)
occurrence_kmers = build_kmer_dict(fichier_fastq, taille_k)
#print(occurrence_kmers)

graph=build_graph(build_kmer_dict(read_fastq("data\eva71_two_reads.fq."), taille_kmer))
nx.draw(graph)

debuts = get_starting_nodes(graph)
fins = get_sink_nodes(graph)
liste_contigs = get_contigs(graph, debuts, fins)
   

    


# In[ ]:





# In[ ]:




