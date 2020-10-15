#!/usr/bin/env python
# coding: utf-8

# In[6]:


def read_fastq(fichier : str):
    with open(fichier) as filin:
        for line in enumerate(filin):
            yield next(filin)[:-1]
            next(filin)
            next(filin)

def cut_kmer(seq : str,k : int):
    for i in range(len(seq)-k+1):
        yield seq[i:i+k]
        
if __name__=="__main__" :
    for i in read_fastq("data\eva71_plus_perfect.fq"):
        print(i)
        for j in cut_kmer(i,k=3):
            print(j,end='')
        break

