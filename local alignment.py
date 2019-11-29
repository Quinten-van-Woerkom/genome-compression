# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 09:12:50 2019

@author: Tygo
"""
from Sparse_SW_for_local_alignment import sparse_SW
from diag_seq_for_local_alignment import diagonal
import numpy as np
#import matplotlib.pyplot as plt

n = 15
m = 15

L = 8
E = 3

seq1 = np.random.randint(0,2,n)
seq2 = np.random.randint(0,2,m)


counter  = 0 


while counter < m + n - 2 * L + 2-1:
    seq = diagonal(seq1,seq2,counter+L-1)
    
    
    
    print(sparse_SW(seq,L,E, False))
    print( [int(i) for i in seq])
    
    
    counter += 1






