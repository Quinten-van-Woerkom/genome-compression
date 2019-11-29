# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 09:17:23 2019

@author: Tygo
"""

def diagonal(seq1, seq2,offset):
    
    n = seq1.size
    m = seq2.size
    
    
    if offset < n:
        cordinate = [n-1-offset,0]
        
    elif offset < m + n - 1:
        cordinate = [0,offset - (n-1)]
        
    else: 
        print('offset is too high')
        return
     
        
    length = min([n-cordinate[0],m-cordinate[1]])
    
    if length <= 0:
        return []
    
    
    diagonal = []
    for i in range(length):
        
        diagonal.append( seq1[i+cordinate[0]] != seq2[i+cordinate[1]] )
  
    return diagonal


