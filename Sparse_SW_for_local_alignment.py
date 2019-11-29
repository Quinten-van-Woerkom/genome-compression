# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 09:53:18 2019

@author: Tygo
"""

import numpy as np

def sparse_SW (seq, L, E, show):
    

    I = len(seq)
    
    
    M = int(np.floor((L+1)/(E+1)))
    
    checks = int(np.ceil(I/M))
    
    interval = (checks) * M
    
   
    
    
    for i in range(interval - I):
        seq.append(True)
    
    
    indexes = []
    
    
    for i in range(checks-E):
        indexes.append(M*i-1+M)
    
    
    
    
    coverage = []
    
    for i in range(E):
        coverage.append(M*i)
    
    
    sumcheckpoints = np.zeros(checks-E, dtype = int)
    print("The stepsize is: "+str(M)+ ", the minimal length L has been reduced to " + str((E+1)*M-1))
    if show:

        print("The interval has been enlarged from " + str(I)+" to "+str(interval))
        print("The number of checks are " + str(checks))
        print("The indexes are " + str(indexes))
        print(sumcheckpoints)
    

    deleted = 0
    Q = 0
    while Q < M:
        J = 0
        if show:
            print("We are shifted by " + str(Q))
        while J < checks-E - deleted:
            
            if show:
                print("The checkpoints give " +str([ seq[M*J-Q+i] for i in coverage ]))
            
            sumcheckpoints[J] += sum([ seq[M*J-Q+i] for i in coverage ])
            
            if show:
                print(sumcheckpoints)
            
            if  sumcheckpoints[J] >= E:
                
                #delete elements if too large
                deleted += 1
                
                if show:
                    print("The element on location "+str(indexes[J])+ ", number "+str(J)+ " in the list, is deleted as it has too many errors.")
                    
                    
                sumcheckpoints = np.delete(sumcheckpoints,J)
                indexes = np.delete(indexes,J)
                
                if show:
                   print("The indexes are now " + str(indexes))
                
                
            else:
                J += 1
                if show:
                    print("J is now " +str(J))
            
            
            
            
        
        
        
            
        Q += 1     
    
    return indexes
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    return 