# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 21:15:11 2019

@author: Tygo
"""

import numpy as np

L = 4
E = 1
I = 100
p = .25
show = True


M = int(np.floor((L+1)/(E+1)))

checks = int(I//M)+1

interval = (int(I//M)+1) * M


match= np.random.uniform(0,1,I)>p # 1 means mismatch, 0 means match
print(len(match))
#match = np.random.randint(0,2,I)
tail = np.zeros([interval - I],dtype = int)
match = np.append(match,tail)


indexes = M*np.linspace(1,checks-E,checks-E,dtype = int)-1
all_points = M*np.linspace(1,checks,checks,dtype = int)-1
checkpoints = match[all_points]

if show:
    print("The stepsize is: "+str(M))
    print("The interval has been enlarged from " + str(I)+" to "+str(interval))
    print("A randomly generated boolean sequence:")
    print(" ")
    print(match)
    print(" ")
    print("The indexes that are checked are: "+str(all_points))
    print(" ")
    print("The boolean values at these points are:")
    print(" ")
    print(checkpoints)
    print(" ")


sumcheckpoints = np.zeros([checks-E],dtype = int)

Range = np.linspace(0,M*E,E+1,dtype = int)
print(Range)
computations = 0 
for Q in range(M):
    n = 0
    i = 0
    looplength = len(indexes)
    while n < looplength: 
        
        sumcheckpoints[i] += sum(match[indexes[i]+Range-Q])
        
        if sumcheckpoints[i] >= E:
            
            sumcheckpoints = np.delete(sumcheckpoints, i)
            indexes = np.delete(indexes, i)
            
        else:
            i += 1
        n += 1
        computations += 1
    
    
    exceed = sumcheckpoints >= E
    if show:
        print(indexes)
        print(sumcheckpoints)
        print(exceed)


cover = np.linspace(0,M-1,M,dtype = int)
for i in indexes:
    print(match(i+cover))

print(computations)

