# -*- coding: utf-8 -*-
"""
Created on Tue Dec  7 15:38:29 2021

@author: mehdi
"""
import numpy as np
from scipy.stats import norm

Archetype1 = [[(3.47,0.14), (3.55,0.13), (3.62,0.12), (3.64,0.14)],
              [(3.86,0.14), (3.93,0.13), (4.01,0.12), (4.05,0.12)],
              [(3.87,0.15), (3.95,0.13), (4.03,0.12), (4.24,0.12)]]


def DS(wind_speed, median, stand_deviation):
    return norm.cdf((np.log(wind_speed)-median)/stand_deviation)

p_Relocation = []
cost = []

v = 60
for s in range(len(Archetype1)):
    N = 1 - DS(v, Archetype1[s][0][0],Archetype1[s][0][1])
    I = DS(v, Archetype1[s][0][0],Archetype1[s][0][1]) - DS(v, Archetype1[s][1][0],Archetype1[s][1][1])
    M = DS(v, Archetype1[s][1][0],Archetype1[s][1][1]) - DS(v, Archetype1[s][2][0],Archetype1[s][2][1])
    H = DS(v, Archetype1[s][2][0],Archetype1[s][2][1]) - DS(v, Archetype1[s][3][0],Archetype1[s][3][1])
    C = DS(v, Archetype1[s][3][0],Archetype1[s][3][1])

    p_Relocation.append([0*N, 0.1*I, 0.6*M, 0.9*H, 1*C])
    cost.append([0*N, 0.005*I, 0.023*M, 0.117*H, 0.234*C])
    
print(p_Relocation)
print(cost)



