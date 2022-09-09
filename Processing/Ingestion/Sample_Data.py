# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 19:27:58 2022

@author: mehdi
"""
import numpy as np
import pandas as pd

class SampleData:
    
    def __init__(self, number_of_clusters=None):
        self.coordinates={}
        self.first_stage_dislocation={}
        self.cost_retrofitting={}
        self.second_stage_dislocation={}
        self.cost_recovery={}
        
'''
#Reading Generated data
instance_name = '50buil_instance5.xlsx'
instance = pd.read_excel(instance_name)
coord_type = [('x_coord', float), ('y_coord', float)]
buildings = np.array([(instance['x_coord'][i], instance['y_coord'][i]) for i in range(len(instance))],  dtype=coord_type)
B = np.sort(buildings, order='x_coord')  

#parameters
S = [1,2,3]
L = [1,2]

num_b = len(B)
num_s = len(S)
num_l = len(L)

w = [[0.1*b + 0.2*s + 1 for s in range(num_s)] for b in range(num_b)]
g = [[[0.3*b + 0.3*s + 0.4*l + 1 for l in range(num_l)] for s in range(num_s)] for b in range(num_b)]
c = [[(b%4+1) + 1/(l+1) for l in range(num_l)] for b in range(num_b)]
d = [[(b%4+1) + 1/(s+1) for s in range(num_s)] for b in range(num_b)]

budget = (np.sum(c, axis=0)[1] + np.sum(d, axis=0)[2] + np.sum(c, axis=0)[0] + np.sum(d, axis=0)[0])/2

delta = 1
length = 7
'''