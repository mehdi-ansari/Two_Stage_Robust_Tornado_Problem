# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 22:33:19 2022

@author: mehdi
"""
import gurobipy as gb
from gurobipy import GRB
import numpy as np

class Subproblem:
    def __init__(self, Param):
        self.Param = Param
        self.model = gb.Model("subproblem_model")
        
        #variables
        self.location_indx = list(Param.InputData.first_stage_dislocation.keys())
        self.retrofit_indx = list(Param.InputData.retrofitting_strategies)
        self.recovery_indx = list(Param.InputData.recovery_strategies)
        
        self.eta = self.model.addVar(lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name = "_eta")
        self.z = self.model.addVars(self.location_indx, vtype=GRB.BINARY, name= "_z")
    