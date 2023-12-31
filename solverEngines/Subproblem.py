# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 22:33:19 2022

@author: mehdi
"""
import gurobipy as gb
from gurobipy import GRB
import numpy as np
from .ValidCuts import valid_cuts_sets
from .UncertaintySet import UncertaintySet

class Subproblem:
    def __init__(self, Param, f_sol):
        self.ValidCuts = valid_cuts_sets(Param)
        self.UncertaintySet = UncertaintySet(Param)
        self.Param = Param
        self.f_sol = f_sol
        self.model = gb.Model("subproblem_model")
        
        #variables
        self.location_indx = list(Param.InputData.first_stage_dislocation.keys())
        self.retrofit_indx = list(Param.InputData.retrofitting_strategies)
        self.recovery_indx = list(Param.InputData.recovery_strategies)
        
        self.eta = self.model.addVar(lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name = "_eta")
        self.z_var = self.model.addVars(self.location_indx, vtype=GRB.BINARY, name= "_z")
        
        self.fix_sol([56,74])
        #objective function
        self.model.setObjective(self.eta, GRB.MAXIMIZE)
        
        
        #add cuts
        self.add_infeasible_pair_cuts()
        self.add_infeasible_triple_cuts()
        # self.add_infeasible_quadruple_cuts()   #It is not ready to use!
        
        self.head = []
        self.tail = []
        
        
    def add_infeasible_pair_cuts(self):
        for pair in self.ValidCuts.infeasible_pair:
            self.model.addConstr(self.z_var[pair[0]] + self.z_var[pair[1]] <= 1)
                
    def add_infeasible_triple_cuts(self):
        numHigherPopulated = int(0.2 * len(self.location_indx))
        sorted_locs = sorted(self.Param.InputData.second_stage_dislocation.items(), key=lambda x: x[1]['Do_nothing']['Do_nothing'], reverse=True)[:numHigherPopulated]
        listLocs = [l[0] for l in sorted_locs]

        for trip in self.ValidCuts.infeasible_tri:
            if (trip[0] in listLocs) or (trip[1] in listLocs) or (trip[2] in listLocs):
                self.model.addConstr(self.z_var[trip[0]] + self.z_var[trip[1]] + self.z_var[trip[2]] <= 2)
        
        '''for trip in self.ValidCuts.infeasible_tri:
            self.model.addConstr(self.z_var[trip[0]] + self.z_var[trip[1]] + self.z_var[trip[2]] <= 2)'''
        
    def add_infeasible_quadruple_cuts(self):
        for quad in self.infeasible_quadruple:
            self.model.addConstr(self.z_var[quad[0]] + self.z_var[quad[1]] + self.z_var[quad[2]] + self.z_var[quad[3]] <= 3)   
    
    
    def generate_constraint(self, r_sol):
        self.model.addConstr(self.eta <= gb.quicksum(self.Param.InputData.second_stage_dislocation[l][s][p] * self.z_var[l] * r_sol[(l,s,p)]
                                            for l in self.location_indx for s in self.retrofit_indx for p in self.recovery_indx))
        
    def fix_sol(self, indxList):
        for i in indxList:
            self.model.addConstr(self.z_var[i] == 1)