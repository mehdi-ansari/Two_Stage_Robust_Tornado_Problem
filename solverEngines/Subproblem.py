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
        
        #objective function
        self.model.setObjective(self.eta, GRB.MAXIMIZE)
        
        
        #add cuts
        self.add_infeasible_pair_cuts()
        self.add_infeasible_triple_cuts()
        # self.add_infeasible_quadruple_cuts()   #It is not ready to use!
        
        self.head = []
        self.tail = []
        
        #self.UncertaintySet.add_uncertainty_constraint(self.model, self.z_var)
        #self.fix_solution()
        
    def add_infeasible_pair_cuts(self):
        for pair in self.ValidCuts.infeasible_pair:
            self.model.addConstr(self.z_var[pair[0]] + self.z_var[pair[1]] <= 1)
                
    def add_infeasible_triple_cuts(self):
        for trip in self.ValidCuts.infeasible_tri:
            self.model.addConstr(self.z_var[trip[0]] + self.z_var[trip[1]] + self.z_var[trip[2]] <= 2)
        
    def add_infeasible_quadruple_cuts(self):
        for quad in self.infeasible_quadruple:
            self.model.addConstr(self.z_var[quad[0]] + self.z_var[quad[1]] + self.z_var[quad[2]] + self.z_var[quad[3]] <= 3)   
    
    
    def generate_constraint(self, r_sol):
        self.model.addConstr(self.eta <= gb.quicksum(self.Param.InputData.second_stage_dislocation[l][s][p] * self.z_var[l] * r_sol[(l,s,p)]
                                            for l in self.location_indx for s in self.retrofit_indx for p in self.recovery_indx))
        
    def fix_solution(self):
        #sol = [32,5,92,21,56,80,55,8,78,34,22,
               #54,75,72,53,2,89,43,9,40,27,46]
   
        #sol = [-1,70,-1,85,71,-1,38,21,41,-1,80,62,-1,10,55,36,44,48,8,-1,61,78,34,-1,-1,64,79,1,22,-1,54,82,11,75,-1,19,99,72,47,49,77,-1,98,53,65,63,14,23,93,-1,2,81,45,88,89,25,-1,43,76,42,7,17,37,26,9,60,16,83,13,67,-1,58,4,68,87,40,30,27,73,84,35,46,39,51,52,6,96,74,-1,29,18,24,94,33,12,57,97,69,3,66]
        
        #for s in sol:
            #if s != -1:
                #self.model.addConstr(self.z_var[s] == 0)
        self.model.addConstr(self.z_var[17] == 1)
        self.model.addConstr(self.z_var[66] == 1)