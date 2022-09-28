# -*- coding: utf-8 -*-
"""
Created on Sun Sep 11 12:36:18 2022

@author: mehdi
"""
import gurobipy as gb
from gurobipy import GRB
import numpy as np

class RecoveryProblem:
    def __init__(self, Param, z_sol, f_sol):
        self.Param = Param
        self.model = gb.Model("recoveryProblem_model")
        
        #variables
        self.location_indx = list(Param.InputData.first_stage_dislocation.keys())
        self.retrofit_indx = list(Param.InputData.retrofitting_strategies)
        self.recovery_indx = list(Param.InputData.recovery_strategies)
        
        self.r_var = self.model.addVars(self.location_indx, self.retrofit_indx, self.recovery_indx, vtype=GRB.BINARY, name = "_r")
        
        #budget constraint
        self.model.addConstr(gb.quicksum(gb.quicksum(self.Param.InputData.cost_recovery[l][s][p]*self.r_var[l,s,p] for p in self.recovery_indx)
                                         + self.Param.InputData.cost_retrofitting[l][s] * f_sol[(l,s)] for l in self.location_indx for s in self.retrofit_indx)
                             <= self.Param.budget)
        
        #corresponding recovery plan can get 1 value
        self.model.addConstrs(gb.quicksum(self.r_var[l,s,p] for p in self.recovery_indx) == f_sol[(l,s)]
                              for l in self.location_indx for s in self.retrofit_indx)
        
        #objective function
        self.model.setObjective(gb.quicksum(self.Param.InputData.second_stage_dislocation[l][s][p] * z_sol[l] * self.r_var[l,s,p]
                                            for l in self.location_indx for s in self.retrofit_indx for p in self.recovery_indx), 
                                GRB.MINIMIZE)
        
        
        self.recovery_sol_dict = {}
        
        
    def solve(self):
        self.model.read(str(self.Param.ROOT_DIR)+'/solverEngines/RecoveryP_parameter_set.prm')
        self.model.write(str(self.Param.ROOT_DIR)+'/Results/recovery.lp')
        self.model.optimize()
        
        
    def get_solutions(self):
        self.model.write(str(self.Param.ROOT_DIR)+'/Results/recovery.sol')
        self.recovery_sol_dict['optimal_value'] = self.model.objVal
        self.recovery_sol_dict['r_sol'] = self.model.getAttr('x', self.r_var)