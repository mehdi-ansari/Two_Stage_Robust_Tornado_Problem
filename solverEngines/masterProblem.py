# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 14:25:17 2022

@author: mehdi
"""
import gurobipy as gb
from gurobipy import GRB
import numpy as np
import random

class masterProblem:
    def __init__(self, Param):
        self.Param = Param
        self.model = gb.Model("master_model")
        
        #variables
        self.location_indx = list(Param.InputData.first_stage_dislocation.keys())
        self.retrofit_indx = list(Param.InputData.retrofitting_strategies)
        self.recovery_indx = list(Param.InputData.recovery_strategies)
        
        self.f_var = self.model.addVars(self.location_indx, self.retrofit_indx, vtype=GRB.BINARY, name = "_f")
        self.theta = self.model.addVar(lb = -GRB.INFINITY, name  = "_Theta")
        
        self.r_var = {}
        
        #only one retrofitting strategy can be selected
        self.model.addConstrs(gb.quicksum(self.f_var[l,s] for s in self.retrofit_indx) == 1 
                              for l in self.location_indx)
        
        #objective function
        obj_expr = gb.quicksum(Param.InputData.first_stage_dislocation[l][s]*self.f_var[l,s]
                               for l in self.location_indx for s in self.retrofit_indx)
        self.model.setObjective(obj_expr + self.theta, GRB.MINIMIZE)
        
        
        self.master_sol_dict = {}
        
        #self.fix_random_solution()        
        self.fix_solution()
        
        
    def generate_column(self, iteration):
        self.r_var[iteration] = self.model.addVars(self.location_indx, self.retrofit_indx, self.recovery_indx, vtype=GRB.BINARY, name= "_r{i}".format(i=iteration))
        
        
        
    def generate_constraints(self, iteration, z_sol):
        self.model.addConstr(self.theta >= gb.quicksum(self.Param.InputData.second_stage_dislocation[l][s][p]*z_sol[iteration][l]*self.r_var[iteration][l,s,p]
                                                       for l in self.location_indx for s in self.retrofit_indx for p in self.recovery_indx))
        
        #budget constraint
        self.model.addConstr(gb.quicksum(gb.quicksum(self.Param.InputData.cost_recovery[l][s][p]*self.r_var[iteration][l,s,p] for p in self.recovery_indx)
                                         + self.Param.InputData.cost_retrofitting[l][s] * self.f_var[l,s] for l in self.location_indx for s in self.retrofit_indx)
                             <= self.Param.budget)
        
        #corresponding recovery plan can get 1 value
        self.model.addConstrs(gb.quicksum(self.r_var[iteration][l,s,p] for p in self.recovery_indx) == self.f_var[l,s]
                              for l in self.location_indx for s in self.retrofit_indx)
        
        
    
    def solve(self):
        self.model.read(str(self.Param.ROOT_DIR)+'/solverEngines/master_parameter_set.prm')
        self.model.write(str(self.Param.ROOT_DIR)+'/Results/master.lp')
        self.model.optimize()
        
    
    
    def get_solutions(self, iteration):
        self.model.write(str(self.Param.ROOT_DIR)+'/Results/master.sol')
        self.master_sol_dict['iteration{}'.format(iteration)] = {}
        self.master_sol_dict['iteration{}'.format(iteration)]['optimal_value'] = self.model.objVal
        self.master_sol_dict['iteration{}'.format(iteration)]['theta'] = self.theta.x
        self.master_sol_dict['iteration{}'.format(iteration)]['f_sol'] = self.model.getAttr('x', self.f_var)
        for key, value in self.r_var.items():
            self.master_sol_dict['iteration{}'.format(iteration)]['r_sol[{}]'.format(key)] = self.model.getAttr('x', value)
            
            
    def fix_random_solution(self):
        budget_to_retrofit = 0.6 * self.Param.budget
        self.model.addConstr(gb.quicksum(self.Param.InputData.cost_retrofitting[l][s] * self.f_var[l,s] for l in self.location_indx for s in self.retrofit_indx)
                             <= budget_to_retrofit)
        

        loc_ret = []
        for loc in self.location_indx:
            for ret in self.retrofit_indx:
                loc_ret.append((loc, ret))
        random.shuffle(loc_ret)
        
        fixed_index = []
        for pair in loc_ret:
            if pair[0] not in fixed_index and pair[1] != 'Do_nothing':
                if self.Param.InputData.cost_retrofitting[pair[0]][pair[1]] < budget_to_retrofit:
                    budget_to_retrofit -= self.Param.InputData.cost_retrofitting[pair[0]][pair[1]]
                    fixed_index.append(pair[0])
                    
                    self.model.addConstr(self.f_var[pair[0], pair[1]] == 1)
    
    
    def fix_solution(self):
        constrList = {l:{s:0 for s in self.retrofit_indx} for l in self.location_indx}

        for l in self.location_indx:
            for s in self.retrofit_indx:
                constrList[l][s] = self.model.addConstr(self.f_var[l, s] <= 0, "fixed"+str(l)+"_"+str(s))
            constrList[l]['Do_nothing'].rhs = 1

        solPair_15M = ((43,'R3'), (9,'R3'), (2,'R3'), (54,'R3'), (21,'R3'), (8,'R3'), (40,'R3'), (46,'R3'), (32,'R3'), (34,'R3'),
                   (55,'R3'), (75,'R3'), (78,'R3'), (72,'R3'), (56,'R3'), (92,'R3'), (80,'R3'), (53,'R3'), (22,'R3'), (27,'R3'),
                   (89,'R3'), (5,'R3'))
        
        for pair in solPair_15M:
            constrList[pair[0]][pair[1]].rhs = 1
            self.model.addConstr(self.f_var[pair[0], pair[1]] == 1)