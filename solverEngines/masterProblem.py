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
        self.fix_optimal_solution()        
        
        
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
        budget_to_retrofit = 1 * self.Param.budget
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
                    
    def fix_optimal_solution(self):
        optimal_solution = {(32, 'Do_nothing'): 0.0, (32, 'R1'): 0.0, (32, 'R2'): 0.0, (32, 'R3'): 1.0, (70, 'Do_nothing'): 1.0, (70, 'R1'): 0.0, (70, 'R2'): 0.0, (70, 'R3'): 0.0, (5, 'Do_nothing'): 0.0, (5, 'R1'): 0.0, (5, 'R2'): 1.0, (5, 'R3'): 0.0, (85, 'Do_nothing'): 1.0, (85, 'R1'): 0.0, (85, 'R2'): 0.0, (85, 'R3'): 0.0, (71, 'Do_nothing'): 1.0, (71, 'R1'): 0.0, (71, 'R2'): 0.0, (71, 'R3'): 0.0, (92, 'Do_nothing'): 0.0, (92, 'R1'): 0.0, (92, 'R2'): 0.0, (92, 'R3'): 1.0, (38, 'Do_nothing'): 1.0, (38, 'R1'): 0.0, (38, 'R2'): 0.0, (38, 'R3'): 0.0, (21, 'Do_nothing'): 0.0, (21, 'R1'): 0.0, (21, 'R2'): 0.0, (21, 'R3'): 1.0, (41, 'Do_nothing'): 1.0, (41, 'R1'): 0.0, (41, 'R2'): 0.0, (41, 'R3'): 0.0, (56, 'Do_nothing'): 0.0, (56, 'R1'): 0.0, (56, 'R2'): 0.0, (56, 'R3'): 1.0, (80, 'Do_nothing'): 0.0, (80, 'R1'): 0.0, (80, 'R2'): 0.0, (80, 'R3'): 1.0, (62, 'Do_nothing'): 1.0, (62, 'R1'): 0.0, (62, 'R2'): 0.0, (62, 'R3'): 0.0, (20, 'Do_nothing'): 1.0, (20, 'R1'): 0.0, (20, 'R2'): 0.0, (20, 'R3'): 0.0, (10, 'Do_nothing'): 1.0, (10, 'R1'): 0.0, (10, 'R2'): 0.0, (10, 'R3'): 0.0, (55, 'Do_nothing'): 0.0, (55, 'R1'): 0.0, (55, 'R2'): 0.0, (55, 'R3'): 1.0, (36, 'Do_nothing'): 1.0, (36, 'R1'): 0.0, (36, 'R2'): 0.0, (36, 'R3'): 0.0, (44, 'Do_nothing'): 1.0, (44, 'R1'): 0.0, (44, 'R2'): 0.0, (44, 'R3'): 0.0, (48, 'Do_nothing'): 1.0, (48, 'R1'): 0.0, (48, 'R2'): 0.0, (48, 'R3'): 0.0, (8, 'Do_nothing'): 0.0, (8, 'R1'): 0.0, (8, 'R2'): 0.0, (8, 'R3'): 1.0, (15, 'Do_nothing'): 1.0, (15, 'R1'): 0.0, (15, 'R2'): 0.0, (15, 'R3'): 0.0, (61, 'Do_nothing'): 1.0, (61, 'R1'): 0.0, (61, 'R2'): 0.0, (61, 'R3'): 0.0, (78, 'Do_nothing'): 0.0, (78, 'R1'): 0.0, (78, 'R2'): 0.0, (78, 'R3'): 1.0, (34, 'Do_nothing'): 0.0, (34, 'R1'): 0.0, (34, 'R2'): 0.0, (34, 'R3'): 1.0, (100, 'Do_nothing'): 1.0, (100, 'R1'): 0.0, (100, 'R2'): 0.0, (100, 'R3'): 0.0, (90, 'Do_nothing'): 1.0, (90, 'R1'): 0.0, (90, 'R2'): 0.0, (90, 'R3'): 0.0, (64, 'Do_nothing'): 1.0, (64, 'R1'): 0.0, (64, 'R2'): 0.0, (64, 'R3'): 0.0, (79, 'Do_nothing'): 1.0, (79, 'R1'): 0.0, (79, 'R2'): 0.0, (79, 'R3'): 0.0, (1, 'Do_nothing'): 1.0, (1, 'R1'): 0.0, (1, 'R2'): 0.0, (1, 'R3'): 0.0, (22, 'Do_nothing'): 0.0, (22, 'R1'): 0.0, (22, 'R2'): 0.0, (22, 'R3'): 1.0, (31, 'Do_nothing'): 1.0, (31, 'R1'): 0.0, (31, 'R2'): 0.0, (31, 'R3'): 0.0, (54, 'Do_nothing'): 0.0, (54, 'R1'): 0.0, (54, 'R2'): 0.0, (54, 'R3'): 1.0, (82, 'Do_nothing'): 1.0, (82, 'R1'): 0.0, (82, 'R2'): 0.0, (82, 'R3'): 0.0, (11, 'Do_nothing'): 1.0, (11, 'R1'): 0.0, (11, 'R2'): 0.0, (11, 'R3'): 0.0, (75, 'Do_nothing'): 0.0, (75, 'R1'): 0.0, (75, 'R2'): 0.0, (75, 'R3'): 1.0, (50, 'Do_nothing'): 1.0, (50, 'R1'): 0.0, (50, 'R2'): 0.0, (50, 'R3'): 0.0, (19, 'Do_nothing'): 1.0, (19, 'R1'): 0.0, (19, 'R2'): 0.0, (19, 'R3'): 0.0, (99, 'Do_nothing'): 1.0, (99, 'R1'): 0.0, (99, 'R2'): 0.0, (99, 'R3'): 0.0, (72, 'Do_nothing'): 0.0, (72, 'R1'): 0.0, (72, 'R2'): 0.0, (72, 'R3'): 1.0, (47, 'Do_nothing'): 1.0, (47, 'R1'): 0.0, (47, 'R2'): 0.0, (47, 'R3'): 0.0, (49, 'Do_nothing'): 1.0, (49, 'R1'): 0.0, (49, 'R2'): 0.0, (49, 'R3'): 0.0, (77, 'Do_nothing'): 1.0, (77, 'R1'): 0.0, (77, 'R2'): 0.0, (77, 'R3'): 0.0, (95, 'Do_nothing'): 1.0, (95, 'R1'): 0.0, (95, 'R2'): 0.0, (95, 'R3'): 0.0, (98, 'Do_nothing'): 1.0, (98, 'R1'): 0.0, (98, 'R2'): 0.0, (98, 'R3'): 0.0, (53, 'Do_nothing'): 0.0, (53, 'R1'): 0.0, (53, 'R2'): 0.0, (53, 'R3'): 1.0, (65, 'Do_nothing'): 1.0, (65, 'R1'): 0.0, (65, 'R2'): 0.0, (65, 'R3'): 0.0, (63, 'Do_nothing'): 1.0, (63, 'R1'): 0.0, (63, 'R2'): 0.0, (63, 'R3'): 0.0, (14, 'Do_nothing'): 1.0, (14, 'R1'): 0.0, (14, 'R2'): 0.0, (14, 'R3'): 0.0, (23, 'Do_nothing'): 1.0, (23, 'R1'): 0.0, (23, 'R2'): 0.0, (23, 'R3'): 0.0, (93, 'Do_nothing'): 1.0, (93, 'R1'): 0.0, (93, 'R2'): 0.0, (93, 'R3'): 0.0, (28, 'Do_nothing'): 1.0, (28, 'R1'): 0.0, (28, 'R2'): 0.0, (28, 'R3'): 0.0, (2, 'Do_nothing'): 0.0, (2, 'R1'): 0.0, (2, 'R2'): 0.0, (2, 'R3'): 1.0, (81, 'Do_nothing'): 1.0, (81, 'R1'): 0.0, (81, 'R2'): 0.0, (81, 'R3'): 0.0, (45, 'Do_nothing'): 1.0, (45, 'R1'): 0.0, (45, 'R2'): 0.0, (45, 'R3'): 0.0, (88, 'Do_nothing'): 1.0, (88, 'R1'): 0.0, (88, 'R2'): 0.0, (88, 'R3'): 0.0, (89, 'Do_nothing'): 0.0, (89, 'R1'): 0.0, (89, 'R2'): 0.0, (89, 'R3'): 1.0, (25, 'Do_nothing'): 1.0, (25, 'R1'): 0.0, (25, 'R2'): 0.0, (25, 'R3'): 0.0, (91, 'Do_nothing'): 1.0, (91, 'R1'): 0.0, (91, 'R2'): 0.0, (91, 'R3'): 0.0, (43, 'Do_nothing'): 0.0, (43, 'R1'): 0.0, (43, 'R2'): 0.0, (43, 'R3'): 1.0, (76, 'Do_nothing'): 1.0, (76, 'R1'): 0.0, (76, 'R2'): 0.0, (76, 'R3'): 0.0, (42, 'Do_nothing'): 1.0, (42, 'R1'): 0.0, (42, 'R2'): 0.0, (42, 'R3'): 0.0, (7, 'Do_nothing'): 1.0, (7, 'R1'): 0.0, (7, 'R2'): 0.0, (7, 'R3'): 0.0, (17, 'Do_nothing'): 1.0, (17, 'R1'): 0.0, (17, 'R2'): 0.0, (17, 'R3'): 0.0, (37, 'Do_nothing'): 1.0, (37, 'R1'): 0.0, (37, 'R2'): 0.0, (37, 'R3'): 0.0, (26, 'Do_nothing'): 1.0, (26, 'R1'): 0.0, (26, 'R2'): 0.0, (26, 'R3'): 0.0, (9, 'Do_nothing'): 0.0, (9, 'R1'): 0.0, (9, 'R2'): 0.0, (9, 'R3'): 1.0, (60, 'Do_nothing'): 1.0, (60, 'R1'): 0.0, (60, 'R2'): 0.0, (60, 'R3'): 0.0, (16, 'Do_nothing'): 1.0, (16, 'R1'): 0.0, (16, 'R2'): 0.0, (16, 'R3'): 0.0, (83, 'Do_nothing'): 1.0, (83, 'R1'): 0.0, (83, 'R2'): 0.0, (83, 'R3'): 0.0, (13, 'Do_nothing'): 1.0, (13, 'R1'): 0.0, (13, 'R2'): 0.0, (13, 'R3'): 0.0, (67, 'Do_nothing'): 1.0, (67, 'R1'): 0.0, (67, 'R2'): 0.0, (67, 'R3'): 0.0, (59, 'Do_nothing'): 1.0, (59, 'R1'): 0.0, (59, 'R2'): 0.0, (59, 'R3'): 0.0, (58, 'Do_nothing'): 1.0, (58, 'R1'): 0.0, (58, 'R2'): 0.0, (58, 'R3'): 0.0, (4, 'Do_nothing'): 1.0, (4, 'R1'): 0.0, (4, 'R2'): 0.0, (4, 'R3'): 0.0, (68, 'Do_nothing'): 1.0, (68, 'R1'): 0.0, (68, 'R2'): 0.0, (68, 'R3'): 0.0, (87, 'Do_nothing'): 1.0, (87, 'R1'): 0.0, (87, 'R2'): 0.0, (87, 'R3'): 0.0, (40, 'Do_nothing'): 0.0, (40, 'R1'): 0.0, (40, 'R2'): 0.0, (40, 'R3'): 1.0, (30, 'Do_nothing'): 1.0, (30, 'R1'): 0.0, (30, 'R2'): 0.0, (30, 'R3'): 0.0, (27, 'Do_nothing'): 0.0, (27, 'R1'): 0.0, (27, 'R2'): 0.0, (27, 'R3'): 1.0, (73, 'Do_nothing'): 1.0, (73, 'R1'): 0.0, (73, 'R2'): 0.0, (73, 'R3'): 0.0, (84, 'Do_nothing'): 1.0, (84, 'R1'): 0.0, (84, 'R2'): 0.0, (84, 'R3'): 0.0, (35, 'Do_nothing'): 1.0, (35, 'R1'): 0.0, (35, 'R2'): 0.0, (35, 'R3'): 0.0, (46, 'Do_nothing'): 0.0, (46, 'R1'): 0.0, (46, 'R2'): 0.0, (46, 'R3'): 1.0, (39, 'Do_nothing'): 1.0, (39, 'R1'): 0.0, (39, 'R2'): 0.0, (39, 'R3'): 0.0, (51, 'Do_nothing'): 1.0, (51, 'R1'): 0.0, (51, 'R2'): 0.0, (51, 'R3'): 0.0, (52, 'Do_nothing'): 1.0, (52, 'R1'): 0.0, (52, 'R2'): 0.0, (52, 'R3'): 0.0, (6, 'Do_nothing'): 0.0, (6, 'R1'): 0.0, (6, 'R2'): 0.0, (6, 'R3'): 1.0, (96, 'Do_nothing'): 0.0, (96, 'R1'): 0.0, (96, 'R2'): 1.0, (96, 'R3'): 0.0, (74, 'Do_nothing'): 1.0, (74, 'R1'): 0.0, (74, 'R2'): 0.0, (74, 'R3'): 0.0, (86, 'Do_nothing'): 1.0, (86, 'R1'): 0.0, (86, 'R2'): 0.0, (86, 'R3'): 0.0, (29, 'Do_nothing'): 1.0, (29, 'R1'): 0.0, (29, 'R2'): 0.0, (29, 'R3'): 0.0, (18, 'Do_nothing'): 1.0, (18, 'R1'): 0.0, (18, 'R2'): 0.0, (18, 'R3'): 0.0, (24, 'Do_nothing'): 1.0, (24, 'R1'): 0.0, (24, 'R2'): 0.0, (24, 'R3'): 0.0, (94, 'Do_nothing'): 1.0, (94, 'R1'): 0.0, (94, 'R2'): 0.0, (94, 'R3'): 0.0, (33, 'Do_nothing'): 1.0, (33, 'R1'): 0.0, (33, 'R2'): 0.0, (33, 'R3'): 0.0, (12, 'Do_nothing'): 1.0, (12, 'R1'): 0.0, (12, 'R2'): 0.0, (12, 'R3'): 0.0, (57, 'Do_nothing'): 1.0, (57, 'R1'): 0.0, (57, 'R2'): 0.0, (57, 'R3'): 0.0, (97, 'Do_nothing'): 1.0, (97, 'R1'): 0.0, (97, 'R2'): 0.0, (97, 'R3'): 0.0, (69, 'Do_nothing'): 1.0, (69, 'R1'): 0.0, (69, 'R2'): 0.0, (69, 'R3'): 0.0, (3, 'Do_nothing'): 1.0, (3, 'R1'): 0.0, (3, 'R2'): 0.0, (3, 'R3'): 0.0, (66, 'Do_nothing'): 1.0, (66, 'R1'): 0.0, (66, 'R2'): 0.0, (66, 'R3'): 0.0}
        
        for i, sol in optimal_solution.items():
            self.model.addConstr(self.f_var[i[0], i[1]] == sol)