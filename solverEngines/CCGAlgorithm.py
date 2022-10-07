# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 14:55:45 2022

@author: mehdi
"""
import numpy as np
import time
from .masterProblem import masterProblem
from .Subproblem_SeparationProblem import SeparationProblem
from .Subproblem_OriginalProblem import OriginalProblem

class CCGAlgorithm:
    def __init__(self, Param):
        self.Param = Param
        
        self.iteration = 0
        self.run_time = 0
        self.subproblem_run_time = 0
        self.number_subproblem_callbacks = 0
        self.subproblem_callback_run_time = 0
        self.number_uncertaintySet_check = 0
        self.uncertaintySet_check_run_time = 0
        
        self.upper_bound = np.inf
        self.lower_bound = -np.inf
        self.z_sol = {self.iteration: {l: 0.0 for l in Param.InputData.coordinates}}
        self.tornado_head = {}
        self.tornado_tail = {}
        
        self.master_problem = masterProblem(Param)
        self.master_problem.generate_column(self.iteration)
        self.master_problem.generate_constraints(self.iteration, self.z_sol)
        
        self.infeasible_collection = []
        self.feasible_collection = []
        
        
    def run(self):
        self.run_time = time.time()
        while self.upper_bound - self.lower_bound > 0.01 * self.lower_bound and time.time() - self.run_time < 3600:
            self.master_problem.solve()
            self.master_problem.get_solutions(self.iteration)
            self.lower_bound = self.master_problem.master_sol_dict['iteration{}'.format(self.iteration)]['optimal_value']
            self.f_sol = self.master_problem.master_sol_dict['iteration{}'.format(self.iteration)]['f_sol']
            
            subproblem_run_time_begin = time.time()
            if self.Param.subproblem_method == 'Decomposed':
                self.subproblem = SeparationProblem(self.Param, self.f_sol, self.infeasible_collection, self.feasible_collection)
                self.subproblem.solve()
                self.subproblem.get_solutions()
                self.upper_bound = min(self.upper_bound, self.subproblem.subproblem_sol_dict['optimal_value'] 
                                       + sum(self.Param.InputData.first_stage_dislocation[l][s]*self.f_sol[(l,s)] for l in self.Param.InputData.first_stage_dislocation.keys() for s in self.Param.InputData.retrofitting_strategies))
                
                self.infeasible_collection = self.subproblem.infeasible_collection
                self.feasible_collection = self.subproblem.feasible_collection

                self.number_uncertaintySet_check += self.subproblem.number_uncertaintySet_check
                self.uncertaintySet_check_run_time += self.subproblem.uncertaintySet_check_run_time
            
            
            if self.Param.subproblem_method == 'Original':
                self.subproblem = OriginalProblem(self.Param, self.f_sol)
                self.subproblem.solve()
                self.subproblem.get_solutions()
                self.upper_bound = min(self.upper_bound, self.subproblem.subproblem_sol_dict['optimal_value'] 
                                       + sum(self.Param.InputData.first_stage_dislocation[l][s]*self.f_sol[(l,s)] for l in self.Param.InputData.first_stage_dislocation.keys() for s in self.Param.InputData.retrofitting_strategies))
            
            
            
            self.subproblem_run_time += time.time() - subproblem_run_time_begin
            self.number_subproblem_callbacks += self.subproblem.number_callbacks
            self.subproblem_callback_run_time += self.subproblem.callback_run_time

            self.iteration += 1
            self.z_sol[self.iteration] = self.subproblem.subproblem_sol_dict['z_sol']
            self.tornado_head[self.iteration] = self.subproblem.head
            self.tornado_tail[self.iteration] = self.subproblem.tail
            
            self.master_problem.generate_column(self.iteration)
            self.master_problem.generate_constraints(self.iteration, self.z_sol)
            
        self.run_time = time.time() - self.run_time