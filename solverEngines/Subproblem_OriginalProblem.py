# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 15:19:44 2022

@author: mehdi
"""
import gurobipy as gb
from gurobipy import GRB
import numpy as np
import time
from .Subproblem import Subproblem
from .RecoveryProblem import RecoveryProblem

class OriginalProblem(Subproblem):
    def __init__(self, Param, f_sol):
        super().__init__(Param, f_sol)
        self.subproblem_sol_dict = {}
        self.number_callbacks = 0
        self.callback_run_time = 0
        
        self.UncertaintySet.add_uncertainty_constraint(self.model, self.z_var)
        
        #intiail
        Recovery0 = RecoveryProblem(self.Param, {l: 0.0 for l in Param.InputData.coordinates}, self.f_sol)
        Recovery0.solve()
        Recovery0.get_solutions()
        r_sol = Recovery0.recovery_sol_dict['r_sol']
        self.generate_constraint(r_sol)
        
    def solve(self):        
        self.model.read(str(self.Param.ROOT_DIR)+'/solverEngines/Subproblem_parameter_set.prm')
        self.model.write(str(self.Param.ROOT_DIR)+'/Results/subprblem.lp')
        self.model.params.lazyConstraints = 1
        self.model.params.NonConvex = 2
        def callback(model, where):
            if where == GRB.Callback.MIPSOL:
                self.update_model()
        self.model.optimize(callback)
        
    def get_solutions(self):
        self.model.write(str(self.Param.ROOT_DIR)+'/Results/subproblem.sol')
        self.subproblem_sol_dict['optimal_value'] = self.model.objVal
        self.subproblem_sol_dict['z_sol'] = self.model.getAttr('x', self.z_var)
        
    def update_model(self):
        '''
        This add lazy cuts to the subproblem whenever it finds new feasible MIP solution.

        we try to find a new upper bound from recovery problem. 

        Returns
        -------
        None.

        '''
        self.number_callbacks += 1
        callback_run_time_begin = time.time()
        
        eta_sol = self.model.cbGetSolution(self.eta)
        z_sol = self.model.cbGetSolution(self.z_var)
        
        Recovery = RecoveryProblem(self.Param, z_sol, self.f_sol)
        Recovery.solve()
        Recovery.get_solutions()
        r_sol = Recovery.recovery_sol_dict['r_sol']
        if eta_sol - Recovery.recovery_sol_dict['optimal_value'] > 0.01 * Recovery.recovery_sol_dict['optimal_value']:
            self.model.cbLazy(self.eta <= gb.quicksum(self.Param.InputData.second_stage_dislocation[l][s][p] * self.z_var[l] * r_sol[(l,s,p)]
                                        for l in self.location_indx for s in self.retrofit_indx for p in self.recovery_indx))
        
        
        self.callback_run_time += time.time() - callback_run_time_begin