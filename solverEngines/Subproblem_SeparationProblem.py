# -*- coding: utf-8 -*-
"""
Created on Sun Sep 11 13:14:48 2022

@author: mehdi
"""
import gurobipy as gb
from gurobipy import GRB
import numpy as np
from .Subproblem import Subproblem
from .UncertaintySet import UncertaintySet
from .StabbingLineAlgorithm import StabbingLine

class SeparationProblem(Subproblem):
    def __init__(self, Param, infeasible_collection, feasible_collection):
        super().__init__(Param)
        
        self.subproblem_sol_dict = {}
        
        self.feasible_collection = feasible_collection
        self.infeasible_collection = infeasible_collection
        self.add_infeasible_collection_cuts()

    
    def add_infeasible_collection_cuts(self):
        '''
        Infeasible collections from previous iterations (user search cuts)

        Returns
        -------
        None.

        '''
        for collection in self.infeasible_collection:
            self.model.addConstr(gb.quicksum(self.z_var[l] for l in collection['locations']) <= collection['maximal_intersections'])
            
        
    
    def solve(self):
        self.model.read(str(self.Param.ROOT_DIR)+'/solverEngines/Subproblem_parameter_set.prm')
        self.model.write(str(self.Param.ROOT_DIR)+'/Results/subprblem.lp')
        self.model.params.lazyConstraints = 1
        def callback(model, where):
            if where == GRB.Callback.MIPSOL:
                self.update_model()
        self.model.optimize(callback)
        
    def get_solutions(self):
        self.model.write(str(self.Param.ROOT_DIR)+'/Results/subproblem.sol')
        self.subproblem_sol_dict['optimal_value'] = self.model.objVal
        self.subproblem_sol_dict['z_sol'] = self.model.getAttr('x', self.z_var)

    def update_model(self):
        eta_sol = self.model.cbGetSolution(self.eta)
        z_sol = self.model.cbGetSolution(self.z_var)
        
        active_index = []
        damaged_location_coordindates = {}
        for l, soluation in z_sol.items():
            if soluation >= 0.5:
                damaged_location_coordindates[l] = self.Param.InputData.coordinates[l]
           
        #We don't need to check the feasibility if already known
        do_check = True
        for collection in self.feasible_collection:
            if set(collection['indices']) == set(active_index):
                do_check = False
        
        if do_check == True:
            stabbing_line_algorithm = StabbingLine(damaged_location_coordindates, np.ones(len(damaged_location_coordindates)), self.Param.width, self.Param.length+2*self.Param.width)
            line_intersecting_maximal_circles = stabbing_line_algorithm.find_line_intersecting_maximal_circles()
            
            #if there is no line crossing through all cirlces whose centers are active locations with radius delta:
            if len(damaged_location_coordindates) > line_intersecting_maximal_circles['max_stabbed_weighted_circles'] + 1:
                self.model.cbLazy(gb.quicksum(self.z_var[loc] for loc in damaged_location_coordindates.keys()) <= line_intersecting_maximal_circles['max_stabbed_weighted_circles'] + 1)
                self.infeasible_collection.append({'locations': damaged_location_coordindates.keys(), 'maximal_intersections': line_intersecting_maximal_circles['max_stabbed_weighted_circles'] + 1})
                print('infeasible_collection', self.infeasible_collection)
                input()
                
            
            print('I stopped here!')
        