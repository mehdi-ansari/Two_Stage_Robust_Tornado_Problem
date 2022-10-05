# -*- coding: utf-8 -*-
"""
Created on Sun Sep 11 13:14:48 2022

@author: mehdi
"""
import gurobipy as gb
from gurobipy import GRB
import numpy as np
import time
from .Subproblem import Subproblem
from .StabbingLineAlgorithm import StabbingLine
from .LineSegment import Segment
from .RecoveryProblem import RecoveryProblem

class SeparationProblem(Subproblem):
    def __init__(self, Param, f_sol, infeasible_collection, feasible_collection):
        super().__init__(Param, f_sol)
        
        self.number_callbacks = 0
        self.callback_run_time = 0
        self.number_uncertaintySet_check = 0
        self.uncertaintySet_check_run_time = 0
        self.subproblem_sol_dict = {}
        
        self.feasible_collection = feasible_collection
        self.infeasible_collection = infeasible_collection
        self.add_infeasible_collection_cuts()
        
        #intiail
        Recovery0 = RecoveryProblem(self.Param, {l: 0.0 for l in Param.InputData.coordinates}, self.f_sol)
        Recovery0.solve()
        Recovery0.get_solutions()
        r_sol = Recovery0.recovery_sol_dict['r_sol']
        self.generate_constraint(r_sol)       
        
        #self.test_this_solution()
    '''def test_this_solution(self):
        for l in [32,56,55,8,78,22,54,75,72,53,2,89,43,9,40,46,5]:
            self.model.addConstr(self.z_var[l] == 1)'''

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
        '''
        This add lazy cuts to the subproblem whenever it finds new feasible MIP solution.
        We check:
            1. If we already know this solution is feasible
            2. If not known, is there a stabbing line intersecting all z_sol = 1, If not cut off
            3. If there is line, is there a line segment; if there is a line segment we are done
            4. If there is not line segment, we call non convex feasible space to investigate more!
        
        If we conclude there is segment to make z_sol, we try to find a new upper bound from recovery problem. 

        Returns
        -------
        None.

        '''
        self.number_callbacks += 1
        callback_run_time_begin = time.time()
        
        eta_sol = self.model.cbGetSolution(self.eta)
        z_sol = self.model.cbGetSolution(self.z_var)
        
        damaged_location_coordindates = {}
        for l, soluation in z_sol.items():
            if soluation >= 0.5:
                damaged_location_coordindates[l] = self.Param.InputData.coordinates[l]
         
        isFeasible = True
        #We don't need to check the feasibility if already known
        do_check = True
        for collection in self.feasible_collection:
            if set(collection['locations']) == set(damaged_location_coordindates.keys()):
                do_check = False
                isFeasible = True
                self.head = collection['tornado_head']
                self.tail = collection['tornado_tail']
                
        
        if do_check == True and len(damaged_location_coordindates)>0:
            stabbing_line_algorithm = StabbingLine(damaged_location_coordindates, np.ones(len(damaged_location_coordindates)), self.Param.width, self.Param.length+2*self.Param.width)
            line_intersecting_maximal_circles = stabbing_line_algorithm.find_line_intersecting_maximal_circles()

            #if there is no line crossing through all cirlces whose centers are active locations with radius delta:
            if len(damaged_location_coordindates) > line_intersecting_maximal_circles['max_stabbed_weighted_circles'] + 1:
                isFeasible = False
                self.model.cbLazy(gb.quicksum(self.z_var[loc] for loc in damaged_location_coordindates.keys()) <= line_intersecting_maximal_circles['max_stabbed_weighted_circles'] + 1)
                self.infeasible_collection.append({'locations': damaged_location_coordindates.keys(), 'maximal_intersections': line_intersecting_maximal_circles['max_stabbed_weighted_circles'] + 1})
            
            #else, if there is a line, check is there a line segment intersecting all circles too?
            elif len(damaged_location_coordindates) == line_intersecting_maximal_circles['max_stabbed_weighted_circles'] + 1:
                '''lineSegment = Segment(line_intersecting_maximal_circles, self.Param, damaged_location_coordindates)

                if lineSegment.isSegment():
                    isFeasible = True
                    self.head = lineSegment.intersection1
                    self.tail = lineSegment.intersection2
                else:'''
                self.number_uncertaintySet_check += 1
                uncertaintySet_check_run_time_begin = time.time()
                check_uncertainty_set = self.UncertaintySet.check_feasibility(damaged_location_coordindates)
                self.uncertaintySet_check_run_time += time.time() - uncertaintySet_check_run_time_begin
                if check_uncertainty_set['status'] == 2:
                    isFeasible = True
                    self.head = check_uncertainty_set['head']
                    self.tail = check_uncertainty_set['tail']
                else:
                    isFeasible = False
                    self.model.cbLazy(gb.quicksum(self.z_var[loc] for loc in damaged_location_coordindates.keys()) <= len(damaged_location_coordindates)-1)
                    self.infeasible_collection.append({'locations': damaged_location_coordindates.keys(), 'maximal_intersections': len(damaged_location_coordindates)-1})
                        
        
        ##if solution is feasible, investigate the new upperbound
        if isFeasible:
            self.feasible_collection.append({'locations': damaged_location_coordindates.keys(), 'tornado_head':self.head, 'tornado_tail': self.tail})
            
            Recovery = RecoveryProblem(self.Param, z_sol, self.f_sol)
            Recovery.solve()
            Recovery.get_solutions()
            r_sol = Recovery.recovery_sol_dict['r_sol']
            if eta_sol - Recovery.recovery_sol_dict['optimal_value'] > 0.01 * Recovery.recovery_sol_dict['optimal_value']:
                self.model.cbLazy(self.eta <= gb.quicksum(self.Param.InputData.second_stage_dislocation[l][s][p] * self.z_var[l] * r_sol[(l,s,p)]
                                            for l in self.location_indx for s in self.retrofit_indx for p in self.recovery_indx))
        
        
        self.callback_run_time += time.time() - callback_run_time_begin