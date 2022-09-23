# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 14:55:45 2022

@author: mehdi
"""
from .masterProblem import masterProblem
from .Subproblem_SeparationProblem import SeparationProblem

class CCGAlgorithm:
    def __init__(self, Param):
        self.master_problem = masterProblem(Param)
        
        infeasible_collection = []
        feasible_collection = []
        self.subproblem = SeparationProblem(Param, infeasible_collection, feasible_collection)
        



