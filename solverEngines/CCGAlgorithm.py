# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 14:55:45 2022

@author: mehdi
"""
from .masterProblem import masterProblem
from .Subproblem import Subproblem

class CCGAlgorithm:
    def __init__(self, Param):
        self.master_problem = masterProblem(Param)
        self.subproblem = Subproblem(Param)
        



