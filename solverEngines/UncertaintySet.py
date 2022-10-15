# -*- coding: utf-8 -*-
"""
Created on Thu Aug 19 17:10:37 2021

@author: mehdi
"""
import numpy as np
import matplotlib.pyplot as plt
import gurobipy as gb
from gurobipy import GRB
import numpy as np

def Euclidean_distance(a,b):
    return np.sqrt((a[0]-b[0])**2+(a[1]-b[1])**2)

class UncertaintySet:
    def __init__(self, Param):
        self.Param = Param
        
    def check_feasibility(self, damaged_location_coordindates):
        model = gb.Model("Uncertainy")    

        active_index = damaged_location_coordindates.keys()
        coordinates = damaged_location_coordindates
            
        
        #variables    
        t = model.addVars(active_index, lb=0, ub=1, vtype=GRB.CONTINUOUS, name = "t")
        v_x = model.addVars(active_index, lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name = "v_x")
        v_y = model.addVars(active_index, lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name = "v_y")
        head_x = model.addVar(lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name = "head_x")
        head_y = model.addVar(lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name = "head_y")
        tail_x = model.addVar(lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name = "tail_x")
        tail_y = model.addVar(lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name = "tail_y")
        
        #constraints
        for l in active_index:
            model.addConstr((1-t[l])*head_x + t[l]*tail_x - coordinates[l][0] == v_x[l])
            model.addConstr((1-t[l])*head_y + t[l]*tail_y - coordinates[l][1] == v_y[l])
            model.addConstr(v_x[l]*v_x[l] + v_y[l]*v_y[l] <= self.Param.width**2)
            
        model.addConstr(head_x*head_x - 2*head_x*tail_x + tail_x*tail_x +
                        head_y*head_y - 2*head_y*tail_y + tail_y*tail_y <= self.Param.length**2)
    
        #Auxiliary Constraints
            #find extreme locations:
        d_sides = 0
        side1 = 0
        side2 = 0
        for p1 in active_index:
            for p2 in active_index:
                if d_sides < Euclidean_distance(coordinates[p1], coordinates[p2]):
                    d_sides = Euclidean_distance(coordinates[p1], coordinates[p2])
                    side1 = p1
                    side2 = p2
        
            #Confined the space
        model.addConstr(head_x >= coordinates[side1][0] - self.Param.width)
        model.addConstr(head_x <= coordinates[side1][0] + self.Param.width)
        model.addConstr(head_y >= coordinates[side1][1] - self.Param.width)
        model.addConstr(head_y <= coordinates[side1][1] + self.Param.width)
        
        model.addConstr(tail_x >= coordinates[side2][0] - self.Param.width)
        model.addConstr(tail_x <= coordinates[side2][0] + self.Param.width)
        model.addConstr(tail_y >= coordinates[side2][1] - self.Param.width)
        model.addConstr(tail_y <= coordinates[side2][1] + self.Param.width)
        
        #solve
        model.params.NonConvex = 2
        model.params.TimeLimit = 900
        model.optimize()
        #model.write('check feasibility.lp')
        #model.write('check feasibility.sol')
        
        
        if model.Status == 2:
            head = [head_x.x, head_y.x]
            tail = [tail_x.x, tail_y.x]
        else: 
            print("model.Status:", model.Status)
            head = None
            tail = None
        
        return {'status': model.Status, 'head': head, 'tail': tail}
        

    def add_uncertainty_constraint(self, model, z_var):
        coordinates = self.Param.InputData.coordinates
        
        BigM = self.compute_BigM(coordinates)
        
        t = model.addVars(coordinates.keys(), lb=0, ub=1, vtype=GRB.CONTINUOUS, name = "t")
        v_x = model.addVars(coordinates.keys(), lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name = "v_x")
        v_y = model.addVars(coordinates.keys(), lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name = "v_y")
        head_x = model.addVar(lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name = "head_x")
        head_y = model.addVar(lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name = "head_y")
        tail_x = model.addVar(lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name = "tail_x")
        tail_y = model.addVar(lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name = "tail_y")
        
        '''x_rand = np.random.uniform(-5164, (-5154-5164)/2)
        y_rand = np.random.uniform(2553, 2563)
        angle = np.random.uniform(-np.pi/2, np.pi/2)
        model.addConstr(head_x == x_rand)
        model.addConstr(head_y == y_rand)
        model.addConstr(tail_x >= x_rand + 4.9*np.cos(angle))
        model.addConstr(tail_y >= y_rand + 4.9*np.sin(angle))
        model.params.NonConvex = 2
        model.addConstr(head_x == -94.567 * 54.6)
        model.addConstr(head_y == 37.06 * 69)
        model.addConstr(tail_x == -94.479 * 54.6)
        model.addConstr(tail_y == 37.078 * 69)
        model.params.NonConvex = 2'''
        #constraints
        for l in coordinates.keys():
            model.addConstr((1-t[l])*head_x + t[l]*tail_x - coordinates[l][0] == v_x[l])
            model.addConstr((1-t[l])*head_y + t[l]*tail_y - coordinates[l][1] == v_y[l])
            model.addConstr(v_x[l]*v_x[l] + v_y[l]*v_y[l] <= self.Param.width**2 + (2*self.Param.width*BigM[l] + BigM[l]**2)*(1-z_var[l]))
            
        model.addConstr(head_x*head_x - 2*head_x*tail_x + tail_x*tail_x +
                        head_y*head_y - 2*head_y*tail_y + tail_y*tail_y <= self.Param.length**2)
        
    
    def compute_BigM(self, coordinates):
        max_x = -np.inf
        min_x = np.inf
        max_y = -np.inf
        min_y = np.inf
        
        for indx, coord in coordinates.items():
            if max_x < coord[0]:
                max_x = coord[0]
            if min_x > coord[0]:
                min_x = coord[0]
            
            if max_y < coord[1]:
                max_y = coord[1]
            if min_y > coord[1]:
                min_y = coord[1]

        
        bigMs = {}
        for i,b in coordinates.items():
            dist1 = Euclidean_distance(b, (max_x, max_y))
            dist2 = Euclidean_distance(b, (max_x, min_y))
            dist3 = Euclidean_distance(b, (min_x, max_y))
            dist4 = Euclidean_distance(b, (min_x, min_y))
            
            bigMs[i] = max(dist1, dist2, dist3, dist4)
        
        return bigMs