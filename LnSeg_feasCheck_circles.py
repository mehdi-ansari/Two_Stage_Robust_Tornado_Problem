# -*- coding: utf-8 -*-
"""
Created on Thu Aug 19 17:10:37 2021

@author: mehdi
"""
import numpy as np
import matplotlib.pyplot as plt
import gurobipy as gb
from gurobipy import GRB
import Stabbing_line as stab_line

def Euclidean_distance(a,b):
    return np.sqrt((a[0]-b[0])**2+(a[1]-b[1])**2)

#Feasibility Check-----------------++++++++++++++++++++------------------------
#Alternative 1
def feasibility_check(B, delta, length, active_index):
    d_sides = 0
    side1 = 0
    side2 = 0
    for p1 in active_index:
        for p2 in active_index:
            if d_sides < Euclidean_distance(B[p1], B[p2]):
                d_sides = Euclidean_distance(B[p1], B[p2])
                side1 = p1
                side2 = p2
    
    
    num_b = len(active_index)
    
    model = gb.Model("feasibility_check")
    
    t_b = model.addVars(num_b, lb=0, ub=1, vtype=GRB.CONTINUOUS, name = "t")
    v_b = model.addVars(num_b, 2, lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name = "v")
    head = model.addVars(2, lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name = "head")
    tail = model.addVars(2, lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name = "tail")
    
    for b in range(num_b):
        model.addConstr((1-t_b[b])*head[0] + t_b[b]*tail[0] - B[active_index[b]][0] == v_b[b,0])
        model.addConstr((1-t_b[b])*head[1] + t_b[b]*tail[1] - B[active_index[b]][1] == v_b[b,1])
        model.addConstr(v_b[b,0]*v_b[b,0] + v_b[b,1]*v_b[b,1] <= delta**2)
        
    model.addConstr(head[0]*head[0] - 2*head[0]*tail[0] + tail[0]*tail[0] +
                    head[1]*head[1] - 2*head[1]*tail[1] + tail[1]*tail[1] <= length**2)
    
    #Confined the space
    model.addConstr(head[0] >= B[side1][0] - delta)
    model.addConstr(head[0] <= B[side1][0] + delta)
    model.addConstr(head[1] >= B[side1][1] - delta)
    model.addConstr(head[1] <= B[side1][1] + delta)
    
    model.addConstr(tail[0] >= B[side2][0] - delta)
    model.addConstr(tail[0] <= B[side2][0] + delta)
    model.addConstr(tail[1] >= B[side2][1] - delta)
    model.addConstr(tail[1] <= B[side2][1] + delta)
    
    '''p1 = stab_line.internal_tangent(B[side1], B[side2], delta)[0]
    p2 = stab_line.internal_tangent(B[side1], B[side2], delta)[1]
    
    q1 = stab_line.angle(p1)
    q2 = stab_line.angle(p2)
    
    angles = np.sort([-1/np.tan(q1), -1/np.tan(q2)])
    model.addConstr(tail[1] - head[1] >= angles[0] * (tail[0] - head[0]))
    model.addConstr(tail[1] - head[1] <= angles[1] * (tail[0] - head[0]))'''
    
    model.params.NonConvex = 2
    model.params.TimeLimit = 1800
    #model.Params.Threads = 16
        
    model.optimize()
    #model.write("feasibility_check.lp")
    head = []
    tail = []
    if model.Status == 2:
        head.append(model.getVarByName("head[0]").x)
        head.append(model.getVarByName("head[1]").x)
        tail.append(model.getVarByName("tail[0]").x)
        tail.append(model.getVarByName("tail[1]").x)
    else: 
        print("model.Status:", model.Status)    
    
    return model.Status, head, tail


#Alternative 2
def line_segement_feasibility(B, delta, length, active_index, side1_index, side2_index):
    indx1 = active_index.index(side1_index)
    indx2 = active_index.index(side2_index)
    
    print(active_index)
    print(side1_index)
    print(side2_index)
    print(indx1)
    print(indx2)
    #input()

    model = gb.Model("feasibility_w/circles")
    
    center = model.addVars(2, len(active_index), lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name = "c")
    slope = model.addVar(lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name = "m")
    bb = model.addVar(lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name = "b")
    
    for b in range(len(active_index)):
        model.addConstr((B[active_index[b]][0] - center[0,b])*(B[active_index[b]][0] - center[0,b])
                        + (B[active_index[b]][1] - center[1,b])*(B[active_index[b]][1] - center[1,b])
                        <= delta**2)
        
        model.addConstr(center[1,b] == slope*(center[0,b] - center[0,0]) + center[1,0])
    
    model.addConstr((center[0,indx1] - center[0,indx2])*(center[0,indx1] - center[0,indx2])
                        + (center[1,indx1] - center[1,indx2])*(center[1,indx1] - center[1,indx2])
                        <= length**2)

    
    model.params.NonConvex = 2
    model.params.TimeLimit = 600
        
    model.optimize()
    #model.write("feasibility_with_circles.lp")
    head = []
    tail = []
    if model.Status == 2:
        head = [model.getVarByName("c[0,{i}]".format(i=indx1)).x, model.getVarByName("c[1,{i}]".format(i=indx1)).x]
        tail = [model.getVarByName("c[0,{i}]".format(i=indx2)).x, model.getVarByName("c[1,{i}]".format(i=indx2)).x]
    else: 
        print("model.Status:", model.Status)

    return model.Status, head, tail
        
        
'''B = [(0,0), (5,0), (5,1), (10,0)]
active_index = [0,2,3]
delta = 1
length = 8
side1_index = 0
side2_index = 3

print(line_segement_feasibility(B, delta, length, active_index, side1_index, side2_index))'''