# -*- coding: utf-8 -*-
"""
Created on Mon May 24 13:24:49 2021

@author: mehdi
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import gurobipy as gb
from gurobipy import GRB
import time
import Stabbing_line as stab_line
import LnSeg_feasCheck_circles as feas_check
from itertools import combinations



#SUBPROBLEM+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++_________
def Euclidean_distance(a,b):
    return np.sqrt((a[0]-b[0])**2+(a[1]-b[1])**2)

def BigM_generator(sorted_Buildings):
    max_x = sorted_Buildings[len(sorted_Buildings)-1][0]
    min_x = sorted_Buildings[0][0]
    
    ys = [sorted_Buildings[i][1] for i in range(len(sorted_Buildings))]
    max_y = max(ys)
    min_y = min(ys)

    bigMs = []
    for b in sorted_Buildings:
        dist1 = Euclidean_distance(b, (max_x, max_y))
        dist2 = Euclidean_distance(b, (max_x, min_y))
        dist3 = Euclidean_distance(b, (min_x, max_y))
        dist4 = Euclidean_distance(b, (min_x, min_y))
        
        bigMs.append(max(dist1, dist2, dist3, dist4))

    return bigMs

    
class subproblem:
    def __init__(self, B, delta, S, L, c, d, f_opt, g, budget, quad_inf, tri_inf, pair_inf, length, inf_collection, fea_collection, max_inter_circ):
        self.B = B
        self.delta = delta
        self.BigM = BigM_generator(self.B)
        self.S = S
        self.L = L
        self.c = c
        self.d = d
        self.f = f_opt
        self.g = g
        self.A = budget
        self.length = length
        self.inf_collection = inf_collection
        self.fea_collection = fea_collection
        
        self.tri_inf = tri_inf
        self.pair_inf = pair_inf
        self.quad_inf = quad_inf
        
        self.max_inter_circ = max_inter_circ
        
        '''self.slope = 0
        self.intercept = 0'''
        self.head = 0
        self.tail = 0
        
        self.new_inf_collection = []
        self.new_fea_collection = []
        
        self.timelim_reached = 0
        
        self.nonConv_feasCheck_time = 0
        self.feasCheck_time = 0
        self.feasCheck_iteration = 0
    #__________________________________________________________________________
    def upper_problem_bigM(self):
        B = self.B
        delta = self.delta
        BigM = self.BigM
        S = self.S
        L = self.L
        c = self.c
        d = self.d
        f_opt = self.f
        g = self.g
        budget = self.A
        tri_inf = self.tri_inf
        pair_inf = self.pair_inf
        length = self.length
        
        num_b = len(B)
        num_s = len(S)
        num_l = len(L)
        upper_model = gb.Model("upper_model_bigM")
        eta = upper_model.addVar(lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name = "eta")        
        
        z = upper_model.addVars(num_b, vtype=GRB.BINARY, name= "_z")
        
        #BIG M constraint------------------------------------------------------
        #Alternative 1:
        '''m = upper_model.addVar(lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name = "m")
        q = upper_model.addVar(lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name = "q")
        
        upper_model.addConstrs((B[b][0]**2 - delta**2) * m*m + q*q - 2*B[b][0]*B[b][1]*m 
                               - 2*B[b][1]*q + 2*B[b][0] * m*q <= delta**2 - B[b][1]**2 + (BigM[b]**2) * (1-z[b]) for b in range(num_b) )'''
    
        #Alternative 2:
        head = upper_model.addVars(2, lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name = "head")
        tail = upper_model.addVars(2, lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name = "tail")
        
        t_b = upper_model.addVars(num_b, lb=0, ub=1, vtype=GRB.CONTINUOUS, name = "t")
        
        v_b = upper_model.addVars(num_b, 2, lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name = "v")
        
        for b in range(num_b):
            upper_model.addConstr((1-t_b[b])*head[0] + t_b[b]*tail[0] - B[b][0] == v_b[b,0])
            upper_model.addConstr((1-t_b[b])*head[1] + t_b[b]*tail[1] - B[b][1] == v_b[b,1])
            upper_model.addConstr(v_b[b,0]*v_b[b,0] + v_b[b,1]*v_b[b,1] <= delta**2 + (2*delta*BigM[b] + BigM[b]**2)*(1-z[b]))
            
        upper_model.addConstr(head[0]*head[0] - 2*head[0]*tail[0] + tail[0]*tail[0] + 
                              head[1]*head[1] - 2*head[1]*tail[1] + tail[1]*tail[1] <= length**2)
        upper_model.params.NonConvex = 2
        
        #Valid cuts------------------------------------------------------------
        for trip in tri_inf:
            upper_model.addConstr(z[trip[0]] + z[trip[1]] + z[trip[2]] <= 2)
        
        #for Alternative 2:
        for i in pair_inf:
            for j in pair_inf[i]:
                upper_model.addConstr(z[i] + z[j] <= 1)
                
        #Objective
        upper_model.setObjective(eta, GRB.MAXIMIZE)
        
        initial_r = subsubproblem(num_b, num_l, num_s, c, d, f_opt, np.zeros(num_b), g, budget)
        r_opt = initial_r[0]
        lw = initial_r[1]
        up = np.sum(g)
        z_opt = np.zeros(num_b)
        iteration = 0
        
        start_t = time.time()
        while up - lw > 0.01 * lw and (time.time() - start_t < 1800):
            upper_model.addConstr(eta <= gb.quicksum(gb.quicksum(gb.quicksum(z[b] * g[b][s][l] * f_opt[b][s] * r_opt[b][l]
                                                                             for l in range(num_l)) for s in range(num_s)) for b in range(num_b)))
            
            upper_model.params.MIPGap = 0.0
            upper_model.params.TimeLimit = 1800
            
            upper_model.optimize()
            #upper_model.write("upp.lp")
            #upper_model.write("upp.sol")
            
            up = upper_model.objVal
            for b in range(num_b):
                 zz = upper_model.getVarByName("_z[{bb}]".format(bb = b))
                 z_opt[b] = zz.x
            
            head = []
            tail = []
            head.append(upper_model.getVarByName("head[0]").x)
            head.append(upper_model.getVarByName("head[1]").x)
            tail.append(upper_model.getVarByName("tail[0]").x)
            tail.append(upper_model.getVarByName("tail[1]").x)
            self.head = head
            self.tail = tail        
            
           #Subsubproblem solving----------------------------------------------
            subsub_opt = subsubproblem(num_b, num_l, num_s, c, d, f_opt, z_opt, g, budget)
            r_opt = subsub_opt[0]
            lw = subsub_opt[1]
            
            iteration += 1
            
        return z_opt, up, iteration
    
    #__________________________________________________________________________    
    def upper_problem_composition(self):
        B = self.B
        delta = self.delta
        S = self.S
        L = self.L
        c = self.c
        d = self.d
        f_opt = self.f
        g = self.g
        budget = self.A
        tri_inf = self.tri_inf
        pair_inf = self.pair_inf
        quad_inf = self.quad_inf
        length = self.length
        
        num_b = len(B)
        num_s = len(S)
        num_l = len(L)
        upper_model = gb.Model("upper_model_composition")
        eta = upper_model.addVar(lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name = "eta")        
        
        z = upper_model.addVars(num_b, vtype=GRB.BINARY, name= "_z")
                
        #Valid cuts (Alternative 3)------------------------------------------------------------
        '''for quad in quad_inf:
            upper_model.addConstr(z[quad[0]] + z[quad[1]] + z[quad[2]] + z[quad[3]] <= 3)   
        
        for trip in tri_inf:
            upper_model.addConstr(z[trip[0]] + z[trip[1]] + z[trip[2]] <= 2)'''
        
        #setOfpair_inf = set()
        for i in pair_inf:
            for j in pair_inf[i]:
                upper_model.addConstr(z[i] + z[j] <= 1)
                #setOfpair_inf.add((i,j))

        #Infeasible collections from previous iterations (search user cuts)
        for collection in self.inf_collection:
            upper_model.addConstr(gb.quicksum(z[buil] for buil in collection) <= len(collection) - 1)
        
        #Objective
        upper_model.setObjective(eta, GRB.MAXIMIZE)
        
        #Initialize
        initial_r = subsubproblem(num_b, num_l, num_s, c, d, f_opt, np.zeros(num_b), g, budget)
        r_opt = initial_r[0]
        lw = initial_r[1]
        up = np.sum(g)
        z_opt = np.zeros(num_b)
        iteration = 0
        
        #Initial cut
        upper_model.addConstr(eta <= gb.quicksum(gb.quicksum(gb.quicksum(z[b] * g[b][s][l] * r_opt[b][s][l]
                                                                             for l in range(num_l)) for s in range(num_s)) for b in range(num_b)))
        
        #Heuristic approach - lower bound
        #upper_model.addConstr(eta >= self.max_inter_circ)
        
        upper_model.params.MIPGap = 0.0
        upper_model.params.TimeLimit = 3600
        #upper_model.Params.Threads = 16
         
        #1.Lazy cut approach
        upper_model._zVar = z
        upper_model._eta = eta
        upper_model.params.lazyConstraints = 1
        
        def mycallback_composition(upper_model, where):
            if where == GRB.Callback.MIPSOL:
                start_feasCheck = time.time()
                eta = upper_model._eta
                z = upper_model._zVar
                z_opt = upper_model.cbGetSolution(z)

                active_index = []
                for b in range(len(B)):
                    if z_opt[b] > 0.5:
                        active_index.append(b)

                #Do we need to check the feasibility if already known?
                model_status = 0
                do_check = True
                for records in self.fea_collection:
                    if records[0] == active_index:  #active index is already known as feasible sol
                        do_check = False
                        self.head = records[1]
                        self.tail = records[2]
                        model_status = 2   

                
                if do_check == True:
                    num_active = len(active_index)
                    #Stabbing line algorithm:
                    max_inter_act = 0
                    max_point = (0, 0, 0)
                    src_b = 0                    
                    
                    for i in active_index:
                        maximum_active = stab_line.stabbing_line(B[i], B[active_index], np.ones(num_active), delta, length+2*delta).max_overlap()
                        inter_circ = maximum_active[1] + 1
                        if inter_circ > max_inter_act:
                            max_inter_act = inter_circ
                            max_point = maximum_active[0]
                            src_b = i
                    
                    
                    if max_inter_act < num_active:
                        model_status = 3
                        upper_model.cbLazy(gb.quicksum(z[bb] for bb in active_index) <= max_inter_act)
                    elif max_inter_act == num_active:
                        segment_finderX = segment_finder(max_point, delta, B, src_b, length, active_index)
                        intersection1 = segment_finderX[0]
                        intersection2 = segment_finderX[1]

                        if Euclidean_distance(intersection1, intersection2) <= length:                            
                            model_status = 2                             
                            self.head = intersection1
                            self.tail = intersection2
                        else:
                            #Alt 1:
                            start_nonConv_feasCheck = time.time()
                            feasibility_check_f = feas_check.feasibility_check(B, delta, length, active_index)
                            #Alt 2:
                            #feasibility_check_f = feas_check.line_segement_feasibility(B, delta, length, active_index, side1_index, side2_index)
                            
                            self.nonConv_feasCheck_time += time.time() - start_nonConv_feasCheck
                            
                            model_status = feasibility_check_f[0]
                            if model_status == 2:
                                self.head = feasibility_check_f[1]
                                self.tail = feasibility_check_f[2]
                    else: 
                        print("Something has gone wrong!")
                        input()
                        
                ####################model status inspection  
                if model_status == 2:
                    #record feasible nodes for next iterations
                    self.new_fea_collection.append([active_index, self.head, self.tail])
                    
                    subsub_opt = subsubproblem(num_b, num_l, num_s, c, d, f_opt, z_opt, g, budget)
                    r_opt = subsub_opt[0]
                    lb = subsub_opt[1]
                    up = upper_model.cbGet(GRB.Callback.MIPSOL_OBJ)
                    if up - lb > 0.01 * lw:
                        upper_model.cbLazy(eta <= gb.quicksum(gb.quicksum(gb.quicksum(z[b] * g[b][s][l] * r_opt[b][s][l]
                                                                                 for l in range(num_l)) for s in range(num_s)) for b in range(num_b)))
                elif model_status == 3:
                    upper_model.cbLazy(gb.quicksum(z[bb] for bb in active_index) <= len(active_index) - 1)
                    self.new_inf_collection.append(active_index)
                    
                    comb3 = set(combinations(active_index, 3))
                    setOfTri_inf = set([tuple(ii) for ii in tri_inf])
                    listof3INF = list(comb3.intersection(setOfTri_inf))
                    for case in listof3INF:
                        upper_model.cbLazy(gb.quicksum(z[bb] for bb in case) <= 2)
                        
                    '''comb2 = set(combinations(active_index, 2))
                    listof2INF = list(comb2.intersection(setOfpair_inf))
                    for case in listof2INF:
                        upper_model.cbLazy(gb.quicksum(z[bb] for bb in case) <= 1)
                    
                    comb4 = set(combinations(active_index, 4))
                    setOfquad_inf = set([tuple(ii) for ii in quad_inf])
                    listof4INF = list(comb4.intersection(setOfquad_inf))
                    for case in listof4INF:
                        upper_model.cbLazy(gb.quicksum(z[bb] for bb in case) <= 3)'''
                    '''if len(active_index) < len(B)/2:
                        ext_p1 = active_index[0]
                        ext_p2 = active_index[len(active_index)-1]
                        max_distanceComb2 = 0
                        comb2 = list(combinations(active_index, 2))
                        for comb_iterator in comb2:
                            if max_distanceComb2 < Euclidean_distance(B[comb_iterator[0]], B[comb_iterator[1]]):
                                max_distanceComb2 = Euclidean_distance(B[comb_iterator[0]], B[comb_iterator[1]])                            
                                ext_p1 = comb_iterator[0]
                                ext_p2 = comb_iterator[1]
                        
                        slopeCenters = (B[ext_p2][1]-B[ext_p1][1]) / (B[ext_p2][0]-B[ext_p1][0])
                        yInterceptCenters = B[ext_p1][1] - slopeCenters*B[ext_p1][0]
                        
                        
                        distanceToCenter = [(distance_point_line(slopeCenters, yInterceptCenters, B[i])) for i in active_index]
                        distanceToCenter1 = sorted(enumerate(distanceToCenter), key=lambda i: i[1])
                        
                        checklist = [ext_p1, ext_p2]
                        for i in range(2, len(distanceToCenter1)):
                            checklist.append(active_index[distanceToCenter1[i][0]])
                            if feas_check.feasibility_check(B, delta, length, checklist)[0] == 3:
                                upper_model.cbLazy(gb.quicksum(z[bb] for bb in checklist) <= len(checklist)-1)
                                self.new_inf_collection.append(checklist)
                                
                                break'''
                        

                    
                else:
                    upper_model.cbLazy(gb.quicksum(z[bb] for bb in range(len(B))) <= 0)
                    self.new_inf_collection.append(active_index)
                    self.timelim_reached = model_status
                    #input()
                
                self.feasCheck_time += time.time() - start_feasCheck
                self.feasCheck_iteration += 1
                
        #SOLVE_________________________________________________________________
        upper_model.optimize(mycallback_composition)
       
        up = upper_model.objVal
        upper_model.write("upp.lp")
        upper_model.write("upp.sol")

        for b in range(num_b):
             zz = upper_model.getVarByName("_z[{bb}]".format(bb = b))
             z_opt[b] = zz.x
        
       
        '''
        #2.Regular cutting generation appraoch:
        start_t = time.time()
        while up - lw > 0.01 * lw and (time.time() - start_t < 1800):
            upper_model.addConstr(eta <= gb.quicksum(gb.quicksum(gb.quicksum(z[b] * g[b][s][l] * f_opt[b][s] * r_opt[b][l]
                                                                             for l in range(num_l)) for s in range(num_s)) for b in range(num_b)))
            
            
            model_status = 3
            while model_status == 3 and (time.time() - start_t < 1800):
                upper_model.optimize()
                #upper_model.write("upp.lp")
                #upper_model.write("upp.sol")

                up = upper_model.objVal
                active_index = []
                for b in range(num_b):
                     zz = upper_model.getVarByName("_z[{bb}]".format(bb = b))
                     z_opt[b] = zz.x                
                     if z_opt[b] > 0.5:
                        active_index.append(b)
                        
                #Do we need to check the feasibility if already known?
                do_check = True
                for records in self.fea_collection:
                    if records[0] == active_index:  #active index is already known as feasible sol
                        do_check = False
                        self.head = records[1]
                        self.tail = records[2]
                        model_status = 2   

                if do_check == True:
                    feasibility_check_f = feasibility_check(B, delta, length, active_index)
                    model_status = feasibility_check_f[0]
                    self.head = feasibility_check_f[1]
                    self.tail = feasibility_check_f[2]
                
                if model_status == 2:
                    #record feasible nodes for next iterations
                    self.new_fea_collection.append([active_index, self.head, self.tail])
                else:
                    upper_model.addConstr(gb.quicksum(z[bb] for bb in active_index) <= len(active_index) - 1)
                    self.new_inf_collection.append(active_index)
           
            
           #Subsubproblem solving----------------------------------------------
            subsub_opt = subsubproblem(num_b, num_l, num_s, c, d, f_opt, z_opt, g, budget)
            r_opt = subsub_opt[0]
            lw = subsub_opt[1]
            
            iteration += 1
            '''

        return z_opt, up, iteration

#Quadratic solution Ax^2+Bx+C = 0
def quad_solution(A,B,C):
    DELTA = B**2 - 4*A*C
    if DELTA >= 0: return (-B+np.sqrt(DELTA))/(2*A), (-B-np.sqrt(DELTA))/(2*A)
    elif -0.01 < DELTA < 0: return (-B+np.sqrt(0))/(2*A), (-B-np.sqrt(0))/(2*A)
    else:
        print("DELTA:", DELTA)
        print("Quadratic equation has no solution!")
        return 1000000, -1000000

def distance_point_line(slope, y_intercept, point):
    return abs(-slope*point[0] + point[1] - y_intercept)/np.sqrt(1+slope**2)

#Finding part of stabbing line that covers all active buildings
def segment_finder(max_point, delta, B, src_b, length, active_index):
    slop_tangent = -1/np.tan(max_point[0])
    x_tangent = delta*np.cos(max_point[0]) + B[src_b][0]
    y_tangent = delta*np.sin(max_point[0]) + B[src_b][1]
    b_tangent = y_tangent - slop_tangent*x_tangent

    intersection1 = [0,0] 
    intersection2 = [0,0]
    side1_index = 0
    side2_index = 0
    max_length = 0
    for i in range(len(active_index)-1):
        for j in range(i+1, len(active_index)):
            if (distance_point_line(slop_tangent, b_tangent, B[active_index[i]]) < 0.99*delta) and (distance_point_line(slop_tangent, b_tangent, B[active_index[j]]) < 0.99*delta):
                #intersection of line y=mx+b and circle (x-x0)^2+(y-y0)^2=delta^2
                A1 = slop_tangent**2 + 1
                B1 = -2 * B[active_index[i]][0] + 2 * slop_tangent * (b_tangent - B[active_index[i]][1])
                C1 =  B[active_index[i]][0]**2 + (b_tangent - B[active_index[i]][1])**2 - delta**2
                SolveEquation1 = quad_solution(A1, B1, C1)
                
                intersections_set1 = [[0,0], [0,0]]
                intersections_set1[0][0] = SolveEquation1[0]
                intersections_set1[0][1] = slop_tangent * intersections_set1[0][0] + b_tangent
                intersections_set1[1][0] = SolveEquation1[1]
                intersections_set1[1][1] = slop_tangent * intersections_set1[1][0] + b_tangent
            
                A2 = slop_tangent**2 + 1
                B2 = -2 * B[active_index[j]][0] + 2 * slop_tangent * (b_tangent - B[active_index[j]][1])
                C2 =  B[active_index[j]][0]**2 + (b_tangent - B[active_index[j]][1])**2 - delta**2
                SolveEquation2 = quad_solution(A2, B2, C2)
                
                intersections_set2 = [[0,0], [0,0]]
                intersections_set2[0][0] = SolveEquation2[0]
                intersections_set2[0][1] = slop_tangent * intersections_set2[0][0] + b_tangent
                intersections_set2[1][0] = SolveEquation2[1]
                intersections_set2[1][1] = slop_tangent * intersections_set2[1][0] + b_tangent
        
                minimum_dist_intersection = np.inf
                intSec1 = [0,0]
                intSec2 = [0,0]
                for k in intersections_set1:
                    for l in intersections_set2:
                        if minimum_dist_intersection > Euclidean_distance(k, l):
                            minimum_dist_intersection = Euclidean_distance(k, l)
                            intSec1 = k
                            intSec2 = l
                
                if max_length < minimum_dist_intersection:
                    max_length = minimum_dist_intersection
                    side1_index = active_index[i]
                    side2_index = active_index[j]
                    intersection1 = intSec1
                    intersection2 = intSec2
                
            
    return intersection1, intersection2, side1_index, side2_index


#SUBSUB+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def subsubproblem(num_b, num_l, num_s, c, d, f_opt, z, g, budget):
    subsubmodel = gb.Model("subsub_model")
    
    r_bsl= subsubmodel.addVars(num_b, num_s, num_l, vtype=GRB.BINARY, name= "_r")
    
    expr_obj = gb.LinExpr()
    cost1 = gb.LinExpr()
    cost2 = 0
        
    for b in range(num_b):
        for s in range(num_s):  
            expr_r = gb.LinExpr()
            for l in range(num_l):
                expr_r.add(r_bsl[b,s,l])
                cost1.add(c[b][s][l]*r_bsl[b,s,l])
                
                expr_obj.add(z[b]*g[b][s][l]*r_bsl[b,s,l])
                 
            cost2 += d[b][s]*f_opt[b][s]                   
            
            subsubmodel.addConstr(expr_r == f_opt[b][s])  
            
    subsubmodel.addConstr(cost1 <= budget - cost2)
    
    subsubmodel.setObjective(expr_obj, GRB.MINIMIZE)
    
    subsubmodel.params.MIPGap = 0.0
    subsubmodel.params.TimeLimit = 900
    subsubmodel.optimize()
    #subsubmodel.write("subsub.lp")
    #subsubmodel.write("subsub.sol")
    optimal = subsubmodel.objVal
    r_opt = np.zeros((num_b, num_s, num_l))
    for b in range(num_b):
        for s in range(num_s):
            for l in range(num_l):
                 rr = subsubmodel.getVarByName("_r[{bb},{ss},{ll}]".format(bb=b, ss=s ,ll=l))
                 r_opt[b][s][l] = rr.x

    return r_opt, optimal