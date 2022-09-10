# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 16:18:40 2021

@author: mehdi
"""
import numpy as np
import matplotlib.pyplot as plt
import LnSeg_feasCheck_circles as feas_check
import time
import Stabbing_line as stab_line
import problems_class as pb

def Euclidean_distance(a,b):
    return np.sqrt((a[0]-b[0])**2+(a[1]-b[1])**2)

#Class of boundaries for tri constraints
class boundary_l2:
    def __init__(self, left_node, right_node, delta):
        self.left_node = left_node
        self.right_node = right_node
        self.delta = delta
        
        self.dist = Euclidean_distance(left_node, right_node) #Distance between a pair of buildings
        if 2*delta < self.dist:
            self.alpha = np.arcsin(2*delta/self.dist)     #Rotating threshhold
        
        epsilon = 0.000001
        self.theta = np.arctan((right_node[1]-left_node[1]+epsilon)/(right_node[0]-left_node[0]+epsilon))      #angle between a pair of buildings
        if abs(self.theta) <= epsilon: self.theta=0
        
    def line1(self, x_input): 
        return np.tan(self.theta) * (x_input - self.left_node[0]) + self.left_node[1] + (2*self.delta/np.cos(self.theta))

    def line2(self, x_input): 
        return np.tan(self.theta) * (x_input - self.left_node[0]) + self.left_node[1] - (2*self.delta/np.cos(self.theta))
    
    def line3(self, x_input):
        return np.tan(self.theta + self.alpha) * (x_input - self.left_node[0]) + self.left_node[1]

    def line4(self, x_input):
        return np.tan(self.theta - self.alpha) * (x_input - self.left_node[0]) + self.left_node[1]
    
    def line5(self, x_input):
        return np.tan(self.theta + self.alpha) * (x_input - self.right_node[0]) + self.right_node[1]

    def line6(self, x_input):
        return np.tan(self.theta - self.alpha) * (x_input - self.right_node[0]) + self.right_node[1]

#class of valid cuts for coupling and tri constraints:
class valid_cuts_sets:
    def __init__(self, Sorted_Buildings, delta, num_buildings, Q_maxLength):
        self.B = Sorted_Buildings
        self.delta = delta
        self.num_buildings = num_buildings
        self.Q_maxLength = Q_maxLength
        
        self.Minus_Gamma = {}
       
    def Coupling_conflict_set(self):    #This set identifies the pair of locations that are in an infeasible distance of 2*delta+Q
        infeasible_forward_pair = {}
        feasible_forward_pair = {}
        BB = self.B
        for i in range(self.num_buildings - 1):
            infeasible_forward_pair[i] = []
            feasible_forward_pair[i] = []
            for j in range(i+1, self.num_buildings):
                if Euclidean_distance(BB[i], BB[j]) > 2*self.delta+self.Q_maxLength:
                    infeasible_forward_pair[i].append(j)
                else:
                    feasible_forward_pair[i].append(j)
                    
        return infeasible_forward_pair, feasible_forward_pair

    def Tri_conflict_set(self):
        infeasible_forward_tri = []
        
        feas_pairs = self.Coupling_conflict_set()[1]
        B = self.B
        for i in range(len(B) - 2):
            for j in range(i + 1, len(B) - 1):
                self.Minus_Gamma[(i,j)] = []
                Bound = boundary_l2(B[i], B[j], self.delta)
                dist = Bound.dist 

                if 2*self.delta < dist:
                    set_intrsectFeas_index = list(set(feas_pairs[i]).intersection(set(feas_pairs[j])))
                    for k in set_intrsectFeas_index:
                        line1 = False 
                        line2 = False
                        line3 = False
                        line4 = False
                        line5 = False
                        line6 = False
                        
                        #The upper parallel line---------------------------------------
                        if B[k][1] > Bound.line1(B[k][0]):
                            line1 = True
                          
                        #The lower parallel line---------------------------------------
                        if B[k][1] < Bound.line2(B[k][0]):
                            line2 = True
                            
                        #Clockwise rotating--------------------------------------------
                        if Bound.theta + Bound.alpha < np.pi/2:
                            if B[k][1] > Bound.line3(B[k][0]):
                                line3 = True
                            
                            if B[k][1] < Bound.line5(B[k][0]):
                                line5 = True
                                
                        else:
                            if B[k][1] < Bound.line3(B[k][0]):
                                line3 = True
                                
                            if B[k][1] > Bound.line5(B[k][0]):
                                line5 = True
                                
                        #Counter-Clockwise rotating------------------------------------
                        if Bound.theta - Bound.alpha > -np.pi/2:
                            if B[k][1] < Bound.line4(B[k][0]):
                                line4 = True
                            
                            if B[k][1] > Bound.line6(B[k][0]):
                                line6 = True
                                
                        else:
                            if B[k][1] > Bound.line4(B[k][0]):
                                line4 = True
                                
                            if B[k][1] < Bound.line6(B[k][0]):
                                line6 = True
                             
                        #Infeasible triple
                        if (line1==True and line3==True and line6==True) or (line2==True and line4==True and line5==True):
                            infeasible_forward_tri.append([i,j,k])
                         
                        #New 3/27/2022 --- quadrant_infeasible: 
                        else:
                            self.Minus_Gamma[(i,j)].append(k)
                                    
                                    
        return infeasible_forward_tri                   
    
    
    def quad_conflict_set(self, timelimit):
        infeasible_quad = []
        Minus_Gamma = self.Minus_Gamma
        B = self.B
        
        start_time = time.time()
        count = 0
        for pairs in Minus_Gamma:
            for third_point in Minus_Gamma[pairs]:
                if third_point != len(B)-1 and (time.time() - start_time < timelimit):
                    set1 = set(Minus_Gamma[pairs])
                    set2 = set(Minus_Gamma[(pairs[0],third_point)])
                    set3 = set(Minus_Gamma[(pairs[1],third_point)])
                    intersec_set = set1.intersection(set2.intersection(set3))
                    
                    for forth_point in list(intersec_set):
                        count+=1
                        pointsToCheck = [pairs[0],pairs[1],third_point,forth_point]                        
                        
                        max_inter_act = 0
                        max_point = (0, 0, 0)
                        src_b = 0                    
                    
                        for i in pointsToCheck:
                            maximum_active = stab_line.stabbing_line(B[i], B[pointsToCheck], np.ones(4), self.delta, self.Q_maxLength+2*self.delta).max_overlap()
                            inter_circ = maximum_active[1] + 1
                            if inter_circ > max_inter_act:
                                max_inter_act = inter_circ
                                max_point = maximum_active[0]
                                src_b = i
                        
                        
                        if max_inter_act < 4:
                            infeasible_quad.append(pointsToCheck)
                        else:
                            segment_finderX = pb.segment_finder(max_point, self.delta, B, src_b, self.Q_maxLength, pointsToCheck)
                            intersection1 = segment_finderX[0]
                            intersection2 = segment_finderX[1]

                            if Euclidean_distance(intersection1, intersection2) > self.Q_maxLength: 
                                feasibilityCheck = feas_check.feasibility_check(B, self.delta, self.Q_maxLength, pointsToCheck)
                                if feasibilityCheck[0] == 3:
                                    infeasible_quad.append(pointsToCheck)
                              
        return infeasible_quad
    
        
'''
B = ((5.50873069, 9.67334179), (6.59114546, 0.75576362), (4.07456233, 5.14895151), (4.36333754, 2.99619805))
delta = 1

for i in range(len(B) - 2):
    for j in range(i + 1, len(B) - 1):
        Bound = boundary_l2(B[i], B[j], delta)
        plt.figure()
        x = np.array(range(int(B[0][0])-10,int(B[len(B)-1][0])+10))
        plt.xlim(-4, 18)
        plt.ylim(-4, 18)
        plt.suptitle(str(B[i])+str(B[j]))
                
        for o in range(len(B)):
            plt.scatter(B[o][0], B[o][1])          
            #plt.scatter(8.5, 12)   
        plt.plot(x, Bound.line1(x))
        plt.plot(x, Bound.line2(x))
        plt.plot(x, Bound.line3(x))
        plt.plot(x, Bound.line4(x))
        plt.plot(x, Bound.line5(x))
        plt.plot(x, Bound.line6(x))
        
setss = valid_cuts_sets(B, delta, len(B), 32)
print(setss.Tri_conflict_set())'''