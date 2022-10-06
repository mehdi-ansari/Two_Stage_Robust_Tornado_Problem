# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 16:18:40 2021

@author: mehdi
"""
import numpy as np
import matplotlib.pyplot as plt
import time
#import Stabbing_line as stab_line

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
    def __init__(self, Param):
        self.coordinates = Param.InputData.coordinates
        self.delta = Param.width
        self.Q_maxLength = Param.length
        
        self.infeasible_pair = []
        self.feasible_pair_dict = {i:[] for i in self.coordinates.keys()}
        self.Coupling_conflict_set()
        
        self.Minus_Gamma = {}
        self.infeasible_tri = []
        self.Tri_conflict_set()
        
        
       
    def Coupling_conflict_set(self):    #This set identifies the pair of locations that are in an infeasible distance of 2*delta+Q      
        for loc1, coord1 in self.coordinates.items():
            for loc2, coord2 in self.coordinates.items():
                if Euclidean_distance(coord1, coord2) > 2*self.delta+self.Q_maxLength:
                    self.infeasible_pair.append((loc1, loc2))
                else:
                    self.feasible_pair_dict[loc1].append(loc2)
                    self.feasible_pair_dict[loc2].append(loc1)


    def Tri_conflict_set(self):
        for loc1, coord1 in self.coordinates.items():
            for loc2, coord2 in self.coordinates.items():
                if self.coordinates[loc1][0] < self.coordinates[loc2][0]:
                    self.Minus_Gamma[(loc1,loc2)] = []
                    
                    Bound = boundary_l2(coord1, coord2, self.delta)
    
                    if 2*self.delta < Bound.dist:
                        set_intrsectFeas_index = list(set(self.feasible_pair_dict[loc1]).intersection(set(self.feasible_pair_dict[loc2])))
                        
                        for k in set_intrsectFeas_index:
                            if self.coordinates[loc2][0] < self.coordinates[k][0]: 
                                line1 = False 
                                line2 = False
                                line3 = False
                                line4 = False
                                line5 = False
                                line6 = False
                                
                                #The upper parallel line---------------------------------------
                                if self.coordinates[k][1] > Bound.line1(self.coordinates[k][0]):
                                    line1 = True
                                  
                                #The lower parallel line---------------------------------------
                                if self.coordinates[k][1] < Bound.line2(self.coordinates[k][0]):
                                    line2 = True
                                    
                                #Clockwise rotating--------------------------------------------
                                if Bound.theta + Bound.alpha < np.pi/2:
                                    if self.coordinates[k][1] > Bound.line3(self.coordinates[k][0]):
                                        line3 = True
                                    
                                    if self.coordinates[k][1] < Bound.line5(self.coordinates[k][0]):
                                        line5 = True
                                        
                                else:
                                    if self.coordinates[k][1] < Bound.line3(self.coordinates[k][0]):
                                        line3 = True
                                        
                                    if self.coordinates[k][1] > Bound.line5(self.coordinates[k][0]):
                                        line5 = True
                                        
                                #Counter-Clockwise rotating------------------------------------
                                if Bound.theta - Bound.alpha > -np.pi/2:
                                    if self.coordinates[k][1] < Bound.line4(self.coordinates[k][0]):
                                        line4 = True
                                    
                                    if self.coordinates[k][1] > Bound.line6(self.coordinates[k][0]):
                                        line6 = True
                                        
                                else:
                                    if self.coordinates[k][1] > Bound.line4(self.coordinates[k][0]):
                                        line4 = True
                                        
                                    if self.coordinates[k][1] < Bound.line6(self.coordinates[k][0]):
                                        line6 = True
                                     
                                #Infeasible triple
                                if (line1==True and line3==True and line6==True) or (line2==True and line4==True and line5==True):
                                    self.infeasible_tri.append((loc1,loc2,k))
                                    
                                 
                                #New 3/27/2022 --- quadrant_infeasible: 
                                else:
                                    self.Minus_Gamma[(loc1,loc2)].append(k)                 
    
    
    '''def quad_conflict_set(self, timelimit):
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
                              
        return infeasible_quad'''
    
        
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