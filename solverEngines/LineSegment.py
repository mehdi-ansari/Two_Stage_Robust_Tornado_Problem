# -*- coding: utf-8 -*-
"""
Created on Sat Sep 24 22:00:56 2022

@author: mehdi
"""
import numpy as np

def Euclidean_distance(a,b):
    return np.sqrt((a[0]-b[0])**2+(a[1]-b[1])**2)

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
class Segment:
    def __init__(self, line_info, Param, damaged_location_coordindates):
        self.circle_center = line_info['center']
        self.angle = line_info['Angle_on_circle']
        self.damaged_location_coordindates = damaged_location_coordindates
        self.delta = Param.width
        self.length = Param.length
        self.coordinate = Param.InputData.coordinates[self.circle_center]
        
        
        self.intersection1 = [0,0] 
        self.intersection2 = [0,0]
        self.side1_index = 0
        self.side2_index = 0
        
        self.segment_finder()
        
        
    def segment_finder(self):
        angle = self.angle        
        delta = self.delta
        coordinate = self.coordinate
        
        
        slop_tangent = -1/np.tan(angle)
        x_tangent = delta*np.cos(angle) + coordinate[0]
        y_tangent = delta*np.sin(angle) + coordinate[1]
        b_tangent = y_tangent - slop_tangent*x_tangent

    
        
        max_length = 0
        for i, coord_i in self.damaged_location_coordindates.items():
            for j, coord_j in self.damaged_location_coordindates.items():

                #if (distance_point_line(slop_tangent, b_tangent, coord_i) < delta) and (distance_point_line(slop_tangent, b_tangent, coord_j) < delta):
                    #intersection of line y=mx+b and circle (x-x0)^2+(y-y0)^2=delta^2
                    A1 = slop_tangent**2 + 1
                    B1 = -2 * coord_i[0] + 2 * slop_tangent * (b_tangent - coord_i[1])
                    C1 =  coord_i[0]**2 + (b_tangent - coord_i[1])**2 - delta**2
                    SolveEquation1 = quad_solution(A1, B1, C1)
                    
                    intersections_set1 = [[0,0], [0,0]]
                    intersections_set1[0][0] = SolveEquation1[0]
                    intersections_set1[0][1] = slop_tangent * intersections_set1[0][0] + b_tangent
                    intersections_set1[1][0] = SolveEquation1[1]
                    intersections_set1[1][1] = slop_tangent * intersections_set1[1][0] + b_tangent
                
                    A2 = slop_tangent**2 + 1
                    B2 = -2 * coord_j[0] + 2 * slop_tangent * (b_tangent - coord_j[1])
                    C2 =  coord_j[0]**2 + (b_tangent - coord_j[1])**2 - delta**2
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
                        self.side1_index = i
                        self.side2_index = j
                        self.intersection1 = intSec1
                        self.intersection2 = intSec2
    
    def isSegment(self):
        if Euclidean_distance(self.intersection1, self.intersection2) <= self.length: 
            return True
        else:
            return False
        