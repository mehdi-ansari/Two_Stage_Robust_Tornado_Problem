# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 21:21:49 2021

@author: mehdi

see link: https://sahandsaba.com/line-intersecting-maximal-number-of-circles.html
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time

epsilon = 10**(-10)

def Euclidean_distance(a,b):
    return np.sqrt((a[0]-b[0])**2+(a[1]-b[1])**2)

def angle(p):
    if p[0] >= 0 and p[1] >= 0:
        return np.arctan((p[1]+epsilon)/(p[0]+epsilon))
    if p[0] < 0:
        return np.pi + np.arctan((p[1]+epsilon)/(p[0]+epsilon))
    if p[0] >= 0 and p[1] < 0:
        return 2*np.pi + np.arctan((p[1]+epsilon)/(p[0]+epsilon))

def external_tangent(p, q, radius):
    if p[1] != q[1]:
        mu = p[0] - q[0]
        eta = p[1] - q[1]
        
        x1 = np.sqrt((mu**2 + eta**2) * (radius**2 * eta**2)) / (mu**2 + eta**2)
        y1 = -mu*x1/eta

        x2 = -np.sqrt((mu**2 + eta**2) * (radius**2 * eta**2)) / (mu**2 + eta**2)
        y2 = -mu*x2/eta
        
        return (x1, y1), (x2, y2)
        
    else:
        return (0, radius), (0, -radius)
        
  
def internal_tangent(p, q, radius):
    mu = p[0] - q[0]
    eta = p[1] - q[1]
    delta = 2 * (radius**2)
        
    if p[1] != q[1]:        
        x1 = (-delta*mu + np.sqrt(delta**2 * mu**2 - (mu**2 + eta**2)*(delta**2 - radius**2 * eta**2))) / (mu**2 + eta**2) 
        y1 = -(mu*x1+delta)/eta

        x2 = (-delta*mu - np.sqrt(delta**2 * mu**2 - (mu**2 + eta**2)*(delta**2 - radius**2 * eta**2))) / (mu**2 + eta**2)
        y2 = -(mu*x2+delta)/eta
        
        return (x1, y1), (x2, y2)
        
    else:
        x = -delta / mu
        y1 = np.sqrt(radius**2 - x**2)
        y2 = -np.sqrt(radius**2 - x**2)
        
        return (x, y1), (x, y2)
    
    
def assign_coefficient(collected_set, c_bldngs):
    Initialization = 0
    #Assign second elements of set regarding the order of external/internal starting angle 0
    Frst_observ = collected_set[0][2]
    Scnd_observ = collected_set[1][2]
    if (Frst_observ == 22 and Scnd_observ == 11) or (Frst_observ  == 11 and Scnd_observ == 22):
        collected_set[0][1] = c_bldngs
        collected_set[1][1] = -c_bldngs
        collected_set[2][1] = c_bldngs
        collected_set[3][1] = -c_bldngs
    else:
        collected_set[0][1] = -c_bldngs
        collected_set[1][1] = c_bldngs
        collected_set[2][1] = -c_bldngs
        collected_set[3][1] = c_bldngs
        
        Initialization += c_bldngs        

    return collected_set, Initialization


class StabbingLine:
    def __init__(self, centers, coefficient, radius, max_dist):
        self.centers = centers
        self.coefficient = coefficient
        self.radius = radius
        self.max_dist = max_dist
    
    def events(self, index):
        C = self.centers[index]
        bldngs = self.centers.values()
        r = self.radius
        coef = self.coefficient
        Q = self.max_dist
        
        event_set_type = [("angle", float), ("counter", float), ("event_type", int)]
        Set_A = np.array([], dtype = event_set_type)
        
        #Event generator
        number_b = 0
        #Initialize
        init = 0
        for B in bldngs:
            collect_set = np.array([(-1, 0, 11), (-1, 0, 11),     #11: external, 22: internal
                                    (-1, 0, 22), (-1, 0, 22)], dtype = event_set_type)
            
            if epsilon < Euclidean_distance(B, C) <= Q:    #related to Length of tornado!! it means the length of tornado is Q-2*Delta
                p1 = external_tangent(C, B, r)[0]
                collect_set[0][0] = angle(p1)
                p2 = external_tangent(C, B, r)[1]
                collect_set[1][0] = angle(p2)
                if Euclidean_distance(B, C) > 2*r:
                    q1 = internal_tangent(C, B, r)[0]
                    collect_set[2][0] = angle(q1)                     
                    q2 = internal_tangent(C, B, r)[1]
                    collect_set[3][0] = angle(q2)
                
                collect_set = np.sort(collect_set, order = "angle")
                
                #-------------------------------------------------------------
                if Euclidean_distance(B, C) <= 2*r:     #The special case of intersection
                    collect_set = collect_set[2:]       #It removes -1 angles
                    if B[0] > C[0]:
                        collect_set[0][1] = -coef[number_b]
                        collect_set[1][1] = coef[number_b]
                        
                        init += coef[number_b]
                    else:
                        collect_set[0][1] = coef[number_b]
                        collect_set[1][1] = -coef[number_b]
                        
                    Set_A = np.r_[Set_A, collect_set]
                
                else:   #Euclidean_distance(B, C) > 2*r
                    assign_call = assign_coefficient(collect_set, coef[number_b])
                    Set_A = np.r_[Set_A, assign_call[0]]
                    init += assign_call[1]
           
            number_b += 1
        
        Set_A = np.sort(Set_A, order = "angle")
        
        return Set_A, init
             
    def max_overlap(self, C):
        events_call = self.events(C)
        A = events_call[0]
        track_sum = events_call[1]
        track = np.zeros(len(A))
        for e in range(len(A)):
            track_sum += A[e][1]
            track[e] = track_sum
        
        if len(A) == 0:
            A = [(0, 0, 0)]
            track = [0]
        return {'Angle_on_circle': round(A[np.argmax(track)][0],3), 'max_stabbed_weighted_circles': max(track)}

    def find_line_intersecting_maximal_circles(self):
        maxx = 0
        max_overlap_point = {}
        for i in self.centers.keys():
            overlaps = self.max_overlap(i)

            if overlaps['max_stabbed_weighted_circles'] > maxx:
                maxx = overlaps['max_stabbed_weighted_circles']
                max_overlap_point = overlaps
                max_overlap_point['center'] = i
                
        return max_overlap_point
                
        
'''
instance_name = 'counterexample1.xlsx'
instance = pd.read_excel(instance_name)
coord_type = [('x_coord', float), ('y_coord', float)]
buildings = np.array([(instance['x_coord'][i], instance['y_coord'][i]) for i in range(len(instance))],  dtype=coord_type)
B = np.sort(buildings, order='x_coord')  

num_buildings=len(B)
delta = 1
max_length = 2  #Q:lower bound, Q+2*delta: upper bound

c_b = [1 for i in range(num_buildings)]

max_inter_circ = 0
max_point = (0, 0, 0)
src_b = 0
start_time = time.time()
for i in range(len(B)):
    maximizer = stabbing_line(B[i], B, c_b, delta, max_length)
    max_overlap_call = maximizer.max_overlap()
    #print(max_overlap_call)
    intersected_circ = max_overlap_call[1] + c_b[i]
    if intersected_circ > max_inter_circ:
        max_inter_circ = intersected_circ
        max_point= max_overlap_call[0]
        src_b = i

print(max_inter_circ, max_point, src_b)
total_time = time.time() - start_time
total_time = "{:.2f}".format(total_time)
print("Running Time: ", total_time, " sec.")

#Plot++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
plt.figure()
plt.suptitle("{name}, max: {max_intrs}, time: {runtime}".format(name=instance_name, max_intrs = max_inter_circ, runtime=total_time))
             
for o in range(len(B)):
    plt.scatter(B[o][0], B[o][1])  

slop = -1/np.tan(max_point[0])
x0 = delta*np.cos(max_point[0]) + B[src_b][0]
y0 = delta*np.sin(max_point[0]) + B[src_b][1]

x_plot = range(-1,2)
plt.plot(x_plot, slop*(x_plot-x0) + y0)
'''