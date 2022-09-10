# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 15:30:39 2022

@author: mehdi
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import gurobipy as gb
from gurobipy import GRB
import problems_class as pb
import time
import ValidCuts_l2_class as ValidCuts
import Stabbing_line as stab_line
import random
import ListExcel_reader as Joplin_data

from multiprocessing import Process



timelimitForQuad = 0

#PROBLEM_________________________________________________________________________
#master_problem
master_model = gb.Model("master_model")
master_model_class = pb.master_problem(master_model, B, S, L, w, g, c, d, budget)   #class of master problem generating constraints and variables

#fixed parts
fixed_part = master_model_class.fixed_part()    #generate fixed constraint of master problem
f_bs = fixed_part[0]    #function returns vars. f
theta = fixed_part[1]   #function returns var. theta


#Valid cuts-------------------------------------------------------------------
cuts_class = ValidCuts.valid_cuts_sets(B, delta, num_b, length)
pair_inf = cuts_class.Coupling_conflict_set()[0]
tri_inf = cuts_class.Tri_conflict_set()
quad_inf = cuts_class.quad_conflict_set(timelimitForQuad)

#Heuristic - lower bound of tornado---------------'------
'''g_b00 = []

for b in range(num_b):
    g_b00.append(g[b][0][0])
max_inter_circ = 0
for i in range(num_b):
    maximum_stabbed = stab_line.stabbing_line(B[i], B, g_b00, delta, length/2).max_overlap()
    intersected_circ = maximum_stabbed[1] + g_b00[i]
    if intersected_circ > max_inter_circ:
        max_inter_circ = intersected_circ
        

print("Heuristic done!")'''
#------------------------------

#main loop++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
begin_time = time.time()
#array of new adding varibales:
ri_bsl = []
iteration = 0
z = [np.zeros(num_b)] #for test
path_record = [([0,0],[0,0])]
overal_subproblem_time = 0
overal_subproblem_iteration = 0
inf_collection = []     #set of infeasible collections that are feasible by valid cuts but infeasible in original formulation
fea_collection = []
total_nonConv_feasCheck_time = 0
total_feasCheck_time = 0
total_feasCheck_iter = 0

slope = 0
intercept = 0

lb = 0
up = np.sum(g)
while up - lb > 0.01 * lb and (time.time()-begin_time<MainRuntimelimit):
    ri_bsl.append(0)
    
    master_model_class.generating_part(f_bs, theta, ri_bsl, z, iteration)    #variable and constraint generation
    
    #Solving master problem
    print("Master Problem: ")
    
    master_model.params.MIPGap = 0.0
    master_model.params.TimeLimit = 3600
    #master_model.Params.Threads = 4
    
    master_model.optimize()
    master_model.write("master_model.sol")
    master_model.write("master_model.lp")
    lb=master_model.objVal

    #Get values
    f_opt = np.zeros((num_b, num_s))
    first_term_sum = 0
    for b in range(num_b):
        for s in range(num_s):
            f_opt1 = master_model.getVarByName("_f[{bb},{ss}]".format(bb=b, ss=s))
            f_opt[b][s] = f_opt1.x
            first_term_sum += f_opt[b][s] * w[b][s]

    #Subproblem++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
    sub_begin = time.time()
    
    subproblem = pb.subproblem(B, delta, S, L, c, d, f_opt, g, budget, quad_inf,
                               tri_inf, pair_inf, length, inf_collection, fea_collection, 0)
    bilevel = subproblem.upper_problem_composition()
    #bilevel = subproblem.upper_problem_bigM()
    z.append(bilevel[0])
    up = bilevel[1] + first_term_sum

    sub_end = time.time()
    overal_subproblem_time += sub_end - sub_begin
    overal_subproblem_iteration += bilevel[2]
    
    head = subproblem.head
    tail = subproblem.tail
    path_record.append((head,tail))
    
    inf_collection = inf_collection + subproblem.new_inf_collection     #update the set of infeasible collection for the next iteration
    fea_collection = fea_collection + subproblem.new_fea_collection
    
    total_nonConv_feasCheck_time += subproblem.nonConv_feasCheck_time
    total_feasCheck_time += subproblem.feasCheck_time
    total_feasCheck_iter += subproblem.feasCheck_iteration
    
    iteration += 1

end_time = time.time()

print("Best objective: {:.2f}".format(lb))
print("Best Bound: {:.2f}".format(up))
print("main iteration: ", iteration)
print("running time: {:.2f} sec.".format(end_time-begin_time))
print("total subproblem run time: {:.2f} sec.".format(overal_subproblem_time))
print("total feasibility check run time: ", total_nonConv_feasCheck_time)
#print("overall subproblem iteration: ", overal_subproblem_iteration)

print("Number of infeasible collections:", len(inf_collection))

print("Time limit reached for feasibility check?", subproblem.timelim_reached)

print("Number of infeasible tiple:", len(tri_inf))
print("Number of infeasible quads:", len(quad_inf))


#PLOT++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#First identify the worst tornado:
summm_max = 0
r_max = 0
z_max = 0
path_max = 0
for iterat in range(iteration):
    summm = 0
    r = np.zeros((num_b,num_s,num_l))
    
    for b in range(num_b):
        for s in S:
            for l in L:
                r[b][s][l] = master_model.getVarByName("_r{iit}[{bb},{ss},{ll}]".format(iit=iterat,
                                                                             bb=b,ss=s,ll=l)).x
                if z[iterat][b] >= 0.5:
                    if r[b][s][l] >= 0.5:
                        summm+=g[b][s][l]
    
    if summm_max <= summm:
        summm_max = summm
        r_max = r
        z_max = z[iterat]
        path_max = path_record[iterat]



'''plt.suptitle(instance_name)
for o in range(len(B)):
    plt.scatter(B[o][0], B[o][1], c=tract[b][2])

plt.plot([head[0], tail[0]], [head[1], tail[1]])'''

#plt.xlim(-94.61806, -94.05767)
#plt.ylim(36.7495, 37.35671)
#plt.suptitle("EF-5 - {:.2f} sec. for {} buildings".format(end_time-begin_time, num_b))
act_indx = []
for i in range(num_b):
    colour = "black"
    
    for s in range(1,num_s):
        if f_opt[i][s] > 0.5:
            colour = 'blue'
            
    if z_max[i] >= 0.5: 
        act_indx.append(i)
        colour = "red"
        for s in range(1,num_s):
            if f_opt[i][s] > 0.5:
                colour = 'darkviolet'
        
    plt.scatter(B[i][0], B[i][1], c = colour)
#plt.colorbar()
point1 = 0
point2 = 0
for collect in fea_collection:
    if collect[0] == act_indx:
        point1 = collect[1]
        point2 = collect[2]
plt.plot([point1[0], point2[0]], [point1[1], point2[1]], c = "red")

plt.xticks([])
plt.yticks([])

#Outputs++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
outfile = open('{}Clusters_Budget{}M_{}Miles.csv'.format(num_b,round(budget/1000000),length),'w')
outfile.write("Budget :" + "," + str(budget) + ',' + "Tornado Length:" + "," + str(length) +'\n')
outfile.write('Best objective (Dislocation):' + ', ' + str(round(lb))+'\n')
outfile.write('Best bound:' + ', ' + str(round(up))+'\n')
outfile.write('Gap:' + ',' + '{:.2f}'.format(abs(up-lb)/up) + '\n')
outfile.write("main iteration: " + ',' + str(iteration)+'\n')
outfile.write("Run time :" + "," + str(round(end_time-begin_time))+'\n')
outfile.write("Subproblem run time:" + "," + str(round(overal_subproblem_time)) + '\n')
outfile.write("Number of feasibility check (CallBack): " + ',' + str(total_feasCheck_iter) + '\n')
outfile.write("Feasibility check run time: " + ',' + str(round(total_feasCheck_time)) + '\n')
outfile.write("Nonconvex to check run time: " + ',' + str(round(total_nonConv_feasCheck_time)) + '\n')
outfile.write("Number of infeasible pairs" + ',' + str(sum([len(pair_inf[i]) for i in range(len(pair_inf))])) + '\n')
outfile.write("Number of infeasible tiples:" + ','+ str(len(tri_inf)) + '\n')
outfile.write("Number of infeasible quadruples:" +',' +  str(len(quad_inf)) + '\n')


outfile.write('\n'+'Block ID' +','+ 
              'Population' + ','+
              'Retrofitting' + ',' + 
              'Retrofitting Cost' + ','
              + 'Damaged?' + ','+
              'Recovery Cost' + ','+
              'Dislocation' +'\n')


#Then, report:
for index,rows in instance.iterrows():
    damaged = "No"
    if z_max[index] >= 0.5:
        damaged = "Yes"
        
    Ret_strategy = "Do Nothing"
    Ret_cost = "0"
    Rec_cost = "0"
    Dislocation = 0.0
    if damaged == "Yes": Dislocation = g[index][0][0]
    if r_max[index][0][1] >= 0.5: 
        Rec_cost = str(c[index][0][1])
        Dislocation = g[index][0][1]
        
    for i in S[1:]:
        if f_opt[index][i] >= 0.5: 
            Ret_strategy = "R{}".format(i)
            Ret_cost = str(rows["Cost "+Ret_strategy])

            if damaged == "Yes": Dislocation = g[index][i][0]

            if r_max[index][i][1] >= 0.5:
                Rec_cost = str(c[index][i][1])
                Dislocation = g[index][i][1]
        
    
    outfile.write(str(rows['Block id']) +','+ 
              str(rows['Population']) + ',' +
              Ret_strategy + ',' + 
              Ret_cost + ',' +
              damaged + ',' +
              Rec_cost + ',' +
              str(Dislocation) +'\n')


outfile.close()