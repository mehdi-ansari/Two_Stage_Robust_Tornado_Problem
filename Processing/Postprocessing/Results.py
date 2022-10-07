# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 22:16:15 2022

@author: mehdi
"""
from datetime import datetime
from solverEngines.UncertaintySet import UncertaintySet

class Results:
    def __init__(self, CCGAlgorithm, ROOT_DIR):
        self.ROOT_DIR = ROOT_DIR
        self.CCGAlgorithm = CCGAlgorithm
        self.Param = CCGAlgorithm.Param
        self.master_sol_dict = self.CCGAlgorithm.master_problem.master_sol_dict
        self.r_sols = {i: self.master_sol_dict['iteration{}'.format(self.CCGAlgorithm.iteration-1)]['r_sol[{}]'.format(i)] for i in range(self.CCGAlgorithm.iteration)}
        self.theta = self.master_sol_dict['iteration{}'.format(self.CCGAlgorithm.iteration-1)]['theta']
        self.f_sol = self.master_sol_dict['iteration{}'.format(self.CCGAlgorithm.iteration-1)]['f_sol']
        self.z_sol = self.CCGAlgorithm.z_sol
        
        self.lower_bound = round(CCGAlgorithm.lower_bound)
        self.upper_bound = round(CCGAlgorithm.upper_bound)
        self.gap = round((self.upper_bound - self.lower_bound)/(self.lower_bound+0.00001), 2)
        self.run_time = round(CCGAlgorithm.run_time, 2)
        self.iteration = CCGAlgorithm.iteration
        self.subproblem_run_time = round(CCGAlgorithm.subproblem_run_time, 2)
        self.number_subproblem_callbacks = CCGAlgorithm.number_subproblem_callbacks
        self.subproblem_callback_run_time = round(CCGAlgorithm.subproblem_callback_run_time, 2)
        self.number_uncertaintySet_check = CCGAlgorithm.number_uncertaintySet_check
        self.uncertaintySet_check_run_time = round(CCGAlgorithm.uncertaintySet_check_run_time, 2)
        
        self.find_worst_case_tornado()
        
    def print_results(self):
        print('Best Bound:', self.lower_bound)
        print('Best Objective:', self.upper_bound)
        print('Gap:', self.gap)
        print('CCG Run time:', self.run_time)
        print('CCG Iteration:', self.iteration)
        print('Subproblem Run time:', self.subproblem_run_time)
        print('Number of Subproblem Callbacks:', self.number_subproblem_callbacks)
        print('Subproblem Callbacks Run Time:', self.subproblem_callback_run_time)
        print('Number of Uncertainty Set Check Call:', self.number_uncertaintySet_check)
        print('Uncertainty Set Check Run Time:', self.uncertaintySet_check_run_time)
        
        
    def find_worst_case_tornado(self):
        self.location_indx = list(self.Param.InputData.first_stage_dislocation.keys())
        self.retrofit_indx = list(self.Param.InputData.retrofitting_strategies)
        self.recovery_indx = list(self.Param.InputData.recovery_strategies)
        self.z_worst={}
        self.head_worst = 0
        self.tail_worst = 0
        for i in self.r_sols.keys():
            if abs(self.theta - sum(self.Param.InputData.second_stage_dislocation[l][s][p]*self.z_sol[i][l]*self.r_sols[i][(l,s,p)]
                                                       for l in self.location_indx for s in self.retrofit_indx for p in self.recovery_indx)) < 0.01:

                self.z_worst = self.z_sol[i]
                self.r_worst = self.r_sols[i]
                
                damaged_location_coordindates = {}
                for l, soluation in self.z_worst.items():
                    if soluation >= 0.5:
                        damaged_location_coordindates[l] = self.Param.InputData.coordinates[l]
                check_for_line = UncertaintySet(self.Param).check_feasibility(damaged_location_coordindates)
                self.head_worst = check_for_line['head']
                self.tail_worst = check_for_line['tail']
                break
        
    def make_file(self):
        now = datetime.now()
        now = str(now.date())+'_'+now.strftime("%H-%M-%S")
        
        with open(str(self.ROOT_DIR)+'/Results/{}_{}_{}Clusters_Budget{}M_{}Miles.csv'.format(now, self.Param.user.input_dict['input_data_name'], len(self.Param.InputData.coordinates) ,round(self.Param.budget/1000000),self.Param.length),'w') as csv_file:
            csv_file.write("Budget :" + "," + str(self.Param.budget) + ',' + "Tornado Length:" + "," + str(self.Param.length) +'\n')
            csv_file.write('Best Bound:'+ "," + str(self.lower_bound)+ '\n')
            csv_file.write('Best Objective:'+ "," + str(self.upper_bound)+ '\n')
            csv_file.write('Gap:'+ "," + str(self.gap)+ '\n')
            csv_file.write('CCG Run time:'+ "," + str(self.run_time)+ '\n')
            csv_file.write('CCG Iteration:'+ "," + str(self.iteration)+ '\n')
            csv_file.write('Subproblem Run time:'+ "," + str(self.subproblem_run_time)+ '\n')
            csv_file.write('Number of Subproblem Callbacks:'+ "," + str(self.number_subproblem_callbacks)+ '\n')
            csv_file.write('Subproblem Callbacks Run Time:'+ "," + str(self.subproblem_callback_run_time)+ '\n')
            csv_file.write('Number of Uncertainty Set Check Call:'+ "," + str(self.number_uncertaintySet_check)+ '\n')
            csv_file.write('Uncertainty Set Check Run Time:'+ "," + str(self.uncertaintySet_check_run_time)+ '\n')
            csv_file.write('Head of worst tornado:'+ "," +str(round(self.head_worst[0]/54.6, 3)) + "," +str(round(self.head_worst[1]/69, 3))+'\n')
            csv_file.write('Tail of worst tornado:'+ "," +str(round(self.tail_worst[0]/54.6, 3)) + "," +str(round(self.tail_worst[1]/69, 3))+'\n')
            
            csv_file.write('\n'+'Block ID' +','+ 
                           'coordinate_x' + ','+
                           'coordinate_y' + ','+
                           'Population' + ','+
                           'Retrofitting' + ',' + 
                           'Retrofitting Cost' + ','
                            + 'Damaged?' + ','+
                            'Recovery Cost' + ','+
                            'Dislocation' +'\n')
            
            for ID in self.location_indx:
                csv_file.write(str(ID) + ',' + 
                               str(round(self.Param.InputData.coordinates[ID][0]/54.6, 3)) + ',' +
                               str(round(self.Param.InputData.coordinates[ID][1]/69, 3)) + ',' +
                               str(self.Param.InputData.population[ID]) + ',')
               
                for s in self.retrofit_indx:
                    if self.f_sol[(ID,s)] > 0.5: 
                        csv_file.write(str(s) + ',')
                        csv_file.write(str(self.Param.InputData.cost_retrofitting[ID][s]) + ',')
                        
                if self.z_worst[ID] > 0.5: 
                    csv_file.write('YES' + ',')
                    for s in self.retrofit_indx:
                        for p in self.recovery_indx:
                            if self.r_worst[(ID,s,p)] > 0.5:
                                csv_file.write(str(self.Param.InputData.cost_recovery[ID][s][p]) + ',')
                                csv_file.write(str(self.Param.InputData.second_stage_dislocation[ID][s][p]) + ',')
                else: 
                    csv_file.write('NO' + ',')
                    csv_file.write(str(0) + ',')
                    csv_file.write(str(0) + ',')
                
                
                            
                csv_file.write('\n')
                        
                    

                
             