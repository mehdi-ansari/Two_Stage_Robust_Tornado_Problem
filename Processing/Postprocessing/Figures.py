# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 23:05:10 2022

@author: mehdi
"""
import matplotlib.pyplot as plt
#from Results import Results

class Figures:
    def __init__(self, Results):
        plt.plot([Results.head_worst[0]/54.6, Results.tail_worst[0]/54.6], 
                 [Results.head_worst[1]/69, Results.tail_worst[1]/69], c = "red")
        
        for ID in Results.location_indx:
            x_coord = round(Results.Param.InputData.coordinates[ID][0]/54.6, 3)
            y_coord = round(Results.Param.InputData.coordinates[ID][1]/69, 3)
            color = 'black'
            
            for s in Results.retrofit_indx:
                if Results.f_sol[(ID,s)] > 0.5 and s != 'Do_nothing': 
                    color = 'blue'
                        
            if Results.z_worst[ID] > 0.5: 
                if color == 'blue':
                    color = 'darkviolet'
                else: color = 'red'
                
                for s in Results.retrofit_indx:
                    for p in Results.recovery_indx:
                        if Results.r_worst[(ID,s,p)] > 0.5 and p != 'Do_nothing':
                            if color == 'darkviolet':
                                color = 'green'
                            else: color = 'yellow'
            
            plt.scatter(x_coord, y_coord, c = color)
            # plt.xlim([-94.58, -94.43])
            # plt.ylim([37, 37.15])