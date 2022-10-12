# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 23:05:10 2022

@author: mehdi
"""
import matplotlib.pyplot as plt
from solverEngines.ValidCuts import boundary_l2
import numpy as np
#from Results import Results

class Figures:
    def __init__(self, Results):
        fig, ax = plt.subplots()
        
        scaled_width = Results.Param.width/54.6
        
        xH = Results.head_worst[0]/54.6
        yH = Results.head_worst[1]/69
        xT = Results.tail_worst[0]/54.6
        yT = Results.tail_worst[1]/69 
        
        
        if Results.head_worst[0] > Results.tail_worst[0]:
            xT = Results.head_worst[0]/54.6
            yT = Results.head_worst[1]/69
            xH = Results.tail_worst[0]/54.6
            yH = Results.tail_worst[1]/69        
        
        lines = boundary_l2((xH, yH), (xT,yT), scaled_width/2)
        
        ax.plot([xH, xT], [yH, yT], c = 'red')
        
        xH_line1 = xH - scaled_width * np.sin(np.arctan((yT-yH)/(xT-xH)))
        xT_line1 = xT - scaled_width * np.sin(np.arctan((yT-yH)/(xT-xH)))
        ax.plot([xH_line1, xT_line1], [lines.line1(xH_line1), lines.line1(xT_line1)], linestyle='dashed', color='red', linewidth=1)
        
        xH_line2 = xH + scaled_width * np.sin(np.arctan((yT-yH)/(xT-xH)))
        xT_line2 = xT + scaled_width * np.sin(np.arctan((yT-yH)/(xT-xH)))
        ax.plot([xH_line2, xT_line2], [lines.line2(xH_line2), lines.line2(xT_line2)], linestyle='dashed', color='red', linewidth=1)
        
        circle1 = plt.Circle((xH, yH), scaled_width, linestyle='dashed', color='red', linewidth=1, fill=False)           
        ax.add_patch(circle1)
        circle2 = plt.Circle((xT, yT), scaled_width, linestyle='dashed', color='red', linewidth=1, fill=False)           
        ax.add_patch(circle2)
        
        #Plot Locations:
        for ID in Results.location_indx:
            x_coord = round(Results.Param.InputData.coordinates[ID][0]/54.6, 3)
            y_coord = round(Results.Param.InputData.coordinates[ID][1]/69, 3)
            symbol = 'o'
            color = 'black'
            
            for s in Results.retrofit_indx:
                if Results.f_sol[(ID,s)] > 0.5 and s != 'Do_nothing': 
                    symbol = 's'
                    color = 'blue'
                        
            if Results.z_worst[ID] > 0.5: 
                if symbol == 's':
                    symbol = 'D'
                    color = 'darkviolet'
                else: 
                    symbol = '*'
                    color = 'red'
                
                for s in Results.retrofit_indx:
                    for p in Results.recovery_indx:
                        if Results.r_worst[(ID,s,p)] > 0.5 and p != 'Do_nothing':
                            if symbol == 'D':
                                symbol = '^'
                                color = 'green'
                            else: 
                                symbol = 'h'
                                color = 'yellow'
            
            plt.scatter(x_coord, y_coord, marker = symbol, c = color)
            # plt.xlim([-94.58, -94.43])
            # plt.ylim([37, 37.15])