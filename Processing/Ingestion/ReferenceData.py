# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 17:56:24 2022

@author: mehdi
"""
import numpy as np
from scipy.stats import lognorm

def lognorm_cdf(x, mu, sigma):
    shape  = sigma
    loc    = 0
    scale  = np.exp(mu)
    return lognorm.cdf(x, shape, loc, scale)

class ReferenceData:
    def __init__(self, **kwargs):
        '''
        This module input data from one of our reference papers:
            Koliou, Maria, and John W. van de Lindt. "Development of building restoration functions for use in community recovery planning to tornadoes." Natural Hazards Review 21.2 (2020): 04020004.

        Returns
        -------
        Input data.

        '''        
        self.cost_per_m2 = 862
        self.cost_per_ft2 = 0.093 * self.cost_per_m2 
        
        #The percentage of total cost to recover the whole area:
        self.costPercentage = [0, 0.005, 0.02, 0.12, 0.23]       #Table 9
        
        #The probability of not being in %100 functionality after x days
        elapsed_days = kwargs['elapsed_days']
        mean_sd_lognormal = ((3.09, 0.51), (3.52, 0.55), (4.62, 0.55), (5.19, 0.52))    #Table 5
        self.P_NotRecoverd = [0.00]    
        for (mu, sigma) in mean_sd_lognormal:
            self.P_NotRecoverd.append(round(1-lognorm_cdf(elapsed_days, mu, sigma), 3))