# -*- coding: utf-8 -*-
"""
Created on Sat Aug 27 18:06:07 2022

@author: mehdi
"""
from .Joplin_Data import JoplinData
from .Sample_Data import SampleData
from .UserInput import UserInput

class Parameter:
    #def __init__(self, InputDataName: str, number_of_clusters=None, **kwargs):
    def __init__(self, ROOT_DIR):
        self.ROOT_DIR = ROOT_DIR
        self.user = UserInput(ROOT_DIR)
        self.budget = float(self.user.input_dict['budget'])
        self.length = float(self.user.input_dict['length'])
        self.width = float(self.user.input_dict['width'])
        self.subproblem_method = self.user.input_dict['subproblem_formulation']
        
        if self.user.input_dict['number_of_clusters'] is not None:
            self.number_of_clusters = int(self.user.input_dict['number_of_clusters'])
        else:
            self.number_of_clusters = None
        
        if self.user.input_dict['input_data_name'] == 'Joplin':
            self.InputData = JoplinData(ROOT_DIR, self.number_of_clusters)
        
        elif self.user.input_dict['input_data_name'][:6] == 'Sample':
            self.InputData = SampleData(ROOT_DIR, self.user.input_dict['input_data_name'])


'''
#Budget:
budget1 = sum(np.mean(d,axis=1))
budget2 = sum(np.mean(np.mean(c,axis=1),axis=1))
budget = budget1 + budget2/10
budget = round(budget)
print(budget)
budget = input("Budget: ")
budget = int(budget)

#scaling delta and length:
delta = (3/8)    #miles - Ref: Spatial Analyses of the 2011 Joplin Tornado Mortality: Deaths by Interpolated Damage Zones and Location of Victims
length = 10    #22.1miles - Ref: The May 22, 2011 Joplin, Missouri EF5 tornado   https://www.ustornadoes.com/2013/05/22/joplin-missouri-ef5-tornado-may-22-2011/
'''