# -*- coding: utf-8 -*-
"""
Created on Fri Sep  2 15:49:26 2022

@author: mehdi
"""
class UserInput:
    def __init__(self, ROOT_DIR):
        #Read text file to input user parameter:
        self.input_dict = {'number_of_clusters': None}
        with open(str(ROOT_DIR)+"/userInputParameter.txt") as user_file:
            for line in user_file:
                (key, value) = line.split()
                self.input_dict[key] = value