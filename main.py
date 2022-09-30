# -*- coding: utf-8 -*-
"""
Created on Mon May 24 14:15:20 2021

@author: mehdi
"""
from multiprocessing import Process
from config.definitions import ROOT_DIR
from Processing.Ingestion.Parameter import Parameter
from solverEngines.CCGAlgorithm import CCGAlgorithm
from Processing.Postprocessing.Results import Results


def main():
    
    modelParameter = Parameter(ROOT_DIR)
    CCG = CCGAlgorithm(modelParameter)
    CCG.run()
    
    results = Results(CCG, ROOT_DIR)
    results.print_results()
    results.make_file()

if __name__ == '__main__':
    
    program = Process(target=main(), args=())
    program.start()
    program.join(timeout=3600)
    program.terminate()
    
    if program.is_alive():
        print("Complete!")
