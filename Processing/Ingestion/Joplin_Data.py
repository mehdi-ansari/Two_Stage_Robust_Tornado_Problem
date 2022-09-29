# -*- coding: utf-8 -*-
"""
Created on Sat Feb 12 16:11:20 2022

@author: mehdi
"""
import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
from .ReferenceData import ReferenceData

class JoplinData:
    
    def __init__(self, ROOT_DIR, number_of_clusters=None):
        self.coordinates={}
        self.first_stage_dislocation={}
        self.cost_retrofitting={}
        self.second_stage_dislocation={}
        self.cost_recovery={}
        
        self.retrofitting_strategies = ('Do_nothing', 'R1', 'R2', 'R3')
        self.recovery_strategies = ('Do_nothing', 'Recover')
        
        self.block_data = pd.read_excel(str(ROOT_DIR)+'/Data/Joplin_Data/ListExcel.xlsx')
        self.INCORE_data = pd.read_csv(str(ROOT_DIR)+'/Data/Joplin_Data/IN-CORE_2ev2_housingunitallocation_1238.csv')
        
        self.repair_time = 60
        
        self.integrate_data()
        self.df_locations = self.df
        
        #If there is number_of_clusters argument, find Kmean centers for locations instead of block ids
        if number_of_clusters is not None:
            self.find_Kmeans(num_clusters=number_of_clusters)
            self.df_locations = self.df_Kmeans
            
        
        self.generate_parameters(self.df_locations)



    def integrate_data(self):
        excel_file = self.block_data
        csv_file = self.INCORE_data
        
        block_area = csv_file.groupby(["blockid"])["gsq_foot"].sum()
        block_popu = csv_file.groupby(["blockid"])["numprec"].sum()


        self.df = pd.DataFrame(columns=['Block id', 'Longitud', 'Latitude', 'Cost nothing', 'Cost R1',
               'Cost R2', 'Cost R3', 'Area', 'Population', 'ProbR0', 'ProbR1', 'ProbR2', 'ProbR3'])
        
        
        for i in range(0,len(excel_file),5):
            
            self.df = self.df.append({'Block id': int(excel_file["Block id"][i]), 
                            'Longitud': excel_file["Longitud"][i], 
                            'Latitude': excel_file["Latitude"][i], 
                            'Cost nothing': round(excel_file["Cost nothing"][i]), 
                            'Cost R1': round(excel_file["Cost R1"][i]),
                            'Cost R2': round(excel_file["Cost R2"][i]), 
                            'Cost R3': round(excel_file["Cost R3"][i]), 
                            'Area': block_area[int(excel_file["Block id"][i])], 
                            'Population': int(block_popu[int(excel_file["Block id"][i])]), 
                            'ProbR0': (1-sum(excel_file["P"][i+k] for k in range(1,5)), excel_file["P"][i+1], excel_file["P"][i+2], excel_file["P"][i+3], excel_file["P"][i+4]), 
                            'ProbR1': (1-sum(excel_file["R1"][i+k] for k in range(1,5)), excel_file["R1"][i+1], excel_file["R1"][i+2], excel_file["R1"][i+3], excel_file["R1"][i+4]), 
                            'ProbR2': (1-sum(excel_file["R2"][i+k] for k in range(1,5)), excel_file["R2"][i+1], excel_file["R2"][i+2], excel_file["R2"][i+3], excel_file["R2"][i+4]), 
                            'ProbR3': (1-sum(excel_file["R3"][i+k] for k in range(1,5)), excel_file["R3"][i+1], excel_file["R3"][i+2], excel_file["R3"][i+3], excel_file["R3"][i+4])}, ignore_index=True)
        
        self.df = self.df.drop_duplicates(subset='Block id')
        self.df = self.df.sort_values('Longitud', ignore_index=True)

    
    def generate_parameters(self, location_dataframe):
        '''
        This function generates the input parameters for the optimization model, namely:
            1. location coordinates dictionary
            2. first stage dislocation dictionary
            3. cost of retrofitting dictionary
            4. second stage dislocation dictionary
            5. cost of recovery dictionary
        
        To see how this parameters are calculated, please review the associated journal paper.

        Parameters
        ----------
        location_dataframe : Pandas Dataframe
            DESCRIPTION. This argument might be either the whole location data or clusters, but must be in the same format

        Returns
        -------
        None.

        '''
        RD = ReferenceData(elapsed_days=self.repair_time)
        self.population = {}
        for index,rows in location_dataframe.iterrows():
            self.coordinates[rows['Block id']] = (round(rows['Longitud']*54.6,3), round(rows['Latitude']*69,3))
            self.first_stage_dislocation[rows['Block id']] = {i: 0 for i in self.retrofitting_strategies}
            self.cost_retrofitting[rows['Block id']] = {'Do_nothing': rows['Cost nothing'], 'R1': rows['Cost R1'], 'R2': rows['Cost R2'], 'R3': rows['Cost R3']}
            self.second_stage_dislocation[rows['Block id']] = {self.retrofitting_strategies[i]: {'Do_nothing': round(((1/np.dot(RD.P_NotRecoverd, rows['ProbR'+str(i)])) + 1)/2 * np.dot(RD.P_NotRecoverd, rows['ProbR'+str(i)])*rows['Population']), 
                                                                                                 'Recover': round(np.dot(RD.P_NotRecoverd, rows['ProbR'+str(i)])*rows['Population'])} 
                                                                                       for i in range(len(self.retrofitting_strategies))}
            '''self.second_stage_dislocation[rows['Block id']] = {self.retrofitting_strategies[i]: {'Do_nothing': round(np.dot(RD.P_NotRecoverd, rows['ProbR'+str(i)])*rows['Population']), 
                                                                                                 'Recover': round(np.dot(RD.P_NotRecoverd, rows['ProbR'+str(i)])*rows['Population'])} 
                                                                                       for i in range(len(self.retrofitting_strategies))}'''
            self.cost_recovery[rows['Block id']] = {self.retrofitting_strategies[i]: {'Do_nothing': 0, 'Recover': round(np.dot(RD.costPercentage, rows['ProbR'+str(i)])*rows['Area']*RD.cost_per_ft2)} 
                                                                                       for i in range(len(self.retrofitting_strategies))}
            
            self.population[rows['Block id']] = rows['Population']
        

    def find_Kmeans(self, **kwargs):
        '''
        This function finds K cluster centers for Joplin blocks and update parameters based on new location centers
        Each center represents the avg/sum of block parameters that are inside the cluster

        Parameters
        ----------
        **kwargs : int
            DESCRIPTION. the value K for kmeans

        Returns
        -------
        A dataframe contaning cluster data

        '''
        locations = np.array(self.df[['Longitud','Latitude']])
        kmeans = KMeans(n_clusters=kwargs['num_clusters'], random_state=0).fit(locations)
        
        self.df["Center Label"] = kmeans.labels_
        cluster_groupby = self.df.groupby(["Center Label"])
        #getAgroup = df1.get_group(50)
        
        
        #New DataFrame with clusters
        self.df_Kmeans = pd.DataFrame(columns=['Block id', 'Longitud', 'Latitude', 'Cost nothing', 'Cost R1',
               'Cost R2', 'Cost R3', 'Area', 'Population', 'ProbR0', 'ProbR1', 'ProbR2', 'ProbR3'])
        
        for region in range(len(cluster_groupby)):
            Agroup = cluster_groupby.get_group(region)
            self.df_Kmeans = self.df_Kmeans.append({'Block id': int(region+1), 
                            'Longitud': round(kmeans.cluster_centers_[region][0],3), 
                            'Latitude': round(kmeans.cluster_centers_[region][1],3), 
                            'Cost nothing': round(sum(Agroup['Cost nothing'])), 
                            'Cost R1': round(sum(Agroup['Cost R1'])),
                            'Cost R2': round(sum(Agroup['Cost R2'])), 
                            'Cost R3': round(sum(Agroup['Cost R3'])), 
                            'Area': int(sum(Agroup['Area'])), 
                            'Population': int(sum(Agroup['Population'])), 
                            'ProbR0': (np.mean([rows["ProbR0"][0] for index,rows in Agroup.iterrows()]), 
                                       np.mean([rows["ProbR0"][1] for index,rows in Agroup.iterrows()]), 
                                       np.mean([rows["ProbR0"][2] for index,rows in Agroup.iterrows()]), 
                                       np.mean([rows["ProbR0"][3] for index,rows in Agroup.iterrows()]), 
                                       np.mean([rows["ProbR0"][4] for index,rows in Agroup.iterrows()])), 
                            'ProbR1': (np.mean([rows["ProbR1"][0] for index,rows in Agroup.iterrows()]), 
                                       np.mean([rows["ProbR1"][1] for index,rows in Agroup.iterrows()]), 
                                       np.mean([rows["ProbR1"][2] for index,rows in Agroup.iterrows()]), 
                                       np.mean([rows["ProbR1"][3] for index,rows in Agroup.iterrows()]), 
                                       np.mean([rows["ProbR1"][4] for index,rows in Agroup.iterrows()])), 
                            'ProbR2': (np.mean([rows["ProbR2"][0] for index,rows in Agroup.iterrows()]), 
                                       np.mean([rows["ProbR2"][1] for index,rows in Agroup.iterrows()]), 
                                       np.mean([rows["ProbR2"][2] for index,rows in Agroup.iterrows()]), 
                                       np.mean([rows["ProbR2"][3] for index,rows in Agroup.iterrows()]), 
                                       np.mean([rows["ProbR2"][4] for index,rows in Agroup.iterrows()])), 
                            'ProbR3': (np.mean([rows["ProbR3"][0] for index,rows in Agroup.iterrows()]), 
                                       np.mean([rows["ProbR3"][1] for index,rows in Agroup.iterrows()]), 
                                       np.mean([rows["ProbR3"][2] for index,rows in Agroup.iterrows()]), 
                                       np.mean([rows["ProbR3"][3] for index,rows in Agroup.iterrows()]), 
                                       np.mean([rows["ProbR3"][4] for index,rows in Agroup.iterrows()]))}, 
                             ignore_index=True)
        
        self.df_Kmeans = self.df_Kmeans.sort_values('Longitud', ignore_index=True)