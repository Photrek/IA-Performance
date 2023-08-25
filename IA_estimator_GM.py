import base64
import requests
import os
clear = lambda: os.system('cls')
clear()
import numpy as np
import random
np.random.seed(10)
import math
import statistics
import pandas as pd
from scipy import stats
from scipy.stats import genpareto
import matplotlib.pyplot as plt
import math as mt 
from scipy.stats.mstats import gmean
from scipy.optimize import fsolve
from scipy.special import digamma
from numpy import euler_gamma
import time 
start_time = time.time()

mu=0
loc=0
excel_file_path = r'F:\IA\data.csv'  # Replace with the path of your actual file
Data_1 = pd.read_csv(excel_file_path)

# Flatten the NumPy array
Data = Data_1.values.flatten()
print('Data',Data)
N=len(Data)


# Calculate the ML parameters
ML= stats.genpareto.fit(Data, loc=0)

data_ML_gene= genpareto.rvs(ML[0], loc=loc, scale=ML[2], size=N)
N_1=len(data_ML_gene)

def average_least_absolute_deviation(data, Data_estimator):
    sorted_data = np.sort(data)
    sorted_Data_estimator = np.sort(Data_estimator)
    absolute_deviations = np.abs(sorted_data - sorted_Data_estimator)
    avg_lad = np.mean(absolute_deviations)
    return avg_lad

avg_lad = average_least_absolute_deviation(Data, data_ML_gene)


geometric_mean1 = gmean(Data)
#=================================
avg_lad_IA_list=[]
kappa_estimator_list=[]
scale_estimator_list=[]
number_list_1=[]
step_1=[]
def fun():
    def process_pairs(data):
        random_pairs = data[np.random.choice(np.random.permutation(len(data)), size=len(data)-len(data)%2, replace=False)]
        pairs = np.reshape(random_pairs, (-1, 2))

        D = np.amax(pairs, axis=1) - np.amin(pairs, axis=1)
        ix = np.argsort(D) 
        minD = D.sort()

        sorted_pairs = pairs[ix, :]
        return sorted_pairs

    sorted_pairs = process_pairs(Data)
#============================
    condition = True
    step = 1
    flag = 1
    nn0=10
    nnf=len(sorted_pairs)
#-----------

    while nn0 < nnf:
        choose_random_pairs=sorted_pairs[0:math.floor(nn0),:]
        pairs_median=abs(np.median(choose_random_pairs, axis=1))
        m=0.0
        for ii in range(nn0):
            m=m+pairs_median[ii];
        first_moment=(1/(nn0))*m;

        #sigma estimate
        scale_estimator=2*first_moment

        # kappa estimate using GM
        def digamma_H(a):
            return digamma(a + 1) + euler_gamma

        def fun(kappa_estimator):
            return (scale_estimator/kappa_estimator)*np.exp(-digamma_H((1/kappa_estimator)-1))-geometric_mean1
        x0 = 0.09

        sol = fsolve(fun, x0)
        kappa_estimator = sol[0]

        data_IA_gene=genpareto.rvs(kappa_estimator, loc=loc, scale=scale_estimator, size=N)

        #Find LAD
        avg_lad_IA = average_least_absolute_deviation(Data, data_IA_gene)

        avg_lad_IA_list.append(avg_lad_IA)
        kappa_estimator_list.append(kappa_estimator)
        scale_estimator_list.append(scale_estimator)
        step_1.append(step)
        print('Iteration-%d, nn0 = %d, kappa_estimator = %0.6f, scale_estimator = %0.6f,avg_lad_IA=%0.8f, first_moment=%0.6f' % (step, nn0, kappa_estimator, scale_estimator,avg_lad_IA,first_moment))
      
        number_list_1.append(nn0)
        nn0 = nn0 + 1
        step = step + 1
          
        if nn0 > nnf:
          flag = 0

    if flag==1: 

       avg_lad_IA_list_min=min(avg_lad_IA_list)
       print('avg_lad_IA_list_min',avg_lad_IA_list_min)
       dic_vals={"scale_estimator":scale_estimator,"nn0":nn0}
       min_val = avg_lad_IA_list[0]
       min_val_idx = 0
       for i in range(len(avg_lad_IA_list)):
          if avg_lad_IA_list[i] < min_val:
             min_val = avg_lad_IA_list[i]
             min_val_idx = i
       result= np.min(avg_lad_IA_list)
      
       print('Iteration-%d, nn0 = %d, scale_estimator = %0.6f, kappa_estimator = %0.6f, avg_lad_IA = %0.6f' % (min_val_idx + 1, min_val_idx+10, scale_estimator_list[min_val_idx], kappa_estimator_list[min_val_idx], result))
       
       print("Average Least Absolute Deviation_ML:", avg_lad)
       print('kappa_esti_ML',ML[0])
       print('scale_esti_ML',ML[2])
       
       return dic_vals
    else:
          print('\nNot Convergent.')
          return False 

#============================ 
sig_list= []
number_list=[]
while len(sig_list) <= 0 :
    value = fun()
    if (value):
        sig_list.append(value['scale_estimator'])
        number_list.append(value["nn0"])
sig_mean= np.mean(sig_list)  
number_mean=np.mean(number_list)

end_time = time.time()

# Calculate the elapsed time in seconds
elapsed_time_seconds = end_time - start_time

# Convert elapsed time to minutes and seconds
elapsed_minutes = int(elapsed_time_seconds // 60)
elapsed_seconds = int(elapsed_time_seconds % 60)

# Print the elapsed time in minutes and seconds
print("Elapsed time:", elapsed_minutes, "minutes and", elapsed_seconds, "seconds")

#print(f"\n value list = {sig_list}")
#print(f"\n  iterations list = {number_list}")   
#print(f"the average  is = {sig_mean} with average iterations = {number_mean}")