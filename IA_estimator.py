import base64
import requests
import os
clear = lambda: os.system('cls')
clear()
import time
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
from scipy.optimize import fsolve
start_time = time.time()

def coupled_exponential(N: int,
                   kappa: float = 0.0,
                   mu: float = 0.0,
                   scale: float = 1.0
                   ) -> np.ndarray:
    u = np.array([(i - 0.4) / (N + 0.2) for i in range(1, N)])
    
    data = mu + (scale / kappa) * ((1 - u) ** (-kappa) - 1)
    return data

N = 1000 
kappa = 0.25 
mu = 0.0  
scale = 0.5 

Data = coupled_exponential(N, kappa, mu, scale)

# Calculate the ML parameters
ML= stats.genpareto.fit(Data, loc=0)

data_ML= coupled_exponential(N, ML[0], mu, ML[2])


def average_least_absolute_deviation(data, Data_estimator):
    sorted_data = np.sort(data)
    sorted_Data_estimator = np.sort(Data_estimator)
    absolute_deviations = np.abs(sorted_data - sorted_Data_estimator)
    avg_lad = np.mean(absolute_deviations)
    return avg_lad

avg_lad = average_least_absolute_deviation(Data, data_ML)

#=================================================
avg_lad_IA_list=[]
scale_estimator_list=[]
kappa_estimator_list=[]
number_list_1=[]
i3_list = []
step_1=[]

def fun():
    def process_and_sort(data, subset_size):
        subset_indices = np.random.choice(np.random.permutation(len(data)), size=len(data) - len(data) % subset_size, replace=False)
        subset = data[subset_indices]
        matrix = np.reshape(subset, (-1, subset_size))
    
        row_distance = np.amax(matrix, axis=1) - np.amin(matrix, axis=1)
        sorted_indices = np.argsort(row_distance)
        matrix.sort(axis=1)
        sorted_matrix = matrix[sorted_indices]
    
        return sorted_matrix

    def calculate_x_median(matrix):
        return np.median(np.abs(matrix), axis=1)

    def process_and_calculate(data, subset_sizes):
        x_results = []
        processed_matrices = []
    
        for size in subset_sizes:
            processed_matrix = process_and_sort(data, size)
            processed_matrices.append(processed_matrix)
        
            x_values = calculate_x_median(processed_matrix)
            x_results.append(x_values)
    
        return processed_matrices, x_results
#subset_sizes of pairs and triplets
    subset_sizes = [2, 3]

    processed_matrices, x_results = process_and_calculate(Data, subset_sizes)

    median_pairs = x_results[0]  # Get x values with subset size 2
    median_triplets=x_results[1] # Get x values with subset size 3
    processed_matrices_pairs=processed_matrices[0]
    processed_matrices_triplets=processed_matrices[1]
#============================
    condition = True
    step = 1
    flag = 1
    nn0=10
    nnf=500
#-----------
    while nn0 < nnf:
      median_pairs_subset=median_pairs[0:math.floor(nn0)]
      m=0.0
      for ii in range(nn0):
          m=m+median_pairs_subset[ii];
      first_moment=(1/(nn0))*m;
      scale_estimator=2*first_moment

      min_avg_lad_IA = float('inf')  # Initialize with a large value
      min_sample_size = None  # Initialize the sample size

      for i3 in range(10, 333):
          if i3 <= len(median_triplets):
            segment = median_triplets[:i3]
            second_moment = stats.moment(segment, moment=2)
              
            kappa_estimator = (2 * scale_estimator ** 2 / (3 * second_moment)) - 3

            data_estimator=coupled_exponential(N, kappa_estimator, mu, scale_estimator)
        
            avg_lad_IA = average_least_absolute_deviation(Data, data_estimator)
            
            avg_lad_IA_list.append(avg_lad_IA)
            scale_estimator_list.append(scale_estimator)
            kappa_estimator_list.append(kappa_estimator)
            i3_list.append(i3)
            step_1.append(step)
            # Track the minimum avg_lad_IA and its corresponding sample size
            if avg_lad_IA < min_avg_lad_IA:
              min_avg_lad_IA = avg_lad_IA
              min_sample_size = i3  # Update the minimum sample size

            print('Iteration-%d, nn0 = %d, scale_estimator = %0.6f, avg_lad_IA = %0.6f, kappa_estimator = %0.6f, i3= %d, first_moment = %0.6f, second_moment = %0.6f' %
              (step, nn0, scale_estimator, avg_lad_IA, kappa_estimator, i3, first_moment, second_moment))

            
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
       min_i3 = i3_list[min_val_idx]

       print('Iteration-%d, nn0 = %d, scale_estimator = %0.6f, kappa_estimator = %0.6f, i3 = %d, avg_lad_IA = %0.6f' % (min_val_idx + 1, min_sample_size, scale_estimator_list[min_val_idx], kappa_estimator_list[min_val_idx], min_i3, result))
       
       print("Average Least Absolute Deviation_ML:", avg_lad)
       print('kappa_esti_ML',ML[0])
       print('scale_esti_ML',ML[2])

       return dic_vals
    else:
          print('\nNot Convergent.')
          return False
#=======================================
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
