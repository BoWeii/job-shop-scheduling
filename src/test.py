import pandas as pd
import numpy as np
import time
import copy
import _jss
''' ================= initialization setting ======================'''


pt_tmp=pd.read_excel("JSP_dataset.xlsx",sheet_name="Processing Time",index_col =[0])
ms_tmp=pd.read_excel("JSP_dataset.xlsx",sheet_name="Machines Sequence",index_col =[0])

dfshape=pt_tmp.shape
num_mc=dfshape[1] # number of machines
num_job=dfshape[0] # number of jobs
num_gene=num_mc*num_job # number of genes in a chromosome
# print('pt1=',pt_tmp)
# print('ms1=',ms_tmp)

pt=[list(map(int, pt_tmp.iloc[i])) for i in range(num_job)]
ms=[list(map(int,ms_tmp.iloc[i])) for i in range(num_job)]


# print('pt2=',pt)
# print('ms2=',ms)

# raw_input is used in python 2
population_size=int(input('Please input the size of population: ') or 30) # default value is 30
crossover_rate=float(input('Please input the size of Crossover Rate: ') or 0.8) # default value is 0.8
mutation_rate=float(input('Please input the size of Mutation Rate: ') or 0.2) # default value is 0.2
mutation_selection_rate=float(input('Please input the mutation selection rate: ') or 0.2)
num_mutation_jobs=round(num_gene*mutation_selection_rate)
num_iteration=int(input('Please input number of iteration: ') or 2000) # default value is 2000



# start_time = time.time()

'''==================== main code ==============================='''
'''----- generate initial population -----'''
Tbest=999999999999999
best_list,best_obj=[],[]
population_list=[]
makespan_record=[]
for i in range(population_size):
    nxm_random_num=list(np.random.permutation(num_gene)) # generate a random permutation of 0 to num_job*num_mc-1
    population_list.append(nxm_random_num) # add to the population_list
    for j in range(num_gene):
        population_list[i][j]=population_list[i][j]%num_job # convert to job number format, every job appears m times

# print("************",population_list)

a=_jss.JSS(population_list,population_size,crossover_rate,mutation_rate,mutation_selection_rate,num_mutation_jobs,num_iteration,num_gene)
print(a.num_iteration)
# print("num_gene",num_gene)
a.gene_algorithm()
# print("population_list->",a.population_list)
# print("parent_list->",a.parent_list)
# print("off->",a.offspring_list)
# print("sss->",a.S)





        
