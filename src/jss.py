import pandas as pd
import numpy as np
import time
import copy
import sys
import random
from chart_studio import plotly as py
import plotly.figure_factory as ff
import datetime


class JSS:
    def read_file(self, data_index):
        file_path = "../data/JSP_dataset" + str(data_index) + ".xlsx"
        self.pt_tmp = pd.read_excel(file_path,
                                    sheet_name="Processing Time",
                                    index_col=[0])
        self.ms_tmp = pd.read_excel(file_path,
                                    sheet_name="Machines Sequence",
                                    index_col=[0])

    def __init__(self, data_index):
        self.read_file(data_index)
        self.population_list =[]
        self.makespan_record=[]
        self.dfshape = self.pt_tmp.shape
        self.num_mc = self.dfshape[1]  # number of machines
        self.num_job = self.dfshape[0]  # number of jobs
        self.num_gene = self.num_mc * self.num_job  # number of genes in a chromosome
        self.Tbest = 9999999
        self.pt = [
            list(map(int, self.pt_tmp.iloc[i])) for i in range(self.num_job)
        ]
        self.ms = [
            list(map(int, self.ms_tmp.iloc[i])) for i in range(self.num_job)
        ]
        print()
        self.population_size = int(
            input('Please input the size of population: ')
            or 30)  # default value is 30
        self.crossover_rate = float(
            input('Please input the size of Crossover Rate: ')
            or 0.8)  # default value is 0.8
        self.mutation_rate = float(
            input('Please input the size of Mutation Rate: ')
            or 0.2)  # default value is 0.2
        self.mutation_selection_rate = float(
            input('Please input the mutation selection rate: ') or 0.2)
        self.num_mutation_jobs = round(self.num_gene *
                                       self.mutation_selection_rate)
        self.num_iteration = int(
            input('Please input number of iteration: ')
            or 2000)  # default value is 2000

    def generate_init_population(self):
        for i in range(self.population_size):
            nxm_random_num = list(
                np.random.permutation(self.num_gene)
            )  # generate a random permutation of 0 to num_job*num_mc-1
            self.population_list.append(
                nxm_random_num)  # add to the population_list
            for j in range(self.num_gene):
                self.population_list[i][j] = self.population_list[i][
                    j] % self.num_job  # convert to job number format, every job appears m times

    def crossover(self):
        # two point crossover

        self.parent_list = copy.deepcopy(self.population_list)
        self.offspring_list = copy.deepcopy(self.population_list)

        # generate a random sequence to select the parent chromosome to crossover
        S = list(np.random.permutation(self.population_size))

        for i in range(int(self.population_size/2)):
            crossover_prob = np.random.rand()
            if self.crossover_rate >= crossover_prob:
                parent_1 = self.population_list[S[2*i]][:]
                parent_2 = self.population_list[S[2*i+1]][:]
                child_1 = parent_1[:]
                child_2 = parent_2[:]
                cutpoint = list(np.random.choice(
                    self.num_gene, 2, replace=False))
                cutpoint.sort()
                child_1[cutpoint[0]:cutpoint[1]
                        ] = parent_2[cutpoint[0]:cutpoint[1]]
                child_2[cutpoint[0]:cutpoint[1]
                        ] = parent_1[cutpoint[0]:cutpoint[1]]
                self.offspring_list[S[2*i]] = child_1[:]
                self.offspring_list[S[2*i+1]] = child_2[:]

    def repairment(self):
        for m in range(self.population_size):
            job_count = {}
            # 'larger' record jobs appear in the chromosome more than m times, and 'less' records less than m times.
            larger, less = [], []
            for i in range(self.num_job):
                if i in self.offspring_list[m]:
                    count = self.offspring_list[m].count(i)
                    pos = self.offspring_list[m].index(i)
                    # store the above two values to the job_count dictionary
                    job_count[i] = [count, pos]
                else:
                    count = 0
                    job_count[i] = [count, 0]
                if count > self.num_mc:
                    larger.append(i)
                elif count < self.num_mc:
                    less.append(i)

            for k in range(len(larger)):
                chg_job = larger[k]
                while job_count[chg_job][0] > self.num_mc:
                    for d in range(len(less)):
                        if job_count[less[d]][0] < self.num_mc:
                            self.offspring_list[m][job_count[chg_job]
                                                   [1]] = less[d]
                            job_count[chg_job][1] = self.offspring_list[m].index(
                                chg_job)
                            job_count[chg_job][0] = job_count[chg_job][0]-1
                            job_count[less[d]][0] = job_count[less[d]][0]+1
                        if job_count[chg_job][0] == self.num_mc:
                            break

    def muation(self):
        for m in range(len(self.offspring_list)):
            mutation_prob = np.random.rand()
            if self.mutation_rate >= mutation_prob:
                # chooses the position to mutation
                m_chg = list(np.random.choice(
                    self.num_gene, self.num_mutation_jobs, replace=False))
                # save the value which is on the first mutation position
                # print(self.offspring_list)
                # print(m_chg)
                t_value_last = self.offspring_list[m][m_chg[0]]
                for i in range(self.num_mutation_jobs-1):
                    # displacement
                    self.offspring_list[m][m_chg[i]] = self.offspring_list[m][m_chg[i+1]]

                # move the value of the first mutation position to the last mutation position
                self.offspring_list[m][m_chg[self.num_mutation_jobs-1]
                                       ] = t_value_last

    def calc_fitness(self):
        # calculate the makespan
        # parent and offspring chromosomes combination

        self.total_chromosome = copy.deepcopy(
            self.parent_list)+copy.deepcopy(self.offspring_list)
        self.chrom_fitness, self.chrom_fit = [], []
        self.total_fitness = 0

        for m in range(self.population_size*2):
            j_keys = [j for j in range(self.num_job)]
            key_count = {key: 0 for key in j_keys}
            j_count = {key: 0 for key in j_keys}
            m_keys = [j+1 for j in range(self.num_mc)]
            m_count = {key: 0 for key in m_keys}

            for i in self.total_chromosome[m]:
                gen_t = int(self.pt[i][key_count[i]])
                gen_m = int(self.ms[i][key_count[i]])
                j_count[i] = j_count[i]+gen_t
                m_count[gen_m] = m_count[gen_m]+gen_t

                if m_count[gen_m] < j_count[i]:
                    m_count[gen_m] = j_count[i]
                elif m_count[gen_m] > j_count[i]:
                    j_count[i] = m_count[gen_m]

                key_count[i] = key_count[i]+1

            makespan = max(j_count.values())
            self.chrom_fitness.append(1/makespan)
            self.chrom_fit.append(makespan)
            self.total_fitness = self.total_fitness+self.chrom_fitness[m]

    def selection(self):
        pk, qk = [], []

        for i in range(self.population_size*2):
            pk.append(self.chrom_fitness[i]/self.total_fitness)
        for i in range(self.population_size*2):
            cumulative = 0
            for j in range(0, i+1):
                cumulative = cumulative+pk[j]
            qk.append(cumulative)

        selection_rand = [np.random.rand()
                          for i in range(self.population_size)]

        for i in range(self.population_size):
            if selection_rand[i] <= qk[0]:
                self.population_list[i] = copy.deepcopy(
                    self.total_chromosome[0])
            else:
                for j in range(0, self.population_size*2-1):
                    if selection_rand[i] > qk[j] and selection_rand[i] <= qk[j+1]:
                        self.population_list[i] = copy.deepcopy(
                            self.total_chromosome[j+1])
                        break

    def comparison(self):
        Tbest_now = 99999999
        for i in range(self.population_size*2):
            if self.chrom_fit[i] < Tbest_now:
                # print("chrom_fit[i]=",chrom_fit[i])
                Tbest_now = self.chrom_fit[i]
                sequence_now = copy.deepcopy(self.total_chromosome[i])
        if Tbest_now < self.Tbest:
            self.Tbest = Tbest_now
            self.sequence_best = copy.deepcopy(sequence_now)
            self.makespan_record.append(Tbest_now)

    def gene_algo(self):
        self.generate_init_population()
        i = 0
        while i < self.num_iteration:
            self.crossover()
            self.repairment()
            self.muation()
            self.calc_fitness()
            self.selection()
            self.comparison()
            i=i+1

        print("====Tbest=====", self.Tbest)
        print("====makespan_record===", self.makespan_record)
        return self.sequence_best

    def plotly(self):

        m_keys = [j + 1 for j in range(self.num_mc)]
        j_keys = [j for j in range(self.num_job)]
        key_count = {key: 0 for key in j_keys}
        j_count = {key: 0 for key in j_keys}
        m_count = {key: 0 for key in m_keys}
        j_record = {}
        for i in self.sequence_best:
            gen_t = int(self.pt[i][key_count[i]])
            gen_m = int(self.ms[i][key_count[i]])
            j_count[i] = j_count[i] + gen_t
            m_count[gen_m] = m_count[gen_m] + gen_t

            if m_count[gen_m] < j_count[i]:
                m_count[gen_m] = j_count[i]
            elif m_count[gen_m] > j_count[i]:
                j_count[i] = m_count[gen_m]

            start_time = str(
                datetime.timedelta(seconds=j_count[i] -
                                   self.pt[i][key_count[i]])
            )  # convert seconds to hours, minutes and seconds
            end_time = str(datetime.timedelta(seconds=j_count[i]))

            j_record[(i, gen_m)] = [start_time, end_time]

            key_count[i] = key_count[i] + 1

        df = []
        def r(): return random.randint(0, 255)
        colors = ['#%02X%02X%02X' % (r(), r(), r())]
        for m in m_keys:
            for j in j_keys:
                colors.append('#%02X%02X%02X' % (r(), r(), r()))
                df.append(
                    dict(Task='Machine %s' % (m),
                         Start='2022-01-03 %s' % (str(j_record[(j, m)][0])),
                         Finish='2022-01-03 %s' % (str(j_record[(j, m)][1])),
                         Resource='Job %s' % (j + 1)))

        fig = ff.create_gantt(df,
                              index_col='Resource',
                              colors=colors if self.num_job > 10 else None,
                              show_colorbar=True,
                              group_tasks=True,
                              showgrid_x=True,
                              title='Job shop Schedule')
        fig.show()


def check_answer(sequence, data_size):
    count = {key: 0 for key in range(data_size)}
    for i in sequence:
        count[i] += 1
    for i in count:
        if count[i] != data_size:
            return False
    return True


def working(data_size):
    jss = JSS(int(data_size))
    best_sequence = jss.gene_algo()
    assert check_answer(best_sequence, data_size)
    jss.plotly()


def main():
    data_index = (sys.argv)[1]
    working(int(data_index))


if __name__ == '__main__':
    main()
