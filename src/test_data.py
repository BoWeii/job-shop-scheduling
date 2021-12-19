import pandas as pd
import numpy as np
import _jss
'''--------plot gantt chart-------'''
import pandas as pd
from chart_studio import plotly as py
import plotly.figure_factory as ff
import datetime
import sys


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
        self.population_list = []
        self.dfshape = self.pt_tmp.shape
        self.num_mc = self.dfshape[1]  # number of machines
        self.num_job = self.dfshape[0]  # number of jobs
        self.num_gene = self.num_mc * self.num_job  # number of genes in a chromosome
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

    def gene_algo(self):
        self.generate_init_population()
        obj = _jss.JSS(self.pt, self.ms, self.population_list,
                       self.population_size, self.crossover_rate,
                       self.mutation_rate, self.mutation_selection_rate,
                       self.num_mutation_jobs, self.num_iteration,
                       self.num_gene, self.num_job, self.num_mc)
        obj.gene_algorithm()
        # print("population_list=",obj.population_list)
        print("S=", obj.S)
        self.sequence_best = obj.sequence_best
        self.Tbest = obj.Tbest
        print("====Tbest=====", self.Tbest)
        print("====makespan_record===", obj.makespan_record)

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
        for m in m_keys:
            for j in j_keys:
                df.append(
                    dict(Task='Machine %s' % (m),
                         Start='2022-01-03 %s' % (str(j_record[(j, m)][0])),
                         Finish='2022-01-03 %s' % (str(j_record[(j, m)][1])),
                         Resource='Job %s' % (j + 1)))

        fig = ff.create_gantt(df,
                              index_col='Resource',
                              show_colorbar=True,
                              group_tasks=True,
                              showgrid_x=True,
                              title='Job shop Schedule')
        fig.show()


def main():
    data_index = (sys.argv)[1]
    jss = JSS(data_index)
    jss.gene_algo()
    jss.plotly()


if __name__ == '__main__':
    main()
