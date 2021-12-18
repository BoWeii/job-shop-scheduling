#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <vector>
#include <map>
#include <time.h>
#include <algorithm>    // std::sort
#include <cmath> 


namespace py = pybind11;

int getNum(std::vector<int>& v)
{

    // Size of the vector
    int n = v.size();

    // Generate a random number
    srand(time(NULL));

    // Make sure the number is within
    // the index range
    int index = rand() % n;

    // Get random number from the vector
    int num = v[index];

    // Remove the number from the vector
    std::swap(v[index], v[n - 1]);
    v.pop_back();

    // Return the removed number
    return num;
}

// Function to generate n non-repeating random numbers
void generateRandom(int n, std::vector<int>& S)
{
    std::vector<int> v(n);
    // Fill the vector with the values
    // 0, 1, 2, ..., n-1
    for (int i = 0; i < n; i++)
        v[i] = i;

    // While vector has elements
    // get a random number from the vector and print it
    while (v.size()) {
        // cout << getNum(v) << " ";
        S.push_back(getNum(v));
    }
}

void gen_random_choice(int size, std::vector<int>& _list, int times) {
    //generate non-duplicated #times number in range=0~size    
    for (int i = 0;i < times;i++) {
        int num = 0;
        unsigned int seed = time(NULL);
        do {
            num = rand_r(&seed) % size;
        } while (std::find(_list.begin(), _list.end(), num) != _list.end());
        _list.push_back(num);
    }
    std::sort(_list.begin(), _list.end());

}

bool is_in_list(int num, std::vector<int>& offspring_list) {
    for (auto i = offspring_list.begin();i != offspring_list.end();i++) {
        if (num == *i)
            return true;
    }
    return false;
}

void get_count(std::vector<int>& offspring_list, int& count, int num) {
    for (auto i = offspring_list.begin();i != offspring_list.end();i++) {
        if (num == *i)
            count++;
    }
}

void get_index(std::vector<int>& offspring_list, int& pos, int num) {
    int index = 0;
    for (auto i = offspring_list.begin();i != offspring_list.end();i++) {
        if (*i == num) {
            pos = index;
            break;
        }
        index++;
    }
}

class JSS
{
    // friend JSS multiply_naive(JSS const &mat1, JSS const &mat2);
    friend void generateRandom(int n);

public:
    JSS(std::vector<std::vector<int>>pt2, std::vector<std::vector<int>>ms2, std::vector<std::vector<int>>population_list2, int population_size2, float crossover_rate2, float mutation_rate2, float mutation_selection_rate2, int num_mutation_jobs2, int num_iteration2, int num_gene2, int num_job2, int num_mc2) :
        pt(pt2), ms(ms2), population_list(population_list2), population_size(population_size2), crossover_rate(crossover_rate2), mutation_rate(mutation_rate2), mutation_selection_rate(mutation_selection_rate2), num_mutation_jobs(num_mutation_jobs2), num_iteration(num_iteration2), num_gene(num_gene2), num_job(num_job2), num_mc(num_mc2)
    {
    }
    void gene_algorithm();
    void crossover();
    void repairment();
    void mutatuon();
    void calculate_makespan();
    void selection();
    void comparison(long Tbest_now);
    long get_Tbest() const { return Tbest; }
    float get_num_iteration() const { return num_iteration; }
    std::vector<std::vector<int>> get_population_list() { return population_list; }
    std::vector<std::vector<int>> get_parent_list() { return parent_list; }
    std::vector<std::vector<int>> get_offspring_list() { return offspring_list; }
    std::vector<int> get_test() { return test; }
    std::vector<int> get_S() { return S; }
    std::vector<int> get_sequence_best() { return sequence_best; }

    long Tbest = 999999999999999;
    std::vector<int> sequence_best, S, test, chrom_fit, makespan_record;
    std::vector<std::vector<int>>pt, ms, population_list, parent_list, offspring_list, total_chromosome;
    std::vector<float> chrom_fitness;
    int population_size;
    float crossover_rate;
    float mutation_rate;
    float mutation_selection_rate;
    int num_mutation_jobs;
    int num_iteration;
    int num_gene;
    int num_job;
    int num_mc;
    int total_fitness = 0;
};

void JSS::gene_algorithm() {
    for (int n = 0;n < num_iteration;n++) {
        long Tbest_now = 99999999999;
        this->crossover();
        this->repairment();
        this->mutatuon();
        this->calculate_makespan();
        this->selection();
        this->comparison(Tbest_now);
    }
}

void JSS::crossover() {
    parent_list.clear();
    offspring_list.clear();
    S.clear();
    for (auto it = population_list.begin(); it != population_list.end(); ++it) {
        parent_list.push_back(*it);
        offspring_list.push_back(*it);
    }
    generateRandom(population_size, S); // generate a random sequence (to select the parent chromosome to crossover)
    for (int m = 0;m < (int)(population_size / 2);m++) {
        unsigned int seed = time(NULL);
        float crossover_prob = ((double)rand_r(&seed) / (RAND_MAX));
        if (crossover_rate >= crossover_prob) {
            std::vector<int> cutpoint;
            gen_random_choice(num_gene, cutpoint, 2);
            for (int i = cutpoint[0];i < cutpoint[1];i++) {
                // std::cout << "S[2 * m]=" << S[2 * m] << " 2 * m" << 2 * m << "\n";
                // printf("i=%d,%d \n", i, offspring_list[S[2 * m]][i]);
                offspring_list[S[2 * m]][i] = population_list[S[2 * m]][i];
                offspring_list[S[2 * m + 1]][i] = population_list[S[2 * m + 1]][i];
            }
        }
    }


}

void JSS::repairment() {
    for (int m = 0;m < population_size;m++) {
        // std::cout << "m=" << m << "\n";
        std::map <int, std::vector<int>> job_count;
        std::vector<int> larger, less;
        for (int i = 0;i < num_job;i++) {
            std::vector<int> temp;
            int count = 0, pos = 0;
            if (is_in_list(i, offspring_list[m])) {
                get_count(offspring_list[m], count, i);
                get_index(offspring_list[m], pos, i);
            }
            // std::cout << "count=" << count << "\n";
            temp.push_back(count);
            temp.push_back(pos);
            job_count[i] = temp;
            if (count > num_mc)
                larger.push_back(i);
            else if (count < num_mc)
                less.push_back(i);
        }
        for (auto k = larger.begin();k != larger.end();k++) {
            int chg_job = *k;
            // std::cout << "yee" << chg_job << "\n";
            while (job_count[chg_job][0] > num_mc) {
                for (int d = 0;d < less.size();d++) {
                    if (job_count[less[d]][0] < num_mc) {
                        offspring_list[m][job_count[chg_job][1]] = less[d];
                        int index;
                        get_index(offspring_list[m], index, chg_job);
                        job_count[chg_job][1] = index;
                        job_count[chg_job][0] = job_count[chg_job][0] - 1;
                        job_count[less[d]][0] = job_count[less[d]][0] + 1;
                    }
                    if (job_count[chg_job][0] == num_mc)
                        break;
                }
            }
        }
    }
}
void JSS::mutatuon() {
    int num_mutation_jobs = std::round(num_gene * mutation_selection_rate);
    for (int m = 0;m < offspring_list.size();m++) {
        unsigned int seed = time(NULL);
        float mutation_prob = ((double)rand_r(&seed) / (RAND_MAX));
        if (mutation_rate >= mutation_prob) {
            // printf("mutation_rate >= mutation_prob");
            std::vector<int> m_chg;
            gen_random_choice(num_gene, m_chg, num_mutation_jobs);
            int t_value_last = offspring_list[m][m_chg[0]];
            for (int i = 0;i < num_mutation_jobs - 1;i++) {
                offspring_list[m][m_chg[i]] = offspring_list[m][m_chg[i + 1]];
            }
            offspring_list[m][m_chg[num_mutation_jobs - 1]] = t_value_last;
        }
    }
}

void JSS::calculate_makespan() {
    total_fitness = 0;
    chrom_fitness.clear();
    chrom_fit.clear();
    total_chromosome.clear();
    total_chromosome.insert(total_chromosome.end(), parent_list.begin(), parent_list.end());
    total_chromosome.insert(total_chromosome.end(), offspring_list.begin(), offspring_list.end());
    for (int m = 0;m < population_size * 2;m++) {
        std::vector<int> j_keys, m_keys;
        std::map<int, int> j_count, key_count, m_count;
        for (int j = 0;j < num_job;j++) {
            j_keys.push_back(j);
        }
        for (auto i = j_keys.begin();i != j_keys.end();i++) {
            j_count[*i] = 0;
            key_count[*i] = 0;
        }
        for (int j = 1;j < num_mc + 1;j++) {
            m_keys.push_back(j);
        }
        for (auto i = m_keys.begin();i != m_keys.end();i++) {
            m_count[*i] = 0;
        }
        for (auto i = total_chromosome[m].begin();i != total_chromosome[m].end();i++) {
            int gen_t = int(pt[*i][key_count[*i]]);
            int gen_m = int(ms[*i][key_count[*i]]);
            j_count[*i] = j_count[*i] + gen_t;
            m_count[gen_m] += gen_t;
            if (m_count[gen_m] < j_count[*i])
                m_count[gen_m] = j_count[*i];
            else if (m_count[gen_m] > j_count[*i])
                j_count[*i] = m_count[gen_m];
            key_count[*i] += 1;
        }
        int makespan = -1;
        for (auto it = j_count.begin();it != j_count.end();it++) {
            if (it->second > makespan)
                makespan = it->second;
        }
        chrom_fitness.push_back(1 / makespan);
        chrom_fit.push_back(makespan);
        total_fitness += chrom_fitness[m];
    }
}

void JSS::selection() {
    std::vector<float> pk, qk;
    for (int i = 0;i < population_size * 2;i++) {
        pk.push_back(chrom_fitness[i] / total_fitness);
    }
    for (int i = 0;i < population_size * 2;i++) {
        int cumulative = 0;
        for (int j = 0;j < i + 1;j++) {
            cumulative += pk[j];
        }
        qk.push_back(cumulative);
    }
    std::vector<float> selection_rand;
    unsigned int seed = time(NULL);
    for (int i = 0;i < population_size;i++) {
        selection_rand.push_back(((double)rand_r(&seed) / (RAND_MAX)));
    }
    for (int i = 0;i < population_size;i++) {
        if (selection_rand[i] <= qk[0]) {
            population_list[i] = total_chromosome[0];
        }
        else {
            for (int j = 0;j < population_size * 2 - 1;j++) {
                if (selection_rand[i] > qk[j] && selection_rand[i] <= qk[j + 1]) {
                    population_list[i] = total_chromosome[j + 1];
                    break;
                }
            }
        }
    }
}

void JSS::comparison(long Tbest_now) {
    std::vector<int> sequence_now;
    for (int i = 0;i < population_size * 2;i++) {
        if (chrom_fit[i] < Tbest_now) {
            Tbest_now = chrom_fit[i];
            sequence_now = total_chromosome[i];
        }
    }
    if (Tbest_now <= Tbest) {
        Tbest = Tbest_now;
        sequence_best = sequence_now;
    }
    makespan_record.push_back(Tbest);
}
PYBIND11_MODULE(_jss, m)
{
    py::class_<JSS>(m, "JSS")
        .def(py::init<std::vector<std::vector<int>>, std::vector<std::vector<int>>, std::vector<std::vector<int>>, int, float, float, float, int, int, int, int, int >())
        .def("gene_algorithm", &JSS::gene_algorithm)
        .def_property("num_iteration", &JSS::get_num_iteration, nullptr)
        .def_property("population_list", &JSS::get_population_list, nullptr)
        .def_property("parent_list", &JSS::get_parent_list, nullptr)
        .def_property("test", &JSS::get_test, nullptr)
        .def_property("S", &JSS::get_S, nullptr)
        .def_property("sequence_best", &JSS::get_sequence_best, nullptr)
        .def_property("Tbest", &JSS::get_Tbest, nullptr)
        .def_property("offspring_list", &JSS::get_offspring_list, nullptr);
    // .def_property("population_list", &JSS::get_population_list, nullptr)

}