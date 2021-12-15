#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <vector>
#include <time.h>


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
        v[i] = i ;

    // While vector has elements
    // get a random number from the vector and print it
    while (v.size()) {
        // cout << getNum(v) << " ";
        S.push_back(getNum(v));
    }
}

void gen_cutpoint(int size, int* cutpoint) {

    unsigned int seed = time(NULL);
    int num1 = rand_r(&seed) % size;
    int num2 = rand_r(&seed) % size;
    if (num2 < num1) {
        cutpoint[0] = num2;
        cutpoint[1] = num1;
    }
    else {
        cutpoint[0] = num1;
        cutpoint[1] = num2;
    }

}

class JSS
{
    // friend JSS multiply_naive(JSS const &mat1, JSS const &mat2);
    friend void generateRandom(int n);

public:
    JSS(std::vector<std::vector<int>>population_list2, int population_size2, float crossover_rate2, float mutation_rate2, float mutation_selection_rate2, int num_mutation_jobs2, int num_iteration2, int num_gene2) :
        population_list(population_list2), population_size(population_size2), crossover_rate(crossover_rate2), mutation_rate(mutation_rate2), mutation_selection_rate(mutation_selection_rate2), num_mutation_jobs(num_mutation_jobs2), num_iteration(num_iteration2), num_gene(num_gene2)
    {
    }
    void gene_algorithm();
    void crossover();
    float get_Tbest() const { return Tbest; }
    float get_num_iteration() const { return num_iteration; }
    std::vector<std::vector<int>> get_population_list() { return population_list; }
    std::vector<std::vector<int>> get_parent_list() { return parent_list; }
    std::vector<std::vector<int>> get_offspring_list() { return offspring_list; }
    std::vector<int> get_test() { return test; }
    std::vector<int> get_S() { return S; }
    std::vector<int> get_sequence_best() { return sequence_best; }

    float Tbest = 0;
    std::vector<int> sequence_best, S, test;
    std::vector<std::vector<int>>population_list, parent_list, offspring_list;
    int population_size;
    float crossover_rate;
    float mutation_rate;
    float mutation_selection_rate;
    int num_mutation_jobs;
    int num_iteration;
    int num_gene;

};

void JSS::gene_algorithm() {
    for (int n = 0;n < 1;n++) {
        // long Tbest_now = 99999999999;
        this->crossover();
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
            int cutpoint[2];
            gen_cutpoint(num_gene, cutpoint);
            for (int i = cutpoint[0];i < cutpoint[1];i++) {
                // std::cout << "S[2 * m]=" << S[2 * m] << " 2 * m" << 2 * m << "\n";
                // printf("i=%d,%d \n", i, offspring_list[S[2 * m]][i]);
                offspring_list[S[2 * m]][i] = population_list[S[2 * m]][i];
                offspring_list[S[2 * m + 1]][i] = population_list[S[2 * m + 1]][i];
            }
        }
    }


}

PYBIND11_MODULE(_jss, m)
{
    py::class_<JSS>(m, "JSS")
        .def(py::init<std::vector<std::vector<int>>, int, float, float, float, int, int, int>())
        .def("gene_algorithm", &JSS::gene_algorithm)
        .def_property("num_iteration", &JSS::get_num_iteration, nullptr)
        .def_property("population_list", &JSS::get_population_list, nullptr)
        .def_property("parent_list", &JSS::get_parent_list, nullptr)
        .def_property("test", &JSS::get_test, nullptr)
        .def_property("S", &JSS::get_S, nullptr)
        .def_property("offspring_list", &JSS::get_offspring_list, nullptr);
    // .def_property("population_list", &JSS::get_population_list, nullptr)
    // .def_property("sequence_best", &JSS::get_sequence_best, nullptr)
    // .def_property("Tbest", &JSS::get_Tbest, nullptr);
}