//
// Created by runeebl on 4/26/22.
//

#include "experiments.h"

#include "graphs.h"
#include "utils.h"
#include "algorithms.h"

#include <cstdlib>
#include <fstream>
#include <sstream>
#include <ctime>
#include <iostream>


int main() {
    std::ofstream csv_stream;

    std::ostringstream filename_stream;

    filename_stream << "benchmark-" << std::time(nullptr) << ".csv";

    csv_stream.open(filename_stream.str());

    csv_stream << "n" << "\t"
               << "m" << "\t"
               << "graph_type" << "\t"
               << "algorithm" << "\t"
               << "preprocessing_time" << "\t"
               << "preprocessing_extract_mins" << "\t"
               << "preprocessing_decrease_keys" << "\t"
               << "avg_query_time" << "\t"
               << "avg_query_extract_mins" << "\t"
               << "avg_query_decrease_keys" << "\t"
               << "mem_usage" << "\t"
               << "size" << "\t"
               << "avg_stretch" << "\t"
               << "max_stretch" << "\t"
               << "theoretical_max_stretch" << "\t"
               << "samples" << "\n";


    uint n_min = 100;
    uint n_max = 3000;
    uint n_step = 100;

    double p = 0.1;
    uint min_weight = 1;
    uint max_weight = 1000000;

    uint k_max = 20;
    uint k_step = 2;

    uint num_samples = 50;
    uint repetitions = 10;

    benchmark_n_p(n_min, n_max, n_step, min_weight, max_weight, p, k_max, k_step, num_samples, repetitions, csv_stream, true);

    uint m = 100000;

    benchmark_n_fixed_m(m, n_min, n_max, n_step, min_weight, max_weight, k_max, k_step, num_samples, repetitions, csv_stream, false);

    uint n = 1000;

    benchmark_m_fixed_n(n, n_min, n_max, n_step, min_weight, max_weight, p, k_max, k_step, num_samples, repetitions, csv_stream, false);

    int max_coord = 500000;

    benchmark_geo(n_min, n_max, n_step, max_coord, p, k_max, k_step, num_samples, repetitions, csv_stream, false);

    benchmark_dimacs(k_max, k_step, num_samples, repetitions, csv_stream, false);

    csv_stream.close();

    return EXIT_SUCCESS;
}
