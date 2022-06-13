//
// Created by runeebl on 4/26/22.
//

#ifndef BACHELOR_REPORT_EXPERIMENTS_H
#define BACHELOR_REPORT_EXPERIMENTS_H

#include "graphs.h"
#include "algorithms.h"

#include <cstdlib>
#include <fstream>
#include <vector>


void benchmark_oracle(Graph& graph, const DistanceOracleFactory& distanceOracleFactory, const std::vector<std::pair<uint, uint>>& samples, const std::vector<uint>& true_distances, uint reps, std::ofstream& csv_stream, std::string graph_type, bool path_query);


void benchmark_oracles(Graph& graph, const std::vector<DistanceOracleFactory>& oracles, uint num_samples, uint reps, std::ofstream& csv_stream, std::string graph_type, bool path_query);


void benchmark_n_p(uint n_min, uint n_max, uint n_step, uint min_weight, uint max_weight, double p, uint k_max, uint k_step, uint num_samples, uint reps, std::ofstream& csv_stream, bool path_query);


void benchmark_n_fixed_m(uint m, uint n_min, uint n_max, uint n_step, uint min_weight, uint max_weight, uint k_max, uint k_step, uint num_samples, uint reps, std::ofstream& csv_stream, bool path_query);


void benchmark_m_fixed_n(uint n, uint m_min, uint m_max, uint m_step, uint min_weight, uint max_weight, double p, uint k_max, uint k_step, uint num_samples, uint reps, std::ofstream& csv_stream, bool path_query);


void benchmark_geo(uint n_min, uint n_max, uint n_step, int max_coord, double p, uint k_max, uint k_step, uint num_samples, uint reps, std::ofstream& csv_stream, bool path_query);


void benchmark_dimacs(uint k_max, uint k_step, uint num_samples, uint reps, std::ofstream& csv_stream, bool path_query);


#endif //BACHELOR_REPORT_EXPERIMENTS_H
