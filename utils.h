//
// Created by runeebl on 4/26/22.
//

#ifndef BACHELOR_REPORT_UTILS_H
#define BACHELOR_REPORT_UTILS_H


#include <boost/random/mersenne_twister.hpp>
#include "graphs.h"


extern boost::random::mt19937 gen;


std::vector<std::pair<uint, uint>> sample_pairs_without_replacement(uint n, uint k);


std::vector<std::pair<uint, uint>> sample_increasing_pairs_without_replacement(uint n, uint k);


uint manhattan_dist(uint a, uint b, const std::vector<std::pair<int, int>>& coords);


double euclidean_dist(uint a, uint b, const std::vector<std::pair<int, int>>& coords);


uint graph_mem_usage(const Graph& graph);


uint graph_size(const Graph& graph);


#endif //BACHELOR_REPORT_UTILS_H
