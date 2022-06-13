//
// Created by runeebl on 4/26/22.
//

#ifndef BACHELOR_REPORT_GRAPHS_H
#define BACHELOR_REPORT_GRAPHS_H

#include <cstdlib>
#include <vector>
#include <string>
#include <functional>


struct arc
{
    uint to;
    uint weight;
};


typedef std::vector<std::vector<arc>> Graph;


void print_graph(const Graph& graph);


Graph parse_graph_c9(const std::string& graph_file_name);


std::vector<std::pair<int, int>> parse_coords_c9(const std::string& graph_file_name);


Graph parse_graph_rg(const std::string& graph_file_name);


Graph random_graph(uint n, double p, uint min_weight, uint max_weight);


Graph random_graph_fixed_m(uint n, uint m, uint min_weight, uint max_weight);


std::pair<Graph, std::vector<std::pair<int, int>>> random_geo_graph(uint n, double p, int max_coord, const std::function<uint(uint, uint, const std::vector<std::pair<int, int>>&)>& dist);


#endif //BACHELOR_REPORT_GRAPHS_H
