//
// Created by runeebl on 4/26/22.
//

#include "graphs.h"

#include "utils.h"

#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <fstream>


#define MAX_LINE_LEN 256


void print_graph(const Graph& graph)
{
    for (uint i = 0; i < graph.size(); ++i)
    {
        std::cout << i + 1 << ": ";
        for (auto arc : graph[i])
        {
            std::cout << arc.to + 1 << " ";
        }
        std::cout << std::endl;
    }
}


Graph parse_graph_c9(const std::string& graph_file_name)
{
    std::ifstream stream;

    stream.open(graph_file_name);

    char first_char;
    char line[MAX_LINE_LEN];

    uint num_nodes = 0;
    uint num_arcs = 0;

    uint from = 0;
    uint to = 0;
    uint weight = 0;

    Graph graph;

    while (stream.get(first_char)) {
//        std::cout << first_char << std::endl;

        stream.getline(line, sizeof(line));

        switch (first_char)
        {
            case 'p':
                // assumes this only happens once
                sscanf(line, " sp %u %u\n", &num_nodes, &num_arcs);

                graph.resize(num_nodes);

                break;
            case 'a':
                sscanf(line, " %u %u %u\n", &from, &to, &weight);

                // input index starts at 1, decrement to start at 0
                from -= 1;
                to -= 1;

                // Dimacs graph file has from->to and to->from as separate edges,
                // so only need to add one edge per line
                graph[from].push_back({to, weight});

                break;
            default:
                break;
        }
    }

    stream.close();

    return graph;
}


std::vector<std::pair<int, int>> parse_coords_c9(const std::string& graph_file_name)
{
    std::ifstream stream;

    stream.open(graph_file_name);

    char first_char;
    char line[MAX_LINE_LEN];

    uint num_nodes = 0;

    uint node = 0;
    int x = 0;
    int y = 0;

    std::vector<std::pair<int, int>> coords;

    while (stream.get(first_char)) {

        stream.getline(line, sizeof(line));

        switch (first_char)
        {
            case 'p':
                // assumes this only happens once
                sscanf(line, " aux sp co %u\n", &num_nodes);

                coords.resize(num_nodes);

                break;
            case 'v':
                sscanf(line, " %u %i %i\n", &node, &x, &y);

                // input index starts at 1, decrement to start at 0
                node -= 1;

                coords[node] = {x, y};

                break;
            default:
                break;
        }
    }

    stream.close();

    return coords;
}


Graph random_graph(uint n, double p, uint min_weight, uint max_weight)
{
    Graph graph(n);

    boost::random::uniform_01<> dist_p;
    boost::random::uniform_int_distribution<uint> dist_weight(min_weight, max_weight);

    for (uint i = 0; i < n - 1; ++i)
    {
        for (uint j = i + 1; j < n; ++j)
        {
            if (dist_p(gen) < p)
            {
                uint weight = dist_weight(gen);
                graph[i].push_back({j, weight});
                graph[j].push_back({i, weight});
            }
        }
    }

    return graph;
}


Graph random_graph_fixed_m(uint n, uint m, uint min_weight, uint max_weight)
{
    Graph graph(n);

    auto pairs = sample_increasing_pairs_without_replacement(n, m);

    boost::random::uniform_int_distribution<uint> dist_weight(min_weight, max_weight);

    for (auto& pair : pairs)
    {
        uint weight = dist_weight(gen);
        graph[pair.first].push_back({pair.second, weight});
        graph[pair.second].push_back({pair.first, weight});
    }

    return graph;
}


std::pair<Graph, std::vector<std::pair<int, int>>> random_geo_graph(uint n, double p, int max_coord, const std::function<uint(uint, uint, const std::vector<std::pair<int, int>>&)>& dist)
{
    Graph graph(n);

    boost::random::uniform_01<> dist_p;
    boost::random::uniform_int_distribution<int> dist_coord(0, max_coord);

    std::vector<std::pair<int, int>> coords(n);

    for (uint i = 0; i < n; ++i)
    {
        int x = dist_coord(gen);
        int y = dist_coord(gen);
        coords[i] = {x, y};
    }

    for (uint i = 0; i < n - 1; ++i)
    {
        for (uint j = i + 1; j < n; ++j)
        {
            if (dist_p(gen) < p)
            {
                uint weight = dist(i, j, coords);
                graph[i].push_back({j, weight});
                graph[j].push_back({i, weight});
            }
        }
    }

    return {graph, coords};
}
