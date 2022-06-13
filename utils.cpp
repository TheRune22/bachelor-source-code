//
// Created by runeebl on 4/26/22.
//

#include "utils.h"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include <set>

boost::random::mt19937 gen(2);


std::vector<std::pair<uint, uint>> sample_pairs_without_replacement(uint n, uint k)
{
    // (n choose 2) * 2
    assert(k <= n * (n - 1));

    boost::random::uniform_int_distribution<uint> dist(0, n - 1);

    std::set<std::pair<uint, uint>> samples;

    while (samples.size() < k)
    {
        uint from = dist(gen);
        uint to = dist(gen);

        if (from != to)
        {
            samples.emplace(from, to);
        }
    }

    std::vector<std::pair<uint, uint>> samples_vec(k);
    std::copy(samples.begin(), samples.end(), samples_vec.begin());

    return samples_vec;
}


std::vector<std::pair<uint, uint>> sample_increasing_pairs_without_replacement(uint n, uint k)
{
    // n choose 2
    assert(k <= n * (n - 1) / 2);

    boost::random::uniform_int_distribution<uint> dist(0, n - 1);

    std::set<std::pair<uint, uint>> samples;

    while (samples.size() < k)
    {
        uint from = dist(gen);
        uint to = dist(gen);

        if (from < to)
        {
            samples.emplace(from, to);
        }
        else if (from > to)
        {
            samples.emplace(to, from);
        }
    }

    std::vector<std::pair<uint, uint>> samples_vec(k);
    std::copy(samples.begin(), samples.end(), samples_vec.begin());

    return samples_vec;
}


uint manhattan_dist(uint a, uint b, const std::vector<std::pair<int, int>>& coords)
{
    uint first_diff = coords[a].first < coords[b].first ? coords[b].first - coords[a].first : coords[a].first - coords[b].first;
    uint second_diff = coords[a].second < coords[b].second ? coords[b].second - coords[a].second : coords[a].second - coords[b].second;

    return first_diff + second_diff;
}


double euclidean_dist(uint a, uint b, const std::vector<std::pair<int, int>>& coords)
{
    double first_diff = coords[a].first < coords[b].first ? coords[b].first - coords[a].first : coords[a].first - coords[b].first;
    double second_diff = coords[a].second < coords[b].second ? coords[b].second - coords[a].second : coords[a].second - coords[b].second;

    return std::sqrt(first_diff * first_diff + second_diff * second_diff);
}


uint graph_mem_usage(const Graph& graph)
{
    uint mem_usage = 0;
    mem_usage += graph.size() * sizeof(graph[0]);

    for (auto & i: graph)
    {
        mem_usage += i.size() * sizeof(i[0]);
    }

    return mem_usage;
}


uint graph_size(const Graph& graph)
{
    uint total = 0;

    for (auto & i: graph)
    {
        total += i.size();
    }

    return total;
}
