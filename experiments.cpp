//
// Created by runeebl on 4/26/22.
//

#include "experiments.h"

#include "graphs.h"
#include "algorithms.h"
#include "utils.h"

#include <cstdlib>
#include <boost/random/uniform_01.hpp>
#include <chrono>


void benchmark_oracle(Graph& graph, const DistanceOracleFactory& distanceOracleFactory, const std::vector<std::pair<uint, uint>>& samples, const std::vector<uint>& true_distances, uint reps, std::ofstream& csv_stream, std::string graph_type, bool path_query)
{
    uint n = graph.size();
    uint num_samples = samples.size();

    for (uint i = 0; i < reps; ++i)
    {
        bool stretch_holds_all = true;
        double max_stretch = 0;
        double stretch_sum = 0;
        double total_query_time = 0.0;
        uint total_query_extract_mins = 0;
        uint total_query_decrease_keys = 0;

        std::cout << "Preprocessing" << std::endl;

        extract_mins = 0;
        decrease_keys = 0;
        auto start = std::chrono::high_resolution_clock::now();
        auto oracle = distanceOracleFactory(graph);
        auto stop = std::chrono::high_resolution_clock::now();
        double preprocessing_time =
                (double) duration_cast<std::chrono::nanoseconds>(stop - start).count() / 1000000000;
        uint preprocessing_extract_mins = extract_mins;
        uint preprocessing_decrease_keys = decrease_keys;

        std::cout << "Preprocessing done" << std::endl;

        for (uint j = 0; j < num_samples; ++j)
        {
            uint from = samples[j].first;
            uint to = samples[j].second;

            uint oracle_dist;

            extract_mins = 0;
            decrease_keys = 0;
            start = std::chrono::high_resolution_clock::now();

            if (path_query)
            {
                oracle_dist = oracle->path_query(from, to);
            } else
            {
                oracle_dist = oracle->distance_query(from, to);
            }
            stop = std::chrono::high_resolution_clock::now();
            double duration =
                    (double) duration_cast<std::chrono::nanoseconds>(stop - start).count() / 1000000000;
            total_query_extract_mins += extract_mins;
            total_query_decrease_keys += decrease_keys;

            total_query_time += duration;

            double stretch;
            if ((true_distances[j] == UINT_MAX && oracle_dist == UINT_MAX) ||
                (true_distances[j] == 0 && oracle_dist == 0))
            {
                stretch = 1;
            }
            else
            {
                stretch = (double) oracle_dist / (double) true_distances[j];
            }
            stretch_sum += stretch;
            max_stretch = std::max(stretch, max_stretch);

            stretch_holds_all &= stretch <= oracle->get_theoretical_max_stretch();
        }

        double avg_stretch = NAN;
        double avg_query_time = NAN;
        uint avg_query_extract_mins = UINT_MAX;
        uint avg_query_decrease_keys = UINT_MAX;

        if (num_samples > 0)
        {
            avg_stretch = stretch_sum / (double) num_samples;
            avg_query_time = total_query_time / (double) num_samples;
            avg_query_extract_mins = total_query_extract_mins / num_samples;
            avg_query_decrease_keys = total_query_decrease_keys / num_samples;
        }

        uint m = graph_size(graph) / 2;

        assert(stretch_holds_all);

        std::string algorithm = oracle->get_name();

        if (path_query)
        {
            algorithm += " (path_query)";
        }

        uint oracle_size = oracle->size();
        uint oracle_mem_usage = oracle->get_mem_usage();
        double oracle_max_stretch = oracle->get_theoretical_max_stretch();


        std::cout << std::endl << std::endl;
        std::cout << "n: " << n << std::endl;
        std::cout << "m: " << m << std::endl;
        std::cout << "graph type: " << graph_type << std::endl;
        std::cout << "algorithm: " << algorithm << std::endl;
        std::cout << "Preprocessing time: " << preprocessing_time << " s" << std::endl;
        std::cout << "Preprocessing extract_mins: " << preprocessing_extract_mins << std::endl;
        std::cout << "Preprocessing decrease_keys: " << preprocessing_decrease_keys << std::endl;
        std::cout << "Average distance_query time: " << avg_query_time << " s" << std::endl;
        std::cout << "Average Query extract_mins: " << avg_query_extract_mins << std::endl;
        std::cout << "Average Query decrease_keys: " << avg_query_decrease_keys << std::endl;
        std::cout << "Memory usage: " << oracle_mem_usage << " B" << std::endl;
        std::cout << "Size: " << oracle_size << std::endl;
        std::cout << "Average stretch: " << avg_stretch << std::endl;
        std::cout << "Max stretch: " << max_stretch << std::endl;
        std::cout << "Theoretical max stretch: " << oracle_max_stretch << std::endl;
        std::cout << "Samples: " << num_samples << std::endl;
        std::cout << "Stretch holds for all: " << stretch_holds_all << std::endl;

        csv_stream << n << "\t"
                   << m << "\t"
                   << graph_type << "\t"
                   << algorithm << "\t"
                   << preprocessing_time << "\t"
                   << preprocessing_extract_mins << "\t"
                   << preprocessing_decrease_keys << "\t"
                   << avg_query_time << "\t"
                   << avg_query_extract_mins << "\t"
                   << avg_query_decrease_keys << "\t"
                   << oracle_mem_usage << "\t"
                   << oracle_size << "\t"
                   << avg_stretch << "\t"
                   << max_stretch << "\t"
                   << oracle_max_stretch << "\t"
                   << num_samples << std::endl;
    }
}


void benchmark_oracles(Graph& graph, const std::vector<DistanceOracleFactory>& oracles, uint num_samples, uint reps, std::ofstream& csv_stream, std::string graph_type, bool path_query)
{
    // Get samples
    std::vector<std::pair<uint, uint>> samples = sample_pairs_without_replacement(graph.size(), num_samples);

    // Calculate true distances
    std::vector<uint> true_distances(samples.size());
    for (uint i = 0; i < samples.size(); ++i)
    {
        true_distances[i] = dijkstra_bidirectional(graph, samples[i].first, samples[i].second);
        std::cout << "True distance calculated for " << i << "\33[2K\r" << std::flush;
    }
    std::cout << "True distances calculated" << std::endl;

    // Run benchmarks for each oracle
    for (const auto & oracle: oracles)
    {
        benchmark_oracle(graph, oracle, samples, true_distances, reps, csv_stream, graph_type, path_query);

        // Also run without path query for comparison
        if (path_query)
        {
            benchmark_oracle(graph, oracle, samples, true_distances, reps, csv_stream, graph_type, false);
        }
    }
}


void benchmark_n_p(uint n_min, uint n_max, uint n_step, uint min_weight, uint max_weight, double p, uint k_max, uint k_step, uint num_samples, uint reps, std::ofstream& csv_stream, bool path_query)
{
    for (uint n = n_min; n <= n_max; n += n_step)
    {
        Graph graph = random_graph(n, p, min_weight, max_weight);

        std::vector<DistanceOracleFactory> oracles;

        const uint log_n = static_cast<uint>(std::log(n));

        oracles.emplace_back([](Graph& graph){return std::make_unique<DijkstraPairOracle>(graph);});
        oracles.emplace_back([](Graph& graph){return std::make_unique<DijkstraPairMapOracle>(graph);});
        oracles.emplace_back([](Graph& graph){return std::make_unique<DijkstraBidirectionalOracle>(graph);});
        oracles.emplace_back([](Graph& graph){return std::make_unique<ThorupZwickOracle>(graph, 2);});
        oracles.emplace_back([log_n](Graph& graph){return std::make_unique<ThorupZwickOracle>(graph, log_n);});
        oracles.emplace_back([log_n](Graph& graph){return std::make_unique<WulffNilsenOracle>(graph, log_n);});


        for (uint k = 3; k <= k_max; k += k_step)
        {
            oracles.emplace_back([k](Graph& graph){return std::make_unique<ThorupZwickOracle>(graph, k);});
            oracles.emplace_back([k](Graph& graph){return std::make_unique<WulffNilsenOracle>(graph, k);});
            oracles.emplace_back([k](Graph& graph){return std::make_unique<AStarOracle<ThorupZwickOracle>>(graph, k, 1.0);});
            oracles.emplace_back([k](Graph& graph){return std::make_unique<AStarOracle<WulffNilsenOracle>>(graph, k, 1.0);});
            oracles.emplace_back([k](Graph& graph){return std::make_unique<AStarOracle<ThorupZwickOracle>>(graph, k, k);});
            oracles.emplace_back([k](Graph& graph){return std::make_unique<AStarOracle<WulffNilsenOracle>>(graph, k, k);});
            oracles.emplace_back([k](Graph& graph){return std::make_unique<AStarOracle<ThorupZwickOracle>>(graph, k, 2*k - 1);});
            oracles.emplace_back([k](Graph& graph){return std::make_unique<AStarOracle<WulffNilsenOracle>>(graph, k, 2*k - 1);});
            oracles.emplace_back([k](Graph& graph){return std::make_unique<SpannerOracle<DijkstraPairMapOracle>>(graph, k);});
            oracles.emplace_back([k](Graph& graph){return std::make_unique<SpannerOracle<DijkstraBidirectionalOracle>>(graph, k);});
        }

        benchmark_oracles(graph, oracles, num_samples, reps, csv_stream, "n_p=" + std::to_string(p), path_query);
    }
}


void benchmark_n_fixed_m(uint m, uint n_min, uint n_max, uint n_step, uint min_weight, uint max_weight, uint k_max, uint k_step, uint num_samples, uint reps, std::ofstream& csv_stream, bool path_query)
{
    for (uint n = n_min; n <= n_max; n += n_step)
    {
        Graph graph = random_graph_fixed_m(n, m, min_weight, max_weight);

        std::vector<DistanceOracleFactory> oracles;

        const uint log_n = static_cast<uint>(std::log(n));

        oracles.emplace_back([](Graph& graph){return std::make_unique<DijkstraPairOracle>(graph);});
        oracles.emplace_back([](Graph& graph){return std::make_unique<DijkstraPairMapOracle>(graph);});
        oracles.emplace_back([](Graph& graph){return std::make_unique<DijkstraBidirectionalOracle>(graph);});
        oracles.emplace_back([](Graph& graph){return std::make_unique<ThorupZwickOracle>(graph, 2);});
        oracles.emplace_back([log_n](Graph& graph){return std::make_unique<ThorupZwickOracle>(graph, log_n);});
        oracles.emplace_back([log_n](Graph& graph){return std::make_unique<WulffNilsenOracle>(graph, log_n);});

        for (uint k = 3; k <= k_max; k += k_step)
        {
            oracles.emplace_back([k](Graph& graph){return std::make_unique<ThorupZwickOracle>(graph, k);});
            oracles.emplace_back([k](Graph& graph){return std::make_unique<WulffNilsenOracle>(graph, k);});
            oracles.emplace_back([k](Graph& graph){return std::make_unique<AStarOracle<ThorupZwickOracle>>(graph, k, 1.0);});
            oracles.emplace_back([k](Graph& graph){return std::make_unique<AStarOracle<WulffNilsenOracle>>(graph, k, 1.0);});
            oracles.emplace_back([k](Graph& graph){return std::make_unique<AStarOracle<ThorupZwickOracle>>(graph, k, k);});
            oracles.emplace_back([k](Graph& graph){return std::make_unique<AStarOracle<WulffNilsenOracle>>(graph, k, k);});
            oracles.emplace_back([k](Graph& graph){return std::make_unique<AStarOracle<ThorupZwickOracle>>(graph, k, 2*k - 1);});
            oracles.emplace_back([k](Graph& graph){return std::make_unique<AStarOracle<WulffNilsenOracle>>(graph, k, 2*k - 1);});
            oracles.emplace_back([k](Graph& graph){return std::make_unique<SpannerOracle<DijkstraPairMapOracle>>(graph, k);});
            oracles.emplace_back([k](Graph& graph){return std::make_unique<SpannerOracle<DijkstraBidirectionalOracle>>(graph, k);});
        }

        benchmark_oracles(graph, oracles, num_samples, reps, csv_stream, "n_fixed_m=" + std::to_string(m), path_query);
    }
}


void benchmark_m_fixed_n(uint n, uint m_min, uint m_max, uint m_step, uint min_weight, uint max_weight, double p, uint k_max, uint k_step, uint num_samples, uint reps, std::ofstream& csv_stream, bool path_query)
{
    for (uint m = m_min; m <= m_max; m += m_step)
    {
        Graph graph = random_graph_fixed_m(n, m, min_weight, max_weight);

        std::vector<DistanceOracleFactory> oracles;

        const uint log_n = static_cast<uint>(std::log(n));

        oracles.emplace_back([](Graph& graph){return std::make_unique<DijkstraPairOracle>(graph);});
        oracles.emplace_back([](Graph& graph){return std::make_unique<DijkstraPairMapOracle>(graph);});
        oracles.emplace_back([](Graph& graph){return std::make_unique<DijkstraBidirectionalOracle>(graph);});
        oracles.emplace_back([](Graph& graph){return std::make_unique<ThorupZwickOracle>(graph, 2);});
        oracles.emplace_back([log_n](Graph& graph){return std::make_unique<ThorupZwickOracle>(graph, log_n);});
        oracles.emplace_back([log_n](Graph& graph){return std::make_unique<WulffNilsenOracle>(graph, log_n);});

        for (uint k = 3; k <= k_max; k += k_step)
        {
            oracles.emplace_back([k](Graph& graph){return std::make_unique<ThorupZwickOracle>(graph, k);});
            oracles.emplace_back([k](Graph& graph){return std::make_unique<WulffNilsenOracle>(graph, k);});
            oracles.emplace_back([k](Graph& graph){return std::make_unique<AStarOracle<ThorupZwickOracle>>(graph, k, 1.0);});
            oracles.emplace_back([k](Graph& graph){return std::make_unique<AStarOracle<WulffNilsenOracle>>(graph, k, 1.0);});
            oracles.emplace_back([k](Graph& graph){return std::make_unique<AStarOracle<ThorupZwickOracle>>(graph, k, k);});
            oracles.emplace_back([k](Graph& graph){return std::make_unique<AStarOracle<WulffNilsenOracle>>(graph, k, k);});
            oracles.emplace_back([k](Graph& graph){return std::make_unique<AStarOracle<ThorupZwickOracle>>(graph, k, 2*k - 1);});
            oracles.emplace_back([k](Graph& graph){return std::make_unique<AStarOracle<WulffNilsenOracle>>(graph, k, 2*k - 1);});
            oracles.emplace_back([k](Graph& graph){return std::make_unique<SpannerOracle<DijkstraPairMapOracle>>(graph, k);});
            oracles.emplace_back([k](Graph& graph){return std::make_unique<SpannerOracle<DijkstraBidirectionalOracle>>(graph, k);});
        }

        benchmark_oracles(graph, oracles, num_samples, reps, csv_stream, "m_fixed_n=" + std::to_string(n), path_query);
    }
}


void benchmark_geo(uint n_min, uint n_max, uint n_step, int max_coord, double p, uint k_max, uint k_step, uint num_samples, uint reps, std::ofstream& csv_stream, bool path_query)
{
    for (uint n = n_min; n <= n_max; n += n_step)
    {
        auto& dist_func = euclidean_dist;

        auto geo_graph = random_geo_graph(n, p, max_coord, [](uint from, uint to, const std::vector<std::pair<int, int>>& coords){return std::ceil(dist_func(from, to, coords));});
        auto& graph = geo_graph.first;
        const auto& coords = geo_graph.second;

        std::vector<DistanceOracleFactory> oracles;

        const uint log_n = static_cast<uint>(std::log(n));

        oracles.emplace_back([](Graph& graph){return std::make_unique<DijkstraPairOracle>(graph);});
        oracles.emplace_back([](Graph& graph){return std::make_unique<DijkstraPairMapOracle>(graph);});
        oracles.emplace_back([](Graph& graph){return std::make_unique<DijkstraBidirectionalOracle>(graph);});
        oracles.emplace_back([log_n](Graph& graph){return std::make_unique<ThorupZwickOracle>(graph, log_n);});
        oracles.emplace_back([log_n](Graph& graph){return std::make_unique<WulffNilsenOracle>(graph, log_n);});

        oracles.emplace_back([&coords](Graph& graph) {return std::make_unique<AStarGeoOracle>(graph, 1.0, [&coords](uint from, uint to){return dist_func(from, to, coords);});});

        for (uint k = 5; k <= k_max; k += k_step)
        {
            oracles.emplace_back([k](Graph& graph){return std::make_unique<ThorupZwickOracle>(graph, k);});
            oracles.emplace_back([k](Graph& graph){return std::make_unique<WulffNilsenOracle>(graph, k);});
            oracles.emplace_back([k](Graph& graph){return std::make_unique<AStarOracle<ThorupZwickOracle>>(graph, k, 1.0);});
            oracles.emplace_back([k](Graph& graph){return std::make_unique<AStarOracle<WulffNilsenOracle>>(graph, k, 1.0);});
            oracles.emplace_back([k](Graph& graph){return std::make_unique<AStarOracle<ThorupZwickOracle>>(graph, k, k);});
            oracles.emplace_back([k](Graph& graph){return std::make_unique<AStarOracle<WulffNilsenOracle>>(graph, k, k);});
            oracles.emplace_back([k](Graph& graph){return std::make_unique<AStarOracle<ThorupZwickOracle>>(graph, k, 2*k - 1);});
            oracles.emplace_back([k](Graph& graph){return std::make_unique<AStarOracle<WulffNilsenOracle>>(graph, k, 2*k - 1);});
            oracles.emplace_back([k](Graph& graph){return std::make_unique<SpannerOracle<DijkstraPairMapOracle>>(graph, k);});
            oracles.emplace_back([k](Graph& graph){return std::make_unique<SpannerOracle<DijkstraBidirectionalOracle>>(graph, k);});

            oracles.emplace_back([&coords, k](Graph& graph) {return std::make_unique<AStarGeoOracle>(graph, 2*k - 1, [&coords](uint from, uint to){return dist_func(from, to, coords);});});
        }

        benchmark_oracles(graph, oracles, num_samples, reps, csv_stream, "geo", path_query);
    }
}


void benchmark_dimacs(uint k_max, uint k_step, uint num_samples, uint reps, std::ofstream& csv_stream, bool path_query)
{
    std::vector<std::string> graph_names = {
            "dimacs/USA-road-d.NY",
            "dimacs/USA-road-d.BAY",
            "dimacs/USA-road-d.COL",
            "dimacs/USA-road-d.FLA",
            "dimacs/USA-road-d.NW",
            "dimacs/USA-road-d.NE",
    };

    auto& dist_func = euclidean_dist;

    for (const auto& graph_name: graph_names)
    {
        Graph graph = parse_graph_c9(graph_name + ".gr");

        const auto coords = parse_coords_c9(graph_name + ".co");

        // Modify weights since distances are not greater or equal to euclidean distance
        for (uint from = 0; from < graph.size(); ++from)
        {
            for (auto & arc: graph[from])
            {
                // Need to ceil here to ensure admissible heuristic
                uint dist = static_cast<uint>(std::ceil(dist_func(from, arc.to, coords)));
                arc.weight = dist;
            }
        }

        std::vector<DistanceOracleFactory> oracles;

        const uint log_n = static_cast<uint>(std::log(graph.size()));

        oracles.emplace_back([](Graph& graph){return std::make_unique<DijkstraPairOracle>(graph);});
        oracles.emplace_back([](Graph& graph){return std::make_unique<DijkstraPairMapOracle>(graph);});
        oracles.emplace_back([](Graph& graph){return std::make_unique<DijkstraBidirectionalOracle>(graph);});

        oracles.emplace_back([&coords](Graph& graph) {return std::make_unique<AStarGeoOracle>(graph, 1.0, [&coords](uint from, uint to){return dist_func(from, to, coords);});});

        for (uint k = 5; k <= k_max; k += k_step)
        {
             oracles.emplace_back([k](Graph& graph){return std::make_unique<ThorupZwickOracle>(graph, k);});
             oracles.emplace_back([k](Graph& graph){return std::make_unique<WulffNilsenOracle>(graph, k);});
             oracles.emplace_back([k](Graph& graph){return std::make_unique<AStarOracle<ThorupZwickOracle>>(graph, k, 1.0);});
             oracles.emplace_back([k](Graph& graph){return std::make_unique<AStarOracle<WulffNilsenOracle>>(graph, k, 1.0);});
             oracles.emplace_back([k](Graph& graph){return std::make_unique<AStarOracle<ThorupZwickOracle>>(graph, k, k);});
             oracles.emplace_back([k](Graph& graph){return std::make_unique<AStarOracle<WulffNilsenOracle>>(graph, k, k);});
             oracles.emplace_back([k](Graph& graph){return std::make_unique<AStarOracle<ThorupZwickOracle>>(graph, k, 2*k - 1);});
             oracles.emplace_back([k](Graph& graph){return std::make_unique<AStarOracle<WulffNilsenOracle>>(graph, k, 2*k - 1);});
             oracles.emplace_back([k](Graph& graph){return std::make_unique<SpannerOracle<DijkstraPairMapOracle>>(graph, k);});
             oracles.emplace_back([k](Graph& graph){return std::make_unique<SpannerOracle<DijkstraBidirectionalOracle>>(graph, k);});

            oracles.emplace_back([&coords, k](Graph& graph) {return std::make_unique<AStarGeoOracle>(graph, 2*k - 1, [&coords](uint from, uint to){return dist_func(from, to, coords);});});
        }

        benchmark_oracles(graph, oracles, num_samples, reps, csv_stream, "dimacs", path_query);
    }
}



