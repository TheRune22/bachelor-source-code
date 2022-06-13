//
// Created by runeebl on 4/26/22.
//

#include "algorithms.h"

#include "graphs.h"
#include "utils.h"

#include <cstdlib>
#include <utility>
#include <vector>
#include <climits>
#include <boost/random/uniform_01.hpp>


uint decrease_keys = 0;
uint extract_mins = 0;


uint IDijkstraOracle::get_mem_usage()
{
    return sizeof(IDijkstraOracle) + graph_mem_usage(graph);
}


uint IDijkstraOracle::size()
{
    return graph_size(graph);
}


std::vector<uint> dijkstra_full(const Graph& graph, uint from)
{
    std::vector<uint> result(graph.size(), UINT_MAX);

    std::vector<heap_node> nodes(graph.size());

    for (uint i = 0; i < nodes.size(); ++i)
    {
        nodes[i] = {i, UINT_MAX};
    }

    nodes[from].key = 0;

    min_heap heap(nodes);

    while (!heap.A.empty())
    {
        heap_node top = heap.extract_min();

        if (top.key == UINT_MAX)
        {
            break;
        }

        result[top.id] = top.key;

        for (arc arc : graph[top.id])
        {
            // Relax
            uint new_d = top.key + arc.weight;

            if (heap.get_key(arc.to) > new_d)
            {
                heap.decrease_key(arc.to, new_d);
            }
        }
    }

    return result;
}


uint dijkstra_pair(const Graph& graph, uint from, uint to)
{
    std::vector<heap_node> nodes(graph.size());

    for (uint i = 0; i < nodes.size(); ++i)
    {
        nodes[i] = {i, UINT_MAX};
    }

    nodes[from].key = 0;

    min_heap heap(nodes);

    while (!heap.A.empty())
    {
        heap_node top = heap.extract_min();

        if (top.id == to || top.key == UINT_MAX)
        {
            return top.key;
        }

        for (arc arc : graph[top.id])
        {
            // Relax
            uint new_d = top.key + arc.weight;

            if (heap.get_key(arc.to) > new_d)
            {
                heap.decrease_key(arc.to, new_d);
            }
        }
    }

    return UINT_MAX;
}


DijkstraPairOracle::DijkstraPairOracle(Graph &graph) : IDijkstraOracle(graph) {}


uint DijkstraPairOracle::distance_query(uint from, uint to)
{
    return dijkstra_pair(graph, from, to);
}


uint dijkstra_pair_map(const Graph& graph, uint from, uint to)
{
    std::vector<heap_node> nodes(1, {from, 0});
    min_heap_map heap(nodes);
// ---------------------------------------------------------------------------------------------------------------------

    while (!heap.A.empty())
    {
        heap_node top = heap.extract_min();

        if (top.id == to || top.key == UINT_MAX)
        {
            return top.key;
        }

        for (arc arc : graph[top.id])
        {
            // Relax
            uint new_d = top.key + arc.weight;

            if (heap.get_key(arc.to) > new_d)
            {
                heap.decrease_key(arc.to, new_d);
            }
        }
    }

    return UINT_MAX;
}


DijkstraPairMapOracle::DijkstraPairMapOracle(Graph &graph) : IDijkstraOracle(graph) {}


uint DijkstraPairMapOracle::distance_query(uint from, uint to)
{
    return dijkstra_pair_map(graph, from, to);
}


uint dijkstra_bidirectional(const Graph& graph, uint from, uint to)
{
    std::vector<heap_node> forward_nodes(1, {from, 0});
    min_heap_map forward_heap(forward_nodes);

    std::vector<heap_node> backward_nodes(1, {to, 0});
    min_heap_map backward_heap(backward_nodes);

    std::unordered_map<uint, uint> forward_results;
    std::unordered_map<uint, uint> backward_results;

    uint min_seen = UINT_MAX;

    if (from == to)
    {
        return 0;
    }

    while (!forward_heap.A.empty() && !backward_heap.A.empty())
    {
        heap_node top_forward = forward_heap.extract_min();
        forward_results[top_forward.id] = top_forward.key;

        heap_node top_backward = backward_heap.extract_min();
        backward_results[top_backward.id] = top_backward.key;

        if (top_forward.key + top_backward.key >= min_seen)
        {
            break;
        }

        // Forward
        if (top_forward.key == UINT_MAX)
        {
            break;
        }

        for (arc arc : graph[top_forward.id])
        {
            // Relax
            uint new_d = top_forward.key + arc.weight;

            if (forward_heap.get_key(arc.to) > new_d)
            {
                forward_heap.decrease_key(arc.to, new_d);
                auto backward_dist = backward_results.find(arc.to);
                uint min_candidate;
                if (backward_dist != backward_results.end() && (min_candidate = backward_dist->second + new_d) < min_seen)
                {
                    min_seen = min_candidate;
                }
            }
        }

        // Backward
        if (top_backward.key == UINT_MAX)
        {
            break;
        }

        for (arc arc : graph[top_backward.id])
        {
            // Relax
            uint new_d = top_backward.key + arc.weight;

            if (backward_heap.get_key(arc.to) > new_d)
            {
                backward_heap.decrease_key(arc.to, new_d);
                auto forward_dist = forward_results.find(arc.to);
                uint min_candidate;
                if (forward_dist != forward_results.end() && (min_candidate = forward_dist->second + new_d) < min_seen)
                {
                    min_seen = min_candidate;
                }
            }
        }
    }

    return min_seen;
}


DijkstraBidirectionalOracle::DijkstraBidirectionalOracle(Graph &graph) : IDijkstraOracle(graph) {}


uint DijkstraBidirectionalOracle::distance_query(uint from, uint to)
{
    return dijkstra_bidirectional(graph, from, to);
}


std::vector<sp_result> ThorupZwickOracle::dijkstra_pi(const Graph& graph, uint from)
{
    std::vector<sp_result> result(graph.size(), {UINT_MAX, from});
    std::vector<heap_node> nodes(graph.size());

    for (uint i = 0; i < nodes.size(); ++i)
    {
        nodes[i] = {i, UINT_MAX};
    }

    nodes[from].key = 0;

    min_heap heap(nodes);

    while (!heap.A.empty())
    {
        heap_node top = heap.extract_min();

        if (top.key == UINT_MAX)
        {
            break;
        }

        result[top.id].d = top.key;

        for (arc arc : graph[top.id])
        {
            // Relax
            uint new_d = top.key + arc.weight;

            if (heap.get_key(arc.to) > new_d)
            {
                heap.decrease_key(arc.to, new_d);
                if (top.id == from) {
                    result[arc.to].p = arc.to;
                } else {
                    result[arc.to].p = result[top.id].p;
                }
            }
        }
    }

    return result;
}


std::unordered_map<uint, shortest_path_tree_node> ThorupZwickOracle::dijkstra_max(const Graph& graph, uint from, uint i)
{
    const std::vector<sp_result>& Adpi = Adp[i];

    std::unordered_map<uint, shortest_path_tree_node> Cw;

    Cw[from] = {from, 0, 0};

    std::vector<heap_node> nodes(1, {from, 0});
    min_heap_map heap(nodes);

    while (!heap.A.empty())
    {
        heap_node top = heap.extract_min();

        if (top.key == UINT_MAX)
        {
            break;
        }

        B[top.id][from] = top.key;

        for (arc arc : graph[top.id])
        {
            // Relax
            uint new_d = top.key + arc.weight;
            if (new_d >= Adpi[arc.to].d)
            {
                continue;
            }

            if (heap.get_key(arc.to) > new_d)
            {
                heap.decrease_key(arc.to, new_d);
                Cw[arc.to] = {top.id, arc.weight, new_d};
            }
        }
    }

    return Cw;
}


uint ThorupZwickOracle::thorup_zwick_query_explicit(uint from, uint to, uint i)
{
    uint u = from;
    uint v = to;
    uint w = Adp[i][u].p;

    std::unordered_map<uint, uint>::const_iterator Bvw;

    while ((Bvw = B[v].find(w)) == B[v].end()) {
        i += 1;

        if (i > k) {
            return UINT_MAX;
        }

        std::swap(u, v);
        w = Adp[i][u].p;
    }

    uint total_dist = 0;
    while (v != u) {
        auto u_node = C[w].find(u)->second;
        auto v_node = C[w].find(v)->second;

        if (u_node.root_dist > v_node.root_dist)
        {
            u = u_node.parent;
            total_dist += u_node.parent_dist;
        }
        else
        {
            v = v_node.parent;
            total_dist += v_node.parent_dist;
        }
    }

    return total_dist;
}


ThorupZwickOracle::ThorupZwickOracle(Graph& graph, uint k) : k(k)
{
    uint n = graph.size();

    // Set up outputs:
    Adp = std::vector<std::vector<sp_result>>(k + 1, std::vector<sp_result>(n));

    C = std::vector<std::unordered_map<uint, shortest_path_tree_node>>(n);

    B = std::vector<std::unordered_map<uint, uint>>(n);

    boost::random::uniform_01<> dist;
    double p = pow((double)n,  (- 1. / (double)k));

    // Construct samples
    std::vector<std::vector<uint>> A(k + 1);

    // init A[0]
    A[0].resize(n);
    for (uint i = 0; i < A[0].size(); ++i)
    {
        A[0][i] = i;
    }

    while (A[k - 1].empty())
    {
        for (uint i = 1; i < k; ++i)
        {
            // Ai
            A[i].clear();
            for (auto elm : A[i - 1])
            {
                if (dist(gen) < p)
                {
                    A[i].push_back(elm);
                }
            }
        }
    }

    Adp[k] = std::vector<sp_result>(n, {UINT_MAX, UINT_MAX});

    for (uint i = k - 1; i < k; --i)
    {
        std::vector<arc> s(A[i].size());

        for (uint j = 0; j < s.size(); ++j)
        {
            s[j] = {A[i][j], 0};
        }

        graph.push_back(s);

        // Calculate shortest distances and witnesses
        Adp[i] = dijkstra_pi(graph, graph.size() - 1);

        // Remove "artificial" vertex
        graph.pop_back();
        Adp[i].pop_back();

        for (uint j = 0; j < Adp[i].size(); ++j)
        {
            if (Adp[i][j].d == Adp[i + 1][j].d)
            {
                // Ensure witness property holds
                Adp[i][j].p = Adp[i + 1][j].p;
            }
        }

        std::cout << "Adp[" << i << "] calculated" << std::endl;
        uint elements_left = A[i].size() - A[i + 1].size();

        // for w in A[i] - A[i + 1]
        auto Ai_plus_1_iterator = A[i + 1].begin();
        for (uint Ai_elm : A[i])
        {
            if (Ai_plus_1_iterator == A[i + 1].end() || Ai_elm != *Ai_plus_1_iterator)
            {
                // Element from A[i] was not in A[i + 1]
                C[Ai_elm] = dijkstra_max(graph, Ai_elm, i + 1);

                std::cout << "Elements left: " << --elements_left << "\33[2K\r" << std::flush;
            } else
            {
                // Element from A[i] was also in A[i + 1]
                ++Ai_plus_1_iterator;
            }
        }
    }
}


uint ThorupZwickOracle::thorup_zwick_query(uint from, uint to, uint i)
{
    uint u = from;
    uint v = to;
    uint w = Adp[i][u].p;

    std::unordered_map<uint, uint>::const_iterator Bvw;

    while ((Bvw = B[v].find(w)) == B[v].end()) {
        i += 1;

        if (i > k) {
            return UINT_MAX;
        }

        std::swap(u, v);
        w = Adp[i][u].p;
    }

    return Adp[i][u].d + Bvw->second;
}


uint ThorupZwickOracle::distance_query(uint from, uint to)
{
    return thorup_zwick_query(from, to, 0);
}


uint ThorupZwickOracle::path_query(uint from, uint to)
{
    return thorup_zwick_query_explicit(from, to, 0);
}


uint ThorupZwickOracle::get_mem_usage()
{
    uint total = 0;

    total += sizeof(*this);

    // all inner vectors are same size, no need to loop
    uint Apd_mem = Adp.size() * Adp[0].size() * sizeof Adp[0][0];

    total += Apd_mem;

    uint C_B_mem = 0;

    for (const auto & i : C)
    {
        C_B_mem += sizeof i;

        // heap allocations by unordered map verified using valgrind massif
        C_B_mem += i.size() * sizeof (std::__detail::_Hash_node<__typeof__(*i.begin()), false>);
        C_B_mem += i.bucket_count() * sizeof (std::__detail::_Hash_node_base*);
    }

    for (const auto & i : B)
    {
        C_B_mem += sizeof i;

        // heap allocations by unordered map verified using valgrind massif
        C_B_mem += i.size() * sizeof (std::__detail::_Hash_node<__typeof__(*i.begin()), false>);
        C_B_mem += i.bucket_count() * sizeof (std::__detail::_Hash_node_base*);
    }

    total += C_B_mem;

    return total;
}

uint ThorupZwickOracle::size()
{
    uint total = Adp.size() * Adp[0].size();

    for (const auto & i : C)
    {
        total += i.size();
    }

    for (const auto & i : B)
    {
        total += i.size();
    }

    return total;
}


inline uint get_i(uint i1, uint i2)
{
    uint i = ((i1 + i2) / 2);
    if (i % 2 == 1)
    {
        i += 1;
    }
    return i;
}


uint WulffNilsenOracle::bdist_explicit(uint from, uint to, uint i1, uint i2, const j_index_tree_node &j_index_tree_node)
{
    uint u = from;
    uint v = to;
    uint k = Adp.size() - 1;

    if ((double)(i2 - i1) <= log2((double)k) || (i2 - i1) < 4)
    {
        return thorup_zwick_query_explicit(from, to, i1);
    }

    uint i = get_i(i1, i2);

    if (B[v].find(Adp[j_index_tree_node.j][u].p) == B[v].end() && B[u].find(Adp[j_index_tree_node.j + 1][v].p) == B[u].end())
    {
        return bdist_explicit(u, v, i, i2, *j_index_tree_node.right);
    } else
    {
        return bdist_explicit(u, v, i1, j_index_tree_node.j, *j_index_tree_node.left);
    }
}


inline void update_arg_max(uint& current_max, uint& current_arg_max, uint candidate_max, uint candidate_arg_max)
{
    if (candidate_max > current_max)
    {
        current_max = candidate_max;
        current_arg_max = candidate_arg_max;
    }
}


void WulffNilsenOracle::create_index_tree(index_tree_node& current_node, const uint u, uint current, uint i1, uint i2)
{
    // Stopping condition
    if (i2 - i1 < 3) {
        // Set i1 as default, then try updating with i2 values
        current_node.dju_max = Adp[i1 + 2][u].d - Adp[i1][u].d;
        current_node.dju_arg_max = i1;
        update_arg_max(current_node.dju_max, current_node.dju_arg_max, Adp[i2 + 2][u].d - Adp[i2][u].d, i2);

        return;
    }

    // Calculate i
    uint i = get_i(i1, i2);

    // Add children to tree
    current_node.left = new index_tree_node(i1, i);
    current_node.right = new index_tree_node(i, i2);

    // Run recursively on children
    create_index_tree(*current_node.left, u, LEFT(current), i1, i);
    create_index_tree(*current_node.right, u, RIGHT(current), i, i2);

    // dju_max and dju_arg_max should now be known for children

    // Set left as default, then try updating with right values
    current_node.dju_max = current_node.left->dju_max;
    current_node.dju_arg_max = current_node.left->dju_arg_max;
    update_arg_max(current_node.dju_max, current_node.dju_arg_max, current_node.right->dju_max, current_node.right->dju_arg_max);
}


uint WulffNilsenOracle::get_j(const index_tree_node* LCA, uint start, uint stop)
{
    if (start == stop)
    {
        return start;
    }
    bool LCA_found = false;

    while (!LCA_found)
    {
        if (LCA->left && LCA->left->stop >= stop)
        {
            // Both start and stop are in left subtree
            LCA = LCA->left;
        } else if (LCA->right && LCA->right->start <= start)
        {
            // Both start and stop are in right subtree
            LCA = LCA->right;
        } else
        {
            // Start and stop are in separate subtrees
            LCA_found = true;
        }
    }

    // Check if LCA matches sequence
    if (LCA->start == start && LCA->stop == stop)
    {
        return LCA->dju_arg_max;
    }

    uint dju_max = 0;
    uint dju_arg_max = start;

    index_tree_node* current = LCA->left;
    while (current->start != start)
    {
        if (current->right->start > start)
        {
            // Start is in left subtree
            current = current->left;
            // Right should be considered for finding j
            update_arg_max(dju_max, dju_arg_max, current->right->dju_max, current->right->dju_arg_max);
        }
        else
        {
            // Start is in right subtree
            current = current->right;
        }
    }
    // Endpoint should be considered for finding j
    update_arg_max(dju_max, dju_arg_max, current->dju_max, current->dju_arg_max);

    current = LCA->right;
    while (current->stop != stop)
    {
        if (current->left->stop < stop)
        {
            // Stop is in right subtree
            current = current->right;
            // Left should be considered for finding j
            update_arg_max(dju_max, dju_arg_max, current->left->dju_max, current->left->dju_arg_max);
        }
        else
        {
            // Start is in left subtree
            current = current->left;
        }
    }
    // Endpoint should be considered for finding j
    update_arg_max(dju_max, dju_arg_max, current->dju_max, current->dju_arg_max);

    return dju_arg_max;
}


void WulffNilsenOracle::create_j_index_tree(const index_tree_node * const index_tree_root, uint k, j_index_tree_node& current_node, uint i1, uint i2)
{
    if ((double)(i2 - i1) <= log2((double)k) || (i2 - i1) < 4)
    {
        return;
    }

    uint i = get_i(i1, i2);

    current_node.j = get_j(index_tree_root, i1, i - 2);

    // Recurse on left subtree
    current_node.left = std::make_unique<j_index_tree_node>();
    create_j_index_tree(index_tree_root, k, *current_node.left, i1, current_node.j);

    // Recurse on right subtree
    current_node.right = std::make_unique<j_index_tree_node>();
    create_j_index_tree(index_tree_root, k, *current_node.right, i, i2);
}


uint WulffNilsenOracle::bdist(uint from, uint to, uint i1, uint i2, const j_index_tree_node& j_index_tree_node)
{
    uint u = from;
    uint v = to;
    uint k = Adp.size() - 1;

    if ((double)(i2 - i1) <= log2((double)k) || (i2 - i1) < 4)
    {
        return thorup_zwick_query(from, to, i1);
    }

    uint i = get_i(i1, i2);

    if (B[v].find(Adp[j_index_tree_node.j][u].p) == B[v].end() && B[u].find(Adp[j_index_tree_node.j + 1][v].p) == B[u].end())
    {
        return bdist(u, v, i, i2, *j_index_tree_node.right);
    } else
    {
        return bdist(u, v, i1, j_index_tree_node.j, *j_index_tree_node.left);
    }
}


WulffNilsenOracle::WulffNilsenOracle(Graph &graph, uint k) : ThorupZwickOracle(graph, k)
{
    uint i_max = (k - 3) % 2 == 0 ? k - 3 : k - 4;
    j_index_trees = std::vector<j_index_tree_node>(graph.size());

    for (uint u = 0; u < graph.size(); ++u)
    {
        index_tree_node index_tree_root(0, i_max);
        create_index_tree(index_tree_root, u, 0, 0, i_max);
        create_j_index_tree(&index_tree_root, k, j_index_trees[u], 0, k - 1);
    }
}


uint WulffNilsenOracle::distance_query(uint from, uint to)
{
    return bdist(from, to, 0, k - 1, j_index_trees[from]);
}


uint WulffNilsenOracle::path_query(uint from, uint to)
{
    return bdist_explicit(from, to, 0, k - 1, j_index_trees[from]);
}


uint WulffNilsenOracle::get_mem_usage()
{
    uint total_mem = ThorupZwickOracle::get_mem_usage();

    total_mem += sizeof(j_index_trees);

    uint j_index_tree_mem = 0;

    for (const auto & i: j_index_trees)
    {
        j_index_tree_mem += sizeof (i) * i.size();
    }

    total_mem += j_index_tree_mem;

    return total_mem;
}

uint WulffNilsenOracle::size()
{
    uint total_size = ThorupZwickOracle::size();

    for (const auto & i: j_index_trees)
    {
        total_size += i.size();
    }

    return total_size;
}


uint a_star(const Graph& graph, uint from, uint to, const std::function<uint(uint, uint)>& h)
{
    std::vector<heap_node> nodes(1, {from, 0 + h(from, to)});
    min_heap_map heap(nodes);

    std::unordered_map<uint, uint> gScores;

    while (!heap.A.empty())
    {
        heap_node top = heap.extract_min();
        uint top_gScore = gScores[top.id];

        if (top.id == to || top.key == UINT_MAX)
        {
            return top_gScore;
        }

        for (arc arc : graph[top.id])
        {
            // Relax
            uint new_gScore = top_gScore + arc.weight;

            auto neighbor_g = gScores.find(arc.to);

            if (neighbor_g == gScores.end())
            {
                // Not seen yet, add to map
                gScores[arc.to] = new_gScore;
                heap.decrease_key(arc.to, new_gScore + h(arc.to, to));
            } else if (new_gScore < neighbor_g->second)
            {
                // Already seen, just update
                neighbor_g->second = new_gScore;
                heap.decrease_key(arc.to, new_gScore + h(arc.to, to));
            }
        }
    }

    return UINT_MAX;
}


AStarGeoOracle::AStarGeoOracle(Graph &graph, double stretch,
                               std::function<uint(uint, uint)> h) : IDijkstraOracle(graph), stretch(stretch), h(std::move(h))
{

}


uint AStarGeoOracle::distance_query(uint from, uint to)
{
    return a_star(graph, from, to, [this](uint from, uint to){return h(from, to) * stretch;});
}


uint AStarGeoOracle::get_mem_usage()
{
    // Adding size of coords vector
    return IDijkstraOracle::get_mem_usage() + graph.size() * sizeof(uint) * 2;
}


uint AStarGeoOracle::size()
{
    return IDijkstraOracle::size();
}


template<class BaseOracle>
AStarOracle<BaseOracle>::AStarOracle(Graph &graph, uint k, double stretch) : IDijkstraOracle(graph), base_oracle(graph, k), stretch(stretch) {}


template<class BaseOracle>
uint AStarOracle<BaseOracle>::distance_query(uint from, uint to)
{
    return a_star(graph, from, to, [=, this](uint from, uint to) { return base_oracle.distance_query(from, to) * (stretch / base_oracle.get_theoretical_max_stretch()); });
}


template<class BaseOracle>
uint AStarOracle<BaseOracle>::get_mem_usage()
{
    return IDijkstraOracle::get_mem_usage() + base_oracle.get_mem_usage() + sizeof(*this);
}


template<class BaseOracle>
uint AStarOracle<BaseOracle>::size()
{
    return IDijkstraOracle::size() + base_oracle.size();
}


template class AStarOracle<ThorupZwickOracle>;
template class AStarOracle<WulffNilsenOracle>;


void oracle_spanner_to_graph(const std::vector<std::unordered_map<uint, shortest_path_tree_node>>& C, Graph& graph)
{
    // Use hashmap initially since there may be duplicate edges in C, but duplicates will have the same weight
    auto graph_map = std::vector<std::unordered_map<uint, uint>>(C.size());

    for (const auto& Cw : C)
    {
        for (const auto& [v, node] : Cw)
        {
            // Don't add the root node w from T(w)
            if (node.parent != v)
            {
                graph_map[v][node.parent] = node.parent_dist;
                graph_map[node.parent][v] = node.parent_dist;
            }
        }
    }

    graph.resize(C.size());

    for (uint i = 0; i < C.size(); i++)
    {
        for (const auto& [v, weight] : graph_map[i])
        {
            graph[i].push_back({v, weight});
        }
    }
}


template <class QueryOracle>
SpannerOracle<QueryOracle>::SpannerOracle(Graph &graph, uint k) : queryOracle(spanner) {
    ThorupZwickOracle oracle = ThorupZwickOracle(graph, k);

    oracle_spanner_to_graph(oracle.C, spanner);

    spannerCreatorName = oracle.get_name();

    spannerCreatorMaxTheoreticalStretch = oracle.get_theoretical_max_stretch();
}


template<class QueryOracle>
uint SpannerOracle<QueryOracle>::distance_query(uint from, uint to)
{
    return queryOracle.distance_query(from, to);
}


template<class QueryOracle>
uint SpannerOracle<QueryOracle>::path_query(uint from, uint to)
{
    return queryOracle.path_query(from, to);
}


template<class QueryOracle>
uint SpannerOracle<QueryOracle>::get_mem_usage()
{
    return queryOracle.get_mem_usage();
}


template<class QueryOracle>
uint SpannerOracle<QueryOracle>::size()
{
    return queryOracle.size();
}


template class SpannerOracle<DijkstraPairOracle>;
template class SpannerOracle<DijkstraPairMapOracle>;
template class SpannerOracle<DijkstraBidirectionalOracle>;
