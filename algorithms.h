//
// Created by runeebl on 4/26/22.
//

#ifndef BACHELOR_REPORT_ALGORITHMS_H
#define BACHELOR_REPORT_ALGORITHMS_H

#include "graphs.h"

#include <vector>
#include <climits>
#include <cstdlib>
#include <memory>


struct sp_result
{
    uint d;
    uint p;
};


#define PARENT(i) (((i) - 1) / 2)

#define LEFT(i) (2 * (i) + 1)

#define RIGHT(i) (2 * (i) + 2)


extern uint decrease_keys;
extern uint extract_mins;

struct heap_node
{
    uint id;
    uint key;
};


struct min_heap
{
    std::vector<heap_node> A;
    std::vector<uint> p;

    uint get_key(uint id)
    {
        uint i = p[id];
        if (i == UINT_MAX)
        {
            // Already popped, don't relax
            return 0;
        } else
        {
            return A[p[id]].key;
        }
    }

    void min_heapify(uint i)
    {
        uint l = LEFT(i);
        uint r = RIGHT(i);

        uint min;

        if (l < A.size() && A[l].key < A[i].key)
        {
            min = l;
        } else
        {
            min = i;
        }
        if (r < A.size() && A[r].key < A[min].key)
        {
            min = r;
        }
        if (min != i)
        {
            // Exchange A[i] with A[min
            std::swap(A[i], A[min]);
            std::swap(p[A[i].id], p[A[min].id]);
            min_heapify(min);
        }
    }

    heap_node extract_min()
    {
        extract_mins += 1;

        heap_node min = A[0];
        A[0] = A[A.size() - 1];
        A.pop_back();

        p[A[0].id] = 0;
        p[min.id] = UINT_MAX;

        min_heapify(0);
        return min;
    }

    void decrease_key(uint id, uint key)
    {
        decrease_keys += 1;

        uint i = p[id];

        A[i].key = key;

        uint parent = PARENT(i);

        while (i > 0 && A[parent].key > A[i].key)
        {
            std::swap(A[i], A[parent]);
            std::swap(p[A[i].id], p[A[parent].id]);
            i = parent;
            parent = PARENT(i);
        }
    }

    explicit min_heap(const std::vector<heap_node>& nodes) : A(nodes), p(nodes.size())
    {
        for (uint i = 0; i < p.size(); ++i)
        {
            p[i] = i;
        }
        for (uint i = (A.size() - 1) / 2; i < A.size(); --i)
        {
            min_heapify(i);
        }
    }
};


struct min_heap_map
{
    std::vector<heap_node> A;
    std::unordered_map<uint, uint> p;

    uint get_key(uint id)
    {
        auto find = p.find(id);

        if (find == p.end())
        {
            insert_inf(id);
            return UINT_MAX;
        }

        uint i = find->second;

        if (i == UINT_MAX)
        {
            // Already popped, don't relax
            return 0;
        } else
        {
            return A[i].key;
        }
    }

    void min_heapify(uint i)
    {
        uint l = LEFT(i);
        uint r = RIGHT(i);

        uint min;

        if (l < A.size() && A[l].key < A[i].key)
        {
            min = l;
        } else
        {
            min = i;
        }
        if (r < A.size() && A[r].key < A[min].key)
        {
            min = r;
        }
        if (min != i)
        {
            // Exchange A[i] with A[min]
            std::swap(A[i], A[min]);
            std::swap(p[A[i].id], p[A[min].id]);
            min_heapify(min);
        }
    }

    heap_node extract_min()
    {
        extract_mins += 1;

        heap_node min = A[0];
        A[0] = A[A.size() - 1];
        A.pop_back();

        p[A[0].id] = 0;
        p[min.id] = UINT_MAX;

        min_heapify(0);
        return min;
    }

    void decrease_key(uint id, uint key)
    {
        decrease_keys += 1;

        auto find = p.find(id);
        uint i;

        if (find == p.end() || find->second == UINT_MAX)
        {
            insert_inf(id);
            i = A.size() - 1;
        } else
        {
            i = find->second;
        }

        A[i].key = key;

        uint parent = PARENT(i);

        while (i > 0 && A[parent].key > A[i].key)
        {
            std::swap(A[i], A[parent]);
            std::swap(p[A[i].id], p[A[parent].id]);
            i = parent;
            parent = PARENT(i);
        }
    }

    void insert_inf(uint id)
    {
        p[id] = A.size();
        A.push_back({id, UINT_MAX});
    }

    void increase_key(uint id, uint key)
    {
        auto find = p.find(id);
        uint i;

        if (find == p.end() || find->second == UINT_MAX)
        {
            insert_inf(id);
            i = A.size() - 1;
        } else
        {
            i = find->second;
        }

        A[i].key = key;

        while (true)
        {
            // Find minimum of parent and children
            uint min = key;
            uint arg_min = i;

            if (LEFT(i) < A.size() && A[LEFT(i)].key < min)
            {
                min = A[LEFT(i)].key;
                arg_min = LEFT(i);
            }
            if (RIGHT(i) < A.size() && A[RIGHT(i)].key < min)
            {
                min = A[RIGHT(i)].key;
                arg_min = RIGHT(i);
            }

            if (arg_min == i)
            {
                // Node is lesser than both children
                return;
            }
            else
            {
                // Min should be the parent, swap and loop from swapped location
                std::swap(A[i], A[arg_min]);
                std::swap(p[A[i].id], p[A[arg_min].id]);
                i = arg_min;
            }
        }
    }

    explicit min_heap_map(const std::vector<heap_node>& nodes) : A(nodes)
    {
        for (uint i = 0; i < nodes.size(); ++i)
        {
            p[nodes[i].id] = i;
        }
        for (uint i = (A.size() - 1) / 2; i < A.size(); --i)
        {
            min_heapify(i);
        }
    }
};


struct index_tree_node
{
    index_tree_node* left = nullptr;
    index_tree_node* right = nullptr;
    uint start;
    uint stop;
    uint dju_max = UINT_MAX;
    uint dju_arg_max = UINT_MAX;

    explicit index_tree_node(uint start, uint stop) : start(start), stop(stop) {};
    ~index_tree_node()
    {
        delete left;
        delete right;
    }
};


struct j_index_tree_node
{
    uint j = UINT_MAX;
    std::unique_ptr<j_index_tree_node> left = nullptr;
    std::unique_ptr<j_index_tree_node> right = nullptr;
    j_index_tree_node() = default;

    uint size() const
    {
        return 1 + (left ? left->size() : 0) + (right ? right->size() : 0);
    }
};


struct shortest_path_tree_node
{
    uint parent;
    uint parent_dist;
    uint root_dist;
};


class IDistanceOracle
{
public:
    virtual uint distance_query(uint from, uint to) = 0;

    virtual uint path_query(uint from, uint to) = 0;

    virtual uint get_mem_usage() = 0;

    virtual uint size() = 0;

    virtual ~IDistanceOracle() = default;

    virtual std::string get_name() = 0;

    virtual double get_theoretical_max_stretch() = 0;
};


typedef std::function<std::unique_ptr<IDistanceOracle>(Graph& graph)> DistanceOracleFactory;


class IDijkstraOracle : public IDistanceOracle
{
protected:
    const Graph& graph;

public:
    explicit IDijkstraOracle(Graph &graph) : graph(graph) {};

    uint path_query(uint from, uint to) override {
        return distance_query(from, to);
    };

    uint get_mem_usage() override;

    uint size() override;

    double get_theoretical_max_stretch() override {
        return 1;
    };
};


class DijkstraPairOracle : public IDijkstraOracle
{
public:
    explicit DijkstraPairOracle(Graph &graph);

    uint distance_query(uint from, uint to) override;

    std::string get_name() override {
        return "DijkstraPairOracle";
    };
};


class DijkstraPairMapOracle : public IDijkstraOracle
{
public:
    explicit DijkstraPairMapOracle(Graph &graph);

    uint distance_query(uint from, uint to) override;

    std::string get_name() override {
        return "DijkstraPairMapOracle";
    };
};


class DijkstraBidirectionalOracle : public IDijkstraOracle
{
public:
    explicit DijkstraBidirectionalOracle(Graph &graph);

    uint distance_query(uint from, uint to) override;

    std::string get_name() override {
        return "DijkstraBidirectionalOracle";
    };
};


class ThorupZwickOracle : public IDistanceOracle
{
private:
    static std::vector<sp_result> dijkstra_pi(const Graph& graph, uint from);

    std::unordered_map<uint, shortest_path_tree_node> dijkstra_max(const Graph& graph, uint from, uint i);

public:
    const uint k;

    std::vector<std::vector<sp_result>> Adp;
    std::vector<std::unordered_map<uint, shortest_path_tree_node>> C;
    std::vector<std::unordered_map<uint, uint>> B;

    uint thorup_zwick_query(uint from, uint to, uint i);

    uint thorup_zwick_query_explicit(uint from, uint to, uint i);

    ThorupZwickOracle(Graph& graph, uint k);

    uint distance_query(uint from, uint to) override;

    uint path_query(uint from, uint to) override;

    std::string get_name() override {
        return "ThorupZwickOracle (k=" + std::to_string(k) + ")";
    };

    double get_theoretical_max_stretch() override {
        return 2 * k - 1;
    };

    uint get_mem_usage() override;

    uint size() override;
};


class WulffNilsenOracle : public ThorupZwickOracle
{
private:
    std::vector<j_index_tree_node> j_index_trees;

    uint bdist(uint from, uint to, uint i1, uint i2, const j_index_tree_node& j_index_tree_node);

    uint bdist_explicit(uint from, uint to, uint i1, uint i2, const j_index_tree_node& j_index_tree_node);

    void create_index_tree(index_tree_node& current_node, const uint u, uint current, uint i1, uint i2);

    static uint get_j(const index_tree_node* LCA, uint start, uint stop);

    static void create_j_index_tree(const index_tree_node * const index_tree_root, uint k, j_index_tree_node& current_node, uint i1, uint i2);

public:
    WulffNilsenOracle(Graph& graph, uint k);

    uint distance_query(uint from, uint to) override;

    uint path_query(uint from, uint to) override;

    std::string get_name() override {
        return "WulffNilsenOracle (k=" + std::to_string(k) + ")";
    };

    uint get_mem_usage() override;

    uint size() override;
};


class AStarGeoOracle : public IDijkstraOracle
{
private:
    const double stretch;
    const std::function<uint(uint from, uint to)> h;

public:
    AStarGeoOracle(Graph& graph, double stretch, std::function<uint(uint from, uint to)>  h);

    uint distance_query(uint from, uint to) override;

    std::string get_name() override {
        return "AStarGeoOracle (stretch=" + std::to_string(stretch) + ")";
    };

    double get_theoretical_max_stretch() override {
        return stretch;
    };

    uint get_mem_usage() override;

    uint size() override;
};


template <class BaseOracle>
class AStarOracle : public IDijkstraOracle
{
private:
    BaseOracle base_oracle;
    const double stretch;

public:
    AStarOracle(Graph& graph, uint k, double stretch);

    uint distance_query(uint from, uint to) override;

    std::string get_name() override {
        return "AStarOracle (stretch=" + std::to_string(stretch) + ", base_oracle=" + base_oracle.get_name() + ")";
    };

    double get_theoretical_max_stretch() override {
        return stretch;
    };

    uint get_mem_usage() override;

    uint size() override;
};


template <class QueryOracle>
class SpannerOracle : public IDistanceOracle
{
private:
    std::string spannerCreatorName;

    double spannerCreatorMaxTheoreticalStretch;

    Graph spanner;

    QueryOracle queryOracle;
public:
    SpannerOracle(Graph& graph, uint k);

    uint distance_query(uint from, uint to) override;

    uint path_query(uint from, uint to) override;

    std::string get_name() override {
        return "SpannerOracle (SpannerCreator=" + spannerCreatorName + ", QueryOracle=" + queryOracle.get_name() + ")";
    };

    double get_theoretical_max_stretch() override {
        return spannerCreatorMaxTheoreticalStretch * queryOracle.get_theoretical_max_stretch();
    };

    uint get_mem_usage() override;

    uint size() override;
};


std::vector<uint> dijkstra_full(const Graph& graph, uint from);


uint dijkstra_pair(const Graph& graph, uint from, uint to);


uint dijkstra_pair_map(const Graph& graph, uint from, uint to);


uint dijkstra_bidirectional(const Graph& graph, uint from, uint to);


std::vector<sp_result> dijkstra_pi(const Graph& graph, uint from);


std::unordered_map<uint, shortest_path_tree_node> dijkstra_max(const Graph& graph, uint from, const std::vector<sp_result>& Adpi, std::vector<std::unordered_map<uint, uint>>& B);


uint thorup_zwick_query_explicit(uint from, uint to, uint i, const std::vector<std::vector<sp_result>>& Adp, const std::vector<std::unordered_map<uint, shortest_path_tree_node>>& C, const std::vector<std::unordered_map<uint, uint>>& B);


uint bdist_explicit(uint from, uint to, uint i1, uint i2, const std::vector<std::vector<sp_result>>& Adp, const std::vector<std::unordered_map<uint, shortest_path_tree_node>>& C, const std::vector<std::unordered_map<uint, uint>>& B, const j_index_tree_node& j_index_tree_node);


std::tuple<std::vector<std::vector<sp_result>>, std::vector<std::unordered_map<uint, shortest_path_tree_node>>, std::vector<std::unordered_map<uint, uint>>, std::vector<j_index_tree_node>, double> a_star_preprocess(Graph& graph, uint k);


uint a_star(const Graph& graph, uint from, uint to, const std::function<uint(uint, uint)>& h);


uint a_star_BPMX(const Graph& graph, uint from, uint to, const std::function<uint(uint, uint)>& h);


void oracle_spanner_to_graph(const std::vector<std::unordered_map<uint, shortest_path_tree_node>>& C, Graph& graph);


#endif //BACHELOR_REPORT_ALGORITHMS_H
