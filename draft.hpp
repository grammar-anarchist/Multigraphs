#include <algorithm>
#include <cmath>
#include <future>
#include <limits>
#include <queue>
#include <vector>
#include <utility>
#include <unordered_set>

#include "homology.h"
#include "testers.h"
#include "edge_inflater.h"
#include "simplex.h"
#include "hash.h"
#include "link_homologies.h"

class CliqueMultiComplex {
    Infestation graph;
    std::vector<std::vector<Simplex>> complex;
    std::vector<Vertex> vertices;
    std::vector<std::vector<Homology>> homologies;
    std::vector<std::unordered_map<std::set<int>, int>> search;

    void find_all_cliques(
            std::set<int> &curr,
            std::vector<int> &candidates,
            std::unordered_set<int> &out,
            double time,
            int a,
            int b,
            size_t prev_ind
    );

    // a < b
    void change_topology(int new_edge, int a, int b, double time, std::set<int> &curr) {
        std::vector<int> candidates;
        std::unordered_set<int> out;
        for (int i = 0; i != a; ++i) {
            if (graph.edge_count(i, a) && graph.edge_count(i, b))
                candidates.push_back(i);
        }
        for (int i = a + 1; i != b; ++i) {
            if (graph.edge_count(a, i) && graph.edge_count(i, b))
                candidates.push_back(i);
        }
        for (int i = b + 1; i != vertices.size(); ++i) {
            if (graph.edge_count(a, i) && graph.edge_count(b, i))
                candidates.push_back(i);
        }
        find_all_cliques(curr, candidates, out, time, a, b, new_edge);
        compute_nmu(a, b, curr, time);
    }

    // a < b
    void compute_nmu(int a, int b, std::set<int> &edge, double time) {
        std::queue<int> order;
        complex[0][search[0][edge]].compute_edge_nmu(complex, order, time, graph(a, b));
        int dim = 1;
        while (!order.empty()) {
            size_t curr_dim_sz = order.size();
            for (size_t i = 0; i != curr_dim_sz; ++i) {
                complex[dim][order.front()].compute_nmu(complex, order, time);
                order.pop();
            }
            ++dim;
        }
    }

    void reduce(size_t ind, int prime) {
        for (size_t i = 0; i != complex[ind].size(); ++i) {
            Homology res = complex[ind][i].reduce(complex[ind], complex[ind - 1], i, prime);
            if (res.is_homology()) {
                homologies[ind].push_back(res);
            }
        }
    }

    void check_for_immortal(size_t ind) {
        std::vector<Simplex> &columns = complex[ind];
        for (size_t i = 0; i != columns.size(); ++i) {
            if (columns[i].no_pair()) {
                homologies[ind + 1].emplace_back(columns[i].birth());
            }
        }
    }

public:
    CliqueMultiComplex(): complex(1), homologies(1), search(1) {}
    void add_vertex(double time) {
        vertices.emplace_back(unpaired, time);
    }
    // two unequal indeces
    void add_edge(int a, int b, double time) {
        // std::cout << a << " " << b << "\n";
        if (a == b) {
            return;
        }
        if (a > b) {
            std::swap(a, b);
        }
        // a < b
        std::set<int> edge = {a, b};
        if (++graph(a, b) == 1) {
            graph.edge_birth(a, b) = time;
            int curr_ind = complex[0].size();
            if (vertices[b].pair == unpaired) {
                vertices[b].pair = curr_ind;
                complex[0].emplace_back(a, b, time, b);
                if (time > vertices[b].time) {
                    homologies[0].emplace_back(vertices[b].time, time);
                }
            }
            else {
                complex[0].emplace_back(a, b, time);
            }
            search[0][edge] = curr_ind;
            change_topology(curr_ind, a, b, time, edge);
        } else {
            compute_nmu(a, b, edge, time);
        }
    }

    std::vector<std::vector<Homology>> &compute_homologies(int prime = 2) {
        homologies.resize(complex.size() + 1);
        for (int i = 0; i != vertices.size(); ++i) {
            if (vertices[i].pair == unpaired) {
                homologies[0].emplace_back(vertices[i].time);
            }
        }
        for (size_t i = 1; i < complex.size(); ++i) {
            reduce(i, prime);
        }
        for (size_t i = 1; i < complex.size(); ++i) {
            check_for_immortal(i);
        }
        for (auto &matrix : complex) {
            for (Simplex &simplex : matrix) {
                int res = simplex.last_nmu();
                if (res > 0) {
                    compute_link_homologies(simplex, graph, homologies, prime);
                }
            }
        }
        return homologies;
    }
};

#include "find_cliques.h"
