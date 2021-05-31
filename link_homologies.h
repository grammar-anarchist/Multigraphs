#include <algorithm>
#include <unordered_set>

#include "homology.h"
#include "testers.h"
#include "edge_inflater.h"
#include "simplex.h"
#include "hash.h"

#ifndef MULTIGRAPHS_LINK_HOMOLOGIES_H
#define MULTIGRAPHS_LINK_HOMOLOGIES_H

class CliqueLinkComplex {
    Infestation graph;
    std::vector<std::vector<Simplex>> complex;
    std::vector<Vertex> vertices;
    std::vector<std::unordered_map<std::set<int>, int>> search;
    // number of vertices in original matrix
    size_t correction;

    void find_all_cliques(
            std::set<int> &curr,
            std::vector<int> &candidates,
            std::unordered_set<int> &out,
            double time,
            int a,
            int b,
            size_t prev_ind
    )
    {
        if (candidates.empty() && out.empty()) {
            return;
        }
        for (int i = 0; i != candidates.size(); ++i) {
            size_t dim = curr.size() - 1;
            size_t dest;
            if (out.find(candidates[i]) != out.end()) {
                curr.insert(candidates[i]);
                dest = search[dim][curr];
                curr.erase(candidates[i]);
                complex[dim][dest].push_sub(prev_ind);
                continue;
            }
            curr.insert(candidates[i]);
            // std::cout << curr << "\n";
            if (dim == complex.size()) {
                complex.emplace_back();
                search.emplace_back();
            }
            unsigned curr_ind = complex[dim].size();
            complex[dim].emplace_back(time);
            Simplex &curr_simplex = complex[dim].back();

            search[dim][curr] = curr_ind;

            curr.erase(a);
            dest = search[dim - 1][curr];
            curr_simplex.push_sub(dest);
            curr.insert(a);

            curr.erase(b);
            dest = search[dim - 1][curr];
            curr_simplex.push_sub(dest);
            curr.insert(b);

            curr_simplex.push_sub(prev_ind);

            std::vector<int> new_candidates;

            for (size_t j = 0; j != candidates.size(); ++j) {
                int c, d;
                c = std::min(candidates[i], candidates[j]);
                d = std::max(candidates[i], candidates[j]);
                if (graph.edge_count(c, d))
                    new_candidates.push_back(candidates[j]);
            }
            std::unordered_set<int> new_out;
            for (auto j : out) {
                int c, d;
                c = std::min(candidates[i], j);
                d = std::max(candidates[i], j);
                if (graph.edge_count(candidates[i], j))
                    new_out.insert(j);
            }
            find_all_cliques(curr, new_candidates, new_out, time, a, b, curr_ind);
            curr.erase(candidates[i]);
            out.insert(candidates[i]);
        }
    }

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
    }

    void reduce(size_t ind, std::vector<Homology> homologies, Simplex &curr, int prime) {
        for (size_t i = 0; i != complex[ind].size(); ++i) {
            Homology res = complex[ind][i].reduce(complex[ind], complex[ind - 1], i, prime);
            if (res.is_homology()) {
                add_homology(res, homologies, curr);
            }
        }
    }

    void check_for_immortal(size_t ind, std::vector<Homology> &homologies, Simplex &curr) {
        std::vector<Simplex> &columns = complex[ind];
        for (size_t i = 0; i != columns.size(); ++i) {
            if (columns[i].no_pair()) {
                double birth = columns[i].birth();
                auto res = Homology(birth);
                add_homology(res, homologies, curr);
            }
        }
    }

    void add_homology(Homology res, std::vector<Homology> &homologies, Simplex &curr) {
        for (int j = 0; j != curr.timeline.size() - 1; ++j) {
            if (res.birth <= curr.timeline[j] && res.death > curr.timeline[j]) {
                homologies.emplace_back(curr.timeline[j],
                                        std::min(res.death, curr.timeline[j + 1]),
                                        curr.nmu_vals[j]);
            }
        }
        if (res.birth <= curr.timeline.back()) {
            homologies.emplace_back(curr.timeline.back(),
                                    res.death,
                                    curr.nmu_vals.back());
        }
    }

public:
    explicit CliqueLinkComplex(size_t correction):
            complex(1), correction(correction), search(1) {}
    void add_vertex(double time) {
        vertices.emplace_back(unpaired, time);
    }
    void add_edge(int a, int b, double time) {
        if (a == b) {
            return;
        }
        if (a > b) {
            std::swap(a, b);
        }
        // a < b
        std::set<int> edge = {a, b};
        graph(a, b);
        int curr_ind = complex[0].size();
        if (vertices[b].pair == unpaired) {
            vertices[b].pair = curr_ind;
            complex[0].emplace_back(a, b, time, b);
        }
        else {
            complex[0].emplace_back(a, b, time);
        }
        search[0][edge] = curr_ind;
        change_topology(curr_ind, a, b, time, edge);
    }

    std::vector<std::vector<Homology>> &compute_homologies(
            std::vector<std::vector<Homology>> &homologies,
            Simplex &curr,
            int prime = 2) {
        for (int i = 0; i != vertices.size(); ++i) {
            if (vertices[i].pair == unpaired) {
                auto res = Homology(vertices[i].time);
                add_homology(res, homologies[correction], curr);
            } else {
                double vertex_time = vertices[i].time;
                double edge_time = complex[0][vertices[i].pair].birth();
                if (edge_time > vertex_time) {
                    auto res = Homology(vertex_time, edge_time);
                    add_homology(res, homologies[correction], curr);
                }
            }
        }
        for (size_t i = 1; i < complex.size(); ++i) {
            reduce(i, homologies[i + correction], curr, prime);
        }
        for (size_t i = 1; i < complex.size(); ++i) {
            check_for_immortal(i, homologies[i + correction], curr);
        }
        return homologies;
    }
};

struct Edge_temp {
public:
    double link_edge_time;
    int a;
    int b;

    bool operator==(const Edge_temp &other) const {
        return this->link_edge_time == other.link_edge_time;
    }
    bool operator!=(const Edge_temp &other) const {
        return this->link_edge_time != other.link_edge_time;
    }
    bool operator<(const Edge_temp &other) const {
        return this->link_edge_time < other.link_edge_time;
    }
    bool operator<=(const Edge_temp &other) const {
        return this->link_edge_time <= other.link_edge_time;
    }
};

void compute_link_homologies(Simplex& curr,
                             Infestation &original_graph,
                             std::vector<std::vector<Homology>> &homologies,
                             int prime) {
    CliqueLinkComplex clc(curr.size());
    std::vector<int> &link_vertices = curr.link();
    std::vector<double> &link_timeline = curr.link_times();
    for (int i = 0; i != link_vertices.size(); ++i) {
        clc.add_vertex(link_timeline[i]);
    }
    std::vector<Edge_temp> edges;
    for (int i = 0; i != link_vertices.size(); ++i) {
        for (int j = i + 1; j != link_vertices.size(); ++j) {
            double edge_time = original_graph.edge_birth(link_vertices[i], link_vertices[j]);
            if (original_graph.edge_count(i, j)) {
                double link_edge_time = std::max(
                        {link_timeline[i],
                            link_timeline[j],
                            edge_time});
                edges.push_back(Edge_temp({link_edge_time, i, j}));
            }

        }
    }
    std::sort(edges.begin(), edges.end());
    for (auto &edge : edges) {
        clc.add_edge(edge.a, edge.b, edge.link_edge_time);
    }
    clc.compute_homologies(homologies, curr, prime);
}

#endif //MULTIGRAPHS_LINK_HOMOLOGIES_H
