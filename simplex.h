#include <map>
#include <vector>
#include <set>

#include "homology.h"
#include "edge_inflater.h"
#include "zp.h"

#ifndef MULTIGRAPHS_SIMPLEX_H
#define MULTIGRAPHS_SIMPLEX_H

enum {
    unpaired = -1
};

class Vertex {
public:
    int pair;
    double time;
    Vertex(int pair, double time): pair(pair), time(time) {}
};

class Simplex {
private:
    int smaller_pair;
    int bigger_pair;
    double time;
    std::vector<int> subsimpleces;
    // это для nmu, хочется удалить
    std::vector<int> parents;
    // for each dim
    std::vector<double> summands = {1};
    std::vector<int> link_vertices;
    std::vector<double> link_timeline;

public:
    std::vector<double> timeline;
    std::vector<unsigned> nmu_vals;

    Simplex() {}
    explicit Simplex(double time):
            time(time),bigger_pair(unpaired), smaller_pair(unpaired) {}
    // edge constructor without known pair
    Simplex(int a, int b, double time):
            time(time), smaller_pair(unpaired), bigger_pair(unpaired), subsimpleces(2) {
        subsimpleces[0] = a;
        subsimpleces[1] = b;
    }
    // edge constructor, with a known vertex pair
    Simplex(int a, int b, double time, int curr_ind):
            time(time), smaller_pair(curr_ind), bigger_pair(unpaired), subsimpleces(2) {
        subsimpleces[0] = a;
        subsimpleces[1] = b;
    }

    bool no_pair() {
        return smaller_pair == unpaired && bigger_pair == unpaired;
    }
    double birth() {
        return time;
    }
    int size() {
        return subsimpleces.size();
    }
    std::vector<int> &link() {
        return link_vertices;
    }
    std::vector<double> &link_times() {
        return link_timeline;
    }

    void push_sub(size_t sub_ind) {
        subsimpleces.push_back(sub_ind);
    }

    void push_parent(size_t sub_ind) {
        parents.push_back(sub_ind);
    }
    void push_vertex(int vertex, double time) {
        link_vertices.push_back(vertex);
        link_timeline.push_back(time);
    }
    void compute_edge_nmu(std::vector<std::vector<Simplex>> &columns,
                          std::queue<int> &order, double change_time, int curr) {
        if (summands[0] != curr) {
            summands[0] = curr;
            if (!timeline.empty() && timeline.back() == change_time) {
                nmu_vals.back() = summands[0] - 1;
            } else {
                timeline.push_back(change_time);
                nmu_vals.push_back(summands[0] - 1);
            }
        }
        for (size_t i = 0; i != parents.size(); ++i) {
            order.push(parents[i]);
        }
    }
    void compute_nmu(std::vector<std::vector<Simplex>> &columns,
                     std::queue<int> &order, double change_time) {
        int sz = subsimpleces.size();
        summands[0] = 0;
        summands.resize(sz - 1, 0);
        summands.back() = 1;
        for (int i = 0; i != sz; ++i) {
            std::vector<double> &nmu_curr = columns[subsimpleces.size() - 3][subsimpleces[i]].summands;
            for (int j = 0; j != nmu_curr.size(); ++j) {
                summands[j] += nmu_curr[j] / (sz - j - 2);
            }
            summands.back() *= std::pow(nmu_curr.back(), 1 / ((double) sz - 2));
        }
        int new_val = sum_nmu();
        if (new_val != 0 && (nmu_vals.empty() || new_val != nmu_vals.back())) {
            if (!timeline.empty() && timeline.back() == change_time) {
                nmu_vals.back() = new_val;
            } else {
                timeline.push_back(change_time);
                nmu_vals.push_back(sum_nmu());
            }
            for (size_t i = 0; i != parents.size(); ++i) {
                order.push(parents[i]);
            }
        }
    }
    int sum_nmu() {
        if (subsimpleces.size() == 2) {
            return (int) summands[0] - 1;
        }
        int ind = (subsimpleces.size() % 2 ? 1 : -1);
        int res = ((int) subsimpleces.size() - 1) * ind;
        for (double to_add : summands) {
            ind *= -1;
            res += (int) std::round(to_add) * ind;
        }
        return res;
    }
    int last_nmu() {
        if (nmu_vals.empty()) {
            return 0;
        } else {
            return nmu_vals.back();
        }
    }

private:
    std::map<int, Zp> coefs;
    void reduce(Simplex &other, Zp red_coef) {
        for (auto j = other.coefs.begin(); j != other.coefs.end(); ++j) {
            coefs[j->first] = coefs[j->first] + red_coef * other.coefs[j->first];
            if (coefs[j->first] == 0) {
                coefs.erase(j->first);
            }
        }
    }

public:
    // coefs should be initialized in reduce
    Zp lowest() {
        return coefs.rbegin()->second;
    }
    Homology reduce(std::vector<Simplex> &columns,
                std::vector<Simplex> &rows,
                size_t curr, int prime) {
        int curr_coef = subsimpleces.size() % 2 ? 1 : -1;
        for (size_t i = 0; i != subsimpleces.size(); ++i) {
            coefs[subsimpleces[i]] = Zp(prime, curr_coef);
            curr_coef *= -1;
        }
        size_t ind = subsimpleces.back();
        while (rows[ind].bigger_pair != unpaired) {
            Simplex &candidate = columns[rows[ind].bigger_pair];
            this->reduce(candidate, coefs[ind] / candidate.lowest());
            if (coefs.empty()) {
                break;
            }
            ind = coefs.rbegin()->first;
        }
        if (rows[ind].bigger_pair == unpaired) {
            rows[ind].bigger_pair = curr;
            columns[curr].smaller_pair = ind;
            if (rows[ind].time < columns[curr].time) {
                return Homology(rows[ind].time, columns[curr].time);
            }
        }
        return Homology();
    }
};



#endif //MULTIGRAPHS_SIMPLEX_H
