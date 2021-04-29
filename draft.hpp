#include <algorithm>
#include <future>
#include <iostream>
#include <limits>
#include <utility>
#include <vector>
#include <unordered_set>
#include <unordered_map>

#include <boost/container_hash/hash.hpp>
#include <cmath>

class Homology {
public:
    Homology(unsigned smaller, unsigned greater, double birth, double death):
            smaller(smaller), greater(greater), birth(birth), death(death) {}

    Homology(unsigned smaller, double birth):
            smaller(smaller), birth(birth), death(INFINITY) {}

    unsigned smaller;
    unsigned greater;
    double birth;
    double death;
};

template<typename Elem>
std::ostream &operator<<(std::ostream &os, std::vector<Elem> curr) {
    for (Elem i : curr) {
        os << i << " ";
    }
    return os;
}

template<typename Elem>
std::ostream &operator<<(std::ostream &os, std::set<Elem> curr) {
    for (auto it = curr.begin(); it != curr.end(); ++it) {
        os << *it << " ";
    }
    return os;
}

namespace std {
    template<>
    struct hash<std::set<unsigned>> {
    public:
        size_t operator()(const std::set<unsigned> &s) const {
            size_t res = 42;
            boost::hash_combine<std::set<unsigned>>(res, s);
            return res;
        }
    };
}

class Infestation {
    // first index is less than second
    std::unordered_map<unsigned, std::unordered_map<unsigned, unsigned>> func;

public:
    // too unequal indeces
    unsigned &operator()(unsigned a, unsigned b) {
        if (a < b)
            return func[a][b];
        return func[b][a];
    }

    unsigned has_edge(unsigned a, unsigned b) {
        if (a == b)
            return false;
        if (a < b) {
            if (func.find(a) != func.end() && func[a].find(b) != func[a].end()) {
                return func[a][b];
            }
        }
        if (func.find(b) != func.end() && func[b].find(a) != func[b].end()) {
            return func[b][a];
        }
        return false;
    }
};
enum {
    unpaired = -1
};

class Vertex {
public:
    ssize_t pair;
    double time;
    Vertex(ssize_t pair, double time): pair(pair), time(time) {}
};

// actually Z2
class Zp {
    int prime;
public:
    Zp(int prime): prime(prime) {}
    int operator()(int a) {
        return 1;
    }
    int del(int a, int b) {
        return a;
    }
    int min(int a, int b) {
        if (a == 0 && b == 1)
            return b;
        return a - b;
    }
    int mul(int a, int b) {
        return a * b;
    }
};

class CliqueLinkComplex;

class Simplex {
public:
    ssize_t smaller_pair;
    ssize_t bigger_pair;
    double time;
    std::vector<unsigned> subsimpleces;
    std::set<unsigned> self;

    Simplex(std::set<unsigned> self, double time):
            time(time),bigger_pair(unpaired), smaller_pair(unpaired), self(std::move(self)){}
    // constructor from the vector of vertices, time, last edge
    // and number which equals zero if a pair should be created
    Simplex(unsigned a, unsigned b, double time):
            time(time), smaller_pair(unpaired), bigger_pair(unpaired), subsimpleces(2) {
        subsimpleces[0] = a;
        subsimpleces[1] = b;
        self.insert(a);
        self.insert(b);
    }
    Simplex(unsigned a, unsigned b, double time, ssize_t curr_ind):
            time(time), smaller_pair(curr_ind), bigger_pair(unpaired), subsimpleces(2) {
        subsimpleces[0] = a;
        subsimpleces[1] = b;
        self.insert(a);
        self.insert(b);
    }
    Simplex() {}

    void push_sub(size_t sub_ind) {
        subsimpleces.push_back(sub_ind);
    }

    std::vector<unsigned> &edge() {
        return subsimpleces;
    }

private:
    void reduce(Simplex &other, int red_coef, Zp zp) {
        for (auto j = other.coefs.begin(); j != other.coefs.end(); ++j) {
            coefs[j->first] = zp.min(coefs[j->first], zp.mul(red_coef, other.coefs[j->first]));
            if (coefs[j->first] == 0) {
                coefs.erase(j->first);
            }
        }
    }

public:
    // it is not necessary to initialize coeffs for ALL simpleces

    std::map<unsigned, int> coefs;
    //coefs initialized
    int lowest() {
        return coefs.rbegin()->second;
    }
    void reduce(std::vector<Simplex> &columns,
                std::vector<Simplex> &rows,
                std::vector<Homology> &curr_homologies,
                size_t curr,
                Zp zp) {
        int curr_coef = subsimpleces.size() % 2 ? 1 : -1;
        for (size_t i = 0; i != subsimpleces.size(); ++i) {
            coefs[subsimpleces[i]] = zp(curr_coef);
            curr_coef *= -1;
        }
        size_t ind = subsimpleces.back();
        while (rows[ind].bigger_pair != unpaired) {
            Simplex &candidate = columns[rows[ind].bigger_pair];
            this->reduce(candidate, zp.del(coefs[ind], candidate.lowest()), zp);
            if (coefs.empty()) {
                break;
            }
            ind = coefs.rbegin()->first;
        }
        if (rows[ind].bigger_pair == unpaired) {
            rows[ind].bigger_pair = curr;
            columns[curr].smaller_pair = ind;
            if (rows[ind].time < columns[curr].time) {
                curr_homologies.emplace_back(ind, curr, rows[ind].time, columns[curr].time);
            }
        }
    }

    std::vector<double> nmu = {1};
    void compute_nmu(Infestation &mu, std::vector<Simplex> &columns) {
        int sz = subsimpleces.size();
        nmu[0] = 0;
        nmu.resize(sz - 1, 0);
        nmu.back() = 1;
        for (int i = 0; i != sz; ++i) {
            std::vector<double> &nmu_curr = columns[subsimpleces[i]].nmu;
            for (int j = 0; j != nmu_curr.size(); ++j) {
                nmu[j] += nmu_curr[j] / (sz - j - 2);
            }
            nmu.back() *= std::pow(nmu_curr.back(), 1 / ((double) sz - 2));
        }
    }
    // only for edge
    void compute_nmu(Infestation &mu) {
        nmu[0] = mu(subsimpleces[0], subsimpleces[1]);
    }
    int sum_nmu() {
        int ind = (subsimpleces.size() % 2 ? 1 : -1);
        int res = ((int) subsimpleces.size() - 1) * ind;
        for (double to_add : nmu) {
            ind *= -1;
            res += (int) std::round(to_add) * ind;
        }
        return res;
    }
};

void compute_link_homologies(Simplex& curr,
                             Infestation &mu,
                             const std::vector<Vertex> &vertices,
                             std::vector<std::vector<Homology>> &homologies);

class CliqueMultiComplex {
    Infestation mu;
    std::vector<std::vector<Simplex>> complex;
    std::vector<Vertex> vertices;
    std::vector<std::vector<Homology>> homologies;
    std::unordered_map<std::set<unsigned>, unsigned> search;

    // what if we trace parents instead of children
    // can be optimized easily as Bron-Kerbosch algorithm
    void find_all_cliques(
            std::set<unsigned> &curr,
            std::vector<unsigned> &candidates,
            std::unordered_set<unsigned> &out,
            double time,
            unsigned a,
            unsigned b,
            size_t prev_ind
    )
    {
        if (candidates.empty() && out.empty()) {
            // std::cout << curr;
            return;
        }
        for (ssize_t i = 0; i != candidates.size(); ++i) {
            size_t sec;
            if (out.find(candidates[i]) != out.end()) {
                sec = curr.size() - 1;
                curr.insert(candidates[i]);
                size_t dest = search[curr];
                curr.erase(candidates[i]);
                complex[sec][dest].push_sub(prev_ind);
                continue;
            }
            curr.insert(candidates[i]);
            // std::cout << curr;
            sec = curr.size() - 2;
            if (sec == complex.size()) {
                complex.emplace_back();
            }
            size_t curr_ind = complex[sec].size();
            complex[sec].emplace_back(curr, time);
            Simplex &curr_simplex = complex[sec].back();

            search[curr] = curr_ind;

            curr.erase(a);
            curr_simplex.push_sub(search[curr]);
            curr.insert(a);

            curr.erase(b);
            curr_simplex.push_sub(search[curr]);
            curr.insert(b);

            complex[sec].back().push_sub(prev_ind);

            std::vector<unsigned> new_candidates;
            for (size_t j = 0; j != candidates.size(); ++j) {
                if (mu.has_edge(candidates[i], candidates[j]))
                    new_candidates.push_back(candidates[j]);
            }
            std::unordered_set<unsigned> new_out;
            for (auto j : out) {
                if (mu.has_edge(candidates[i], j))
                    new_out.insert(j);
            }
            find_all_cliques(curr, new_candidates, new_out, time, a, b, curr_ind);
            curr.erase(candidates[i]);
            out.insert(candidates[i]);
        }
    }

    void change_topology(size_t new_edge, unsigned a, unsigned b, double time) {
        std::vector<unsigned> candidates;
        std::unordered_set<unsigned> out;
        for (ssize_t i = vertices.size() - 1; i != -1; --i) {
            if (i != a && i != b && mu.has_edge(a, i) && mu.has_edge(b, i))
                candidates.push_back(i);
        }
        std::set<unsigned> curr = {a, b};
        find_all_cliques(curr, candidates, out, time, a, b, new_edge);
    }

    static void reduce(std::vector<Simplex> &columns,
                       std::vector<Simplex> &rows,
                       std::vector<Homology> &curr_homologies,
                       Zp &zp) {
        for (size_t i = 0; i != columns.size(); ++i) {
            columns[i].reduce(columns, rows, curr_homologies, i, zp);
        }
    }

    static void check_for_immortal(std::vector<Simplex> &columns,
                                   std::vector<Homology> &curr_homologies) {
        for (size_t i = 0; i != columns.size(); ++i) {
            if (columns[i].smaller_pair == unpaired && columns[i].bigger_pair == unpaired) {
                curr_homologies.emplace_back(i, columns[i].time);
            }
        }
    }

public:
    CliqueMultiComplex(): complex(1), homologies(1) {}
    void add_vertex(double time) {
        vertices.emplace_back(unpaired, time);
    }
    void add_edge(unsigned a, unsigned b, double time) {
        // std::cout << a << " " << b << "\n";
        if (a == b) {
            return;
        }
        if (a > b) {
            std::swap(a, b);
        }
        if (++mu(a, b) == 1) {
            size_t curr_ind = complex[0].size();
            if (vertices[b].pair == unpaired) {
                vertices[b].pair = curr_ind;
                complex[0].emplace_back(a, b, time, b);
                if (time > vertices[b].time) {
                    homologies[0].emplace_back(b, curr_ind, vertices[b].time, time);
                }
            }
            else {
                complex[0].emplace_back(a, b, time);
            }
            search[complex[0].back().self] = curr_ind;
            change_topology(curr_ind, a, b, time);
        }
    }

    // add_parallelism
    std::vector<std::vector<Homology>> &compute_homologies(int prime = 2) {
        class Zp zp(prime);
        homologies.resize(complex.size() + 1);
        for (unsigned i = 0; i != vertices.size(); ++i) {
            if (vertices[i].pair == unpaired) {
                homologies[0].emplace_back(i, vertices[i].time);
            }
        }
        /*printf("%lu dims\n", complex.size());
        for (size_t i = 0; i < complex.size(); ++i) {
            printf("%lu ", complex[i].size());
        }
        printf("\n");*/
        for (size_t i = 1; i < complex.size(); ++i) {
            // printf("reducing %lu\n", i);
            reduce(complex[i], complex[i - 1], homologies[i], zp);
        }
        for (size_t i = 0; i < complex.size(); ++i) {
            check_for_immortal(complex[i], homologies[i + 1]);
        }
        for (Simplex &simplex : complex[0]) {
            simplex.compute_nmu(mu);
        }
        for (size_t i = 1; i < complex.size(); ++i) {
            for (Simplex &simplex : complex[i]) {
                simplex.compute_nmu(mu, complex[i - 1]);
            }
        }
        int link = 0;
        for (auto &matrix : complex) {
            for (Simplex &simplex : matrix) {
                int res = simplex.sum_nmu();
                if (res > 0) {
                    /*std::cout << simplex.self << "\t";
                    std::cout << res << "\n";*/
                    ++link;
                    compute_link_homologies(simplex, mu, vertices, homologies);
                }
            }
        }
        printf("link: %d\n", link);
        return homologies;
    }
};

class CliqueLinkComplex {
    Infestation mu;
    std::vector<std::vector<Simplex>> complex;
    std::vector<Vertex> vertices;
    std::unordered_map<std::set<unsigned>, unsigned> search;
    std::vector<std::vector<Homology>> &homologies;
    // number of vertices in original matrix
    size_t correction;

    void find_all_cliques(
            std::set<unsigned> &curr,
            std::vector<unsigned> &candidates,
            std::unordered_set<unsigned> &out,
            double time,
            unsigned a,
            unsigned b,
            size_t prev_ind
    )
    {
        if (candidates.empty() && out.empty()) {
            // std::cout << curr;
            return;
        }
        for (ssize_t i = 0; i != candidates.size(); ++i) {
            size_t sec;
            if (out.find(candidates[i]) != out.end()) {
                sec = curr.size() - 1;
                curr.insert(candidates[i]);
                size_t dest = search[curr];
                curr.erase(candidates[i]);
                complex[sec][dest].push_sub(prev_ind);
                continue;
            }
            curr.insert(candidates[i]);
            // std::cout << curr;
            sec = curr.size() - 2;
            if (sec == complex.size()) {
                complex.emplace_back();
            }
            size_t curr_ind = complex[sec].size();
            complex[sec].emplace_back(curr, time);
            Simplex &curr_simplex = complex[sec].back();

            search[curr] = curr_ind;

            curr.erase(a);
            curr_simplex.push_sub(search[curr]);
            curr.insert(a);

            curr.erase(b);
            curr_simplex.push_sub(search[curr]);
            curr.insert(b);

            complex[sec].back().push_sub(prev_ind);

            std::vector<unsigned> new_candidates;
            for (size_t j = 0; j != candidates.size(); ++j) {
                if (mu.has_edge(candidates[i], candidates[j]))
                    new_candidates.push_back(candidates[j]);
            }
            std::unordered_set<unsigned> new_out;
            for (auto j : out) {
                if (mu.has_edge(candidates[i], j))
                    new_out.insert(j);
            }
            find_all_cliques(curr, new_candidates, new_out, time, a, b, curr_ind);
            curr.erase(candidates[i]);
            out.insert(candidates[i]);
        }
    }

    void change_topology(size_t new_edge, unsigned a, unsigned b, double time) {
        std::vector<unsigned> candidates;
        std::unordered_set<unsigned> out;
        for (ssize_t i = vertices.size() - 1; i != -1; --i) {
            if (i != a && i != b && mu.has_edge(a, i) && mu.has_edge(b, i))
                candidates.push_back(i);
        }
        std::set<unsigned> curr = {a, b};
        find_all_cliques(curr, candidates, out, time, a, b, new_edge);
    }

    static void reduce(std::vector<Simplex> &columns,
                       std::vector<Simplex> &rows,
                       std::vector<Homology> &curr_homologies,
                       Zp &zp) {
        for (size_t i = 0; i != columns.size(); ++i) {
            columns[i].reduce(columns, rows, curr_homologies, i, zp);
        }
    }

    static void check_for_immortal(std::vector<Simplex> &columns,
                                   std::vector<Homology> &curr_homologies) {
        for (size_t i = 0; i != columns.size(); ++i) {
            if (columns[i].smaller_pair == unpaired && columns[i].bigger_pair == unpaired) {
                curr_homologies.emplace_back(i, columns[i].time);
            }
        }
    }

public:
    CliqueLinkComplex(std::vector<std::vector<Homology>> &outer, size_t correction):
        complex(1), homologies(outer), correction(correction) {}
    void add_bundle(double time) {
        add_vertex(time);
        unsigned curr = vertices.size() - 1;
        for (size_t i = 0; i != curr; ++i) {
            add_edge(curr, i, time);
        }
    }
    void add_vertex(double time) {
        vertices.emplace_back(unpaired, time);
    }
    void add_edge(unsigned a, unsigned b, double time) {
        // std::cout << a << " " << b << "\n";
        if (a > b) {
            std::swap(a, b);
        }
        size_t curr_ind = complex[0].size();
        if (vertices[b].pair == unpaired) {
            vertices[b].pair = curr_ind;
            complex[0].emplace_back(a, b, time, b);
            if (time > vertices[b].time) {
                homologies[0 + correction].emplace_back(b, curr_ind, vertices[b].time, time);
            }
        }
        else {
            complex[0].emplace_back(a, b, time);
        }
        search[complex[0].back().self] = curr_ind;
        change_topology(curr_ind, a, b, time);
    }

    std::vector<std::vector<Homology>> &compute_homologies(int prime = 2) {
        class Zp zp(prime);
        for (unsigned i = 0; i != vertices.size(); ++i) {
            if (vertices[i].pair == unpaired) {
                homologies[correction].emplace_back(i, vertices[i].time);
            }
        }
        /*printf("%lu dims\n", complex.size());
        for (size_t i = 0; i < complex.size(); ++i) {
            printf("%lu ", complex[i].size());
        }
        printf("\n");*/
        for (size_t i = 1; i < complex.size(); ++i) {
            reduce(complex[i], complex[i - 1], homologies[i + correction], zp);
        }
        for (size_t i = 0; i < complex.size(); ++i) {
            check_for_immortal(complex[i], homologies[i + 1 + correction]);
        }
        return homologies;
    }
};

void compute_link_homologies(Simplex& curr,
                             Infestation &mu,
                             const std::vector<Vertex> &vertices,
                             std::vector<std::vector<Homology>> &homologies) {
    CliqueLinkComplex clc(homologies, curr.self.size());
    int added = 0;
    for (size_t i = 0; i != vertices.size(); ++i) {
        double time_vertex = -1;
        for (unsigned v : curr.self) {
            if (!mu.has_edge(v, i)) {
                time_vertex = -1;
                break;
            } else {
                // здесь должно быть время для данного ребра
                time_vertex = -2;
            }
        }
        if (time_vertex != -1) {
            clc.add_bundle(time_vertex);
            ++added;
        }
    }
    if (added == 0) {
        // pending correction
        homologies[curr.self.size() - 1].emplace_back(0, INFINITY);
    }
    // clc.compute_homologies();
}