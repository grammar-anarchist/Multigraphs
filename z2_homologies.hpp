#include <algorithm>
#include <future>
#include <iostream>
#include <limits>
#include <vector>

struct Simplex {
    std::vector<int> vertices;
    double time;
    int dim;

    bool operator==(const Simplex &other) const {
        if (vertices == other.vertices) {
            return 1;
        }
        return 0;
    }
    bool operator!=(const Simplex &other) const {
        return !((*this) == other);
    }

    friend std::ostream& operator<<(std::ostream& os, const Simplex &curr) {
        os << "[ ";
        for (size_t i = 0; i != curr.vertices.size(); ++i) {
            os << curr.vertices[i] << " ";
        }
        os << "] – ";
        os << curr.time;
        return os;
    }
};

struct Homology {
    size_t ind_birth;
    size_t ind_death;
    double birth;
    double death;

    Homology(size_t ind_birth, size_t ind_death, double birth, double death):
        ind_birth(ind_birth), ind_death(ind_death), birth(birth), death(death) {}
};

class Column {
public:
    std::vector<int> column;
    size_t lowest = 0;

    Column(
        const Simplex &curr,
        const std::vector<size_t> &sub_inds, 
        const std::vector<Simplex> &simpleces) 
    {
        for (size_t i = 0, gathered = 0; i != sub_inds.size() && gathered != curr.dim + 1; ++i) {
            if (std::includes(
                    curr.vertices.begin(), curr.vertices.end(),
                    simpleces[sub_inds[i]].vertices.begin(), 
                    simpleces[sub_inds[i]].vertices.end())) {
                column.push_back(1);
                ++gathered;
                lowest = i + 1;
            } else {
                column.push_back(0);
            }
        }
    }

    void reduce(Column &reducer) {
        size_t new_lowest = 0;
        for (size_t i = 0; i != lowest; ++i) {
            if (reducer[i] == 1 && (*this)[i] == 1) {
                (*this)[i] = 0;
            } else if (reducer[i] == 1 && (*this)[i] == 0) {
                (*this)[i] = 1;
                new_lowest = i + 1;
            } else if ((*this)[i] == 1) {
                new_lowest = i + 1;
            }
        }
        lowest = new_lowest;
    }

    int& operator[](size_t i) {
        return column[i];
    }

    void print() {
        for (size_t i = 0; i != column.size(); ++ i) {
            std::cout << column[i] << " ";
        }
    }
};

class Matrix {
public:
    std::vector<size_t> sub_inds;
    std::vector<size_t> higher_inds;
    std::vector<Column> matrix;

    void push_smaller(size_t simplex_ind) {
        sub_inds.push_back(simplex_ind);
    }

    void push_bigger(size_t simplex_ind, 
        const Simplex &curr,
        const std::vector<Simplex> &simpleces) 
    {
        higher_inds.push_back(simplex_ind);
        matrix.emplace_back(curr, sub_inds, simpleces);
    }

    void reduce(
            std::vector<Homology> &curr_dim_homologies, 
            const std::vector<Simplex> &simpleces,
            std::vector<int> &paired)
    {
        for (size_t curr_reduced = sub_inds.size(); curr_reduced != 0; --curr_reduced) {
            /*std::cout << "curr_reduced" << curr_reduced << "\n";
            for (size_t i = 0; i != matrix.size(); ++ i) {
                matrix[i].print();
                std::cout << "\n";
            }*/
            size_t ind = 0;
            while (ind < matrix.size() && matrix[ind].lowest != curr_reduced) {
                ++ind;
            }
            if (ind == matrix.size()) {
                continue;
            }
            Column &reducer = matrix[ind];
            size_t bearer = sub_inds[curr_reduced - 1];
            size_t murderer = higher_inds[ind];
            double birth = simpleces[bearer].time;
            double death = simpleces[murderer].time;
            if (birth != death) {
                curr_dim_homologies.push_back(Homology(
                    bearer,
                    murderer,
                    simpleces[bearer].time, 
                    simpleces[murderer].time
                ));
            }
            paired[sub_inds[curr_reduced - 1]] = 1;
            paired[higher_inds[ind]] = 1;
            ++ind;
            while (ind < matrix.size()) {
                if (matrix[ind].lowest == curr_reduced) {
                    matrix[ind].reduce(reducer);
                }
                ++ind;
            }
        }
    }

    void reduce_link(
            std::vector<Homology> &curr_dim_homologies, 
            const std::vector<Simplex> &simpleces,
            std::vector<int> &paired,
            std::vector<int> &ignored)
    {
        for (size_t curr_reduced = sub_inds.size(); curr_reduced != 0; --curr_reduced) {
            if (ignored[sub_inds[curr_reduced - 1]]) {
                continue;
            }
            /*std::cout << "curr_reduced" << curr_reduced << "\n";
            for (size_t i = 0; i != matrix.size(); ++ i) {
                matrix[i].print();
                std::cout << "\n";
            }*/
            size_t ind = 0;
            while (ind < matrix.size() && 
                (matrix[ind].lowest < curr_reduced || 
                matrix[ind][curr_reduced - 1] == 0 || 
                ignored[higher_inds[ind]])) 
            {
                ++ind;
            }
            if (ind == matrix.size()) {
                continue;
            }
            Column &reducer = matrix[ind];
            size_t bearer = sub_inds[curr_reduced - 1];
            size_t murderer = higher_inds[ind];
            double birth = simpleces[bearer].time;
            double death = simpleces[murderer].time;
            if (birth != death) {
                curr_dim_homologies.push_back(Homology(
                    bearer,
                    murderer,
                    simpleces[bearer].time, 
                    simpleces[murderer].time
                ));
            }
            paired[sub_inds[curr_reduced - 1]] = 1;
            paired[higher_inds[ind]] = 1;
            ignored[higher_inds[ind]] = 1;
            ++ind;
            while (ind < matrix.size()) {
                if (matrix[ind].lowest == curr_reduced && !ignored[higher_inds[ind]]) {
                    matrix[ind].reduce(reducer);
                }
                ++ind;
            }
        }
    }
};

class Filtration {
private:
    std::vector<Matrix> matrices;
    std::vector<Simplex> simpleces;
    size_t sz = 0;
    int dim = -1;

public:
    //push Simplex unsafely:
    // all younger Simpleces should be present
    // otherwise Filtration will be incorrect
    void push(const Simplex &simplex) {
        simpleces.push_back(simplex);
        dim = std::max(dim, simplex.dim);
        if (simplex.dim >= matrices.size()) {
            matrices.emplace_back();
        }
        if (simplex.dim > 0) {
            matrices[simplex.dim - 1].push_bigger(sz, simplex, simpleces);
        }
        matrices[simplex.dim].push_smaller(sz);
        ++sz;
    }

    // finds persistent homologies over Z2
    std::vector<std::vector<Homology>> persistent_homologies() {
        std::vector<std::vector<Homology>> ans;
        if (dim == -1) {
            return ans;
        }
        ans.resize(dim + 1);
        std::vector<int> paired(sz, 0);
        std::vector<std::future<void>> futures;
        int future_num = 8;
        for (int i = 0; i != future_num; ++i) {
            futures.push_back(std::async(
                [i, future_num, this, &ans, &paired] {
                for (size_t k = i; k < this->matrices.size() - 1; k += future_num) {
                    this->matrices[k].reduce(ans[k], this->simpleces, paired);
                }
            }
            ));
        }
        for (int i = 0; i != future_num; ++i) {
            futures[i].get();
        }
        double inf = std::numeric_limits<double>::infinity();
        for (size_t i = 0; i != sz; ++i) {
            if (!paired[i]) {
                ans[simpleces[i].dim].emplace_back(i, i, simpleces[i].time, inf);
            }
        }
        return ans;
    }

    // finds persistent homologies of link, which index is in to_ignore_ind,over Z2
    std::vector<std::vector<Homology>> persistent_homologies_link(size_t to_ignore_ind) {
        std::vector<std::vector<Homology>> ans;
        if (dim == -1) {
            return ans;
        }
        const Simplex &to_ignore = simpleces[to_ignore_ind];
        ans.resize(dim - to_ignore.dim);
        std::vector<int> paired(sz, 0);
        std::vector<int> ignored(sz, 0);
        for (size_t i = 0; i != sz; ++i) {
            const Simplex &curr = simpleces[i];
            if (!std::includes(
                    curr.vertices.begin(), 
                    curr.vertices.end(),
                    to_ignore.vertices.begin(), 
                    to_ignore.vertices.end()
                    ) || i == to_ignore_ind) 
            {
                ignored[i] = 1;
            }
        }
        std::vector<std::future<void>> futures;
        int future_num = 8;
        int margin = to_ignore.dim + 1;
        for (int i = margin; i != margin + future_num; ++i) {
            futures.push_back(std::async(
                [i, future_num, this, &ans, &ignored, &paired, margin] {
                for (size_t k = i; k < this->matrices.size() - 1; k += future_num) {
                    this->matrices[k].reduce_link(ans[k - margin], 
                        this->simpleces, paired, ignored);
                }
            }
            ));
        }
        for (int i = 0; i != future_num; ++i) {
            futures[i].get();
        }
        //this->matrices[2].reduce_link(ans[1], this->simpleces, paired, ignored);
        double inf = std::numeric_limits<double>::infinity();
        for (size_t i = 0; i != sz; ++i) {
            if (!paired[i] && !ignored[i]) {
                ans[simpleces[i].dim - margin].emplace_back(i, i, simpleces[i].time, inf);
            }
        }
        return ans;
    }

    // get the ith added simplex
    Simplex& operator[](size_t i) {
        return simpleces[i];
    }
};
