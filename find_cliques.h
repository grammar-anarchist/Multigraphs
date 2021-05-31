#ifndef MULTIGRAPHS_FIND_CLIQUES_H
#define MULTIGRAPHS_FIND_CLIQUES_H

// Bron-Kerbosch algorithm
void CliqueMultiComplex::find_all_cliques(
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
            complex[dim - 1][prev_ind].push_parent(dest);
            complex[dim - 1][prev_ind].push_vertex(candidates[i], time);
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
        complex[dim - 1].back().push_parent(curr_ind);
        complex[dim - 1].back().push_vertex(candidates[i], time);
        Simplex &curr_simplex = complex[dim].back();

        search[dim][curr] = curr_ind;

        curr.erase(a);
        dest = search[dim - 1][curr];
        complex[dim - 1][dest].push_parent(curr_ind);
        complex[dim - 1][dest].push_vertex(a, time);
        curr_simplex.push_sub(dest);
        curr.insert(a);

        curr.erase(b);
        dest = search[dim - 1][curr];
        complex[dim - 1][dest].push_parent(curr_ind);
        complex[dim - 1][dest].push_vertex(b, time);
        curr_simplex.push_sub(dest);
        curr.insert(b);

        curr_simplex.push_sub(prev_ind);

        std::vector<int> new_candidates;

        // можно вынести в функцию
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

#endif //MULTIGRAPHS_FIND_CLIQUES_H
