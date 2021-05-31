#include <unordered_map>

#ifndef MULTIGRAPHS_EDGE_INFLATER_H
#define MULTIGRAPHS_EDGE_INFLATER_H

class Infestation {
    // first index is less than second
    std::unordered_map<int, std::unordered_map<int, int>> func;
    std::unordered_map<int, std::unordered_map<int, double>> time;

public:
    // (creates) and returns multiplicity of edge
    int &operator()(int a, int b) {
        if (a < b)
            return func[a][b];
        return func[b][a];
    }

    double &edge_birth(int a, int b) {
        if (a < b)
            return time[a][b];
        return time[b][a];
    }

    // number of edges for two vertices a and b
    int edge_count(int a, int b) {
        if (a == b)
            return 0;
        if (a < b) {
            if (func.find(a) != func.end() && func[a].find(b) != func[a].end()) {
                return func[a][b];
            }
        }
        if (func.find(b) != func.end() && func[b].find(a) != func[b].end()) {
            return func[b][a];
        }
        return 0;
    }
};

#endif //MULTIGRAPHS_EDGE_INFLATER_H
