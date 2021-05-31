#include <iostream>
#include <set>
#include <vector>

#ifndef MULTIGRAPHS_TESTERS_H
#define MULTIGRAPHS_TESTERS_H

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

#endif //MULTIGRAPHS_TESTERS_H
