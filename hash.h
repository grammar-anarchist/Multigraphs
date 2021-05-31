#include <boost/functional/hash.hpp>

#ifndef MULTIGRAPHS_HASH_H
#define MULTIGRAPHS_HASH_H

namespace std {
    template<>
    struct hash<std::set<int>> {
    public:
        size_t operator()(const std::set<int> &s) const {
            size_t res = 42;
            boost::hash_combine<std::set<int>>(res, s);
            return res;
        }
    };
}

#endif //MULTIGRAPHS_HASH_H
