#include <cmath>


#ifndef MULTIGRAPHS_HOMOLOGY_H
#define MULTIGRAPHS_HOMOLOGY_H

class Homology {
public:
    double birth;
    double death;
    unsigned coef = 1;

    Homology(): birth(INFINITY) {}

    Homology(double birth, double death):
        birth(birth), death(death) {}
    Homology(double birth, double death, unsigned coef):
            birth(birth), death(death), coef(coef) {}

    explicit Homology(double birth):
        birth(birth), death(INFINITY) {}
    Homology(double birth, unsigned multiplier):
            birth(birth), death(INFINITY), coef(1) {}

    bool is_homology() {
        return birth != INFINITY;
    }
};

#endif //MULTIGRAPHS_HOMOLOGY_H
