#ifndef THERMAL_SOLVER_H
#define THERMAL_SOLVER_H

#include "Grid.h"

class ThermalSolver
{
public:
    ThermalSolver(Grid2D& grid, double max_error = 1e-5, int max_iter = 10000)
        : _grid(grid), _max_error(max_error), _max_iter(max_iter)
    {}
    void solve();
private:
    Grid2D& _grid;
    double _max_error;
    int _max_iter;
};

#endif