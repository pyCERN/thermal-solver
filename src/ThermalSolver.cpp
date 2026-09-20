#include "ThermalSolver.h"
#include <cstdio>

void
ThermalSolver::solve()
{
    int nx = _grid.nx();
    int ny = _grid.ny();
    double dx = _grid.dx();
    double dy = _grid.dy();
    const std::vector<std::vector<Node>>& grid = _grid.grid();

    int iter = 0;
    double global_error = 0.0f;

    printf("Before iteration\n");
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            printf("%lf ", grid[j][i].T);
        }
        printf("\n");
    }

    while (global_error > _max_error && iter < _max_iter) {
        for (int j = 1; j < ny - 1; j++) {
            for (int i = 1; i < nx - 1; i++) {
                double T_old = grid[j][i].T;
                double k = grid[j][i].k;
                double Q = grid[j][i].Q;
                double T_new = 0.25 * (
                    grid[j][i+1].T + grid[j][i-1].T + grid[j+1][i].T + grid[j-1][i].T +
                    (Q * dx * dy) / k
                );

                global_error = std::max(global_error, std::abs(T_new - T_old));
            }
        }
    }

    printf("After iteration\n");
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            printf("%lf ", grid[j][i].T);
        }
        printf("\n");
    }
}