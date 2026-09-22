#include "ThermalSolver.h"
#include <cmath>
#include <cstdio>
#include <omp.h>

void
ThermalSolver::solve()
{
    int nx = _grid.nx();
    int ny = _grid.ny();
    double dx = _grid.dx();
    double dy = _grid.dy();
    std::vector<std::vector<Node>>& grid = _grid.grid();

    int iter = 0;
    double global_error = 1.0;
    std::vector<std::vector<Node>> grid_old = grid;

    while (global_error > _max_error && iter < _max_iter) {
        global_error = 0.0;
        iter++;

        for (int j = 1; j < ny - 1; j++) {
            for (int i = 1; i < nx - 1; i++) {
                if (grid_old[j][i].is_fixed_temp) {
                    continue;
                }
                double T_old = grid_old[j][i].T;
                double k = grid_old[j][i].k;
                double Q = grid_old[j][i].Q;
                double T_new = 0.25 * (
                    grid_old[j][i+1].T + grid_old[j][i-1].T + grid_old[j+1][i].T + grid_old[j-1][i].T +
                    (Q * dx * dy) / k
                );

                global_error = std::max(global_error, std::abs(T_new - T_old));
                grid[j][i].T = T_new;
            }
        }

        std::swap(grid, grid_old);
    }

    printf("Reached convergence after %d iterations with global error %lf\n",
        iter, global_error);
}

void
ThermalSolver::solve_mt()
{
#pragma omp parallel
{
    #pragma omp single
    {
        printf("Number of threads: %d\n", omp_get_num_threads());
    }
}

    int nx = _grid.nx();
    int ny = _grid.ny();
    double dx = _grid.dx();
    double dy = _grid.dy();
    std::vector<std::vector<Node>>& grid = _grid.grid();

    int iter = 0;
    double global_error = 1.0;
    std::vector<std::vector<Node>> grid_old = grid;

    while (global_error > _max_error && iter < _max_iter) {
        global_error = 0.0;
        iter++;

        double max_err_local = 0.0; // local max error per threads

        #pragma omp parallel for reduction(max:max_err_local) collapse(2)
        for (int j = 1; j < ny - 1; j++) {
            for (int i = 1; i < nx - 1; i++) {
                if (grid_old[j][i].is_fixed_temp) {
                    continue;
                }
                double T_old = grid_old[j][i].T;
                double k = grid_old[j][i].k;
                double Q = grid_old[j][i].Q;
                double T_new = 0.25 * (
                    grid_old[j][i+1].T + grid_old[j][i-1].T + grid_old[j+1][i].T + grid_old[j-1][i].T +
                    (Q * dx * dy) / k
                );

                double err = std::abs(T_new - T_old);
                if (err >  max_err_local) {
                    max_err_local = err;
                }
                grid[j][i].T = T_new;
            }
        }

        global_error = max_err_local;
        std::swap(grid, grid_old);
    }

    printf("[MT] Reached convergence after %d iterations with global error %lf\n",
        iter, global_error);
}