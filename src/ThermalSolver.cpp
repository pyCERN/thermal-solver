#include "ThermalSolver.h"
#include <cmath>
#include <cstdio>
#include <omp.h>

void
ThermalSolver::solve_jacobi()
{
    int nx = _grid.nx();
    int ny = _grid.ny();
    double dx = _grid.dx();
    double dy = _grid.dy();
    std::vector<Node>& grid = _grid.grid();

    int iter = 0;
    double global_error = 1.0;
    std::vector<Node> grid_old = grid;

    while (global_error > _max_error && iter < _max_iter) {
        global_error = 0.0;
        iter++;

        for (int j = 1; j < ny - 1; j++) {
            for (int i = 1; i < nx - 1; i++) {
                if (grid_old[j*nx+i].is_fixed_temp) {
                    continue;
                }
                double T_old = grid_old[j*nx+i].T;
                double k = grid_old[j*nx+i].k;
                double Q = grid_old[j*nx+i].Q;
                double T_new = 0.25 * (
                    grid_old[j*nx+i+1].T + grid_old[j*nx+i-1].T + grid_old[(j+1)*nx+i].T + grid_old[(j-1)*nx+i].T +
                    (Q * dx * dy) / k
                );

                global_error = std::max(global_error, std::abs(T_new - T_old));
                grid[j*nx+i].T = T_new;
            }
        }

        std::swap(grid, grid_old);
    }

    printf("[Gauss-Seidel] Reached convergence after %d iterations with global error %lf\n",
        iter, global_error);
}

void
ThermalSolver::solve_gauss_seidel()
{
    int nx = _grid.nx();
    int ny = _grid.ny();
    double dx = _grid.dx();
    double dy = _grid.dy();
    std::vector<Node>& grid = _grid.grid();

    int iter = 0;
    double global_error = 1.0;

    while (global_error > _max_error && iter < _max_iter) {
        global_error = 0.0;
        iter++;

        for (int j = 1; j < ny - 1; j++) {
            for (int i = 1; i < nx - 1; i++) {
                if (grid[j*nx+i].is_fixed_temp) {
                    continue;
                }
                double T_old = grid[j*nx+i].T;
                double k = grid[j*nx+i].k;
                double Q = grid[j*nx+i].Q;
                double T_new = 0.25 * (
                    grid[j*nx+i+1].T + grid[j*nx+i-1].T + grid[(j+1)*nx+i].T + grid[(j-1)*nx+i].T +
                    (Q * dx * dy) / k
                );

                global_error = std::max(global_error, std::abs(T_new - T_old));
                grid[j*nx+i].T = T_new;
            }
        }
    }

    printf("[Jacobi] Reached convergence after %d iterations with global error %lf\n",
        iter, global_error);
}

void
ThermalSolver::solve_sor()
{
    int nx = _grid.nx();
    int ny = _grid.ny();
    double dx = _grid.dx();
    double dy = _grid.dy();
    std::vector<Node>& grid = _grid.grid();

    int iter = 0;
    double global_error = 1.0;

    double relax_factor = 1.5;

    while (global_error > _max_error && iter < _max_iter) {
        global_error = 0.0;
        iter++;

        for (int j = 1; j < ny - 1; j++) {
            for (int i = 1; i < nx - 1; i++) {
                if (grid[j*nx+i].is_fixed_temp) {
                    continue;
                }
                double T_old = grid[j*nx+i].T;
                double k = grid[j*nx+i].k;
                double Q = grid[j*nx+i].Q;
                double T_new = 0.25 * (
                    grid[j*nx+i+1].T + grid[j*nx+i-1].T + grid[(j+1)*nx+i].T + grid[(j-1)*nx+i].T +
                    (Q * dx * dy) / k
                );

                T_new = T_old + relax_factor * (T_new - T_old);

                global_error = std::max(global_error, std::abs(T_new - T_old));
                grid[j*nx+i].T = T_new;
            }
        }
    }

    printf("[SOR] Reached convergence after %d iterations with global error %lf\n",
        iter, global_error);
}

void
ThermalSolver::solve_mt()
{
    int nx = _grid.nx();
    int ny = _grid.ny();
    double dx = _grid.dx();
    double dy = _grid.dy();
    std::vector<Node>& grid = _grid.grid();

    int iter = 0;
    double global_error = 1.0;
    std::vector<Node> grid_old = grid;
    bool terminate = false;
    double max_err_local = 0.0;

    #pragma omp parallel
    {
        #pragma omp single
        {
            printf("Number of threads: %d\n", omp_get_num_threads());
        }

        while (true) {
            #pragma omp single
            {
                global_error = 0.0;
                iter++;
                max_err_local = 0.0;
                terminate = false;
            }
            #pragma omp barrier

            #pragma omp for reduction(max:max_err_local) collapse(2)
            for (int j = 1; j < ny - 1; j++) {
                for (int i = 1; i < nx - 1; i++) {
                    if (grid_old[j*nx+i].is_fixed_temp) {
                        continue;
                    }
                    double T_old = grid_old[j*nx+i].T;
                    double k = grid_old[j*nx+i].k;
                    double Q = grid_old[j*nx+i].Q;
                    double T_new = 0.25 * (
                        grid_old[j*nx+i+1].T + grid_old[j*nx+i-1].T + grid_old[(j+1)*nx+i].T + grid_old[(j-1)*nx+i].T +
                        (Q * dx * dy) / k
                    );

                    double err = std::abs(T_new - T_old);
                    if (err > max_err_local) {
                        max_err_local = err;
                    }
                    grid[j*nx+i].T = T_new;
                }
            }

            #pragma omp single
            {
                std::swap(grid, grid_old);
                global_error = max_err_local;
                if (global_error <= _max_error || iter >= _max_iter) {
                    terminate = true;
                }
            }
            #pragma omp barrier

            if (terminate) {
                break;
            }
        }
    }

    printf("[MT] Reached convergence after %d iterations with global error %lf\n",
        iter, global_error);
}