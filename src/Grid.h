#ifndef GRID_H
#define GRID_H

#include <vector>

enum class MaterialType
{
    Silicon,
    SiO2,
    HeatSource
};

struct Node
{
    double T = 300.0; // temperature in Kelvin
    double k = 148.0; // thermal conductivity (W/mK)
    double Q = 0.0;
    bool is_fixed_temp = false;
};

class Grid2D
{
public:
    Grid2D(int nx, int ny, int nz, double dx, double dy, double dz)
        : _nx(nx), _ny(ny), _nz(nz),
          _dx(dx), _dy(dy), _dz(dz),
          _grid(std::vector<std::vector<Node>>(ny))
    {
        for (int i = 0; i < ny; i++) {
            _grid[i].assign(nx, Node());
        }
    }

    int nx() const { return _nx; }
    int ny() const { return _ny; }
    int nz() const { return _nz; }
    double dx() const { return _dx; }
    double dy() const { return _dy; }
    double dz() const { return _dz; }
    std::vector<std::vector<Node>>& grid() { return _grid; }
private:
    int _nx, _ny, _nz;
    double _dx, _dy, _dz;
    std::vector<std::vector<Node>> _grid;
};

#endif