#ifndef GRID_H
#define GRID_H

#include <vector>

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
    Grid2D(int nx, int ny, double dx, double dy)
        : _nx(nx), _ny(ny),
          _dx(dx), _dy(dy),
          _grid(std::vector<Node>(ny * nx))
    {}

    int nx() const { return _nx; }
    int ny() const { return _ny; }
    double dx() const { return _dx; }
    double dy() const { return _dy; }
    std::vector<Node>& grid() { return _grid; }
private:
    int _nx, _ny;
    double _dx, _dy;
    std::vector<Node> _grid;
};

#endif