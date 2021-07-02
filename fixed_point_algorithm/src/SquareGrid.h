#ifndef BIOMIMETICS_SQUAREGRID_H
#define BIOMIMETICS_SQUAREGRID_H

#include <vector>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>

class SquareGrid {
public:
    SquareGrid(unsigned int N, double R, double constant);
    double evaluate(const dealii::Point<2>& pt) const;
    void attach_mesh(const dealii::DoFHandler<2>* dof_handler);
    void attach_solution(const dealii::Vector<double>& sol);
    void set_constant(double c);
    void save_grid(std::string filename) const;

private:
    unsigned int N;
    double R;
    double constant;
    double h_x;
    const dealii::DoFHandler<2>* handler = nullptr;
    std::vector<std::vector<double>> values;
};


#endif //BIOMIMETICS_SQUAREGRID_H
