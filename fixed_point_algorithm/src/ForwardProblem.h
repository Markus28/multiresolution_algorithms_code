#ifndef MAIN_FORWARDPROBLEM_H
#define MAIN_FORWARDPROBLEM_H

#include <deal.II/base/point.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/numerics/data_out.h>

#include "Absorption.h"
#include <iostream>
#include <fstream>

using namespace dealii;

class ForwardProblem{
public:
    ForwardProblem(const Absorption& absorption_model, double (*g)(const Point<2>&));
    void run();
    void output_results(std::string file_name) const;
    Vector<double> get_solution() const;
    const DoFHandler<2>& get_dof_handler() const;

private:
    void make_grid();
    void setup_system();
    void assemble_system();
    void solve();

    Triangulation<2> triangulation;
    FE_Q<2> fe;
    DoFHandler<2> dof_handler;
    SparsityPattern sparsity_pattern;
    SparseMatrix<double> system_matrix;
    Vector<double> solution;
    Vector<double> system_rhs;

    const Absorption& absorption_model;
    double (*g)(const Point<2>&);
};

#endif //MAIN_FORWARDPROBLEM_H
