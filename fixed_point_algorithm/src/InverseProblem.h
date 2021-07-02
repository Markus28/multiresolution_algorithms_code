#ifndef BIOMIMETICS_INVERSEPROBLEM_H
#define BIOMIMETICS_INVERSEPROBLEM_H

#include <vector>
#include <deal.II/grid/tria.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/lac/sparse_matrix.h>
#include "InverseRadonTransform.h"
#include "ForwardProblem.h"
#include "MeshAbsorption.h"
#include "MeshEvaluator.h"
#include "SquareGrid.h"

class InverseProblem {
public:
    InverseProblem(std::vector<std::vector<double> > measurements, double eta, double q_star, double minimum,
            double maximum, double (*g)(const Point<2>&), unsigned int N, double w_l1);
    void output_results(std::string file_name) const;
    void run_iteration();

private:
    void make_grid();
    void setup_system();
    void solve();
    void assemble_system();

    std::vector<std::vector<double>> measurements;
    InverseRadonTransform psi;
    double q_star;
    double mini, maxi;
    double eta;
    double (*g)(const Point<2>&);

    MeshAbsorption approximate_absorption;
    ForwardProblem forward_problem;
    SquareGrid light_field;

    dealii::Triangulation<2> triangulation;
    dealii::FE_Q<2> fe;
    dealii::DoFHandler<2> dof_handler;
    dealii::SparsityPattern sparsity_pattern;
    dealii::SparseMatrix<double> system_matrix;
    dealii::Vector<double> system_rhs;
    dealii::Vector<double> solution;
};


#endif //BIOMIMETICS_INVERSEPROBLEM_H
