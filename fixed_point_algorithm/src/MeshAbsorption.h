#ifndef BIOMIMETICS_MESHABSORPTION_H
#define BIOMIMETICS_MESHABSORPTION_H


#include <deal.II/base/point.h>
#include "Absorption.h"
#include "MeshEvaluator.h"
#include "SquareGrid.h"

class MeshAbsorption: public Absorption {
public:
    MeshAbsorption(unsigned int N, double R0, double c): evaluator(N, R0, c) {};
    double absorption(const dealii::Point<2>& p) const override;
    void attach_mesh(const dealii::DoFHandler<2>* dof_handler);
    void attach_solution(const dealii::Vector<double>& sol);
    void set_constant(double val);
    void enable_constant();
    void set_bounds(double minimum, double maximum);
    void save_grid(std::string filename) const;

private:
    double mini, maxi;
    bool use_constant = true;
    double constant = 0;
    SquareGrid evaluator;
};


#endif //BIOMIMETICS_MESHABSORPTION_H
