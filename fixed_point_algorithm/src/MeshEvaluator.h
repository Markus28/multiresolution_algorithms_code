#ifndef BIOMIMETICS_MESHEVALUATOR_H
#define BIOMIMETICS_MESHEVALUATOR_H


#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>

class MeshEvaluator {
public:
    MeshEvaluator(){};
    double evaluate(const dealii::Point<2>& pt) const;
    void attach_mesh(const dealii::DoFHandler<2>* dof_handler);
    void attach_solution(const dealii::Vector<double>& sol);
    void set_constant(double c);

private:
    double constant = -1;
    const dealii::DoFHandler<2>* handler = nullptr;
    std::unique_ptr<dealii::GridTools::Cache<2, 2>> search_structure = nullptr;
    dealii::Vector<double> solution;
};


#endif //BIOMIMETICS_MESHEVALUATOR_H
