#include "MeshEvaluator.h"


void MeshEvaluator::attach_mesh(const dealii::DoFHandler<2>* dof_handler){
    handler = dof_handler;
    search_structure = std::make_unique<dealii::GridTools::Cache<2,2>>(handler->get_triangulation());
}

void MeshEvaluator::attach_solution(const dealii::Vector<double> &sol) {
    assert(handler!=nullptr);
    solution = sol;
}

void MeshEvaluator::set_constant(double c){
    constant = c;
}

double MeshEvaluator::evaluate(const dealii::Point<2> &p) const {
    dealii::DoFHandler<2>::active_cell_iterator surrounding_cell = handler->begin_active();
    dealii::Point<2> unit_point;
    double total = 0;

    //surrounding_cell->vertex()

    try {
        std::tie(surrounding_cell, unit_point) = dealii::GridTools::find_active_cell_around_point(*search_structure, p);
    }

    catch(...){             //TODO: Only catch no point found exception
        return constant;
    }

    assert(surrounding_cell->point_inside(p));

    unit_point = dealii::GeometryInfo<2>::project_to_unit_cell(unit_point);

    const dealii::FiniteElement<2>& fe = handler->get_fe();

    std::vector<unsigned int> dof_indices(fe.dofs_per_cell, 0);

    surrounding_cell->get_dof_indices(dof_indices);

    for(unsigned int i = 0; i<dof_indices.size(); ++i){
        total += solution(dof_indices[i])*fe.shape_value(i, unit_point);
    }

    return total;
}