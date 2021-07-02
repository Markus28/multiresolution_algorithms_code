#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/lac/vector.h>
#include "measure.h"


double make_measurement(const dealii::DoFHandler<2>& dof_handler, double (*g)(const dealii::Point<2>&),
        const dealii::Vector<double>& solution0, const dealii::Vector<double>& solution1){

    const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;

    dealii::Vector<double> solution0_values(dofs_per_cell);
    dealii::Vector<double> solution1_values(dofs_per_cell);

    dealii::QGauss<1> quadrature_formula_boundary(dof_handler.get_fe().degree+1);
    dealii::FEFaceValues<2> fe_face_values(dof_handler.get_fe(), quadrature_formula_boundary, dealii::update_values | dealii::update_quadrature_points | dealii::update_gradients | dealii::update_normal_vectors | dealii::update_JxW_values);
    double total = 0;

    for (const auto &cell : dof_handler.active_cell_iterators()){
        for(const auto &face: cell->face_iterators()){
            if(face->at_boundary()){
                fe_face_values.reinit(cell, face);
                cell->get_dof_values(solution0, solution0_values);
                cell->get_dof_values(solution1, solution1_values);

                //for (const unsigned int i : fe_face_values.dof_indices()) {
                for(unsigned int i = 0; i<dofs_per_cell; ++i){
                    for (const unsigned int q_index: fe_face_values.quadrature_point_indices()) {
                        total += g(fe_face_values.quadrature_point(q_index))*(solution0_values(i)-solution1_values(i))
                                *fe_face_values.normal_vector(q_index)*fe_face_values.shape_grad(i, q_index)*fe_face_values.JxW(q_index);
                    }
                }
            }
        }
    }

    return total;
}


InverseRadonTransform generate_psi(const std::vector<std::vector<double>>& measurements, unsigned int N, double eta, double w_l1){
    unsigned int N_phi = measurements.size() - 1;
    unsigned int N_r = measurements[0].size() - 1;
    double h = 2.0/N_r;
    InverseRadonTransform psi(N_phi, N_r, 1);
    std::vector<std::vector<double>> integrals(N_phi+1, std::vector<double>(N_r+1, 0));

    for(unsigned int i = 0; i<=N_phi; ++i){
        for(unsigned int j = 1; j<=N_r; ++j){
            integrals[i][j] = integrals[i][j-1]+0.5*h*(measurements[i][j-1]+measurements[i][j]);      //TODO: STEP SIZE
        }
    }

    for(unsigned int i = 0; i<=N_phi; ++i){
        for(unsigned int j = 0; j<=N_r; ++j){
            integrals[i][j] *= -1.0/(pow(eta, 2)*w_l1*2.0*M_PI);
        }
    }

    psi.transform(integrals, N);
    return psi;
}