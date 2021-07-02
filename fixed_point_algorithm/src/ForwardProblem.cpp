#include "ForwardProblem.h"

ForwardProblem::ForwardProblem(const Absorption& absorption_model, double g(const Point<2>&))
        :fe(1), dof_handler(triangulation), absorption_model(absorption_model), g(g){
    make_grid();
    setup_system();
}

Vector<double> ForwardProblem::get_solution() const {
    return solution;
}


void ForwardProblem::make_grid()
{
    GridGenerator::hyper_cube(triangulation, -1, 1);
    //triangulation.refine_global(5);
    triangulation.refine_global(7);
    std::cout << "Number of cells: " << triangulation.n_active_cells() << std::endl;
}


void ForwardProblem::setup_system()
{
    dof_handler.distribute_dofs(fe);
    std::cout << "Number of degrees of freedom: " << dof_handler.n_dofs()
              << std::endl;
    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp);
    sparsity_pattern.copy_from(dsp);
    system_matrix.reinit(sparsity_pattern);
    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());
}


void ForwardProblem::assemble_system()
{
    /*
     * May be wrong:
     */

    solution = 0;
    system_matrix = 0;
    system_rhs = 0;

    QGauss<2> quadrature_formula(fe.degree + 1);
    QGauss<1> quadrature_formula_boundary(fe.degree+1);

    FEValues<2> fe_values(fe, quadrature_formula, update_values | update_gradients | update_quadrature_points | update_JxW_values);
    FEFaceValues<2> fe_face_values(fe, quadrature_formula_boundary, update_values | update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double> cell_rhs(dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    double absorption_value;

    for (const auto &cell : dof_handler.active_cell_iterators()){
        fe_values.reinit(cell);
        cell_matrix = 0;
        cell_rhs    = 0;

        for (const unsigned int q_index : fe_values.quadrature_point_indices()){
            absorption_value = absorption_model.absorption(fe_values.quadrature_point(q_index));
            for (const unsigned int i : fe_values.dof_indices()){
                for (const unsigned int j : fe_values.dof_indices()) {
                    cell_matrix(i, j) +=
                            (fe_values.shape_grad(i, q_index) * fe_values.shape_grad(j, q_index)
                             +(absorption_value)*fe_values.shape_value(i, q_index)*fe_values.shape_value(j, q_index))*
                            fe_values.JxW(q_index);

                }
            }
        }

        /* Added :*/
        for(const auto &face: cell->face_iterators()) {
            if (face->at_boundary()) {
                fe_face_values.reinit(cell, face);
                for(unsigned int i: fe_values.dof_indices()){
                    for(unsigned int j: fe_values.dof_indices()){
                        for (const unsigned int q_index: fe_face_values.quadrature_point_indices()) {
                            cell_matrix(i, j) += fe_face_values.shape_value(i, q_index) *
                                    fe_face_values.shape_value(j, q_index)
                                    * fe_face_values.JxW(q_index);
                        }
                    }
                }
            }
        }
        /* Added */

        for(const auto &face: cell->face_iterators()){
            if(face->at_boundary()){
                fe_face_values.reinit(cell, face);

                for (const unsigned int i : fe_values.dof_indices()){
                    for(const unsigned int q_index: fe_face_values.quadrature_point_indices()){
                        cell_rhs(i) += fe_face_values.shape_value(i, q_index)*g(fe_face_values.quadrature_point(q_index))
                                       *fe_face_values.JxW(q_index);
                    }
                }
            }
        }

        cell->get_dof_indices(local_dof_indices);

        for (const unsigned int i : fe_values.dof_indices()) {
            for (const unsigned int j : fe_values.dof_indices()) {
                system_matrix.add(local_dof_indices[i],
                                  local_dof_indices[j],
                                  cell_matrix(i, j));
            }
        }

        for (const unsigned int i : fe_values.dof_indices()) {
            system_rhs(local_dof_indices[i]) += cell_rhs(i);
        }
    }
}


void ForwardProblem::solve()
{
    SolverControl solver_control(1000, 1e-12);
    SolverCG<Vector<double>> solver(solver_control);
    solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());
}

void ForwardProblem::output_results(std::string file_name) const
{
    DataOut<2> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution, "solution");
    data_out.build_patches();
    std::ofstream output(file_name);
    data_out.write_vtk(output);
}

void ForwardProblem::run()
{
    //make_grid();
    //setup_system();
    std::cout << "Assembling forward problem"<<std::endl;
    assemble_system();
    std::cout << "Solving forward problem"<<std::endl;
    solve();
    //output_results(file_name);
}

const DoFHandler<2>& ForwardProblem::get_dof_handler() const{
    return dof_handler;
}