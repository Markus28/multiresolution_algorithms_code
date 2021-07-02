#include <deal.II/grid/grid_generator.h>

#include "InverseProblem.h"
#include "measure.h"


InverseProblem::InverseProblem(std::vector<std::vector<double> > measurements, double eta, double q_star, double mini, double maxi, double (*g)(const Point<2>&), unsigned int N, double w_l1):
measurements(measurements), eta(eta), q_star(q_star), g(g), forward_problem(approximate_absorption, g), fe(1), dof_handler(triangulation),
psi(generate_psi(measurements, N, eta, w_l1)), mini(mini), maxi(maxi), light_field(500, 1, 0), approximate_absorption(500, 1, q_star){
    make_grid();
    setup_system();
    approximate_absorption.set_constant(q_star);        //TODO: REMOVE
    approximate_absorption.attach_mesh(&dof_handler);
    approximate_absorption.set_bounds(mini, maxi);
    light_field.attach_mesh(&forward_problem.get_dof_handler());
}

void InverseProblem::make_grid(){
    dealii::GridGenerator::hyper_ball_balanced(triangulation, {0,0}, 0.9);
    triangulation.refine_global(7);
}

void InverseProblem::setup_system(){
    dof_handler.distribute_dofs(fe);
    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp);
    sparsity_pattern.copy_from(dsp);
    system_matrix.reinit(sparsity_pattern);
    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());
}

void InverseProblem::solve()
{
    SolverControl solver_control(18000, 1e-12);
    SolverCG<Vector<double>> solver(solver_control);
    solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());
}

void InverseProblem::assemble_system(){
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

    Point<2> quad_point;
    std::pair<double, double> value;
    dealii::Tensor<1, 2> shape_grad;

    for (const auto &cell : dof_handler.active_cell_iterators()){
        fe_values.reinit(cell);
        cell_matrix = 0;
        cell_rhs    = 0;

        for (const unsigned int q_index : fe_values.quadrature_point_indices()){
            for (const unsigned int i : fe_values.dof_indices()){
                for (const unsigned int j : fe_values.dof_indices()) {
                    cell_matrix(i, j) +=
                            fe_values.shape_grad(i, q_index) * fe_values.shape_grad(j, q_index) * pow(light_field.evaluate(fe_values.quadrature_point(q_index)),2)*fe_values.JxW(q_index);
                }       //TODO: SQUARE light_field!

                quad_point = fe_values.quadrature_point(q_index);
                value = psi.interpolate_gradient(quad_point(0), quad_point(1));
                shape_grad = fe_values.shape_grad(i, q_index);
                cell_rhs(i) += - (shape_grad[0]*value.first + shape_grad[1]*value.second)
                               *fe_values.JxW(q_index);
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


    std::map<types::global_dof_index, double> boundary_values;
    VectorTools::interpolate_boundary_values(dof_handler,
                                             0,
                                             Functions::ConstantFunction<2>(q_star),
                                             boundary_values);

    MatrixTools::apply_boundary_values(boundary_values,
                                       system_matrix,
                                       solution,
                                       system_rhs);
}

void InverseProblem::run_iteration(){
    forward_problem.run();
    std::cout<<"Attaching solution";
    light_field.attach_solution(forward_problem.get_solution());
    //light_field.save_grid("light_field.txt");
    std::cout<<"Assembling inverse problem"<<std::endl;
    assemble_system();
    std::cout<<"Solving inverse problem"<<std::endl;
    solve();
    std::cout<<"Attaching solution";
    approximate_absorption.attach_solution(solution);
    approximate_absorption.save_grid("absorption.txt");
}

void InverseProblem::output_results(std::string file_name) const
{
    DataOut<2> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution, "solution");
    data_out.build_patches();
    std::ofstream output(file_name);
    data_out.write_vtk(output);
}