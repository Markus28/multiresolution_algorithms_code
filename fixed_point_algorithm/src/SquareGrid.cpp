#include "SquareGrid.h"
#include <deal.II/fe/mapping.h>
#include <fstream>


SquareGrid::SquareGrid(unsigned int N, double R, double constant): N(N), R(R), constant(constant), values(N, std::vector<double>(N, constant)) {
    h_x = 2.0*R/(N-1);
}

void SquareGrid::set_constant(double c){
    constant = c;
    for(unsigned int i = 0; i<N; ++i){
        for(unsigned int j = 0; j<N; ++j){
            values[i][j] = constant;
        }
    }
}

double SquareGrid::evaluate(const dealii::Point<2>& pt) const{
    double x = pt[0];
    double y = pt[1];

    if(x< -R || x>R || y< -R || y > R){
        throw;
    }

    double estimated_i1 = (x+R)/h_x;
    double estimated_i2 = (y+R)/h_x;

    double x_to_corner = (estimated_i1-floor(estimated_i1))*h_x;
    double y_to_corner = (estimated_i2-floor(estimated_i2))*h_x;

    double val_1 = values[floor(estimated_i1)][floor(estimated_i2)];
    double val_2 = values[floor(estimated_i1)][ceil(estimated_i2)];
    double val_3 = values[ceil(estimated_i1)][floor(estimated_i2)];
    double val_4 = values[ceil(estimated_i1)][ceil(estimated_i2)];

    return 1/pow(h_x, 2)*(val_1*(h_x-x_to_corner)*(h_x-y_to_corner) + val_2*(h_x-x_to_corner)*y_to_corner + val_3*x_to_corner*(h_x-y_to_corner) + val_4*x_to_corner*y_to_corner);
}

void SquareGrid::attach_mesh(const dealii::DoFHandler<2>* dof_handler){
    handler = dof_handler;
}

void SquareGrid::attach_solution(const dealii::Vector<double>& sol){
    assert(handler != nullptr);
    unsigned int i_min, i_max, j_min, j_max;
    double x, y;
    double total;
    dealii::BoundingBox<2> bounds;

    const dealii::FiniteElement<2>& fe = handler->get_fe();

    dealii::MappingQ1<2> mapping;

    std::vector<unsigned int> dof_indices(fe.dofs_per_cell, 0);

    for(unsigned int i = 0; i<N; ++i){
        for(unsigned int j = 0; j<N; ++j){
            values[i][j] = constant;
        }
    }

    for (const auto &cell : handler->active_cell_iterators()){
        bounds = cell->bounding_box();
        i_min = std::max(0, (int) ceil((bounds.lower_bound(0)+R)/h_x));
        i_max = std::min(N-1, (unsigned int) floor((bounds.upper_bound(0)+R)/h_x));
        j_min = std::max(0, (int) ceil((bounds.lower_bound(1)+R)/h_x));
        j_max = std::min(N-1, (unsigned int) floor((bounds.upper_bound(1)+R)/h_x));


        cell->get_dof_indices(dof_indices);

        for(unsigned int i = i_min; i<= i_max; ++i){
            for(unsigned int j = j_min; j<=j_max; ++j){
                x = -R + h_x*i;
                y = -R + h_x*j;

                if(cell->point_inside({x, y})){
                    total = 0;

                    for(unsigned int k = 0; k<dof_indices.size(); ++k){
                        total += sol(dof_indices[k])*fe.shape_value(k, mapping.transform_real_to_unit_cell(cell, {x, y}));
                        //total += sol(dof_indices[k])*fe.shape_value(k, dealii::Mapping<2>::transform_real_to_unit_cell(cell, {x,y}));
                    }

                    values[i][j] = total;
                }
            }
        }
    }
}

void SquareGrid::save_grid(std::string file_name) const{
    std::ofstream my_file;
    my_file.open (file_name);
    my_file.precision(15);
    my_file<<values.size()<<","<<values[0].size()<<"\n";

    for(unsigned int i = 0; i<values.size(); ++i){
        assert(values[i].size() == values[0].size());
        for(unsigned int j = 0; j<values[0].size(); ++j){
            my_file<< i << "," << j << "," << values[i][j] << "\n";
        }
    }

    my_file.close();
}