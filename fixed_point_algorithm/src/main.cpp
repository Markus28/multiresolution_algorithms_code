#include "BoxAbsorption.h"
#include "PerturbedAbsorption.h"
#include "ExponentialWave.h"
#include "ForwardProblem.h"
#include "InverseRadonTransform.h"
#include "SheppLoganAbsorption.h"
#include "measure.h"
#include "InverseProblem.h"
#include "FuzzyCircleAbsorption.h"

#include <vector>

void save_grid(const std::vector<std::vector<double>>& data, std::string file_name, unsigned int precision){
    std::ofstream my_file;
    my_file.open (file_name);
    my_file.precision(precision);
    my_file<<data.size()<<","<<data[0].size()<<"\n";

    for(unsigned int i = 0; i<data.size(); ++i){
        assert(data[i].size() == data[0].size());
        for(unsigned int j = 0; j<data[0].size(); ++j){
            my_file<< i << "," << j << "," << data[i][j] << "\n";
        }
    }

    my_file.close();
}

void save_grid(const std::vector<std::vector<double>>& data, std::string file_name){
    save_grid(data, file_name, 15);
}

std::vector<std::vector<double>> read_grid(std::string file_name){
    std::ifstream my_file;
    my_file.open(file_name);
    char comma;
    unsigned int i, j;
    double value;
    unsigned int n, m;
    my_file >> n >> comma >> m;
    std::vector<std::vector<double>> result(n, std::vector<double>(m, 0));
    while(my_file >> i >> comma >> j >> comma >> value){
        result[i][j] = value;
    }

    my_file.close();

    return result;
}


double exponential_g(const Point<2>& p){
    double x = p(0);
    double y = p(1);
    double x0 = 1;
    double y0 = 0;
    double eps = 0.05;
    return exp(-0.5*(pow((x-x0)/eps,2)+pow((y-y0)/eps,2)))/(pow(eps, 2)*2*M_PI);
}


double cut_g(const Point<2>& p){
    if(p(0)<=-0.9999999)
        return 1;
    return 0;
}


void save_psi(){
    double R_0 = 1;
    unsigned int N_phi = 64;
    unsigned int N_r = 255;
    unsigned int N = 400;

    double h_phi = 2.0*M_PI/(N_phi+1);
    double h_r = 2.0*R_0/N_r;
    double c = 1;

    std::vector<std::vector<double>> measurements(N_phi+1, std::vector<double>(N_r+1, 0));

    BoxAbsorption normal_absorption;
    PerturbedAbsorption<BoxAbsorption, ExponentialWave> my_absorption(Point<2>(1, 0), c, 0.02);

    ForwardProblem problem(my_absorption, cut_g);
    ForwardProblem problem_star(normal_absorption, cut_g);

    const DoFHandler<2>& my_handler = problem.get_dof_handler();

    problem_star.run();

    Vector<double> solution0 = problem_star.get_solution();

    for(unsigned int i = 0; i<=N_phi; ++i){
        my_absorption.set_center({R_0*cos(i*h_phi), R_0*sin(i*h_phi)});
        for(unsigned int j = 1; j<=N_r; ++j){
            my_absorption.set_time(j*h_r/c);
            problem.run();
            measurements[i][j] = make_measurement(my_handler, cut_g, solution0, problem.get_solution());
        }
    }

    save_grid(measurements, "measurements.txt");

    InverseRadonTransform psi = generate_psi(measurements, N, 0.02, ExponentialWave::w_l1);
    save_grid(psi.get_raw_solution(), "psi.txt");
}


void test_inverse_radon(unsigned int N_phi, unsigned int N_r, unsigned int N_quadrature, unsigned int subdivision){
    SheppLoganAbsorption my_phantom;
    double h_phi = 2*M_PI/(N_phi + 1);
    double h_subdivsion_phi = h_phi/subdivision;
    double h_r = 2.0/N_r;
    double h_phi_quad = 2*M_PI/N_quadrature;
    double approximate_integral;
    std::vector<std::vector<double>> F(N_phi+1, std::vector<double>(N_r+1, 0));

    std::cout<< "Calculating Integral";

    double val1;
    double val2;

    for(unsigned int k = 0; k<=N_phi; ++k){
        for(unsigned int m=1; m<=N_r; ++m){
            approximate_integral = 0;

            for(unsigned int i = 0; i<N_quadrature; ++i){
                val1 = my_phantom.absorption({cos(k*h_phi) + m*h_r*cos(i*h_phi_quad), sin(k*h_phi) + m*h_r*sin(i*h_phi_quad)});
                val2 = my_phantom.absorption({cos(k*h_phi) + m*h_r*cos((i+1)*h_phi_quad), sin(k*h_phi) + m*h_r*sin((i+1)*h_phi_quad)});

                if(val1==val2) {
                    approximate_integral += 2 * M_PI * h_r * m / N_quadrature * (0.5 * val1 + 0.5 * val2);
                }

                else{       //Subdivide
                    for(unsigned int j = 0; j<subdivision;++j){
                        val1 = my_phantom.absorption({cos(k*h_phi) + m*h_r*cos(i*h_phi_quad+j*h_subdivsion_phi), sin(k*h_phi) + m*h_r*sin(i*h_phi_quad + j*h_subdivsion_phi)});
                        val2 = my_phantom.absorption({cos(k*h_phi) + m*h_r*cos(i*h_phi_quad+(j+1)*h_subdivsion_phi), sin(k*h_phi) + m*h_r*sin(i*h_phi_quad+(j+1)*h_subdivsion_phi)});

                        approximate_integral += 2*M_PI*h_r*m/(N_quadrature*subdivision)*(0.5*val1 + 0.5*val2);
                    }
                }
            }

            F[k][m] = approximate_integral/(2*M_PI*h_r*m);
        }
        F[k][0]= my_phantom.absorption({cos(k*h_phi), sin(k*h_phi)});
    }

    std::cout<<"Computing inverse";

    InverseRadonTransform transformer(N_phi, N_r, 1);

    transformer.transform(F, 800);

    save_grid(transformer.get_raw_solution(), "inverted_radon_transform.txt");
    save_grid(F, "radon_transform.txt");
}

void visualize_forward(){
    unsigned int N_r = 200;
    double c = 1;
    double h_r = 2.0/N_r;

    BoxAbsorption normal_absorption;
    PerturbedAbsorption<BoxAbsorption, ExponentialWave> my_absorption(Point<2>(1, 0), c, 0.02);
    my_absorption.set_center({cos(2.3), sin(2.3)});

    ForwardProblem problem(my_absorption, cut_g);
    ForwardProblem problem_star(normal_absorption, cut_g);

    SquareGrid solution_grid(400, 1, 0);
    const DoFHandler<2>& my_handler = problem.get_dof_handler();
    solution_grid.attach_mesh(&my_handler);

    problem_star.run();

    Vector<double> solution0 = problem_star.get_solution();
    solution_grid.attach_solution(solution0);
    solution_grid.save_grid("solution_0.txt");

    for(unsigned int j = 1; j<=N_r; ++j) {
        my_absorption.set_time(j * h_r / c);
        problem.run();
        solution_grid.attach_solution(problem.get_solution());
        solution_grid.save_grid("solution_"+std::to_string(j)+".txt");
    }
}


void reconstruct(){
    std::vector<std::vector<double>> measurements = read_grid("measurements.txt");
    unsigned int N = 500;
    double q_star = 1;
    double mini = 0.95;
    double maxi = 1.05;

    InverseProblem inv_problem(measurements, 0.02, q_star, mini, maxi, cut_g, N, ExponentialWave::w_l1);
    std::cout<<"First iteration"<<std::endl;
    inv_problem.run_iteration();
    std::cout<<"Second iteration"<<std::endl;
    inv_problem.run_iteration();
    inv_problem.output_results("solution.vtk");
}



int main()
{
    std::cout<<"Creating measurements..."<<std::endl;
    save_psi();

    std::cout<<std::endl<<"Reconstructing absorption..."<<std::endl;
    reconstruct();
    visualize_forward();
    return 0;
}