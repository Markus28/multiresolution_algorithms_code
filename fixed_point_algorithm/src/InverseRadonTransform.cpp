#include <cassert>
#include <math.h>
#include "InverseRadonTransform.h"


InverseRadonTransform::InverseRadonTransform(unsigned int N_phi, unsigned int N_r, double r0): a(N_r+1, std::vector<double>(N_r+1, 0)), b(N_r+1, std::vector<double>(N_r+1, 0)){
    R_0 = r0;
    this->N_phi = N_phi;
    this->N_r = N_r;
    h_phi = 2.0*M_PI/(N_phi+1);
    h_r = 2.0*R_0/N_r;
    for(unsigned int m=0; m<=N_r; ++m){
        for(unsigned int m_prime = 0; m_prime<=N_r; ++m_prime){
            a[m][m_prime] = a_value(h_r, m, m_prime+1) - a_value(h_r, m, m_prime);
            b[m][m_prime] = -h_r*m_prime*a[m][m_prime] + b_value(h_r, m, m_prime+1) - b_value(h_r, m, m_prime);
        }
    }
}

double InverseRadonTransform::a_value(double h_r, unsigned int m, unsigned int m_eval){
    if(m==m_eval && m_eval!=0)
        return 2*m_eval*h_r*log(2*m_eval*h_r)-2*h_r*m_eval;
    else if(m==m_eval)
        return 0;
    return (h_r*m_eval - h_r*m)*log(abs(h_r*m - h_r*m_eval)) + (h_r*m + h_r*m_eval)*log(abs(h_r*m + h_r*m_eval)) - 2*h_r*m_eval;
}

double InverseRadonTransform::b_value(double h_r, unsigned int m, unsigned int m_eval){
    if(m==m_eval)
        return -0.5*pow(h_r*m_eval, 2);
    return 0.5*((pow(h_r*m_eval,2) - pow(h_r*m,2))*log(abs(pow(h_r*m,2) - pow(h_r*m_eval,2))) - pow(h_r*m_eval, 2));
}

void InverseRadonTransform::transform(std::vector<std::vector<double>> F, unsigned int N){
    assert(F.size() == N_phi+1 && F[0].size() == N_r+1);
    h_x = 2.0*R_0/N;
    solution = std::vector<std::vector<double>>(N+1, std::vector<double>(N+1, 0));

    std::vector<double> temp_data(N_r+1, 0);
    double distance_to_k = 0;
    double x_i = 0;
    double y_i = 0;

    for(unsigned int k = 0; k<=N_phi; ++k){
        for(unsigned int m = 1; m<N_r; ++m){
            temp_data[m] = (m+0.5)*F[k][m+1] + (m-0.5)*F[k][m-1] - 2*m*F[k][m];
        }

        //TODO: Divide by h_r?
        F[k][0] = 0.5*F[k][1]/h_r;
        F[k][N_r] = ((N_r-0.5)*F[k][N_r-1] - 2*N_r*F[k][N_r])/h_r;

        for(unsigned int j=1; j<N_r; ++j){
            F[k][j] = temp_data[j]/h_r;
        }

        for(unsigned int m = 0; m<=N_r; ++m){
            temp_data[m] = 0;
            for(unsigned int m_prime = 0; m_prime<N_r; ++m_prime){
                temp_data[m] += a[m][m_prime]*F[k][m_prime] + b[m][m_prime]*(F[k][m_prime+1]-F[k][m_prime])/h_r;
            }
        }

        for(unsigned int m = 0; m<=N_r; ++m){
            F[k][m] = temp_data[m];
        }
    }

    unsigned int m = 0;
    double T = 0;

    for(unsigned int i1 = 0; i1<=N; ++i1){
        for(unsigned int i2=0; i2<=N; ++i2){
            x_i = -R_0 + i1* h_x;
            y_i = -R_0 + i2*h_x;

            if(sqrt(pow(x_i, 2)+pow(y_i, 2))<R_0) {
                for (unsigned int k = 0; k <= N_phi; ++k) {
                    distance_to_k = sqrt(pow(x_i - R_0 * cos(k * h_phi), 2) + pow(y_i - R_0 * sin(k * h_phi), 2));
                    m = floor(distance_to_k / h_r);
                    assert(m < N_r);
                    T = F[k][m] + (distance_to_k - m * h_r) * (F[k][m + 1] - F[k][m]) / h_r;        //IS distance_to_k right?
                    solution[i1][i2] += T / (N_phi + 1);
                }
            }
        }
    }
}


double InverseRadonTransform::interpolate_value(double x, double y) const{
    if(x< -R_0 || x>R_0 || y< -R_0 || y > R_0){
        throw;
    }

    double estimated_i1 = (x+R_0)/h_x;
    double estimated_i2 = (y+R_0)/h_x;

    double x_to_corner = (estimated_i1-floor(estimated_i1))*h_x;
    double y_to_corner = (estimated_i2-floor(estimated_i2))*h_x;

    double val_1 = solution[floor(estimated_i1)][floor(estimated_i2)];
    double val_2 = solution[floor(estimated_i1)][ceil(estimated_i2)];
    double val_3 = solution[ceil(estimated_i1)][floor(estimated_i2)];
    double val_4 = solution[ceil(estimated_i1)][ceil(estimated_i2)];

    return 1/pow(h_x, 2)*(val_1*(h_x-x_to_corner)*(h_x-y_to_corner) + val_2*(h_x-x_to_corner)*y_to_corner + val_3*x_to_corner*(h_x-y_to_corner) + val_4*x_to_corner*y_to_corner);
}

std::pair<double, double> InverseRadonTransform::interpolate_gradient(double x, double y) const{
    return std::make_pair((interpolate_value(x+h_x, y)-interpolate_value(x-h_x, y))/(2*h_x),
                          (interpolate_value(x, y+h_x)-interpolate_value(x, y-h_x))/(2*h_x));
}

std::vector<std::vector<double>> InverseRadonTransform::get_raw_solution() const{
    return solution;
}