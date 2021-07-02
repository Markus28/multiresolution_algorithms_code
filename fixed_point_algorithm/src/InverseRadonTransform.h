#ifndef MAIN_INVERSERADONTRANSFORM_H
#define MAIN_INVERSERADONTRANSFORM_H


#include <utility>
#include <vector>

class InverseRadonTransform{
public:
    InverseRadonTransform(unsigned int N_phi, unsigned int N_r, double r0);
    void transform(std::vector<std::vector<double>> data, unsigned int N);
    double interpolate_value(double x, double y) const;
    std::pair<double, double> interpolate_gradient(double x, double y) const;
    std::vector<std::vector<double>> get_raw_solution() const;

private:
    static double a_value(double r_phi, unsigned int m, unsigned int m_eval);
    static double b_value(double r_phi, unsigned int m, unsigned int m_eval);
    double h_phi;
    double h_r;
    unsigned int N_phi;
    unsigned int N_r;
    double R_0;
    double h_x;
    std::vector<std::vector<double>> a;
    std::vector<std::vector<double>> b;
    std::vector<std::vector<double>> solution;
};


#endif //MAIN_INVERSERADONTRANSFORM_H
