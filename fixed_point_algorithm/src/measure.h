#ifndef BIOMIMETICS_MEASURE_H
#define BIOMIMETICS_MEASURE_H

#include "InverseRadonTransform.h"
#include <deal.II/base/point.h>
#include <deal.II/dofs/dof_handler.h>


double make_measurement(const dealii::DoFHandler<2>& dof_handler, double (*g)(const dealii::Point<2>&), const dealii::Vector<double>& solution0, const dealii::Vector<double>& solution1);
InverseRadonTransform generate_psi(const std::vector<std::vector<double>>& measurements, unsigned int N, double eta, double w_l1);

#endif //BIOMIMETICS_MEASURE_H
