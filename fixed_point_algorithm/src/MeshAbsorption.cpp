#include "MeshAbsorption.h"
#include "MeshEvaluator.h"


void MeshAbsorption::set_constant(double val){
    constant = val;
    evaluator.set_constant(constant);
}

void MeshAbsorption::set_bounds(double minimum, double maximum){
    mini = minimum;
    maxi = maximum;
}

void MeshAbsorption::enable_constant(){
    use_constant = true;
}

void MeshAbsorption::attach_mesh(const dealii::DoFHandler<2>* dof_handler){
    evaluator.attach_mesh(dof_handler);
}

void MeshAbsorption::attach_solution(const dealii::Vector<double> &sol) {
    use_constant = false;
    evaluator.attach_solution(sol);
}

double MeshAbsorption::absorption(const dealii::Point<2> &p) const {
    if(use_constant)
        return constant;

    return std::max(std::min(evaluator.evaluate(p), maxi), mini);
}


void MeshAbsorption::save_grid(std::string filename) const{
    evaluator.save_grid(filename);
}
