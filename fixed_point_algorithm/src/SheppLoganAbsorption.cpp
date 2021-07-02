#include "SheppLoganAbsorption.h"


double SheppLoganAbsorption::absorption(const dealii::Point<2>& p) const{
    double total = 0;
    std::pair<double, double> pt = std::make_pair(p(0), p(1));

    for(unsigned int i=0; i<ellipses.size(); ++i){
        if(ellipses[i].contains_point(pt)){
            total += weights[i];
        }
    }

    return total;
}

const std::vector<double> SheppLoganAbsorption::weights = {2, -0.98, -0.02, -0.02, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};
const std::vector<Ellipse> SheppLoganAbsorption::ellipses = {Ellipse({0,0}, 0.69, 0.92, 0), Ellipse({0,-0.0184}, 0.6624, 0.874, 0),
                                                             Ellipse({0.22,0}, 0.11, 0.31, -18), Ellipse({-0.22,0}, 0.16, 0.41, 18),
                                                             Ellipse({0,0.35}, 0.21, 0.25, 0), Ellipse({0,0.1}, 0.046, 0.046, 0),
                                                             Ellipse({0,-0.1}, 0.046, 0.046, 0), Ellipse({-0.08,-0.605}, 0.046, 0.023, 0),
                                                             Ellipse({0,-0.605}, 0.023, 0.023, 0), Ellipse({0.06,-0.605}, 0.023, 0.046, 0)};
