#include "FuzzyCircleAbsorption.h"


double FuzzyCircleAbsorption::absorption(const dealii::Point<2>& p) const{
    double r = pow(pow(p(0), 2)+pow(p(1), 2), 0.5);

    if(r <= 0.25) {
        return 2;
    }

    else if(r<=0.3){
        double diff = (r-0.25)*20.0;
        return 2.0-3*pow(diff, 2)+2*pow(diff, 3);
    }

    return 1;
}