#include "BoxAbsorption.h"

double BoxAbsorption::absorption(const dealii::Point<2>& p) const{
    double lower = -0.25;
    double upper = 0.25;

    if(p(0)<=upper && p(0)>=lower && p(1)<=upper && p(1)>=lower) {
        return 3;
    }

    return 1;
}