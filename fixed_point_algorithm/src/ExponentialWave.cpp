#include <cmath>
#include "ExponentialWave.h"


double ExponentialWave::f(double x) const{
    if(fabs(x)>=1.0) {
        return 0.0;
    }

    return exp(1.0/(pow(x, 2) - 1.0));
}

const double ExponentialWave::w_l1 = 0.44399;