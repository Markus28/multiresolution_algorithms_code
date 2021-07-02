#ifndef MAIN_ABSORPTION_H
#define MAIN_ABSORPTION_H

#include <deal.II/base/point.h>

class Absorption{
public:
    virtual double absorption(const dealii::Point<2>& p) const = 0;
    virtual ~Absorption(){};
};


#endif //MAIN_ABSORPTION_H
