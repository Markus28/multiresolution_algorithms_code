#ifndef MAIN_BOXABSORPTION_H
#define MAIN_BOXABSORPTION_H

#include <deal.II/base/point.h>
#include "Absorption.h"

class BoxAbsorption: public Absorption{
public:
    virtual double absorption(const dealii::Point<2>& p) const override;
};


#endif //MAIN_BOXABSORPTION_H
