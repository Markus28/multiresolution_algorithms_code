#ifndef BIOMIMETICS_FUZZYCIRCLEABSORPTION_H
#define BIOMIMETICS_FUZZYCIRCLEABSORPTION_H

#include <deal.II/base/point.h>
#include "Absorption.h"

class FuzzyCircleAbsorption: public Absorption{
public:
    virtual double absorption(const dealii::Point<2>& p) const override;
};


#endif //BIOMIMETICS_FUZZYCIRCLEABSORPTION_H
