#ifndef MAIN_SHEPPLOGANABSORPTION_H
#define MAIN_SHEPPLOGANABSORPTION_H

#include <vector>
#include <deal.II/base/point.h>
#include "Ellipse.h"
#include "Absorption.h"

class SheppLoganAbsorption: public Absorption{
public:
    SheppLoganAbsorption() = default;
    virtual double absorption(const dealii::Point<2>& p) const override;

private:
    static const std::vector<Ellipse> ellipses;
    static const std::vector<double> weights;
};


#endif