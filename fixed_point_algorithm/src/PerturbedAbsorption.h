#ifndef MAIN_PERTURBEDABSORPTION_H
#define MAIN_PERTURBEDABSORPTION_H

#include "Absorption.h"


template<class AbsorptionPolicy, class WavePolicy>
class PerturbedAbsorption: public AbsorptionPolicy, private WavePolicy{
public:
    PerturbedAbsorption(const dealii::Point<2>& center, double c, double eta);
    void set_time(double t);
    void set_center(const dealii::Point<2>& center);
    virtual double absorption(const dealii::Point<2>& p) const override;

private:
    using WavePolicy::f;

    dealii::Point<2> center;
    double T;
    double c;
    double eta;
};


template<class AbsorptionPolicy, class WavePolicy>
PerturbedAbsorption<AbsorptionPolicy, WavePolicy>::PerturbedAbsorption(const dealii::Point<2>& center, double c, double eta)
        :T(0), center(center), c(c),eta(eta) {}

template<class AbsorptionPolicy, class WavePolicy>
void PerturbedAbsorption<AbsorptionPolicy, WavePolicy>::set_time(double t) {T=t;}

template<class AbsorptionPolicy, class WavePolicy>
void PerturbedAbsorption<AbsorptionPolicy, WavePolicy>::set_center(const dealii::Point<2>& ct) {
    center=ct;
}

template<class AbsorptionPolicy, class WavePolicy>
double PerturbedAbsorption<AbsorptionPolicy, WavePolicy>::absorption(const dealii::Point<2>& p) const{
    return AbsorptionPolicy::absorption(p-eta/(c*T)*f((p.distance(center)-c*T)/eta)*(p-center)/p.distance(center));
}


#endif //MAIN_PERTURBEDABSORPTION_H
