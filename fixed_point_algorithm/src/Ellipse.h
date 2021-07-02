#ifndef MAIN_ELLIPSE_H
#define MAIN_ELLIPSE_H

#include <utility>

class Ellipse{
public:
    Ellipse(std::pair<double, double> center, double major_axis, double minor_axis, double tilt);
    bool contains_point(std::pair<double, double> pt) const;

private:
    std::pair<double, double> transform_point(std::pair<double, double> pt) const;
    std::pair<double, double> center;
    double major_axis;
    double minor_axis;
    double tilt;
    double cos_tilt;
    double sin_tilt;
};




#endif //MAIN_ELLIPSE_H
