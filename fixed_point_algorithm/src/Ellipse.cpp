#include <cmath>
#include "Ellipse.h"

Ellipse::Ellipse(std::pair<double, double> center, double major_axis, double minor_axis, double tilt): center(center), major_axis(major_axis), minor_axis(minor_axis), tilt(tilt) {
    cos_tilt = cos(-2*M_PI*tilt/360);
    sin_tilt = sin(-2*M_PI*tilt/360);
}

bool Ellipse::contains_point(std::pair<double, double> pt) const{
    std::pair<double, double> transformed_point = transform_point(pt);
    return pow(transformed_point.first/major_axis, 2) + pow(transformed_point.second/minor_axis, 2) <= 1;
}

std::pair<double, double> Ellipse::transform_point(std::pair<double, double> pt) const{
    pt.first -= center.first;
    pt.second -= center.second;
    return std::make_pair(pt.first*cos_tilt - pt.second*sin_tilt, pt.first*sin_tilt + pt.second*cos_tilt);
}