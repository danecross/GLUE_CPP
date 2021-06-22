
#include <vector>

#ifndef ELLIPSE_FUNTIONS_H
#define ELLIPSE_FUNCTIONS_H

std::vector<double> polar_to_cartesian(double r, double theta);
double eccentricity (double a, double b);
double r(double theta, double a, double b, double alpha=0);

// TODO: add center as default argument
std::vector<std::vector<double> > ellipse_cut ( std::vector<std::vector<double> > p, std::vector<double> lower_ellipse, std::vector<double> upper_ellipse);



#endif

