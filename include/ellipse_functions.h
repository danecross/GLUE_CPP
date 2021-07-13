
#include <vector>

#ifndef ELLIPSE_FUNTIONS_H
#define ELLIPSE_FUNCTIONS_H

std::vector<double> polar_to_cartesian(double r, double theta);
std::vector<double> cartesian_to_polar(double x, double y);
double eccentricity (double a, double b);
double r(double theta, double a, double b, double alpha=0);

// TODO: add center as default argument
void ellipse_cut ( std::vector<std::vector<double> >* res, std::vector<std::vector<double> >* p, std::vector<double> lower_ellipse, std::vector<double> upper_ellipse);



#endif


