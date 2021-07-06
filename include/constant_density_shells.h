

#ifndef CDS_H
#define CDS_H

#include <vector>

std::vector<std::vector<double> > create_const_box ( double w, double rho, double offset=0);
std::vector<std::vector<double> > create_const_shell (std::vector<double> lower_ellipse, std::vector<double> upper_ellipse, int num_stars);

#endif



