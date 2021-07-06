
#include <algorithm>
#include <vector>
#include <math.h>

#include "constant_density_shells.h"
#include "helpers.h"
#include "ellipse_functions.h"


std::vector<std::vector<double> > create_const_shell (std::vector<double> lower_ellipse, std::vector<double> upper_ellipse, int num_stars){

	double A, rho; std::vector<std::vector<double> > box, shell;

	A = M_PI * (upper_ellipse[0] * upper_ellipse[1] - lower_ellipse[0] * lower_ellipse[1]);
	rho = (double)(num_stars)/A; // <-- number density (~mass density for equal mass system)
	box = create_const_box(std::max(2*upper_ellipse[0], 2*upper_ellipse[1]), rho);

	return ellipse_cut(box, lower_ellipse, upper_ellipse);
}


std::vector<std::vector<double> > create_const_box ( double w, double rho, double offset){
	
	float spacing = w/((int)(w*sqrt(rho)));
	offset = offset*spacing;

	std::vector<double> x = linspace(-w/2-offset, w/2-offset+spacing, (int)(w*sqrt(rho)));
	std::vector<double> y = linspace(-w/2-offset, w/2-offset+spacing, (int)(w*sqrt(rho)));

	return cartesian_product(x, y);
}


