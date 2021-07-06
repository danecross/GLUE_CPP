
#include <math.h>

#include "ellipse_functions.h"




std::vector<double> polar_to_cartesian(double r, double theta){

	std::vector<double> coords {r*cos(theta), r*sin(theta)};
	return coords;

}

std::vector<double> cartesian_to_polar(double x, double y){
	
	std::vector<double> coords {sqrt(x*x + y*y), atan(y/x)};
	return coords;
}

double eccentricity (double a, double b){ return sqrt(1-pow(a/b, 2));}

double r(double theta, double a, double b, double alpha){return a*b/sqrt(pow((b*cos(theta-alpha)),2) + pow((a*sin(theta-alpha)),2));}


std::vector<std::vector<double> > ellipse_cut ( std::vector<std::vector<double> > p, std::vector<double> lower_ellipse, std::vector<double> upper_ellipse){

	std::vector<std::vector<double> > cut_pos;
	double a_lower, b_lower, alpha_lower, a_upper, b_upper, alpha_upper;
	
	a_lower = lower_ellipse[0] ; b_lower = lower_ellipse[1] ; alpha_lower=0;
	a_upper = upper_ellipse[0] ; b_upper = upper_ellipse[1] ; alpha_upper=0;
	
	if (lower_ellipse.size() == 3){alpha_lower = lower_ellipse[2];}
	if (upper_ellipse.size() == 3){alpha_upper = upper_ellipse[2];}

	double radius, theta;
	for ( auto coord : p ){
		radius = sqrt(pow(coord[0], 2) + pow(coord[1], 2));
		theta = atan(coord[1]/coord[0]);
		if (r(theta, a_lower, b_lower, alpha_lower) < radius
			&& r(theta, a_upper, b_upper, alpha_upper) > radius)
			cut_pos.push_back(coord);
	}

	return cut_pos;

}





