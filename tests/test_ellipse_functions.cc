#include <string>
#include <vector>
#include <iostream>
#include <math.h>

#define _USE_MATH_DEFINES

#include "ellipse_functions.h"
#include "helpers.h"

int main(){

	/* test conversion */
	double rad, theta; std::vector<double> conv_result;
	
	rad = 10 ; theta = M_PI/3;
	conv_result = polar_to_cartesian(rad, theta);
	if (int(conv_result[0]) != 5) {print_error_message("polar_to_cartesian conversion", 5, conv_result[0]);}

	/* test eccentricity */
	double a, b, ecc, exp_ecc;

	a = 3 ; b = 5 ; ecc = eccentricity(a, b); exp_ecc = sqrt(1-pow((a/b),2));
	if (ecc != exp_ecc){
		print_error_message("eccentricity calculation", exp_ecc, ecc);
	}

	/* test radius calculation */
	double exp_r, r_calc;

	// test circle
	theta = M_PI/4 ; a = 3 ; b = 3 ; exp_r = 3;
	r_calc = r(theta, a, b);
	if (r_calc != exp_r){print_error_message("radius calculation: circle", exp_r, r_calc);}

	// test ellipse
	theta = 0 ; a = 4 ; b = 10 ; exp_r = a;
	r_calc = r(theta, a, b);
        if (r_calc != exp_r){print_error_message("radius calculation: general", exp_r, r_calc);}

	/* test ellipse cut */
	std::vector<std::vector<double> > pos, cut_result; std::vector<double> lower_ellipse, upper_ellipse;	
	pos = {{0, 0}, {0, 1}, {0, 2}, {0, 3}, {0, 4}}; lower_ellipse = {1, 0.5}; upper_ellipse = {3, 5};

	cut_result = ellipse_cut(pos, lower_ellipse, upper_ellipse);
	if(cut_result.size() != 4){print_error_message("ellipse cut", 4, cut_result.size());}

	return 0;

}





