
#include <fstream>
#include <iostream>

#include "ellipse_functions.h"
#include "helpers.h"
#include "constant_density_shells.h"

int test_const_box(){

	std::vector<std::vector<double> > box = create_const_box(15.0, 100.0);

	auto num = *max_element(std::begin(box), std::end(box));
	
	if (abs(num[0]-7.5) < 0.5 && abs(num[1]-7.5) < 0.5){
		std::ofstream saveCoords("cds_data_and_plots/const_box_coords.csv");
		for ( auto coord : box ){
			saveCoords << coord[0] << "," << coord[1] << std::endl;
		}
		saveCoords.close();
		return 1;
	}
	else{
		print_error_message("create constant box", std::to_string(7.5), std::to_string(num[0])+", "+std::to_string(num[1]));
		return 0;
	}
}

int test_const_shell(){
	
	std::vector<std::vector<double> > shell; std::vector<double> lower_ellipse, upper_ellipse;
	int num_stars=1000;

	lower_ellipse = {20, 10};
	upper_ellipse = {30, 20};

	shell = create_const_shell(lower_ellipse, upper_ellipse, num_stars);

	bool EXCEEDS_UPPER_RADIUS = false; bool EXCEEDS_LOWER_RADIUS = false;
	double r_l, theta_l, r_exp_l, theta_exp_l, r_u, theta_u, r_exp_u, theta_exp_u;
	std::vector<double> pol_coords;
	for ( std::vector<double> c_coords : shell ){
		pol_coords = cartesian_to_polar(c_coords[0], c_coords[1]);
		if      ( r(pol_coords[1], lower_ellipse[0], lower_ellipse[1]) > pol_coords[0]){
			EXCEEDS_LOWER_RADIUS = true;
			r_l = pol_coords[0];
			r_exp_l = r(pol_coords[1], lower_ellipse[0], lower_ellipse[1]);
		}
		else if ( r(pol_coords[1], upper_ellipse[0], upper_ellipse[1]) < pol_coords[0]){
			EXCEEDS_UPPER_RADIUS = true;
			r_u = pol_coords[0];
                        r_exp_u = r(pol_coords[1], upper_ellipse[0], upper_ellipse[1]);
		}

		if (EXCEEDS_LOWER_RADIUS && EXCEEDS_UPPER_RADIUS)
			break;
	}

	if (EXCEEDS_LOWER_RADIUS)
		print_error_message("constant shell: lower radius", "above "+std::to_string(r_exp_l), std::to_string(r_l));
	if (EXCEEDS_UPPER_RADIUS)
		print_error_message("constant shell: upper radius", "below "+std::to_string(r_exp_u), std::to_string(r_u));
	
	std::ofstream saveCoords("cds_data_and_plots/const_shell_coords.csv");
        for ( auto coord : shell ){
        	saveCoords << coord[0] << "," << coord[1] << std::endl;
        }
        saveCoords.close();


	if (!EXCEEDS_LOWER_RADIUS && !EXCEEDS_UPPER_RADIUS)
		return 1;
	return 0;
}

int main(){

	int res = test_const_box() + test_const_shell();


	if (res < 1){
		return EXIT_SUCCESS;
	}
	else{
		return EXIT_FAILURE;
	}

}



