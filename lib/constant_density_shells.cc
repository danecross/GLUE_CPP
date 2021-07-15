
//#define DEBUG

#include <chrono>
using namespace std::chrono;

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <algorithm>
#include <vector>
#include <math.h>
#include <iostream>
#include <thread>
#include <utility>

#include "constant_density_shells.h"
#include "helpers.h"
#include "ellipse_functions.h"
#include "axis_convergence.h"

std::vector<double> shell_fit(std::vector<std::vector<double> >* p, std::vector<double> lower_ellipse, double sep, int num_cpus, int guesses_per_branch, int numiter, double alpha_res){
	std::vector<std::vector<double> > guesses ; 
	std::vector<double> best_guess {lower_ellipse[0]+sep, (lower_ellipse[0]+sep)*(lower_ellipse[1]/lower_ellipse[0]), lower_ellipse[2]};

	double b_res = sep; 
	for ( int i = 0 ; i < numiter ; ++i ){
#ifdef DEBUG
		std::cout << "iteration " << std::to_string(i) << std::endl
		   	  << "\tb resolution: \t\t" << std::to_string(b_res) << std::endl
			  << "\talpha_resolution: \t" << std::to_string(alpha_res)
			  << std::endl << std::endl;
#endif
		generate_guesses(&guesses, best_guess, lower_ellipse, sep, b_res, alpha_res, guesses_per_branch);
		best_guess = test_guesses(p, lower_ellipse, sep, guesses, num_cpus);
		
		//update resolution
		b_res = 3*b_res/double(guesses_per_branch) ; alpha_res = 3*alpha_res/double(guesses_per_branch);
		guesses.clear();
	}

	std::vector<double> best_upper_ellipse {lower_ellipse[0]+sep, best_guess[0], best_guess[1]};

#ifdef DEBUG
	std::cout << "best fit: a = " << lower_ellipse[0]+sep << ", b = " << best_guess[0] << ", c = " << best_guess[1]  << std::endl;
#endif
	return best_upper_ellipse;
}

void generate_guesses( std::vector<std::vector<double> >* guesses, std::vector<double> best_guess, std::vector<double> lower_ellipse, double sep, double b_res, double alpha_res, int guesses_per_branch){

	auto b_grid = linspace(best_guess[0]-b_res, best_guess[0]+b_res, guesses_per_branch);
	auto alpha_grid = linspace(best_guess[1]-alpha_res, best_guess[1]+alpha_res, guesses_per_branch);
	
	cartesian_product(guesses, b_grid, alpha_grid);
}

std::vector<double> test_guesses (std::vector<std::vector<double> >* p, std::vector<double> lower_ellipse, double sep, std::vector<std::vector<double> > guesses, int num_cpus){

	auto const_ratios = get_const_ratios(p->size(), lower_ellipse, sep, guesses, num_cpus);
	auto data_ratios = get_data_ratios(p, lower_ellipse, sep, guesses, num_cpus);

	sort(data_ratios.begin(), data_ratios.end());
	sort(const_ratios.begin(), const_ratios.end());

	double diff, min_diff; int best_i; min_diff=INFINITY;
	for ( int i = 0 ; i < data_ratios.size() ; ++i){
		diff = abs(data_ratios[i].second - const_ratios[i].second);
		if ( diff < min_diff ){
			min_diff = diff;
			best_i = i;
		}
	}

	return data_ratios[best_i].first;
}

void get_ratio(std::vector<std::vector<double> > cut, std::vector<double>* upper_ellipse, std::vector<std::pair< std::vector<double>, double> > *ratios){
	
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver; Eigen::VectorXd evals;
	double M1(0), M2(0);

	solver = iterate(&cut);
	evals = solver.eigenvalues();
	M1 = evals(0); M2 = evals(1);

	std::pair< std::vector<double>, double> result;
	result.first = *upper_ellipse;
	result.second = sqrt(M1/M2);

	ratios->push_back(result);
}

std::vector<std::pair< std::vector<double>, double> > get_data_ratios(std::vector<std::vector<double> >* p, std::vector<double> lower_ellipse, double sep, std::vector<std::vector<double> > guesses, int num_cpus){

	std::vector<std::pair< std::vector<double>, double> > ratios;
	std::vector<std::thread> workers; 
	std::vector<std::vector<std::vector<double> > > cuts;

	int i = -1; int numiter = 0;
	for ( auto guess : guesses ){
		
		std::vector<double> upper_ellipse = {lower_ellipse[0]+sep, guess[0], guess[1]};
		std::vector<std::vector<double> > cut;
		ellipse_cut(&cut, p, lower_ellipse, upper_ellipse);
		if (cut.size() == 0){continue;}
		cuts.push_back(cut); i++;

		workers.push_back(std::thread(get_ratio, cuts[i], &upper_ellipse, &ratios));
		if ( workers.size() >= num_cpus ){
			for ( std::thread &t : workers){ 
				t.join(); numiter++;
#ifdef DEBUG
				if ( numiter%100 == 0 ){
					std::cout << "\ttest data iteration " << numiter << " of " << guesses.size() << " done" << std::endl;
				}
#endif
			}
			workers.clear();
		}
	}

	for ( std::thread &t : workers){ t.join(); }

	return ratios;

}


std::vector<std::pair< std::vector<double>, double> > get_const_ratios(int num_particles, std::vector<double> lower_ellipse, double sep, std::vector<std::vector<double> > guesses, int num_cpus){

	std::vector<std::pair< std::vector<double>, double> > ratios;
	std::vector<std::thread> workers; 
	std::vector<std::vector<std::vector<double> > > shells;

	int i = -1; int numiter = 0;
	for ( auto guess : guesses ){
	
		std::vector<double> upper_ellipse = {lower_ellipse[0]+sep, guess[0], guess[1]};
		std::vector<std::vector<double> > shell;
		create_const_shell(&shell, lower_ellipse, upper_ellipse, num_particles);
		if (shell.size() == 0){continue;}
		shells.push_back(shell); i++;

		workers.push_back(std::thread(get_ratio, shells[i], &upper_ellipse, &ratios));
		if ( workers.size() >= num_cpus ){
                        for ( std::thread &t : workers){ 
				t.join(); numiter++;
#ifdef DEBUG
                                if ( numiter%100 == 0 ){
                                        std::cout << "\ttest const iteration " << numiter << " of " << guesses.size() << " done" << std::endl;
                                }
#endif
			}
                        workers.clear();
                } 
	}

	for ( std::thread &t : workers ){ t.join(); }

	return ratios;
}

void create_const_shell (std::vector<std::vector<double> >* res, std::vector<double> lower_ellipse, std::vector<double> upper_ellipse, int num_stars){

	double A, rho; std::vector<std::vector<double> > box;
	res->clear();

	A = M_PI * (upper_ellipse[0] * upper_ellipse[1] - lower_ellipse[0] * lower_ellipse[1]);
	rho = (double)(num_stars)/A; // <-- number density (~mass density for equal mass system)

	if ( A < 0 || rho < 0 ){ return;}
	
	create_const_box(&box, std::max(2*upper_ellipse[0], 2*upper_ellipse[1]), rho);
	ellipse_cut(res, &box, lower_ellipse, upper_ellipse);
	
}


void create_const_box ( std::vector<std::vector<double> >* res, double w, double rho, double offset){
	
	float spacing = w/((int)(w*sqrt(rho)));
	offset = offset*spacing;

	std::vector<double> x = linspace(-w/2-offset, w/2-offset+spacing, (int)(w*sqrt(rho)));
	std::vector<double> y = linspace(-w/2-offset, w/2-offset+spacing, (int)(w*sqrt(rho)));

	cartesian_product(res, x, y);
}


