
#include <chrono>
using namespace std::chrono;

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <algorithm>
#include <vector>
#include <math.h>
#include <iostream>

#include "constant_density_shells.h"
#include "helpers.h"
#include "ellipse_functions.h"
#include "axis_convergence.h"

std::vector<double> shell_fit(std::vector<std::vector<double> >* p, std::vector<double> lower_ellipse, double sep, int guesses_per_branch, int numiter, double alpha_res){
	std::vector<std::vector<double> > guesses ; 
	std::vector<double> best_guess {lower_ellipse[0]+sep, (lower_ellipse[0]+sep)*(lower_ellipse[1]/lower_ellipse[0]), lower_ellipse[2]};

	double b_res = sep; 
	for ( int i = 0 ; i < numiter ; ++i ){
		std::cout << "iteration " << std::to_string(i) << std::endl
		   	  << "\tb resolution: \t\t" << std::to_string(b_res) << std::endl
			  << "\talpha_resolution: \t" << std::to_string(alpha_res)
			  << std::endl << std::endl;

		generate_guesses(&guesses, best_guess, lower_ellipse, sep, b_res, alpha_res, guesses_per_branch);
		best_guess = test_guesses(p, lower_ellipse, sep, guesses);
		
		//update resolution
		b_res = 3*b_res/double(guesses_per_branch) ; alpha_res = 3*alpha_res/double(guesses_per_branch);
		guesses.clear();
	}

	std::vector<double> best_upper_ellipse {lower_ellipse[0]+sep, best_guess[0], best_guess[1]};
	
	std::cout << "best fit: a = " << lower_ellipse[0]+sep << ", b = " << best_guess[0] << ", c = " << best_guess[1]  << std::endl;
	
	return best_upper_ellipse;
}

void generate_guesses( std::vector<std::vector<double> >* guesses, std::vector<double> best_guess, std::vector<double> lower_ellipse, double sep, double b_res, double alpha_res, int guesses_per_branch){

	auto b_grid = linspace(best_guess[0]-b_res, best_guess[0]+b_res, guesses_per_branch);
	auto alpha_grid = linspace(best_guess[1]-alpha_res, best_guess[1]+alpha_res, guesses_per_branch);
	
	cartesian_product(guesses, b_grid, alpha_grid);
}

std::vector<double> test_guesses (std::vector<std::vector<double> >* p, std::vector<double> lower_ellipse, double sep, std::vector<std::vector<double> > guesses){

	std::vector<double> data_ratios = get_data_ratios(p, lower_ellipse, sep, guesses);
	std::vector<double> const_ratios = get_const_ratios(p->size(), lower_ellipse, sep, guesses);

	double diff, min_diff; int best_i; min_diff=INFINITY;
	for ( int i = 0 ; i < data_ratios.size() ; ++i){
		diff = abs(data_ratios[i] - const_ratios[i]);
		if ( diff < min_diff ){
			min_diff = diff;
			best_i = i;
		}
	}

	return guesses[best_i];
}

std::vector<double> get_data_ratios(std::vector<std::vector<double> >* p, std::vector<double> lower_ellipse, double sep, std::vector<std::vector<double> > guesses){

	std::vector<double> ratios;
	std::vector<std::vector<double> > cut; std::vector<double> upper_ellipse;
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver; Eigen::VectorXd evals;
	double M1(0), M2(0);

	int i = 0;
	for ( auto guess : guesses ){
		
		cut.clear();
		auto start = high_resolution_clock::now();

		upper_ellipse = {lower_ellipse[0]+sep, guess[0], guess[1]};
		ellipse_cut(&cut, p, lower_ellipse, upper_ellipse);
		solver = iterate(&cut);
		evals = solver.eigenvalues();
        	M1 = evals(0); M2 = evals(1);
		ratios.push_back(sqrt(M1/M2));
		
		if ( i%100 == 0 ){
			auto stop = high_resolution_clock::now();
        		auto duration = duration_cast<milliseconds>(stop - start);
			std::cout << "\ttest data iteration " << i << " of " << guesses.size() << " done in " << duration.count() << " milliseconds" << std::endl;
		}
		i++; 

	}

	return ratios;

}


std::vector<double> get_const_ratios(int num_particles, std::vector<double> lower_ellipse, double sep, std::vector<std::vector<double> > guesses){

	std::vector<double> ratios;
	std::vector<double> upper_ellipse; std::vector<std::vector<double> > shell;
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver; Eigen::VectorXd evals;
        double M1(0), M2(0);
	
	int i = 0;
	for ( auto guess : guesses ){
	
		auto start = high_resolution_clock::now();
		shell.clear();
		
		upper_ellipse = {lower_ellipse[0]+sep, guess[0], guess[1]};
		create_const_shell(&shell, lower_ellipse, upper_ellipse, num_particles);
		if (shell.size() == 0){ratios.push_back(INFINITY); continue;}
		solver = iterate(&shell);
		
                evals = solver.eigenvalues();
                M1 = evals(0); M2 = evals(1);
                ratios.push_back(sqrt(M1/M2));

		if ( i%100 == 0 ){
			auto stop = high_resolution_clock::now();
        		auto duration = duration_cast<millliseconds>(stop - start);
			std::cout << "\tconst shell iteration " << i << " of " << guesses.size() << " done in " << duration.count() << " milliseconds" << std::endl;
		}
		i++;

	}

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


