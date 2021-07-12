

#ifndef CDS_H
#define CDS_H

#include <vector>

std::vector<double> shell_fit(std::vector<std::vector<double> >* p, std::vector<double> lower_ellipse, double sep, int guesses_per_branch=30, int numiter=2, double alpha_resolution=M_PI/3);

std::vector<std::vector<double> > generate_guesses( std::vector<double> best_guess, std::vector<double> lower_ellipse, double sep, double b_res, double alpha_res, int guesses_per_branch);
std::vector<double> test_guesses (std::vector<std::vector<double> >* p, std::vector<double> lower_ellipse, double sep, std::vector<std::vector<double> > guesses);

std::vector<double> get_data_ratios(std::vector<std::vector<double> >* p, std::vector<double> lower_ellipse, double sep, std::vector<std::vector<double> > guesses);
std::vector<double> get_const_ratios(std::vector<double> lower_ellipse, double sep, std::vector<std::vector<double> > guesses);

std::vector<std::vector<double> >* create_const_box ( double w, double rho, double offset=0);
std::vector<std::vector<double> >* create_const_shell (std::vector<double> lower_ellipse, std::vector<double> upper_ellipse, int num_stars);

#endif



