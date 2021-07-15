
#ifndef HELPER_H
#define HELPER_H

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include <vector>
#include <string>

// numpy-esque list functions
std::vector<double> linspace(double start_in, double end_in, int num_in);
void cartesian_product( std::vector<std::vector<double> >*res, std::vector<double> x, std::vector<double> y);

// debugging
void save_coords(std::vector<std::vector<double> >* p, Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver, std::string c_name, std::string e_name);
void save_coords(std::vector<std::vector<double> >* p, std::vector<double> lower_ellipse, std::vector<double> upper_ellipse, std::string name);

// error message printing
void print_error_message(std::string test_name, double expected_value, double actual_value);
void print_error_message(std::string test_name, std::string expected_value, std::string actual_value);

#endif



