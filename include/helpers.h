
#ifndef HELPER_H
#define HELPER_H

#include <vector>
#include <string>

// numpy-esque list functions
std::vector<double> linspace(double start_in, double end_in, int num_in);
std::vector<std::vector<double> > cartesian_product( std::vector<double> x, std::vector<double> y);

// error message printing
void print_error_message(std::string test_name, double expected_value, double actual_value);
void print_error_message(std::string test_name, std::string expected_value, std::string actual_value);

#endif



