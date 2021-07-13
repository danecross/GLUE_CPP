
#include <Eigen/Dense>
#include <Eigen/Geometry>

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>

std::vector<double> linspace(double start, double end, int num)
{

  std::vector<double> linspaced;

  if (num == 0) { return linspaced; }
  if (num == 1)
    {
      linspaced.push_back(start);
      return linspaced;
    }

  double delta = (end - start) / (num - 1);

  for(int i=0; i < num-1; ++i)
    {
      linspaced.push_back(start + delta * i);
    }
  linspaced.push_back(end); // I want to ensure that start and end
                            // are exactly the same as the input
  return linspaced;
}


void cartesian_product( std::vector<std::vector<double> >* res, std::vector<double> x, std::vector<double> y){

	res->clear();

        std::vector<std::vector<double> > p ; std::vector<double> coord;
        for ( int i = 0 ; i < x.size() ; ++i){
                for (int j = 0 ; j < y.size() ; ++j){
                        coord.push_back(x[i]) ; coord.push_back(y[j]);
                        p.push_back(coord);
                        coord.clear();
                }
        }

	for ( auto c : p){
		res->push_back(c);
	}

}


void save_coords(std::vector<std::vector<double> >* p, Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver, std::string name){

	std::ofstream saveCoords("basic_data_and_plots/data/coords_"+name+".csv");
        for ( auto coord : *p ){
                saveCoords << coord[0] << ", " << coord[1] << std::endl;
        }   
        saveCoords.close();

        std::ofstream saveEig("basic_data_and_plots/data/eigenresults_"+name+".csv");
        saveEig << "eigenvalues: " << solver.eigenvalues() << std::endl;
        saveEig << "eigenvectors: " << solver.eigenvectors() << std::endl;
        saveEig.close();

}

void save_coords(std::vector<std::vector<double> >* p, std::vector<double> lower_ellipse, std::vector<double> upper_ellipse, std::string name){

	std::ofstream saveCoords("cds_data_and_plots/data/coords_"+name+".csv");
        for ( auto coord : *p ){
                saveCoords << coord[0] << ", " << coord[1] << std::endl;
        }
        saveCoords.close();

	std::ofstream saveEllipses("cds_data_and_plots/data/ellipses_"+name+".csv");
	saveEllipses << lower_ellipse[0] << " " << lower_ellipse[1] << " " << lower_ellipse[2] << std::endl;
	saveEllipses << upper_ellipse[0] << " " << upper_ellipse[1] << " " << upper_ellipse[2] << std::endl;
	saveEllipses.close();	

}

void print_error_message(std::string test_name, double expected_value, double actual_value){

	std::cout << std::endl;
        std::cout << "WARNING: " << test_name << " FAILED" << std::endl;
        std::cout << "\tEXPECTED: " << expected_value << std::endl;
        std::cout << "\tRETURNED: " << actual_value << std::endl;

}

void print_error_message(std::string test_name, std::string expected_value, std::string actual_value){

	std::cout << std::endl;
        std::cout << "WARNING: " << test_name << " FAILED" << std::endl;
        std::cout << "\tEXPECTED: " << expected_value << std::endl;
        std::cout << "\tRETURNED: " << actual_value << std::endl;

}


