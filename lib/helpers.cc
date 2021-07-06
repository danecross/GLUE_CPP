
#include <vector>
#include <string>
#include <iostream>


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


std::vector<std::vector<double> > cartesian_product( std::vector<double> x, std::vector<double> y){

        std::vector<std::vector<double> > p ; std::vector<double> coord;
        for ( int i = 0 ; i < x.size() ; ++i){
                for (int j = 0 ; j < y.size() ; ++j){
                        coord.push_back(x[i]) ; coord.push_back(y[j]);
                        p.push_back(coord);
                        coord.clear();
                }
        }

        return p;

}




void print_error_message(std::string test_name, double expected_value, double actual_value){

        std::cout << "WARNING: " << test_name << " FAILED" << std::endl;
        std::cout << "\tEXPECTED: " << expected_value << std::endl;
        std::cout << "\tRETURNED: " << actual_value << std::endl;

}

void print_error_message(std::string test_name, std::string expected_value, std::string actual_value){

        std::cout << "WARNING: " << test_name << " FAILED" << std::endl;
        std::cout << "\tEXPECTED: " << expected_value << std::endl;
        std::cout << "\tRETURNED: " << actual_value << std::endl;

}


