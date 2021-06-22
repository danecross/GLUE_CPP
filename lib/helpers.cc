
#include <string>
#include <iostream>


void print_error_message(std::string test_name, int expected_value, int actual_value){

        std::cout << "WARNING: " << test_name << " FAILED" << std::endl;
        std::cout << "\tEXPECTED: " << expected_value << std::endl;
        std::cout << "\tRETURNED: " << actual_value << std::endl;

}



