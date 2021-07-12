
#include <vector>
#include <iostream>
#include <random>
#include <map>
#include <fstream>

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include "ellipse_functions.h"
#include "axis_convergence.h"
#include "constant_density_shells.h"
#include "helpers.h"


int test_get_stars(){

	std::vector<std::vector<double> > p = {{0, 1}, {0, 2}, {0, 3}, {0, 4}};
        std::vector<std::vector<double> >* p_new;
	auto p_ptr = &p;

        // test null selection
        p_new = get_stars(p_ptr);
        if (p_new->size() != 4 ){print_error_message("get_stars null selection", 4, p_new->size());}

        // test shell selection
        p_new = get_stars(p_ptr, 0.0, 0.4);
        if (p_new->size() != 2) { print_error_message("get_stars, shell selection (basic)", 2, p_new->size());}

        // test shell selection with half-mass-radius
        p_new = get_stars(p_ptr, 0.0, 0.5, 10);
        if (p_new->size() != 4){print_error_message("get_stars, shell selection w/ hmr", 4, p_new->size()); return 0;}

	return 1;

}

int test_rotate_vector(){

	Eigen::VectorXd rotated_vector(2), v(2);
        double alpha = M_PI/3; v << 1, 0;
        rotated_vector = rotate_vector(v, alpha);
        double exp_x = 0.5; double rec_x = rotated_vector[0];
        if (float(rec_x) != float(exp_x)) {print_error_message("vector rotation", exp_x, rec_x); return 0;}

	return 1;
}

int test_rotate_evecs(){

	Eigen::MatrixXd M1(2,2); M1 << 1,1,1,1;
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es1(M1);
        auto evals1 = es1.eigenvalues();
        auto evecs1 = es1.eigenvectors();

        Eigen::MatrixXd M2(2,2); M2 << sqrt(2)/2, sqrt(2)/2, -sqrt(2)/2, sqrt(2)/2;
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es2(M2);
        auto evals2 = es2.eigenvalues();
        auto evecs2 = es2.eigenvectors();

        auto evecs_rot = rotate_evecs(evecs1.col(0).real(), evecs1.col(1).real(),
                                      evecs2.col(0).real(), evecs2.col(1).real());
        Eigen::MatrixXd exp_vecs(2,2) ; exp_vecs << 0,1,0.5,0.5;
        std::vector<double> diffs = {evecs_rot(0)-exp_vecs(0), evecs_rot(1)-exp_vecs(1)};
        if (diffs[0] > 1e-5 || diffs[1] > 1e-5){
                std::stringstream exp ; std::stringstream got;
                exp << exp_vecs; got << evecs_rot;
                print_error_message("eigenvector rotation", exp.str(), got.str());
		return 0;
        }

	return 1;

}

int test_rotate_coords(){

	std::vector<std::vector<double> > p = {{0, 1}, {0, 2}, {0, 3}, {0, 4}};
        Eigen::MatrixXd m(2,2) ; m << 1,0,0,1;
        std::vector<std::vector<double> >* p_rot = rotate_coords(m, &p);

        std::vector<std::vector<double> > p_exp = {{-1,0},{-2,0},{-3,0},{-4,0}};

        if (p_exp[0][0]-(*p_rot)[0][0] > 1e-5 || p_exp[2][1]-(*p_rot)[2][1] > 1e-5){
                std::cout << (p_exp[0][0]-(*p_rot)[0][0]) << std::endl;
                print_error_message("coordinate rotation", "??", "??");
		return 0;
        }
	
	return 1;
}

int test_iterate(){

	// generate random distribution
	double a = 10.0 ; double b = 5.0 ;
	std::vector<double> low, up; low = {0.1, 0.1} ; up = {a, b};
	std::vector< std::vector<double> >* p = create_const_shell(low, up, 5000);

	// run iterate
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver = iterate(p);
	auto evals = solver.eigenvalues();
	double M1(evals(0)), M2(evals(1));

	save_coords(p, solver, "final");

	// expected ratio
	if (abs(sqrt(M1/M2)-b/a) > 0.01){
		print_error_message("iterate test", b/a, sqrt(M1/M2));
	}

	return 1;
}

int main() {
	

	int res = 0;
	res += test_get_stars();
	res += test_rotate_vector();
	res += test_rotate_evecs();
	res += test_rotate_coords();
	res += test_iterate();
	
	if ( res != 5 ){
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;

}

