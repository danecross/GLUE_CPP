
//#define DEBUG

#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <math.h>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>

#include "axis_convergence.h"
#include "helpers.h"


Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> iterate(std::vector<std::vector<double> >* p, int maxiter, double converge_radius){


	Eigen::VectorXd M(2), M_last(2); Eigen::MatrixXd evecs, evecs_new; std::vector<double> q, q_last; 
	M << 1.0,0.1 ; M_last << 1.0, 1.0; evecs = Eigen::MatrixXd::Identity(2,2); 
	q = q_calc(p, M_last);

	// iterate
	int i = 0; Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver; Eigen::VectorXd x(2); x << 1,0;
	while ( abs(M(0)/M(1) - M_last(0)/M_last(1)) > converge_radius && i < maxiter ){
	
		M_last = M; 
		solver = M_calc( p, q );
		M = solver.eigenvalues().real(); evecs_new = solver.eigenvectors().real(); 
		evecs = rotate_evecs(evecs.col(0), evecs.col(1), evecs_new.col(0), evecs_new.col(1));

		q_last = q;
		p = rotate_coords(evecs_new, p);
		q = q_calc(p, M);
		i++;
#ifdef DEBUG
		save_coords(p, solver, std::to_string(i));
#endif
	}

#ifdef DEBUG
	std::cout << "exited on iteration number: " << i << std::endl;
	std::cout << "radius of convergence: " << abs(M(0)/M(1) - M_last(0)/M_last(1)) << std::endl; 
#endif

	return solver;

}

std::vector<std::vector<double> >* rotate_coords( Eigen::MatrixXd evecs, std::vector<std::vector<double> >* p){

	Eigen::VectorXd x(2); x << 1,0;
	auto v1 = evecs.col(0) ; auto v2 = evecs.col(1);
	double alpha = std::min(acos(x.dot(v1)/(x.norm()*v1.norm())), acos(x.dot(v2)/(x.norm()*v2.norm())));

	Eigen::VectorXd v(2); std::vector<double> coord;
	for ( int i = 0 ; i < p->size() ; ++i ){
		
		coord = (*p)[i];
		v << coord[0], coord[1];
		v = rotate_vector(v, alpha);
		(*p)[i][0] = v(0) ; (*p)[i][1] = v(1);
	}

	return p;

}

Eigen::MatrixXd rotate_evecs(Eigen::VectorXd x1, Eigen::VectorXd y1, Eigen::VectorXd x2, Eigen::VectorXd y2){

	// get angles between evecs
	double alpha = std::min(acos(x1.dot(x2)/(x1.norm()*x2.norm())), acos(x1.dot(y2)/(x1.norm()*y2.norm())));

	Eigen::MatrixXd evecs_f(2,2);
	auto x1p = rotate_vector(x1, alpha); auto x2p = rotate_vector(x2, alpha);
	evecs_f(0,0) = x1p(0); evecs_f(0,1) = x1p(1);
	evecs_f(1,0) = x2p(0); evecs_f(1,1) = x2p(1);	

	return evecs_f;
}

Eigen::VectorXd rotate_vector(Eigen::VectorXd v, double alpha){

        Eigen::VectorXd rotated_vector;
	Eigen::Rotation2D<double> rotation(alpha);
           
	rotated_vector = rotation*v;

        return rotated_vector;

}

std::vector<double> q_calc( std::vector< std::vector<double> >* p, Eigen::VectorXd M){

	std::vector<double> q;
	double ab = sqrt(M[0]/M[1]);
	
	for ( auto coords : *p ){
		q.push_back(sqrt(pow(coords[0], 2) + pow(coords[1]/ab, 2)));
	}

	return q;
}

Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> M_calc( std::vector<std::vector<double> >* p, std::vector<double> q){

	Eigen::MatrixXd M(2,2); M << 0,0,0,0; 
	
	for ( int i = 0 ; i < q.size() ; ++i){
	
		M(0,0) += (*p)[i][0]*(*p)[i][0]/pow(q[i], 2);
		M(1,0) += (*p)[i][1]*(*p)[i][0]/pow(q[i], 2);
		M(0,1) += (*p)[i][0]*(*p)[i][1]/pow(q[i], 2);
		M(1,1) += (*p)[i][1]*(*p)[i][1]/pow(q[i], 2);

	}

	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(M);

	return es;


}









