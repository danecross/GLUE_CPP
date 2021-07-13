

#ifndef AXIS_CONVERGENCE_H
#define AXIS_CONVERGENCE_H

#include <vector>
#include <map>
#include <string>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>


Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> iterate(std::vector<std::vector<double> >* p, int maxiter=40, double converge_radius=10e-4);

std::vector<double> q_calc( std::vector< std::vector<double> >* p, Eigen::VectorXd M);
Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> M_calc( std::vector<std::vector<double> >* p, std::vector<double> q);
Eigen::VectorXd rotate_vector(Eigen::VectorXd v, double alpha);
Eigen::MatrixXd rotate_evecs(Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd);
std::vector<std::vector<double> >* rotate_coords( Eigen::MatrixXd evecs, std::vector<std::vector<double> >* p);



#endif

