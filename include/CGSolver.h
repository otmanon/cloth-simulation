#pragma once
#include <eigen\Dense>
#include <eigen\Sparse>

struct CGSolver
{
private:
	Eigen::SparseMatrix<double> P;	//preconditioner
	Eigen::SparseMatrix<double> Pinv;	//preconditioner
	Eigen::SparseMatrix<double> A;	//lhs 
	Eigen::VectorXd b;				//rhs
	Eigen::VectorXd x;				//output, what we are solving for
	Eigen::SparseMatrix<double> Si; //multiplies vectors to filter out constraints

	Eigen::VectorXd z;				//desired motion for particles that are constrained... should be zero for constrained particles if we dont want motion Right?

	float epsilon = 1e-6;			//threshold
	Eigen::VectorXd r;				//residual

public:
	/*
	Solves system using modified conjugate gradient scheme outlined in Baraff/Witkin 5.3
	*/
	Eigen::VectorXd solve(Eigen::SparseMatrix<double>& lhs, Eigen::VectorXd& rhs);

	//Need list of vertex indices to initialize Si.
	void setUpSi(Eigen::MatrixXd V);

	void setPreconditionerInfo();
	

	Eigen::VectorXd filter(Eigen::VectorXd a);
};
