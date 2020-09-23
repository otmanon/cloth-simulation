#pragma once
#include <eigen\Dense>
#include <eigen\Sparse>

struct CGSolver
{
	Eigen::SparseMatrix<double> P;	//preconditioner
	Eigen::SparseMatrix<double> Pinv;	//preconditioner
	Eigen::SparseMatrix<double> A;	//lhs 
	Eigen::VectorXd b;				//rhs
	Eigen::VectorXd x;				//output, what we are solving for
	Eigen::SparseMatrix<double> Si; //multiplies vectors to filter out constraints

	Eigen::VectorXd z;				//desired motion for particles that are constrained... should be zero for constrained particles if we dont want motion Right?

	float epsilon = 1e-6;			//threshold
	Eigen::VectorXd r;				//residual

	/*
	Solves system using modified conjugate gradient scheme outlined in Baraff/Witkin 5.3
	*/
	Eigen::VectorXd solve(Eigen::SparseMatrix<double>& lhs, Eigen::VectorXd& rhs)
	{
		A = lhs;
		b = rhs;

		setPreconditionerInfo();
		x.resize(b.rows(), 1);
		x.setZero(); //set vertices that cant move to zero
		float  delta0 = filter(b).dot(P * filter(b));
		r = filter(b - A * x);
		Eigen::VectorXd c = filter(Pinv *  r);
		float deltaNew = r.transpose()*c;

		Eigen::VectorXd q;
		Eigen::VectorXd s;
		while (deltaNew > epsilon*epsilon*delta0)
		{
			q = filter(A*c);
			float alpha = deltaNew / (c.transpose()*q);
			x += alpha * c;
			r = r - alpha*q;
			s = Pinv * r;
			float deltaOld = deltaNew;
			deltaNew = r.transpose() * s;
			c = filter(s + c * deltaNew / deltaOld);

		}

		//x vector will hold the final solution
		return x;
	}

	//Need list of vertex indices to initialize Si.
	void setUpSi(Eigen::MatrixXd V)
	{

		//Clamp Top most vertices
		Eigen::VectorXd indices(V.rows(), 1);
		indices.setZero();
		int maxZ = V.col(2).maxCoeff();
		int index = 0;
		for (int i = 0; i < V.rows(); i++)
		{
			if (V(i, 2) == maxZ)
			{
				indices(index) = i;
				index++;
			}
		}
		indices.conservativeResize(index, 1);

		Si.resize(V.rows()*3, V.rows()*3);
		Si.setIdentity();
		for (int i = 0; i < indices.rows(); i++)
		{
			Si.coeffRef(3 * indices(i) + 0, 3 * indices(i) + 0) = 0;
			Si.coeffRef(3 * indices(i) + 1, 3 * indices(i) + 1) = 0;
			Si.coeffRef(3 * indices(i) + 2, 3 * indices(i) + 2) = 0;
		}
	}

	void setPreconditionerInfo()
	{
		P.resize(A.rows(), A.cols());
		Pinv.resize(A.rows(), A.cols());
		for (int i = 0; i < A.outerSize(); ++i)
		{
			for (Eigen::SparseMatrix<double>::InnerIterator it(A, i); it; ++it)
			{
				if (it.row() == it.col())
				{
					P.coeffRef(it.row(), it.col()) = 1/it.value();
					Pinv.coeffRef(it.row(), it.col()) = it.value();

				}
			}
		}
	}

	Eigen::VectorXd filter(Eigen::VectorXd a)
	{
		return Si * a;
	}
};
