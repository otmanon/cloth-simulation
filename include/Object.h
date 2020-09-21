#pragma once
#include <eigen/Dense>
#include <eigen/Sparse>
#include "Triangle.h"
#include <Eigen/IterativeLinearSolvers>
struct Object
{
	Eigen::MatrixXd V;					//Vertices
	Eigen::MatrixXi F;					//Faces
	Eigen::MatrixXd E;					//Edges

};
struct Cloth : public Object
{
	Eigen::VectorXd energy;				//sum of stretch, shear, and bending energies
	Eigen::MatrixXd f;					// |V| x 3 force matrix. Each element is the force vector of the corresponding vertex
	Eigen::SparseMatrix<double> K;		// stiffness matrix (df/dx)
	Eigen::MatrixXd y;					//displacement vector to avoid mesh tangling
	Eigen::SparseMatrix<double> M;		//Mass matrix
	Eigen::SparseMatrix<double> W;		//inverse mass matrix used for constraint handling
	Eigen::MatrixXd UV;					//Contains uv coordinates for each vertex in cloth
	
	std::vector<Triangle> triangles;	//individual triangles found in cloth. Used to compute locally force and stiffness values

	/*
	Assume cloth is straight and is a square with evenly spaced vertices in both directions
	*/
	void buildUVCoords()
	{
		
		Eigen::Vector3d tlV, trV, blV, brV;		// corner edges of square cloth
		trV << V.col(0).maxCoeff(), V.col(1).maxCoeff(), 0.0f;
		blV << V.col(0).minCoeff(), V.col(1).minCoeff(), 0.0f;

		tlV << V.col(1).minCoeff(), V.col(1).maxCoeff(), 0.0f;
		brV << V.col(0).maxCoeff(), V.col(1).minCoeff(), 0.0f;
		
		float w, h;
		w = trV.x() - blV.x();
		h = trV.y() - blV.y();
		
		Eigen::Vector2d tr_uv, bl_uv;			//topright most corner has uv coordinates (1, 1), bottomleft has uv coords (0, 0)
		tr_uv << 1.0f, 1.0f;
		bl_uv << 0.0f, 0.0f;

		//go through vertices and interpolate UV coords
		UV.resize(V.rows(), 2);
		Eigen::Vector2d uv;
		float x, y, s_x, s_y;					//x, y, and fraction coordinates
		for (int i = 0; i < V.rows(); i++)
		{
			x = V.row(i).x();
			y = V.row(i).y();

			s_x = (x - blV.x()) / w;
			s_y = (y - blV.y()) / h;
			
			uv(0) = 0.0f * (1 - s_x) + 1.0f * s_x;
			uv(1) = 0.0f * (1 - s_y) + 1.0f * s_y;
			UV.row(i) = uv;
		}

		UV *= 11;


	}

	/*
	Initializes values of cloth
	*/
	void initCloth()
	{
		f.resize(V.rows(), 3);
		K.resize(3 * V.rows(), 3 * V.rows());
		M.resize(3 * V.rows(), 3 * V.rows());
		M.setIdentity();

		triangles = std::vector<Triangle>(F.rows());
		Triangle t;
		float massDensity = 1;
		float onethird = 1 / 3;
		float massthird;
		for (int i = 0; i < F.rows(); i++)
		{
			//fill triangle info
			t.init(V, F, UV, i);
			triangles[i] = t;

			massthird = onethird * t.a;

			//Fill mass matrix
			M.coeffRef(t.indexi + 0, t.indexi + 0) += massthird; //add triangle area times 1/3
			M.coeffRef(t.indexi + 1, t.indexi + 1) += massthird;
			M.coeffRef(t.indexi + 2, t.indexi + 2) += massthird;

			M.coeffRef(t.indexj + 0, t.indexj + 0) += massthird; //add triangle area times 1/3
			M.coeffRef(t.indexj + 1, t.indexj + 1) += massthird;
			M.coeffRef(t.indexj + 2, t.indexj + 2) += massthird;

			M.coeffRef(t.indexk + 0, t.indexk + 0) += massthird; //add triangle area times 1/3
			M.coeffRef(t.indexk + 1, t.indexk + 1) += massthird;
			M.coeffRef(t.indexk + 2, t.indexk + 2) += massthird;

		}


	};

	/*
	Calculates the stretching energy
	*/
	Eigen::VectorXd stretchingEnergy(Triangle& t) 
	{
		
		for (int i = 0; i < F.rows(); i++)
		{
			t = triangles[i];
			//compute triangle contributions to force/stiffness
			t.computeLocalValues(V);

			//distribute triangle contributions to global stiffness matrix
			t.distributeLocalToGlobal(f, K);
		}

	
	};

	/*
	Updates the vertex positions of the cloth by solving linear system with Conjugate Gradient method
	Assemble into Ax = b form followin Eq. 15 in Baraff/Witkin with df/dv = 0
	*/
	void updateCloth(float dt)
	{
		Triangle t;
		Eigen::SparseMatrix<double> A;
		Eigen::MatrixXd f0, v0, b;
		f0 = f;
		v0 = V;
		for (int i = 0; i < F.rows(); i++)
		{
			t = triangles[i];
			//compute triangle contributions to force/stiffness
			t.computeLocalValues(V);

			//distribute triangle contributions to global stiffness matrix
			t.distributeLocalToGlobal(f, K);
		}

		A = (M - dt * dt * K);
		b = (f0 + dt * K*v0);

		//Set up solver
		Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> cg;
		cg.compute(A);
		cg.setMaxIterations(10);

		Eigen::MatrixXd x;
		x = cg.solve(b);
		

	}
	Eigen::VectorXd bendingEnergy() {};
	Eigen::VectorXd shearingEnergy() {};
	
	Eigen::MatrixXd stretchingStiffness() {};

	void removeSelfIntersections() {};
	void solve() {};
	void resolveContactForces() {};
	void collisionDetection() {};

	
};
