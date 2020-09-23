#pragma once
#include "Triangle.h"
#include "CGSolver.h"
struct Object
{
	Eigen::MatrixXd V;					//Vertices
	Eigen::MatrixXi F;					//Faces
	Eigen::MatrixXd E;					//Edges

};
struct Cloth : public Object
{


	Eigen::VectorXd x;

	Eigen::VectorXd energy;				//sum of stretch, shear, and bending energies
	Eigen::MatrixXd fMat;					// |V| x 3 force matrix. Each element is the force vector of the corresponding vertex
	Eigen::VectorXd f;					//Same as above but flattened to vector
	Eigen::VectorXd v;					// |V| x 3 velocity matrix. Each element is the velocity vector of the corresponding vertex
	Eigen::SparseMatrix<double> K;		// stiffness matrix (df/dx)
	Eigen::MatrixXd y;					//displacement vector to avoid mesh tangling
	Eigen::SparseMatrix<double> M;		//Mass matrix
	Eigen::SparseMatrix<double> W;		//inverse mass matrix used for constraint handling
	Eigen::MatrixXd UV;					//Contains uv coordinates for each vertex in cloth
	
	std::vector<Triangle> triangles;	//individual triangles found in cloth. Used to compute locally force and stiffness values

	float massDensity = 1.0f;
	float gravity = 0.1f;
	CGSolver solver;

	/*
	Initializes values of cloth
	*/
	void initCloth();
	

	/*
	Assume cloth is straight and is a square with evenly spaced vertices in both directions
	*/
	void buildUVCoords();


	/*
	fills mass matrix
	*/
	void calculateMassMatrix();

	/*
	Updates the vertex positions of the cloth by solving linear system with Conjugate Gradient method
	Assemble into Ax = b form followin Eq. 15 in Baraff/Witkin with df/dv = 0
	*/
	void updateCloth(float dt, float k);


	/*
	Applies gravity force to fMat, the |V| x 3 Force matrix
	*/
	void applyGravity(Eigen::MatrixXd& fMat);

	/*
	Rotates vertices about x by a certain degrees
	*/
	void rotateAboutX(float degrees);


	/*
	Helper function converts Matrix to flat vector
	Inpout matrix, output vector
	*/
	void convertMatToVector(Eigen::MatrixXd& mat, Eigen::VectorXd& vec);

	/*
	Helpfer function converts flat vector to matrix
	Inputs vector, and vector size, outputs matrix
	*/
	void convertVecToMat(Eigen::VectorXd& vec, int dim, Eigen::MatrixXd& mat);

};
