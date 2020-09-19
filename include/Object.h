#pragma once
#include <eigen/Dense>
#include <eigen/Sparse>

struct Triangle
{
	int indexi, indexj, indexk;			//indices of the i, j and k vertices in our triangle
	
	Eigen::Vector3d xi, xj, xk;			//3D world coordinates of i, j, k vertex

	float ui, uj, uk;					//u component of i, j, k vertex
	float vi, vj, vk;					//v component of i, j, k vertex

	float Cu;							//Condition on wu component of stretch
	float Cv;							//Condition on wv component of stretch

	Eigen::Vector3d dCudi, dCudj, dCudk;  //Force vector contribution from the u stretch acting on vertices i, j, k
	Eigen::Vector3d dCvdi, dCvdj, dCvdk;  //Force vector contribution from the v stretch acting on vertices i, j, k

	Eigen::Matrix3d	d2Cudidi, d2Cudidj, d2Cudidk,	//Stiffness Contribution from u stretch... how the u stretch force on i changes by moving i, j, and k appropriately
					d2Cudjdi, d2Cudjdj, d2Cudjdk,	//Stiffness Contribution from u stretch... how the u strecth force on j changes by moving i, j, and k appropriately
					d2Cudkdi, d2Cudkdj, d2Cudkdk;	//Stiffness Contribution from u stretch... how the u stretch force on k changes by moving i, j, and k appropriately

	Eigen::Matrix3d	d2Cvdidi, d2Cvdidj, d2Cvdidk,	//Stiffness Contribution from v stretch... how the v stretch force on i changes by moving i, j, and k appropriately
					d2Cvdjdi, d2Cvdjdj, d2Cvdjdk,	//Stiffness Contribution from v stretch... how the v strecth force on j changes by moving i, j, and k appropriately
					d2Cvdkdi, d2Cvdkdj, d2Cvdkdk;	//Stiffness Contribution from v stretch... how the v stretch force on k changes by moving i, j, and k appropriately


};
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
	Eigen::MatrixXd K;					// stiffness matrix (df/dx)
	Eigen::MatrixXd y;					//displacement vector to avoid mesh tangling
	Eigen::SparseMatrix<double> M;		//Mass matrix
	Eigen::SparseMatrix<double> W;		//inverse mass matrix used for constraint handling
	Eigen::MatrixXd UV;					//Contains uv coordinates for each vertex in cloth
	
	/*
	Assume cloth is straight and is a square with evenly spaced vertices in both directions
	*/
	void buildUVCoords()
	{
		
		Eigen::Vector3d tlV, trV, blV, brV;		// corner edges of square cloth
		trV = V.rowwise().maxCoeff();
		blV = V.rowwise().minCoeff();
		tlV << blV.x(), trV.y(), 0;
		brV << trV.x(), blV.y(), 0;
		
		float w, h;
		w = trV.x() - blV.x();
		h = trV.y() - blV.y();
		
		Eigen::Vector2d tr_uv, bl_uv;			//topright most corner has uv coordinates (1, 1), bottomleft has uv coords (0, 0)
		tr_uv << 1, 1;
		bl_uv << 0, 0;

		//go through vertices and interpolate UV coords
		UV.resize(V.rows(), 2);
		Eigen::Vector2d uv;
		float x, y, s_x, s_y;					//x, y, and fraction coordinates
		for (int i = 0; i < V.rows(); i++)
		{
			x = V.row(i).x();
			y = V.row(i).y();

			s_x = (x - bl_uv.x()) / w;
			s_y = (y - bl_uv.y()) / h;
			
			uv(0) = 0 * (1 - s_x) + 1 * s_x;
			uv(1) = 0 * (1 - s_y) + 1 * s_y;
			UV.row(i) = uv;
		}



	}

	void calculateEnergy() {};

	/*
	Calculates the stretching energy
	*/
	Eigen::VectorXd stretchingEnergy() {

		for (int i = 0; i < F.rows(); i++)
		{
			//Assemble f, and K matrix block by block, C by C.
			
		
		}

		
	
	};

	Eigen::VectorXd bendingEnergy() {};
	Eigen::VectorXd shearingEnergy() {};
	
	Eigen::MatrixXd stretchingStiffness() {};

	void removeSelfIntersections() {};
	void solve() {};
	void resolveContactForces() {};
	void collisionDetection() {};

	
};
