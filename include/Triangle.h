#pragma once
#include <eigen/Dense>
#include <eigen/Sparse>
struct Triangle
{
public:
	//This stuff remains fixed throughout simulation
	int indexi, indexj, indexk;						//indices of the i, j and k vertices in our triangle
	float a;										//triangle area in uv coordinates
private:
	Eigen::Vector2d uvi, uvj, uvk;
	float ui, uj, uk,								//u component of i, j, k vertex
		vi, vj, vk;								//v component of i, j, k vertex



	Eigen::Matrix2d D;								// Inverted Matrix used to calculate dw/du and dw/dv , as given by Eq9 in baraff/witkin

	Eigen::Matrix3d dwudi, dwudj, dwudk,			//derivatives of wu and wv; THESE CAN BE PRECOMPUTED as they dont depend on x
		dwvdi, dwvdj, dwvdk;

	float k;										//spring stiffness? idk the term here but its a constant that makes cloth stiffer

	float du1;										//Stores u spacing between j and i
	float dv1;										//Stores v spacing between j and i
	float du2;										//Stores u spacing between k and i
	float dv2;										//Stores v spacing between k and i
	
	// This stuff changes and needs to be recomputed throughout sim
	Eigen::Vector3d xi, xj, xk;						//3D world coordinates of i, j, k vertex

	Eigen::Vector2d C;								//Condition vector who's components are Cu, and Cv
	float Cu;										//Condition on wu component of stretch
	float Cv;										//Condition on wv component of stretch
	float Cs;										//Condition on shear 
	float uEnergy, vEnergy, sEnergy, energy;					// stretching Energy = 0.5*Cu*Cu + 0.5 *Cv*Cv; 
	
	Eigen::Vector3d fui, fuj, fuk,
					fvi, fvj, fvk,
					fsi, fsj, fsk;

	Eigen::Matrix3d Kuij, Kuik, Kujk,					//Stiffness matrices u comp
					Kuji, Kuki, Kukj,
					Kuii, Kujj, Kukk;

	Eigen::Matrix3d Kvij, Kvik, Kvjk,					//Stiffness matrices v comp
					Kvji, Kvki, Kvkj,
					Kvii, Kvjj, Kvkk;

	Eigen::Matrix3d Ksij, Ksik, Ksjk,					//Stiffness matrices shear
					Ksji, Kski, Kskj,
					Ksii, Ksjj, Kskk;

	Eigen::Matrix3d Kij, Kik, Kjk,					//Stiffness matrices	total
					Kji, Kki, Kkj,
					Kii, Kjj, Kkk;

	Eigen::MatrixXd w;								//dw/duv row vector as described by Baraff/Witkin eq 9; 3x2 matrix
	Eigen::VectorXd wu, wv;							//Individual collumns of dwduv.

	Eigen::Vector3d dCudi, dCudj, dCudk,			//Force vector contribution from the u stretch acting on vertices i, j, k
					dCvdi, dCvdj, dCvdk,			//Force vector contribution from the v stretch acting on vertices i, j, k
					dCsdi, dCsdj, dCsdk;

	Eigen::Matrix3d	d2Cudidi, d2Cudidj, d2Cudidk,	//Stiffness Contribution from u stretch... how the u stretch force on i changes by moving i, j, and k appropriately
					d2Cudjdi, d2Cudjdj, d2Cudjdk,				//Stiffness Contribution from u stretch... how the u strecth force on j changes by moving i, j, and k appropriately
					d2Cudkdi, d2Cudkdj, d2Cudkdk;				//Stiffness Contribution from u stretch... how the u stretch force on k changes by moving i, j, and k appropriately

	Eigen::Matrix3d	d2Cvdidi, d2Cvdidj, d2Cvdidk,	//Stiffness Contribution from v stretch... how the v stretch force on i changes by moving i, j, and k appropriately
					d2Cvdjdi, d2Cvdjdj, d2Cvdjdk,		//Stiffness Contribution from v stretch... how the v strecth force on j changes by moving i, j, and k appropriately
					d2Cvdkdi, d2Cvdkdj, d2Cvdkdk;		//Stiffness Contribution from v stretch... how the v stretch force on k changes by moving i, j, and k appropriately

	Eigen::Matrix3d	d2Csdidi, d2Csdidj, d2Csdidk,	//Stiffness Contribution from v stretch... how the v stretch force on i changes by moving i, j, and k appropriately
					d2Csdjdi, d2Csdjdj, d2Csdjdk,		//Stiffness Contribution from v stretch... how the v strecth force on j changes by moving i, j, and k appropriately
					d2Csdkdi, d2Csdkdj, d2Csdkdk;		//Stiffness Contribution from v stretch... how the v stretch force on k changes by moving i, j, and k appropriately



	float dwudiScalar;								// Stores scalar factor we will need to calculate various C derivatives
	float dwudjScalar;								// Stores scalar factor we will need to calculate various C derivatives
	float dwudkScalar;								// Stores scalar factor we will need to calculate various C derivatives
			  		 								
	float dwvdiScalar;								// Stores scalar factor we will need to calculate various C derivatives
	float dwvdjScalar;								// Stores scalar factor we will need to calculate various C derivatives
	float dwvdkScalar;								// Stores scalar factor we will need to calculate various C derivatives

public:

	/*
	fills out initialized fields for triangle. the indices, u values and vector positions
	Inputs:
		V : |V|x Ndim Matrix of Vertices in mesh
		F : |F| x 3 matrix of faces (3 vertex indices per face)
		UV: |V| x 2 matrix of uv coordinates (one uv coord per vertex)
		index: index of face we are interested in

	*/
	void init(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& UV, int index);

	/*
	Computes dw/duv, a row vector that describes how stretched the triangle is compared to its reference uv coordinates
	dw/duv is given by equation 9
	*/
	void computeDwduv();

	/*
	Computes C, the stretch condition associated with stretching the vertices in the u and v direction.
	The cloth system's dynamics moves the cloth vertices in a direction as to minimize this energy
	C is given by eq 10 in Baraff/Witkin
	*/
	void computeC();

	/*
	Computes local energy, both the u stretch energy and the v stretch energy, as well as their sum
	*/
	void computeLocalEnergy();


	/*
	Fills force vectors for each index. Make sure computeLocalEnergy was called before this
	*/
	void computeLocalForces();

	/*
	Fills Kij, Kik and Kjk... the other Kji, Kki and Kkj are transposes of the respective matrices from before.
	Make sure computeLocalEnergy was called before calling this as we need some variables set in that method

	*/
	void computeLocalStiffnesses();

	/*
	Called every timestep for the cloth... updartes x position values for the three vertices and recalculates all the stiffness/Conditions/force vectors locally
	*/
	void computeLocalValues(Eigen::MatrixXd V, float kConst);

	/*
	Now we have our local force vectors and stiffness matrices.... We need to figure out where these matrices and vectors lie
	on their global counterpart using the triangle's indices.
	*/
	void distributeLocalToGlobal(Eigen::MatrixXd& Force, Eigen::SparseMatrix<double>& K);

};