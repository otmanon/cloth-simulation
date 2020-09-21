#pragma once
#include <eigen/Dense>
#include <eigen/Sparse>

struct Triangle
{
	//This stuff remains fixed throughout simulation
	int indexi, indexj, indexk;						//indices of the i, j and k vertices in our triangle
	
	Eigen::Vector2d uvi, uvj, uvk;
	float ui, uj, uk,								//u component of i, j, k vertex
		  vi, vj, vk;								//v component of i, j, k vertex

	float a;										//triangle area in uv coordinates

	Eigen::Matrix2d D;								// Inverted Matrix used to calculate dw/du and dw/dv , as given by Eq9 in baraff/witkin

	Eigen::Matrix3d dwudi, dwudj, dwudk,			//derivatives of wu and wv; THESE CAN BE PRECOMPUTED as they dont depend on x
					dwvdi, dwvdj, dwvdk;

	float k;
	// This stuff changes and needs to be recomputed throughout sim

	Eigen::Vector3d xi, xj, xk;						//3D world coordinates of i, j, k vertex

	Eigen::Vector2d C;								//Condition vector who's components are Cu, and Cv
	float Cu;										//Condition on wu component of stretch
	float Cv;										//Condition on wv component of stretch
	float uEnergy, vEnergy, energy;					// stretching Energy = 0.5*Cu*Cu + 0.5 *Cv*Cv; 
	Eigen::Vector3d fui, fuj, fuk,
					fvi, fvj, fvk;

	Eigen::Matrix3d Kij, Kik, Kjk,					//Stiffness matrices
					Kji, Kki, Kkj;

	Eigen::MatrixXd w;								//dw/duv row vector as described by Baraff/Witkin eq 9; 3x2 matrix
	Eigen::VectorXd wu, wv;							//Individual collumns of dwduv.

	Eigen::Vector3d dCudi, dCudj, dCudk,			//Force vector contribution from the u stretch acting on vertices i, j, k
					dCvdi, dCvdj, dCvdk;			//Force vector contribution from the v stretch acting on vertices i, j, k

	Eigen::Matrix3d	d2Cudidi, d2Cudidj, d2Cudidk,	//Stiffness Contribution from u stretch... how the u stretch force on i changes by moving i, j, and k appropriately
					d2Cudjdi, d2Cudjdj, d2Cudjdk,	//Stiffness Contribution from u stretch... how the u strecth force on j changes by moving i, j, and k appropriately
					d2Cudkdi, d2Cudkdj, d2Cudkdk;	//Stiffness Contribution from u stretch... how the u stretch force on k changes by moving i, j, and k appropriately

	Eigen::Matrix3d	d2Cvdidi, d2Cvdidj, d2Cvdidk,	//Stiffness Contribution from v stretch... how the v stretch force on i changes by moving i, j, and k appropriately
					d2Cvdjdi, d2Cvdjdj, d2Cvdjdk,	//Stiffness Contribution from v stretch... how the v strecth force on j changes by moving i, j, and k appropriately
					d2Cvdkdi, d2Cvdkdj, d2Cvdkdk;	//Stiffness Contribution from v stretch... how the v stretch force on k changes by moving i, j, and k appropriately

	/*
	fills out initialized fields for triangle. the indices, u values and vector positions
	Inputs:
		V : |V|x Ndim Matrix of Vertices in mesh
		F : |F| x 3 matrix of faces (3 vertex indices per face)
		UV: |V| x 2 matrix of uv coordinates (one uv coord per vertex)
		index: index of face we are interested in

	*/
	void init(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& UV, int index)
	{
		Eigen::Vector3i face = F.row(index);
		indexi = face(0);
		indexj = face(1);
		indexk = face(2);

		xi = V.row(indexi);
		xj = V.row(indexj);
		xk = V.row(indexk);

		ui = UV.row(indexi)(0); vi = UV.row(indexi)(1);
		uj = UV.row(indexj)(0); vj = UV.row(indexj)(1);
		uk = UV.row(indexk)(0); vk = UV.row(indexk)(1);
		uvi << ui, vi;
		uvj << uj, vj;
		uvk << uk, vk;

		//area = 1/2(duv1 X duv2);
		Eigen::Vector2d duv1 = uvk - uvi;
		Eigen::Vector2d duv2 = uvj - uvi;
		a = 0.5f*(duv1.cross(duv2).norm());
		k = 1;
		//calculate D as given by Eq. 9 in Baraff/Witkin
		Eigen::Matrix2d Di;
		Di   << uj - ui, uk - ui,
				vj - vi, vk - vi;
		D = Di.inverse();

		//Gotta do some vector calculus to derive this from baraff witkin
		//Get derivative of w with respect to xi, xj xk... (How the deformation map from uv coordinates to world coordinates changes by moving uv coordinates and xi/j/k coord
		dwudi = -Eigen::Matrix3d::Identity()*(D(0, 0) + D(1, 0));
		dwudj = Eigen::Matrix3d::Identity()*D(0, 0);
		dwudk = Eigen::Matrix3d::Identity()*D(1, 0);

		dwvdi = -Eigen::Matrix3d::Identity()*(D(0, 1) + D(1, 1));
		dwvdj = Eigen::Matrix3d::Identity()*D(0, 1);
		dwvdk = Eigen::Matrix3d::Identity()*D(1, 1);


		//Calculate first derivative of condition quantities... not dependent on x
		//Use information of dwudi/j/k to calculate C derivatives... gotta derive it yourself using vector calc :) ihatemylife
		float wuNorm = wu.norm();
		dCudi = a * dwudi * wu / wuNorm;
		dCudj = a * dwudj * wu / wuNorm;
		dCudk = a * dwudk * wu / wuNorm;

		float wvNorm = wv.norm();
		dCvdi = a * dwvdi * wv / wvNorm;
		dCvdj = a * dwvdj * wv / wvNorm;
		dCvdk = a * dwvdk * wv / wvNorm;


		//Calculate second derivative of Condition Quantities
		//Do u component
		float wuNorm = wu.norm();
		float wuNorm3 = wuNorm * wuNorm * wuNorm;
		Eigen::Matrix3d matrixTerm = Eigen::Matrix3d::Identity() / wuNorm - wu * wu.transpose() / wuNorm3;
		d2Cudidi = a * (dwudi*dwudi*matrixTerm);
		d2Cudidj = a * (dwudi*dwudj*matrixTerm);
		d2Cudidk = a * (dwudi*dwudk*matrixTerm);

		d2Cudjdi = a * (dwudj*dwudi*matrixTerm);
		d2Cudjdj = a * (dwudj*dwudj*matrixTerm);
		d2Cudjdk = a * (dwudj*dwudk*matrixTerm);

		d2Cudkdi = a * (dwudk*dwudi*matrixTerm);
		d2Cudkdj = a * (dwudk*dwudj*matrixTerm);
		d2Cudkdk = a * (dwudk*dwudk*matrixTerm);

		//Do v component
		float wvNorm = wv.norm();
		float wvNorm3 = wvNorm * wvNorm * wvNorm;
		Eigen::Matrix3d matrixTerm = Eigen::Matrix3d::Identity() / wvNorm - wv * wv.transpose() / wvNorm3;
		d2Cvdidi = a * (dwvdi*dwvdi*matrixTerm);
		d2Cvdidj = a * (dwvdi*dwvdj*matrixTerm);
		d2Cvdidk = a * (dwvdi*dwvdk*matrixTerm);

		d2Cvdjdi = a * (dwvdj*dwvdi*matrixTerm);
		d2Cvdjdj = a * (dwvdj*dwvdj*matrixTerm);
		d2Cvdjdk = a * (dwvdj*dwvdk*matrixTerm);

		d2Cvdkdi = a * (dwvdk*dwvdi*matrixTerm);
		d2Cvdkdj = a * (dwvdk*dwvdj*matrixTerm);
		d2Cvdkdk = a * (dwvdk*dwvdk*matrixTerm);

	}

	/*
	Computes dw/duv, a row vector that describes how stretched the triangle is compared to its reference uv coordinates
	dw/duv is given by equation 9
	*/
	void computeDwduv()
	{
		Eigen::MatrixXd dx;
		dx.resize(3, 2);
		dx.col(0) = xj - xi;
		dx.col(1) = xk - xi;
		w = dx * D;
		wu = w.col(0);
		wv = w.col(1);
	}

	/*
	Computes C, the stretch condition associated with stretching the vertices in the u and v direction. 
	The cloth system's dynamics moves the cloth vertices in a direction as to minimize this energy
	C is given by eq 10 in Baraff/Witkin
	*/
	float computeC()
	{
		computeDwduv();
		C << wu.norm() - 1, wv.norm - 1; //TODO: replace with bu, bv
		C *= a;
		Cu = C(0);
		Cv = C(1);
	}


	/*
	Computes local energy, both the u stretch energy and the v stretch energy, as well as their sum
	*/
	void computeLocalEnergy()
	{
		computeC();
		energy = 0.5f*C.transpose()*C;
		uEnergy = 0.5f*Cu*Cu;
		vEnergy = 0.5f*Cv*Cv;
	}


	/*
	Fills force vectors for each index. Make sure computeLocalEnergy was called before this
	*/
	void computeLocalForces()
	{


		//Compute force info using info above
		fui = -k * dCudi * Cu;
		fuj = -k * dCudj * Cu;
		fuk = -k * dCudk * Cu;
						   
		fvi = -k * dCvdi * Cv;
		fvj = -k * dCvdj * Cv;
		fvk = -k * dCvdk * Cv;
	}

	/*
	Fills Kij, Kik and Kjk... the other Kji, Kki and Kkj are transposes of the respective matrices from before.
	Make sure computeLocalEnergy was called before calling this as we need some variables set in that method
	I do not wish this vector calc upon my worst enemy
	*/
	void computeLocalStiffnesses()
	{
		

	}
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

	void calculateEnergy() {};

	/*
	Calculates the stretching energy
	*/
	Eigen::VectorXd stretchingEnergy() {
		std::vector<Triangle> triangles(F.rows());
		Triangle triangle;
		for (int i = 0; i < F.rows(); i++)
		{
			//fill triangle info
			triangle.init(V, F, UV, i);
			
			//compute triangle contributions to force/stiffness
			triangle.computeLocalForces();
			triangle.computeLocalStiffnesses();

			//distribute triangle contributions to global stiffness matrix
		
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
