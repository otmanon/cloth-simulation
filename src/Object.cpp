#include "Object.h"

void Cloth::initCloth()
{
	fMat.resize(V.rows(), 3);
	f.resize(3 * V.rows());
	v.resize(3 * V.rows());
	x.resize(3 * V.rows());

	f.setZero();
	v.setZero();
	x.setZero();
	convertMatToVector(V, x);
	fMat.setZero();

	K.resize(3 * V.rows(), 3 * V.rows());
	M.resize(3 * V.rows(), 3 * V.rows());
	M.setIdentity();

	triangles = std::vector<Triangle>(F.rows());
	Triangle t;


	for (int i = 0; i < F.rows(); i++)
	{
		//fill triangle info
		t.init(V, F, UV, i);
		triangles[i] = t;

	}
	calculateMassMatrix();
	solver.setUpSi(V);

}

void Cloth::buildUVCoords()
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

	UV *= 2;

}



void Cloth::calculateMassMatrix()
{
	Triangle t;

	float onethird = 1.0f / 3.0f;
	float massthird;
	M.setZero();
	for (int i = 0; i < F.rows(); i++)
	{
		t = triangles[i];

		massthird = massDensity * onethird * t.a;

		//Fill mass matrix						 1 ;//
		M.coeffRef(3 * t.indexi + 0, 3 * t.indexi + 0) += massthird; //add triangle area times 1/3
		M.coeffRef(3 * t.indexi + 1, 3 * t.indexi + 1) += massthird;
		M.coeffRef(3 * t.indexi + 2, 3 * t.indexi + 2) += massthird;
							 							 
		M.coeffRef(3 * t.indexj + 0, 3 * t.indexj + 0) += massthird; //add triangle area times 1/3
		M.coeffRef(3 * t.indexj + 1, 3 * t.indexj + 1) += massthird;
		M.coeffRef(3 * t.indexj + 2, 3 * t.indexj + 2) += massthird;
					 									 
		M.coeffRef(3 * t.indexk + 0, 3 * t.indexk + 0) += massthird; //add triangle area times 1/3
		M.coeffRef(3 * t.indexk + 1, 3 * t.indexk + 1) += massthird;
		M.coeffRef(3 * t.indexk + 2, 3 * t.indexk + 2) += massthird;

	}
	
}



void Cloth::updateCloth(float dt, float k)
{
	Triangle t;
	Eigen::SparseMatrix<double> A;
	Eigen::VectorXd v0, b;

	v0 = v;
	calculateMassMatrix();
	for (int i = 0; i < F.rows(); i++)
	{
		t = triangles[i];
		//compute triangle contributions to force/stiffness
		t.computeLocalValues(V, k);

		//distribute triangle contributions to global stiffness matrix
		t.distributeLocalToGlobal(fMat, K);
	}

	applyGravity(fMat);
	convertMatToVector(fMat, f);


	Eigen::VectorXd y(V.rows() * 3);
	y.setZero();


	A = (M - dt * dt * K);
	b = (dt * f) + dt * dt * K*v0 + dt * K*y;


	Eigen::VectorXd dv = solver.solve(A, b);

	v = v0 + dv;	//update velocity
	x += dt * v;
	convertVecToMat(x, V.cols(), V);
	K.setZero();
	f.setZero();
	fMat.setZero();
}


void Cloth::applyGravity(Eigen::MatrixXd& fMat)
{
	fMat.rowwise() += Eigen::RowVector3d(0.0, -gravity, 0.0);
}

void Cloth::rotateAboutX(float degrees)
{
	float pi = 3.14159;
	float deg2rad = pi / 180;
	float theta = degrees * deg2rad;
	Eigen::Matrix3d rot;
	rot << 1, 0, 0,
		0, cos(theta), -sin(theta),
		0, sin(theta), sin(theta);
	V = (rot * V.transpose()).transpose();
}


void Cloth::convertMatToVector(Eigen::MatrixXd& mat, Eigen::VectorXd& vec)
{
	vec.resize(mat.rows()*mat.cols());
	int index = 0;
	for (int i = 0; i < mat.rows(); i++)
	{
		for (int j = 0; j < mat.cols(); j++)
		{
			vec(index) = mat(i, j);
			index++;
		}
	}
}


void Cloth::convertVecToMat(Eigen::VectorXd& vec, int dim, Eigen::MatrixXd& mat)
{
	int index;
	assert(vec.rows() % dim == 0); //make sure vector and dims are approriate and compatible
	mat.resize(vec.rows() / dim, dim);
	int i, j;
	for (int index = 0; index < vec.rows(); index++)
	{
		i = (int)(index / dim);
		j = index % dim;
		mat(i, j) = vec(index);
	}
}