#include "PhysicsHook.h"
#include "Object.h"
#include <igl/readOBJ.h>
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/harmonic.h>

class ClothHook : public PhysicsHook
{
private:
	
	double dt;							//time step size

	//std::vector<Object> objects;		//list of all objects in scene that are not cloth
	Object floor;						//object representing floor
	Cloth cloth;						//cloth object
		
	Eigen::MatrixXd renderV;
	Eigen::MatrixXi renderF;


	Eigen::MatrixXd V_uv;

public:
	ClothHook() : PhysicsHook() {}

	virtual bool simulateOneStep()
	{
		
		return false;
	}

	virtual void drawGUI(igl::opengl::glfw::imgui::ImGuiMenu &menu)
	{

	}

	virtual void initSimulation()
	{
		dt = 1e-3;
	
		igl::readOBJ("data/floor.obj", floor.V, floor.F);
		floor.V.col(1).array() -= 1.5;
		
		igl::readOBJ("data/cloth.obj", cloth.V, cloth.F);
		cloth.V *= 0.25f; //scale to reasonable
		cloth.V.col(1).array() += 1.5; //translate cloth upwards
	
		std::cout << cloth.F << std::endl;
			

		// Find the open boundary
		Eigen::VectorXi bnd;
		igl::boundary_loop(cloth.F, bnd);

		// Map the boundary to a circle, preserving edge proportions
		Eigen::MatrixXd bnd_uv;
		igl::map_vertices_to_circle(cloth.V, bnd, bnd_uv);

		// Harmonic parametrization for the internal vertices
		igl::harmonic(cloth.V, cloth.F, bnd, bnd_uv, 1, V_uv);

		// Scale UV to make the texture more clear
		V_uv *= 5;
	}

	virtual void updateRenderGeometry()
	{
		renderV.resize(floor.V.rows() + cloth.V.rows(), 3);
		renderF.resize(floor.F.rows() + cloth.F.rows(), 3);

		renderV << floor.V, cloth.V;
	
		renderF << floor.F, cloth.F.array() + floor.V.rows();
	}


	virtual void renderRenderGeometry(igl::opengl::glfw::Viewer &viewer)
	{
		viewer.data().set_mesh(cloth.V, cloth.F);
		viewer.data().set_uv(V_uv);
		viewer.data().show_texture = true;
	}

};