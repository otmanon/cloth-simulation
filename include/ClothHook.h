#include "PhysicsHook.h"
#include "Object.h"
#include <igl/readOBJ.h>
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/harmonic.h>

class ClothHook : public PhysicsHook
{
private:
	
	float dt;							//time step size
	float k = 1000;							//stiffness constant

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
		cloth.updateCloth(dt, k);
		return false;
	}

	virtual void drawGUI(igl::opengl::glfw::imgui::ImGuiMenu &menu)
	{
		if (ImGui::CollapsingHeader("Simulation Parameters"))
		{
			ImGui::SliderFloat("TimeStep", &dt, 0.001, 1.0f, "%.3f", 10.0f)  ;
			ImGui::SliderFloat("Stiffness", &k, 100.0f, 100000.0f, "%.1f", 10.0f);
			ImGui::SliderFloat("Mass Density", &cloth.massDensity, 0.1f, 5.0f, "%.1f", 10.0f);
			ImGui::SliderFloat("Gravity", &cloth.gravity, 0.05f, 10.0f, "%.1f", 10.0f);
		}
	
	}

	virtual void initSimulation()
	{
		dt = 1e-1;
	
		igl::readOBJ("data/floor.obj", floor.V, floor.F);
		floor.V.col(1).array() -= 1.5;
		
		igl::readOBJ("data/cloth4.obj", cloth.V, cloth.F);
		//cloth.V *= 0.25f; //scale to reasonable
		
	
		cloth.buildUVCoords();
		cloth.rotateAboutX(90);
		cloth.V.col(1).array() += 1.5; //translate cloth upwards
		cloth.initCloth();

		// Scale UV to make the texture more clear
		
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
	//	viewer.data().set_uv(cloth.UV);
	//	viewer.data().show_texture = true;
	}

};