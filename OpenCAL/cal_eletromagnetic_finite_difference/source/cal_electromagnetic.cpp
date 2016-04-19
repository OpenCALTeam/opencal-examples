//============================================================================
// Name        : cal_elecromagnetic.cpp
// Author      : Davide Spataro
// Version     :
// Copyright   : 
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>

#include <OpenCAL-GL/calgl2D.h>
#include <OpenCAL-GL/calgl2DWindow.h>

using namespace std;

#define SIZE (100)
#define ROWS (SIZE)
#define COLS (SIZE)
#define STEPS (20000)
#define EPSILON (0.01)

//model&materials parameters
#define DELTA_X (0.001)
#define DELTA_Y (0.001)

// declare CA, substate and simulation objects
struct CALModel2D* elmagnModel;							//the cellular automaton
struct CALSubstate2Dr *Q_electric_field;							//the substate Q
struct CALSubstate2Dr *Q_magnetic_field;							//the substate Q
struct CALRun2D* elMagn_simulation;

// The cell's transition function (first and only elementary process)
void elmagnModel_TransitionFunction(struct CALModel2D* heatModel, int i, int j)
{

}

// Callback unction called just before program termination
void exitFunction(void)
{
	//finalizations
	calRunFinalize2D(heat_simulation);
	calFinalize2D(heatModel);
}

int main(int argc, char** argv) {
	// Declare a viewer object
	struct CALGLDrawModel3D* drawModel;
	atexit(exitFunction);
	//cadef and rundef
	elmagnModel = calCADef2D(ROWS, COLS, CAL_MOORE_NEIGHBORHOOD_2D, CAL_SPACE_FLAT, CAL_OPT_ACTIVE_CELLS);
	elMagn_simulation = calRunDef2D(elmagnModel, 1, STEPS, CAL_UPDATE_IMPLICIT);
	//add substates
	Q_electric_field = calAddSubstate2Dr(elmagnModel);
	Q_magnetic_field = calAddSubstate2Dr(elmagnModel);
	//add transition function's elementary processes
	calAddElementaryProcess2D(elmagnModel, elmagnModel_TransitionFunction);

	//simulation run setup
	calRunAddInitFunc2D(elMagn_simulation, elMagn_simulation_SimulationInit);
	calRunInitSimulation2D(elMagn_simulation);	//It is required in the case the simulation main loop is explicitated; similarly for calRunFinalizeSimulation3D
	calRunAddStopConditionFunc2D(elMagn_simulation, elMagn_simulation_SimulationStopCondition);



	// Initialize the viewer
	calglInitViewer("mod2 3D CA viewer", 1.0f, 400, 400, 40, 40, CAL_TRUE, 0);
	//drawModel definition
	drawModel = calglDefDrawModel3D(CALGL_DRAW_MODE_FLAT, "3D view", elMagn_simulation, elMagn_simulation);
	calglAdd3Dr(drawModel, NULL, &Q_electric_field, CALGL_TYPE_INFO_VERTEX_DATA, CALGL_TYPE_INFO_USE_RED_SCALE, CALGL_DATA_TYPE_DYNAMIC);
	calglColor3D(drawModel, 0.9f, 0.5f, 0.5f, 0.2f);
	calglAdd3Dr(drawModel, Q_electric_field, &Q_electric_field, CALGL_TYPE_INFO_COLOR_DATA, CALGL_TYPE_INFO_USE_RED_SCALE , CALGL_DATA_TYPE_DYNAMIC);
	calglAdd3Dr(drawModel, Q_electric_field, &Q_electric_field, CALGL_TYPE_INFO_NORMAL_DATA, CALGL_TYPE_INFO_USE_NO_COLOR, CALGL_DATA_TYPE_DYNAMIC);
//calglDisableLights();
	calglInfoBar2Dr(drawModel, Q_electric_field, "Temperature", CALGL_TYPE_INFO_USE_RED_SCALE, 20, 120, 300, 80);
//calglDisableLights();
	puts("asdasdasd"); /* prints  */
	calglMainLoop3D(argc, argv);
	return EXIT_SUCCESS;
}

