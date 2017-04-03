#include<stdio.h>
#include<stdlib.h>

#include<string>

#include <model.h>
#include <initialise.h>
#include<start.h>
#include<timestep.h>
#include<PdV.h>
using std::string;
FILE* g_in, *g_out;
string in_file, out_file;

/******************************************************************************
* Initialize Global constants and variables
******************************************************************************/

/**----------Cloverleaf Vars/Consts--------------**/

CALModel2D* clover_model;
CLOVER_SBS clover_sbs;


float   g_version = 1.0;
int     g_ibig = 640000;
double  g_small = 1.0e-16;
double  g_big  = 1.0e+21;
int     g_name_len_max = 255 ,
g_xdir = 1,
g_ydir = 2;

int number_of_states;
int test_problem;
int profiler_on;
int state_max;
bool complete; //logical = boolean

std::array<int, FIELD::NUM_FIELDS> fields = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

double end_time;
int step;
int end_step;
int visit_frequency;
int summary_frequency;
bool use_vector_loops;
bool advect_x;
double dtold, dt, clover_time, dtinit, dtmin, dtmax, dtrise, dtu_safe, dtv_safe, dtc_safe,
dtdiv_safe, dtc, dtu, dtv, dtdiv;

grid_type grid;
state_type * states; //global variable holding state info
field_type field; //global variable holding info of fields

int     g_rect=1,
g_circ=2,
g_point=3;


int jdt, kdt;


/******************************************************************************
* Main program
******************************************************************************/

int main(int argc, char **argv)
{
    CLOVER_MODEL* clover_model;
    printf(" Clover version %f\n", g_version);
    /**-------------------------- OPENCAL Initialisation --------------------------**/



    double ct0, ct1, et0, et1;
    /**---------------------initialize and generate chunk----------------------**/
    in_file = "../data/clover.in";
    out_file = "../data/clover.out";
    initialise(clover_model,g_in,g_out,in_file,out_file);

    start(g_out);

    /***************************************************************************
    **-----------------------------hydro loop---------------------------------**
    ***************************************************************************/
    while(true) {
        step = step + 1;
        timestep(g_out);

        PdV(true);
    }



    printf("\nTotal Wall time %lf\n",et1-et0);
    return 0;
}
