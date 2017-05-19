#ifndef MODEL_H_
#define MODEL_H_


extern "C"{
#include<OpenCAL/cal2D.h>
#include<OpenCAL/cal2DBuffer.h>
}
#include<cstdint>
#include<array>

#define FTNREF1D(i_index,i_lb) ((i_index)-(i_lb))
#define FTNREF2D(i_index,j_index,i_size,i_lb,j_lb) ((i_size)*(j_index-(j_lb))+(i_index)-(i_lb))




using CALBuffer = CALreal*;


typedef struct grid_type
{
    double  xmin, ymin, xmax, ymax;
    int x_cells, y_cells;
} grid_type;


typedef struct
{
    int defined;  //logical
    double density,
    energy,
    xvel,
    yvel;
    int geometry;
    double xmin,
    xmax,
    ymin,
    ymax,
    radius;
} state_type;

typedef struct field_type
{
    int left, right, bottom, top ,left_boundary, right_boundary,
    bottom_boundary, top_boundary;
    int x_min, y_min, x_max ,y_max;
} field_type;



/**----------OPENCAL SUBSTATES--------------**/
typedef struct{

    struct CALSubstate2Dr * density0;
    struct CALSubstate2Dr * density1;
    struct CALSubstate2Dr * energy0;
    struct CALSubstate2Dr * energy1;
    struct CALSubstate2Dr * pressure;
    struct CALSubstate2Dr * viscosity;
    struct CALSubstate2Dr * soundspeed;

    struct CALSubstate2Dr * xvel0;
    struct CALSubstate2Dr * xvel1;
    struct CALSubstate2Dr * yvel0;
    struct CALSubstate2Dr * yvel1;
    struct CALSubstate2Dr * vol_flux_x;
    struct CALSubstate2Dr * mass_flux_x;
    struct CALSubstate2Dr * vol_flux_y;
    struct CALSubstate2Dr * mass_flux_y;

    struct CALSubstate2Dr * volume;

    //work arrays
    struct CALSubstate2Dr * work_array1;
    struct CALSubstate2Dr * work_array2;
    struct CALSubstate2Dr * work_array3;
    struct CALSubstate2Dr * work_array4;
    struct CALSubstate2Dr * work_array5;
    struct CALSubstate2Dr * work_array6;
    struct CALSubstate2Dr * work_array7;

    struct CALSubstate2Dr * xarea;
    struct CALSubstate2Dr * yarea;

    CALBuffer cellx;
    CALBuffer celldx;
    CALBuffer celly;
    CALBuffer celldy;


    CALBuffer vertexx;
    CALBuffer vertexy;
    CALBuffer vertexdx;
    CALBuffer vertexdy;

} CLOVER_SBS ;

/** OpenCAL cloverleaf model **/
extern CALModel2D* clover_model;

extern CLOVER_SBS clover_sbs;


extern float   g_version;
extern int     g_ibig ;
extern double  g_small ;
extern double  g_big  ;
extern int     g_name_len_max, g_xdir, g_ydir ;

extern int test_problem;
extern int profiler_on;
extern int state_max;
extern bool complete;
extern int number_of_states;

extern int step;
extern double end_time;
extern int end_step;
extern int visit_frequency;
extern int summary_frequency;
extern bool use_vector_loops;
extern bool advect_x; //logical

extern double dtold, dt, clover_time, dtinit, dtmin, dtmax, dtrise, dtu_safe, dtv_safe, dtc_safe,
dtdiv_safe, dtc, dtu, dtv, dtdiv;


extern grid_type grid;
extern state_type * states;
extern field_type field;
extern int g_rect, g_circ, g_point;

extern int jdt, kdt;

typedef struct {
    //    /**----------Cloverleaf Vars/Consts--------------**/
    //     const float   g_version = 1.0;

    //    int     g_ibig = 640000;

    //    double  g_small = 1.0e-16;

    //    double  g_big  = 1.0e+21;

    //    int     g_name_len_max = 255 ,
    //            g_xdir = 1,
    //            g_ydir = 2;

    //    int     number_of_states;

    //    double  dtold,
    //            dt,
    //            clover_time,
    //            dtinit,
    //            dtmin,
    //            dtmax,
    //            dtrise,
    //            dtu_safe,
    //            dtv_safe,
    //            dtc_safe,
    //            dtdiv_safe,
    //            dtc,
    //            dtu,
    //            dtv,
    //            dtdiv;

    //    int x_min,
    //        y_min,
    //        x_max,
    //        y_max,
    //        x_cells,
    //        y_cells;

    //    grid_type grid;



} CLOVER_MODEL ;



enum  FIELD : std::int8_t  {
    DENSITY0 = 0
    ,  DENSITY1
    ,  ENERGY0
    ,  ENERGY1
    ,  PRESSURE
    ,  VISCOSITY
    ,  SOUNDSPEED
    ,  XVEL0
    ,  XVEL1
    ,  YVEL0
    ,  YVEL1
    ,  VOL_FLUX_X
    ,  VOL_FLUX_Y
    ,  MASS_FLUX_X
    ,  MASS_FLUX_Y
    , NUM_FIELDS
};

extern std::array<int, FIELD::NUM_FIELDS> fields;
extern std::array<CALSubstate2Dr*,FIELD::NUM_FIELDS> FIELDS_SUBSTATE;
extern unsigned int FIELD_DEPTH;


auto get = [&](auto sbs, auto x, auto y){return calGet2Dr(clover_model,sbs,x,y); };
auto set = [&](auto sbs, auto x, auto y, auto val){return calSet2Dr(clover_model,sbs,x,y,val); };
auto setCurr = [&](auto sbs, auto x, auto y, auto val){return calSetCurrent2Dr(clover_model,sbs,x,y,val); };
#endif // MODEL_H_
