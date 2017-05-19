#ifndef _TIMESTEP_
#define _TIMESTEP_

#include<ideal_gas.h>
#include<viscosity.h>
#include<calc_dt.h>

void timestep( FILE* g_out){
    int jldt, kldt;
    double dtlp;
    double x_pos, y_pos, xl_pos, yl_pos;
    //char dt_control[8];
    char dtl_control[8];

    int small = 0;

    x_pos  = 0.0;
    y_pos  = 0.0;
    xl_pos = 0.0;
    yl_pos = 0.0;

    ideal_gas(false);

    //---updateHalo----
    setFields(0);
    fields[FIELD::DENSITY0] = 1;
    fields[FIELD::ENERGY0] = 1;
    fields[FIELD::PRESSURE] = 1;
    fields[FIELD::XVEL0] = 1;
    fields[FIELD::YVEL0] = 1;
    FIELD_DEPTH = 1;
    update_halo();
    //---updateHalo----

    viscosity();

    calc_dt(xl_pos,yl_pos,jldt,kldt);
dtlp = dt_min_val;
    if(dtlp <= dt) {
        dt=dtlp;
        //dt_control=dtl_control;
        x_pos=xl_pos;
        y_pos=yl_pos;
        jdt=jldt;
        kdt=kldt;
    }

    dt = std::min( std::min(dt, (dtold * dtrise)), dtmax);
    if(dt < dtmin) small=1;
    printf(   " Step %d time %11.7lf control %s timestep  %3.2E  %d, %d x  %E  y %E\n",
                step,   clover_time,    dtl_control,dt,          jdt, kdt,  x_pos,y_pos);
    fprintf(g_out,
            " Step %d time %11.7lf control %s timestep  %3.2E  %d, %d x  %E  y %E\n",
            step,   clover_time,    dtl_control,dt,          jdt, kdt,  x_pos,y_pos);

    if(small == 1) {
        printf("timestep :small timestep\n");
        exit(-2);
    }

    dtold = dt;

}

#endif //_TIMESTEP_
