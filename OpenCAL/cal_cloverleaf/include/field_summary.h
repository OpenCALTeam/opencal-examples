#ifndef _FIELD_SUMMARY_
#define _FIELD_SUMMARY_

#include<model.h>
#include<ideal_gas.h>
double vol= 0.0 , mass = 0.0, ie = 0.0, ke = 0.0, press = 0.0;

void field_summary_kernel(struct CALModel2D  * clover_model,  int i,  int j){
    int cx = FTNREF1D(i,field.x_min-2);
    int cy = FTNREF1D(j, field.y_min-2);
    if(!checkBoundaries(cx,cy,field.x_min,field.x_max,field.x_min,field.x_max))
        return;

    double vsqrd, cell_vol, cell_mass;

    vsqrd = 0.0;


    //vsqrd = vsqrd + 0.25 * ( xvel0[OPS_ACC4(0,0)] * xvel0[OPS_ACC4(0,0)] + yvel0[OPS_ACC5(0,0)] * yvel0[OPS_ACC5(0,0)]);
    vsqrd = vsqrd + 0.25 * ( get(clover_sbs.xvel0,cx,cy)*get(clover_sbs.xvel0,cx,cy) +
                             get(clover_sbs.yvel0,cx,cy)* get(clover_sbs.yvel0,cx,cy)
                             );

    //vsqrd = vsqrd + 0.25 * ( xvel0[OPS_ACC4(1,0)] * xvel0[OPS_ACC4(1,0)] + yvel0[OPS_ACC5(1,0)] * yvel0[OPS_ACC5(1,0)]);
    vsqrd = vsqrd + 0.25 * ( get(clover_sbs.xvel0,cx+1,cy)*get(clover_sbs.xvel0,cx+1,cy) +
                             get(clover_sbs.yvel0,cx+1,cy)* get(clover_sbs.yvel0,cx+1,cy)
                             );

    //vsqrd = vsqrd + 0.25 * ( xvel0[OPS_ACC4(0,1)] * xvel0[OPS_ACC4(0,1)] + yvel0[OPS_ACC5(0,1)] * yvel0[OPS_ACC5(0,1)]);
    vsqrd = vsqrd + 0.25 * ( get(clover_sbs.xvel0,cx,cy+1)*get(clover_sbs.xvel0,cx,cy+1) +
                             get(clover_sbs.yvel0,cx,cy+1)* get(clover_sbs.yvel0,cx,cy+1)
                             );


    //vsqrd = vsqrd + 0.25 * ( xvel0[OPS_ACC4(1,1)] * xvel0[OPS_ACC4(1,1)] + yvel0[OPS_ACC5(1,1)] * yvel0[OPS_ACC5(1,1)]);
    vsqrd = vsqrd + 0.25 * ( get(clover_sbs.xvel0,cx+1,cy+1)*get(clover_sbs.xvel0,cx+1,cy+1) +
                             get(clover_sbs.yvel0,cx+1,cy+1)* get(clover_sbs.yvel0,cx+1,cy+1)
                             );


    cell_vol    = get(clover_sbs.volume,cx,cy);
    cell_mass   = cell_vol * get(clover_sbs.density0,cx,cy);

    vol = vol + cell_vol;
    mass = mass + cell_mass;
    ie = ie + cell_mass * get(clover_sbs.energy0,cx,cy);
    ke = ke + cell_mass * 0.5 * vsqrd;
    press = press + cell_vol * get(clover_sbs.pressure,cx,cy);
}

void field_summary(){
    vol= 0.0 , mass = 0.0, ie = 0.0, ke = 0.0, press = 0.0;
    ideal_gas(false);
    calApplyElementaryProcess2D(clover_model, field_summary_kernel);


    fprintf(stdout,"\n");
    fprintf(stdout,"\n Time %lf\n",clover_time);
    fprintf(stdout,"              %-10s  %-10s  %-10s  %-10s  %-15s  %-15s  %-s\n",
    " Volume"," Mass"," Density"," Pressure"," Internal Energy","Kinetic Energy","Total Energy");
    fprintf(stdout," step:   %3d   %-10.3E  %-10.3E  %-10.3E  %-10.3E  %-15.3E  %-15.3E  %-.3E",
            step, vol, mass, mass/vol, press/vol, ie, ke, ie+ke);

//    fprintf(g_out,"\n");
//    fprintf(g_out,"\n Time %lf\n",clover_time);
//    fprintf(g_out,"              %-10s  %-10s  %-10s  %-10s  %-15s  %-15s  %-s\n",
//    " Volume"," Mass"," Density"," Pressure"," Internal Energy","Kinetic Energy","Total Energy");
//    fprintf(g_out," step:   %3d   %-10.3E  %-10.3E  %-10.3E  %-10.3E  %-15.3E  %-15.3E  %-.3E",
//            step, vol, mass, mass/vol, press/vol, ie, ke, ie+ke);


}


#endif //_FIELD_SUMMARY_
