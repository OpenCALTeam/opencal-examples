#ifndef CALC_DT_H_
#define CALC_DT_H_

#include<model.h>
#include<utils.h>
#include<algorithm>

double dt_min_val=0.0;

void calc_dt_kernel(struct CALModel2D  * clover_model,  int i,  int j){
    int cx = FTNREF1D(i,field.x_min-2);
    int cy = FTNREF1D(j, field.y_min-2);
    if(!checkBoundaries(cx,cy,field.x_min,field.x_max,field.y_min,field.y_max))
        return;

    double div, dsx, dsy, dtut, dtvt, dtct, dtdivt, cc, dv1, dv2;//, jk_control;
    //dsx=celldx(j)
    dsx= calGetMatrixElement(clover_sbs.celldx , 1, cx , 0 );
    //dsy=celldy(k)
    dsy= calGetMatrixElement(clover_sbs.celldx , 1, cy , 0 );

    //cc=soundspeed(j,k)*soundspeed(j,k)
    const auto sndspd = get(clover_sbs.soundspeed,cx,cy);
    cc = sndspd*sndspd;

    const auto visc = get(clover_sbs.viscosity,cx,cy);
    const auto dens0 = get(clover_sbs.density0,cx,cy);
    cc=cc+2.0*visc/dens0;
    cc=std::max(sqrt(cc),g_small);


    dtct=dtc_safe*std::min(dsx,dsy)/cc;

    div=0.0;

    //dv1=(xvel0(j  ,k)+xvel0(j  ,k+1))*xarea(j  ,k);
    dv1= (get(clover_sbs.xvel0,cx,cy)+get(clover_sbs.xvel0,cx,cy+1))*get(clover_sbs.xarea,cx,cy);
    //dv2=(xvel0(j+1,k)+xvel0(j+1,k+1))*xarea(j+1,k)
    dv2= (get(clover_sbs.xvel0,cx+1,cy) +get(clover_sbs.xvel0,cx+1,cy+1))*get(clover_sbs.xarea,cx+1,cy);
    div=div+dv2-dv1;

    g_small=1.0000000168623835e-16;
    const auto vol = get(clover_sbs.volume,cx,cy);
    // dtut=dtu_safe*2.0*volume(j,k)/MAX(ABS(dv1),ABS(dv2),g_small*volume(j,k))
    double fdv1 = fabs(dv1);
    dtut = dtu_safe*2.0*vol/std::max(
                                    std::max(fabs(dv1),fabs(dv2) ) ,
                                    g_small*vol);

    //dv1=(yvel0(j,k  )+yvel0(j+1,k  ))*yarea(j,k  )
    dv1 = (get(clover_sbs.yvel0,cx,cy) + get(clover_sbs.yvel0,cx+1,cy))*get(clover_sbs.yarea,cx,cy);

    //dv2=(yvel0(j,k+1)+yvel0(j+1,k+1))*yarea(j,k+1)
    dv2 = (get(clover_sbs.yvel0,cx,cy+1)+get(clover_sbs.yvel0,cx+1,cy+1))*get(clover_sbs.yarea,cx,cy+1);

    div=div+dv2-dv1;

    dtvt=dtv_safe*2.0*vol/std::max(std::max(fabs(dv1),fabs(dv2)), g_small * vol);

    div = div/(2.0 * vol);
    if(div < -g_small)
        dtdivt = dtdiv_safe * (-1.0/div);
    else
        dtdivt = g_big;

    std::array<double,5>  mincandidates = {dt_min_val,dtct,dtut,dtvt,dtdivt};
    dt_min_val = *std::min_element(begin(mincandidates),end(mincandidates));
    //dt_min_val=std::min(dt_min_val,std::min(dtct,dtut),                 std::min(dtvt,dtdivt)

}

void calc_dt(double& xl_pos, double& yl_pos,int& jldt, int& kldt) {

    int small;
    double jk_control = 1.1;
    int dtl_control;

    small = 0;
    dt_min_val = g_big;
    jk_control = 1.1;
    calApplyElementaryProcess2D(clover_model, calc_dt_kernel);

    //initialize sizes using global values
    int x_min = field.x_min;
    int x_max = field.x_max;
    int y_min = field.y_min;
    int y_max = field.y_max;


    //Extract the mimimum timestep information
    dtl_control = 10.01 * (jk_control - (int)(jk_control));
    jk_control = jk_control - (jk_control - (int)(jk_control));
    //*jldt = ((int)jk_control)%x_max;
    //*kldt = 1 + (jk_control/x_max);
    jldt = ((int)jk_control)%(x_max);
    kldt = 1 + (jk_control/(x_max));

    //serve a rimappare indici calcolati in fortrain in C.1.0000000168623835e-16
    jldt+=1;
    kldt+=1;

    xl_pos= calGetMatrixElement(clover_sbs.cellx , 1, jldt , 0 );
    yl_pos= calGetMatrixElement(clover_sbs.celly , 1, kldt , 0 );

    if(dt_min_val < dtmin) small=1;

    if(small != 0){
        printf("Timestep information:");
        printf("\nj %i, k %i                 : ",jldt,kldt);
        printf("\nx %i, y %i                 : ",calGetMatrixElement(clover_sbs.cellx , 1, jldt , 0 ),calGetMatrixElement(clover_sbs.celly , 1, kldt , 0 ));
        printf("\ntimestep %f: ",dt_min_val);
        printf("Cells velocities:");
        //          WRITE(0,*) xvel0(jldt  ,kldt  ),yvel0(jldt  ,kldt  )
        //          WRITE(0,*) xvel0(jldt+1,kldt  ),yvel0(jldt+1,kldt  )
        //          WRITE(0,*) xvel0(jldt+1,kldt+1),yvel0(jldt+1,kldt+1)
        //          WRITE(0,*) xvel0(jldt  ,kldt+1),yvel0(jldt  ,kldt+1)
        //          WRITE(0,*) 'density, energy, pressure, soundspeed '
        //          WRITE(0,*) density0(jldt,kldt),energy0(jldt,kldt),pressure(jldt,kldt),soundspeed(jldt,kldt)
    }


}

#endif //CALC_DT_H_
