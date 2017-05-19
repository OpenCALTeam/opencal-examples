#ifndef _CALC_FLUX_H_
#define _CALC_FLUX_H_

#include<iostream>
#include<model.h>





void calc_flux_kernel(struct CALModel2D  * clover_model,  int i,  int j){
    int cx = FTNREF1D(i,field.x_min-2);
    int cy = FTNREF1D(j, field.y_min-2);
    if(!checkBoundaries(cx,cy,field.x_min,field.x_max+1,field.y_min,field.y_max+1))
        return;



    //    vol_flux_x(j,k)=0.25_8*dt*xarea(j,k)                  &
    //    *(xvel0(j,k)+xvel0(j,k+1)+xvel1(j,k)+xvel1(j,k+1))
    const auto mult = 0.25*dt;
    const auto vol_flux_x_j_k = mult * get(clover_sbs.xarea,cx,cy) *
            (
                get(clover_sbs.xvel0,cx,cy)+get(clover_sbs.xvel0,cx,cy+1)+
                get(clover_sbs.xvel1,cx,cy)+get(clover_sbs.xvel1,cx,cy+1)
                );

    //        vol_flux_y(j,k)=0.25_8*dt*yarea(j,k)                  &
    //          *(yvel0(j,k)+yvel0(j+1,k)+yvel1(j,k)+yvel1(j+1,k))
    const auto vol_flux_y_j_k = mult * get(clover_sbs.yarea,cx,cy) *
            (
                get(clover_sbs.yvel0,cx,cy)+get(clover_sbs.yvel0,cx+1,cy)+
                get(clover_sbs.yvel1,cx,cy)+get(clover_sbs.yvel1,cx+1,cy)
                );



    set(clover_sbs.vol_flux_x,cx,cy,vol_flux_x_j_k);
    set(clover_sbs.vol_flux_y,cx,cy,vol_flux_y_j_k);

}

void calc_flux(){


    calApplyElementaryProcess2D(clover_model, calc_flux_kernel);

    calUpdateSubstate2Dr(clover_model,clover_sbs.vol_flux_x);
    calUpdateSubstate2Dr(clover_model,clover_sbs.vol_flux_y);


}


#endif //_CALC_FLUX_H_
