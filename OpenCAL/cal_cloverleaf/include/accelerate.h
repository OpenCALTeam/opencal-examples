#ifndef _ACCELERATE_H_
#define _ACCELERATE_H_

#include<iostream>
#include<model.h>
void accelerate_kernel(struct CALModel2D  * clover_model,  int i,  int j){
    int cx = FTNREF1D(i,field.x_min-2);
    int cy = FTNREF1D(j, field.y_min-2);
    if(!checkBoundaries(cx,cy,field.x_min,field.x_max+1,field.y_min,field.y_max+1))
        return;



    double nodal_mass, stepbymass_s;
    const double halfdt=0.5*dt;

    //{0,0, -1,0, 0,-1, -1,-1};
    stepbymass_s = halfdt/((
                               get(clover_sbs.density0,cx-1,cy-1)  * get(clover_sbs.volume,cx-1,cy-1)
                               +
                               get(clover_sbs.density0,cx,cy-1)    * get(clover_sbs.volume,cx,cy-1)
                               +
                               get(clover_sbs.density0,cx,cy)     * get(clover_sbs.volume,cx,cy)
                               +
                               get(clover_sbs.density0,cx-1,cy)   * get(clover_sbs.volume,cx-1,cy)
                               )*0.25);


    printf("\t%d\t%d\t%.16e\t%.16e\n",cx,cy,stepbymass_s,halfdt);
    //xvel1(j,k)=xvel0(j,k)-stepbymass_s*(xarea(j  ,k  )*(pressure(j  ,k  )-pressure(j-1,k  ))    &
    //        +xarea(j  ,k-1)*(pressure(j  ,k-1)-pressure(j-1,k-1)))

    auto pressure_diff1 = get(clover_sbs.pressure,cx,cy) - get(clover_sbs.pressure,cx-1,cy);
    auto pressure_diff2 = get(clover_sbs.pressure,cx,cy-1) - get(clover_sbs.pressure,cx-1,cy-1);

    auto new_xvel1 = get(clover_sbs.xvel0,cx,cy) - stepbymass_s *    (
                get(clover_sbs.xarea,cx,cy)* pressure_diff1 + get(clover_sbs.xarea,cx,cy-1)* pressure_diff2
                );


    //yvel1(j,k)=yvel0(j,k)-stepbymass_s*(yarea(j  ,k  )*(pressure(j  ,k  )-pressure(j  ,k-1))    &
    //+yarea(j-1,k  )*(pressure(j-1,k  )-pressure(j-1,k-1)))
    pressure_diff1 = get(clover_sbs.pressure,cx,cy) - get(clover_sbs.pressure,cx,cy-1);
    pressure_diff2 = get(clover_sbs.pressure,cx-1,cy) - get(clover_sbs.pressure,cx-1,cy-1);

    auto new_yvel1 = get(clover_sbs.yvel0,cx,cy) - stepbymass_s *    (
                get(clover_sbs.yarea,cx,cy)* pressure_diff1 + get(clover_sbs.yarea,cx-1,cy)* pressure_diff2
                );

    //printf("\t%d\t%d\t%.16e\n",cx,cy,pressure_diff2);

    //xvel1(j,k)=xvel1(j,k)-stepbymass_s*(xarea(j  ,k  )*(viscosity(j  ,k  )-viscosity(j-1,k  )) &
    //+xarea(j  ,k-1)*(viscosity(j  ,k-1)-viscosity(j-1,k-1)))

    auto viscosity_diff1 = get(clover_sbs.viscosity,cx,cy) - get(clover_sbs.viscosity,cx-1,cy);
    auto viscosity_diff2 = get(clover_sbs.viscosity,cx,cy-1) - get(clover_sbs.viscosity,cx-1,cy-1);

    new_xvel1 = new_xvel1 - stepbymass_s *    (
                get(clover_sbs.xarea,cx,cy)* viscosity_diff1 + get(clover_sbs.xarea,cx,cy-1)* viscosity_diff2
                );


    //       vel1(j,k)=yvel1(j,k)-stepbymass_s*(yarea(j  ,k  )*(viscosity(j  ,k  )-viscosity(j  ,k-1)) &
    //       +yarea(j-1,k  )*(viscosity(j-1,k  )-viscosity(j-1,k-1)))

    viscosity_diff1 = get(clover_sbs.viscosity,cx,cy) - get(clover_sbs.viscosity,cx,cy-1);
    viscosity_diff2 = get(clover_sbs.viscosity,cx-1,cy) - get(clover_sbs.viscosity,cx-1,cy-1);

    new_yvel1 = new_yvel1 - stepbymass_s *    (
                get(clover_sbs.yarea,cx,cy)* viscosity_diff1 + get(clover_sbs.yarea,cx-1,cy)* viscosity_diff2
                );


    //printf("\t%d\t%d\t%.16e\t%.16e\n",cx,cy,new_xvel1,new_yvel1);
    set(clover_sbs.xvel1,cx,cy, new_xvel1);
    set(clover_sbs.yvel1,cx,cy,new_yvel1);

}

void accelerate(){


    calApplyElementaryProcess2D(clover_model, accelerate_kernel);

    //aggiorna i sottostati ceh hai upppppppddddattato
    calUpdate2D(clover_model);

}


#endif //_ACCELERATE_H_
