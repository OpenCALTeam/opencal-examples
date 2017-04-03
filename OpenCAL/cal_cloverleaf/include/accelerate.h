#ifndef _ACCELERATE_H_
#define _ACCELERATE_H_

#include<iostream>
#include<model.h>
void accelerate_kernel(struct CALModel2D  * clover_model,  int i,  int j){
    int cx = FTNREF1D(i,field.x_min-2);
    int cy = FTNREF1D(j, field.y_min-2);
    if(!checkBoundaries(cx,cy,field.x_min,field.x_max,field.y_min,field.y_max))
        return;

    double nodal_mass;
    const double halfdt=0.5*dt;

    //{0,0, -1,0, 0,-1, -1,-1};
    stepbymass = halfdt/((get(clover_sbs.density0,cx-1,cy-1)*get(clover_sbs.density0,cx-1,cy-1)*
                          get(clover_sbs.density0,cx,cy-1)*get(clover_sbs.density0,cx,cy-1)*
                          get(clover_sbs.density0,cx,cy)*get(clover_sbs.density0,cx,cy)*
                          get(clover_sbs.density0,cx-1,cy)*get(clover_sbs.density0,cx-1,cy)
                          )*0.25);


    xvel1[OPS_ACC4(0,0)] = xvel0[OPS_ACC3(0,0)] - stepbymass[OPS_ACC2(0,0)] *
            ( xarea[OPS_ACC5(0,0)]  * ( pressure[OPS_ACC6(0,0)] - pressure[OPS_ACC6(-1,0)] ) +
            xarea[OPS_ACC5(0,-1)] * ( pressure[OPS_ACC6(0,-1)] - pressure[OPS_ACC6(-1,-1)] ) );

}

void accelerate(bool predict){


    calApplyElementaryProcess2D(clover_model, accelerate_kernel);

    //aggiorna i sottostati ceh hai upppppppddddattato

}


#endif //_ACCELERATE_H_
