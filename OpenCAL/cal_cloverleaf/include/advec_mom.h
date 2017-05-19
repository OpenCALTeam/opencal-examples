#ifndef _ADVECTION_MOM_H_
#define _ADVECTION_MOM_H_

#include<iostream>
#include<model.h>



int advec_mom_direction=0;


//probabilmente devi fare +4
void advec_cell_kernel_x1(struct CALModel2D  * clover_model,  int i,  int j){
    int cx = FTNREF1D(i,field.x_min-2);
    int cy = FTNREF1D(j, field.y_min-2);
    if(!checkBoundaries(cx,cy,field.x_min,field.x_max+2,field.y_min,field.y_max+2))
        return;

    auto pre_vol = clover_sbs.work_array5;
    auto post_vol = clover_sbs.work_array6;

    const auto vol_j_k           = get(clover_sbs.volume,cx,cy);
    const auto vol_flux_x_j1_k   = get(clover_sbs.vol_flux_x,cx+1,cy);
    const auto vol_flux_x_j_k   = get(clover_sbs.vol_flux_x,cx,cy);
    const auto vol_flux_y_j_k   = get(clover_sbs.vol_flux_y,cx,cy);
    const auto vol_flux_y_j_k1  = get(clover_sbs.vol_flux_y,cx,cy+1);

    const auto post_vol_j_k = vol_j_k + vol_flux_y_j_k1 - vol_flux_y_j_k;
    const auto pre_vol_j_k =  post_vol_j_k +  vol_flux_x_j1_k - vol_flux_x_j_k ;


    set(pre_vol, cx, cy, pre_vol_j_k);
    set(post_vol, cx, cy, post_vol_j_k);

}

void advec_cell_kernel_y1(struct CALModel2D  * clover_model,  int i,  int j){
    int cx = FTNREF1D(i,field.x_min-2);
    int cy = FTNREF1D(j, field.y_min-2);
    if(!checkBoundaries(cx,cy,field.x_min,field.x_max+2,field.y_min,field.y_max+2))
        return;

    auto pre_vol = clover_sbs.work_array5;
    auto post_vol = clover_sbs.work_array6;

    const auto vol_j_k           = get(clover_sbs.volume,cx,cy);
    const auto vol_flux_x_j1_k   = get(clover_sbs.vol_flux_x,cx+1,cy);
    const auto vol_flux_x_j_k   = get(clover_sbs.vol_flux_x,cx,cy);
    const auto vol_flux_y_j_k   = get(clover_sbs.vol_flux_y,cx,cy);
    const auto vol_flux_y_j_k1  = get(clover_sbs.vol_flux_y,cx,cy+1);

    const auto post_vol_j_k = vol_j_k + vol_flux_x_j1_k - vol_flux_x_j_k;
    const auto pre_vol_j_k =  post_vol_j_k + ( vol_flux_y_j_k1 - vol_flux_y_j_k );


    set(pre_vol, cx, cy, pre_vol_j_k);
    set(post_vol, cx, cy, post_vol_j_k);
}


void advec_cell_kernel_x2(struct CALModel2D  * clover_model,  int i,  int j){
    int cx = FTNREF1D(i,field.x_min-2);
    int cy = FTNREF1D(j, field.y_min-2);
    if(!checkBoundaries(cx,cy,field.x_min,field.x_max+2,field.y_min,field.y_max+2))
        return;

    auto pre_vol = clover_sbs.work_array5;
    auto post_vol = clover_sbs.work_array6;


    const auto vol_j_k          = get(clover_sbs.volume,cx,cy);
    const auto vol_flux_y_j_k   = get(clover_sbs.vol_flux_y,cx,cy);
    const auto vol_flux_y_j_k1  = get(clover_sbs.vol_flux_y,cx,cy+1);

    const auto post_vol_j_k = vol_j_k;
    const auto pre_vol_j_k =  post_vol_j_k + ( vol_flux_y_j_k1 - vol_flux_y_j_k );


    set(pre_vol, cx, cy, pre_vol_j_k);
    set(post_vol, cx, cy, post_vol_j_k);

}

void advec_cell_kernel_y2(struct CALModel2D  * clover_model,  int i,  int j){
    int cx = FTNREF1D(i,field.x_min-2);
    int cy = FTNREF1D(j, field.y_min-2);
    if(!checkBoundaries(cx,cy,field.x_min,field.x_max+2,field.y_min,field.y_max+2))
        return;


    auto pre_vol = clover_sbs.work_array5;
    auto post_vol = clover_sbs.work_array6;


    const auto vol_j_k          = get(clover_sbs.volume,cx,cy);
    const auto vol_flux_x_j_k   = get(clover_sbs.vol_flux_x,cx,cy);
    const auto vol_flux_x_j1_k  = get(clover_sbs.vol_flux_x,cx+1,cy);

    const auto post_vol_j_k = vol_j_k;
    const auto pre_vol_j_k =  post_vol_j_k + ( vol_flux_x_j1_k - vol_flux_x_j_k );


    set(pre_vol, cx, cy, pre_vol_j_k);
    set(post_vol, cx, cy, post_vol_j_k);


}

//-----------------------------
//IF(direction.EQ.1)THEN
//  IF(which_vel.EQ.1)THEN
void advec_cell_kernel_dir1_vel_1(struct CALModel2D  * clover_model,  int i,  int j){
    int cx = FTNREF1D(i,field.x_min-2);
    int cy = FTNREF1D(j, field.y_min-2);
    if(!checkBoundaries(cx,cy,field.x_min,field.x_max+2,field.y_min,field.y_max+1))
        return;

    auto node_flux = clover_sbs.work_array1;


    const auto mass_flux_x_j_km1  = get(clover_sbs.mass_flux_x,cx,cy-1);
    const auto mass_flux_x_j_k  = get(clover_sbs.mass_flux_x,cx,cy);
    const auto mass_flux_x_j1_km1  = get(clover_sbs.mass_flux_x,cx+1,cy-1);
    const auto mass_flux_x_j1_k  = get(clover_sbs.mass_flux_x,cx+1,cy);

    const auto new_node_flux_j_k = 0.25*
            ( mass_flux_x_j_km1 + mass_flux_x_j_k + mass_flux_x_j1_km1 + mass_flux_x_j1_k  );

    set(node_flux, cx, cy, new_node_flux_j_k);

}

void advec_cell_kernel_dir1_vel_2(struct CALModel2D  * clover_model,  int i,  int j){
    int cx = FTNREF1D(i,field.x_min-2);
    int cy = FTNREF1D(j, field.y_min-2);
    if(!checkBoundaries(cx,cy,field.x_min,field.x_max+2,field.y_min,field.y_max+1))
        return;
    auto node_mass_post = clover_sbs.work_array2;
    auto node_mass_pre  = clover_sbs.work_array3;
    auto post_vol = clover_sbs.work_array6;
    auto node_flux = clover_sbs.work_array1;

    const auto density1_j_km1       = get(clover_sbs.density1,cx,cy-1);
    const auto density1_j_k         = get(clover_sbs.density1,cx,cy);
    const auto density1_jm1_km1     = get(clover_sbs.density1,cx-1,cy-1);
    const auto density1_jm1_k       = get(clover_sbs.density1,cx-1,cy);

    const auto post_vol_j_km1       = get(post_vol,cx,cy-1);
    const auto post_vol_j_k         = get(post_vol,cx,cy);
    const auto post_vol_jm1_km1     = get(post_vol,cx-1,cy-1);
    const auto post_vol_jm1_k       = get(post_vol,cx-1,cy);

    const auto new_node_mass_post_j_k = 0.25* (
                density1_j_km1 * post_vol_j_km1 +
                density1_j_k*  post_vol_j_k+
                density1_jm1_km1*post_vol_jm1_km1+
                density1_jm1_k*post_vol_jm1_k
                );

    const auto node_flux_jm1_k = get(node_flux, cx-1,cy);
    const auto node_flux_j_k = get(node_flux, cx,cy);

    const auto new_node_mass_pre_j_k = new_node_mass_post_j_k - node_flux_jm1_k + node_flux_j_k;

    set(node_mass_post, cx, cy , new_node_mass_post_j_k);
    set(node_mass_pre, cx, cy, new_node_mass_pre_j_k);
}

void advec_cell_kernel_dir1_vel_3(struct CALModel2D  * clover_model,  int i,  int j){
    int cx = FTNREF1D(i,field.x_min-2);
    int cy = FTNREF1D(j, field.y_min-2);
    if(!checkBoundaries(cx,cy,field.x_min,field.x_max+1,field.y_min,field.y_max+2))
        return;
    auto node_mass_post = clover_sbs.work_array2;
    auto node_mass_pre  = clover_sbs.work_array3;
    auto post_vol = clover_sbs.work_array6;
    auto node_flux = clover_sbs.work_array1;

    const auto density1_j_km1       = get(clover_sbs.density1,cx,cy-1);
    const auto density1_j_k         = get(clover_sbs.density1,cx,cy);
    const auto density1_jm1_km1     = get(clover_sbs.density1,cx-1,cy-1);
    const auto density1_jm1_k       = get(clover_sbs.density1,cx-1,cy);

    const auto post_vol_j_km1       = get(post_vol,cx,cy-1);
    const auto post_vol_j_k         = get(post_vol,cx,cy);
    const auto post_vol_jm1_km1     = get(post_vol,cx-1,cy-1);
    const auto post_vol_jm1_k       = get(post_vol,cx-1,cy);

    const auto new_node_mass_post_j_k = 0.25* (
                density1_j_km1 * post_vol_j_km1 +
                density1_j_k*  post_vol_j_k+
                density1_jm1_km1*post_vol_jm1_km1+
                density1_jm1_k*post_vol_jm1_k
                );

    const auto node_flux_j_km1 = get(node_flux, cx,cy-1);
    const auto node_flux_j_k = get(node_flux, cx,cy);

    const auto new_node_mass_pre_j_k = new_node_mass_post_j_k - node_flux_j_km1 + node_flux_j_k;

    set(node_mass_post, cx, cy , new_node_mass_post_j_k);
    set(node_mass_pre, cx, cy, new_node_mass_pre_j_k);

}



void advec_cell_kernel_dir1_vel_4(struct CALModel2D  * clover_model,  int i,  int j){
    int cx = FTNREF1D(i,field.x_min-2);
    int cy = FTNREF1D(j, field.y_min-2);
    if(!checkBoundaries(cx,cy,field.x_min,field.x_max+1,field.y_min,field.y_max+1))
        return;

    auto node_flux = clover_sbs.work_array1;
    auto node_mass_pre  = clover_sbs.work_array3;
    auto mom_flux = clover_sbs.work_array4;

    auto vel1 = advec_mom_direction == 1 ? clover_sbs.xvel1 : clover_sbs.yvel1 ;

    int upwind,donor,downwind,dif;
    double limiter,sigma,wind, width,vdiffuw,vdiffdw,auw,adw,advec_vel_s;



    const auto node_flux_j_k = get(node_flux, cx, cy);
    if(node_flux_j_k < 0.0)
    {
        upwind      = cx+2; //j+2
        donor       = cx+1; //j+1
        downwind    = cx; //j
        dif         = donor;
    }else
    {
        upwind      = cx-1; //j-1
        donor       = cx; //j
        downwind    = cx+1; //j+1
        dif         = upwind;
    }


    const auto node_mass_pre_donor_k = get(node_mass_pre, donor, cy);

    sigma = std::fabs(node_flux_j_k/node_mass_pre_donor_k);
    width =calGetMatrixElement(clover_sbs.celldx,1,cx,0);

    const auto vel1_donor_k = get(vel1, donor,cy);
    const auto vel1_upwind_k = get(vel1, upwind,cy);
    const auto vel1_downwind_k = get(vel1, downwind,cy);


    vdiffuw=vel1_donor_k-vel1_upwind_k;
    vdiffdw=vel1_downwind_k-vel1_donor_k;
    limiter=0.0;

    if(vdiffuw*vdiffuw > 0.0){
        auw= std::fabs(vdiffuw);
        adw=std::fabs(vdiffdw);
        wind=1.0;
        if(vdiffdw <= 0.0)
            wind=-1.0;

        const auto tmp = width*((2.0-sigma)*adw/width+(1.0+sigma)*auw/calGetMatrixElement(clover_sbs.celldx,1,dif,0))/6.0;
        limiter = wind * (std::min( tmp ,std::min(auw,adw)));
    }


    advec_vel_s = vel1_donor_k+(1.0-sigma)*limiter;
    set(mom_flux , cx, cy, advec_vel_s*node_flux_j_k);

}

void advec_cell_kernel_dir1_vel_5(struct CALModel2D  * clover_model,  int i,  int j){
    int cx = FTNREF1D(i,field.x_min-2);
    int cy = FTNREF1D(j, field.y_min-2);
    if(!checkBoundaries(cx,cy,field.x_min,field.x_max+1,field.y_min,field.y_max+2))
        return;


    auto vel1 = advec_mom_direction == 1 ? clover_sbs.xvel1 : clover_sbs.yvel1 ;

    const  auto node_mass_post = clover_sbs.work_array2;
    const auto node_mass_pre  = clover_sbs.work_array3;
    const auto mom_flux = clover_sbs.work_array4;

    const  auto vel1_j_k = get(vel1,cx,cy);
    const auto node_mass_pre_j_k = get(node_mass_pre,cx,cy);
    const auto node_mass_post_j_k = get(node_mass_post,cx,cy);
    const auto mom_flux_j_k = get(mom_flux,cx,cy);
    const auto mom_flux_jm1_k = get(mom_flux,cx-1,cy);;



    const auto new_vel1_j_k = (vel1_j_k * node_mass_pre_j_k+mom_flux_jm1_k-mom_flux_j_k)/node_mass_post_j_k;
    set(vel1,cx,cy, new_vel1_j_k);

}


//-----------------------------



//-----------------------------
//ELSEIF(direction.EQ.2)THEN
//     IF(which_vel.EQ.1)THEN
void advec_cell_kernel_dir2_vel_1(struct CALModel2D  * clover_model,  int i,  int j){
    int cx = FTNREF1D(i,field.x_min-2);
    int cy = FTNREF1D(j, field.y_min-2);
    if(!checkBoundaries(cx,cy,field.x_min,field.x_max+1,field.y_min,field.y_max+2))
        return;

    auto node_flux = clover_sbs.work_array1;


    const auto mass_flux_y_jm1_k    = get(clover_sbs.mass_flux_y,cx-1,cy);
    const auto mass_flux_y_j_k      = get(clover_sbs.mass_flux_y,cx,cy);
    const auto mass_flux_y_jm1_k1   = get(clover_sbs.mass_flux_y,cx-1,cy+1);
    const auto mass_flux_y_j_k1     = get(clover_sbs.mass_flux_y,cx,cy+1);

    const auto new_node_flux_j_k = 0.25*
            ( mass_flux_y_jm1_k + mass_flux_y_j_k + mass_flux_y_jm1_k1 + mass_flux_y_j_k1  );

    set(node_flux, cx, cy, new_node_flux_j_k);

}

void advec_cell_kernel_dir2_vel_2(struct CALModel2D  * clover_model,  int i,  int j){

    int cx = FTNREF1D(i,field.x_min-2);
    int cy = FTNREF1D(j, field.y_min-2);
    if(!checkBoundaries(cx,cy,field.x_min,field.x_max+1,field.y_min,field.y_max+2))
        return;

    auto node_mass_post = clover_sbs.work_array2;
    auto node_mass_pre  = clover_sbs.work_array3;
    auto post_vol = clover_sbs.work_array6;
    auto node_flux = clover_sbs.work_array1;

    const auto density1_j_km1       = get(clover_sbs.density1,cx,cy-1);
    const auto density1_j_k         = get(clover_sbs.density1,cx,cy);
    const auto density1_jm1_km1     = get(clover_sbs.density1,cx-1,cy-1);
    const auto density1_jm1_k       = get(clover_sbs.density1,cx-1,cy);

    const auto post_vol_j_km1       = get(post_vol,cx,cy-1);
    const auto post_vol_j_k         = get(post_vol,cx,cy);
    const auto post_vol_jm1_km1     = get(post_vol,cx-1,cy-1);
    const auto post_vol_jm1_k       = get(post_vol,cx-1,cy);

    const auto new_node_mass_post_j_k = 0.25* (
                density1_j_km1 * post_vol_j_km1 +
                density1_j_k*  post_vol_j_k+
                density1_jm1_km1*post_vol_jm1_km1+
                density1_jm1_k*post_vol_jm1_k
                );

    const auto node_flux_j_km1 = get(node_flux, cx,cy-1);
    const auto node_flux_j_k = get(node_flux, cx,cy);

    const auto new_node_mass_pre_j_k = new_node_mass_post_j_k - node_flux_j_km1 + node_flux_j_k;

    set(node_mass_post, cx, cy , new_node_mass_post_j_k);
    set(node_mass_pre, cx, cy, new_node_mass_pre_j_k);

}

void advec_cell_kernel_dir2_vel_3(struct CALModel2D  * clover_model,  int i,  int j){
    int cx = FTNREF1D(i,field.x_min-2);
    int cy = FTNREF1D(j, field.y_min-2);
    if(!checkBoundaries(cx,cy,field.x_min,field.x_max+1,field.y_min,field.y_max+1))
        return;

    auto node_flux = clover_sbs.work_array1;
    auto node_mass_pre  = clover_sbs.work_array3;
    auto mom_flux = clover_sbs.work_array4;

    auto vel1 = advec_mom_direction == 1 ? clover_sbs.xvel1 : clover_sbs.yvel1 ;

    int upwind,donor,downwind,dif;
    double limiter,sigma,wind, width,vdiffuw,vdiffdw,auw,adw,advec_vel_s;



    const auto node_flux_j_k = get(node_flux, cx, cy);
    if(node_flux_j_k < 0.0)
    {
        upwind      = cy+2; //k+2
        donor       = cy+1; //k+1
        downwind    = cy; //j
        dif         = donor;
    }else
    {
        upwind      = cy-1; //k-1
        donor       = cy; //k
        downwind    = cy+1; //k+1
        dif         = upwind;
    }


    const auto node_mass_pre_donor_k = get(node_mass_pre, donor, cy);

    sigma = std::fabs(node_flux_j_k/node_mass_pre_donor_k);
    width =calGetMatrixElement(clover_sbs.celldx,1,cx,0);

    const auto vel1_donor_k = get(vel1, donor,cy);
    const auto vel1_upwind_k = get(vel1, upwind,cy);
    const auto vel1_downwind_k = get(vel1, downwind,cy);


    vdiffuw=vel1_donor_k-vel1_upwind_k;
    vdiffdw=vel1_downwind_k-vel1_donor_k;
    limiter=0.0;

    if(vdiffuw*vdiffuw > 0.0){
        auw= std::fabs(vdiffuw);
        adw=std::fabs(vdiffdw);
        wind=1.0;
        if(vdiffdw <= 0.0)
            wind=-1.0;

        const auto tmp = width*((2.0-sigma)*adw/width+(1.0+sigma)*auw/calGetMatrixElement(clover_sbs.celldx,1,dif,0))/6.0;
        limiter = wind * (std::min( tmp ,std::min(auw,adw)));
    }


    advec_vel_s = vel1_donor_k+(1.0-sigma)*limiter;
    set(mom_flux , cx, cy, advec_vel_s*node_flux_j_k);

}

void advec_cell_kernel_dir2_vel_4(struct CALModel2D  * clover_model,  int i,  int j){
    int cx = FTNREF1D(i,field.x_min-2);
    int cy = FTNREF1D(j, field.y_min-2);
    if(!checkBoundaries(cx,cy,field.x_min,field.x_max+1,field.y_min,field.y_max+1))
        return;


    auto vel1 = advec_mom_direction == 1 ? clover_sbs.xvel1 : clover_sbs.yvel1 ;

    const  auto node_mass_post = clover_sbs.work_array2;
    const auto node_mass_pre  = clover_sbs.work_array3;
    const auto mom_flux = clover_sbs.work_array4;

    const  auto vel1_j_k = get(vel1,cx,cy);
    const auto node_mass_pre_j_k = get(node_mass_pre,cx,cy);
    const auto node_mass_post_j_k = get(node_mass_post,cx,cy);
    const auto mom_flux_j_k = get(mom_flux,cx,cy);
    const auto mom_flux_j_km1 = get(mom_flux,cx,cy-1);



    const auto new_vel1_j_k = (vel1_j_k * node_mass_pre_j_k+mom_flux_j_km1-mom_flux_j_k)/node_mass_post_j_k;
    set(vel1,cx,cy, new_vel1_j_k);

}


//-----------------------------

void advection_mom(const int which_vel, const int direction, const int sweep_number){

    int mom_sweep = direction+2*(sweep_number-1);

    switch (mom_sweep) {
    case 1:
        calApplyElementaryProcess2D(clover_model, advec_cell_kernel_x1);
        break;
    case 2:
        calApplyElementaryProcess2D(clover_model, advec_cell_kernel_y1);
        break;
    case 3:
        calApplyElementaryProcess2D(clover_model, advec_cell_kernel_x2);
        break;
    case 4:

        calApplyElementaryProcess2D(clover_model, advec_cell_kernel_y2);
        break;

    }

    if(direction==1){
        if(which_vel ==1){
            calApplyElementaryProcess2D(clover_model,advec_cell_kernel_dir1_vel_1);
            calApplyElementaryProcess2D(clover_model,advec_cell_kernel_dir1_vel_2);
            calApplyElementaryProcess2D(clover_model,advec_cell_kernel_dir1_vel_3);
            calApplyElementaryProcess2D(clover_model,advec_cell_kernel_dir1_vel_4);
            calApplyElementaryProcess2D(clover_model,advec_cell_kernel_dir1_vel_5);



        }
    }
    if(direction==2){
        if(which_vel ==1){
            calApplyElementaryProcess2D(clover_model,advec_cell_kernel_dir2_vel_1);
            calApplyElementaryProcess2D(clover_model,advec_cell_kernel_dir2_vel_2);
            calApplyElementaryProcess2D(clover_model,advec_cell_kernel_dir2_vel_3);
            calApplyElementaryProcess2D(clover_model,advec_cell_kernel_dir2_vel_4);

        }
    }


}

#endif //_ADVECTION_MOM_H_

