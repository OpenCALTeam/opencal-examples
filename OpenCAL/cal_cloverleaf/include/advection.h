#ifndef _ADVECTION_H_
#define _ADVECTION_H_

#include<iostream>
#include<model.h>
#include<advec_mom.h>


//check loop boundaries
void advec_cell_kernel_xdir_sweep0(struct CALModel2D  * clover_model,  int i,  int j){
    int cx = FTNREF1D(i,field.x_min-2);
    int cy = FTNREF1D(j, field.y_min-2);
    if(!checkBoundaries(cx,cy,field.x_min-2,field.x_max+2,field.y_min-2,field.y_max+2))
        return;
    auto pre_vol = clover_sbs.work_array1;
    auto post_vol = clover_sbs.work_array2;

    const auto vol_j_k              = get(clover_sbs.volume,cx,cy);
    const auto vol_flux_x_j_k   = get(clover_sbs.vol_flux_x,cx,cy);
    const auto vol_flux_x_j1_k  = get(clover_sbs.vol_flux_x,cx+1,cy);

    const auto pre_vol_j_k =  vol_j_k + ( vol_flux_x_j1_k - vol_flux_x_j_k );
    const auto post_vol_j_k = vol_j_k;



    set(pre_vol, cx, cy, pre_vol_j_k);
    set(post_vol, cx, cy, post_vol_j_k);

}


//check loop boundaries
void advec_cell_kernel_xdir_sweep1(struct CALModel2D  * clover_model,  int i,  int j){
    int cx = FTNREF1D(i,field.x_min-2);
    int cy = FTNREF1D(j, field.y_min-2);
    if(!checkBoundaries(cx,cy,field.x_min-2,field.x_max+2,field.y_min-2,field.y_max+2))
        return;
    auto pre_vol = clover_sbs.work_array1;
    auto post_vol = clover_sbs.work_array2;

    const auto vol_j_k          = get(clover_sbs.volume,cx,cy);
    const auto vol_flux_x_j_k   = get(clover_sbs.vol_flux_x,cx,cy);
    const auto vol_flux_y_j_k   = get(clover_sbs.vol_flux_y,cx,cy);
    const auto vol_flux_x_j1_k  = get(clover_sbs.vol_flux_x,cx+1,cy);
    const auto vol_flux_y_j_k1  = get(clover_sbs.vol_flux_y,cx,cy+1);

    auto pre_vol_j_k =  vol_j_k + ( vol_flux_x_j1_k - vol_flux_x_j_k +  vol_flux_y_j_k1 - vol_flux_y_j_k  );

    auto post_vol_j_k =  pre_vol_j_k - (vol_flux_x_j1_k - vol_flux_x_j_k);




    //sets
    set(pre_vol, cx, cy, pre_vol_j_k);
    set(post_vol, cx, cy, post_vol_j_k);

}

//ricontrolla i set che vanno fatti a partire dalle variabili automatiche e const
//check loop boundaries
void advec_cell_kernel_xdir2(struct CALModel2D  * clover_model,  int i,  int j){
    int cx = FTNREF1D(i,field.x_min-2);
    int cy = FTNREF1D(j, field.y_min-2);

    if(!checkBoundaries(cx,cy,field.x_min,field.x_max+2,field.y_min,field.y_max))
        return;

    double sigma,sigmat, sigmav, sigmam, sigma3, sigma4, wind;
    double diffuw, diffdw, limiter;
    double one_by_six = 1.0/6.0;

    int x_max=field.x_max;

    int upwind,donor,downwind,dif;

    const auto vol_flux_x_j_k   = get(clover_sbs.vol_flux_x,cx,cy);

    if(vol_flux_x_j_k > 0.0) {
        upwind      = cx-2; //j-2
        donor       = cx-1; //j-1
        downwind    = cx; //j
        dif         = donor;
    }
    else{
        upwind      = std::min(cx+1,x_max+2);
        donor       =cx;
        downwind    = cx-1;
        dif         =upwind;
    }
    auto pre_vol = clover_sbs.work_array1;
    auto energ_flux= clover_sbs.work_array7;

    // sigmat=ABS(vol_flux_x(j,k))/pre_vol(donor,k)
    sigmat = fabs(vol_flux_x_j_k)/get(pre_vol,donor,cy);
    // sigmat=ABS(vol_flux_x(j,k))/pre_vol(donor,k)

    sigma3 = (1.0+sigmat)* (calGetMatrixElement(clover_sbs.vertexdx,1,cx,0) /  calGetMatrixElement(clover_sbs.vertexdx,1,dif,0));
    //sigma3 = (1.0 + sigmat)*(vertexdx[OPS_ACC3(0,0)]/vertexdx[OPS_ACC3(dif,0)]);
    sigma4=2.0-sigmat;

    sigma=sigmat;
    sigmav=sigmat;

    //diffuw=density1(donor,k)-density1(upwind,k)
    const auto dd = get(clover_sbs.density1,donor,cy);
    const auto cc = get(clover_sbs.density1,upwind,cy);
    auto ddd = dd+cc;
    diffuw = get(clover_sbs.density1,donor,cy) - get(clover_sbs.density1,upwind,cy);
    //diffdw=density1(downwind,k)-density1(donor,k)
    diffdw = get(clover_sbs.density1,downwind,cy) - get(clover_sbs.density1,donor,cy);
    wind=1.0;

    if(diffdw <= 0.0)
        wind=-1.0;
    if(diffuw*diffdw > 0.0){
        limiter = (1.0-sigmav)*wind* std::min( std::min(fabs(diffuw), fabs(diffdw)) , one_by_six*(sigma3*fabs(diffuw) + sigma4*fabs(diffdw)));
    }else
        limiter=0.0;

    //mass_flux_x(j,k)=vol_flux_x(j,k)*(density1(donor,k)+limiter)
    const auto density1_donor_k = get(clover_sbs.density1, donor, cy);


    auto new_mass_flux_x_j_k = vol_flux_x_j_k * (density1_donor_k+limiter);
   // printf("\t %d \t %d \t %.16e \t %.16e \t %d \t %.16e\n",cx,cy,vol_flux_x_j_k,density1_donor_k,donor,limiter);
    //std::cout<<cx<<"\t"<<cy<<"\t"<<vol_flux_x_j_k<<"\t"<<density1_donor_k<<"\t"<<donor<<"\t"<<limiter<<std::endl;


    const auto pre_vol_donor_k = get(pre_vol, donor, cy);
    //sigmam=ABS(mass_flux_x(j,k))/(density1(donor,k)*pre_vol(donor,k))
    sigmam = fabs(new_mass_flux_x_j_k)/(density1_donor_k*pre_vol_donor_k);

    // diffuw=energy1(donor,k)-energy1(upwind,k)
    diffuw=get(clover_sbs.energy1,donor,cy)-get(clover_sbs.energy1,upwind,cy);
    //diffdw=energy1(downwind,k)-energy1(donor,k)
    diffdw=get(clover_sbs.energy1,downwind,cy)-get(clover_sbs.energy1,donor,cy);
    wind=1.0;
    if(diffdw <= 0.0)
        wind=-1.0;
    if(diffuw*diffdw > 0.0)
        limiter = (1.0-sigmam)*wind* std::min( std::min(fabs(diffuw), fabs(diffdw)) , one_by_six*(sigma3*fabs(diffuw) + sigma4*fabs(diffdw)));
    else
        limiter=0.0;




    //sets here
    set(clover_sbs.mass_flux_x, cx, cy , new_mass_flux_x_j_k);
    //ener_flux(j,k)=mass_flux_x(j,k)*(energy1(donor,k)+limiter)
    set(energ_flux, cx, cy , new_mass_flux_x_j_k * (get(clover_sbs.energy1,donor,cy)+limiter) );

}




//check loop boundaries
void advec_cell_kernel_xdir3(struct CALModel2D  * clover_model,  int i,  int j){
    int cx = FTNREF1D(i,field.x_min-2);
    int cy = FTNREF1D(j, field.y_min-2);
    if(!checkBoundaries(cx,cy,field.x_min,field.x_max,field.y_min,field.y_max))
        return;

    auto pre_vol = clover_sbs.work_array1;
    auto ener_flux = clover_sbs.work_array7;

    //!!!in fortran sets are only present for density1 and energy1 matrices!!!!
    //pre_mass_s=density1(j,k)*pre_vol(j,k)
    auto new_pre_mass_j_k = get(clover_sbs.density1,cx,cy) * get(pre_vol,cx,cy);

    //post_mass_s=pre_mass_s+mass_flux_x(j,k)-mass_flux_x(j+1,k)
    auto new_post_mass_j_k = new_pre_mass_j_k + get(clover_sbs.mass_flux_x,cx,cy) - get(clover_sbs.mass_flux_x,cx+1,cy);

    //post_ener_s=(energy1(j,k)*pre_mass_s+ener_flux(j,k)-ener_flux(j+1,k))/post_mass_s
    auto new_post_energ_j_k = (get(clover_sbs.energy1,cx,cy) * new_pre_mass_j_k + get(ener_flux,cx,cy)-get(ener_flux,cx+1,cy))/new_post_mass_j_k;

    //advec_vol_s=pre_vol(j,k)+vol_flux_x(j,k)-vol_flux_x(j+1,k)
    auto new_advec_vol_s = get(pre_vol, cx, cy) + get(clover_sbs.vol_flux_x, cx, cy) - get(clover_sbs.vol_flux_x, cx+1, cy);

    //              density1(j,k)=post_mass_s/advec_vol_s
    set(clover_sbs.density1,cx,cy , new_post_mass_j_k/new_advec_vol_s);
    //              energy1(j,k)=post_ener_s
    set(clover_sbs.energy1,cx,cy , new_post_energ_j_k);


}


//vedere la ipsilon




//dir Y ------------------------------
//check loop boundaries
void advec_cell_kernel_ydir_sweep0(struct CALModel2D  * clover_model,  int i,  int j){
    int cx = FTNREF1D(i,field.x_min-2);
    int cy = FTNREF1D(j, field.y_min-2);
    if(!checkBoundaries(cx,cy,field.x_min-2,field.x_max+2,field.y_min-2,field.y_max+2))
        return;
    auto pre_vol = clover_sbs.work_array1;
    auto post_vol = clover_sbs.work_array2;

    const auto vol              = get(clover_sbs.volume,cx,cy);
    const auto vol_flux_y_j_k   = get(clover_sbs.vol_flux_y,cx,cy);
    const auto vol_flux_y_j_k1  = get(clover_sbs.vol_flux_y,cx,cy+1);

    const auto pre_vol_j_k =  vol + vol_flux_y_j_k1 - vol_flux_y_j_k;
    const auto post_vol_j_k = vol;




    set(pre_vol, cx, cy, pre_vol_j_k);
    set(post_vol, cx, cy, post_vol_j_k);

}

void advec_cell_kernel_ydir_sweep1(struct CALModel2D  * clover_model,  int i,  int j){
    int cx = FTNREF1D(i,field.x_min-2);
    int cy = FTNREF1D(j, field.y_min-2);
    if(!checkBoundaries(cx,cy,field.x_min-2,field.x_max+2,field.y_min-2,field.y_max+2))
        return;
    auto pre_vol = clover_sbs.work_array1;
    auto post_vol = clover_sbs.work_array2;

    const auto vol              = get(clover_sbs.volume,cx,cy);
    const auto vol_flux_x_j_k   = get(clover_sbs.vol_flux_x,cx,cy);
    const auto vol_flux_y_j_k   = get(clover_sbs.vol_flux_y,cx,cy);
    const auto vol_flux_y_j_k1  = get(clover_sbs.vol_flux_y,cx,cy+1);

    const auto pre_vol_j_k =  vol + ( vol_flux_y_j_k1 - vol_flux_y_j_k +  vol_flux_y_j_k1 - vol_flux_x_j_k  );
    const auto post_vol_j_k =  pre_vol_j_k - (vol_flux_y_j_k1 - vol_flux_y_j_k);

    //sets
    set(pre_vol, cx, cy, pre_vol_j_k);
    set(post_vol, cx, cy, post_vol_j_k);

}



//check loop boundaries
void advec_cell_kernel_ydir2(struct CALModel2D  * clover_model,  int i,  int j){
    int cx = FTNREF1D(i,field.x_min-2);
    int cy = FTNREF1D(j, field.y_min-2);
    if(!checkBoundaries(cx,cy,field.x_min,field.x_max,field.y_min,field.y_max+2))
        return;

    double sigma,sigmat, sigmav, sigmam, sigma3, sigma4, wind;
    double diffuw, diffdw, limiter;
    double one_by_six = 1.0/6.0;

    int y_max=field.y_max;

    int upwind,donor,downwind,dif;

    const auto vol_flux_y_j_k   = get(clover_sbs.vol_flux_y,cx,cy);

    if(vol_flux_y_j_k > 0.0) {
        upwind      = cy-2; //j-2
        donor       = cy-1; //j-1
        downwind    = cy; //j
        dif         = donor;
    }
    else{
        upwind      = std::min(cx+1,y_max+2);
        donor       = cy;
        downwind    = cy-1;
        dif         = upwind;
    }
    auto pre_vol = clover_sbs.work_array1;
    auto energ_flux= clover_sbs.work_array7;

    //  sigmat=ABS(vol_flux_y(j,k))/pre_vol(j,donor)
    sigmat = fabs(vol_flux_y_j_k)/get(pre_vol,cx,donor);
    // sigmat=ABS(vol_flux_x(j,k))/pre_vol(donor,k)

    sigma3 = fabs(1.0+sigmat)* (calGetMatrixElement(clover_sbs.vertexdy,1,cy,0) /  calGetMatrixElement(clover_sbs.vertexdy,1,dif,0));
    //sigma3 = (1.0 + sigmat)*(vertexdx[OPS_ACC3(0,0)]/vertexdx[OPS_ACC3(dif,0)]);
    sigma4=2.0-sigmat;

    sigma=sigmat;
    sigmav=sigmat;

    //diffuw=density1(donor,k)-density1(upwind,k)
    diffuw = get(clover_sbs.density1,cx,donor) - get(clover_sbs.density1,cx,upwind);
    //diffdw=density1(downwind,k)-density1(donor,k)
    diffdw = get(clover_sbs.density1,cx,downwind) - get(clover_sbs.density1,cx,donor);
    wind=1.0;

    if(diffdw <= 0.0)
        wind=-1.0;
    if(diffuw*diffdw > 0.0){
        limiter = (1.0-sigmav)*wind* std::min( std::min(fabs(diffuw), fabs(diffdw)) , one_by_six*(sigma3*fabs(diffuw) + sigma4*fabs(diffdw)));
    }else
        limiter=0.0;

    //mass_flux_x(j,k)=vol_flux_x(j,k)*(density1(donor,k)+limiter)
    const auto density1_j_donor = get(clover_sbs.density1, cx,donor);
    auto new_mass_flux_y = vol_flux_y_j_k * (density1_j_donor+limiter);

    const auto pre_vol_j_donor = get(pre_vol, cx, donor);
    //sigmam=ABS(mass_flux_x(j,k))/(density1(donor,k)*pre_vol(donor,k))
    sigmam = fabs(new_mass_flux_y/(density1_j_donor*pre_vol_j_donor));

    // diffuw=energy1(j,donor)-energy1(j,upwind)
    diffuw=get(clover_sbs.energy1,cx,donor)-get(clover_sbs.energy1,cx,upwind);
    //diffdw=energy1(j,downwind)-energy1(j,donor)
    diffdw=get(clover_sbs.energy1,downwind,cy)-get(clover_sbs.energy1,donor,cy);
    wind=1.0;
    if(diffdw <= 0.0)
        wind=-1.0;
    if(diffuw*diffdw > 0.0)
        limiter = (1.0-sigmam)*wind* std::min( std::min(fabs(diffuw), fabs(diffdw)) , one_by_six*(sigma3*fabs(diffuw) + sigma4*fabs(diffdw)));
    else
        limiter=0.0;

    //sets here
    set(clover_sbs.mass_flux_y, cx, cy , new_mass_flux_y);
    //ener_flux(j,k)=mass_flux_x(j,k)*(energy1(donor,k)+limiter)
    set(energ_flux, cx, cy , new_mass_flux_y * (get(clover_sbs.energy1,cx,donor)+limiter) );

}

//check loop boundaries
void advec_cell_kernel_ydir3(struct CALModel2D  * clover_model,  int i,  int j){
    int cx = FTNREF1D(i,field.x_min-2);
    int cy = FTNREF1D(j, field.y_min-2);
    if(!checkBoundaries(cx,cy,field.x_min,field.x_max,field.y_min,field.y_max))
        return;

    auto pre_vol = clover_sbs.work_array1;
    //auto post_vol = clover_sbs.work_array2;
    //auto pre_mass = clover_sbs.work_array3;
    //auto post_mass = clover_sbs.work_array4;
    //auto advec_vol = clover_sbs.work_array5;
    //auto post_ener= clover_sbs.work_array6;
    auto ener_flux = clover_sbs.work_array7;

    //!!!in fortran sets are only present for density1 and energy1 matrices!!!!
    //pre_mass_s=density1(j,k)*pre_vol(j,k)
    auto new_pre_mass_j_k = get(clover_sbs.density1,cx,cy) * get(pre_vol,cx,cy);

    //post_mass_s=pre_mass_s+mass_flux_x(j,k)-mass_flux_x(j+1,k)
    auto new_post_mass_j_k = new_pre_mass_j_k + get(clover_sbs.mass_flux_y,cx,cy) - get(clover_sbs.mass_flux_y,cx,cy+1);

    //post_ener_s=(energy1(j,k)*pre_mass_s+ener_flux(j,k)-ener_flux(j+1,k))/post_mass_s
    auto new_post_energ_j_k = (get(clover_sbs.energy1,cx,cy) * new_pre_mass_j_k + get(ener_flux,cx,cy)-get(ener_flux,cx,cy+1))/new_post_mass_j_k;

    //advec_vol_s=pre_vol(j,k)+vol_flux_x(j,k)-vol_flux_x(j+1,k)
    auto new_advec_vol_s = get(pre_vol, cx, cy) + get(clover_sbs.vol_flux_y, cx, cy) - get(clover_sbs.vol_flux_y, cx, cy+1);



    //              density1(j,k)=post_mass_s/advec_vol_s
    set(clover_sbs.density1,cx,cy , new_post_mass_j_k/new_advec_vol_s);
    //              energy1(j,k)=post_ener_s
    set(clover_sbs.energy1,cx,cy , new_post_energ_j_k);


}

void advectDriver(int direction, int sweep_number){
    int xvel,yvel;
    xvel = g_xdir;
    yvel = g_ydir;
    if(direction == g_xdir){

        if(sweep_number==1){
            calApplyElementaryProcess2D(clover_model, advec_cell_kernel_xdir_sweep1);
            calUpdate2D(clover_model);
        }
        else{
            calApplyElementaryProcess2D(clover_model, advec_cell_kernel_xdir_sweep0);
            calUpdate2D(clover_model);
        }

        calApplyElementaryProcess2D(clover_model, advec_cell_kernel_xdir2);
        calUpdate2D(clover_model);
        calApplyElementaryProcess2D(clover_model, advec_cell_kernel_xdir3);
        calUpdate2D(clover_model);
    }else{
        if(sweep_number==1)
            calApplyElementaryProcess2D(clover_model, advec_cell_kernel_ydir_sweep1);
        else
            calApplyElementaryProcess2D(clover_model, advec_cell_kernel_ydir_sweep0);

        calApplyElementaryProcess2D(clover_model, advec_cell_kernel_ydir2);
        calApplyElementaryProcess2D(clover_model, advec_cell_kernel_ydir3);
    }


}

void advection(){

    int sweep_number, direction;


    sweep_number = 1;
    if(advect_x) direction = g_xdir;
    else direction = g_ydir;

    int xvel = g_xdir;
    int yvel = g_ydir;

    //----updateHalo-----
    setFields(0);
    fields[FIELD::ENERGY1] =  fields[FIELD::DENSITY1] = fields[FIELD::VOL_FLUX_X] = fields[FIELD::VOL_FLUX_Y] = 1;
    FIELD_DEPTH=2;
    update_halo();
    //----updateHalo-----


    advectDriver(direction, sweep_number);

    exit(-1);
    //----updateHalo-----
    setFields(0);
    fields[FIELD::ENERGY1] =  fields[FIELD::DENSITY1] = fields[FIELD::XVEL1] = fields[FIELD::YVEL1] = fields[FIELD::MASS_FLUX_X] =fields[FIELD::MASS_FLUX_Y] = 1;
    FIELD_DEPTH=2;
    update_halo();
    //----updateHalo-----
    advection_mom(xvel,direction,sweep_number);
    advection_mom(yvel,direction,sweep_number);


    sweep_number = 2;
    if(advect_x) direction = g_xdir;
    else direction = g_ydir;

    advectDriver(direction,sweep_number);

    //----updateHalo-----
    setFields(0);
    fields[FIELD::ENERGY1] =  fields[FIELD::DENSITY1] = fields[FIELD::XVEL1] = fields[FIELD::YVEL1] = fields[FIELD::MASS_FLUX_X] =fields[FIELD::MASS_FLUX_Y] = 1;
    FIELD_DEPTH=2;
    update_halo();
    //----updateHalo-----

    advection_mom(xvel,direction,sweep_number);
    advection_mom(yvel,direction,sweep_number);

    calUpdate2D(clover_model);


}

#endif //_ADVECTION_H_

