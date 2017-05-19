#ifndef _PDV_H_
#define _PDV_H_

#include<iostream>
#include<model.h>
#include<functional>
#include<update_halo.h>
#include<ideal_gas.h>

bool predict = false;
/*void PdV_kernel_PREDICT(struct CALModel2D  * clover_model,  int i,  int j){
    int cx = FTNREF1D(i,field.x_min-2);
    int cy = FTNREF1D(j, field.y_min-2);
    if(!checkBoundaries(cx,cy,field.x_min,field.x_max,field.y_min,field.y_max))
        return;

    double recip_volume, energy_change;//, min_cell_volume;
    double right_flux, left_flux, top_flux, bottom_flux, total_flux, volume_change_s;


    if(cx==481 && cy==2){
        printf("blablabla\n");
    }
    //          left_flux=  (xarea(j  ,k  )*(xvel0(j  ,k  )+xvel0(j  ,k+1)                     &
    //+xvel0(j  ,k  )+xvel0(j  ,k+1)))*0.25_8*dt*0.5
    left_flux= (
                get(clover_sbs.xarea,cx,cy)* (
                    get(clover_sbs.xvel0,cx,cy) + get(clover_sbs.xvel0,cx,cy+1) + get(clover_sbs.xvel0,cx,cy) + get(clover_sbs.xvel0,cx,cy+1))
                )*0.25*dt*0.5;

    //right_flux= (xarea(j+1,k  )*(xvel0(j+1,k  )+xvel0(j+1,k+1)                     &
    //   +xvel0(j+1,k  )+xvel0(j+1,k+1)))*0.25_8*dt*0.5
    right_flux= (
                get(clover_sbs.xarea,cx+1,cy)* (
                    get(clover_sbs.xvel0,cx+1,cy) + get(clover_sbs.xvel0,cx+1,cy+1) + get(clover_sbs.xvel0,cx+1,cy) + get(clover_sbs.xvel0,cx+1,cy+1))
                )*0.25*dt*0.5;

    //   bottom_flux=(yarea(j  ,k  )*(yvel0(j  ,k  )+yvel0(j+1,k  )                     &
    //+yvel0(j  ,k  )+yvel0(j+1,k  )))*0.25_8*dt*0.5
    bottom_flux = (
                get(clover_sbs.yarea,cx,cy)* (
                    get(clover_sbs.yvel0,cx,cy) + get(clover_sbs.yvel0,cx+1,cy) + get(clover_sbs.yvel0,cx,cy) + get(clover_sbs.yvel0,cx+1,cy))
                )*0.25*dt*0.5;
    //  top_flux=   (yarea(j  ,k+1)*(yvel0(j  ,k+1)+yvel0(j+1,k+1)                     &
    //+yvel0(j  ,k+1)+yvel0(j+1,k+1)))*0.25_8*dt*0.5
    bottom_flux = (
                get(clover_sbs.yarea,cx,cy+1)* (
                    get(clover_sbs.yvel0,cx,cy+1) + get(clover_sbs.yvel0,cx+1,cy+1) + get(clover_sbs.yvel0,cx,cy+1) + get(clover_sbs.yvel0,cx+1,cy+1))
                )*0.25*dt*0.5;

    total_flux=right_flux-left_flux+top_flux-bottom_flux;


    const auto _vol = get(clover_sbs.volume,cx,cy);

    volume_change_s = _vol/(_vol+total_flux);



    //     min_cell_volume=MIN(volume(j,k)+right_flux-left_flux+top_flux-bottom_flux &
    //       ,volume(j,k)+right_flux-left_flux                      &
    //       ,volume(j,k)+top_flux-bottom_flux)

    recip_volume = 1.0/_vol;

    const auto _press = get(clover_sbs.pressure,cx,cy);
    const auto _den = get(clover_sbs.density0,cx,cy);
    const auto _visc = get(clover_sbs.viscosity,cx,cy);

    energy_change = (_press/_den + _visc/_den)*total_flux*recip_volume;

    const auto _enrg = get(clover_sbs.energy0,cx,cy);
    set(clover_sbs.energy1,cx,cy,_enrg-energy_change);

    set(clover_sbs.density1,cx,cy,_den*volume_change_s);

}*/


void PdV_kernel(struct CALModel2D  * clover_model,  int i,  int j){
    int cx = FTNREF1D(i,field.x_min-2);
    int cy = FTNREF1D(j, field.y_min-2);
    if(!checkBoundaries(cx,cy,field.x_min,field.x_max,field.y_min,field.y_max))
        return;

    //default is no predict
    double mult=0.25;
    auto xvel0_sbs= clover_sbs.xvel0;
    auto yvel0_sbs= clover_sbs.yvel0;
    auto xvel1_sbs= clover_sbs.xvel1;
    auto yvel1_sbs= clover_sbs.yvel1;

    if(predict){
        mult*=0.5 ;
        xvel1_sbs = clover_sbs.xvel0;
        yvel1_sbs = clover_sbs.yvel0;
    }


    double recip_volume, energy_change;//, min_cell_volume;
    double right_flux, left_flux, top_flux, bottom_flux, total_flux, volume_change_s;



    //          left_flux=  (xarea(j  ,k  )*(xvel0(j  ,k  )+xvel0(j  ,k+1)                     &
    //+xvel0(j  ,k  )+xvel0(j  ,k+1)))*0.25_8*dt*0.5
    left_flux= (
                get(clover_sbs.xarea,cx,cy)* (
                    get(xvel0_sbs,cx,cy) + get(xvel0_sbs,cx,cy+1) + get(xvel1_sbs,cx,cy) + get(xvel1_sbs,cx,cy+1))
                )*dt*mult;

    //right_flux= (xarea(j+1,k  )*(xvel0(j+1,k  )+xvel0(j+1,k+1)                     &
    //   +xvel0(j+1,k  )+xvel0(j+1,k+1)))*0.25_8*dt*0.5
    right_flux= (
                get(clover_sbs.xarea,cx+1,cy)* (
                    get(xvel0_sbs,cx+1,cy) + get(xvel0_sbs,cx+1,cy+1) + get(xvel1_sbs,cx+1,cy) + get(xvel1_sbs,cx+1,cy+1))
                )*dt*mult;

    //   bottom_flux=(yarea(j  ,k  )*(yvel0(j  ,k  )+yvel0(j+1,k  )                     &
    //+yvel0(j  ,k  )+yvel0(j+1,k  )))*0.25_8*dt*0.5
    bottom_flux = (
                get(clover_sbs.yarea,cx,cy)* (
                    get(yvel0_sbs,cx,cy) + get(yvel0_sbs,cx+1,cy) + get(yvel1_sbs,cx,cy) + get(yvel1_sbs,cx+1,cy))
                )*dt*mult;
    //  top_flux=   (yarea(j  ,k+1)*(yvel0(j  ,k+1)+yvel0(j+1,k+1)                     &
    //+yvel0(j  ,k+1)+yvel0(j+1,k+1)))*0.25_8*dt*0.5
    bottom_flux = (
                get(clover_sbs.yarea,cx,cy+1)* (
                    get(yvel0_sbs,cx,cy+1) + get(yvel0_sbs,cx+1,cy+1) + get(yvel1_sbs,cx,cy+1) + get(yvel1_sbs,cx+1,cy+1))
                )*dt*mult;

    total_flux=right_flux-left_flux+top_flux-bottom_flux;


    const auto _vol = get(clover_sbs.volume,cx,cy);

    volume_change_s = _vol/(_vol+total_flux);



    //     min_cell_volume=MIN(volume(j,k)+right_flux-left_flux+top_flux-bottom_flux &
    //       ,volume(j,k)+right_flux-left_flux                      &
    //       ,volume(j,k)+top_flux-bottom_flux)

    recip_volume = 1.0/_vol;

    const auto _press = get(clover_sbs.pressure,cx,cy);
    const auto _den = get(clover_sbs.density0,cx,cy);
    const auto _visc = get(clover_sbs.viscosity,cx,cy);

    energy_change = (_press/_den + _visc/_den)*total_flux*recip_volume;

    const auto _enrg = get(clover_sbs.energy0,cx,cy);
    set(clover_sbs.energy1,cx,cy,_enrg-energy_change);


    //if(!predict)printf("\t %d \t %d \t %.16e \t %.16e \n",cx,cy,_den*volume_change_s,_enrg-energy_change);

    set(clover_sbs.density1,cx,cy,_den*volume_change_s);
}

void revert_kernel(struct CALModel2D  * clover_model,  int i,  int j){
    int cx = FTNREF1D(i,field.x_min-2);
    int cy = FTNREF1D(j, field.y_min-2);
    if(!checkBoundaries(cx,cy,field.x_min,field.x_max,field.y_min,field.y_max))
        return;


   // printf("\t %d \t %d \t %.16e \t %.16e \n",cx,cy, get(clover_sbs.density0,cx,cy),get(clover_sbs.energy0,cx,cy));
    setCurr(clover_sbs.density1,cx,cy, get(clover_sbs.density0,cx,cy));
    setCurr(clover_sbs.energy1,cx,cy, get(clover_sbs.energy0,cx,cy));

}

void PdV(bool _predict){

    predict = _predict;
    calApplyElementaryProcess2D(clover_model, PdV_kernel);

    //aggiorna i sottostati ceh hai upppppppddddattato
   // calUpdateSubstate2Dr(clover_model, clover_sbs.density1);
 //   calUpdateSubstate2Dr(clover_model, clover_sbs.energy1);
    calUpdate2D(clover_model);

    if(predict){
        ideal_gas(true); //predict = true

        setFields(0);
        fields[FIELD::PRESSURE] = 1;
        FIELD_DEPTH=1;
        update_halo();

    }

    if(predict)
         calApplyElementaryProcess2D(clover_model,revert_kernel);

    calUpdate2D(clover_model);





}


#endif //_PDV_H_
