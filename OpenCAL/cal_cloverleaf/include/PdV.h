#ifndef _PDV_H_
#define _PDV_H_

#include<iostream>
#include<model.h>
void PdV_kernel_PREDICT(struct CALModel2D  * clover_model,  int i,  int j){
    int cx = FTNREF1D(i,field.x_min-2);
    int cy = FTNREF1D(j, field.y_min-2);
    if(!checkBoundaries(cx,cy,field.x_min,field.x_max,field.y_min,field.y_max))
        return;

    double recip_volume, energy_change;//, min_cell_volume;
    double right_flux, left_flux, top_flux, bottom_flux, total_flux, volume_change_s;



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

}

void PdV(bool predict){

    if(predict)
        calApplyElementaryProcess2D(clover_model, PdV_kernel_PREDICT);
    else{
        std::cerr<<"PDVKERNEL _ PREDICT=FALSE NOT IMPLEMENTED YET\n EXITING...\n";
        exit(-1);
    }
    //aggiorna i sottostati ceh hai upppppppddddattato
    calUpdateSubstate2Dr(clover_model, clover_sbs.density1);
    calUpdateSubstate2Dr(clover_model, clover_sbs.energy1);

}


#endif //_PDV_H_
