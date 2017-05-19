#ifndef _IDEAL_GAS
#define _IDEAL_GAS

#include<model.h>
#include<math.h>



void ideal_gas_kernel0(struct CALModel2D  * clover_model,  int i,  int j){
    int cx = FTNREF1D(i,field.x_min-2);
    int cy = FTNREF1D(j, field.y_min-2);
    if(!checkBoundaries(cx,cy,field.x_min,field.x_max,field.y_min,field.y_max))
        return;

    double sound_speed_squared, v, pressurebyenergy, pressurebyvolume;
    const auto den     = calGet2Dr(clover_model, clover_sbs.density0, cx, cy);
    const auto enrg    = calGet2Dr(clover_model, clover_sbs.energy0, cx, cy);
    const auto press = (1.4 - 1.0) * den * enrg;

    v = 1.0 / den;

    pressurebyenergy = (1.4 - 1.0) * den;
    pressurebyvolume = -1*den * press;
    sound_speed_squared = v*v*(press * pressurebyenergy-pressurebyvolume);

    const auto sp2 =  sqrt(sound_speed_squared);
    calInit2Dr(clover_model, clover_sbs.soundspeed , cx, cy, sp2);
    calInit2Dr(clover_model, clover_sbs.pressure , cx, cy, press);
}


void ideal_gas_kernel1(struct CALModel2D  * clover_model,  int i,  int j){
    int cx = FTNREF1D(i,field.x_min-2);
    int cy = FTNREF1D(j, field.y_min-2);
    if(!checkBoundaries(cx,cy,field.x_min,field.x_max,field.y_min,field.y_max))
        return;

    double sound_speed_squared, v, pressurebyenergy, pressurebyvolume;
    const auto den     = calGet2Dr(clover_model, clover_sbs.density1, cx, cy);
    const auto enrg    = calGet2Dr(clover_model, clover_sbs.energy1, cx, cy);
    const auto press = (1.4 - 1.0) * den * enrg;

    v = 1.0 / den;

    pressurebyenergy = (1.4 - 1.0) * den;
    pressurebyvolume = -1*den * press;
    sound_speed_squared = v*v*(press * pressurebyenergy-pressurebyvolume);

    const auto sp2 =  sqrt(sound_speed_squared);
    calInit2Dr(clover_model, clover_sbs.soundspeed , cx, cy, sp2);
    calInit2Dr(clover_model, clover_sbs.pressure , cx, cy, press);
}

void ideal_gas(const bool predict){

    if(!predict)
        calApplyElementaryProcess2D(clover_model, ideal_gas_kernel0);
    else
        calApplyElementaryProcess2D(clover_model, ideal_gas_kernel1);

     calUpdate2D(clover_model);

}


#endif // _IDEAL_GAS
