#ifndef _UPDATE_HALO_H_
#define _UPDATE_HALO_H_

#include<iostream>
#include<model.h>
#include<array>




void update_halo_kernel_left_right(struct CALModel2D  * clover_model,  int i,  int j){
    int cx = FTNREF1D(i,field.x_min-2);
    int cy = FTNREF1D(j, field.y_min-2);

    if(!checkBoundaries(cx,cy,field.x_min-2,field.x_min,field.y_min-2,field.y_max+2))
        return;
    for(int f=0;f<FIELD::NUM_FIELDS ; f++){
        if(fields[f]){
            auto target_field = FIELDS_SUBSTATE[f];
            for(int i=0 ; i < FIELD_DEPTH ; i++){
                setCurr(target_field,i, cy , get(target_field, FIELD_DEPTH+1-i, cy ));
                setCurr(target_field, i+field.y_max , cy, get(target_field, field.y_max-FIELD_DEPTH-i+1, cy));
            }
        }
    }
}


void update_halo_kernel_top_bottom(struct CALModel2D  * clover_model,  int i,  int j){
    int cx = FTNREF1D(i,field.x_min-2);
    int cy = FTNREF1D(j, field.y_min-2);

    if(!checkBoundaries(cx,cy,field.x_min-2,field.x_max+2,field.y_min-2,field.y_min))
        return;

    for(int f=0;f<FIELD::NUM_FIELDS ; f++){
        if(fields[f]){
            auto target_field = FIELDS_SUBSTATE[f];
            for(int i=0 ; i < FIELD_DEPTH ; i++){
                setCurr(target_field,cx, i , get(target_field, cx , FIELD_DEPTH+1-i));
                setCurr(target_field,cx, i+field.y_max , get(target_field, cx , field.y_max-FIELD_DEPTH-i+1));
            }
        }
    }
}

void update_halo(){


    calApplyElementaryProcess2D(clover_model, update_halo_kernel_top_bottom);
    calApplyElementaryProcess2D(clover_model, update_halo_kernel_left_right);


}


#endif //_UPDATE_HALO_H_
