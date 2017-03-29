#ifndef GENERATE_H_
#define GENERATE_H_

#include<model.h>
#include<utils.h>

void generate_chunk_kernel(struct CALModel2D  * clover_model,  int i,  int j){

    int cx = FTNREF1D(i,field.x_min-2);
    int cy = FTNREF1D(j, field.y_min-2);
    if(!checkBoundaries(cx,cy,field.x_min-2,field.x_max+2,field.y_min-2,field.y_max+2))
        return;



    double radius, x_cent, y_cent;
    int is_in = 0;
    int is_in2 = 0;

    //State 0(aka 1) is always the background state
    //     energy0[cx]= states[0].energy;
    calInit2Dr(clover_model, clover_sbs.energy0, cx , cy , states[0].energy);
    //    density0[OPS_ACC3(0,0)]= states[0].density;
    calInit2Dr(clover_model, clover_sbs.density0, cx , cy , states[0].density);
    //    xvel0[OPS_ACC4(0,0)]=states[0].xvel;
    calInit2Dr(clover_model, clover_sbs.xvel0, cx , cy , states[0].xvel);
    //    yvel0[OPS_ACC5(0,0)]=states[0].yvel;
    calInit2Dr(clover_model, clover_sbs.yvel0, cx , cy , states[0].yvel);



    for(int i = 1; i<number_of_states; i++) {
        x_cent=states[i].xmin;
        y_cent=states[i].ymin;
        is_in = 0;
        is_in2 = 0;

        if (states[i].geometry == g_rect) {
            if(calGetMatrixElement(clover_sbs.vertexx,1 , cx+1, 0 ) >= states[i].xmin &&
                    calGetMatrixElement(clover_sbs.vertexx,1 , cx, 0 )< states[i].xmax )


                if(calGetMatrixElement(clover_sbs.vertexy,1 , cy+1, 0 ) >= states[i].ymin &&
                        calGetMatrixElement(clover_sbs.vertexy,1 , cy, 0 ) < states[i].ymax)
                {
                    calInit2Dr(clover_model, clover_sbs.energy0, cx , cy , states[i].energy);
                    calInit2Dr(clover_model, clover_sbs.density0, cx , cy , states[i].density);

                    for(int cxt=cx ; cxt<=cx+1 ;cxt++){
                        for(int cyt = cy ; cyt <= cy+1 ; cyt++){
                            calInit2Dr(clover_model, clover_sbs.xvel0, cxt, cyt , states[i].xvel);
                            calInit2Dr(clover_model, clover_sbs.yvel0, cxt , cyt , states[i].yvel);

                        }
                    }

                }

        }//geometry rect

        /**TO-DO
         * Missing g_circ and g_point geometry initializer here!
         * */

    }
}

void generate() {

    calApplyElementaryProcess2D(clover_model, generate_chunk_kernel);

}
#endif //GENERATE_H_
