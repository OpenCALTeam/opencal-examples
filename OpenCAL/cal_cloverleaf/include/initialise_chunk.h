#ifndef INITIALISE_CHUNK_H_
#define INITIALISE_CHUNK_H_

#include<model.h>
#include<utils.h>

//int rangex[] = {x_min-2, x_max+3, y_min-2, y_max+3};
//int rangey[] = {x_min-2, x_max+3, y_min-2, y_max+3};
//int rangefull[] = {-2, x_cells+8, -2, y_cells+8};




void initialise_chunk_kernel_x(struct CALModel2D  * clover_model,  int i,  int j){
    int cx = FTNREF1D(i,field.x_min-2);
    int cy = FTNREF1D(j, field.y_min-2);
    if(!checkBoundaries(cx,cy,field.x_min-2,field.x_max+3,0,1))
        return;

    int fx = i - field.x_min; //2
    const int x_min=field.x_min-2;
    const double d_x = (grid.xmax - grid.xmin)/(double)grid.x_cells;
    const double min_x=grid.xmin+d_x*field.left;

    //vertexx[FTNREF1D(j,x_min-2)]=min_x+d_x*(double)(j-x_min);
    //    vertexx[OPS_ACC0(0,0)] = min_x + d_x * (xx[OPS_ACC1(0,0)] - x_min);

    calSetMatrixElement(clover_sbs.vertexx,1,cx,0,min_x+d_x*(double)(fx-x_min));
    //      vertexdy[OPS_ACC2(0,0)] = (double)d_y;

    //vertexdy[FTNREF1D(k,y_min-2)]=d_y;
    calSetMatrixElement(clover_sbs.vertexdx,1,cx,0,d_x);


}

/**
 * @brief works on 1d buffers so J = 0 has to be enforced
 * @param clover_model
 * @param i
 * @param j always zero
 */
void initialise_chunk_kernel_y(struct CALModel2D  * clover_model,  int i,  int j){

    int cx = FTNREF1D(i,field.x_min-2);
    int cy = FTNREF1D(j, field.y_min-2);



    if( !checkBoundaries(cx,cy,field.x_min-2,field.x_max+3,0,1))
        return;

    int fy = i - field.y_min; //-2
    const int y_min=field.y_min-2;
    const double d_y = (grid.ymax - grid.ymin)/(double)grid.y_cells;
    const double min_y=grid.ymin+d_y*field.bottom;

    //    vertexy[FTNREF1D(k,y_min-2)]=min_y+d_y*(double)(k-y_min);
    //    vertexy[OPS_ACC0(0,0)] = min_y + d_y * (yy[OPS_ACC1(0,0)] - y_min);

    calSetMatrixElement(clover_sbs.vertexy,1,cx,0,min_y+d_y*(double)(fy-y_min));

    //vertexdy[OPS_ACC2(0,0)] = (double)d_y;
    //vertexdy[FTNREF1D(k,y_min-2)]=d_y;
    calSetMatrixElement(clover_sbs.vertexdy,1,cx,0,d_y);


}



void initialise_chunk_kernel_celly(struct CALModel2D  * clover_model,  int i,  int j){
    int cx = FTNREF1D(i,field.x_min-2);
    int cy = FTNREF1D(j, field.y_min-2);
    if(!checkBoundaries(cx,cy,field.x_min-2,field.x_max+2,0,1))
        return;

    const double d_y = (grid.ymax - grid.ymin)/(double)grid.y_cells;

    //celly[FTNREF1D(k,y_min-2)]=0.5*(vertexy[FTNREF1D(k,y_min-2)]+vertexy[FTNREF1D(k+1,x_min-2)]);
    //celly[OPS_ACC1(0,0)] = 0.5*( vertexy[OPS_ACC0(0,0)]+ vertexy[OPS_ACC0(0,1)] );

    calSetMatrixElement(clover_sbs.celly,1,cx,0,
                        0.5*(calGetMatrixElement(clover_sbs.vertexy,1,cx,0)  + calGetMatrixElement(clover_sbs.vertexy,1,cx+1,0))
                        );
    // celldy[OPS_ACC2(0,0)] = d_y;
    calSetMatrixElement(clover_sbs.celldy,1,cx,0,d_y);


}

void initialise_chunk_kernel_volume(struct CALModel2D  * clover_model,  int i,  int j){

    int cx = FTNREF1D(i,field.x_min-2);
    int cy = FTNREF1D(j, field.y_min-2);
    if(!checkBoundaries(cx,cy,field.x_min-2,field.x_max+2,field.y_min-2,field.y_max+2))
        return;


    const double d_x = (grid.xmax - grid.xmin)/(double)grid.x_cells;
    const double d_y = (grid.ymax - grid.ymin)/(double)grid.y_cells;

    //ricontrolla dimensioni
    calInit2Dr(clover_model, clover_sbs.volume, cx , cy , d_x*d_y);

    //  xarea[FTNREF2D(j,k,x_max+5,x_min-2,y_min-2)]=celldy[FTNREF1D(k,y_min-2)];
    calInit2Dr(clover_model, clover_sbs.xarea, cx , cy ,
               calGetMatrixElement(clover_sbs.celldy , 1, cx , 0 )
               );

    //yarea[OPS_ACC4(0,0)] = celldx[OPS_ACC3(0,0)];
    calInit2Dr(clover_model, clover_sbs.yarea, cx , cy ,
               calGetMatrixElement(clover_sbs.celldx , 1, cx , 0 )
               );


}


void initialise_chunk_kernel_cellx(struct CALModel2D  * clover_model,  int i,  int j){
    int cx = FTNREF1D(i,field.x_min-2);
    int cy = FTNREF1D(j, field.y_min-2);
    if(!checkBoundaries(cx,cy,field.x_min-2,field.x_max+2,0,1))
        return;

    const double d_x = (grid.xmax - grid.xmin)/(double)grid.x_cells;

    //cellx[OPS_ACC1(0,0)]  = 0.5*( vertexx[OPS_ACC0(0,0)] + vertexx[OPS_ACC0(1,0)] );
    //cellx[FTNREF1D(j,x_min-2)]=0.5*(vertexx[FTNREF1D(j,x_min-2)]+vertexx[FTNREF1D(j+1,x_min-2)]);

    calSetMatrixElement(clover_sbs.cellx,1,cx,0,
                        0.5 * (calGetMatrixElement(clover_sbs.vertexx,1,cx,0)  +calGetMatrixElement(clover_sbs.vertexx,1,cx+1,0))
                        );

    //celldx[OPS_ACC2(0,0)]  = d_x;
    //celldx(j)=d_x
    calSetMatrixElement(clover_sbs.celldx,1,cx,0,d_x);

}


void initialise_chunk(){

    calApplyElementaryProcess2D(clover_model, initialise_chunk_kernel_x);
    calApplyElementaryProcess2D(clover_model, initialise_chunk_kernel_y);
    calApplyElementaryProcess2D(clover_model, initialise_chunk_kernel_cellx);
    calApplyElementaryProcess2D(clover_model, initialise_chunk_kernel_celly);
    calApplyElementaryProcess2D(clover_model, initialise_chunk_kernel_volume);

}

#endif //INITIALISE_CHUNK
