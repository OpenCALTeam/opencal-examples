#ifndef VISCOSITY_H_
#define VISCOSITY_H_
#include<model.h>
#include<cmath>
void viscosity_kernel(struct CALModel2D  * clover_model,  int i,  int j){

    int cx = FTNREF1D(i,field.x_min-2);
    int cy = FTNREF1D(j, field.y_min-2);
    if(!checkBoundaries(cx,cy,field.x_min,field.x_max,field.y_min,field.y_max))
        return;

    double ugrad, vgrad,
            grad2,
            pgradx,pgrady,
            pgradx2,pgrady2,
            grad,
            ygrad, xgrad,
            div,
            strain2,
            limiter,
            pgrad,
            dirx,diry;

    //ugrad=(xvel0(j+1,k  )+xvel0(j+1,k+1))-(xvel0(j  ,k  )+xvel0(j  ,k+1))
    ugrad = (get(clover_sbs.xvel0,cx+1,cy) + get(clover_sbs.xvel0,cx+1,cy+1))
            -
            (get(clover_sbs.xvel0,cx,cy)+get(clover_sbs.xvel0,cx,cy+1));

    //vgrad=(yvel0(j  ,k+1)+yvel0(j+1,k+1))-(yvel0(j  ,k  )+yvel0(j+1,k  ))
    vgrad = (get(clover_sbs.yvel0,cx,cy+1) + get(clover_sbs.yvel0,cx+1,cy+1))
            -
            (get(clover_sbs.yvel0,cx,cy) + get(clover_sbs.yvel0,cx+1,cy) );

    //div = (celldx(j)*(ugrad)+  celldy(k)*(vgrad))
    //calGetMatrixElement(clover_sbs.vertexy,1,cx,0)
    div =   calGetMatrixElement(clover_sbs.celldx,1,cx,0)*ugrad
            +
            calGetMatrixElement(clover_sbs.celldy,1,cy,0)*vgrad;

    //strain2 = 0.5_8*(xvel0(j,  k+1) + xvel0(j+1,k+1)-xvel0(j  ,k  )-xvel0(j+1,k  ))/celldy(k) &
    //         + 0.5_8*(yvel0(j+1,k  ) + yvel0(j+1,k+1)-yvel0(j  ,k  )-yvel0(j  ,k+1))/celldx(j)

    strain2 = 0.5* ( get(clover_sbs.xvel0,cx,cy+1) + get(clover_sbs.xvel0,cx+1,cy+1) - get(clover_sbs.xvel0,cx,cy) -get(clover_sbs.xvel0,cx+1,cy)  )/ calGetMatrixElement(clover_sbs.celldy,1,cy,0)
            +
            0.5* ( get(clover_sbs.yvel0,cx+1,cy) + get(clover_sbs.yvel0,cx+1,cy+1) - get(clover_sbs.yvel0,cx,cy) -get(clover_sbs.yvel0,cx+1,cy)  )/ calGetMatrixElement(clover_sbs.celldx,1,cx,0);

    // pgradx=(pressure(j+1,k)-pressure(j-1,k))/(celldx(j)+celldx(j+1))
    pgradx = (get(clover_sbs.pressure,cx+1,cy) - get(clover_sbs.pressure,cx-1,cy))/(calGetMatrixElement(clover_sbs.celldx,1,cx,0)+calGetMatrixElement(clover_sbs.celldx,1,cx+1,0));
    //pgrady=(pressure(j,k+1)-pressure(j,k-1))/(celldy(k)+celldy(k+1))
    pgrady =  (get(clover_sbs.pressure,cx,cy+1) - get(clover_sbs.pressure,cx,cy-1))/(calGetMatrixElement(clover_sbs.celldy,1,cy,0)+calGetMatrixElement(clover_sbs.celldy,1,cy+1,0));

    pgradx2 = pgradx*pgradx;
    pgrady2 = pgrady*pgrady;

    //  limiter = ((0.5_8*(ugrad)/celldx(j))*pgradx2+(0.5_8*(vgrad)/celldy(k))*pgrady2+strain2*pgradx*pgrady)  &
    //                /MAX(pgradx2+pgrady2,1.0e-16_8)

    limiter =   (
                (0.5*ugrad/calGetMatrixElement(clover_sbs.celldx,1,cx,0)) * pgradx2 +
                (0.5*vgrad/calGetMatrixElement(clover_sbs.celldy,1,cy,0))*pgrady2   +
                strain2*pgradx*pgrady
                )/std::max(pgradx2 + pgrady2 , 1.0e-16);

    if(limiter > 0.0 || div >= 0.0)
        set(clover_sbs.viscosity,cx,cy,0.0);

    else{
        dirx=1.0;
        if(pgradx < 0.0) dirx=-1.0;

        pgradx = dirx*std::max(1.0e-16,std::fabs(pgradx));
        diry=1.0;
        if(pgradx < 0.0)
            diry=-1.0;
        pgrady = diry*std::max(1.0e-16,std::fabs(pgrady));
        pgrad = std::sqrt(pgradx*pgradx+pgrady*pgrady);

        xgrad = std::fabs(calGetMatrixElement(clover_sbs.celldx,1,cx,0)*pgrad/pgradx);

        ygrad = std::fabs(calGetMatrixElement(clover_sbs.celldy,1,cy,0)*pgrad/pgrady);
        grad  = std::min(xgrad,ygrad);
        grad2 = grad*grad;

        //viscosity[OPS_ACC6(0,0)] = 2.0 * (density0[OPS_ACC5(0,0)]) * grad2 * limiter * limiter;
        const auto den = get(clover_sbs.density0,cx,cy);
        set(clover_sbs.viscosity,cx,cy,
            2.0 * ( den ) * grad2 * limiter * limiter
                );

    }

}

void viscosity(){
    calApplyElementaryProcess2D(clover_model, viscosity_kernel);
    calUpdateSubstate2Dr(clover_model,clover_sbs.viscosity);
}



#endif //VISCOSITY_H_
