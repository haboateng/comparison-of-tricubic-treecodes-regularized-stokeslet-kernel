#include <cmath>
#include "kernel_utils.h"

//**********************************************************//

double MxMatrix(double x,double y,double z)
{
  switch(fid){
     case 0:
       return -4.0*x;
     case 1:
       return -y;
     case 2:
       return -z;
     case 3:
       return -2.0*x;
     case 4:
       return 0.0;
     case 5:
       return -2.0*x;
     default:
       return  0.0;
   }
}

double MyMatrix(double x,double y,double z)
{
  switch(fid){
     case 0:
       return -2.0*y;
     case 1:
       return -x;
     case 2:
       return  0.0;
     case 3:
       return -4.0*y;
     case 4:
       return -z;
     case 5:
       return -2.0*y;
     default:
       return  0.0;
   }
}

double MzMatrix(double x,double y,double z)
{
  switch(fid){
     case 0:
       return -2.0*z;
     case 1:
       return  0.0;
     case 2:
       return -x;
     case 3:
       return -2.0*z;
     case 4:
       return -y;
     case 5:
       return -4.0*z;
     default:
       return  0.0;
   }
}


//**********************************************************//

double phiEval(double gr,double x,double y,double z)
{
  switch(fid){
     case 0:
       return gr + x*x;
     case 1:
       return x*y;
     case 2:
       return x*z;
     case 3:
       return gr + y*y;
     case 4:
       return y*z;
     case 5:
       return gr + z*z;
     default:
       return  0.0;
   }
}

//**********************************************************//

double gEval(double r2,double x,double y,double z)
{
  double gvalue;
  switch(fid){
     case 0:
       return r2 + epsq + x*x;
     case 1:
       return x*y;
     case 2:
       return x*z;
     case 3:
       return r2 + epsq + y*y;
     case 4:
       return y*z;
     case 5:
       return r2 + epsq + z*z;
     default:
       return  0.0;
   }
}
//*****************************************************************************//

double fEval_RS(double r2, double x, double y, double z) 
{
  double repsq = r2 + epsq;
  return r8pi*gEval(repsq,x,y,z)/(repsq*sqrt(repsq));
}
//*****************************************************************************//

//*****************************************************************************//
vec_3d f1derivOrd2_3d(double r2h, double x, double y, double z)
{ 
  vec_3d grad;
  double h = hgrid;
  double tth;
   
  tth = x*twoh;
  grad.val[0] = rtwoh*(fEval_RS(r2h - tth,x-h,y,z) - fEval_RS(r2h + tth,x+h,y,z));
    
  tth = y*twoh;
  grad.val[1] = rtwoh*(fEval_RS(r2h - tth,x,y-h,z) - fEval_RS(r2h + tth,x,y+h,z));
    
  tth = z*twoh;
  grad.val[2] = rtwoh*(fEval_RS(r2h - tth,x,y,z-h) - fEval_RS(r2h + tth,x,y,z+h));
 
  return grad;
// the input is  (x_b - x) and we want the derivative wrt to x, so it's the opposite:
// i.e. (f(x-h) - f(x+h))/2h

}
//*****************************************************************************//

double fEval(double rrr) //Use 3-point interpolation to compute function value
{
  int ll     = static_cast<int>(rrr*rdr);
  int l1     = ll+1;
  int l2     = ll+2;
  double ppp = rrr*rdr-static_cast<double>(ll); 

  double vk0 = fval[ll];
  double vk1 = fval[l1];
  double vk2 = fval[l2];
  double t1  = vk0+(vk1-vk0)*ppp;
  double t2  = vk1+(vk2-vk1)*(ppp-1.0);

  double fvalue = (t1+(t2-t1)*ppp*0.5);
  return (t1+(t2-t1)*ppp*0.5);

}
//*****************************************************************************//
double f1derivOrd2(double r2h, double t)
{
  double tth = t*twoh;
  return rtwoh*(fEval(sqrt(r2h - tth)) - fEval(sqrt(r2h + tth))); // the input is
// (x_b - x) and we want the derivative wrt to x, so it's the opposite:
// i.e. (f(x-h) - f(x+h))/2h
}
//*****************************************************************************//
double f1derivOrd4(double r2h, double r4h, double t)
{
  double tth = t*twoh; double fth = 2.0*tth;
  return rtwoh*(  fEval(sqrt(r4h + fth)) 
                + 8.0*(fEval(sqrt(r2h - tth)) - fEval(sqrt(r2h + tth)))
                - fEval(sqrt(r4h - fth)) )/6.0 ; // the input is
// (x_b - x) and we want the derivative wrt to x, so it's the opposite:
// i.e. (f(x+2h) - 8f(x+h)+8f(x-h)-f(x-2h))/2h
}
//*****************************************************************************//

void kerEval_B2_2ndOrder(double x, double y, double z, double b[])
{
  double r2 = x*x + y*y + z*z; 
  double r2h = r2 + hsq;
  double r4h = r2 + 4.0*hsq;

  b[0] = fEval(sqrt(r2));
 
  b[1] = f1derivOrd2(r2h, x);
  b[2] = f1derivOrd2(r2h, y);
  b[3] = f1derivOrd2(r2h, z);

  // b[1] = f1derivOrd4(r2h, r4h, x); b[2] = f1derivOrd4(r2h, r4h, y); b[3] = f1derivOrd4(r2h, r4h, z);
}
//*****************************************************************************//

void kerEval_B2_2ndOrder_RS(double x, double y, double z, double b[])
{
  double r2 = x*x + y*y + z*z;
  double r2h = r2 + hsq;
  double r4h = r2 + 4.0*hsq;

  b[0] = fEval_RS(r2,x,y,z);

  vec_3d grad = f1derivOrd2_3d(r2h, x, y, z);
  b[1] = grad.val[0]; b[2] = grad.val[1]; b[3] = grad.val[2];
}

void kerEval_B2_RS(double x, double y, double z, double b[])
{
  double r2     = x*x + y*y + z*z; double repsq = r2 + epsq;
  double h2r    = r8pi/(repsq*sqrt(repsq));
  double gr     = repsq + epsq; double qr = 3.0/repsq;
  double pval   = phiEval(gr,x,y,z);
  double qrpval = qr*pval;

  b[0] = pval*h2r;
  b[1] = (MxMatrix(x,y,z) + x*qrpval)*h2r; 
  b[2] = (MyMatrix(x,y,z) + y*qrpval)*h2r; 
  b[3] = (MzMatrix(x,y,z) + z*qrpval)*h2r;
}

void kerEval_LM_RS(double x, double y, double z, double b[])
{
  double r2     = x*x + y*y + z*z; double repsq = r2 + epsq;
  double h2r    = r8pi/(repsq*sqrt(repsq));
  double gr     = repsq + epsq; double qr = 3.0/repsq;
  double pval   = phiEval(gr,x,y,z);
  double qrpval = qr*pval;

  double Mx = MxMatrix(x,y,z); double My = MyMatrix(x,y,z); double Mz = MzMatrix(x,y,z);

  b[0] = pval*h2r;
  b[1] = (Mx + x*qrpval)*h2r;
  b[2] = (My + y*qrpval)*h2r;
  b[3] = (Mz + z*qrpval)*h2r;

  double Ixy    = 0.0; double Ixz = 0.0; double Iyz = 0.0;
  if (fid==1) Ixy = h2r; if (fid==2) Ixz = h2r; if (fid==4) Iyz = h2r;                
  double qrh2r = qr*h2r; double xy = x*y; double xz = x*z; double yz = y*z;

  Mx = Mx * qrh2r; My = My * qrh2r; Mz = Mz * qrh2r;
  double ftqr2pvh2r = 5.0/3.0*qrpval*qrh2r;
 
  b[4] = Ixy + x*My + y*Mx + xy*ftqr2pvh2r;
  b[5] = Ixz + x*Mz + z*Mx + xz*ftqr2pvh2r;
  b[6] = Iyz + y*Mz + z*My + yz*ftqr2pvh2r;

  double Ixyz = x*Iyz + y*Ixz + z*Ixy;
  double Mxyz = 5.0/3.0*(yz*Mx + xz*My + xy*Mz);
  b[7] = qr*(Ixyz + Mxyz + 7.0/3.0*x*yz*ftqr2pvh2r);
}

