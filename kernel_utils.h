
#include "utilities.h"

//*******************************************************************//

double fEval(double rrr);
double f1derivOrd2(double r2h, double t);
double f1derivOrd4(double r2h, double r4h, double t);
void kerEval_B2_2ndOrder(double x, double y, double z, double b[]);

double MxMatrix(double x, double y, double z);
double MyMatrix(double x, double y, double z);
double MzMatrix(double x, double y, double z);
double phiEval(double gr,double x,double y,double z);

double gEval(double r2,double x,double y,double z);
double fEval_RS(double r2, double x, double y, double z);
vec_3d f1derivOrd2_3d(double r2h, double x, double y, double z);
void kerEval_B2_2ndOrder_RS(double x, double y, double z, double b[]);
void kerEval_B2_RS(double x, double y, double z, double b[]);
void kerEval_LM_RS(double x, double y, double z, double b[]);
