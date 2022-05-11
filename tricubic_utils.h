
#include "tree_utils.h"

//*******************************************************************//

void Panel_Moment_Tricubic(size_t panel_index,
			   double *lambda[3],
			   struct xyz &particles,
			   double m1[],
                           double m1x[], double m1y[], double m1z[],
                           double m2[],
                           double m2x[], double m2y[], double m2z[],
                           double m3[],
                           double m3x[], double m3y[], double m3z[]);

void getClusterPoints_B2(size_t panel_index,
			  double &dx, double &dy, double &dz,
			  double ccx[], double ccy[], double ccz[]);

void getClusterPoints_B6(size_t panel_index,
                          double &dx, double &dy, double &dz,
                          double ccx[], double ccy[], double ccz[]);

void getClusterCorners_LM(size_t panel_index,
                          double &dx, double &dy, double &dz,
                          double ccx[], double ccy[], double ccz[]);

void scaleDerivatives_B2(double b[], double dx, double dy, double dz);

void scaleDerivatives_LM(double b[], double dx, double dy, double dz);

void Compute_bvec_B2_fdiff(double bvec[][64],
		     double x, double y, double z,
		     size_t panel_index);

void Compute_bvec_B2(double bvec[][64],
                     double x, double y, double z,
                     size_t panel_index);

void Compute_bvec_B6(double bvec[][64],
                     double x, double y, double z,
                     size_t panel_index);

void Compute_bvec_LM(double bvec[][64],
                     double x, double y, double z,
                     size_t panel_index);

void Compute_Coeff_B2_fdiff(double x, double y, double z,
                     size_t panel_index);

void Compute_Coeff_B2(double x, double y, double z,
                     size_t panel_index);
                    
void Compute_Coeff_B6(double x, double y, double z,
                     size_t panel_index); 

void Compute_Coeff_LM(double x, double y, double z,
                     size_t panel_index);

vec_4d assemble_velocity(const double vec[][64], int panel_index);

vec_4d Call_Tricubic_B2_fdiff(double px, double py, double pz,
			int panel_index);

vec_4d Call_Tricubic_B2(double px, double py, double pz,
                        int panel_index);

vec_4d Call_Tricubic_B6(double px, double py, double pz,
                        int panel_index);

vec_4d Call_Tricubic_LM(double px, double py, double pz,
                        int panel_index);

void BinvMultiply_Transpose_B2(const double b[], double a[]);

void BinvMultiply_B2(const double b[], double a[]);

void BinvMultiply_Transpose_LM(const double b[], double a[]);

void BinvMultiply_LM(const double b[], double a[]);

void BinvMultiply_Transpose_B6(const double b[], double a[]);

void BinvMultiply_B6(const double b[], double a[]);

