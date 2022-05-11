#include <vector>
#include <string>

#include "kernel_utils.h"

using namespace std;

extern size_t node_count;

extern int N0; // leaf size
extern double mtheta;
extern double sq_theta; // theta^2
extern int max_level;
extern vec_4d *cpvelo; // cp treecode approximation of the potential

extern int MACflag;
extern double maxCradius;

//**************//

struct panel
{
  size_t members[2];
  double xinterval[2];
  double yinterval[2];
  double zinterval[2];
  double xc; // panel center x coordinate
  double yc; // panel center y coordinate
  double zc; // panel center z coordinate
  double rxl; // inverse of panel length parallel to x axis
  double ryl; // inverse of panel length parallel to y axis
  double rzl; // inverse of panel length parallel to z axis
  vector<size_t> children;
  double radius;
  double MAC; // r^2 / theta^2
  double moments1[Pflat];
  double dxmoments1[Pflat];
  double dymoments1[Pflat];
  double dzmoments1[Pflat];
  double moments2[Pflat];
  double dxmoments2[Pflat];
  double dymoments2[Pflat];
  double dzmoments2[Pflat];
  double moments3[Pflat];
  double dxmoments3[Pflat];
  double dymoments3[Pflat];
  double dzmoments3[Pflat];


  int moment_flag;
  panel() // initialization
  {
    moment_flag = 0;
    members[0] = 0;
    members[1] = -1;
    for (size_t kk = 0; kk < Pflat + 1; kk++){
      moments1[kk] = 0.0;
      dxmoments1[kk] = 0.0;
      dymoments1[kk] = 0.0;
      dzmoments1[kk] = 0.0;
      moments2[kk] = 0.0;
      dxmoments2[kk] = 0.0;
      dymoments2[kk] = 0.0;
      dzmoments2[kk] = 0.0;
      moments3[kk] = 0.0;
      dxmoments3[kk] = 0.0;
      dymoments3[kk] = 0.0;
      dzmoments3[kk] = 0.0;

   }
  }
};

extern vector<panel> tree;
extern vector<size_t> leaf;
//********************************************************************//

void read_tree_params(string &treeMethod,
		      int &N_cube,
		      int &N0,
		      double &theta,
                      double &epsilon);

void read_direct_sum(const string& kernelName,
		     int N_cube,
		     vec_3d v_true[],
		     long &ds_cpu_time);

void build_tree_init(int newMAC, struct xyz &particles);

void Swap(size_t i, size_t j, double *lambda, struct xyz &s);

void split_tree_node(int newMAC,
		     size_t panel_index,
		     double *lambda[3],
		     struct xyz &particles);

void build_tree_3D_Recursive(int newMAC,
			     size_t panel_index,
			     double *lambda[3],
			     struct xyz &particles,
			     int level);

vec_3d Call_Ds_PC(int limit_1, int limit_2, 
	       double p_x, double p_y, double p_z,
	       struct xyz &particles,
	       double *lambda[3]);

void Call_Ds_CP(int limit_1, int limit_2, 
               double p_x, double p_y, double p_z,
               struct xyz &particles);

void write_treecode_sum(const string& kernelName,
			const string &type,
			int N_cube,
			long ds_cpu_time,
			long treecode_cpu_time,
			const vec_3d v_true[],
			const vec_4d velo[]);

void compute_error(const vec_3d v_true[],
		   const vec_4d velo[],
		   double &e_n,
		   double &e_d,
		   double &e_d_ex,
		   double &E);
