// Tricubic LM version 

#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "tricubic_utils.h"

using namespace std;

//*****************************************************************************//

void Compute_CP2(size_t panel_index,
                           struct xyz &particles)
{
 if (tree[panel_index].moment_flag == 1)
 {

   double tm1[Pflat];
   double tm2[Pflat];
   double tm3[Pflat];

   double dm[Pflat];

// Multiply coefficient by B inverse

   BinvMultiply_LM(tree[panel_index].moments1, tm1);
   BinvMultiply_LM(tree[panel_index].moments2, tm2);
   BinvMultiply_LM(tree[panel_index].moments3, tm3);

//============================================================
  double xmin = tree[panel_index].xinterval[0];
  double ymin = tree[panel_index].yinterval[0];
  double zmin = tree[panel_index].zinterval[0];

  double rdx = tree[panel_index].rxl; 
  double rdy = tree[panel_index].ryl; 
  double rdz = tree[panel_index].rzl;

  double tp0 = tree[panel_index].members[0];
  double tp1 = tree[panel_index].members[1];

  for (size_t tp_j = tp0; tp_j <= tp1; tp_j++)
    {
      size_t old_j = particles.old_index[tp_j];

      double x = (particles.x[tp_j] - xmin) * rdx;
      double y = (particles.y[tp_j] - ymin) * rdy;
      double z = (particles.z[tp_j] - zmin) * rdz;

      double peng1 = 0.0; double peng2 = 0.0; double peng3 = 0.0;
      double si   = 1.0;
      for (int i = 0; i < P + 1; i++)
       {
          double sij = si;
          for (int j = 0; j < P + 1; j++)
            {
              int i4j = i + 4*j;

              double s = sij;
              for (int k = 0; k < P + 1; k++)
                {
                  int ii = i4j + 16*k;

                  peng1  += tm1[ii] * s;
                  peng2  += tm2[ii] * s;
                  peng3  += tm3[ii] * s;

                  dm[ii] = s;

                  s *= z;
                }
                sij *= y;
            }
            si *= x;
        }
         cpvelo[old_j].val[0] += peng1;
         cpvelo[old_j].val[1] += peng2;
         cpvelo[old_j].val[2] += peng3;

      double tvx = 0.0; double tvy = 0.0; double tvz = 0.0;
      for (int i = 1; i < P + 1; i++)
        { 

          double di = static_cast<double>(i);

          int fourI = 4*i; int sixtI = 16*i;
          int fII = 4*(i-1); int sII = 16*(i-1);

          for (int j = 0; j < P + 1; j++)
            { 

              int fourJ = 4*j;
              int i4j   = i + fourJ;
              int j4i   = j + fourI;

              int fJsI  = fourJ + sixtI;
              int jfII  = j     + fII;
              int fJsII = fourJ + sII;

              for (int k = 0; k < P + 1; k++)
                {
                  int sk = 16*k;
                  int ii = i4j + sk;
                  int jj = j4i + sk;
                  int kk = k   + fJsI;
 
                  int njj = jfII + sk;
                  int nkk = k    + fJsII;

                  tvx   += tm1[ii] * di * dm[ii-1];
                  tvy   += tm2[jj] * di * dm[njj];
                  tvz   += tm3[kk] * di * dm[nkk];

                }
            }
        }
         cpvelo[old_j].val[3] += rdx*tvx + rdy*tvy + rdz*tvz;

    }
 }

// Loop over children
   size_t length = tree[panel_index].children.size();
   for (size_t i = 0; i < length; i++)
     {
       size_t index = tree[panel_index].children[i];
       Compute_CP2(index, particles);
     }
}

//*****************************************************************************//

void Compute_CP1(struct xyz &particles,
                 double p_x, double p_y, double p_z,
                 size_t panel_index)
{
    size_t limit_1 = tree[panel_index].members[0];
    size_t limit_2 = tree[panel_index].members[1];   
        
    double tpx = p_x - tree[panel_index].xc; 
    double tpy = p_y - tree[panel_index].yc; 
    double tpz = p_z - tree[panel_index].zc;
        
    double R_sq = tpx * tpx + tpy * tpy + tpz * tpz;
    
    if (tree[panel_index].MAC < R_sq)
      {
        tree[panel_index].moment_flag = 1;
        Compute_Coeff_LM(p_x, p_y, p_z, panel_index);

        maxCradius = (tree[panel_index].radius > maxCradius) ? tree[panel_index].radius : maxCradius;

      }
    
    else
      {
        if (limit_2 - limit_1 < N0) //  otherwise, if cluster is a leaf, use direct sum
          {
            Call_Ds_CP(limit_1, limit_2, p_x,  p_y, p_z, particles);
          }
        else // othervise, if cluster is not a leaf, look at children
          {
            size_t length = tree[panel_index].children.size();
            for (size_t i = 0; i < length; i++)
              {
                size_t index = tree[panel_index].children[i];
                Compute_CP1(particles, p_x,  p_y, p_z, index);
              }            
          }
      }
}

//*****************************************************************************//

vec_4d Compute_Velocity(double *lambda[3],
			struct xyz &particles,
			double p_x, double p_y, double p_z,
			size_t panel_index)
{
    size_t limit_1 = tree[panel_index].members[0];
    size_t limit_2 = tree[panel_index].members[1];   
	
    double tpx = p_x - tree[panel_index].xc; 
    double tpy = p_y - tree[panel_index].yc; 
    double tpz = p_z - tree[panel_index].zc;
	
    double R_sq = tpx * tpx + tpy * tpy + tpz * tpz;
    
    vec_4d velocity;
    velocity.val[0] = 0.0; velocity.val[1] = 0.0; velocity.val[2] = 0.0; 
    velocity.val[3] = 0.0;

    if (tree[panel_index].MAC < R_sq)
      {
        vec_4d tree_result = Call_Tricubic_LM(p_x, p_y, p_z,
					      panel_index);
					      
        velocity.val[0] += tree_result.val[0];
        velocity.val[1] += tree_result.val[1];
        velocity.val[2] += tree_result.val[2];
        velocity.val[3] += tree_result.val[3];

        maxCradius = (tree[panel_index].radius > maxCradius) ? tree[panel_index].radius : maxCradius;

      }
    
    else
      {
	if (limit_2 - limit_1 < N0) //if cluster is a leaf, use direct sum
	  {
            vec_3d DS_result = Call_Ds_PC(limit_1, limit_2,
				       p_x,  p_y, p_z,
				       particles,
				       lambda);
	    
            velocity.val[0] += DS_result.val[0] ;
            velocity.val[1] += DS_result.val[1] ;
            velocity.val[2] += DS_result.val[2] ;

            
	  }
	else //if cluster is not a leaf, look at children
	  {
	    velocity.val[0] = 0.0; 
            velocity.val[1] = 0.0;
            velocity.val[2] = 0.0;
            velocity.val[3] = 0.0;

	    size_t length = tree[panel_index].children.size();
	    for (size_t i = 0; i < length; i++)
	      {
		size_t index = tree[panel_index].children[i];
                vec_4d temp_result = Compute_Velocity(lambda,
						      particles,
						      p_x,  p_y, p_z,
						      index);

                velocity.val[0] += temp_result.val[0];
                velocity.val[1] += temp_result.val[1];
                velocity.val[2] += temp_result.val[2];
                velocity.val[3] += temp_result.val[3];
	      }            
	  }
      }
    return velocity;
}

//*****************************************************************************//

int main()
{
    tree.reserve(5000);
    leaf.reserve(5000);

    double theta;
    string kernelName;
    string treeMethod;
    read_tree_params(treeMethod, N_cube, N0, theta, epsilon);

    cout<<" "<<endl;
    cout<<"Using "<<treeMethod<<" treecode"<<endl;
    cout<<" "<<endl;

    epsq = epsilon*epsilon;
    kernelName = "RegStokeslet";

    mtheta = theta;
    sq_theta = theta*theta; // theta^2

    struct xyz particles(N_cube);
    double *lambda[3];
    lambda[0] = new double[N_cube];
    lambda[1] = new double[N_cube];
    lambda[2] = new double[N_cube];

    cout << "===== Tricubic LM version =====" << endl;
    if (fdiff==1) cout << "===== with finite differences =====" << endl;
    cout << "Kernel = Regularized Stokeslet" << endl;
    cout << "P = " << P << endl;
    cout << "N = " << N_cube << endl;
    cout << "theta = " << theta << endl;
    cout << "N0 = " << N0 << endl;

    // read particle coordinates and strengths from files
    read_particle_data(N_cube, particles, lambda);

    
    if (fdiff==1) {
      twoh = 2.0*hgrid; rtwoh = 1.0/twoh; hsq = hgrid*hgrid; 
      rfhsq = 0.25/hsq; r8hcb = rtwoh*rfhsq;
   
      // grid size
      
      double xlmax = particles.x[0];
      double xlmin = particles.x[0];
      double ylmax = particles.y[0];
      double ylmin = particles.y[0];
      double zlmax = particles.z[0];
      double zlmin = particles.z[0];
      
      for (int i = 1; i < N_cube; i++)
	{
	  double xi, yi, zi;
	  xi=particles.x[i]; yi=particles.y[i]; zi=particles.z[i];
	  
	  if (xi > xlmax) xlmax = xi;
	  if (xi < xlmin) xlmin = xi;
	  if (yi > ylmax) ylmax = yi;
	  if (yi < ylmin) ylmin = yi;
	  if (zi > zlmax) zlmax = zi;
	  if (zi < zlmin) zlmin = zi;
	}
      
      double maxlen = max(xlmax-xlmin,ylmax-ylmin);
      maxlen = sqrt(3.0)*max(maxlen,zlmax-zlmin);
      double dr = maxlen/static_cast<double>(mgrid);
      rdr = 1.0/dr;
      
      // Call look-up table to compute function values at grid points
      fval = new double[mgrid]; //Allocate mgrid doubles and save ptr in fval
    }
    
    cout << "Starting treecode" << endl;
	
    //***************** Set up tree *******************************
    long Start_total, Start_btree;
    long End_total, End_btree;
    
    Start_total = getTickCount(); // Get currenct CPU time
    Start_btree = getTickCount();

    int newMAC = MACflag;
    build_tree_init(newMAC, particles);
    build_tree_3D_Recursive(newMAC, 0, lambda, particles, 0);

    if(treeMethod=="Particle-Cluster"){    

    //***************** Compute moment for each panel **************
       size_t size = tree.size();
    
       for (size_t k = 1; k < size; k++) // skip root
         {

           double tm1[Pflat], tm1x[Pflat], tm1y[Pflat], tm1z[Pflat];
           double tm2[Pflat], tm2x[Pflat], tm2y[Pflat], tm2z[Pflat];
           double tm3[Pflat], tm3x[Pflat], tm3y[Pflat], tm3z[Pflat];

         
           Panel_Moment_Tricubic(k, lambda, particles,
                                    tm1, tm1x, tm1y, tm1z,
                                    tm2, tm2x, tm2y, tm2z,
                                    tm3, tm3x, tm3y, tm3z);

           BinvMultiply_Transpose_LM(tm1, tree[k].moments1);
           BinvMultiply_Transpose_LM(tm1x, tree[k].dxmoments1);
           BinvMultiply_Transpose_LM(tm1y, tree[k].dymoments1);
           BinvMultiply_Transpose_LM(tm1z, tree[k].dzmoments1);

           BinvMultiply_Transpose_LM(tm2, tree[k].moments2);
           BinvMultiply_Transpose_LM(tm2x, tree[k].dxmoments2);
           BinvMultiply_Transpose_LM(tm2y, tree[k].dymoments2);
           BinvMultiply_Transpose_LM(tm2z, tree[k].dzmoments2);

           BinvMultiply_Transpose_LM(tm3, tree[k].moments3);
           BinvMultiply_Transpose_LM(tm3x, tree[k].dxmoments3);
           BinvMultiply_Transpose_LM(tm3y, tree[k].dymoments3);
           BinvMultiply_Transpose_LM(tm3z, tree[k].dzmoments3);

         }
    
    }
    End_btree = getTickCount();
    
    cout << "build tree time (incl. moments) = "
	 << End_btree - Start_btree << endl;

    //***************** Run Treecode *****************

    divergence = 0.0;

    vec_4d *velo = new vec_4d[N_cube];
    cpvelo = new vec_4d[N_cube];

    if(treeMethod == "Particle-Cluster"){ 
     
      //********Particle-Cluster Treecode**************   
	
       for (int i = 0; i < N_cube; i++)
         {
	   double p_x = particles.x[i];
	   double p_y = particles.y[i];
	   double p_z = particles.z[i];
	
	   int old_i = particles.old_index[i];
	
	   velo[old_i] = Compute_Velocity(lambda,
	 			         particles,
				         p_x, p_y, p_z,
				         0);
 
           divergence += abs(velo[old_i].val[3]);	
//           divergence += velo[old_i].val[3];
         }
    }else if(treeMethod == "Cluster-Particle"){

      //********Cluster-Particle Treecode**************
    
    for (int i = 0; i < N_cube; i++)
      {
        cpvelo[i].val[0] = 0.0;
        cpvelo[i].val[1] = 0.0;
        cpvelo[i].val[2] = 0.0;
        cpvelo[i].val[3] = 0.0;
      }

    double bfac = 1.0;
    for (int i = 0; i < N_cube; i++)
      {
        double p_x = particles.x[i];
        double p_y = particles.y[i];
        double p_z = particles.z[i];

        int old_i = particles.old_index[i];
        sweight1  = lambda[0][i];
        sfweight1 = sweight1;
        sweight2  = lambda[1][i];
        sfweight2 = sweight2;
        sweight3  = lambda[2][i];
        sfweight3 = sweight3;

        sID       = i;

        Compute_CP1(particles, p_x, p_y, p_z, 0);

      }

      Compute_CP2(0, particles);

      for (int i = 0; i < N_cube; i++)
        {
          divergence += abs(cpvelo[i].val[3]);
//          divergence += cpvelo[i].val[3];

        }
     }else{
      // No tree method specified
      cout << "Please specify a tree method" <<endl;
      exit (EXIT_FAILURE);
    }

    End_total = getTickCount(); // Time for all treecode computing
    long treecode_cpu_time;
    treecode_cpu_time = End_total - Start_total;
    
    cout << "treecode_cpu_time = " << treecode_cpu_time << endl;
    
    //***************** End Run Treecode *****************
    
    // ********* read exact data from a file *********

    vec_3d *v_true = new vec_3d[N_cube];

    long ds_cpu_time;
    read_direct_sum(kernelName,
		    N_cube,
		    v_true,
		    ds_cpu_time);
    cout << "ds time = " << ds_cpu_time << endl;

    // ********* write treecode data to a file *********
    string method = "RS_C1";

   if(treeMethod == "Particle-Cluster"){    
       write_treecode_sum(kernelName,
		       method,
		       N_cube,
		       ds_cpu_time,
		       treecode_cpu_time,
		       v_true,
		       velo);
    }

   if(treeMethod == "Cluster-Particle"){
       write_treecode_sum(kernelName,
                       method,
                       N_cube,
                       ds_cpu_time,
                       treecode_cpu_time,
                       v_true,
                       cpvelo);
    }

    // compute error
    double e_n, e_d, e_d_ex, E, fd,fe, FE;

    if(treeMethod == "Particle-Cluster"){
       compute_error(v_true, velo, e_n, e_d, e_d_ex, E);
    }

    if(treeMethod == "Cluster-Particle"){
       compute_error(v_true, cpvelo, e_n, e_d, e_d_ex, E);
    }

    cout << "L2 velocity Error = " << E << endl;
    cout << "Numerator of L2 error = " << e_n << endl;
    cout << "Denominator of L2 error (uses exact value) = " << e_d_ex << endl;
    cout << " " << endl;
    cout << "Divergence = "<< divergence <<endl;
    cout << " " << endl;
    cout << "Divergence per particle = "<< divergence/N_cube <<endl;
    cout << " " << endl;

// Output to file

    ofstream output_file;
    output_file.open("output_RS_C1.txt", ios::out | ios::ate | ios::app);
    output_file << N_cube <<"  "<< N0 << "  " << theta << "  "<< setprecision(16)
    << E <<"  "<< divergence << "  "<< divergence/N_cube <<"  "<<ds_cpu_time <<"  "<<treecode_cpu_time << endl; 

    output_file.close();

   //**************************************************
    
    delete [] lambda[0];
    delete [] lambda[1];
    delete [] lambda[2];
    delete [] v_true;
   
    delete [] velo;

    delete [] cpvelo;

//    if (treeMethod == "Cluster-Particle") delete *Binv;
    if (fdiff==1) delete [] fval;  

    return 0;
}
