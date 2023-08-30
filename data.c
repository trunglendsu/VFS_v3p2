static char help[] = "Testing programming!";

#include <vector>
#include "petscdmda.h"
#include "petscts.h"
#include "petscpc.h"
#include "petscsnes.h"
#include "petscsys.h"   
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#define NEWMETRIC

#ifdef TECIO
#include "TECIO.h"
#endif
PetscBool subgrid_flag=PETSC_FALSE;
char subgrid[256];
int ti, block_number, Flux_in;
PetscInt binary_input=0;
PetscInt xyz_input=0;
PetscInt tis, tie, tsteps=5;
PetscReal angle;
PetscInt nv_once=0;
PetscInt onlyV=0;
PetscInt vorticity_budget=0;
PetscInt tke_budget=0;
PetscInt mke_budget=0;
PetscInt zm_budget=0;
PetscInt k_average=0;
PetscInt j_average=0;
PetscInt i_average=0;
PetscInt ik_average=0;
PetscInt ikc_average=0;	// conditional spatial averaging in ik directions (channel flow)
PetscInt reynolds=0;	// 1: contravariant reynolds stress

PetscInt i_begin, i_end;
PetscInt j_begin, j_end;
PetscInt k_begin, k_end;

PetscInt pcr=0;
PetscInt avg=0, rans=0, rans_output=0, levelset=0;
PetscInt vc = 1;

PetscInt les=0;
PetscInt i_periodic=0;
PetscInt j_periodic=0;
PetscInt k_periodic=0;
PetscInt kk_periodic=0;
PetscInt averaging_option=0;
PetscInt phase=0;
PetscInt ti_phase=0;
PetscInt pi=-1, pk=-1;
PetscInt shear=0;

char prefix[256];
double gravity_x=0, gravity_y=0, gravity_z=0;


//int l, m, n;
/* Symmetric matrix A -> eigenvectors in columns of V, corresponding eigenvalues in d. */
void eigen_decomposition(double A[3][3], double V[3][3], double d[3]);




typedef struct {
  PetscReal t, f;
} FlowWave;

typedef struct {
  PassiveScalar u, v, w, p;
} PassiveField;

typedef struct {
  PetscScalar u, v, w;
} Field;

typedef struct {
  PetscScalar x, y, z;
} Cmpnts;

typedef struct {
  PetscScalar x, y;
} Cmpnts2;

typedef struct {
  PassiveScalar csi[3], eta[3], zet[3], aj;
} Metrics;

typedef struct {
  Vec	Ubcs; // An ugly hack, waste of memory
} BCS;


typedef struct {
  int	nbnumber;
  int	n_v, n_elmt;	// number of vertices and number of elements
  int	*nv1, *nv2, *nv3;	// Node index of each triangle
  PetscReal	*nf_x, *nf_y, *nf_z;	// Normal direction
  PetscReal	*x_bp, *y_bp, *z_bp;	// Coordinates of IBM surface nodes
  PetscReal	*x_bp0, *y_bp0, *z_bp0;
  PetscReal	*x_bp_o, *y_bp_o, *z_bp_o;
/*   PetscReal	x_bp_in[101][3270], y_bp_in[101][3270], z_bp_in[101][3270]; */
  Cmpnts	*u, *uold;

  // xiaolei add
  PetscReal 	*uu, *vv, *ww, *uv, *vw, *wu, *U, *V, *W, *val; 
  PetscReal	*cent_x, *cent_y, *cent_z, *dA;
  int        	*i_min, *i_max, *j_min, *j_max, *k_min, *k_max;
  PetscReal 	*dhx, *dhy, *dhz; 
  // end add
  
} IBMNodes;

typedef struct {
  int IM, JM, KM; // dimensions of grid
  DM da;	/* Data structure for scalars (include the grid geometry
		   informaion, to obtain the grid information, 
		   use DMDAGetCoordinates) */
  DM fda, fda2;	// Data Structure for vectors
  DMDALocalInfo info;

  Vec	Cent;	// Coordinates of cell centers
  Vec 	Csi, Eta, Zet, Aj;
  Vec 	ICsi, IEta, IZet, IAj;
  Vec 	JCsi, JEta, JZet, JAj;
  Vec 	KCsi, KEta, KZet, KAj;
  Vec 	Ucont;	// Contravariant velocity components
  Vec 	Ucat;	// Cartesian velocity components
  Vec	Ucat_o;
  Vec 	P;
  Vec	Phi;
  Vec	GridSpace;
  Vec	Nvert;
  Vec	Nvert_o;
  BCS	Bcs;

  PetscInt	*nvert;//ody property

  PetscReal	ren;	// Reynolds number
  PetscReal	dt; 	// time step
  PetscReal	st;	// Strouhal number

  PetscReal	r[101], tin[101], uinr[101][1001];

  Vec	lUcont, lUcat, lP, lPhi;
  Vec	lCsi, lEta, lZet, lAj;
  Vec	lICsi, lIEta, lIZet, lIAj;
  Vec	lJCsi, lJEta, lJZet, lJAj;
  Vec	lKCsi, lKEta, lKZet, lKAj;
  Vec	lGridSpace;
  Vec	lNvert, lNvert_o;
  Vec	lCent;
  
  Vec Ucat_sum;		// u, v, w
  Vec Ucat_cross_sum;		// uv, vw, wu
  Vec Ucat_square_sum;	// u^2, v^2, w^2

  int _this;

  FlowWave *inflow, *kinematics;
  PetscInt number_flowwave, number_kinematics;

} UserCtx;

PetscErrorCode ReadCoordinates(UserCtx *user);
PetscErrorCode QCriteria(UserCtx *user);
PetscErrorCode Velocity_Magnitude(UserCtx *user);
PetscErrorCode Lambda2(UserCtx *user);
PetscErrorCode FormMetrics(UserCtx *user);
void Calc_avg_shear_stress(UserCtx *user);

// xiaolei add

int Itpwidth=4;
PetscInt NumberOfObjects=1;
PetscReal reflength_wt=1.0;
PetscReal z1;
PetscReal z2;
PetscInt  Nz12;
PetscInt  DiskAvg;
PetscReal xc_tb=0.0, yc_tb=0.0, r_tb=0.25;

PetscErrorCode disk_read_ucd(IBMNodes *ibm, int ibi, char fname[80], PetscReal zc_tb);

PetscErrorCode Pre_process(UserCtx *user, IBMNodes *ibm);

PetscErrorCode TECIOOut_AveragCirc(UserCtx *user, IBMNodes *ibm);
PetscErrorCode TECIOOut_MKE_Budget1(UserCtx *user);
PetscErrorCode TECIOOut_zMomentum_Budget(UserCtx *user);
// end add



int file_exist(char *str)
{
  int r=0;

  /*if(!my_rank)*/ {
    FILE *fp=fopen(str, "r");
    if(!fp) {
      r=0;
      printf("\nFILE !!! %s does not exist !!!\n", str);
    }
    else {
      fclose(fp);
      r=1;
    }
  }
  MPI_Bcast(&r, 1, MPI_INT, 0, PETSC_COMM_WORLD);
  return r;
};

void Calculate_Covariant_metrics(double g[3][3], double G[3][3])
{
	/*
		| csi.x  csi.y csi.z |-1		| x.csi  x.eta x.zet | 
		| eta.x eta.y eta.z |	 =	| y.csi   y.eta  y.zet |
		| zet.x zet.y zet.z |		| z.csi  z.eta z.zet |
	
	*/
	const double a11=g[0][0], a12=g[0][1], a13=g[0][2];
	const double a21=g[1][0], a22=g[1][1], a23=g[1][2];
	const double a31=g[2][0], a32=g[2][1], a33=g[2][2];

	double det= a11*(a33*a22-a32*a23)-a21*(a33*a12-a32*a13)+a31*(a23*a12-a22*a13);
	
	G[0][0] = (a33*a22-a32*a23)/det,	G[0][1] = - (a33*a12-a32*a13)/det, 	G[0][2] = (a23*a12-a22*a13)/det;
	G[1][0] = -(a33*a21-a31*a23)/det, G[1][1] = (a33*a11-a31*a13)/det,	G[1][2] = - (a23*a11-a21*a13)/det;
	G[2][0] = (a32*a21-a31*a22)/det,	G[2][1] = - (a32*a11-a31*a12)/det,	G[2][2] = (a22*a11-a21*a12)/det;
};

void Calculate_normal(Cmpnts csi, Cmpnts eta, Cmpnts zet, double ni[3], double nj[3], double nk[3])
{
	double g[3][3];
	double G[3][3];
	
	g[0][0]=csi.x, g[0][1]=csi.y, g[0][2]=csi.z;
	g[1][0]=eta.x, g[1][1]=eta.y, g[1][2]=eta.z;
	g[2][0]=zet.x, g[2][1]=zet.y, g[2][2]=zet.z;
	
	Calculate_Covariant_metrics(g, G);
	double xcsi=G[0][0], ycsi=G[1][0], zcsi=G[2][0];
	double xeta=G[0][1], yeta=G[1][1], zeta=G[2][1];
	double xzet=G[0][2], yzet=G[1][2], zzet=G[2][2];
	      
	double nx_i = xcsi, ny_i = ycsi, nz_i = zcsi;
	double nx_j = xeta, ny_j = yeta, nz_j = zeta;
	double nx_k = xzet, ny_k = yzet, nz_k = zzet;
	      
	double sum_i=sqrt(nx_i*nx_i+ny_i*ny_i+nz_i*nz_i);
	double sum_j=sqrt(nx_j*nx_j+ny_j*ny_j+nz_j*nz_j);
	double sum_k=sqrt(nx_k*nx_k+ny_k*ny_k+nz_k*nz_k);

	nx_i /= sum_i, ny_i /= sum_i, nz_i /= sum_i;
	nx_j /= sum_j, ny_j /= sum_j, nz_j /= sum_j;
	nx_k /= sum_k, ny_k /= sum_k, nz_k /= sum_k;

	ni[0] = nx_i, ni[1] = ny_i, ni[2] = nz_i;
	nj[0] = nx_j, nj[1] = ny_j, nj[2] = nz_j;
	nk[0] = nx_k, nk[1] = ny_k, nk[2] = nz_k;
}

double Contravariant_Reynolds_stress(double uu, double uv, double uw, double vv, double vw, double ww,
	double csi0, double csi1, double csi2, double eta0, double eta1, double eta2)
{
	
	double A = uu*csi0*eta0 + vv*csi1*eta1 + ww*csi2*eta2 + uv * (csi0*eta1+csi1*eta0)	+ uw * (csi0*eta2+csi2*eta0) + vw * (csi1*eta2+csi2*eta1);
	double B = sqrt(csi0*csi0+csi1*csi1+csi2*csi2)*sqrt(eta0*eta0+eta1*eta1+eta2*eta2);
	
	return A/B;
}

void Compute_du_center (int i, int j, int k,  int mx, int my, int mz, Cmpnts ***ucat, PetscReal ***nvert, 
				double *dudc, double *dvdc, double *dwdc, 
				double *dude, double *dvde, double *dwde,
				double *dudz, double *dvdz, double *dwdz)
{
	if ((nvert[k][j][i+1])> 0.1 || (i==mx-2 && !i_periodic)) {
		*dudc = ( ucat[k][j][i].x - ucat[k][j][i-1].x );
		*dvdc = ( ucat[k][j][i].y - ucat[k][j][i-1].y );
		*dwdc = ( ucat[k][j][i].z - ucat[k][j][i-1].z );
	}
	else if ((nvert[k][j][i-1])> 0.1 || (i==1 && !i_periodic)) {
		*dudc = ( ucat[k][j][i+1].x - ucat[k][j][i].x );
		*dvdc = ( ucat[k][j][i+1].y - ucat[k][j][i].y );
		*dwdc = ( ucat[k][j][i+1].z - ucat[k][j][i].z );
	}
	else {
		if(i_periodic && i==1) {
			*dudc = ( ucat[k][j][i+1].x - ucat[k][j][mx-2].x ) * 0.5;
			*dvdc = ( ucat[k][j][i+1].y - ucat[k][j][mx-2].y ) * 0.5;
			*dwdc = ( ucat[k][j][i+1].z - ucat[k][j][mx-2].z ) * 0.5;
		}
		else if(i_periodic && i==mx-2) {
			*dudc = ( ucat[k][j][1].x - ucat[k][j][i-1].x ) * 0.5;
			*dvdc = ( ucat[k][j][1].y - ucat[k][j][i-1].y ) * 0.5;
			*dwdc = ( ucat[k][j][1].z - ucat[k][j][i-1].z ) * 0.5;
		}
		else {
			*dudc = ( ucat[k][j][i+1].x - ucat[k][j][i-1].x ) * 0.5;
			*dvdc = ( ucat[k][j][i+1].y - ucat[k][j][i-1].y ) * 0.5;
			*dwdc = ( ucat[k][j][i+1].z - ucat[k][j][i-1].z ) * 0.5;
		}
	}

	if ((nvert[k][j+1][i])> 0.1 || (j==my-2 && !j_periodic)) {
		*dude = ( ucat[k][j][i].x - ucat[k][j-1][i].x );
		*dvde = ( ucat[k][j][i].y - ucat[k][j-1][i].y );
		*dwde = ( ucat[k][j][i].z - ucat[k][j-1][i].z );
	}
	else if ((nvert[k][j-1][i])> 0.1 || (j==1 && !j_periodic)) {
		*dude = ( ucat[k][j+1][i].x - ucat[k][j][i].x );
		*dvde = ( ucat[k][j+1][i].y - ucat[k][j][i].y );
		*dwde = ( ucat[k][j+1][i].z - ucat[k][j][i].z );
	}
	else {
		if(j_periodic && j==1) {
			*dude = ( ucat[k][j+1][i].x - ucat[k][my-2][i].x ) * 0.5;
			*dvde = ( ucat[k][j+1][i].y - ucat[k][my-2][i].y ) * 0.5;
			*dwde = ( ucat[k][j+1][i].z - ucat[k][my-2][i].z ) * 0.5;
		}
		else if(j_periodic && j==my-2) {
			*dude = ( ucat[k][1][i].x - ucat[k][j-1][i].x ) * 0.5;
			*dvde = ( ucat[k][1][i].y - ucat[k][j-1][i].y ) * 0.5;
			*dwde = ( ucat[k][1][i].z - ucat[k][j-1][i].z ) * 0.5;
		}
		else {
			*dude = ( ucat[k][j+1][i].x - ucat[k][j-1][i].x ) * 0.5;
			*dvde = ( ucat[k][j+1][i].y - ucat[k][j-1][i].y ) * 0.5;
			*dwde = ( ucat[k][j+1][i].z - ucat[k][j-1][i].z ) * 0.5;
		}
	}

	if ((nvert[k+1][j][i])> 0.1 || (k==mz-2 && !k_periodic)) {
		*dudz = ( ucat[k][j][i].x - ucat[k-1][j][i].x );
		*dvdz = ( ucat[k][j][i].y - ucat[k-1][j][i].y );
		*dwdz = ( ucat[k][j][i].z - ucat[k-1][j][i].z );
	}
	else if ((nvert[k-1][j][i])> 0.1 || (k==1 && !k_periodic)) {
		*dudz = ( ucat[k+1][j][i].x - ucat[k][j][i].x );
		*dvdz = ( ucat[k+1][j][i].y - ucat[k][j][i].y );
		*dwdz = ( ucat[k+1][j][i].z - ucat[k][j][i].z );
	}
	else {
		if(k_periodic && k==1) {
			*dudz = ( ucat[k+1][j][i].x - ucat[mz-2][j][i].x ) * 0.5;
			*dvdz = ( ucat[k+1][j][i].y - ucat[mz-2][j][i].y ) * 0.5;
			*dwdz = ( ucat[k+1][j][i].z - ucat[mz-2][j][i].z ) * 0.5;
		}
		else if(k_periodic && k==mz-2) {
			*dudz = ( ucat[1][j][i].x - ucat[k-1][j][i].x ) * 0.5;
			*dvdz = ( ucat[1][j][i].y - ucat[k-1][j][i].y ) * 0.5;
			*dwdz = ( ucat[1][j][i].z - ucat[k-1][j][i].z ) * 0.5;
		}
		else {
			*dudz = ( ucat[k+1][j][i].x - ucat[k-1][j][i].x ) * 0.5;
			*dvdz = ( ucat[k+1][j][i].y - ucat[k-1][j][i].y ) * 0.5;
			*dwdz = ( ucat[k+1][j][i].z - ucat[k-1][j][i].z ) * 0.5;
		}
	}
}

void Compute_dscalar_center (int i, int j, int k,  int mx, int my, int mz, PetscReal ***p, PetscReal ***nvert, double *dpdc, double *dpde, double *dpdz)
{
	if ((nvert[k][j][i+1])> 0.1 || (i==mx-2 && !i_periodic)) {
		*dpdc = ( p[k][j][i] - p[k][j][i-1] );
	}
	else if ((nvert[k][j][i-1])> 0.1 || (i==1 && !i_periodic)) {
		*dpdc = ( p[k][j][i+1] - p[k][j][i] );
	}
	else {
		if(i_periodic && i==1) {
			*dpdc = ( p[k][j][i+1] - p[k][j][mx-2] ) * 0.5;
		}
		else if(i_periodic && i==mx-2) {
			*dpdc = ( p[k][j][1] - p[k][j][i-1] ) * 0.5;
		}
		else {
			*dpdc = ( p[k][j][i+1] - p[k][j][i-1] ) * 0.5;
		}
	}

	if ((nvert[k][j+1][i])> 0.1 || (j==my-2 && !j_periodic)) {
		*dpde = ( p[k][j][i] - p[k][j-1][i] );
	}
	else if ((nvert[k][j-1][i])> 0.1 || (j==1 && !j_periodic)) {
		*dpde = ( p[k][j+1][i] - p[k][j][i] );
	}
	else {
		if(j_periodic && j==1) {
			*dpde = ( p[k][j+1][i] - p[k][my-2][i] ) * 0.5;
		}
		else if(j_periodic && j==my-2) {
			*dpde = ( p[k][1][i] - p[k][j-1][i] ) * 0.5;
		}
		else {
			*dpde = ( p[k][j+1][i] - p[k][j-1][i] ) * 0.5;
		}
	}

	if ((nvert[k+1][j][i])> 0.1 || (k==mz-2 && !k_periodic)) {
		*dpdz = ( p[k][j][i] - p[k-1][j][i] );
	}
	else if ((nvert[k-1][j][i])> 0.1 || (k==1 && !k_periodic)) {
		*dpdz = ( p[k+1][j][i] - p[k][j][i] );
	}
	else {
		if(k_periodic && k==1) {
			*dpdz = ( p[k+1][j][i] - p[mz-2][j][i] ) * 0.5;
		}
		else if(k_periodic && k==mz-2) {
			*dpdz = ( p[1][j][i] - p[k-1][j][i] ) * 0.5;
		}
		else {
			*dpdz = ( p[k+1][j][i] - p[k-1][j][i] ) * 0.5;
		}
	}
}

void Compute_du_dxyz (	double csi0, double csi1, double csi2, double eta0, double eta1, double eta2, double zet0, double zet1, double zet2, double ajc,
					double dudc, double dvdc, double dwdc, double dude, double dvde, double dwde, double dudz, double dvdz, double dwdz,
					double *du_dx, double *dv_dx, double *dw_dx, double *du_dy, double *dv_dy, double *dw_dy, double *du_dz, double *dv_dz, double *dw_dz )
{
	*du_dx = (dudc * csi0 + dude * eta0 + dudz * zet0) * ajc;
	*du_dy = (dudc * csi1 + dude * eta1 + dudz * zet1) * ajc;
	*du_dz = (dudc * csi2 + dude * eta2 + dudz * zet2) * ajc;
	*dv_dx = (dvdc * csi0 + dvde * eta0 + dvdz * zet0) * ajc;
	*dv_dy = (dvdc * csi1 + dvde * eta1 + dvdz * zet1) * ajc;
	*dv_dz = (dvdc * csi2 + dvde * eta2 + dvdz * zet2) * ajc;
	*dw_dx = (dwdc * csi0 + dwde * eta0 + dwdz * zet0) * ajc;
	*dw_dy = (dwdc * csi1 + dwde * eta1 + dwdz * zet1) * ajc;	
	*dw_dz = (dwdc * csi2 + dwde * eta2 + dwdz * zet2) * ajc;
};

void Compute_dscalar_dxyz (double csi0, double csi1, double csi2, double eta0, double eta1, double eta2, double zet0, double zet1, double zet2, double ajc,
					double dpdc, double dpde, double dpdz, double *dp_dx, double *dp_dy, double *dp_dz)
{
	*dp_dx = (dpdc * csi0 + dpde * eta0 + dpdz * zet0) * ajc;
	*dp_dy = (dpdc * csi1 + dpde * eta1 + dpdz * zet1) * ajc;
	*dp_dz = (dpdc * csi2 + dpde * eta2 + dpdz * zet2) * ajc;
};


void IKavg(float *x, int xs, int xe, int ys, int ye, int zs, int ze, int mx, int my, int mz)
{
	int i, j, k;

	for (j=ys; j<ye-2; j++) {
		double iksum=0;
		int count=0;
		for (i=xs; i<xe-2; i++)
		for (k=zs; k<ze-2; k++) {
			iksum += x[k * (mx-2)*(my-2) + j*(mx-2) + i];
			count++;
		}
		double ikavg = iksum/(double)count;
		for (i=xs; i<xe-2; i++)
		for (k=zs; k<ze-2; k++) 	x[k * (mx-2)*(my-2) + j*(mx-2) + i] = ikavg;
	}
}

/*
	pi, pk : # of grid points correcsponding to the period
	conditional averaging
*/
void IKavg_c (float *x, int xs, int xe, int ys, int ye, int zs, int ze, int mx, int my, int mz)
{
  //int i, j, k;
	
  if(pi<=0) pi = (xe-xs-2); // no averaging
  if(pk<=0) pk = (ze-zs-2); // no averaging 
	
	int ni = (xe-xs-2) / pi;
	int nk = (ze-zs-2) / pk;

	std::vector< std::vector<float> > iksum (pk);
	
	for(int k=0; k<pk; k++) iksum[k].resize(pi);
	
	for (int j=ys; j<ye-2; j++) {
		
		for(int k=0; k<pk; k++) std::fill( iksum[k].begin(), iksum[k].end(), 0.0 );
		
		for (int i=xs; i<xe-2; i++)
		for (int k=zs; k<ze-2; k++) {
			iksum [ ( k - zs ) % pk ] [ ( i - xs ) % pi] += x [k * (mx-2)*(my-2) + j*(mx-2) + i];
		}

		for (int i=xs; i<xe-2; i++)
		for (int k=zs; k<ze-2; k++) x [k * (mx-2)*(my-2) + j*(mx-2) + i] = iksum [ ( k - zs ) % pk ] [ ( i - xs ) % pi ] / (ni*nk);
	}
}


void Kavg(float *x, int xs, int xe, int ys, int ye, int zs, int ze, int mx, int my, int mz)
{
	int i, j, k;

	for (j=ys; j<ye-2; j++) 
	for (i=xs; i<xe-2; i++) {
	  double ksum=0;
	  int count=0;
	  for (k=zs; k<ze-2; k++) {
		ksum += x[k * (mx-2)*(my-2) + j*(mx-2) + i];
		count++;
	  }
	  double kavg = ksum/(double)count;
	  for (k=zs; k<ze-2; k++) {
		x[k * (mx-2)*(my-2) + j*(mx-2) + i] = kavg;
	  }
	}
      
}

void Javg(float *x, int xs, int xe, int ys, int ye, int zs, int ze, int mx, int my, int mz)
{
	int i, j, k;
	
	for (k=zs; k<ze-2; k++)
	for (i=xs; i<xe-2; i++) {
	  double jsum=0;
	  int count=0;
	  for (j=ys; j<ye-2; j++) {
		jsum += x[k * (mx-2)*(my-2) + j*(mx-2) + i];
		count++;
	  }
	  double javg = jsum/(double)count;
	  for (j=ys; j<ye-2; j++) {
		x[k * (mx-2)*(my-2) + j*(mx-2) + i] = javg;
	  }
	}
      
}

void Iavg(float *x, int xs, int xe, int ys, int ye, int zs, int ze, int mx, int my, int mz)
{
	int i, j, k;
	
	for (k=zs; k<ze-2; k++) {
	  for (j=ys; j<ye-2; j++) {
	    double isum=0;
	    int count=0;
	    for (i=xs; i<xe-2; i++) {
	      isum += x[k * (mx-2)*(my-2) + j*(mx-2) + i];
	      count++;
	    }
	    double iavg = isum/(double)count;
	    for (i=xs; i<xe-2; i++) {
	      x[k * (mx-2)*(my-2) + j*(mx-2) + i] = iavg;
	    }
	  }
	}
}

void Iavg(Cmpnts ***x, int xs, int xe, int ys, int ye, int zs, int ze, int mx, int my, int mz)
{
	int i, j, k;
	
	for (k=zs; k<ze-2; k++) 
	for (j=ys; j<ye-2; j++) {
		Cmpnts isum, iavg;
		isum.x = isum.y = isum.z = 0;
	    
	  int count=0;
		for (i=xs; i<xe-2; i++) {
	     isum.x += x[k+1][j+1][i+1].x;
	     isum.y += x[k+1][j+1][i+1].y;
	     isum.z += x[k+1][j+1][i+1].z;
	     count++;
	  }
	  
	  iavg.x = isum.x /(double)count;
	  iavg.y = isum.y /(double)count;
	  iavg.z = isum.z /(double)count;
	    
		for (i=xs; i<xe-2; i++) {
	  	x[k+1][j+1][i+1] = iavg;
		}
	}
}

void Iavg(PetscReal ***x, int xs, int xe, int ys, int ye, int zs, int ze, int mx, int my, int mz)
{
	int i, j, k;
	
	for (k=zs; k<ze-2; k++) 
	for (j=ys; j<ye-2; j++) {
		double isum, iavg;
		isum = 0;
	    
	  int count=0;
		for (i=xs; i<xe-2; i++) {
	     isum += x[k+1][j+1][i+1];
	     count++;
	  }
	  iavg = isum /(double)count;
	    
		for (i=xs; i<xe-2; i++) {
	  	x[k+1][j+1][i+1] = iavg;
		}
	}
}

PetscErrorCode TECIOOut_V(UserCtx *user, int only_V)	// seokkoo
{
	PetscInt bi;

	char filen[80];
	sprintf(filen, "%sResult%06d.plt", prefix, ti);

	INTEGER4 I, Debug, VIsDouble, DIsDouble, III, IMax, JMax, KMax;
	VIsDouble = 0;
	DIsDouble = 0;
	Debug = 0;
	
	
	/*
	if(only_V)   {
		if(les) I = TECINI100((char*)"Flow", (char*)"X Y Z UU Cs", filen, (char*)".", &Debug, &VIsDouble);
		else if(only_V==2) I = TECINI100((char*)"Flow", (char*)"X Y Z UU", filen, (char*)".", &Debug, &VIsDouble);
		else I = TECINI100((char*)"Flow", (char*)"X Y Z UU Nv", filen, (char*)".", &Debug, &VIsDouble);
	}*/
	only_V=0;

	/*
	X Y Z U V W UU P Nv 
	K Omega Nut || Cs
	Level Fr
	*/
	
	char str_rans[256];
	char str_les[256];
	char str_levelset[256];
	char str_all[512];
	
	sprintf(str_rans, "");
	sprintf(str_les, "");
	sprintf(str_levelset, "");
	
	if(rans) sprintf(str_rans, "K Omega Nut ");
	else if(les) sprintf(str_les, "Cs ");
	
	if(levelset) sprintf(str_levelset, "Level Fr ");
	
	
	sprintf(str_all, "X Y Z U V W UU P Nv %s%s%s", str_rans, str_les, str_levelset);
	
	printf("%s\n", str_all);
	
	I = TECINI100((char*)"Flow", str_all, filen, (char*)".", &Debug, &VIsDouble);
	
	for (bi=0; bi<block_number; bi++) {
		DM da = user[bi].da, fda = user[bi].fda;
		DMDALocalInfo info = user[bi].info;

		PetscInt	xs = info.xs, xe = info.xs + info.xm;
		PetscInt  	ys = info.ys, ye = info.ys + info.ym;
		PetscInt	zs = info.zs, ze = info.zs + info.zm;
		PetscInt	mx = info.mx, my = info.my, mz = info.mz;
    
		PetscInt	lxs, lys, lzs, lxe, lye, lze;
		PetscInt	i, j, k;
		PetscReal	***aj;
		Cmpnts	***ucat, ***coor, ***ucat_o, ***csi, ***eta, ***zet;
		PetscReal	***p, ***nvert, ***level;
		Vec		Coor, zCoor, nCoor;
		VecScatter	ctx;
		Vec K_Omega;
		
		DMDAVecGetArray(user[bi].da, user[bi].Nvert, &nvert);
		DMDAVecGetArray(user[bi].fda, user[bi].Ucat, &ucat);

		INTEGER4	ZoneType=0, ICellMax=0, JCellMax=0, KCellMax=0;
		INTEGER4	IsBlock=1, NumFaceConnections=0, FaceNeighborMode=0;
		INTEGER4    ShareConnectivityFromZone=0;
		INTEGER4	LOC[40] = {1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0}; /* 1 is cell-centered 0 is node centered */
		
		
		/*************************/
		printf("mi=%d, mj=%d, mk=%d\n", (int)mx, (int)my, (int)mz);
		printf("xs=%d, xe=%d\n", (int)xs, (int)xe);
		printf("ys=%d, ye=%d\n", (int)ys, (int)ye);
		printf("zs=%d, ze=%d\n", (int)zs, (int)ze);
		//exit(0);
		
		i_begin = 1, i_end = mx-1;	// cross section in tecplot
		j_begin = 1, j_end = my-1;
		k_begin = 1, k_end = mz-1;
		
		
		
		xs = i_begin - 1, xe = i_end+1;
		ys = j_begin - 1, ye = j_end+1;
		zs = k_begin - 1, ze = k_end+1;
		
		printf("xs=%d, xe=%d\n", (int)xs, (int)xe);
		printf("ys=%d, ye=%d\n", (int)ys, (int)ye);
		printf("zs=%d, ze=%d\n", (int)zs, (int)ze);
		//exit(0);
		//xs=0, xe=nsection+1;
		/*************************/
		
		if (vc) {
			LOC[3]=0; LOC[4]=0; LOC[5]=0; LOC[6]=0;
		}
		else if(only_V) {
			LOC[4]=0; LOC[5]=0; LOC[6]=0;
		}
		/*
		IMax = mx-1;
		JMax = my-1;
		KMax = mz-1;
		*/
		IMax = i_end - i_begin + 1;
		JMax = j_end - j_begin + 1;
		KMax = k_end - k_begin + 1;

		I = TECZNE100((char*)"Block 1",
			&ZoneType, 	/* Ordered zone */
			&IMax,
			&JMax,
			&KMax,
			&ICellMax,
			&JCellMax,
			&KCellMax,
			&IsBlock,	/* ISBLOCK  BLOCK format */
			&NumFaceConnections,
			&FaceNeighborMode,
			LOC,
			NULL,
			&ShareConnectivityFromZone); /* No connectivity sharing */

		//III = (mx-1) * (my-1) * (mz-1);
		III = IMax*JMax*KMax;
		
		DMDAGetCoordinates(da, &Coor);
		DMDAVecGetArray(fda, Coor, &coor);

		float *x;
		x = new float [III];
		
		for (k=zs; k<ze-1; k++)
		for (j=ys; j<ye-1; j++)
		for (i=xs; i<xe-1; i++) {
			x[(k-zs) * IMax*JMax + (j-ys)*IMax + (i-xs)] = coor[k][j][i].x;
		}
		I = TECDAT100(&III, &x[0], &DIsDouble);

		for (k=zs; k<ze-1; k++)
		for (j=ys; j<ye-1; j++)
		for (i=xs; i<xe-1; i++) {
			x[(k-zs) * IMax*JMax + (j-ys)*IMax + (i-xs)] = coor[k][j][i].y;
		}
		I = TECDAT100(&III, &x[0], &DIsDouble);

		for (k=zs; k<ze-1; k++)
		for (j=ys; j<ye-1; j++)
		for (i=xs; i<xe-1; i++) {
			x[(k-zs) * IMax*JMax + (j-ys)*IMax + (i-xs)] = coor[k][j][i].z;
		}
		I = TECDAT100(&III, &x[0], &DIsDouble);
		
		DMDAVecRestoreArray(fda, Coor, &coor);
		delete []x;
    		
		if(!vc) {
			x = new float [(mx-1)*(my-1)*(mz-1)];
			DMDAVecGetArray(user[bi].fda, user[bi].Ucat_o, &ucat_o);
			for (k=0; k<mz-1; k++)
			for (j=0; j<my-1; j++)
			for (i=0; i<mx-1; i++) {
				ucat_o[k][j][i].x = 0.125 *
					(ucat[k][j][i].x + ucat[k][j][i+1].x +
					ucat[k][j+1][i].x + ucat[k][j+1][i+1].x +
					ucat[k+1][j][i].x + ucat[k+1][j][i+1].x +
					ucat[k+1][j+1][i].x + ucat[k+1][j+1][i+1].x);
				ucat_o[k][j][i].y = 0.125 *
					(ucat[k][j][i].y + ucat[k][j][i+1].y +
					ucat[k][j+1][i].y + ucat[k][j+1][i+1].y +
					ucat[k+1][j][i].y + ucat[k+1][j][i+1].y +
					ucat[k+1][j+1][i].y + ucat[k+1][j+1][i+1].y);
				ucat_o[k][j][i].z = 0.125 *
					(ucat[k][j][i].z + ucat[k][j][i+1].z +
					ucat[k][j+1][i].z + ucat[k][j+1][i+1].z +
					ucat[k+1][j][i].z + ucat[k+1][j][i+1].z +
					ucat[k+1][j+1][i].z + ucat[k+1][j+1][i+1].z);
			}
	      
			for (k=zs; k<ze-1; k++)
			for (j=ys; j<ye-1; j++)
			for (i=xs; i<xe-1; i++) {
				x[k * (mx-1)*(my-1) + j*(mx-1) + i] = ucat_o[k][j][i].x;
			}
			if(!only_V) I = TECDAT100(&III, &x[0], &DIsDouble);

			for (k=zs; k<ze-1; k++)
			for (j=ys; j<ye-1; j++)
			for (i=xs; i<xe-1; i++) {
				x[k * (mx-1)*(my-1) + j*(mx-1) + i] = ucat_o[k][j][i].y;
			}
			if(!only_V) I = TECDAT100(&III, &x[0], &DIsDouble);

			for (k=zs; k<ze-1; k++)
			for (j=ys; j<ye-1; j++)
			for (i=xs; i<xe-1; i++) {
				x[k * (mx-1)*(my-1) + j*(mx-1) + i] = ucat_o[k][j][i].z;
			}
		
			if(!only_V) I = TECDAT100(&III, &x[0], &DIsDouble);
	      
			for (k=zs; k<ze-1; k++)
			for (j=ys; j<ye-1; j++)
			for (i=xs; i<xe-1; i++) {
				x[k * (mx-1)*(my-1) + j*(mx-1) + i] = sqrt( ucat_o[k][j][i].x*ucat_o[k][j][i].x + ucat_o[k][j][i].y*ucat_o[k][j][i].y + ucat_o[k][j][i].z*ucat_o[k][j][i].z );
			}
			I = TECDAT100(&III, &x[0], &DIsDouble);

			DMDAVecRestoreArray(user[bi].fda, user[bi].Ucat_o, &ucat_o);
			delete []x;
		}
		else {
		  //x = new float [(mx-2)*(my-2)*(mz-2)];
			//III = (mx-2) * (my-2) * (mz-2);
			III = (IMax-1)*(JMax-1)*(KMax-1);
			x = new float [III];

			if(!only_V)  {
				for (k=zs; k<ze-2; k++)
				for (j=ys; j<ye-2; j++)
				for (i=xs; i<xe-2; i++) {
					x[ (k-zs) * (IMax-1)*(JMax-1) + (j-ys) * (IMax-1) + (i-xs)] = ucat[k+1][j+1][i+1].x;
				}
				I = TECDAT100(&III, &x[0], &DIsDouble);

				for (k=zs; k<ze-2; k++)
				for (j=ys; j<ye-2; j++)
				for (i=xs; i<xe-2; i++) {
					x[ (k-zs) * (IMax-1)*(JMax-1) + (j-ys) * (IMax-1) + (i-xs)] = ucat[k+1][j+1][i+1].y;
				}
				I = TECDAT100(&III, &x[0], &DIsDouble);

				for (k=zs; k<ze-2; k++)
				for (j=ys; j<ye-2; j++)
				for (i=xs; i<xe-2; i++) {
					x[ (k-zs) * (IMax-1)*(JMax-1) + (j-ys) * (IMax-1) + (i-xs)] = ucat[k+1][j+1][i+1].z;
				}
				I = TECDAT100(&III, &x[0], &DIsDouble);
			}
			
			for (k=zs; k<ze-2; k++)
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) {
				x[ (k-zs) * (IMax-1)*(JMax-1) + (j-ys) * (IMax-1) + (i-xs)] = sqrt(ucat[k+1][j+1][i+1].x*ucat[k+1][j+1][i+1].x+ucat[k+1][j+1][i+1].y*ucat[k+1][j+1][i+1].y+ucat[k+1][j+1][i+1].z*ucat[k+1][j+1][i+1].z);
			}
			I = TECDAT100(&III, &x[0], &DIsDouble);
			delete []x;
		}

		
		III = (IMax-1)*(JMax-1)*(KMax-1);
   		x = new float [III];
		
		DMDAVecGetArray(user[bi].da, user[bi].P, &p);
		for (k=zs; k<ze-2; k++)
		for (j=ys; j<ye-2; j++)
		for (i=xs; i<xe-2; i++) {
			x[ (k-zs) * (IMax-1)*(JMax-1) + (j-ys) * (IMax-1) + (i-xs)] = p[k+1][j+1][i+1];
		}
		I = TECDAT100(&III, &x[0], &DIsDouble);
		DMDAVecRestoreArray(user[bi].da, user[bi].P, &p);
		
		
		for (k=zs; k<ze-2; k++)	// Nv
		for (j=ys; j<ye-2; j++)
		for (i=xs; i<xe-2; i++) {
			x[ (k-zs) * (IMax-1)*(JMax-1) + (j-ys) * (IMax-1) + (i-xs)] = nvert[k+1][j+1][i+1];
		}
		I = TECDAT100(&III, &x[0], &DIsDouble);
			
		if(rans) {
			Cmpnts2 ***komega;
			DMCreateGlobalVector(user->fda2, &K_Omega);
			PetscViewer	viewer;
			sprintf(filen, "%skfield%06d_%1d.dat", prefix, ti, user->_this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoad(K_Omega,viewer);
			PetscViewerDestroy(&viewer);
			DMDAVecGetArray(user[bi].fda2, K_Omega, &komega);
			
			for (k=zs; k<ze-2; k++)
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) x[ (k-zs) * (IMax-1)*(JMax-1) + (j-ys) * (IMax-1) + (i-xs)] = komega[k+1][j+1][i+1].x;
			I = TECDAT100(&III, &x[0], &DIsDouble);
			
			for (k=zs; k<ze-2; k++)
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) x[ (k-zs) * (IMax-1)*(JMax-1) + (j-ys) * (IMax-1) + (i-xs)] = komega[k+1][j+1][i+1].y;
			I = TECDAT100(&III, &x[0], &DIsDouble);
			
			for (k=zs; k<ze-2; k++)
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) x[ (k-zs) * (IMax-1)*(JMax-1) + (j-ys) * (IMax-1) + (i-xs)] = komega[k+1][j+1][i+1].x/(komega[k+1][j+1][i+1].y+1.e-20);
			I = TECDAT100(&III, &x[0], &DIsDouble);
			
			DMDAVecRestoreArray(user[bi].fda2, K_Omega, &komega);
			VecDestroy(&K_Omega);
		}
		else if (les) {
			PetscReal ***cs;
			Vec Cs;
			
			sprintf(filen, "%scs_%06d_%1d.dat", prefix, ti, user->_this);
			
			DMCreateGlobalVector(user->da, &Cs);
			
			if(file_exist(filen)) {
				PetscViewer viewer;
				PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
				VecLoad(Cs,viewer);
				PetscViewerDestroy(&viewer);
			}
			else {
				printf("Cannot open %s! Setting to zero.\n", filen);
				VecSet(Cs, 0);
			}
			
			DMDAVecGetArray(user[bi].da, Cs, &cs);
			
			for (k=zs; k<ze-2; k++)
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) x[ (k-zs) * (IMax-1)*(JMax-1) + (j-ys) * (IMax-1) + (i-xs)] = cs[k+1][j+1][i+1];
			I = TECDAT100(&III, &x[0], &DIsDouble);
			
			DMDAVecRestoreArray(user[bi].da, Cs, &cs);
			VecDestroy(&Cs);
		}
		
		if(levelset) {
			PetscReal ***level;
			Vec Levelset;
			DMCreateGlobalVector(user->da, &Levelset);
			PetscViewer	viewer;
			sprintf(filen, "%slfield%06d_%1d.dat", prefix, ti, user->_this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoad(Levelset,viewer);
			PetscViewerDestroy(&viewer);
			DMDAVecGetArray(user[bi].da, Levelset, &level);
			
			// level
			for (k=zs; k<ze-2; k++)
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) x[ (k-zs) * (IMax-1)*(JMax-1) + (j-ys) * (IMax-1) + (i-xs)] = level[k+1][j+1][i+1];
			I = TECDAT100(&III, &x[0], &DIsDouble);
		
			// Froude number
			for (k=zs; k<ze-2; k++)
			for (i=xs; i<xe-2; i++) {
				double Fr=0; // Froude number
				double G=0;
				double min_phi=1.e10, max_phi=-1.e10;
				double sum_velocity_magnitude = 0;
				int count = 0;
					
				if( fabs(gravity_x)>1.e-10 ) G=fabs(gravity_x);
				else if( fabs(gravity_y)>1.e-10 ) G=fabs(gravity_y);
				else if( fabs(gravity_z)>1.e-10 ) G=fabs(gravity_z);

				for (j=ys; j<ye-2; j++) {
					if(level[k+1][j+1][i+1]>0 && nvert[k+1][j+1][i+1]<0.1) {
						min_phi = std::min (min_phi, level[k+1][j+1][i+1]);
						max_phi = std::max (max_phi, level[k+1][j+1][i+1]);
						sum_velocity_magnitude += sqrt(ucat[k+1][j+1][i+1].x*ucat[k+1][j+1][i+1].x+ucat[k+1][j+1][i+1].y*ucat[k+1][j+1][i+1].y+ucat[k+1][j+1][i+1].z*ucat[k+1][j+1][i+1].z);
						count += 1;
					}
				}
				
				if(count) {
					double UU = sum_velocity_magnitude / (double)count;
					double depth = max_phi - min_phi;
					Fr = UU / sqrt ( G * depth );
					if(count<=2 || Fr>1000) Fr=0;
				}
				
				for (j=ys; j<ye-2; j++) x[ (k-zs) * (IMax-1)*(JMax-1) + (j-ys) * (IMax-1) + (i-xs)] = Fr;
			}
			I = TECDAT100(&III, &x[0], &DIsDouble);
			
			DMDAVecRestoreArray(user[bi].da, Levelset, &level);
			VecDestroy(&Levelset);
		}
		
		delete []x;
		
		DMDAVecRestoreArray(user[bi].fda, user[bi].Ucat, &ucat);
		DMDAVecRestoreArray(user[bi].da, user[bi].Nvert, &nvert);
	}
	I = TECEND100();
	
	return 0;
}


PetscErrorCode TECIOOut_Averaging(UserCtx *user)	// seokkoo
{
	PetscInt bi;

	char filen[80];
	
	if(phase) sprintf(filen, "%sResult_Phase%03d-avg.plt", prefix, (int)phase);
	else if(pcr) sprintf(filen, "%sResult_PCR%06d-avg.plt", prefix, (int)ti);
	else sprintf(filen, "%sResult%06d-avg.plt", prefix, (int)ti);

	INTEGER4 I, Debug, VIsDouble, DIsDouble, III, IMax, JMax, KMax;
	VIsDouble = 0;
	DIsDouble = 0;
	Debug = 0;

	char str_avg3[256];
	char str_rans[256];
	char str_les[256];
	char str_levelset[256];
	//char str_pcr[256];
	char str_all[512];
	
	sprintf(str_avg3, "");
	sprintf(str_rans, "");
	sprintf(str_les, "");
	sprintf(str_levelset, "");
	//sprintf(str_pcr, "");
	
	/*
	if(rans) sprintf(str_rans, "K Omega Nut ");
	else if(les) sprintf(str_les, "Cs ");
	*/
	if(levelset) sprintf(str_levelset, "Level Level_RMS ");
	if(averaging_option==3 && !phase) sprintf(str_avg3, "P pp Vortx Vorty Vortz vortx2 vorty2 vortz2 Anisotropy ");
	if(phase) sprintf(str_avg3, "P pp ");
	//if(pcr) sprintf(str_pcr, "P_fluct");
	//if(ucr) sprintf(str_pcr, "u_fluct v_fluct w_fluct");
	
	sprintf(str_all, "X Y Z U V W uu vv ww uv vw uw Nv %s%s%s%s", str_avg3, str_rans, str_les, str_levelset);
	
	printf("%s\n", str_all);
	
	if(pcr) I = TECINI100((char*)"Averaging", "X Y Z P Velocity_Magnitude Nv",  filen, (char*)".",  &Debug,  &VIsDouble);
	else I = TECINI100((char*)"Averaging", str_all, filen, (char*)".", &Debug, &VIsDouble);

	for (bi=0; bi<block_number; bi++) {
		DM da = user[bi].da, fda = user[bi].fda;
		DMDALocalInfo info = user[bi].info;

		PetscInt	xs = info.xs, xe = info.xs + info.xm;
		PetscInt  	ys = info.ys, ye = info.ys + info.ym;
		PetscInt	zs = info.zs, ze = info.zs + info.zm;
		PetscInt	mx = info.mx, my = info.my, mz = info.mz;
    
		PetscInt	lxs, lys, lzs, lxe, lye, lze;
		PetscInt	i, j, k;
		PetscReal	***aj;
		Cmpnts	***ucat, ***coor, ***ucat_o;
		Cmpnts	***u2sum, ***u1sum,  ***usum;
		PetscReal	***p, ***nvert;
		Vec		Coor, zCoor, nCoor;
		//VecScatter ctx;

		Vec X, Y, Z, U, V, W;

		INTEGER4	ZoneType=0, ICellMax=0, JCellMax=0, KCellMax=0;
		INTEGER4	IsBlock=1, NumFaceConnections=0, FaceNeighborMode=0;
		INTEGER4    ShareConnectivityFromZone=0;
		INTEGER4	LOC[100] = {1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; 
		/* 1 is cell-centered   0 is node centered */

		IMax = mx-1; JMax = my-1; KMax = mz-1;

		I = TECZNE100((char*)"Block 1",
			&ZoneType, 	/* Ordered zone */
			&IMax,
			&JMax,
			&KMax,
			&ICellMax,
			&JCellMax,
			&KCellMax,
			&IsBlock,	/* ISBLOCK  BLOCK format */
			&NumFaceConnections,
			&FaceNeighborMode,
			LOC,
			NULL,
			&ShareConnectivityFromZone); /* No connectivity sharing */

		float *x;
	
		x = new float [(mx-1)*(my-1)*(mz-1)];
		III = (mx-1) * (my-1) * (mz-1);

		DMDAGetCoordinates(da, &Coor);
		DMDAVecGetArray(fda, Coor, &coor);

		// X
		for (k=zs; k<ze-1; k++)
		for (j=ys; j<ye-1; j++)
		for (i=xs; i<xe-1; i++) x[k * (mx-1)*(my-1) + j*(mx-1) + i] = coor[k][j][i].x;
		I = TECDAT100(&III, x, &DIsDouble);

		// Y
		for (k=zs; k<ze-1; k++)
		for (j=ys; j<ye-1; j++)
		for (i=xs; i<xe-1; i++)	x[k * (mx-1)*(my-1) + j*(mx-1) + i] = coor[k][j][i].y;
		I = TECDAT100(&III, x, &DIsDouble);

		// Z
		for (k=zs; k<ze-1; k++)
		for (j=ys; j<ye-1; j++)
		for (i=xs; i<xe-1; i++)	x[k * (mx-1)*(my-1) + j*(mx-1) + i] = coor[k][j][i].z;
		I = TECDAT100(&III, x, &DIsDouble);
	
		DMDAVecRestoreArray(fda, Coor, &coor);
    
		//delete []x;
		double N=(double)tis+1.0;
		if(phase) N = phase;
		//x = new float [(mx-2)*(my-2)*(mz-2)];

		III = (mx-2) * (my-2) * (mz-2);
	
		if(pcr)  {
			DMDAVecGetArray(user[bi].da, user[bi].P, &p);
			for (k=zs; k<ze-2; k++) 
			for (j=ys; j<ye-2; j++) 
			for (i=xs; i<xe-2; i++) x[k * (mx-2)*(my-2) + j*(mx-2) + i] = p[k+1][j+1][i+1];
			I = TECDAT100(&III, x, &DIsDouble);
			DMDAVecRestoreArray(user[bi].da, user[bi].P, &p);

			// Load ucat
			PetscViewer	viewer;
			sprintf(filen, "%sufield%06d_%1d.dat", prefix, ti, user->_this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoad(user[bi].Ucat,viewer);
			PetscViewerDestroy(&viewer);
			
			DMDAVecGetArray(user[bi].fda, user[bi].Ucat, &ucat);
			for (k=zs; k<ze-2; k++)
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) x[k * (mx-2)*(my-2) + j*(mx-2) + i] 
								= sqrt ( ucat[k+1][j+1][i+1].x*ucat[k+1][j+1][i+1].x + ucat[k+1][j+1][i+1].y*ucat[k+1][j+1][i+1].y + ucat[k+1][j+1][i+1].z*ucat[k+1][j+1][i+1].z );
			I = TECDAT100(&III, x, &DIsDouble);
			DMDAVecRestoreArray(user[bi].fda, user[bi].Ucat, &ucat);
		}
		else {
			PetscViewer viewer;
			char filen[128];
			
			DMCreateGlobalVector(user[bi].fda, &user[bi].Ucat_sum);
			DMCreateGlobalVector(user[bi].fda, &user[bi].Ucat_cross_sum);
			DMCreateGlobalVector(user[bi].fda, &user[bi].Ucat_square_sum);

			if(phase) {
				sprintf(filen, "%sphase%03d_su0_%06d_%1d.dat", prefix, (int)phase, (int)ti_phase, (int)user[bi]._this);
				PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
				VecLoad(user[bi].Ucat_sum,viewer);
				PetscViewerDestroy(&viewer);
				
				PetscPrintf(PETSC_COMM_WORLD, "Read %s\n", filen);

				sprintf(filen, "%sphase%03d_su1_%06d_%1d.dat", prefix, (int)phase, (int)ti_phase, (int)user[bi]._this);
				PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
				VecLoad(user[bi].Ucat_cross_sum,viewer);
				PetscViewerDestroy(&viewer);
				
				PetscPrintf(PETSC_COMM_WORLD, "Read %s\n", filen);
			
				sprintf(filen, "%sphase%03d_su2_%06d_%1d.dat", prefix, (int)phase, (int)ti_phase, (int)user[bi]._this);
				PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
				VecLoad(user[bi].Ucat_square_sum,viewer);
				PetscViewerDestroy(&viewer);
				
				PetscPrintf(PETSC_COMM_WORLD, "Read %s\n", filen);
			}
			else {
				sprintf(filen, "%ssu0_%06d_%1d.dat", prefix, ti, user[bi]._this);
				PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
				VecLoad(user[bi].Ucat_sum,viewer);
				PetscViewerDestroy(&viewer);

				sprintf(filen, "%ssu1_%06d_%1d.dat", prefix, ti, user[bi]._this);
				PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
				VecLoad(user[bi].Ucat_cross_sum,viewer);
				PetscViewerDestroy(&viewer);
			
				sprintf(filen, "%ssu2_%06d_%1d.dat", prefix, ti, user[bi]._this);
				PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
				VecLoad(user[bi].Ucat_square_sum,viewer);
				PetscViewerDestroy(&viewer);
			}
			
			
			DMDAVecGetArray(user[bi].fda, user[bi].Ucat_sum, &usum);
			DMDAVecGetArray(user[bi].fda, user[bi].Ucat_cross_sum, &u1sum);
			DMDAVecGetArray(user[bi].fda, user[bi].Ucat_square_sum, &u2sum);
			
			// U
			for (k=zs; k<ze-2; k++)
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) x[k * (mx-2)*(my-2) + j*(mx-2) + i] = usum[k+1][j+1][i+1].x/N;
			if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			I = TECDAT100(&III, x, &DIsDouble);

			// V
			for (k=zs; k<ze-2; k++)
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) x[k * (mx-2)*(my-2) + j*(mx-2) + i] = usum[k+1][j+1][i+1].y/N;
			if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			I = TECDAT100(&III, x, &DIsDouble);

			// W
			for (k=zs; k<ze-2; k++) 
			for (j=ys; j<ye-2; j++) 
			for (i=xs; i<xe-2; i++) x[k * (mx-2)*(my-2) + j*(mx-2) + i] = usum[k+1][j+1][i+1].z/N;
			if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			I = TECDAT100(&III, x, &DIsDouble);
	      
			// uu, u rms
			for (k=zs; k<ze-2; k++)
			for (j=ys; j<ye-2; j++)
			for(i=xs; i<xe-2; i++) {
				double U = usum[k+1][j+1][i+1].x/N;
				x[k * (mx-2)*(my-2) + j*(mx-2) + i] = ( u2sum[k+1][j+1][i+1].x/N - U*U );
			}
			if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			I = TECDAT100(&III, x, &DIsDouble);
	      
			// vv, v rms
			for (k=zs; k<ze-2; k++)
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) {
				double V = usum[k+1][j+1][i+1].y/N;
				x[k * (mx-2)*(my-2) + j*(mx-2) + i] = ( u2sum[k+1][j+1][i+1].y/N - V*V );
			}
			if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			I = TECDAT100(&III, x, &DIsDouble);
	      
			// ww, w rms
			for (k=zs; k<ze-2; k++) 
			for (j=ys; j<ye-2; j++) 
			for (i=xs; i<xe-2; i++) {
				double W = usum[k+1][j+1][i+1].z/N;
				x[k * (mx-2)*(my-2) + j*(mx-2) + i] = ( u2sum[k+1][j+1][i+1].z/N - W*W );
			}
			if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			I = TECDAT100(&III, x, &DIsDouble);
	      
			// uv
			for (k=zs; k<ze-2; k++) 
			for (j=ys; j<ye-2; j++) 
			for (i=xs; i<xe-2; i++) {
				double UV = usum[k+1][j+1][i+1].x*usum[k+1][j+1][i+1].y / (N*N);
				x[k * (mx-2)*(my-2) + j*(mx-2) + i] = u1sum[k+1][j+1][i+1].x/N - UV;
			}
			if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			I = TECDAT100(&III, x, &DIsDouble);
	      
			// vw
			for (k=zs; k<ze-2; k++)
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) {
				double VW = usum[k+1][j+1][i+1].y*usum[k+1][j+1][i+1].z / (N*N);
				x[k * (mx-2)*(my-2) + j*(mx-2) + i] = u1sum[k+1][j+1][i+1].y/N - VW;
			}
			if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			I = TECDAT100(&III, x, &DIsDouble);     
	      
			// wu
			for (k=zs; k<ze-2; k++) 
			for (j=ys; j<ye-2; j++) 
			for (i=xs; i<xe-2; i++) {
				double WU = usum[k+1][j+1][i+1].z*usum[k+1][j+1][i+1].x / (N*N);
				x[k * (mx-2)*(my-2) + j*(mx-2) + i] = u1sum[k+1][j+1][i+1].z/N - WU;
			}
			if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
			I = TECDAT100(&III, x, &DIsDouble);
			
			// Nv
			DMDAVecGetArray(user[bi].da, user[bi].Nvert, &nvert);
			for (k=zs; k<ze-2; k++)
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) x[k * (mx-2)*(my-2) + j*(mx-2) + i] = nvert[k+1][j+1][i+1];
			I = TECDAT100(&III, x, &DIsDouble);
			DMDAVecRestoreArray(user[bi].da, user[bi].Nvert, &nvert);
			
			if(averaging_option==3 || phase) // P pp
			{
				Vec P_sum, P_square_sum;
				PetscReal ***psum, ***p2sum;
				
				DMCreateGlobalVector(user[bi].da, &P_sum);
				DMCreateGlobalVector(user[bi].da, &P_square_sum);
				
				if(phase) {
					sprintf(filen, "%sphase%03d_sp_%06d_%1d.dat", prefix, (int)phase, (int)ti_phase, (int)user[bi]._this);
					PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
					VecLoad(P_sum,viewer);
					PetscViewerDestroy(&viewer);
					
					PetscPrintf(PETSC_COMM_WORLD, "Read %s\n", filen);
					
					sprintf(filen, "%sphase%03d_sp2_%06d_%1d.dat", prefix, (int)phase, (int)ti_phase, (int)user[bi]._this);
					PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
					VecLoad(P_square_sum,viewer);
					PetscViewerDestroy(&viewer);
					
					PetscPrintf(PETSC_COMM_WORLD, "Read %s\n", filen);
				}
				else {
					sprintf(filen, "%ssp_%06d_%1d.dat", prefix, ti, user[bi]._this);
					PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
					VecLoad(P_sum,viewer);
					PetscViewerDestroy(&viewer);

					sprintf(filen, "%ssp2_%06d_%1d.dat", prefix, ti, user[bi]._this);
					PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
					VecLoad(P_square_sum,viewer);
					PetscViewerDestroy(&viewer);
				}
				
				DMDAVecGetArray(user[bi].da, P_sum, &psum);
				DMDAVecGetArray(user[bi].da, P_square_sum, &p2sum);
				
				for (k=zs; k<ze-2; k++)
				for (j=ys; j<ye-2; j++)
				for (i=xs; i<xe-2; i++) {
					double P = psum[k+1][j+1][i+1]/N;
					x[k * (mx-2)*(my-2) + j*(mx-2) + i] = P;
				}
				if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				I = TECDAT100(&III, x, &DIsDouble);

				
				for (k=zs; k<ze-2; k++)
				for (j=ys; j<ye-2; j++)
				for (i=xs; i<xe-2; i++) {
					double P = psum[k+1][j+1][i+1]/N;
					x[k * (mx-2)*(my-2) + j*(mx-2) + i] = ( p2sum[k+1][j+1][i+1]/N - P*P );
				}
				if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				I = TECDAT100(&III, x, &DIsDouble);

				DMDAVecRestoreArray(user[bi].da, P_sum, &psum);
				DMDAVecRestoreArray(user[bi].da, P_square_sum, &p2sum);
				
				VecDestroy(&P_sum);
				VecDestroy(&P_square_sum);
			}

			if(averaging_option==3 && !phase) { // Vortx Vorty Vortz vortx2 vorty2 vortz2 Anisotropy
				Vec Vort_sum, Vort_square_sum;
				Cmpnts ***vortsum, ***vort2sum;

				DMCreateGlobalVector(user[bi].fda, &Vort_sum);
				DMCreateGlobalVector(user[bi].fda, &Vort_square_sum);

				sprintf(filen, "%ssvo_%06d_%1d.dat", prefix, ti, user[bi]._this);
				PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
				VecLoad(Vort_sum,viewer);
				PetscViewerDestroy(&viewer);

				sprintf(filen, "%ssvo2_%06d_%1d.dat", prefix, ti, user[bi]._this);
				PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
				VecLoad(Vort_square_sum,viewer);
				PetscViewerDestroy(&viewer);
								
				DMDAVecGetArray(user[bi].fda, Vort_sum, &vortsum);
				DMDAVecGetArray(user[bi].fda, Vort_square_sum, &vort2sum);
                        
				for (k=zs; k<ze-2; k++)
				for (j=ys; j<ye-2; j++)
				for (i=xs; i<xe-2; i++) {
					double vortx = vortsum[k+1][j+1][i+1].x/N;
					x[k * (mx-2)*(my-2) + j*(mx-2) + i] = vortx;
				}
				if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				I = TECDAT100(&III, x, &DIsDouble);
				
				for (k=zs; k<ze-2; k++)
				for (j=ys; j<ye-2; j++)
				for (i=xs; i<xe-2; i++) {
					double vorty = vortsum[k+1][j+1][i+1].y/N;
					x[k * (mx-2)*(my-2) + j*(mx-2) + i] = vorty;
				}
				if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				I = TECDAT100(&III, x, &DIsDouble);
				
				for (k=zs; k<ze-2; k++)
				for (j=ys; j<ye-2; j++)
				for (i=xs; i<xe-2; i++) {
					double vortz = vortsum[k+1][j+1][i+1].z/N;
					x[k * (mx-2)*(my-2) + j*(mx-2) + i] = vortz;
				}
				if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				I = TECDAT100(&III, x, &DIsDouble);

				for (k=zs; k<ze-2; k++)
				for (j=ys; j<ye-2; j++)
				for (i=xs; i<xe-2; i++) {
					double vortx = vortsum[k+1][j+1][i+1].x/N;
					x[k * (mx-2)*(my-2) + j*(mx-2) + i] = ( vort2sum[k+1][j+1][i+1].x/N - vortx*vortx );
				}
				if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				I = TECDAT100(&III, x, &DIsDouble);

				for (k=zs; k<ze-2; k++)
				for (j=ys; j<ye-2; j++)
				for (i=xs; i<xe-2; i++) {
					double vorty = vortsum[k+1][j+1][i+1].y/N;
					x[k * (mx-2)*(my-2) + j*(mx-2) + i] = ( vort2sum[k+1][j+1][i+1].y/N - vorty*vorty );
				}
				if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				I = TECDAT100(&III, x, &DIsDouble);
				
				for (k=zs; k<ze-2; k++)
				for (j=ys; j<ye-2; j++)
				for (i=xs; i<xe-2; i++) {
					double vortz = vortsum[k+1][j+1][i+1].z/N;
					x[k * (mx-2)*(my-2) + j*(mx-2) + i] = ( vort2sum[k+1][j+1][i+1].z/N - vortz*vortz );
				}
				if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				I = TECDAT100(&III, x, &DIsDouble);
				
				
				// anisotropy 
				DMDAVecGetArray(user[bi].da, user[bi].Nvert, &nvert);
				for (k=zs; k<ze-2; k++)
				for (j=ys; j<ye-2; j++)
				for (i=xs; i<xe-2; i++) {
					double U = usum[k+1][j+1][i+1].x/N;
					double V = usum[k+1][j+1][i+1].y/N;
					double W = usum[k+1][j+1][i+1].z/N;
					double UV = usum[k+1][j+1][i+1].x*usum[k+1][j+1][i+1].y / (N*N);
					double VW = usum[k+1][j+1][i+1].y*usum[k+1][j+1][i+1].z / (N*N);
					double WU = usum[k+1][j+1][i+1].z*usum[k+1][j+1][i+1].x / (N*N);

					double uu = ( u2sum[k+1][j+1][i+1].x/N - U*U );
					double vv = ( u2sum[k+1][j+1][i+1].y/N - V*V );
					double ww = ( u2sum[k+1][j+1][i+1].z/N - W*W );
					double uv = u1sum[k+1][j+1][i+1].x/N - UV;
					double vw = u1sum[k+1][j+1][i+1].y/N - VW;
					double wu = u1sum[k+1][j+1][i+1].z/N - WU;
					double vu = uv, wv = vw, uw = wu;
					
					double TKE = 0.5 * (uu  + vv + ww);
					double b11 = 0.5 * uu/TKE - 1./3.;
					double b12 = 0.5 * uv/TKE;
					double b13 = 0.5 * uw/TKE;
					double b21 = b12;
					double b22 = 0.5 * vv/TKE - 1./3.;
					double b23 = 0.5 * vw/TKE;
					double b31 = b13;
					double b32 = b23;
					double b33 = 0.5 * ww/TKE - 1./3.;
					
					double PI = 0;
					PI += b11 * b11;
					PI += b12 * b12;
					PI += b13 * b13;
					PI += b21 * b21;
					PI += b22 * b22;
					PI += b23 * b23;
					PI += b31 * b31;
					PI += b32 * b32;
					PI += b33 * b33;
					PI *= -0.5;
					
					if( nvert[k+1][j+1][i+1] > 0.1 ) PI = 0;
					x[k * (mx-2)*(my-2) + j*(mx-2) + i] = PI;
				}
				if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				I = TECDAT100(&III, x, &DIsDouble);
				DMDAVecRestoreArray(user[bi].da, user[bi].Nvert, &nvert);

				DMDAVecRestoreArray(user[bi].fda, Vort_sum, &vortsum);
				DMDAVecRestoreArray(user[bi].fda, Vort_square_sum, &vort2sum);

				VecDestroy(&Vort_sum);
				VecDestroy(&Vort_square_sum);
			}
			
			DMDAVecRestoreArray(user[bi].fda, user[bi].Ucat_sum, &usum);
			DMDAVecRestoreArray(user[bi].fda, user[bi].Ucat_cross_sum, &u1sum);
			DMDAVecRestoreArray(user[bi].fda, user[bi].Ucat_square_sum, &u2sum);
			
			VecDestroy(&user[bi].Ucat_sum);
			VecDestroy(&user[bi].Ucat_cross_sum);
			VecDestroy(&user[bi].Ucat_square_sum);
			
			if(rans) {
			}
			
			if(les) {
			}
			
			if(levelset) // L Lrms
			{
				Vec L_sum, L_square_sum;
				PetscReal ***lsum, ***l2sum;
				
				DMCreateGlobalVector(user[bi].da, &L_sum);
				DMCreateGlobalVector(user[bi].da, &L_square_sum);
				
				if(phase) {
					sprintf(filen, "%sphase%03d_sl_%06d_%1d.dat", prefix, (int)phase, (int)ti_phase, (int)user[bi]._this);
					PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
					VecLoad(L_sum,viewer);
					PetscViewerDestroy(&viewer);
					
					PetscPrintf(PETSC_COMM_WORLD, "Read %s\n", filen);
					
					sprintf(filen, "%sphase%03d_sl2_%06d_%1d.dat", prefix, (int)phase, (int)ti_phase, (int)user[bi]._this);
					PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
					VecLoad(L_square_sum,viewer);
					PetscViewerDestroy(&viewer);
					
					PetscPrintf(PETSC_COMM_WORLD, "Read %s\n", filen);
				}
				else {
					sprintf(filen, "%ssl_%06d_%1d.dat", prefix, ti, user[bi]._this);
					PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
					VecLoad(L_sum,viewer);
					PetscViewerDestroy(&viewer);

					sprintf(filen, "%ssl2_%06d_%1d.dat", prefix, ti, user[bi]._this);
					PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
					VecLoad(L_square_sum,viewer);
					PetscViewerDestroy(&viewer);
				}
				
				DMDAVecGetArray(user[bi].da, L_sum, &lsum);
				DMDAVecGetArray(user[bi].da, L_square_sum, &l2sum);
				
				for (k=zs; k<ze-2; k++)
				for (j=ys; j<ye-2; j++)
				for (i=xs; i<xe-2; i++) {
					double L = lsum[k+1][j+1][i+1]/N;
					x[k * (mx-2)*(my-2) + j*(mx-2) + i] = L;
				}
				if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				I = TECDAT100(&III, x, &DIsDouble);

				
				for (k=zs; k<ze-2; k++)
				for (j=ys; j<ye-2; j++)
				for (i=xs; i<xe-2; i++) {
					double L = lsum[k+1][j+1][i+1]/N;
					x[k * (mx-2)*(my-2) + j*(mx-2) + i] = sqrt( l2sum[k+1][j+1][i+1]/N - L*L );	// rms
				}
				if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
				I = TECDAT100(&III, x, &DIsDouble);

				DMDAVecRestoreArray(user[bi].da, L_sum, &lsum);
				DMDAVecRestoreArray(user[bi].da, L_square_sum, &l2sum);
				
				VecDestroy(&L_sum);
				VecDestroy(&L_square_sum);
			}
		}
		delete []x;
	}
	I = TECEND100();
	return 0;
}


PetscErrorCode TECIOOut_TKE_Budget(UserCtx *user)
{
	PetscInt bi;

	char filen[80];
	
	if(phase) sprintf(filen, "%sTKEBudget_Phase%03d-avg.plt", prefix, (int)phase);
	else sprintf(filen, "%sTKEBudget%06d-avg.plt", prefix, (int)ti);

	printf("TKE Budget !\n");
	
	INTEGER4 I, Debug, VIsDouble, DIsDouble, III, IMax, JMax, KMax;
	VIsDouble = 0;
	DIsDouble = 0;
	Debug = 0;

	I = TECINI100((char*)"Averaging", "X Y Z TKE Production Dissipation ViscousDiffusion PressureTransport TurbulentConvection MeanConvection",  filen, (char*)".",  &Debug,  &VIsDouble);
	
	for (bi=0; bi<block_number; bi++) {
		DM da = user[bi].da, fda = user[bi].fda;
		DMDALocalInfo info = user[bi].info;

		PetscInt	xs = info.xs, xe = info.xs + info.xm;
		PetscInt  	ys = info.ys, ye = info.ys + info.ym;
		PetscInt	zs = info.zs, ze = info.zs + info.zm;
		PetscInt	mx = info.mx, my = info.my, mz = info.mz;
    
		PetscInt	lxs, lys, lzs, lxe, lye, lze;
		PetscInt	i, j, k;
		Cmpnts	***coor, ***u2sum, ***u1sum,  ***usum;
		Cmpnts	***ucat, ***csi, ***eta, ***zet; 
		PetscReal	***nvert, ***uu, ***vv, ***ww, ***uv, ***vw, ***wu, ***aj;
		Vec Coor;
		Vec Ucat, UU, VV, WW, UV, VW, WU;	// define vectors

		INTEGER4	ZoneType=0, ICellMax=0, JCellMax=0, KCellMax=0;
		INTEGER4	IsBlock=1, NumFaceConnections=0, FaceNeighborMode=0;
		INTEGER4    ShareConnectivityFromZone=0;
		INTEGER4	LOC[100] = {1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; 
		/* 1 is cell-centered   0 is node centered */

		IMax = mx-1; JMax = my-1; KMax = mz-1;

		I = TECZNE100((char*)"Block 1",
			&ZoneType, 	/* Ordered zone */
			&IMax,
			&JMax,
			&KMax,
			&ICellMax,
			&JCellMax,
			&KCellMax,
			&IsBlock,	/* ISBLOCK  BLOCK format */
			&NumFaceConnections,
			&FaceNeighborMode,
			LOC,
			NULL,
			&ShareConnectivityFromZone); /* No connectivity sharing */

		float *x;
	
		x = new float [(mx-1)*(my-1)*(mz-1)];
		III = (mx-1) * (my-1) * (mz-1);

		DMDAVecGetArray(user[bi].fda, user[bi].Csi, &csi);
		DMDAVecGetArray(user[bi].fda, user[bi].Eta, &eta);
		DMDAVecGetArray(user[bi].fda, user[bi].Zet, &zet);
		DMDAVecGetArray(user[bi].da, user[bi].Aj, &aj);
		DMDAVecGetArray(user[bi].da, user[bi].Nvert, &nvert);
				
		DMDAGetCoordinates(da, &Coor);
		DMDAVecGetArray(fda, Coor, &coor);
//
		// X
		for (k=zs; k<ze-1; k++)
		for (j=ys; j<ye-1; j++)
		for (i=xs; i<xe-1; i++) x[k * (mx-1)*(my-1) + j*(mx-1) + i] = coor[k][j][i].x;
		I = TECDAT100(&III, x, &DIsDouble);

		// Y
		for (k=zs; k<ze-1; k++)
		for (j=ys; j<ye-1; j++)
		for (i=xs; i<xe-1; i++)	x[k * (mx-1)*(my-1) + j*(mx-1) + i] = coor[k][j][i].y;
		I = TECDAT100(&III, x, &DIsDouble);

		// Z
		for (k=zs; k<ze-1; k++)
		for (j=ys; j<ye-1; j++)
		for (i=xs; i<xe-1; i++)	x[k * (mx-1)*(my-1) + j*(mx-1) + i] = coor[k][j][i].z;
		I = TECDAT100(&III, x, &DIsDouble);
	
		DMDAVecRestoreArray(fda, Coor, &coor);
    
		double N=(double)tis+1.0;
		if(phase) N = phase;
	
		III = (mx-2) * (my-2) * (mz-2);

		PetscViewer viewer;
		char filen[128];
		
		DMCreateGlobalVector(user[bi].fda, &user[bi].Ucat_sum);
		DMCreateGlobalVector(user[bi].fda, &user[bi].Ucat_cross_sum);
		DMCreateGlobalVector(user[bi].fda, &user[bi].Ucat_square_sum);
		
		DMCreateGlobalVector(user[bi].fda, &Ucat);
		DMCreateGlobalVector(user[bi].da, &UU);
		DMCreateGlobalVector(user[bi].da, &VV);
		DMCreateGlobalVector(user[bi].da, &WW);
		DMCreateGlobalVector(user[bi].da, &UV);
		DMCreateGlobalVector(user[bi].da, &VW);
		DMCreateGlobalVector(user[bi].da, &WU);
		/*
		if(phase) {
			sprintf(filen, "phase%03d_su0_%06d_%1d.dat", phase, ti_phase, user[bi]._this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoad(viewer, (user[bi].Ucat_sum));
			PetscViewerDestroy(&viewer);
				
			PetscPrintf(PETSC_COMM_WORLD, "Read %s\n", filen);

			sprintf(filen, "phase%03d_su1_%06d_%1d.dat", phase, ti_phase, user[bi]._this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoad(viewer, (user[bi].Ucat_cross_sum));
			PetscViewerDestroy(&viewer);
				
			PetscPrintf(PETSC_COMM_WORLD, "Read %s\n", filen);
			
			sprintf(filen, "phase%03d_su2_%06d_%1d.dat", phase, ti_phase, user[bi]._this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoad(viewer, (user[bi].Ucat_square_sum));
			PetscViewerDestroy(&viewer);
				
			PetscPrintf(PETSC_COMM_WORLD, "Read %s\n", filen);
		}
		else*/ 
		{
			sprintf(filen, "%ssu0_%06d_%1d.dat", prefix, ti, user[bi]._this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoad(user[bi].Ucat_sum,viewer);
			PetscViewerDestroy(&viewer);

			sprintf(filen, "%ssu1_%06d_%1d.dat", prefix, ti, user[bi]._this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoad(user[bi].Ucat_cross_sum,viewer);
			PetscViewerDestroy(&viewer);
			
			sprintf(filen, "%ssu2_%06d_%1d.dat", prefix, ti, user[bi]._this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoad(user[bi].Ucat_square_sum,viewer);
			PetscViewerDestroy(&viewer);
		}
			
		
		DMDAVecGetArray(user[bi].fda, user[bi].Ucat_cross_sum, &u1sum);
		DMDAVecGetArray(user[bi].fda, user[bi].Ucat_square_sum, &u2sum);
		
		DMDAVecGetArray(user[bi].da, UU, &uu);
		DMDAVecGetArray(user[bi].da, VV, &vv);
		DMDAVecGetArray(user[bi].da, WW, &ww);
		DMDAVecGetArray(user[bi].da, UV, &uv);
		DMDAVecGetArray(user[bi].da, VW, &vw);
		DMDAVecGetArray(user[bi].da, WU, &wu);
		
		DMDAVecGetArray(user[bi].fda, user[bi].Ucat_sum, &usum);
		// Reynolds stresses
		for (k=zs; k<ze-2; k++)
		for (j=ys; j<ye-2; j++)
		for (i=xs; i<xe-2; i++) {
			double _U = usum[k+1][j+1][i+1].x/N;
			double _V = usum[k+1][j+1][i+1].y/N;
			double _W = usum[k+1][j+1][i+1].z/N;
			double _UV = usum[k+1][j+1][i+1].x*usum[k+1][j+1][i+1].y / (N*N);
			double _VW = usum[k+1][j+1][i+1].y*usum[k+1][j+1][i+1].z / (N*N);
			double _WU = usum[k+1][j+1][i+1].z*usum[k+1][j+1][i+1].x / (N*N);
			
			uu[k+1][j+1][i+1] = ( u2sum[k+1][j+1][i+1].x/N - _U*_U );
			vv[k+1][j+1][i+1] = ( u2sum[k+1][j+1][i+1].y/N - _V*_V );
			ww[k+1][j+1][i+1] = ( u2sum[k+1][j+1][i+1].z/N - _W*_W );
			uv[k+1][j+1][i+1] = u1sum[k+1][j+1][i+1].x/N - _UV;
			vw[k+1][j+1][i+1] = u1sum[k+1][j+1][i+1].y/N - _VW;
			wu[k+1][j+1][i+1] = u1sum[k+1][j+1][i+1].z/N - _WU;
		}
		DMDAVecRestoreArray(user[bi].fda, user[bi].Ucat_sum, &usum);
		
					
		// Mean velocity
		VecAXPY(Ucat, 1./N, user[bi].Ucat_sum);
		
		// Destroy
		VecDestroy(&user[bi].Ucat_sum);
		
		
		DMDAVecGetArray(user[bi].fda, Ucat, &ucat);
		
		/* TKE */
		for (k=zs; k<ze-2; k++)
		for (j=ys; j<ye-2; j++)
		for (i=xs; i<xe-2; i++) {
			x[k * (mx-2)*(my-2) + j*(mx-2) + i] = 0.5 * ( uu[k+1][j+1][i+1] + vv[k+1][j+1][i+1] + ww[k+1][j+1][i+1] );
		}
		if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		I = TECDAT100(&III, x, &DIsDouble);
		
		/* Production */
		for (k=zs; k<ze-2; k++)
		for (j=ys; j<ye-2; j++)
		for (i=xs; i<xe-2; i++) {
			double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
			double eta0= eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
			double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
			double ajc = aj[k+1][j+1][i+1];
			double dudc, dude, dudz;
			double dvdc, dvde, dvdz;
			double dwdc, dwde, dwdz;
			double du_dx, du_dy, du_dz;
			double dv_dx, dv_dy, dv_dz;
			double dw_dx, dw_dy, dw_dz;
			
			Compute_du_center (i+1, j+1, k+1, mx, my, mz, ucat, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
			Compute_du_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
			
			double Prod = 0;
			Prod -= uu[k+1][j+1][i+1] * du_dx;
			Prod -= uv[k+1][j+1][i+1] * du_dy;
			Prod -= wu[k+1][j+1][i+1] * du_dz;
			Prod -= uv[k+1][j+1][i+1] * dv_dx;
			Prod -= vv[k+1][j+1][i+1] * dv_dy;
			Prod -= vw[k+1][j+1][i+1] * dv_dz;
			Prod -= wu[k+1][j+1][i+1] * dw_dx;
			Prod -= vw[k+1][j+1][i+1] * dw_dy;
			Prod -= ww[k+1][j+1][i+1] * dw_dz;

			x[k * (mx-2)*(my-2) + j*(mx-2) + i] = Prod;
		}
		if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		I = TECDAT100(&III, x, &DIsDouble);
		
		/* Dissipation */
		Vec dU2_sum;
		PetscReal ***du2sum;
				
		DMCreateGlobalVector(user[bi].da, &dU2_sum);
				
		sprintf(filen, "%ssu4_%06d_%1d.dat", prefix, ti, user[bi]._this);
		PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
		VecLoad(dU2_sum,viewer);
		PetscViewerDestroy(&viewer);
				
		DMDAVecGetArray(user[bi].da, dU2_sum, &du2sum);
				
		for (k=zs; k<ze-2; k++)
		for (j=ys; j<ye-2; j++)
		for (i=xs; i<xe-2; i++) {
			double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
			double eta0= eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
			double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
			double ajc = aj[k+1][j+1][i+1];
			double dudc, dude, dudz;
			double dvdc, dvde, dvdz;
			double dwdc, dwde, dwdz;
			double du_dx, du_dy, du_dz;
			double dv_dx, dv_dy, dv_dz;
			double dw_dx, dw_dy, dw_dz;
			
			Compute_du_center (i+1, j+1, k+1, mx, my, mz, ucat, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
			Compute_du_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
								
			double Diss = 0;
			Diss -= pow(du_dx, 2.) + pow(du_dy, 2.) + pow(du_dz, 2.);  // xiaolei change + to -
			Diss -= pow(dv_dx, 2.) + pow(dv_dy, 2.) + pow(dv_dz, 2.);
			Diss -= pow(dw_dx, 2.) + pow(dw_dy, 2.) + pow(dw_dz, 2.);
			Diss += du2sum[k+1][j+1][i+1] / N;
					
			Diss *= 1./user[bi].ren;
			Diss *= -1.;
			x[k * (mx-2)*(my-2) + j*(mx-2) + i] = Diss;
		}
		if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		I = TECDAT100(&III, x, &DIsDouble);
				
		DMDAVecRestoreArray(user[bi].da, dU2_sum, &du2sum);
		VecDestroy(&dU2_sum);
		
		/* Viscous Diffusion */
		
		// Make dk/dxj vector
		Vec GradK;
		Cmpnts ***gradk;
		
		DMCreateGlobalVector(user[bi].fda, &GradK);
		DMDAVecGetArray(user[bi].fda, GradK, &gradk);
		
		for (k=zs; k<ze-2; k++)
		for (j=ys; j<ye-2; j++)
		for (i=xs; i<xe-2; i++) {
			double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
			double eta0= eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
			double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
			double ajc = aj[k+1][j+1][i+1];
			double duudc, duude, duudz;
			double dvvdc, dvvde, dvvdz;
			double dwwdc, dwwde, dwwdz;
			double duu_dx, duu_dy, duu_dz;
			double dvv_dx, dvv_dy, dvv_dz;
			double dww_dx, dww_dy, dww_dz;
			double dk_dx, dk_dy, dk_dz;
			
			Compute_dscalar_center (i+1, j+1, k+1, mx, my, mz, uu, nvert, &duudc, &duude, &duudz);
			Compute_dscalar_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, duudc, duude, duudz, &duu_dx, &duu_dy, &duu_dz);
			
			Compute_dscalar_center (i+1, j+1, k+1, mx, my, mz, vv, nvert, &dvvdc, &dvvde, &dvvdz);
			Compute_dscalar_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dvvdc, dvvde, dvvdz, &dvv_dx, &dvv_dy, &dvv_dz);
			
			Compute_dscalar_center (i+1, j+1, k+1, mx, my, mz, ww, nvert, &dwwdc, &dwwde, &dwwdz);
			Compute_dscalar_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dwwdc, dwwde, dwwdz, &dww_dx, &dww_dy, &dww_dz);
			
			dk_dx = 0;
			dk_dx += 0.5 * duu_dx;
			dk_dx += 0.5 * dvv_dx;
			dk_dx += 0.5 * dww_dx;
			
			dk_dy = 0;
			dk_dy += 0.5 * duu_dy;
			dk_dy += 0.5 * dvv_dy;
			dk_dy += 0.5 * dww_dy;
			
			dk_dz = 0;
			dk_dz += 0.5 * duu_dz;
			dk_dz += 0.5 * dvv_dz;
			dk_dz += 0.5 * dww_dz;
			
			gradk[k+1][j+1][i+1].x = dk_dx;
			gradk[k+1][j+1][i+1].y = dk_dy;
			gradk[k+1][j+1][i+1].z = dk_dz;
		}
		
		// Viscous Diffusion
		for (k=zs; k<ze-2; k++)
		for (j=ys; j<ye-2; j++)
		for (i=xs; i<xe-2; i++) {
			double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
			double eta0= eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
			double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
			double ajc = aj[k+1][j+1][i+1];
			double dudc, dude, dudz;
			double dvdc, dvde, dvdz;
			double dwdc, dwde, dwdz;
			double du_dx, du_dy, du_dz;
			double dv_dx, dv_dy, dv_dz;
			double dw_dx, dw_dy, dw_dz;
			
			Compute_du_center (i+1, j+1, k+1, mx, my, mz, gradk, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
			Compute_du_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
			
			double ddk_dxdx = du_dx;
			double ddk_dydy = dv_dy;
			double ddk_dzdz = dw_dz;
			double nu = 1./user[bi].ren;
			
			double ViscousDiffusion = 0;
			
			ViscousDiffusion = nu * (ddk_dxdx);
			ViscousDiffusion += nu * (ddk_dydy);
			ViscousDiffusion += nu * (ddk_dzdz);
			
			x[k * (mx-2)*(my-2) + j*(mx-2) + i] = ViscousDiffusion;
		}
		if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		I = TECDAT100(&III, x, &DIsDouble);
		
		DMDAVecRestoreArray(user[bi].fda, GradK, &gradk);
		VecDestroy(&GradK);
		
		/* Pressure Transport */
		Vec Udp_sum;
		Vec P_sum;
				
		PetscReal ***psum;
		PetscReal ***udpsum;
				
		DMCreateGlobalVector(user[bi].da, &P_sum);
		DMCreateGlobalVector(user[bi].da, &Udp_sum);
				
		sprintf(filen, "%ssp_%06d_%1d.dat", prefix, ti, user[bi]._this);
		PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
		VecLoad(P_sum,viewer);
		PetscViewerDestroy(&viewer);
				
		sprintf(filen, "%ssu3_%06d_%1d.dat", prefix, ti, user[bi]._this);
		PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
		VecLoad(Udp_sum,viewer);
		PetscViewerDestroy(&viewer);
				
		DMDAVecGetArray(user[bi].da, P_sum, &psum);
		DMDAVecGetArray(user[bi].da, Udp_sum, &udpsum);
				
		for (k=zs; k<ze-2; k++)
		for (j=ys; j<ye-2; j++)
		for (i=xs; i<xe-2; i++) {
			double dpdc, dpde, dpdz;
			double dp_dx, dp_dy, dp_dz;
			double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
			double eta0= eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
			double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
			double ajc = aj[k+1][j+1][i+1];
					
			Compute_dscalar_center (i+1, j+1, k+1, mx, my, mz, psum, nvert, &dpdc, &dpde, &dpdz);
			Compute_dscalar_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dpdc, dpde, dpdz, &dp_dx, &dp_dy, &dp_dz);

			dp_dx /= N;
			dp_dy /= N;
			dp_dz /= N;
			
			double U = ucat[k+1][j+1][i+1].x;
			double V = ucat[k+1][j+1][i+1].y;
			double W = ucat[k+1][j+1][i+1].z;
			
			double PressureTrans = 0;
			PressureTrans -= U * dp_dx;
			PressureTrans -= V * dp_dy;
			PressureTrans -= W * dp_dz;
			PressureTrans += udpsum[k+1][j+1][i+1] / N;
					
			PressureTrans *= -1.;
			x[k * (mx-2)*(my-2) + j*(mx-2) + i] = PressureTrans;
		}
		if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		I = TECDAT100(&III, x, &DIsDouble);

		DMDAVecRestoreArray(user[bi].da, Udp_sum, &udpsum);
		DMDAVecRestoreArray(user[bi].da, P_sum, &psum);
				
		VecDestroy(&Udp_sum);
		VecDestroy(&P_sum);
		
		/* Turbulent Convection */
		Vec UUU_sum, UUU;
		Cmpnts ***uuusum;
		Cmpnts ***uuu;
		
		DMCreateGlobalVector(user[bi].fda, &UUU_sum);
		DMCreateGlobalVector(user[bi].fda, &UUU);
		
		sprintf(filen, "%ssu5_%06d_%1d.dat", prefix, ti, user[bi]._this);
		PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
		VecLoad(UUU_sum,viewer);
		PetscViewerDestroy(&viewer);
		
		DMDAVecGetArray(user[bi].fda, UUU_sum, &uuusum);
		DMDAVecGetArray(user[bi].fda, UUU, &uuu);
		
		// a new vector for <ui'ui'uk'>
		for (k=zs; k<ze-2; k++)
		for (j=ys; j<ye-2; j++)
		for (i=xs; i<xe-2; i++) {
			double UU = u2sum[k+1][j+1][i+1].x/N;
			double VV = u2sum[k+1][j+1][i+1].y/N;
			double WW = u2sum[k+1][j+1][i+1].z/N;
			double UV = u1sum[k+1][j+1][i+1].x/N;
			double VW = u1sum[k+1][j+1][i+1].y/N;
			double WU = u1sum[k+1][j+1][i+1].z/N;
			
			double _U = ucat[k+1][j+1][i+1].x;
			double _V = ucat[k+1][j+1][i+1].y;
			double _W = ucat[k+1][j+1][i+1].z;
			
			uuu[k+1][j+1][i+1].x = ( UU + VV + WW ) * _U - 2. * ( _U*_U + _V*_V + _W*_W ) * _U - uuusum[k+1][j+1][i+1].x/N + 2. * ( UU * _U + UV * _V + WU * _W );
			uuu[k+1][j+1][i+1].y = ( UU + VV + WW ) * _V - 2. * ( _U*_U + _V*_V + _W*_W ) * _V - uuusum[k+1][j+1][i+1].y/N + 2. * ( UV * _U + VV * _V + VW * _W );
			uuu[k+1][j+1][i+1].z = ( UU + VV + WW ) * _W - 2. * ( _U*_U + _V*_V + _W*_W ) * _W - uuusum[k+1][j+1][i+1].z/N + 2. * ( WU * _U + VW * _V + WW * _W );
		}
		
		// calculate the divergence of <uuu>
		for (k=zs; k<ze-2; k++)
		for (j=ys; j<ye-2; j++)
		for (i=xs; i<xe-2; i++) {
			double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
			double eta0= eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
			double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
			double ajc = aj[k+1][j+1][i+1];
			double dudc, dude, dudz;
			double dvdc, dvde, dvdz;
			double dwdc, dwde, dwdz;
			double du_dx, du_dy, du_dz;
			double dv_dx, dv_dy, dv_dz;
			double dw_dx, dw_dy, dw_dz;
			
			Compute_du_center (i+1, j+1, k+1, mx, my, mz, uuu, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
			Compute_du_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
			
			double TurbulentConvection = du_dx + dv_dy + dw_dz;
			TurbulentConvection *= -1.0 * 0.5;
			
			x[k * (mx-2)*(my-2) + j*(mx-2) + i] = TurbulentConvection;
		}
		if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		I = TECDAT100(&III, x, &DIsDouble);
		
		DMDAVecRestoreArray(user[bi].fda, UUU_sum, &uuusum);
		DMDAVecRestoreArray(user[bi].fda, UUU, &uuu);
		VecDestroy(&UUU_sum);
		VecDestroy(&UUU);
		
		/* Mean Convection */
		for (k=zs; k<ze-2; k++)
		for (j=ys; j<ye-2; j++)
		for (i=xs; i<xe-2; i++) {
			double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
			double eta0= eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
			double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
			double ajc = aj[k+1][j+1][i+1];
			double duudc, duude, duudz;
			double dvvdc, dvvde, dvvdz;
			double dwwdc, dwwde, dwwdz;
			double duu_dx, duu_dy, duu_dz;
			double dvv_dx, dvv_dy, dvv_dz;
			double dww_dx, dww_dy, dww_dz;
			double dk_dx, dk_dy, dk_dz;
			
			Compute_dscalar_center (i+1, j+1, k+1, mx, my, mz, uu, nvert, &duudc, &duude, &duudz);
			Compute_dscalar_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, duudc, duude, duudz, &duu_dx, &duu_dy, &duu_dz);
			
			Compute_dscalar_center (i+1, j+1, k+1, mx, my, mz, vv, nvert, &dvvdc, &dvvde, &dvvdz);
			Compute_dscalar_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dvvdc, dvvde, dvvdz, &dvv_dx, &dvv_dy, &dvv_dz);
			
			Compute_dscalar_center (i+1, j+1, k+1, mx, my, mz, ww, nvert, &dwwdc, &dwwde, &dwwdz);
			Compute_dscalar_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dwwdc, dwwde, dwwdz, &dww_dx, &dww_dy, &dww_dz);
			
			dk_dx = 0;
			dk_dx += 0.5 * duu_dx;
			dk_dx += 0.5 * dvv_dx;
			dk_dx += 0.5 * dww_dx;
			
			dk_dy = 0;
			dk_dy += 0.5 * duu_dy;
			dk_dy += 0.5 * dvv_dy;
			dk_dy += 0.5 * dww_dy;
			
			dk_dz = 0;
			dk_dz += 0.5 * duu_dz;
			dk_dz += 0.5 * dvv_dz;
			dk_dz += 0.5 * dww_dz;
			
			double MeanConvection = 0;
			
			MeanConvection += ucat[k][j][i].x * dk_dx;
			MeanConvection += ucat[k][j][i].y * dk_dy;
			MeanConvection += ucat[k][j][i].z * dk_dz;
			MeanConvection *= -1.0;

			x[k * (mx-2)*(my-2) + j*(mx-2) + i] = MeanConvection;
		}
		if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		I = TECDAT100(&III, x, &DIsDouble);
		/*
		// Write Nv
		for (k=zs; k<ze-2; k++)
		for (j=ys; j<ye-2; j++)
		for (i=xs; i<xe-2; i++) x[k * (mx-2)*(my-2) + j*(mx-2) + i] = nvert[k+1][j+1][i+1];
		I = TECDAT100(&III, x, &DIsDouble);
		*/
		
		DMDAVecRestoreArray(user[bi].fda, Ucat, &ucat);
		
		DMDAVecRestoreArray(user[bi].da, UU, &uu);
		DMDAVecRestoreArray(user[bi].da, VV, &vv);
		DMDAVecRestoreArray(user[bi].da, WW, &ww);
		DMDAVecRestoreArray(user[bi].da, UV, &uv);
		DMDAVecRestoreArray(user[bi].da, VW, &vw);
		DMDAVecRestoreArray(user[bi].da, WU, &wu);
				
		DMDAVecRestoreArray(user[bi].fda, user[bi].Csi, &csi);
		DMDAVecRestoreArray(user[bi].fda, user[bi].Eta, &eta);
		DMDAVecRestoreArray(user[bi].fda, user[bi].Zet, &zet);
		DMDAVecRestoreArray(user[bi].da, user[bi].Aj, &aj);
		DMDAVecRestoreArray(user[bi].da, user[bi].Nvert, &nvert);
		
		DMDAVecRestoreArray(user[bi].fda, user[bi].Ucat_cross_sum, &u1sum);
		DMDAVecRestoreArray(user[bi].fda, user[bi].Ucat_square_sum, &u2sum);
		
		delete []x;
		
		VecDestroy(&user[bi].Ucat_cross_sum);
		VecDestroy(&user[bi].Ucat_square_sum);
		VecDestroy(&Ucat);
		VecDestroy(&UU);
		VecDestroy(&VV);
		VecDestroy(&WW);
		VecDestroy(&UV);
		VecDestroy(&VW);
		VecDestroy(&WU);
	}
	
	I = TECEND100();
	
	return 0;
}


PetscErrorCode TECIOOut_Vorticity_Budget(UserCtx *user)	// x-direction vorticity
{
	PetscInt bi;

	char filen[80];
	
	if(phase) sprintf(filen, "%sVorticityBudget_Phase%03d-avg.plt", prefix, (int)phase);
	else sprintf(filen, "%sVorticityBudget%06d-avg.plt", prefix, (int)ti);

	printf("X-Vorticity Budget !\n");
	
	INTEGER4 I, Debug, VIsDouble, DIsDouble, III, IMax, JMax, KMax;
	VIsDouble = 0;
	DIsDouble = 0;
	Debug = 0;

	
	//I = TECINI100((char*)"Averaging", "X Y Z Vortx Convec Diff Stretch P1 P2 P3 P4 Nv",  filen, (char*)".",  &Debug,  &VIsDouble);
	I = TECINI100((char*)"Averaging", "X Y Z Vortx Convection-x Convection-y Convection-z Diffusion Stretching Skewing-y Skewing-z P2 P3 P4",  filen, (char*)".",  &Debug,  &VIsDouble);
	
	
	for (bi=0; bi<block_number; bi++) {
		DM da = user[bi].da, fda = user[bi].fda;
		DMDALocalInfo info = user[bi].info;

		PetscInt	xs = info.xs, xe = info.xs + info.xm;
		PetscInt  	ys = info.ys, ye = info.ys + info.ym;
		PetscInt	zs = info.zs, ze = info.zs + info.zm;
		PetscInt	mx = info.mx, my = info.my, mz = info.mz;
    
		PetscInt	lxs, lys, lzs, lxe, lye, lze;
		PetscInt	i, j, k;
		Cmpnts	***coor, ***u2sum, ***u1sum,  ***usum;
		Cmpnts	***ucat, ***csi, ***eta, ***zet; 
		Cmpnts 	***vort;
		PetscReal	***nvert, ***uu, ***vv, ***ww, ***uv, ***vw, ***wu, ***aj;
		
		Vec Coor;
		Vec Ucat, UU, VV, WW, UV, VW, WU;	// define vectors
		Vec Vort;

		INTEGER4	ZoneType=0, ICellMax=0, JCellMax=0, KCellMax=0;
		INTEGER4	IsBlock=1, NumFaceConnections=0, FaceNeighborMode=0;
		INTEGER4    ShareConnectivityFromZone=0;
		INTEGER4	LOC[100] = {1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; 
		/* 1 is cell-centered   0 is node centered */

		IMax = mx-1; JMax = my-1; KMax = mz-1;

		I = TECZNE100((char*)"Block 1",
			&ZoneType, 	/* Ordered zone */
			&IMax,
			&JMax,
			&KMax,
			&ICellMax,
			&JCellMax,
			&KCellMax,
			&IsBlock,	/* ISBLOCK  BLOCK format */
			&NumFaceConnections,
			&FaceNeighborMode,
			LOC,
			NULL,
			&ShareConnectivityFromZone); /* No connectivity sharing */

		float *x;
	
		x = new float [(mx-1)*(my-1)*(mz-1)];
		III = (mx-1) * (my-1) * (mz-1);

		DMDAVecGetArray(user[bi].fda, user[bi].Csi, &csi);
		DMDAVecGetArray(user[bi].fda, user[bi].Eta, &eta);
		DMDAVecGetArray(user[bi].fda, user[bi].Zet, &zet);
		DMDAVecGetArray(user[bi].da, user[bi].Aj, &aj);
				
		DMDAGetCoordinates(da, &Coor);
		DMDAVecGetArray(fda, Coor, &coor);
//
		// X
		for (k=zs; k<ze-1; k++)
		for (j=ys; j<ye-1; j++)
		for (i=xs; i<xe-1; i++) x[k * (mx-1)*(my-1) + j*(mx-1) + i] = coor[k][j][i].x;
		I = TECDAT100(&III, x, &DIsDouble);

		// Y
		for (k=zs; k<ze-1; k++)
		for (j=ys; j<ye-1; j++)
		for (i=xs; i<xe-1; i++)	x[k * (mx-1)*(my-1) + j*(mx-1) + i] = coor[k][j][i].y;
		I = TECDAT100(&III, x, &DIsDouble);

		// Z
		for (k=zs; k<ze-1; k++)
		for (j=ys; j<ye-1; j++)
		for (i=xs; i<xe-1; i++)	x[k * (mx-1)*(my-1) + j*(mx-1) + i] = coor[k][j][i].z;
		I = TECDAT100(&III, x, &DIsDouble);
	
		DMDAVecRestoreArray(fda, Coor, &coor);
    
		//delete []x;
		double N=(double)tis+1.0;
		if(phase) N = phase;
		//x = new float [(mx-2)*(my-2)*(mz-2)];

		III = (mx-2) * (my-2) * (mz-2);

		PetscViewer viewer;
		char filen[128];
		
		DMCreateGlobalVector(user[bi].fda, &user[bi].Ucat_sum);
		DMCreateGlobalVector(user[bi].fda, &user[bi].Ucat_cross_sum);
		DMCreateGlobalVector(user[bi].fda, &user[bi].Ucat_square_sum);
		
		DMCreateGlobalVector(user[bi].fda, &Vort);
		
		DMCreateGlobalVector(user[bi].fda, &Ucat);
		DMCreateGlobalVector(user[bi].da, &UU);
		DMCreateGlobalVector(user[bi].da, &VV);
		DMCreateGlobalVector(user[bi].da, &WW);
		DMCreateGlobalVector(user[bi].da, &UV);
		DMCreateGlobalVector(user[bi].da, &VW);
		DMCreateGlobalVector(user[bi].da, &WU);

		/*
		if(phase) {
			sprintf(filen, "phase%03d_su0_%06d_%1d.dat", phase, ti_phase, user[bi]._this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoad(viewer, (user[bi].Ucat_sum));
			PetscViewerDestroy(&viewer);
				
			PetscPrintf(PETSC_COMM_WORLD, "Read %s\n", filen);

			sprintf(filen, "phase%03d_su1_%06d_%1d.dat", phase, ti_phase, user[bi]._this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoad(viewer, (user[bi].Ucat_cross_sum));
			PetscViewerDestroy(&viewer);
				
			PetscPrintf(PETSC_COMM_WORLD, "Read %s\n", filen);
			
			sprintf(filen, "phase%03d_su2_%06d_%1d.dat", phase, ti_phase, user[bi]._this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoad(viewer, (user[bi].Ucat_square_sum));
			PetscViewerDestroy(&viewer);
				
			PetscPrintf(PETSC_COMM_WORLD, "Read %s\n", filen);
		}
		else */
		{
			sprintf(filen, "%ssu0_%06d_%1d.dat", prefix, ti, user[bi]._this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoad(user[bi].Ucat_sum,viewer);
			PetscViewerDestroy(&viewer);

			sprintf(filen, "%ssu1_%06d_%1d.dat", prefix, ti, user[bi]._this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoad(user[bi].Ucat_cross_sum,viewer);
			PetscViewerDestroy(&viewer);
			
			sprintf(filen, "%ssu2_%06d_%1d.dat", prefix, ti, user[bi]._this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoad(user[bi].Ucat_square_sum,viewer);
			PetscViewerDestroy(&viewer);
			
			sprintf(filen, "%ssvo_%06d_%1d.dat", prefix, ti, user[bi]._this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoad(Vort,viewer);
			PetscViewerDestroy(&viewer);
			VecScale(Vort, 1./N);
		}
			
		DMDAVecGetArray(user[bi].fda, user[bi].Ucat_sum, &usum);
		DMDAVecGetArray(user[bi].fda, user[bi].Ucat_cross_sum, &u1sum);
		DMDAVecGetArray(user[bi].fda, user[bi].Ucat_square_sum, &u2sum);
		
		DMDAVecGetArray(user[bi].fda, Vort, &vort);
		DMDAVecGetArray(user[bi].da, user[bi].Nvert, &nvert);
		DMDAVecGetArray(user[bi].da, UU, &uu);
		DMDAVecGetArray(user[bi].da, VV, &vv);
		DMDAVecGetArray(user[bi].da, WW, &ww);
		DMDAVecGetArray(user[bi].da, UV, &uv);
		DMDAVecGetArray(user[bi].da, VW, &vw);
		DMDAVecGetArray(user[bi].da, WU, &wu);
			
		// Reynolds stresses
		for (k=zs; k<ze-2; k++)
		for (j=ys; j<ye-2; j++)
		for (i=xs; i<xe-2; i++) {
			double _U = usum[k+1][j+1][i+1].x/N;
			double _V = usum[k+1][j+1][i+1].y/N;
			double _W = usum[k+1][j+1][i+1].z/N;
			double _UV = usum[k+1][j+1][i+1].x*usum[k+1][j+1][i+1].y / (N*N);
			double _VW = usum[k+1][j+1][i+1].y*usum[k+1][j+1][i+1].z / (N*N);
			double _WU = usum[k+1][j+1][i+1].z*usum[k+1][j+1][i+1].x / (N*N);
			
			uu[k+1][j+1][i+1] = ( u2sum[k+1][j+1][i+1].x/N - _U*_U );
			vv[k+1][j+1][i+1] = ( u2sum[k+1][j+1][i+1].y/N - _V*_V );
			ww[k+1][j+1][i+1] = ( u2sum[k+1][j+1][i+1].z/N - _W*_W );
			uv[k+1][j+1][i+1] = u1sum[k+1][j+1][i+1].x/N - _UV;
			vw[k+1][j+1][i+1] = u1sum[k+1][j+1][i+1].y/N - _VW;
			wu[k+1][j+1][i+1] = u1sum[k+1][j+1][i+1].z/N - _WU;
		}
		
		DMDAVecRestoreArray(user[bi].fda, user[bi].Ucat_sum, &usum);
		DMDAVecRestoreArray(user[bi].fda, user[bi].Ucat_cross_sum, &u1sum);
		DMDAVecRestoreArray(user[bi].fda, user[bi].Ucat_square_sum, &u2sum);
					
		// Mean velocity
		VecAXPY(Ucat, 1./N, user[bi].Ucat_sum);
		
		// Destroy
		VecDestroy(&user[bi].Ucat_sum);
		VecDestroy(&user[bi].Ucat_cross_sum);
		VecDestroy(&user[bi].Ucat_square_sum);
		
		DMDAVecGetArray(user[bi].fda, Ucat, &ucat);
				
		// Write X-Vorticity
		for (k=zs; k<ze-2; k++)
		for (j=ys; j<ye-2; j++)
		for (i=xs; i<xe-2; i++) {
			x[k * (mx-2)*(my-2) + j*(mx-2) + i] =vort[k+1][j+1][i+1].x;
		}
		if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		I = TECDAT100(&III, x, &DIsDouble);
		
		/* Convection: u dot dvort */
		/* Convection-1 */
		for (k=zs; k<ze-2; k++)
		for (j=ys; j<ye-2; j++)
		for (i=xs; i<xe-2; i++) {
			double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
			double eta0= eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
			double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
			double ajc = aj[k+1][j+1][i+1];
			double dudc, dude, dudz;
			double dvdc, dvde, dvdz;
			double dwdc, dwde, dwdz;
			double du_dx, du_dy, du_dz;
			double dv_dx, dv_dy, dv_dz;
			double dw_dx, dw_dy, dw_dz;
			
			Compute_du_center (i+1, j+1, k+1, mx, my, mz, vort, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
			Compute_du_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
			
			double Convection = 0;
			
			Convection += ucat[k+1][j+1][i+1].x * du_dx;
			//Convection += ucat[k+1][j+1][i+1].y * du_dy;
			//Convection += ucat[k+1][j+1][i+1].z * du_dz;
			Convection *= -1.;
			
			x[k * (mx-2)*(my-2) + j*(mx-2) + i] = Convection;
		}
		if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		I = TECDAT100(&III, x, &DIsDouble);
		
		/* Convection-2 */
		for (k=zs; k<ze-2; k++)
		for (j=ys; j<ye-2; j++)
		for (i=xs; i<xe-2; i++) {
			double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
			double eta0= eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
			double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
			double ajc = aj[k+1][j+1][i+1];
			double dudc, dude, dudz;
			double dvdc, dvde, dvdz;
			double dwdc, dwde, dwdz;
			double du_dx, du_dy, du_dz;
			double dv_dx, dv_dy, dv_dz;
			double dw_dx, dw_dy, dw_dz;
			
			Compute_du_center (i+1, j+1, k+1, mx, my, mz, vort, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
			Compute_du_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
			
			double Convection = 0;
			
			//Convection += ucat[k+1][j+1][i+1].x * du_dx;
			Convection += ucat[k+1][j+1][i+1].y * du_dy;
			//Convection += ucat[k+1][j+1][i+1].z * du_dz;
			Convection *= -1.;
			
			x[k * (mx-2)*(my-2) + j*(mx-2) + i] = Convection;
		}
		if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		I = TECDAT100(&III, x, &DIsDouble);
		
		/* Convection-3 */
		for (k=zs; k<ze-2; k++)
		for (j=ys; j<ye-2; j++)
		for (i=xs; i<xe-2; i++) {
			double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
			double eta0= eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
			double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
			double ajc = aj[k+1][j+1][i+1];
			double dudc, dude, dudz;
			double dvdc, dvde, dvdz;
			double dwdc, dwde, dwdz;
			double du_dx, du_dy, du_dz;
			double dv_dx, dv_dy, dv_dz;
			double dw_dx, dw_dy, dw_dz;
			
			Compute_du_center (i+1, j+1, k+1, mx, my, mz, vort, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
			Compute_du_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
			
			double Convection = 0;
			
			//Convection += ucat[k+1][j+1][i+1].x * du_dx;
			//Convection += ucat[k+1][j+1][i+1].y * du_dy;
			Convection += ucat[k+1][j+1][i+1].z * du_dz;
			Convection *= -1.;
			
			x[k * (mx-2)*(my-2) + j*(mx-2) + i] = Convection;
		}
		if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		I = TECDAT100(&III, x, &DIsDouble);
		
		/* Diffusion */
		// make a vector for the gradient of the X-vorticity
		Vec Grad_Vortx;
		Cmpnts ***grad_vortx;
		
		DMCreateGlobalVector(user[bi].fda, &Grad_Vortx);
		DMDAVecGetArray(user[bi].fda, Grad_Vortx, &grad_vortx);
		
		for (k=zs; k<ze-2; k++)
		for (j=ys; j<ye-2; j++)
		for (i=xs; i<xe-2; i++) {
			double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
			double eta0= eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
			double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
			double ajc = aj[k+1][j+1][i+1];
			double dudc, dude, dudz;
			double dvdc, dvde, dvdz;
			double dwdc, dwde, dwdz;
			double du_dx, du_dy, du_dz;
			double dv_dx, dv_dy, dv_dz;
			double dw_dx, dw_dy, dw_dz;
			
			Compute_du_center (i+1, j+1, k+1, mx, my, mz, vort, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
			Compute_du_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
			
			grad_vortx[k+1][j+1][i+1].x = du_dx;
			grad_vortx[k+1][j+1][i+1].y = du_dy;
			grad_vortx[k+1][j+1][i+1].z = du_dz;
		}
		
		// diffusion
		for (k=zs; k<ze-2; k++)
		for (j=ys; j<ye-2; j++)
		for (i=xs; i<xe-2; i++) {
			double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
			double eta0= eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
			double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
			double ajc = aj[k+1][j+1][i+1];
			double dudc, dude, dudz;
			double dvdc, dvde, dvdz;
			double dwdc, dwde, dwdz;
			double du_dx, du_dy, du_dz;
			double dv_dx, dv_dy, dv_dz;
			double dw_dx, dw_dy, dw_dz;
			
			Compute_du_center (i+1, j+1, k+1, mx, my, mz, grad_vortx, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
			Compute_du_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
			
			double Diffusion = 0;
			
			Diffusion += du_dx;	// Laplacian = Div.Grad
			Diffusion += dv_dy;
			Diffusion += dw_dz;
			Diffusion *= 1./user[bi].ren;
					
			x[k * (mx-2)*(my-2) + j*(mx-2) + i] = Diffusion;
		}
		if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		I = TECDAT100(&III, x, &DIsDouble);
		
		DMDAVecRestoreArray(user[bi].fda, Grad_Vortx, &grad_vortx);
		VecDestroy(&Grad_Vortx);
		
		
		/* Vortex Stretching*/
		for (k=zs; k<ze-2; k++)
		for (j=ys; j<ye-2; j++)
		for (i=xs; i<xe-2; i++) {
			double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
			double eta0= eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
			double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
			double ajc = aj[k+1][j+1][i+1];
			double dudc, dude, dudz;
			double dvdc, dvde, dvdz;
			double dwdc, dwde, dwdz;
			double du_dx, du_dy, du_dz;
			double dv_dx, dv_dy, dv_dz;
			double dw_dx, dw_dy, dw_dz;
			
			Compute_du_center (i+1, j+1, k+1, mx, my, mz, ucat, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
			Compute_du_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
			
			double Stretching = vort[k+1][j+1][i+1].x * du_dx;
			x[k * (mx-2)*(my-2) + j*(mx-2) + i] = Stretching;
		}
		if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		I = TECDAT100(&III, x, &DIsDouble);
		
		/* Vortex Skewing*/
		/* Vortex Skewing-y*/
		for (k=zs; k<ze-2; k++)
		for (j=ys; j<ye-2; j++)
		for (i=xs; i<xe-2; i++) {
			double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
			double eta0= eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
			double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
			double ajc = aj[k+1][j+1][i+1];
			double dudc, dude, dudz;
			double dvdc, dvde, dvdz;
			double dwdc, dwde, dwdz;
			double du_dx, du_dy, du_dz;
			double dv_dx, dv_dy, dv_dz;
			double dw_dx, dw_dy, dw_dz;
			
			Compute_du_center (i+1, j+1, k+1, mx, my, mz, ucat, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
			Compute_du_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
			
			double Skewing = 0;
			Skewing += vort[k+1][j+1][i+1].y * du_dy;
			//Skewing += vort[k+1][j+1][i+1].z * du_dz;
			x[k * (mx-2)*(my-2) + j*(mx-2) + i] = Skewing;
		}
		if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		I = TECDAT100(&III, x, &DIsDouble);
		
		/* Vortex Skewing-z*/
		for (k=zs; k<ze-2; k++)
		for (j=ys; j<ye-2; j++)
		for (i=xs; i<xe-2; i++) {
			double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
			double eta0= eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
			double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
			double ajc = aj[k+1][j+1][i+1];
			double dudc, dude, dudz;
			double dvdc, dvde, dvdz;
			double dwdc, dwde, dwdz;
			double du_dx, du_dy, du_dz;
			double dv_dx, dv_dy, dv_dz;
			double dw_dx, dw_dy, dw_dz;
			
			Compute_du_center (i+1, j+1, k+1, mx, my, mz, ucat, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
			Compute_du_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
			
			double Skewing = 0;
			//Skewing += vort[k+1][j+1][i+1].y * du_dy;
			Skewing += vort[k+1][j+1][i+1].z * du_dz;
			x[k * (mx-2)*(my-2) + j*(mx-2) + i] = Skewing;
		}
		if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		I = TECDAT100(&III, x, &DIsDouble);
		
		/* P2 term of H. J. Perkins */
		// make a vector for the gradient
		Vec Grad2;
		PetscReal ***grad2;
		
		DMCreateGlobalVector(user[bi].da, &Grad2);
		DMDAVecGetArray(user[bi].da, Grad2, &grad2);
		
		for (k=zs; k<ze-2; k++)
		for (j=ys; j<ye-2; j++)
		for (i=xs; i<xe-2; i++) {
			double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
			double eta0= eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
			double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
			double ajc = aj[k+1][j+1][i+1];
			double duvdc, duvde, duvdz;
			double dwudc, dwude, dwudz;
			double duv_dx, duv_dy, duv_dz;
			double dwu_dx, dwu_dy, dwu_dz;
			
			Compute_dscalar_center (i+1, j+1, k+1, mx, my, mz, uv, nvert, &duvdc, &duvde, &duvdz);
			Compute_dscalar_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, duvdc, duvde, duvdz, &duv_dx, &duv_dy, &duv_dz);
			
			Compute_dscalar_center (i+1, j+1, k+1, mx, my, mz, wu, nvert, &dwudc, &dwude, &dwudz);
			Compute_dscalar_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dwudc, dwude, dwudz, &dwu_dx, &dwu_dy, &dwu_dz);
			
			grad2[k+1][j+1][i+1] = duv_dz - dwu_dy;
		}
		
		// P2
		for (k=zs; k<ze-2; k++)
		for (j=ys; j<ye-2; j++)
		for (i=xs; i<xe-2; i++) {
			double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
			double eta0= eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
			double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
			double ajc = aj[k+1][j+1][i+1];
			double dpdc, dpde, dpdz;
			double dp_dx, dp_dy, dp_dz;
			
			Compute_dscalar_center (i+1, j+1, k+1, mx, my, mz, grad2, nvert, &dpdc, &dpde, &dpdz);
			Compute_dscalar_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dpdc, dpde, dpdz, &dp_dx, &dp_dy, &dp_dz);
					
			x[k * (mx-2)*(my-2) + j*(mx-2) + i] = dp_dx;
		}
		if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		I = TECDAT100(&III, x, &DIsDouble);
		
		DMDAVecRestoreArray(user[bi].da, Grad2, &grad2);
		VecDestroy(&Grad2);
		
		/* P3 term of H. J. Perkins */
		// make a vector for the gradient
		Vec Grad3;
		PetscReal ***grad3;
		
		DMCreateGlobalVector(user[bi].da, &Grad3);
		DMDAVecGetArray(user[bi].da, Grad3, &grad3);
		
		for (k=zs; k<ze-2; k++)
		for (j=ys; j<ye-2; j++)
		for (i=xs; i<xe-2; i++) {
			double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
			double eta0= eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
			double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
			double ajc = aj[k+1][j+1][i+1];
			double dvvdc, dvvde, dvvdz;
			double dwwdc, dwwde, dwwdz;
			double dvv_dx, dvv_dy, dvv_dz;
			double dww_dx, dww_dy, dww_dz;
			
			Compute_dscalar_center (i+1, j+1, k+1, mx, my, mz, vv, nvert, &dvvdc, &dvvde, &dvvdz);
			Compute_dscalar_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dvvdc, dvvde, dvvdz, &dvv_dx, &dvv_dy, &dvv_dz);
			
			Compute_dscalar_center (i+1, j+1, k+1, mx, my, mz, ww, nvert, &dwwdc, &dwwde, &dwwdz);
			Compute_dscalar_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dwwdc, dwwde, dwwdz, &dww_dx, &dww_dy, &dww_dz);
			
			grad3[k+1][j+1][i+1] = dvv_dz - dww_dz;
		}
		
		// P3
		for (k=zs; k<ze-2; k++)
		for (j=ys; j<ye-2; j++)
		for (i=xs; i<xe-2; i++) {
			double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
			double eta0= eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
			double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
			double ajc = aj[k+1][j+1][i+1];
			double dpdc, dpde, dpdz;
			double dp_dx, dp_dy, dp_dz;
			
			Compute_dscalar_center (i+1, j+1, k+1, mx, my, mz, grad3, nvert, &dpdc, &dpde, &dpdz);
			Compute_dscalar_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dpdc, dpde, dpdz, &dp_dx, &dp_dy, &dp_dz);
					
			x[k * (mx-2)*(my-2) + j*(mx-2) + i] = dp_dy;
		}
		if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		I = TECDAT100(&III, x, &DIsDouble);
		
		DMDAVecRestoreArray(user[bi].da, Grad3, &grad3);
		VecDestroy(&Grad3);
		
		/* P4 term of H. J. Perkins */
		// make a vector for the gradient
		Vec Grad4;
		Cmpnts ***grad4;
		
		DMCreateGlobalVector(user[bi].fda, &Grad4);
		DMDAVecGetArray(user[bi].fda, Grad4, &grad4);
		
		for (k=zs; k<ze-2; k++)
		for (j=ys; j<ye-2; j++)
		for (i=xs; i<xe-2; i++) {
			double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
			double eta0= eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
			double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
			double ajc = aj[k+1][j+1][i+1];
			double dvwdc, dvwde, dvwdz;
			double dvw_dx, dvw_dy, dvw_dz;
			
			Compute_dscalar_center (i+1, j+1, k+1, mx, my, mz, vw, nvert, &dvwdc, &dvwde, &dvwdz);
			Compute_dscalar_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dvwdc, dvwde, dvwdz, &dvw_dx, &dvw_dy, &dvw_dz);
			
			grad4[k+1][j+1][i+1].x = dvw_dx;
			grad4[k+1][j+1][i+1].y = dvw_dy;
			grad4[k+1][j+1][i+1].z = dvw_dz;
		}
		
		// P4
		for (k=zs; k<ze-2; k++)
		for (j=ys; j<ye-2; j++)
		for (i=xs; i<xe-2; i++) {
			double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
			double eta0= eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
			double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
			double ajc = aj[k+1][j+1][i+1];
			double dudc, dude, dudz;
			double dvdc, dvde, dvdz;
			double dwdc, dwde, dwdz;
			double du_dx, du_dy, du_dz;
			double dv_dx, dv_dy, dv_dz;
			double dw_dx, dw_dy, dw_dz;
			
			Compute_du_center (i+1, j+1, k+1, mx, my, mz, grad4, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
			Compute_du_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
			
			x[k * (mx-2)*(my-2) + j*(mx-2) + i] = dw_dz - dv_dy;
		}
		if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		I = TECDAT100(&III, x, &DIsDouble);
		
		DMDAVecRestoreArray(user[bi].fda, Grad4, &grad4);
		VecDestroy(&Grad4);
		/*
		// Write Nv
		
		for (k=zs; k<ze-2; k++)
		for (j=ys; j<ye-2; j++)
		for (i=xs; i<xe-2; i++) x[k * (mx-2)*(my-2) + j*(mx-2) + i] = nvert[k+1][j+1][i+1];
		I = TECDAT100(&III, x, &DIsDouble);
		
		*/
		delete []x;
		
		DMDAVecRestoreArray(user[bi].fda, Ucat, &ucat);
		DMDAVecRestoreArray(user[bi].da, user[bi].Nvert, &nvert);
		DMDAVecRestoreArray(user[bi].fda, Vort, &vort);
		DMDAVecRestoreArray(user[bi].da, UU, &uu);
		DMDAVecRestoreArray(user[bi].da, VV, &vv);
		DMDAVecRestoreArray(user[bi].da, WW, &ww);
		DMDAVecRestoreArray(user[bi].da, UV, &uv);
		DMDAVecRestoreArray(user[bi].da, VW, &vw);
		DMDAVecRestoreArray(user[bi].da, WU, &wu);
				
		DMDAVecRestoreArray(user[bi].fda, user[bi].Csi, &csi);
		DMDAVecRestoreArray(user[bi].fda, user[bi].Eta, &eta);
		DMDAVecRestoreArray(user[bi].fda, user[bi].Zet, &zet);
		DMDAVecRestoreArray(user[bi].da, user[bi].Aj, &aj);
		
		VecDestroy(&Vort);
		VecDestroy(&Ucat);
		VecDestroy(&UU);
		VecDestroy(&VV);
		VecDestroy(&WW);
		VecDestroy(&UV);
		VecDestroy(&VW);
		VecDestroy(&WU);
	}
	
	I = TECEND100();
	
	return 0;
}



PetscErrorCode TECIOOutQ(UserCtx *user, int Q)
{
	PetscInt bi;

	char filen[80];
	sprintf(filen, "%sQCriteria%06d.plt", prefix, ti);

	INTEGER4 I, Debug, VIsDouble, DIsDouble, III, IMax, JMax, KMax;
	VIsDouble = 0;
	DIsDouble = 0;
	Debug = 0;

	if(Q==1) {
		printf("qcr=%d, Q-Criterion\n", Q);
		//I = TECINI100((char*)"Result", (char*)"X Y Z Q Velocity_Magnitude Nv",   filen,  (char*)".",   &Debug,  &VIsDouble);
		I = TECINI100((char*)"Result", (char*)"X Y Z Q Velocity_Magnitude",   filen,  (char*)".",   &Debug,  &VIsDouble);
	}
	else if(Q==2) {
		printf("Lambda2-Criterion\n");
		//I = TECINI100((char*)"Result", (char*)"X Y Z Lambda2 Velocity_Magnitude Nv",   filen,  (char*)".",   &Debug,  &VIsDouble);
		I = TECINI100((char*)"Result", (char*)"X Y Z Lambda2 Velocity_Magnitude",   filen,  (char*)".",   &Debug,  &VIsDouble);
	}
	else if(Q==3) {
		printf("Q-Criterion from saved file\n");
		I = TECINI100((char*)"Result", (char*)"X Y Z Q Velocity_Magnitude",   filen,  (char*)".",   &Debug,  &VIsDouble);
	}

	for (bi=0; bi<block_number; bi++) {
		DM da = user[bi].da, fda = user[bi].fda;
		DMDALocalInfo info = user[bi].info;

		PetscInt	xs = info.xs, xe = info.xs + info.xm;
		PetscInt  	ys = info.ys, ye = info.ys + info.ym;
		PetscInt	zs = info.zs, ze = info.zs + info.zm;
		PetscInt	mx = info.mx, my = info.my, mz = info.mz;
    
		PetscInt	lxs, lys, lzs, lxe, lye, lze;
		PetscInt	i, j, k;
		PetscReal	***aj;
		Cmpnts	***ucat, ***coor, ***ucat_o;
		PetscReal	***p, ***nvert;
		Vec		Coor, zCoor, nCoor;
		VecScatter	ctx;

		Vec X, Y, Z, U, V, W;

		INTEGER4	ZoneType=0, ICellMax=0, JCellMax=0, KCellMax=0;
		INTEGER4	IsBlock=1, NumFaceConnections=0, FaceNeighborMode=0;
		INTEGER4    ShareConnectivityFromZone=0;
		INTEGER4	LOC[8] = {1, 1, 1, 0, 0, 0, 0}; /* 1 is cell-centered 0 is node centered */

		IMax = mx-1; JMax = my-1; KMax = mz-1;

		I = TECZNE100((char*)"Block 1",
			&ZoneType, 	/* Ordered zone */
			&IMax,
			&JMax,
			&KMax,
			&ICellMax,
			&JCellMax,
			&KCellMax,
			&IsBlock,	/* ISBLOCK  BLOCK format */
			&NumFaceConnections,
			&FaceNeighborMode,
			LOC,
			NULL,
			&ShareConnectivityFromZone); /* No connectivity sharing */

		III = (mx-1) * (my-1) * (mz-1);
		float	*x;
		PetscMalloc(mx*my*mz*sizeof(float), &x);	// seokkoo
	
		DMDAGetCoordinates(da, &Coor);
		DMDAVecGetArray(fda, Coor, &coor);

		for (k=zs; k<ze-1; k++)
		for (j=ys; j<ye-1; j++)
		for (i=xs; i<xe-1; i++) {
			x[k * (mx-1)*(my-1) + j*(mx-1) + i] = coor[k][j][i].x;
		}
		I = TECDAT100(&III, x, &DIsDouble);

		for (k=zs; k<ze-1; k++)
		for (j=ys; j<ye-1; j++)
		for (i=xs; i<xe-1; i++) {
			x[k * (mx-1)*(my-1) + j*(mx-1) + i] = coor[k][j][i].y;
		}
		I = TECDAT100(&III, x, &DIsDouble);

		for (k=zs; k<ze-1; k++)
		for (j=ys; j<ye-1; j++)
		for (i=xs; i<xe-1; i++) {
			x[k * (mx-1)*(my-1) + j*(mx-1) + i] = coor[k][j][i].z;
		}
      
		I = TECDAT100(&III, x, &DIsDouble);
		DMDAVecRestoreArray(fda, Coor, &coor);
    
		III = (mx-2) * (my-2) * (mz-2);

		if(Q==1) {
			QCriteria(user);
			DMDAVecGetArray(user[bi].da, user[bi].P, &p);
			for (k=zs; k<ze-2; k++)
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) {
				x[k * (mx-2)*(my-2) + j*(mx-2) + i] =p[k+1][j+1][i+1];
			}
			I = TECDAT100(&III, x, &DIsDouble);
			DMDAVecRestoreArray(user[bi].da, user[bi].P, &p);
		}
		else if(Q==2) {
			Lambda2(user);
			DMDAVecGetArray(user[bi].da, user[bi].P, &p);
			for (k=zs; k<ze-2; k++)
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) {
				x[k * (mx-2)*(my-2) + j*(mx-2) + i] = p[k+1][j+1][i+1];
			}
			I = TECDAT100(&III, x, &DIsDouble);
			DMDAVecRestoreArray(user[bi].da, user[bi].P, &p);
		}
		else if(Q==3) {
			char filen2[128];
			PetscViewer	viewer;
			
			sprintf(filen2, "%sqfield%06d_%1d.dat", prefix, ti, user->_this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen2, FILE_MODE_READ, &viewer);
			VecLoad(user[bi].P,viewer);
			PetscViewerDestroy(&viewer);  
			
			DMDAVecGetArray(user[bi].da, user[bi].P, &p);
			for (k=zs; k<ze-2; k++)
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) {
				x[k * (mx-2)*(my-2) + j*(mx-2) + i] = p[k+1][j+1][i+1];
			}
			I = TECDAT100(&III, x, &DIsDouble);
			DMDAVecRestoreArray(user[bi].da, user[bi].P, &p);
		}
		
		Velocity_Magnitude(user);
		DMDAVecGetArray(user[bi].da, user[bi].P, &p);
		for (k=zs; k<ze-2; k++)
		for (j=ys; j<ye-2; j++)
		for (i=xs; i<xe-2; i++) {
			x[k * (mx-2)*(my-2) + j*(mx-2) + i] = p[k+1][j+1][i+1];
		}
		I = TECDAT100(&III, x, &DIsDouble);
		DMDAVecRestoreArray(user[bi].da, user[bi].P, &p);
    
		/*
		DMDAVecGetArray(user[bi].da, user[bi].Nvert, &nvert);
		for (k=zs; k<ze-2; k++)
		for (j=ys; j<ye-2; j++)
		for (i=xs; i<xe-2; i++) {
			x[k * (mx-2)*(my-2) + j*(mx-2) + i] = nvert[k+1][j+1][i+1];
		}
		I = TECDAT100(&III, x, &DIsDouble);
		DMDAVecRestoreArray(user[bi].da, user[bi].Nvert, &nvert);
		*/
		
		PetscFree(x);
	}
	I = TECEND100();
	
	return 0;
}

PetscErrorCode FormMetrics(UserCtx *user)
{
  DM		cda;
  Cmpnts	***csi, ***eta, ***zet;
  PetscScalar	***aj;
  Vec		coords;
  Cmpnts	***coor;

  DM		da = user->da, fda = user->fda;
  Vec		Csi = user->Csi, Eta = user->Eta, Zet = user->Zet;
  Vec		Aj = user->Aj;
  Vec		ICsi = user->ICsi, IEta = user->IEta, IZet = user->IZet;
  Vec		JCsi = user->JCsi, JEta = user->JEta, JZet = user->JZet;
  Vec		KCsi = user->KCsi, KEta = user->KEta, KZet = user->KZet;
  Vec		IAj = user->IAj, JAj = user->JAj, KAj = user->KAj;

  
  Cmpnts	***icsi, ***ieta, ***izet;
  Cmpnts	***jcsi, ***jeta, ***jzet;
  Cmpnts	***kcsi, ***keta, ***kzet;
  Cmpnts	***gs;
  PetscReal	***iaj, ***jaj, ***kaj;

  Vec		Cent = user->Cent; //local working array for storing cell center geometry

  Vec		Centx, Centy, Centz, lCoor;
  Cmpnts	***cent, ***centx, ***centy, ***centz;

  PetscInt	xs, ys, zs, xe, ye, ze;
  DMDALocalInfo	info;

  PetscInt	mx, my, mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  PetscScalar	dxdc, dydc, dzdc, dxde, dyde, dzde, dxdz, dydz, dzdz;
  PetscInt	i, j, k, ia, ja, ka, ib, jb, kb;
  PetscInt	gxs, gxe, gys, gye, gzs, gze;
  PetscErrorCode	ierr;

  PetscReal	xcp, ycp, zcp, xcm, ycm, zcm;
  DMDAGetLocalInfo(da, &info);
  mx = info.mx; my = info.my; mz = info.mz;
  xs = info.xs; xe = xs + info.xm;
  ys = info.ys; ye = ys + info.ym;
  zs = info.zs; ze = zs + info.zm;

  gxs = info.gxs; gxe = gxs + info.gxm;
  gys = info.gys; gye = gys + info.gym;
  gzs = info.gzs; gze = gzs + info.gzm;

  DMDAGetCoordinateDA(da, &cda);
  DMDAVecGetArray(cda, Csi, &csi);
  DMDAVecGetArray(cda, Eta, &eta);
  DMDAVecGetArray(cda, Zet, &zet);
  ierr = DMDAVecGetArray(da, Aj,  &aj); CHKERRQ(ierr);

  DMDAGetGhostedCoordinates(da, &coords);
  DMDAVecGetArray(fda, coords, &coor);


  //  VecDuplicate(coords, &Cent);
  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  /* Calculating transformation metrics in i direction */
  for (k=lzs; k<lze; k++){
    for (j=lys; j<lye; j++){
      for (i=xs; i<lxe; i++){
	/* csi = de X dz */
	dxde = 0.5 * (coor[k  ][j  ][i  ].x + coor[k-1][j  ][i  ].x -
		      coor[k  ][j-1][i  ].x - coor[k-1][j-1][i  ].x);
	dyde = 0.5 * (coor[k  ][j  ][i  ].y + coor[k-1][j  ][i  ].y -
		      coor[k  ][j-1][i  ].y - coor[k-1][j-1][i  ].y);
	dzde = 0.5 * (coor[k  ][j  ][i  ].z + coor[k-1][j  ][i  ].z -
		      coor[k  ][j-1][i  ].z - coor[k-1][j-1][i  ].z);
				       		    	    	 
	dxdz = 0.5 * (coor[k  ][j-1][i  ].x + coor[k  ][j  ][i  ].x -
		      coor[k-1][j-1][i  ].x - coor[k-1][j  ][i  ].x);
	dydz = 0.5 * (coor[k  ][j-1][i  ].y + coor[k  ][j  ][i  ].y -
		      coor[k-1][j-1][i  ].y - coor[k-1][j  ][i  ].y);
	dzdz = 0.5 * (coor[k  ][j-1][i  ].z + coor[k  ][j  ][i  ].z -
		      coor[k-1][j-1][i  ].z - coor[k-1][j  ][i  ].z);
	  
	csi[k][j][i].x = dyde * dzdz - dzde * dydz;
	csi[k][j][i].y =-dxde * dzdz + dzde * dxdz;
	csi[k][j][i].z = dxde * dydz - dyde * dxdz;

	
      }
    }
  }

  // Need more work -- lg65
  /* calculating j direction metrics */
  for (k=lzs; k<lze; k++){
    for (j=ys; j<lye; j++){
      for (i=lxs; i<lxe; i++){

	/* eta = dz X de */
	dxdc = 0.5 * (coor[k  ][j  ][i  ].x + coor[k-1][j  ][i  ].x -
		      coor[k  ][j  ][i-1].x - coor[k-1][j  ][i-1].x);
	dydc = 0.5 * (coor[k  ][j  ][i  ].y + coor[k-1][j  ][i  ].y -
		      coor[k  ][j  ][i-1].y - coor[k-1][j  ][i-1].y);
	dzdc = 0.5 * (coor[k  ][j  ][i  ].z + coor[k-1][j  ][i  ].z -
		      coor[k  ][j  ][i-1].z - coor[k-1][j  ][i-1].z);
			    		         	 		   	 
	dxdz = 0.5 * (coor[k  ][j  ][i  ].x + coor[k  ][j  ][i-1].x -
		      coor[k-1][j  ][i  ].x - coor[k-1][j  ][i-1].x);
	dydz = 0.5 * (coor[k  ][j  ][i  ].y + coor[k  ][j  ][i-1].y -
		      coor[k-1][j  ][i  ].y - coor[k-1][j  ][i-1].y);
	dzdz = 0.5 * (coor[k  ][j  ][i  ].z + coor[k  ][j  ][i-1].z -
		      coor[k-1][j  ][i  ].z - coor[k-1][j  ][i-1].z);
	  
	eta[k][j][i].x = dydz * dzdc - dzdz * dydc;
	eta[k][j][i].y =-dxdz * dzdc + dzdz * dxdc;
	eta[k][j][i].z = dxdz * dydc - dydz * dxdc;

      }
    }
  }

  /* calculating k direction metrics */
  for (k=zs; k<lze; k++){
    for (j=lys; j<lye; j++){
      for (i=lxs; i<lxe; i++){
	dxdc = 0.5 * (coor[k  ][j  ][i  ].x + coor[k  ][j-1][i  ].x -
		      coor[k  ][j  ][i-1].x - coor[k  ][j-1][i-1].x);
	dydc = 0.5 * (coor[k  ][j  ][i  ].y + coor[k  ][j-1][i  ].y -
		      coor[k  ][j  ][i-1].y - coor[k  ][j-1][i-1].y);
	dzdc = 0.5 * (coor[k  ][j  ][i  ].z + coor[k  ][j-1][i  ].z -
		      coor[k  ][j  ][i-1].z - coor[k  ][j-1][i-1].z);
			    		    	     	     	 
	dxde = 0.5 * (coor[k  ][j  ][i  ].x + coor[k  ][j  ][i-1].x -
		      coor[k  ][j-1][i  ].x - coor[k  ][j-1][i-1].x);
	dyde = 0.5 * (coor[k  ][j  ][i  ].y + coor[k  ][j  ][i-1].y -
		      coor[k  ][j-1][i  ].y - coor[k  ][j-1][i-1].y);
	dzde = 0.5 * (coor[k  ][j  ][i  ].z + coor[k  ][j  ][i-1].z -
		      coor[k  ][j-1][i  ].z - coor[k  ][j-1][i-1].z);
	  
	zet[k][j][i].x = dydc * dzde - dzdc * dyde;
	zet[k][j][i].y =-dxdc * dzde + dzdc * dxde;
	zet[k][j][i].z = dxdc * dyde - dydc * dxde;

	/*	if ((i==1 || i==mx-2) && j==1 && (k==1 || k==0)) {
	  PetscPrintf(PETSC_COMM_WORLD, "%e %e %e\n", dxdc * dyde, dydc * dxde, dzdc);
	  PetscPrintf(PETSC_COMM_WORLD, "%e %e %e\n", dxde, dyde, dzde);
	  PetscPrintf(PETSC_COMM_WORLD, "Met %e %e %e\n", zet[k][j][i].x, zet[k][j][i].y, zet[k][j][i].z);
	  }*/
	
      }
    }
  }

  /* calculating Jacobian of the transformation */
  for (k=lzs; k<lze; k++){
    for (j=lys; j<lye; j++){
      for (i=lxs; i<lxe; i++){

	if (i>0 && j>0 && k>0) {
	  dxdc = 0.25 * (coor[k  ][j  ][i  ].x + coor[k  ][j-1][i  ].x +
			 coor[k-1][j  ][i  ].x + coor[k-1][j-1][i  ].x -
			 coor[k  ][j  ][i-1].x - coor[k  ][j-1][i-1].x -
			 coor[k-1][j  ][i-1].x - coor[k-1][j-1][i-1].x);
	  dydc = 0.25 * (coor[k  ][j  ][i  ].y + coor[k  ][j-1][i  ].y +
			 coor[k-1][j  ][i  ].y + coor[k-1][j-1][i  ].y -
			 coor[k  ][j  ][i-1].y - coor[k  ][j-1][i-1].y -
			 coor[k-1][j  ][i-1].y - coor[k-1][j-1][i-1].y);
	  dzdc = 0.25 * (coor[k  ][j  ][i  ].z + coor[k  ][j-1][i  ].z +
			 coor[k-1][j  ][i  ].z + coor[k-1][j-1][i  ].z -
			 coor[k  ][j  ][i-1].z - coor[k  ][j-1][i-1].z -
			 coor[k-1][j  ][i-1].z - coor[k-1][j-1][i-1].z);

	  dxde = 0.25 * (coor[k  ][j  ][i  ].x + coor[k  ][j  ][i-1].x +
			 coor[k-1][j  ][i  ].x + coor[k-1][j  ][i-1].x - 
			 coor[k  ][j-1][i  ].x - coor[k  ][j-1][i-1].x -
			 coor[k-1][j-1][i  ].x - coor[k-1][j-1][i-1].x);
	  dyde = 0.25 * (coor[k  ][j  ][i  ].y + coor[k  ][j  ][i-1].y +
			 coor[k-1][j  ][i  ].y + coor[k-1][j  ][i-1].y - 
			 coor[k  ][j-1][i  ].y - coor[k  ][j-1][i-1].y -
			 coor[k-1][j-1][i  ].y - coor[k-1][j-1][i-1].y);
	  dzde = 0.25 * (coor[k  ][j  ][i  ].z + coor[k  ][j  ][i-1].z +
			 coor[k-1][j  ][i  ].z + coor[k-1][j  ][i-1].z - 
			 coor[k  ][j-1][i  ].z - coor[k  ][j-1][i-1].z -
			 coor[k-1][j-1][i  ].z - coor[k-1][j-1][i-1].z);

	  dxdz = 0.25 * (coor[k  ][j  ][i  ].x + coor[k  ][j-1][i  ].x +
			 coor[k  ][j  ][i-1].x + coor[k  ][j-1][i-1].x -
			 coor[k-1][j  ][i  ].x - coor[k-1][j-1][i  ].x -
			 coor[k-1][j  ][i-1].x - coor[k-1][j-1][i-1].x);
	  dydz = 0.25 * (coor[k  ][j  ][i  ].y + coor[k  ][j-1][i  ].y +
			 coor[k  ][j  ][i-1].y + coor[k  ][j-1][i-1].y -
			 coor[k-1][j  ][i  ].y - coor[k-1][j-1][i  ].y -
			 coor[k-1][j  ][i-1].y - coor[k-1][j-1][i-1].y);
	  dzdz = 0.25 * (coor[k  ][j  ][i  ].z + coor[k  ][j-1][i  ].z +
			 coor[k  ][j  ][i-1].z + coor[k  ][j-1][i-1].z -
			 coor[k-1][j  ][i  ].z - coor[k-1][j-1][i  ].z -
			 coor[k-1][j  ][i-1].z - coor[k-1][j-1][i-1].z);
		
	  aj[k][j][i] = dxdc * (dyde * dzdz - dzde * dydz) -
	    dydc * (dxde * dzdz - dzde * dxdz) +
	    dzdc * (dxde * dydz - dyde * dxdz);
	  aj[k][j][i] = 1./aj[k][j][i];
	  
		#ifdef NEWMETRIC
		csi[k][j][i].x = dyde * dzdz - dzde * dydz;
		csi[k][j][i].y =-dxde * dzdz + dzde * dxdz;
		csi[k][j][i].z = dxde * dydz - dyde * dxdz;
		
		eta[k][j][i].x = dydz * dzdc - dzdz * dydc;
		eta[k][j][i].y =-dxdz * dzdc + dzdz * dxdc;
		eta[k][j][i].z = dxdz * dydc - dydz * dxdc;
	  
		zet[k][j][i].x = dydc * dzde - dzdc * dyde;
		zet[k][j][i].y =-dxdc * dzde + dzdc * dxde;
		zet[k][j][i].z = dxdc * dyde - dydc * dxde;
		#endif
	}
      }
    }
  }

  // mirror grid outside the boundary
if (xs==0) {
		i = xs;
		for (k=zs; k<ze; k++) 
		for (j=ys; j<ye; j++) {
			#ifdef NEWMETRIC
			csi[k][j][i] = csi[k][j][i+1];
			#endif
			eta[k][j][i] = eta[k][j][i+1];
			zet[k][j][i] = zet[k][j][i+1];
			aj[k][j][i] = aj[k][j][i+1];
		}
	}

	if (xe==mx) {
		i = xe-1;
		for (k=zs; k<ze; k++)
		for (j=ys; j<ye; j++) {
			#ifdef NEWMETRIC
			csi[k][j][i] = csi[k][j][i-1];
			#endif
			eta[k][j][i] = eta[k][j][i-1];
			zet[k][j][i] = zet[k][j][i-1];
			aj[k][j][i] = aj[k][j][i-1];
		}
	}
  

	if (ys==0) {
		j = ys;
		for (k=zs; k<ze; k++)
		for (i=xs; i<xe; i++) {
			#ifdef NEWMETRIC
			eta[k][j][i] = eta[k][j+1][i];
			#endif
			csi[k][j][i] = csi[k][j+1][i];
			zet[k][j][i] = zet[k][j+1][i];
			aj[k][j][i] = aj[k][j+1][i];
		}
	}
  

	if (ye==my) {
		j = ye-1;
		for (k=zs; k<ze; k++)
		for (i=xs; i<xe; i++) {
			#ifdef NEWMETRIC
			eta[k][j][i] = eta[k][j-1][i];
			#endif
			csi[k][j][i] = csi[k][j-1][i];
			zet[k][j][i] = zet[k][j-1][i];
			aj[k][j][i] = aj[k][j-1][i];
		}
	}
  

	if (zs==0) {
		k = zs;
		for (j=ys; j<ye; j++)
		for (i=xs; i<xe; i++) {
			#ifdef NEWMETRIC
			zet[k][j][i] = zet[k+1][j][i];
			#endif
			eta[k][j][i] = eta[k+1][j][i];
			csi[k][j][i] = csi[k+1][j][i];
			aj[k][j][i] = aj[k+1][j][i];
		}
	}
	

	if (ze==mz) {
		k = ze-1;
		for (j=ys; j<ye; j++)
		for (i=xs; i<xe; i++) {
			#ifdef NEWMETRIC
			zet[k][j][i] = zet[k-1][j][i];
			#endif
			eta[k][j][i] = eta[k-1][j][i];
			csi[k][j][i] = csi[k-1][j][i];
			aj[k][j][i] = aj[k-1][j][i];
		}
	}


  //  PetscPrintf(PETSC_COMM_WORLD, "Local info: %d", info.mx);



  DMDAVecRestoreArray(cda, Csi, &csi);
  DMDAVecRestoreArray(cda, Eta, &eta);
  DMDAVecRestoreArray(cda, Zet, &zet);
  DMDAVecRestoreArray(da, Aj,  &aj);


  DMDAVecRestoreArray(cda, coords, &coor);


  VecAssemblyBegin(Csi);
  VecAssemblyEnd(Csi);
  VecAssemblyBegin(Eta);
  VecAssemblyEnd(Eta);
  VecAssemblyBegin(Zet);
  VecAssemblyEnd(Zet);
  VecAssemblyBegin(Aj);
  VecAssemblyEnd(Aj);

  PetscBarrier(PETSC_NULL);
  return 0;
}

PetscErrorCode Ucont_P_Binary_Input(UserCtx *user)
{
  PetscViewer	viewer;

  char filen2[128];

  PetscOptionsClearValue("-vecload_block_size");
  sprintf(filen2, "%spfield%06d_%1d.dat", prefix, ti, user->_this);
  
  PetscViewer	pviewer;
  //Vec temp;
  int rank;
  PetscReal norm;
	
  if(file_exist(filen2))
  if(!onlyV) {
    //DACreateNaturalVector(user->da, &temp);
	PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen2, FILE_MODE_READ, &pviewer);
	VecLoad(user->P,pviewer);
	VecNorm(user->P, NORM_INFINITY, &norm);
	PetscPrintf(PETSC_COMM_WORLD, "PIn %le\n", norm);
	PetscViewerDestroy(&pviewer);
	//VecDestroy(&temp);
  }

  if(nv_once) sprintf(filen2, "%snvfield%06d_%1d.dat", prefix, 0, user->_this);
  else sprintf(filen2, "%snvfield%06d_%1d.dat", prefix, ti, user->_this);
    
  if( !nv_once || (nv_once && ti==tis) )
  {
	PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen2, FILE_MODE_READ, &pviewer);
	VecLoad(user->Nvert,pviewer);
	PetscViewerDestroy(&pviewer);  
  }
	return 0;
}

PetscErrorCode Ucont_P_Binary_Input1(UserCtx *user)
{
  PetscViewer viewer;
  char filen[128];
  
  sprintf(filen, "%sufield%06d_%1d.dat", prefix, ti, user->_this);
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);

  PetscInt N;

  VecGetSize(user->Ucat, &N);
  PetscPrintf(PETSC_COMM_WORLD, "PPP %d\n", N);
  VecLoad(user->Ucat,viewer);
  PetscViewerDestroy(&viewer);

  PetscBarrier(PETSC_NULL);
	return 0;
}

PetscErrorCode Ucat_Binary_Input_Phase_Averaged(UserCtx *user)
{
	PetscViewer viewer;
	char filen[128];
	int N = phase;
	
	sprintf(filen, "%sphase%03d_su0_%06d_%1d.dat", prefix, (int)phase, (int)ti_phase, user->_this);
	if(!file_exist(filen)) PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s !\n", filen), exit(0);
	
	PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
	VecLoad(user->Ucat,viewer);
	PetscViewerDestroy(&viewer);
	
	VecScale(user->Ucat, 1./N);
	
	PetscPrintf(PETSC_COMM_WORLD, "Finished Reading %s\n", filen);
	PetscBarrier(PETSC_NULL);
	
	return 0;
}

PetscErrorCode Ucat_Binary_Input_Averaged(UserCtx *user)
{
	PetscViewer viewer;
	char filen[128];
	int N = tis+1;
	
	sprintf(filen, "%ssu0_%06d_%1d.dat", prefix, ti, user->_this);
	if(!file_exist(filen)) PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s !\n", filen), exit(0);
	
	PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
	VecLoad(user->Ucat,viewer);
	PetscViewerDestroy(&viewer);
	
	VecScale(user->Ucat, 1./N);
	
	PetscPrintf(PETSC_COMM_WORLD, "Finished Reading %s\n", filen);
	PetscBarrier(PETSC_NULL);
		
	if(nv_once) sprintf(filen, "%snvfield%06d_%1d.dat", prefix, 0, user->_this);
	else sprintf(filen, "%snvfield%06d_%1d.dat", prefix, ti, user->_this);

	if( !nv_once || (nv_once && ti==tis) ) {
		PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
		VecLoad(user->Nvert,viewer);
		PetscViewerDestroy(&viewer);
	}
	  
	return 0;
}


PetscErrorCode Ucont_P_Binary_Input_Averaging(UserCtx *user)
{
	PetscViewer viewer;
	char filen[128];
	/*
	sprintf(filen, "su0_%06d_%1d.dat", ti, user->_this);
	PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
	VecLoad(viewer, (user->Ucat_sum));
	PetscViewerDestroy(&viewer);

	sprintf(filen, "su1_%06d_%1d.dat", ti, user->_this);
	PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
	VecLoad(viewer, (user->Ucat_cross_sum));
	PetscViewerDestroy(&viewer);
	
	sprintf(filen, "su2_%06d_%1d.dat", ti, user->_this);
	PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
	VecLoad(viewer, (user->Ucat_square_sum));
	PetscViewerDestroy(&viewer);
	*/
	/*
	sprintf(filen, "sp_%06d_%1d.dat", ti, user->_this);
	PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
	VecLoad(viewer, (user->P));
	PetscViewerDestroy(&viewer);
	*/
  
	if(pcr) {
		Vec Ptmp;
		VecDuplicate(user->P, &Ptmp);
		
		sprintf(filen, "%spfield%06d_%1d.dat", prefix, ti, user->_this);
		PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
		VecLoad(user->P,viewer);
		PetscViewerDestroy(&viewer);
		
		sprintf(filen, "%ssp_%06d_%1d.dat", prefix, ti, user->_this);
		PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
		VecLoad(Ptmp,viewer);
		PetscViewerDestroy(&viewer);
		
		
		VecScale(Ptmp, -1./((double)tis+1.0));
		VecAXPY(user->P, 1., Ptmp);
		
		VecDestroy(&Ptmp);
	}
	
	if(nv_once) sprintf(filen, "%snvfield%06d_%1d.dat", prefix, 0, user->_this);
	else sprintf(filen, "%snvfield%06d_%1d.dat", prefix, ti, user->_this);


	if( !nv_once || (nv_once && ti==tis) )
	  {
	    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
	    VecLoad(user->Nvert,viewer);
	    PetscViewerDestroy(&viewer);
	  }
	/*
	if( !nv_once || (nv_once && ti==tis) ) {
		sprintf(filen, "nvfield%06d_%1d.dat", ti, user->_this);
		PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
		VecLoad(viewer, (user->Nvert));
		PetscViewerDestroy(&viewer);
	}
	*/
	PetscBarrier(PETSC_NULL);
	return 0;
}

#undef __FUNCT__
#define __FUNCT__ "main"

int main(int argc, char **argv)
{
	PetscBool flag;
	
	
	DM	da, fda;
	Vec	qn, qnm;
	Vec	c;
	UserCtx	*user;

	PetscErrorCode ierr;
		
	IBMNodes	*ibm, *ibm0, *ibm1;

	PetscInitialize(&argc, &argv, (char *)0, help);

	
	PetscOptionsInsertFile(PETSC_COMM_WORLD, "pcontrol.dat", PETSC_TRUE);
	
	
	char tmp_str[256];
	//PetscOptionsGetString(PETSC_NULL, "-prefix", tmp_str, 256, &flag);
	
	PetscInt n=3;
	
	PetscOptionsGetReal(PETSC_NULL, "-gx", &gravity_x, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-gy", &gravity_y, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-gz", &gravity_z, PETSC_NULL);

	PetscOptionsGetString(PETSC_NULL, "-subgrid", subgrid, 255, &subgrid_flag);
	

	PetscOptionsGetInt(PETSC_NULL, "-i_begin", &i_begin, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-i_end", &i_end, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-j_begin", &j_begin, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-j_end", &j_end, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-k_begin", &k_begin, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-k_end", &k_end, PETSC_NULL);
		
	PetscOptionsGetInt(PETSC_NULL, "-vc", &vc, PETSC_NULL);	
	PetscOptionsGetInt(PETSC_NULL, "-binary", &binary_input, &flag);
	PetscOptionsGetInt(PETSC_NULL, "-xyz", &xyz_input, &flag);
	PetscOptionsGetInt(PETSC_NULL, "-rans", &rans, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-ransout", &rans_output, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-levelset", &levelset, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-avg", &avg, &flag);
	PetscOptionsGetInt(PETSC_NULL, "-shear", &shear, &flag);
	PetscOptionsGetInt(PETSC_NULL, "-averaging", &averaging_option, &flag);	// from control.dat
	PetscOptionsGetInt(PETSC_NULL, "-phase", &phase, &flag);	// from control.dat
	PetscOptionsGetInt(PETSC_NULL, "-ti_phase", &ti_phase, &flag);	// from control.dat
	
	PetscOptionsGetInt(PETSC_NULL, "-les", &les, &flag);
	PetscOptionsGetInt(PETSC_NULL, "-i_periodic", &i_periodic, &flag);
	PetscOptionsGetInt(PETSC_NULL, "-j_periodic", &j_periodic, &flag);
	PetscOptionsGetInt(PETSC_NULL, "-k_periodic", &k_periodic, &flag);
	
	PetscOptionsGetInt(PETSC_NULL, "-ii_periodic", &i_periodic, &flag);
	PetscOptionsGetInt(PETSC_NULL, "-jj_periodic", &j_periodic, &flag);
	PetscOptionsGetInt(PETSC_NULL, "-kk_periodic", &k_periodic, &flag);
	
	PetscOptionsGetInt(PETSC_NULL, "-vorticity_budget", &vorticity_budget, &flag);
	PetscOptionsGetInt(PETSC_NULL, "-tke_budget", &tke_budget, &flag);
	PetscOptionsGetInt(PETSC_NULL, "-mke_budget", &mke_budget, &flag);
	PetscOptionsGetInt(PETSC_NULL, "-zm_budget", &zm_budget, &flag);
	
	PetscOptionsGetInt(PETSC_NULL, "-nv", &nv_once, &flag);
	printf("nv_once=%d\n", (int)nv_once);
	
	PetscInt  QCR = 0;
	PetscOptionsGetInt(PETSC_NULL, "-qcr", &QCR, PETSC_NULL);
	
	PetscOptionsGetInt(PETSC_NULL, "-tis", &tis, &flag);
	if (!flag) PetscPrintf(PETSC_COMM_WORLD, "Need the starting number!\n");    
	
	PetscOptionsGetInt(PETSC_NULL, "-tie", &tie, &flag);
	if (!flag) tie = tis;
    
	PetscOptionsGetInt(PETSC_NULL, "-ts", &tsteps, &flag);
	if (!flag) tsteps = 5; /* Default increasement is 5 */
    
	PetscOptionsGetInt(PETSC_NULL, "-onlyV", &onlyV, &flag);
	PetscOptionsGetInt(PETSC_NULL, "-iavg", &i_average, &flag);
	PetscOptionsGetInt(PETSC_NULL, "-javg", &j_average, &flag);
	PetscOptionsGetInt(PETSC_NULL, "-kavg", &k_average, &flag);
	PetscOptionsGetInt(PETSC_NULL, "-ikavg", &ik_average, &flag);
	PetscOptionsGetInt(PETSC_NULL, "-pcr", &pcr, &flag);
	PetscOptionsGetInt(PETSC_NULL, "-reynolds", &reynolds, &flag);
	
	PetscOptionsGetInt(PETSC_NULL, "-ikcavg", &ikc_average, &flag);
	if(flag) {
		PetscBool flag1, flag2;
		PetscOptionsGetInt(PETSC_NULL, "-pi", &pi, &flag1);
		PetscOptionsGetInt(PETSC_NULL, "-pk", &pk, &flag2);
		
		if(!flag1 || !flag2) {
			printf("To use -ikcavg you must set -pi and -pk, which are number of points in i- and k- directions.\n");
			exit(0);
		}
	}
	
  
	if(pcr) avg=1;
	if(i_average) avg=1;
	if(j_average) avg=1;
	if(k_average) avg=1;
	if(ik_average) avg=1;
	if(ikc_average) avg=1;
 

	PetscOptionsGetReal(PETSC_NULL, "-r_tb", &r_tb, PETSC_NULL); // xiaolei add
	PetscOptionsGetReal(PETSC_NULL, "-xc_tb", &xc_tb, PETSC_NULL); // xiaolei add
	PetscOptionsGetReal(PETSC_NULL, "-yc_tb", &yc_tb, PETSC_NULL); // xiaolei add
	PetscOptionsGetInt(PETSC_NULL, "-DiskAvg", &DiskAvg, PETSC_NULL);  // xiaolei add 

  
	if(i_average + j_average + k_average >1) PetscPrintf(PETSC_COMM_WORLD, "Iavg and Javg cannot be set together !! !\n"), exit(0);
      
	int rank, bi;

	PetscMalloc(sizeof(IBMNodes), &ibm); // xiaolei add 
	PetscMalloc(sizeof(IBMNodes), &ibm0);
	PetscMalloc(sizeof(IBMNodes), &ibm1);

	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	
	if(subgrid_flag) {
		sprintf(prefix, "subgrid_%s_", subgrid);
		if(!xyz_input) binary_input=1;
	}
	else sprintf(prefix, "");

	if(xyz_input) {block_number=1;}
	else {
		FILE *fd;
		char str[256];
		
		if(subgrid_flag) sprintf(str, "subgrid_%s.dat", subgrid);
		else sprintf(str, "grid.dat");
		
		fd = fopen(str, "r");
		if(fd==NULL) PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!! !\n", str), exit(0);
		
		if(binary_input) fread(&block_number, sizeof(int), 1, fd);
		else fscanf(fd, "%i\n", &block_number);
		MPI_Bcast(&block_number, 1, MPI_INT, 0, PETSC_COMM_WORLD);
		fclose(fd);
	}
	
	
	
	PetscMalloc(block_number*sizeof(UserCtx), &user);
	PetscOptionsGetReal(PETSC_NULL, "-ren", &user->ren, PETSC_NULL);

	ReadCoordinates(user);

	PetscPrintf(PETSC_COMM_WORLD, "read coord!\n");

	for (bi=0; bi<block_number; bi++) {
		DMCreateGlobalVector(user[bi].da, &user[bi].Nvert);
		if(shear) {
			DMCreateGlobalVector(user[bi].fda, &user[bi].Csi);
			DMCreateGlobalVector(user[bi].fda, &user[bi].Eta);
			DMCreateGlobalVector(user[bi].fda, &user[bi].Zet);
			DMCreateGlobalVector(user[bi].da, &user[bi].Aj);
			FormMetrics(&(user[bi]));
			
			Calc_avg_shear_stress(&(user[bi]));
						
			VecDestroy(&user[bi].Csi);
			VecDestroy(&user[bi].Eta);
			VecDestroy(&user[bi].Zet);
			VecDestroy(&user[bi].Aj);
			exit(0);
		}
		else if(!avg || QCR) {
			DMCreateGlobalVector(user[bi].da, &user[bi].P);
			DMCreateGlobalVector(user[bi].fda, &user[bi].Ucat);
			if(!vc) DMCreateGlobalVector(user[bi].fda, &user[bi].Ucat_o);
			
			if(QCR) {
				if(QCR==1 || QCR==2) {
					DMCreateGlobalVector(user[bi].fda, &user[bi].Csi);
					DMCreateGlobalVector(user[bi].fda, &user[bi].Eta);
					DMCreateGlobalVector(user[bi].fda, &user[bi].Zet);
					DMCreateGlobalVector(user[bi].da, &user[bi].Aj);
					FormMetrics(&(user[bi]));
				}
			}
		}
		else {
			if(pcr) {
				DMCreateGlobalVector(user[bi].da, &user[bi].P);
				DMCreateGlobalVector(user[bi].fda, &user[bi].Ucat);
			}
			else if(avg==1) {
				if(averaging_option==3 || vorticity_budget || tke_budget || DiskAvg || mke_budget || zm_budget) {
					DMCreateGlobalVector(user[bi].fda, &user[bi].Csi);
					DMCreateGlobalVector(user[bi].fda, &user[bi].Eta);
					DMCreateGlobalVector(user[bi].fda, &user[bi].Zet);
					DMCreateGlobalVector(user[bi].da, &user[bi].Aj);
					FormMetrics(&(user[bi]));
				}
			}
		}
		
	}


  
	if(avg) {
		if(i_average) PetscPrintf(PETSC_COMM_WORLD, "Averaging in I direction!\n");
		else if(j_average) PetscPrintf(PETSC_COMM_WORLD, "Averaging in J direction!\n");
		else if(k_average) PetscPrintf(PETSC_COMM_WORLD, "Averaging in K direction!\n");
		else if(ik_average) PetscPrintf(PETSC_COMM_WORLD, "Averaging in IK direction!\n");
		else PetscPrintf(PETSC_COMM_WORLD, "Averaging !\n");
		/*
		DMCreateGlobalVector(user[bi].fda, &user->Ucat_sum);
		DMCreateGlobalVector(user[bi].fda, &user->Ucat_cross_sum);
		DMCreateGlobalVector(user[bi].fda, &user->Ucat_square_sum);
		*/
		
	}
  

 
	for (ti=tis; ti<=tie; ti+=tsteps) {
		for (bi=0; bi<block_number; bi++) {
			if(phase && QCR) {
				printf("Begin Ucat_Binary_Input_Phase_Averaged\n");
				Ucat_Binary_Input_Phase_Averaged(&user[bi]);
				VecSet(user[bi].Nvert, 0);
				printf("End Ucat_Binary_Input_Phase_Averaged\n");
			}
			else if(avg && QCR) {
				Ucat_Binary_Input_Averaged(&user[bi]);
			}
			else if(avg) {
				Ucont_P_Binary_Input_Averaging(&user[bi]);
			}
			else {
				Ucont_P_Binary_Input(&user[bi]);
				Ucont_P_Binary_Input1(&user[bi]);
			}
		}


		if(QCR) {
			TECIOOutQ(user, QCR);
		}
		else {
			if(avg && vorticity_budget) TECIOOut_Vorticity_Budget(user);
			else if(avg && tke_budget) TECIOOut_TKE_Budget(user);
			else if(avg && mke_budget) TECIOOut_MKE_Budget1(user);
			else if(avg && zm_budget) TECIOOut_zMomentum_Budget(user);	
			else if(avg && DiskAvg) TECIOOut_AveragCirc(user, &ibm[0]);
			else if(avg) TECIOOut_Averaging(user);
			else TECIOOut_V(user, onlyV);
			//TECIOOut(user);
		}
		
	}


	PetscFinalize();
}



PetscErrorCode ReadCoordinates(UserCtx *user)
{
	Cmpnts ***coor;

	Vec Coor;
	int bi, i, j, k, IM, JM, KM;
	PetscReal *gc;
	FILE *fd;
	PetscReal	d0 = 1.;
	int rank;
	
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	PetscReal	cl = 1.;
	PetscOptionsGetReal(PETSC_NULL, "-cl", &cl, PETSC_NULL);

	char str[256];
	
	
	if(xyz_input) {
		if(subgrid_flag) sprintf(str, "subgrid_xyz_%s.dat", subgrid);
		else sprintf(str, "xyz.dat");
	}
	else {
		if(subgrid_flag) sprintf(str, "subgrid_%s.dat", subgrid);
		else sprintf(str, "grid.dat");
	}
	
	fd = fopen(str, "r");
	
	if(fd==NULL) printf("Cannot open %s !\n", str),exit(0);

	printf("Begin Reading %s !\n", str);
	  
	if(xyz_input) {i=1;}
	else if(binary_input) {
		fread(&i, sizeof(int), 1, fd);
		if(i!=1) PetscPrintf(PETSC_COMM_WORLD, "This seems to be a text file !\n"),exit(0);
	}
	else {
		fscanf(fd, "%i\n", &i);
		if(i!=1) PetscPrintf(PETSC_COMM_WORLD, "This seems to be a binary file !\n"),exit(0);
	}
	

	for (bi=block_number-1; bi>=0; bi--) {
	  
		std::vector<double> X, Y,Z;
		double tmp;
		
		if(xyz_input) {
			fscanf(fd, "%i %i %i\n", &(user[bi].IM), &(user[bi].JM), &(user[bi].KM));
			X.resize(user[bi].IM);
			Y.resize(user[bi].JM);
			Z.resize(user[bi].KM);
			
			for (i=0; i<user[bi].IM; i++) fscanf(fd, "%le %le %le\n", &X[i], &tmp, &tmp);
			for (j=0; j<user[bi].JM; j++) fscanf(fd, "%le %le %le\n", &tmp, &Y[j], &tmp);
			for (k=0; k<user[bi].KM; k++) fscanf(fd, "%le %le %le\n", &tmp, &tmp, &Z[k]);
		}
		else if(binary_input) {
			fread(&(user[bi].IM), sizeof(int), 1, fd);
			fread(&(user[bi].JM), sizeof(int), 1, fd);
			fread(&(user[bi].KM), sizeof(int), 1, fd);
		}
		else fscanf(fd, "%i %i %i\n", &(user[bi].IM), &(user[bi].JM), &(user[bi].KM));
		
		IM = user[bi].IM; JM = user[bi].JM; KM = user[bi].KM;
    
    
		PetscInt m = PETSC_DECIDE, n = PETSC_DECIDE, p = PETSC_DECIDE, s=2;
		DMDABoundaryType bx=DMDA_BOUNDARY_NONE, by=DMDA_BOUNDARY_NONE, bz=DMDA_BOUNDARY_NONE;
		
		DMDACreate3d(PETSC_COMM_WORLD, bx, by, bz, DMDA_STENCIL_BOX, user[bi].IM+1, user[bi].JM+1, user[bi].KM+1, m, n, p, 1, s, PETSC_NULL, PETSC_NULL, PETSC_NULL, &(user[bi].da));
		if(rans) {
			DMDACreate3d(PETSC_COMM_WORLD, bx, by, bz, DMDA_STENCIL_BOX, user[bi].IM+1, user[bi].JM+1, user[bi].KM+1, m, n, p, 2, s, PETSC_NULL, PETSC_NULL, PETSC_NULL, &(user[bi].fda2));
		}
		/*
		DMDACreate3d(PETSC_COMM_WORLD, DA_NONPERIODIC, DA_STENCIL_BOX,
				user[bi].IM+1, user[bi].JM+1, user[bi].KM+1, 1,1,
				PETSC_DECIDE, 1, 2, PETSC_NULL, PETSC_NULL, PETSC_NULL,
				&(user[bi].da));
		if(rans) {
			DMDACreate3d(PETSC_COMM_WORLD, DA_NONPERIODIC, DA_STENCIL_BOX,
				user[bi].IM+1, user[bi].JM+1, user[bi].KM+1, 1,1,
				PETSC_DECIDE, 2, 2, PETSC_NULL, PETSC_NULL, PETSC_NULL,
				&(user[bi].fda2));
		}
		*/
		
		DMDASetUniformCoordinates(user[bi].da, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
		DMDAGetCoordinateDA(user[bi].da, &(user[bi].fda));
	
		DMDAGetLocalInfo(user[bi].da, &(user[bi].info));

		DMDALocalInfo	info = user[bi].info;
		PetscInt	xs = info.xs, xe = info.xs + info.xm;
		PetscInt  	ys = info.ys, ye = info.ys + info.ym;
		PetscInt	zs = info.zs, ze = info.zs + info.zm;
		PetscInt	mx = info.mx, my = info.my, mz = info.mz;
		
		DMDAGetGhostedCoordinates(user[bi].da, &Coor);
		DMDAVecGetArray(user[bi].fda, Coor, &coor);
		
		double buffer;
		
		for (k=0; k<KM; k++)
		for (j=0; j<JM; j++)
		for (i=0; i<IM; i++) {
			if(xyz_input) {}
			else if(binary_input) fread(&buffer, sizeof(double), 1, fd);
			else fscanf(fd, "%le", &buffer);
				
			if( k>=zs && k<=ze && j>=ys && j<ye && i>=xs && i<xe ) {
				if(xyz_input) coor[k][j][i].x = X[i]/cl;
				else coor[k][j][i].x = buffer/cl;
			}
		}
			
		for (k=0; k<KM; k++)
		for (j=0; j<JM; j++)
		for (i=0; i<IM; i++) {
			if(xyz_input) {}
			else if(binary_input) fread(&buffer, sizeof(double), 1, fd);
			else fscanf(fd, "%le", &buffer);
				
			if( k>=zs && k<=ze && j>=ys && j<ye && i>=xs && i<xe ) {
				if(xyz_input) coor[k][j][i].y = Y[j]/cl;
				else coor[k][j][i].y = buffer/cl;
			}
		}
	
		for (k=0; k<KM; k++)
		for (j=0; j<JM; j++)
		for (i=0; i<IM; i++) {
			if(xyz_input) {}
			else if(binary_input) fread(&buffer, sizeof(double), 1, fd);
			else fscanf(fd, "%le", &buffer);
				
			if( k>=zs && k<=ze && j>=ys && j<ye && i>=xs && i<xe ) {
				if(xyz_input) coor[k][j][i].z = Z[k]/cl;
				else coor[k][j][i].z = buffer/cl;
			}
		}

		DMDAVecRestoreArray(user[bi].fda, Coor, &coor);

		Vec	gCoor;
		DMDAGetCoordinates(user[bi].da, &gCoor);
		DMLocalToGlobalBegin(user[bi].fda, Coor, INSERT_VALUES, gCoor);
		DMLocalToGlobalEnd(user[bi].fda, Coor, INSERT_VALUES, gCoor);
		DMGlobalToLocalBegin(user[bi].fda, gCoor, INSERT_VALUES, Coor);
		DMGlobalToLocalEnd(user[bi].fda, gCoor, INSERT_VALUES, Coor);

	}
	
	fclose(fd);
	
	printf("Finish Reading %s !\n", str);
  
	for (bi=0; bi<block_number; bi++) {
		user[bi]._this = bi;
	}
	return(0);
}

void Calc_avg_shear_stress(UserCtx *user)
{
	double N=(double)tis+1.0;
  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  PetscInt i, j, k;
  lxs = xs; lxe = xe; lys = ys; lye = ye; lzs = zs; lze = ze;

  Cmpnts ***usum, ***csi, ***eta, ***zet;
  PetscReal ***aj, ***psum, ***nvert;
  if (lxs==0) lxs++;
  if (lxe==mx) lxe--;
  if (lys==0) lys++;
  if (lye==my) lye--;
  if (lzs==0) lzs++;
  if (lze==mz) lze--;

	char filen[256];
	PetscViewer	viewer;
		
	Vec P_sum;
	DMCreateGlobalVector(user->da, &P_sum);
	DMCreateGlobalVector(user->fda, &user->Ucat_sum);
  	  
  ti=tis;
  sprintf(filen, "%ssu0_%06d_%1d.dat", prefix, ti, user->_this);
	PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
	VecLoad(user->Ucat_sum,viewer);
	PetscViewerDestroy(&viewer);
	
	ti=tis;		
	sprintf(filen, "%ssp_%06d_%1d.dat", prefix, ti, user->_this);
	PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
	VecLoad(P_sum,viewer);
	PetscViewerDestroy(&viewer);

  DMDAVecGetArray(user->fda, user->Csi, &csi);
  DMDAVecGetArray(user->fda, user->Eta, &eta);
  DMDAVecGetArray(user->fda, user->Zet, &zet);
  DMDAVecGetArray(user->da, user->Aj, &aj);
  DMDAVecGetArray(user->da, user->Nvert, &nvert);
  DMDAVecGetArray(user->fda, user->Ucat_sum, &usum);
  DMDAVecGetArray(user->da, P_sum, &psum);


	double force_skin_bottom = 0;
	double force_pressure_bottom = 0;
	double force_bottom = 0;
	double area_bottom = 0;
	
	double force_skin_top = 0;
	double force_pressure_top = 0;
	double force_top = 0;
	double area_top = 0;
			
	j=0;
	for (k=lzs; k<lze; k++)
	for (i=lxs; i<lxe; i++) {
		if (nvert[k][j+1][i] < 0.1) {
			double du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz;
			double dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz;
			
			dudc=0, dvdc=0, dwdc=0;
			
			dude=usum[k][j+1][i].x * 2.0 / N;
			dvde=usum[k][j+1][i].y * 2.0 / N;
			dwde=usum[k][j+1][i].z * 2.0 / N;
			
			dudz=0, dvdz=0, dwdz=0;
			
			double ajc = aj[k][j+1][i];
			double csi0 = csi[k][j+1][i].x, csi1 = csi[k][j+1][i].y, csi2 = csi[k][j+1][i].z;
			double eta0 = eta[k][j+1][i].x, eta1 = eta[k][j+1][i].y, eta2 = eta[k][j+1][i].z;
			double zet0 = zet[k][j+1][i].x, zet1 = zet[k][j+1][i].y, zet2 = zet[k][j+1][i].z;

			Compute_du_dxyz (csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, 
					dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz,
					&du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );

			double j_area = sqrt( eta[k][j+1][i].x*eta[k][j+1][i].x + eta[k][j+1][i].y*eta[k][j+1][i].y + eta[k][j+1][i].z*eta[k][j+1][i].z );
			double ni[3], nj[3], nk[3];
			double nx, ny, nz;
			Calculate_normal(csi[k][j+1][i], eta[k][j+1][i], zet[k][j+1][i], ni, nj, nk);
			nx = nj[0]; //inward normal
			ny = nj[1]; //inward normal
			nz = nj[2]; //inward normal
			
			
			double Fp = - psum[k][j+1][i] * eta2 / N;
			double Fs = (dw_dx * nx + dw_dy * ny + dw_dz * nz) / user->ren * j_area;
			//double Fs = (du_dx * nx + du_dy * ny + du_dz * nz) / user->ren * j_area;
			
			force_skin_bottom += Fs;
			force_pressure_bottom += Fp;
			force_bottom += Fs + Fp;
			area_bottom += fabs(eta1);	// projected area
		}
	}

	j=my-2;
	for (k=lzs; k<lze; k++)
	for (i=lxs; i<lxe; i++) {
		if (nvert[k][j][i] < 0.1) {
			double du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz;
			double dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz;
			
			dudc=0, dvdc=0, dwdc=0;
			
			dude = -usum[k][j][i].x * 2.0 / N;
			dvde = -usum[k][j][i].y * 2.0 / N;
			dwde = -usum[k][j][i].z * 2.0 / N;
			
			dudz=0, dvdz=0, dwdz=0;
			
			double ajc = aj[k][j][i];
			double csi0 = csi[k][j][i].x, csi1 = csi[k][j][i].y, csi2 = csi[k][j][i].z;
			double eta0 = eta[k][j][i].x, eta1 = eta[k][j][i].y, eta2 = eta[k][j][i].z;
			double zet0 = zet[k][j][i].x, zet1 = zet[k][j][i].y, zet2 = zet[k][j][i].z;

			Compute_du_dxyz (csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, 
					dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz,
					&du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );

			double j_area = sqrt( eta[k][j][i].x*eta[k][j][i].x + eta[k][j][i].y*eta[k][j][i].y + eta[k][j][i].z*eta[k][j][i].z );
			double ni[3], nj[3], nk[3];
			double nx, ny, nz;
			Calculate_normal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);
			nx = -nj[0]; //inward normal
			ny = -nj[1]; //inward normal
			nz = -nj[2]; //inward normal
			
			
			double Fp = - psum[k][j][i] * eta2 / N;
			double Fs = (dw_dx * nx + dw_dy * ny + dw_dz * nz) / user->ren * j_area;
			//double Fs = (du_dx * nx + du_dy * ny + du_dz * nz) / user->ren * j_area;
			
			force_skin_top += Fs;
			force_pressure_top += Fp;
			force_top += Fs + Fp;
			area_top += fabs(eta1);	// projected area
		}
	}
	
  DMDAVecRestoreArray(user->fda, user->Csi, &csi);
  DMDAVecRestoreArray(user->fda, user->Eta, &eta);
  DMDAVecRestoreArray(user->fda, user->Zet, &zet);
  DMDAVecRestoreArray(user->da, user->Aj, &aj);
  DMDAVecRestoreArray(user->da, user->Nvert, &nvert);
  DMDAVecRestoreArray(user->fda, user->Ucat_sum, &usum);
  DMDAVecRestoreArray(user->da, P_sum, &psum);

	VecDestroy(&P_sum);
  VecDestroy(&user->Ucat_sum);
  
  printf("Top:\tarea=%f, force=%f, skin force=%f, pressure force=%f\n",
  			area_top, force_top, force_skin_top, force_pressure_top);
  			
	printf("\tstress=%f, skin stress=%f, pressure stress=%f, u*=%f, Re*=%f\n",
  			force_top/area_top, force_skin_top/area_top, force_pressure_top/area_top,
  			sqrt(fabs(force_top/area_top)), sqrt(fabs(force_top/area_top))*user->ren);
  			
  printf("\n");
  
  printf("Bottom:\tarea=%f, force=%f, skin force=%f, pressure force=%f\n",
  			area_bottom, force_bottom, force_skin_bottom, force_pressure_bottom);
  			
	printf("\tstress=%f, skin stress=%f, pressure stress=%f, u*=%f, Re*=%f\n",
  			force_bottom/area_bottom, force_skin_bottom/area_bottom, force_pressure_bottom/area_bottom,
  			sqrt(fabs(force_bottom/area_bottom)), sqrt(fabs(force_bottom/area_bottom))*user->ren);
}

PetscErrorCode Lambda2(UserCtx *user)
{
  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  PetscInt i, j, k;
  lxs = xs; lxe = xe; lys = ys; lye = ye; lzs = zs; lze = ze;

  Cmpnts ***ucat, ***csi, ***eta, ***zet;
  PetscReal ***aj, ***q, ***nvert;
  if (lxs==0) lxs++;
  if (lxe==mx) lxe--;
  if (lys==0) lys++;
  if (lye==my) lye--;
  if (lzs==0) lzs++;
  if (lze==mz) lze--;

  DMDAVecGetArray(user->fda, user->Ucat, &ucat);
  DMDAVecGetArray(user->fda, user->Csi, &csi);
  DMDAVecGetArray(user->fda, user->Eta, &eta);
  DMDAVecGetArray(user->fda, user->Zet, &zet);

  DMDAVecGetArray(user->da, user->Aj, &aj);
  DMDAVecGetArray(user->da, user->Nvert, &nvert);
  DMDAVecGetArray(user->da, user->P, &q);

  PetscReal uc, vc, wc, ue, ve, we, uz, vz, wz;

  PetscReal s11, s12, s13, s21, s22, s23, s31, s32, s33;
  PetscReal d11, d12, d13, d21, d22, d23, d31, d32, d33;

  PetscReal w11, w12, w13, w21, w22, w23, w31, w32, w33;
  //PetscReal so, wo;
  PetscReal csi1, csi2, csi3, eta1, eta2, eta3, zet1, zet2, zet3;
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	      
	double dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz;
	double du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz;
	double csi0 = csi[k][j][i].x, csi1 = csi[k][j][i].y, csi2 = csi[k][j][i].z;
	double eta0= eta[k][j][i].x, eta1 = eta[k][j][i].y, eta2 = eta[k][j][i].z;
	double zet0 = zet[k][j][i].x, zet1 = zet[k][j][i].y, zet2 = zet[k][j][i].z;
	double ajc = aj[k][j][i];
	
	Compute_du_center (i, j, k, mx, my, mz, ucat, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
	Compute_du_dxyz (	csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz,
							&du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
							
	double Sxx = 0.5*( du_dx + du_dx ), Sxy = 0.5*(du_dy + dv_dx), Sxz = 0.5*(du_dz + dw_dx);
	double Syx = Sxy, Syy = 0.5*(dv_dy + dv_dy),	Syz = 0.5*(dv_dz + dw_dy);
	double Szx = Sxz, Szy=Syz, Szz = 0.5*(dw_dz + dw_dz);

	
	w11 = 0;
	w12 = 0.5*(du_dy - dv_dx);
	w13 = 0.5*(du_dz - dw_dx);
	w21 = -w12;
	w22 = 0.;
	w23 = 0.5*(dv_dz - dw_dy);
	w31 = -w13;
	w32 = -w23;
	w33 = 0.;
	
	
	double S[3][3], W[3][3], D[3][3];
	
	D[0][0] = du_dx, D[0][1] = du_dy, D[0][2] = du_dz;
	D[1][0] = dv_dx, D[1][1] = dv_dy, D[1][2] = dv_dz;
	D[2][0] = dw_dx, D[2][1] = dw_dy, D[2][2] = dw_dz;
	
	S[0][0] = Sxx;
	S[0][1] = Sxy;
	S[0][2] = Sxz;

	S[1][0] = Syx;
	S[1][1] = Syy;
	S[1][2] = Syz;

	S[2][0] = Szx;
	S[2][1] = Szy;
	S[2][2] = Szz;

	W[0][0] = w11;
	W[0][1] = w12;
	W[0][2] = w13;
	W[1][0] = w21;
	W[1][1] = w22;
	W[1][2] = w23;
	W[2][0] = w31;
	W[2][1] = w32;
	W[2][2] = w33;
	
	// lambda-2
	double A[3][3], V[3][3], d[3];
	
	for(int row=0; row<3; row++)
	for(int col=0; col<3; col++) A[row][col]=0;
	
	for(int row=0; row<3; row++)
	for(int col=0; col<3; col++) {
		A[row][col] += S[row][0] * S[0][col];
		A[row][col] += S[row][1] * S[1][col];
		A[row][col] += S[row][2] * S[2][col];
	}
	
	for(int row=0; row<3; row++)
	for(int col=0; col<3; col++) {
		A[row][col] += W[row][0] * W[0][col];
		A[row][col] += W[row][1] * W[1][col];
		A[row][col] += W[row][2] * W[2][col];
	}
	
	if(nvert[k][j][i]<0.1) {
		eigen_decomposition(A, V, d);
		q[k][j][i] = d[1];
	}
	else q[k][j][i] = 1000.0;
/*	
	// delta criterion
	double DD[3][3];
	for(int row=0; row<3; row++)
	for(int col=0; col<3; col++) DD[row][col]=0;
	
	for(int row=0; row<3; row++)
	for(int col=0; col<3; col++) {
		DD[row][col] += D[row][0] * D[0][col];
		DD[row][col] += D[row][1] * D[1][col];
		DD[row][col] += D[row][2] * D[2][col];
	}
	double tr_DD = DD[0][0] + DD[1][1] + DD[2][2];
	double det_D = D[0][0]*(D[2][2]*D[1][1]-D[2][1]*D[1][2])-D[1][0]*(D[2][2]*D[0][1]-D[2][1]*D[0][2])+D[2][0]*(D[1][2]*D[0][1]-D[1][1]*D[0][2]);
	
	//double Q = -0.5*tr_DD;
	
	double SS=0, WW=0;
	for(int row=0; row<3; row++)
	for(int col=0; col<3; col++) {
		SS+=S[row][col]*S[row][col];
		WW+=W[row][col]*W[row][col];
	}
	double Q = 0.5*(WW - SS);
	
	double R = - det_D;
	if(nvert[k][j][i]<0.1) {
		q[k][j][i] = pow( 0.5*R, 2. )  + pow( Q/3., 3.);
	}
	else q[k][j][i] = -10;
	if(q[k][j][i]<0) q[k][j][i]=-10;
*/	
      }
    }
  }

  DMDAVecRestoreArray(user->fda, user->Ucat, &ucat);
  DMDAVecRestoreArray(user->fda, user->Csi, &csi);
  DMDAVecRestoreArray(user->fda, user->Eta, &eta);
  DMDAVecRestoreArray(user->fda, user->Zet, &zet);

  DMDAVecRestoreArray(user->da, user->Aj, &aj);
  DMDAVecRestoreArray(user->da, user->Nvert, &nvert);
  DMDAVecRestoreArray(user->da, user->P, &q);

  return 0;
}

PetscErrorCode QCriteria(UserCtx *user)
{

  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  PetscInt i, j, k;
  lxs = xs; lxe = xe; lys = ys; lye = ye; lzs = zs; lze = ze;

  Cmpnts ***ucat, ***csi, ***eta, ***zet;
  PetscReal ***aj, ***q, ***nvert;
  if (lxs==0) lxs++;
  if (lxe==mx) lxe--;
  if (lys==0) lys++;
  if (lye==my) lye--;
  if (lzs==0) lzs++;
  if (lze==mz) lze--;
	
  DMDAVecGetArray(user->fda, user->Ucat, &ucat);
  DMDAVecGetArray(user->fda, user->Csi, &csi);
  DMDAVecGetArray(user->fda, user->Eta, &eta);
  DMDAVecGetArray(user->fda, user->Zet, &zet);

  DMDAVecGetArray(user->da, user->Aj, &aj);
  DMDAVecGetArray(user->da, user->Nvert, &nvert);
  DMDAVecGetArray(user->da, user->P, &q);

  PetscReal uc, vc, wc, ue, ve, we, uz, vz, wz;

  PetscReal s11, s12, s13, s21, s22, s23, s31, s32, s33;
  PetscReal d11, d12, d13, d21, d22, d23, d31, d32, d33;

  PetscReal w11, w12, w13, w21, w22, w23, w31, w32, w33;
  PetscReal so, wo;
  PetscReal csi1, csi2, csi3, eta1, eta2, eta3, zet1, zet2, zet3;
  
  
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	double dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz;
	double du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz;
	double csi0 = csi[k][j][i].x, csi1 = csi[k][j][i].y, csi2 = csi[k][j][i].z;
	double eta0= eta[k][j][i].x, eta1 = eta[k][j][i].y, eta2 = eta[k][j][i].z;
	double zet0 = zet[k][j][i].x, zet1 = zet[k][j][i].y, zet2 = zet[k][j][i].z;
	double ajc = aj[k][j][i];
	
	Compute_du_center (i, j, k, mx, my, mz, ucat, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
	Compute_du_dxyz (	csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz,
							&du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
							
	double Sxx = 0.5*( du_dx + du_dx ), Sxy = 0.5*(du_dy + dv_dx), Sxz = 0.5*(du_dz + dw_dx);
	double Syx = Sxy, Syy = 0.5*(dv_dy + dv_dy),	Syz = 0.5*(dv_dz + dw_dy);
	double Szx = Sxz, Szy=Syz, Szz = 0.5*(dw_dz + dw_dz);
	so = Sxx*Sxx + Sxy*Sxy + Sxz*Sxz + Syx*Syx + Syy*Syy + Syz*Syz + Szx*Szx + Szy*Szy + Szz*Szz;
	
	w11 = 0;
	w12 = 0.5*(du_dy - dv_dx);
	w13 = 0.5*(du_dz - dw_dx);
	w21 = -w12;
	w22 = 0.;
	w23 = 0.5*(dv_dz - dw_dy);
	w31 = -w13;
	w32 = -w23;
	w33 = 0.;
	
	wo = w11*w11 + w12*w12 + w13*w13 + w21*w21 + w22*w22 + w23*w23 + w31*w31 + w32*w32 + w33*w33;

/*
	so = ( d11 *  d11 + d22 * d22 + d33 * d33) + 0.5* ( (d12 + d21) * (d12 + d21) + (d13 + d31) * (d13 + d31) + (d23 + d32) * (d23 + d32) );
	wo = 0.5 * ( (d12 - d21)*(d12 - d21) + (d13 - d31) * (d13 - d31) + (d23 - d32) * (d23 - d32) );
	V19=0.5 * ( (V13 - V11)*(V13 - V11) + (V16 - V12) * (V16 - V12) + (V17 - V15) * (V17 - V15) ) - 0.5 * ( V10 *  V10 + V14 * V14 + V18 * V18) - 0.25* ( (V13 + V11) * (V13 + V11) + (V16 + V12) * (V16 + V12) + (V17 + V15) * (V17 + V15) )
*/
	
	if( nvert[k][j][i]>0.1 ) q[k][j][i] = -100;
	else q[k][j][i] = (wo - so) / 2.;
      }
    }
  }

	DMDAVecRestoreArray(user->fda, user->Ucat, &ucat);
	DMDAVecRestoreArray(user->fda, user->Csi, &csi);
	DMDAVecRestoreArray(user->fda, user->Eta, &eta);
	DMDAVecRestoreArray(user->fda, user->Zet, &zet);

	DMDAVecRestoreArray(user->da, user->Aj, &aj);
	DMDAVecRestoreArray(user->da, user->Nvert, &nvert);
	DMDAVecRestoreArray(user->da, user->P, &q);

  return 0;
}

PetscErrorCode Velocity_Magnitude(UserCtx *user)	// store at P
{
  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  PetscInt i, j, k;
  lxs = xs; lxe = xe; lys = ys; lye = ye; lzs = zs; lze = ze;

  Cmpnts ***ucat;
  PetscReal ***aj, ***q, ***nvert;
  if (lxs==0) lxs++;
  if (lxe==mx) lxe--;
  if (lys==0) lys++;
  if (lye==my) lye--;
  if (lzs==0) lzs++;
  if (lze==mz) lze--;

  DMDAVecGetArray(user->fda, user->Ucat, &ucat);
  DMDAVecGetArray(user->da, user->P, &q);

  PetscReal uc, vc, wc, ue, ve, we, uz, vz, wz;
  PetscReal d11, d12, d13, d21, d22, d23, d31, d32, d33;

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	q[k][j][i] = sqrt( ucat[k][j][i].x*ucat[k][j][i].x + ucat[k][j][i].y*ucat[k][j][i].y + ucat[k][j][i].z*ucat[k][j][i].z );
      }
    }
  }

  DMDAVecRestoreArray(user->fda, user->Ucat, &ucat);
  DMDAVecRestoreArray(user->da, user->P, &q);

  return 0;
}


#ifdef MAX
#undef MAX
#endif

#define MAX(a, b) ((a)>(b)?(a):(b))

#define n 3

static double hypot2(double x, double y) {
  return sqrt(x*x+y*y);
}

// Symmetric Householder reduction to tridiagonal form.

static void tred2(double V[n][n], double d[n], double e[n]) {

//  This is derived from the Algol procedures tred2 by
//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
//  Fortran subroutine in EISPACK.

  for (int j = 0; j < n; j++) {
    d[j] = V[n-1][j];
  }

  // Householder reduction to tridiagonal form.

  for (int i = n-1; i > 0; i--) {

    // Scale to avoid under/overflow.

    double scale = 0.0;
    double h = 0.0;
    for (int k = 0; k < i; k++) {
      scale = scale + fabs(d[k]);
    }
    if (scale == 0.0) {
      e[i] = d[i-1];
      for (int j = 0; j < i; j++) {
        d[j] = V[i-1][j];
        V[i][j] = 0.0;
        V[j][i] = 0.0;
      }
    } else {

      // Generate Householder vector.

      for (int k = 0; k < i; k++) {
        d[k] /= scale;
        h += d[k] * d[k];
      }
      double f = d[i-1];
      double g = sqrt(h);
      if (f > 0) {
        g = -g;
      }
      e[i] = scale * g;
      h = h - f * g;
      d[i-1] = f - g;
      for (int j = 0; j < i; j++) {
        e[j] = 0.0;
      }

      // Apply similarity transformation to remaining columns.

      for (int j = 0; j < i; j++) {
        f = d[j];
        V[j][i] = f;
        g = e[j] + V[j][j] * f;
        for (int k = j+1; k <= i-1; k++) {
          g += V[k][j] * d[k];
          e[k] += V[k][j] * f;
        }
        e[j] = g;
      }
      f = 0.0;
      for (int j = 0; j < i; j++) {
        e[j] /= h;
        f += e[j] * d[j];
      }
      double hh = f / (h + h);
      for (int j = 0; j < i; j++) {
        e[j] -= hh * d[j];
      }
      for (int j = 0; j < i; j++) {
        f = d[j];
        g = e[j];
        for (int k = j; k <= i-1; k++) {
          V[k][j] -= (f * e[k] + g * d[k]);
        }
        d[j] = V[i-1][j];
        V[i][j] = 0.0;
      }
    }
    d[i] = h;
  }

  // Accumulate transformations.

  for (int i = 0; i < n-1; i++) {
    V[n-1][i] = V[i][i];
    V[i][i] = 1.0;
    double h = d[i+1];
    if (h != 0.0) {
      for (int k = 0; k <= i; k++) {
        d[k] = V[k][i+1] / h;
      }
      for (int j = 0; j <= i; j++) {
        double g = 0.0;
        for (int k = 0; k <= i; k++) {
          g += V[k][i+1] * V[k][j];
        }
        for (int k = 0; k <= i; k++) {
          V[k][j] -= g * d[k];
        }
      }
    }
    for (int k = 0; k <= i; k++) {
      V[k][i+1] = 0.0;
    }
  }
  for (int j = 0; j < n; j++) {
    d[j] = V[n-1][j];
    V[n-1][j] = 0.0;
  }
  V[n-1][n-1] = 1.0;
  e[0] = 0.0;
} 

// Symmetric tridiagonal QL algorithm.

static void tql2(double V[n][n], double d[n], double e[n]) {

//  This is derived from the Algol procedures tql2, by
//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
//  Fortran subroutine in EISPACK.

  for (int i = 1; i < n; i++) {
    e[i-1] = e[i];
  }
  e[n-1] = 0.0;

  double f = 0.0;
  double tst1 = 0.0;
  double eps = pow(2.0,-52.0);
  for (int l = 0; l < n; l++) {

    // Find small subdiagonal element

    tst1 = MAX(tst1,fabs(d[l]) + fabs(e[l]));
    int m = l;
    while (m < n) {
      if (fabs(e[m]) <= eps*tst1) {
        break;
      }
      m++;
    }

    // If m == l, d[l] is an eigenvalue,
    // otherwise, iterate.

    if (m > l) {
      int iter = 0;
      do {
        iter = iter + 1;  // (Could check iteration count here.)

        // Compute implicit shift

        double g = d[l];
        double p = (d[l+1] - g) / (2.0 * e[l]);
        double r = hypot2(p,1.0);
        if (p < 0) {
          r = -r;
        }
        d[l] = e[l] / (p + r);
        d[l+1] = e[l] * (p + r);
        double dl1 = d[l+1];
        double h = g - d[l];
        for (int i = l+2; i < n; i++) {
          d[i] -= h;
        }
        f = f + h;

        // Implicit QL transformation.

        p = d[m];
        double c = 1.0;
        double c2 = c;
        double c3 = c;
        double el1 = e[l+1];
        double s = 0.0;
        double s2 = 0.0;
        for (int i = m-1; i >= l; i--) {
          c3 = c2;
          c2 = c;
          s2 = s;
          g = c * e[i];
          h = c * p;
          r = hypot2(p,e[i]);
          e[i+1] = s * r;
          s = e[i] / r;
          c = p / r;
          p = c * d[i] - s * g;
          d[i+1] = h + s * (c * g + s * d[i]);

          // Accumulate transformation.

          for (int k = 0; k < n; k++) {
            h = V[k][i+1];
            V[k][i+1] = s * V[k][i] + c * h;
            V[k][i] = c * V[k][i] - s * h;
          }
        }
        p = -s * s2 * c3 * el1 * e[l] / dl1;
        e[l] = s * p;
        d[l] = c * p;

        // Check for convergence.

      } while (fabs(e[l]) > eps*tst1);
    }
    d[l] = d[l] + f;
    e[l] = 0.0;
  }
  
  // Sort eigenvalues and corresponding vectors.

  for (int i = 0; i < n-1; i++) {
    int k = i;
    double p = d[i];
    for (int j = i+1; j < n; j++) {
      if (d[j] < p) {
        k = j;
        p = d[j];
      }
    }
    if (k != i) {
      d[k] = d[i];
      d[i] = p;
      for (int j = 0; j < n; j++) {
        p = V[j][i];
        V[j][i] = V[j][k];
        V[j][k] = p;
      }
    }
  }
}

void eigen_decomposition(double A[n][n], double V[n][n], double d[n]) {
  double e[n];
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      V[i][j] = A[i][j];
    }
  }
  tred2(V, d, e);
  tql2(V, d, e);
}

// xiaolei add 
PetscErrorCode disk_read_ucd(IBMNodes *ibm, int ibi, char fname[80], PetscReal zc_tb)
{
  	int	rank;
  	int	n_v , n_elmt ;
  	PetscReal	*x_bp , *y_bp , *z_bp ;
  	int	*nv1 , *nv2 , *nv3 ;
  	PetscReal	*nf_x, *nf_y, *nf_z;
  	int	i,ii;
  	int	n1e, n2e, n3e;
  	PetscReal	dx12, dy12, dz12, dx13, dy13, dz13;
  	PetscReal     dr;

  	char   ss[20];
  	char string[128];

  	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  	if(!rank) { // root processor read in the data
    		FILE *fd;
    		PetscPrintf(PETSC_COMM_SELF, "READ %s\n", fname);
    		char filen[160];  
    		sprintf(filen,"./%s" ,fname);
 
    		fd = fopen(filen, "r"); 
    		if (!fd) printf("Cannot open %s !!", filen),exit(0);
    		else printf("Opened %s !\n", filen);
    		n_v =0;

    		if (fd) {
      			fgets(string, 128, fd);
      			fgets(string, 128, fd);
      			fgets(string, 128, fd);

      			fscanf(fd, "%i %i %i %i %i",&n_v,&n_elmt,&ii,&ii,&ii);
      			PetscPrintf(PETSC_COMM_SELF, "number of nodes & elements %d %d\n",n_v, n_elmt);
      
      			ibm->n_v = n_v;
      			ibm->n_elmt = n_elmt;      
      
      			MPI_Bcast(&(ibm->n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
      
      			PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp));	// added by seokkoo 03.04.2009
      			PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp));
      			PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp));
      
      			x_bp = ibm->x_bp;	// seokkoo
      			y_bp = ibm->y_bp;
      			z_bp = ibm->z_bp;

      			for (i=0; i<n_v; i++) {
				fscanf(fd, "%i %le %le %le", &ii, &x_bp[i], &y_bp[i], &z_bp[i]);//, &t, &t, &t);
	
				x_bp[i]=x_bp[i]/reflength_wt;
				y_bp[i]=y_bp[i]/reflength_wt;
				z_bp[i]=z_bp[i]/reflength_wt;

				ibm->x_bp[i] = x_bp[i]+xc_tb;
				ibm->y_bp[i] = y_bp[i]+yc_tb;
				ibm->z_bp[i] = z_bp[i]+zc_tb;

			}

      			MPI_Bcast(ibm->x_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      			MPI_Bcast(ibm->y_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      			MPI_Bcast(ibm->z_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      
      			MPI_Bcast(&(ibm->n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);

      			PetscMalloc(n_elmt*sizeof(int), &nv1);
      			PetscMalloc(n_elmt*sizeof(int), &nv2);
      			PetscMalloc(n_elmt*sizeof(int), &nv3);
      
      			PetscMalloc(n_elmt*sizeof(PetscReal), &nf_x);
      			PetscMalloc(n_elmt*sizeof(PetscReal), &nf_y);
      			PetscMalloc(n_elmt*sizeof(PetscReal), &nf_z);
      
      			PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->dA); //Area

      			PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_x));
      			PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_y));
      			PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_z));

      			PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dhx));
      			PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dhy));
      			PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dhz));


		        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->uu));
		        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->vv));
		        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->ww));
		        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->uv));
		        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->vw));
		        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->wu));
		        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->val));

		        PetscMalloc(n_elmt*sizeof(int), &(ibm->i_min));
		        PetscMalloc(n_elmt*sizeof(int), &(ibm->j_min));
		        PetscMalloc(n_elmt*sizeof(int), &(ibm->k_min));

		        PetscMalloc(n_elmt*sizeof(int), &(ibm->i_max));
		        PetscMalloc(n_elmt*sizeof(int), &(ibm->j_max));
		        PetscMalloc(n_elmt*sizeof(int), &(ibm->k_max));

      			for (i=0; i<n_elmt; i++) {
				fscanf(fd, "%i %i %s %i %i %i\n", &ii,&ii, &ss, nv1+i, nv2+i, nv3+i);
				nv1[i] = nv1[i] - 1; nv2[i] = nv2[i]-1; nv3[i] = nv3[i] - 1;
		      	}
		        ibm->nv1 = nv1; ibm->nv2 = nv2; ibm->nv3 = nv3;


		        fclose(fd);
 	      	}
      
	      	for (i=0; i<n_elmt; i++) {
      
			n1e = ibm->nv1[i]; n2e =ibm->nv2[i]; n3e = ibm->nv3[i];
		        dx12 = x_bp[n2e] - x_bp[n1e];
		        dy12 = y_bp[n2e] - y_bp[n1e];
		        dz12 = z_bp[n2e] - z_bp[n1e];
      
		        dx13 = x_bp[n3e] - x_bp[n1e];
		        dy13 = y_bp[n3e] - y_bp[n1e];
		        dz13 = z_bp[n3e] - z_bp[n1e];
      
		        nf_x[i] = dy12 * dz13 - dz12 * dy13;
		        nf_y[i] = -dx12 * dz13 + dz12 * dx13;
		        nf_z[i] = dx12 * dy13 - dy12 * dx13;
      
		        dr = sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i] + nf_z[i]*nf_z[i]);
      
		        nf_x[i] /=dr; nf_y[i]/=dr; nf_z[i]/=dr;
      
		        ibm->dA[i] = dr/2.; 
      
      			ibm->cent_x[i]= (x_bp[n1e]+x_bp[n2e]+x_bp[n3e])/3.;
		        ibm->cent_y[i]= (y_bp[n1e]+y_bp[n2e]+y_bp[n3e])/3.;
		        ibm->cent_z[i]= (z_bp[n1e]+z_bp[n2e]+z_bp[n3e])/3.;	
		}
    
    
    		ibm->nf_x = nf_x; ibm->nf_y = nf_y;  ibm->nf_z = nf_z;
    
	        MPI_Bcast(ibm->nv1, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
	        MPI_Bcast(ibm->nv2, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
	        MPI_Bcast(ibm->nv3, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    
    	        MPI_Bcast(ibm->nf_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
	        MPI_Bcast(ibm->nf_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
	        MPI_Bcast(ibm->nf_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
      
    		MPI_Bcast(ibm->dA, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
    		MPI_Bcast(ibm->cent_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    		MPI_Bcast(ibm->cent_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    		MPI_Bcast(ibm->cent_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

  	}
  	else if (rank) {
    		n_v =0;

      		MPI_Bcast(&(ibm->n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
     
		n_v=ibm->n_v; 
      		PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp));	// added by seokkoo 03.04.2009
      		PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp));
      		PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp));
      
      		MPI_Bcast(ibm->x_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      		MPI_Bcast(ibm->y_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      		MPI_Bcast(ibm->z_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      
      		MPI_Bcast(&(ibm->n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);
		n_elmt=ibm->n_elmt;
		
      		PetscMalloc(n_elmt*sizeof(int), &ibm->nv1);
      		PetscMalloc(n_elmt*sizeof(int), &ibm->nv2);
      		PetscMalloc(n_elmt*sizeof(int), &ibm->nv3);
      
      		PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->nf_x);
      		PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->nf_y);
      		PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->nf_z);
      
      		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dA)); //Area

      		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_x));
      		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_y));
      		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_z));

      		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dhx));
      		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dhy));
      		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dhz));


		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->uu));
	        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->vv));
		PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->ww));
	        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->uv));
	        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->vw));
	        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->wu));
	        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->val));

	        PetscMalloc(n_elmt*sizeof(int), &(ibm->i_min));
	        PetscMalloc(n_elmt*sizeof(int), &(ibm->j_min));
	        PetscMalloc(n_elmt*sizeof(int), &(ibm->k_min));

	        PetscMalloc(n_elmt*sizeof(int), &(ibm->i_max));
	        PetscMalloc(n_elmt*sizeof(int), &(ibm->j_max));
	        PetscMalloc(n_elmt*sizeof(int), &(ibm->k_max));

    
	        MPI_Bcast(ibm->nv1, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
	        MPI_Bcast(ibm->nv2, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
	        MPI_Bcast(ibm->nv3, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    
    	        MPI_Bcast(ibm->nf_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
	        MPI_Bcast(ibm->nf_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
	        MPI_Bcast(ibm->nf_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
      
    		MPI_Bcast(ibm->dA, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
    		MPI_Bcast(ibm->cent_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    		MPI_Bcast(ibm->cent_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    		MPI_Bcast(ibm->cent_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);



	}
	PetscPrintf(PETSC_COMM_WORLD, "Read %s !\n", fname);

  	return(0);
}



// Call in main.c, after reading the grid for turbines
PetscErrorCode Pre_process(UserCtx *user, IBMNodes *ibm)
{
  	DM              da = user->da, fda = user->fda;
  	DMDALocalInfo  	info;
  	PetscInt        xs, xe, ys, ye, zs, ze; // Local grid information
  	PetscInt        mx, my, mz; // Dimensions in three directions
  	PetscInt        i, j, k, l, ibi;
  	PetscInt	lxs, lxe, lys, lye, lzs, lze;

  	Cmpnts		***coor, ***csi, ***eta, ***zet;

  	PetscReal 	***aj;

  	Vec		Coor;


	printf("here 1\n");
  	DMDAGetLocalInfo(da, &info);
  	mx = info.mx; my = info.my; mz = info.mz;
  	xs = info.xs; xe = xs + info.xm;
  	ys = info.ys; ye = ys + info.ym;
  	zs = info.zs; ze = zs + info.zm;

  	lxs = xs; lxe = xe;
  	lys = ys; lye = ye;
  	lzs = zs; lze = ze;

  	if (xs==0) lxs = xs+1;
  	if (ys==0) lys = ys+1;
  	if (zs==0) lzs = zs+1;

  	if (xe==mx) lxe = xe-1;
  	if (ye==my) lye = ye-1;
  	if (ze==mz) lze = ze-1;


	double *dhx_, *dhy_, *dhz_;
	int *iclose, *jclose, *kclose;

	double *sum_dhx_, *sum_dhy_, *sum_dhz_;
	int *sum_iclose, *sum_jclose, *sum_kclose;


  	DMDAGetGhostedCoordinates(da, &Coor);
  	DMDAVecGetArray(fda, Coor, &coor);

  	DMDAVecGetArray(fda, user->Csi,  &csi);
  	DMDAVecGetArray(fda, user->Eta,  &eta);
  	DMDAVecGetArray(fda, user->Zet,  &zet);
  	DMDAVecGetArray(da,  user->Aj,  &aj);


	printf("here 2\n");
  	int n1e, n2e, n3e;

  	for (ibi=0; ibi<NumberOfObjects; ibi++) {

		dhx_= (double *) malloc(ibm[ibi].n_elmt*sizeof(double));
		dhy_= (double *) malloc(ibm[ibi].n_elmt*sizeof(double));
		dhz_= (double *) malloc(ibm[ibi].n_elmt*sizeof(double));

		iclose= (int *) malloc(ibm[ibi].n_elmt*sizeof(int));
		jclose= (int *) malloc(ibm[ibi].n_elmt*sizeof(int));
		kclose= (int *) malloc(ibm[ibi].n_elmt*sizeof(int));

		sum_dhx_= (double *) malloc(ibm[ibi].n_elmt*sizeof(double));
		sum_dhy_= (double *) malloc(ibm[ibi].n_elmt*sizeof(double));
		sum_dhz_= (double *) malloc(ibm[ibi].n_elmt*sizeof(double));

		sum_iclose= (int *) malloc(ibm[ibi].n_elmt*sizeof(int));
		sum_jclose= (int *) malloc(ibm[ibi].n_elmt*sizeof(int));
		sum_kclose= (int *) malloc(ibm[ibi].n_elmt*sizeof(int));

	printf("here 3\n");
	        for (l=0; l<ibm[ibi].n_elmt; l++) {

	               	dhx_[l] = 0.0;
	                dhy_[l] = 0.0;
	                dhz_[l] = 0.0;

	               	sum_dhx_[l] = 0.0;
        	        sum_dhy_[l] = 0.0;
                	sum_dhz_[l] = 0.0;

			iclose[l]=0;
			jclose[l]=0;
			kclose[l]=0;

			sum_iclose[l]=0;
			sum_jclose[l]=0;
			sum_kclose[l]=0;

		}

	printf("here 3-1\n");

	        for (l=0; l<ibm[ibi].n_elmt; l++) {
			int imark=0;
			double dmin=1.e6;
			int ic,jc,kc;
	                for (i=lxs; i<lxe; i++)
        	        for (j=lys; j<lye; j++) 
                	for (k=lzs; k<lze; k++){ 
			
				double r1=ibm[ibi].cent_x[l]-coor[k][j][i].x, r2=ibm[ibi].cent_y[l]-coor[k][j][i].y, r3=ibm[ibi].cent_z[l]-coor[k][j][i].z; 
				double d1=sqrt(r1*r1+r2*r2+r3*r3);
				if (d1<dmin) {
					dmin=d1;
					ic=i; jc=j; kc=k;
				}
                	}
	
			double dmin_global;
			MPI_Allreduce (&dmin, &dmin_global, 1, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD);
			double diff=fabs(dmin-dmin_global);
			if (diff>1.e-6) {
				ic=0; jc=0; kc=0;
				iclose[l]=0; jclose[l]=0; kclose[l]=0;
				dhx_[l]=0.0;
				dhy_[l]=0.0;
				dhz_[l]=0.0;
			} else {

				iclose[l]=ic; jclose[l]=jc; kclose[l]=kc;
				i=ic; j=jc; k=kc;
		                double area = sqrt( csi[k][j][i].x*csi[k][j][i].x + csi[k][j][i].y*csi[k][j][i].y + csi[k][j][i].z*csi[k][j][i].z );
				dhx_[l]=1.0/aj[k][j][i]/area;

			       	area = sqrt( eta[k][j][i].x*eta[k][j][i].x + eta[k][j][i].y*eta[k][j][i].y + eta[k][j][i].z*eta[k][j][i].z );
				dhy_[l]=1.0/aj[k][j][i]/area;

			       	area = sqrt( zet[k][j][i].x*zet[k][j][i].x + zet[k][j][i].y*zet[k][j][i].y + zet[k][j][i].z*zet[k][j][i].z );
				dhz_[l]=1.0/aj[k][j][i]/area;	
			}

	        }


	printf("here 4\n");

	  	PetscBarrier(PETSC_NULL);

		MPI_Allreduce (&dhx_[0], &sum_dhx_[0], ibm[ibi].n_elmt, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
		MPI_Allreduce (&dhy_[0], &sum_dhy_[0], ibm[ibi].n_elmt, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
		MPI_Allreduce (&dhz_[0], &sum_dhz_[0], ibm[ibi].n_elmt, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

		MPI_Allreduce (&iclose[0], &sum_iclose[0], ibm[ibi].n_elmt, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);
		MPI_Allreduce (&jclose[0], &sum_jclose[0], ibm[ibi].n_elmt, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);
		MPI_Allreduce (&kclose[0], &sum_kclose[0], ibm[ibi].n_elmt, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);


	        for (l=0; l<ibm[ibi].n_elmt; l++) {

        	       	ibm[ibi].dhx[l]=sum_dhx_[l];
	               	ibm[ibi].dhy[l]=sum_dhy_[l];
        	       	ibm[ibi].dhz[l]=sum_dhz_[l];

	               	dhx_[l]=sum_dhx_[l];
	               	dhy_[l]=sum_dhy_[l];
	               	dhz_[l]=sum_dhz_[l];

			iclose[l]=sum_iclose[l];
			jclose[l]=sum_jclose[l];
			kclose[l]=sum_kclose[l];

		}

	printf("here 5\n");

	        int iii1 = 0;
		int iii2 = 0;
	        int iii = 0;

	        for (l=0; l<ibm[ibi].n_elmt; l++) {

                        int n1e = ibm[ibi].nv1[i], n2e = ibm[ibi].nv2[i], n3e = ibm[ibi].nv3[i];

			/*
			double d12=sqrt(pow(ibm[ibi].x_bp[n1e]-ibm[ibi].x_bp[n2e],2)+pow(ibm[ibi].y_bp[n1e]-ibm[ibi].y_bp[n2e],2)+pow(ibm[ibi].z_bp[n1e]-ibm[ibi].z_bp[n2e],2));
			double d23=sqrt(pow(ibm[ibi].x_bp[n3e]-ibm[ibi].x_bp[n2e],2)+pow(ibm[ibi].y_bp[n3e]-ibm[ibi].y_bp[n2e],2)+pow(ibm[ibi].z_bp[n3e]-ibm[ibi].z_bp[n2e],2));
			double d31=sqrt(pow(ibm[ibi].x_bp[n1e]-ibm[ibi].x_bp[n3e],2)+pow(ibm[ibi].y_bp[n1e]-ibm[ibi].y_bp[n3e],2)+pow(ibm[ibi].z_bp[n1e]-ibm[ibi].z_bp[n3e],2));

			double d_max=max(d12,d23);
			d_max=max(d_max,d31);
			
			Itpwidth=(int)d_max/ibm[ibi].dhx[l]+4;
			*/

			int ic=iclose[l];
			int jc=jclose[l];
			int kc=kclose[l];

			int ic1=ic-Itpwidth, ic2=ic+Itpwidth; 
		
			if (ic1>lxe||ic2<lxs) {
				ibm[ibi].i_min[l]=lxs;
				ibm[ibi].i_max[l]=lxs;
			} else {
				ibm[ibi].i_min[l]=PetscMax(ic1, lxs);
				ibm[ibi].i_max[l]=PetscMin(ic2, lxe);
			}

			// Itpwidth=(int)d_max/ibm[ibi].dhy[l]+4;
			int jc1=jc-Itpwidth, jc2=jc+Itpwidth; 

	       		if (jc1>lye||jc2<lys) {
				ibm[ibi].j_min[l]=lys;
				ibm[ibi].j_max[l]=lys;
			} else {
				ibm[ibi].j_min[l]=PetscMax(jc1, lys);
				ibm[ibi].j_max[l]=PetscMin(jc2, lye);
			}

			// Itpwidth=(int)d_max/ibm[ibi].dhz[l]+4;
			int kc1=kc-Itpwidth, kc2=kc+Itpwidth; 

	       		if (kc1>lze||kc2<lzs) {
				ibm[ibi].k_min[l]=lzs;
				ibm[ibi].k_max[l]=lzs;
			} else {
				ibm[ibi].k_min[l]=PetscMax(kc1, lzs);
				ibm[ibi].k_max[l]=PetscMin(kc2, lze);
			}


		}

	printf("here 6\n");
/*
	        for (l=0; l<ibm[ibi].n_elmt; l++) {

			ibm[ibi].i_min[l]=lxs;
			ibm[ibi].i_max[l]=lxe;

			ibm[ibi].j_min[l]=lys;
			ibm[ibi].j_max[l]=lye;


			ibm[ibi].k_min[l]=lzs;
			ibm[ibi].k_max[l]=lze;


		}	
*/

		/*
		if (rotor_model == 2) {

		int imin_g=1000000, imax_g=-1000000, jmin_g=1000000, jmax_g=-1000000, kmin_g=1000000, kmax_g=-1000000;

	        for (l=0; l<ibm[ibi].n_elmt; l++) {
			if (imin_g>ibm[ibi].i_min[l]) imin_g=ibm[ibi].i_min[l];	
			if (imax_g<ibm[ibi].i_max[l]) imax_g=ibm[ibi].i_max[l];	

			if (jmin_g>ibm[ibi].j_min[l]) jmin_g=ibm[ibi].j_min[l];	
			if (jmax_g<ibm[ibi].j_max[l]) jmax_g=ibm[ibi].j_max[l];	

			if (kmin_g>ibm[ibi].k_min[l]) kmin_g=ibm[ibi].k_min[l];	
			if (kmax_g<ibm[ibi].k_max[l]) kmax_g=ibm[ibi].k_max[l];	
		}

	        for (l=0; l<ibm[ibi].n_elmt; l++) {

			ibm[ibi].i_min[l]=imin_g;	
			ibm[ibi].j_min[l]=jmin_g;	
			ibm[ibi].k_min[l]=kmin_g;	

			ibm[ibi].i_max[l]=imax_g;	
			ibm[ibi].j_max[l]=jmax_g;	
			ibm[ibi].k_max[l]=kmax_g;	

		}	

		}
		*/

		free(dhx_);
		free(dhy_); 
		free(dhz_);
	        free(iclose);
		free(jclose);
		free(kclose);

		free(sum_dhx_);
		free(sum_dhy_); 
		free(sum_dhz_);
	        free(sum_iclose);
		free(sum_jclose);
		free(sum_kclose);

	}

  	DMDAVecRestoreArray(fda, Coor, &coor);

	DMDAVecRestoreArray(fda, user->Csi,  &csi);
  	DMDAVecRestoreArray(fda, user->Eta,  &eta);
  	DMDAVecRestoreArray(fda, user->Zet,  &zet);
  	DMDAVecRestoreArray(da,  user->Aj,  &aj);

  	return(0);

}



PetscErrorCode TECIOOut_AveragCirc(UserCtx *user, IBMNodes *ibm)
{
	PetscInt bi;

	int i,j,k;

	double N=(double)tis+1.0;
	
	for (bi=0; bi<block_number; bi++) {

		DM da = user[bi].da, fda = user[bi].fda;
		DMDALocalInfo info = user[bi].info;

		PetscInt	xs = info.xs, xe = info.xs + info.xm;
		PetscInt  	ys = info.ys, ye = info.ys + info.ym;
		PetscInt	zs = info.zs, ze = info.zs + info.zm;
		PetscInt	mx = info.mx, my = info.my, mz = info.mz;
    
		Cmpnts		***u2sum, ***u1sum,  ***usum, ***coor, ***csi, ***eta, ***zet;
		PetscReal	***aj;
		//VecScatter ctx;

		double uu[mz], vv[mz], ww[mz], uv[mz], vw[mz], wu[mz], U[mz], V[mz], W[mz], SN[mz];
		PetscViewer viewer;
		char filen[128];
			
		DMCreateGlobalVector(user[bi].fda, &user[bi].Ucat_sum);
		DMCreateGlobalVector(user[bi].fda, &user[bi].Ucat_cross_sum);
		DMCreateGlobalVector(user[bi].fda, &user[bi].Ucat_square_sum);

		sprintf(filen, "%ssu0_%06d_%1d.dat", prefix, ti, user[bi]._this);
		PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
		VecLoad(user[bi].Ucat_sum,viewer);
		PetscViewerDestroy(&viewer);

		sprintf(filen, "%ssu1_%06d_%1d.dat", prefix, ti, user[bi]._this);
		PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
		VecLoad(user[bi].Ucat_cross_sum,viewer);
		PetscViewerDestroy(&viewer);
			
		sprintf(filen, "%ssu2_%06d_%1d.dat", prefix, ti, user[bi]._this);
		PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
		VecLoad(user[bi].Ucat_square_sum,viewer);
		PetscViewerDestroy(&viewer);

		Vec Coor;
		DMDAGetCoordinates(da, &Coor);
		DMDAVecGetArray(fda, Coor, &coor);
		DMDAVecGetArray(user[bi].fda, user[bi].Ucat_sum, &usum);
		DMDAVecGetArray(user[bi].fda, user[bi].Ucat_cross_sum, &u1sum);
		DMDAVecGetArray(user[bi].fda, user[bi].Ucat_square_sum, &u2sum);
		DMDAVecGetArray(user[bi].fda, user[bi].Csi, &csi);
		DMDAVecGetArray(user[bi].fda, user[bi].Eta, &eta);
		DMDAVecGetArray(user[bi].fda, user[bi].Zet, &zet);
		DMDAVecGetArray(user[bi].da, user[bi].Aj, &aj);
		
	
		for (k=zs; k<ze-2; k++) {
			double A_tot=0.0, uu_tot=0.0, vv_tot=0.0, ww_tot=0.0, U_tot=0.0, V_tot=0.0, W_tot=0.0, uv_tot=0.0, vw_tot=0.0, wu_tot=0.0, Gt=0.0, Gn=0.0;
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) {
				double rr=sqrt(pow(coor[k+1][j+1][i+1].x-xc_tb,2)+pow(coor[k+1][j+1][i+1].y-yc_tb,2));
			       	double area=sqrt( zet[k+1][j+1][i+1].x*zet[k+1][j+1][i+1].x + zet[k+1][j+1][i+1].y*zet[k+1][j+1][i+1].y + zet[k+1][j+1][i+1].z*zet[k+1][j+1][i+1].z );
				double nt_y=(coor[k+1][j+1][i+1].x-xc_tb)/rr;
				double nt_x=-(coor[k+1][j+1][i+1].y-yc_tb)/rr;

				if (rr<r_tb) {
					A_tot+=area;
					U_tot += area*usum[k+1][j+1][i+1].x/N;
					V_tot += area*usum[k+1][j+1][i+1].y/N;
					W_tot += area*usum[k+1][j+1][i+1].z/N;
					double _U = usum[k+1][j+1][i+1].x/N;
					uu_tot += area*( u2sum[k+1][j+1][i+1].x/N - _U*_U );
					double _V = usum[k+1][j+1][i+1].y/N;
					vv_tot += area*( u2sum[k+1][j+1][i+1].y/N - _V*_V );
					double _W = usum[k+1][j+1][i+1].z/N;
					ww_tot += area*( u2sum[k+1][j+1][i+1].z/N - _W*_W );
					double UV = usum[k+1][j+1][i+1].x*usum[k+1][j+1][i+1].y / (N*N);
					uv_tot += area*u1sum[k+1][j+1][i+1].x/N - UV;
					double VW = usum[k+1][j+1][i+1].y*usum[k+1][j+1][i+1].z / (N*N);
					vw_tot += area*u1sum[k+1][j+1][i+1].y/N - VW;
					double WU = usum[k+1][j+1][i+1].z*usum[k+1][j+1][i+1].x / (N*N);
					wu_tot += area*u1sum[k+1][j+1][i+1].z/N - WU;

					double Ut=_U*nt_x+_V*nt_y;
					Gt+=rr*Ut*_W*area;
					Gn+=_W*_W*area;
				}
			}
			U[k+1]=U_tot/A_tot;
			V[k+1]=V_tot/A_tot;
			W[k+1]=W_tot/A_tot;
			uu[k+1]=uu_tot/A_tot;
			vv[k+1]=vv_tot/A_tot;
			ww[k+1]=ww_tot/A_tot;
			uv[k+1]=uv_tot/A_tot;
			vw[k+1]=vw_tot/A_tot;
			wu[k+1]=wu_tot/A_tot;
			SN[k+1]=Gt/(r_tb*Gn+1.e-20);
			
		}

	
			
		FILE *f;
		char fname[128];
	   	sprintf(fname, "DiskAveraged.plt");
		f = fopen(fname, "w");
	      	PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=\"z\", \"U\", \"V\", \"W\", \"uu\", \"vv\", \"ww\", \"uv\", \"vw\", \"wu\", \"Sn\"\n");

		i=1;j=1;
		for (k=zs+1; k<ze-1; k++) {
      			PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le %le %le %le %le %le %le %le %le\n", coor[k][j][i].z, U[k], V[k], W[k], uu[k], vv[k], ww[k], uv[k], vw[k], wu[k], SN[k]);
		}
	      	fclose(f);

		DMDAVecRestoreArray(fda, Coor, &coor);
		DMDAVecRestoreArray(user[bi].fda, user[bi].Ucat_sum, &usum);
		DMDAVecRestoreArray(user[bi].fda, user[bi].Ucat_cross_sum, &u1sum);
		DMDAVecRestoreArray(user[bi].fda, user[bi].Ucat_square_sum, &u2sum);
		DMDAVecRestoreArray(user[bi].fda, user[bi].Csi, &csi);
		DMDAVecRestoreArray(user[bi].fda, user[bi].Eta, &eta);
		DMDAVecRestoreArray(user[bi].fda, user[bi].Zet, &zet);
		DMDAVecRestoreArray(user[bi].da, user[bi].Aj, &aj);
			
		VecDestroy(&user[bi].Ucat_sum);
		VecDestroy(&user[bi].Ucat_cross_sum);
		VecDestroy(&user[bi].Ucat_square_sum);
			

	}

	return 0;
}



PetscErrorCode TECIOOut_MKE_Budget(UserCtx *user)
{
	PetscInt bi;

	char filen[80];
	
	if(phase) sprintf(filen, "%sMKEBudget_Phase%03d-avg.plt", prefix, (int)phase);
	else sprintf(filen, "%sMKEBudget%06d-avg.plt", prefix, (int)ti);

	printf("MKE Budget !\n");
	
	INTEGER4 I, Debug, VIsDouble, DIsDouble, III, IMax, JMax, KMax;
	VIsDouble = 0;
	DIsDouble = 0;
	Debug = 0;

	I = TECINI100((char*)"Averaging", "X Y Z MKE TurbulenceProduction Dissipation ViscousDiffusion PressureTransport TurbulentConvection MeanConvection",  filen, (char*)".",  &Debug,  &VIsDouble);
	
	for (bi=0; bi<block_number; bi++) {
		DM da = user[bi].da, fda = user[bi].fda;
		DMDALocalInfo info = user[bi].info;

		PetscInt	xs = info.xs, xe = info.xs + info.xm;
		PetscInt  	ys = info.ys, ye = info.ys + info.ym;
		PetscInt	zs = info.zs, ze = info.zs + info.zm;
		PetscInt	mx = info.mx, my = info.my, mz = info.mz;
    
		PetscInt	lxs, lys, lzs, lxe, lye, lze;
		PetscInt	i, j, k;
		Cmpnts	***coor, ***u2sum, ***u1sum,  ***usum;
		Cmpnts	***ucat, ***csi, ***eta, ***zet; 
		PetscReal	***nvert, ***uu, ***vv, ***ww, ***uv, ***vw, ***wu, ***aj, ***uu_mean, ***vv_mean, ***ww_mean;
		Vec Coor;
		Vec Ucat, UU, VV, WW, UV, VW, WU, UU_Mean, VV_Mean, WW_Mean;	// define vectors

		INTEGER4	ZoneType=0, ICellMax=0, JCellMax=0, KCellMax=0;
		INTEGER4	IsBlock=1, NumFaceConnections=0, FaceNeighborMode=0;
		INTEGER4    ShareConnectivityFromZone=0;
		INTEGER4	LOC[100] = {1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; 
		/* 1 is cell-centered   0 is node centered */

		IMax = mx-1; JMax = my-1; KMax = mz-1;

		I = TECZNE100((char*)"Block 1",
			&ZoneType, 	/* Ordered zone */
			&IMax,
			&JMax,
			&KMax,
			&ICellMax,
			&JCellMax,
			&KCellMax,
			&IsBlock,	/* ISBLOCK  BLOCK format */
			&NumFaceConnections,
			&FaceNeighborMode,
			LOC,
			NULL,
			&ShareConnectivityFromZone); /* No connectivity sharing */

		float *x;
	
		x = new float [(mx-1)*(my-1)*(mz-1)];
		III = (mx-1) * (my-1) * (mz-1);

		DMDAVecGetArray(user[bi].fda, user[bi].Csi, &csi);
		DMDAVecGetArray(user[bi].fda, user[bi].Eta, &eta);
		DMDAVecGetArray(user[bi].fda, user[bi].Zet, &zet);
		DMDAVecGetArray(user[bi].da, user[bi].Aj, &aj);
		DMDAVecGetArray(user[bi].da, user[bi].Nvert, &nvert);
				
		DMDAGetCoordinates(da, &Coor);
		DMDAVecGetArray(fda, Coor, &coor);


//
		// X
		for (k=zs; k<ze-1; k++)
		for (j=ys; j<ye-1; j++)
		for (i=xs; i<xe-1; i++) x[k * (mx-1)*(my-1) + j*(mx-1) + i] = coor[k][j][i].x;
		I = TECDAT100(&III, x, &DIsDouble);

		// Y
		for (k=zs; k<ze-1; k++)
		for (j=ys; j<ye-1; j++)
		for (i=xs; i<xe-1; i++)	x[k * (mx-1)*(my-1) + j*(mx-1) + i] = coor[k][j][i].y;
		I = TECDAT100(&III, x, &DIsDouble);

		// Z
		for (k=zs; k<ze-1; k++)
		for (j=ys; j<ye-1; j++)
		for (i=xs; i<xe-1; i++)	x[k * (mx-1)*(my-1) + j*(mx-1) + i] = coor[k][j][i].z;
		I = TECDAT100(&III, x, &DIsDouble);
	
    
		double N=(double)tis+1.0;
		if(phase) N = phase;
	
		III = (mx-2) * (my-2) * (mz-2);

		PetscViewer viewer;
		char filen[128];
		
		DMCreateGlobalVector(user[bi].fda, &user[bi].Ucat_sum);
		DMCreateGlobalVector(user[bi].fda, &user[bi].Ucat_cross_sum);
		DMCreateGlobalVector(user[bi].fda, &user[bi].Ucat_square_sum);
		
		DMCreateGlobalVector(user[bi].fda, &Ucat);
		DMCreateGlobalVector(user[bi].da, &UU);
		DMCreateGlobalVector(user[bi].da, &VV);
		DMCreateGlobalVector(user[bi].da, &WW);
		DMCreateGlobalVector(user[bi].da, &UV);
		DMCreateGlobalVector(user[bi].da, &VW);
		DMCreateGlobalVector(user[bi].da, &WU);
		DMCreateGlobalVector(user[bi].da, &UU_Mean);
		DMCreateGlobalVector(user[bi].da, &VV_Mean);
		DMCreateGlobalVector(user[bi].da, &WW_Mean);
		/*
		if(phase) {
			sprintf(filen, "phase%03d_su0_%06d_%1d.dat", phase, ti_phase, user[bi]._this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoad(viewer, (user[bi].Ucat_sum));
			PetscViewerDestroy(&viewer);
				
			PetscPrintf(PETSC_COMM_WORLD, "Read %s\n", filen);

			sprintf(filen, "phase%03d_su1_%06d_%1d.dat", phase, ti_phase, user[bi]._this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoad(viewer, (user[bi].Ucat_cross_sum));
			PetscViewerDestroy(&viewer);
				
			PetscPrintf(PETSC_COMM_WORLD, "Read %s\n", filen);
			
			sprintf(filen, "phase%03d_su2_%06d_%1d.dat", phase, ti_phase, user[bi]._this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoad(viewer, (user[bi].Ucat_square_sum));
			PetscViewerDestroy(&viewer);
				
			PetscPrintf(PETSC_COMM_WORLD, "Read %s\n", filen);
		}
		else*/ 
		{
			sprintf(filen, "%ssu0_%06d_%1d.dat", prefix, ti, user[bi]._this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoad(user[bi].Ucat_sum,viewer);
			PetscViewerDestroy(&viewer);

			sprintf(filen, "%ssu1_%06d_%1d.dat", prefix, ti, user[bi]._this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoad(user[bi].Ucat_cross_sum,viewer);
			PetscViewerDestroy(&viewer);
			
			sprintf(filen, "%ssu2_%06d_%1d.dat", prefix, ti, user[bi]._this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoad(user[bi].Ucat_square_sum,viewer);
			PetscViewerDestroy(&viewer);
		}
			
		
		DMDAVecGetArray(user[bi].fda, user[bi].Ucat_cross_sum, &u1sum);
		DMDAVecGetArray(user[bi].fda, user[bi].Ucat_square_sum, &u2sum);
		
		DMDAVecGetArray(user[bi].da, UU, &uu);
		DMDAVecGetArray(user[bi].da, VV, &vv);
		DMDAVecGetArray(user[bi].da, WW, &ww);
		DMDAVecGetArray(user[bi].da, UV, &uv);
		DMDAVecGetArray(user[bi].da, VW, &vw);
		DMDAVecGetArray(user[bi].da, WU, &wu);
		DMDAVecGetArray(user[bi].da, UU_Mean, &uu_mean);
		DMDAVecGetArray(user[bi].da, VV_Mean, &vv_mean);
		DMDAVecGetArray(user[bi].da, WW_Mean, &ww_mean);

		DMDAVecGetArray(user[bi].fda, user[bi].Ucat_sum, &usum);


		// Reynolds stresses
		for (k=zs; k<ze-2; k++)
		for (j=ys; j<ye-2; j++)
		for (i=xs; i<xe-2; i++) {
			double _U = usum[k+1][j+1][i+1].x/N;
			double _V = usum[k+1][j+1][i+1].y/N;
			double _W = usum[k+1][j+1][i+1].z/N;
			double _UV = usum[k+1][j+1][i+1].x*usum[k+1][j+1][i+1].y / (N*N);
			double _VW = usum[k+1][j+1][i+1].y*usum[k+1][j+1][i+1].z / (N*N);
			double _WU = usum[k+1][j+1][i+1].z*usum[k+1][j+1][i+1].x / (N*N);
			
			uu[k+1][j+1][i+1] = ( u2sum[k+1][j+1][i+1].x/N - _U*_U );
			vv[k+1][j+1][i+1] = ( u2sum[k+1][j+1][i+1].y/N - _V*_V );
			ww[k+1][j+1][i+1] = ( u2sum[k+1][j+1][i+1].z/N - _W*_W );
			uv[k+1][j+1][i+1] = u1sum[k+1][j+1][i+1].x/N - _UV;
			vw[k+1][j+1][i+1] = u1sum[k+1][j+1][i+1].y/N - _VW;
			wu[k+1][j+1][i+1] = u1sum[k+1][j+1][i+1].z/N - _WU;
			uu_mean[k+1][j+1][i+1] = _U*_U;
			vv_mean[k+1][j+1][i+1] = _V*_V;
			ww_mean[k+1][j+1][i+1] = _W*_W;
		}
		DMDAVecRestoreArray(user[bi].fda, user[bi].Ucat_sum, &usum);
		
					
		// Mean velocity
		VecAXPY(Ucat, 1./N, user[bi].Ucat_sum);
		
		// Destroy
		VecDestroy(&user[bi].Ucat_sum);
		
		
		DMDAVecGetArray(user[bi].fda, Ucat, &ucat);
	
		double mke[mz], MeanConv[mz], PresTrans[mz], ViscTrans[mz], ViscDisp[mz], TurbTrans[mz], TurbProd[mz]; 
		
	
		/* MKE */
		for (k=zs; k<ze-2; k++) {
			double A_sum=0.0,Val_sum=0.0;
	
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) {

				double rr=sqrt(pow(coor[k+1][j+1][i+1].x-xc_tb,2)+pow(coor[k+1][j+1][i+1].y-yc_tb,2));
			       	double area=sqrt( zet[k+1][j+1][i+1].x*zet[k+1][j+1][i+1].x + zet[k+1][j+1][i+1].y*zet[k+1][j+1][i+1].y + zet[k+1][j+1][i+1].z*zet[k+1][j+1][i+1].z );
				double _U = ucat[k+1][j+1][i+1].x;
				double _V = ucat[k+1][j+1][i+1].y;
				double _W = ucat[k+1][j+1][i+1].z;

				x[k * (mx-2)*(my-2) + j*(mx-2) + i] = 0.5 * ( _U*_U+_V*_V+_W*_W );

				if (rr<r_tb) {
					A_sum += area;
					Val_sum += 0.5 * ( _U*_U+_V*_V+_W*_W )*area;
				}	
		
			}

			mke[k+1]=Val_sum/(A_sum+1.e-19);
		}

		if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		I = TECDAT100(&III, x, &DIsDouble);
	

	
		/* Turbulence Production */
		for (k=zs; k<ze-2; k++) {
			double A_sum=0.0,Val_sum=0.0;
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) {

				double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
				double eta0= eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
				double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
				double ajc = aj[k+1][j+1][i+1];
				double dudc, dude, dudz;
				double dvdc, dvde, dvdz;
				double dwdc, dwde, dwdz;
				double du_dx, du_dy, du_dz;
				double dv_dx, dv_dy, dv_dz;
				double dw_dx, dw_dy, dw_dz;
			
				Compute_du_center (i+1, j+1, k+1, mx, my, mz, ucat, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
				Compute_du_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
			
				double Prod = 0;
				Prod -= uu[k+1][j+1][i+1] * du_dx;
				Prod -= uv[k+1][j+1][i+1] * du_dy;
				Prod -= wu[k+1][j+1][i+1] * du_dz;
				Prod -= uv[k+1][j+1][i+1] * dv_dx;
				Prod -= vv[k+1][j+1][i+1] * dv_dy;
				Prod -= vw[k+1][j+1][i+1] * dv_dz;
				Prod -= wu[k+1][j+1][i+1] * dw_dx;
				Prod -= vw[k+1][j+1][i+1] * dw_dy;
				Prod -= ww[k+1][j+1][i+1] * dw_dz;

				x[k * (mx-2)*(my-2) + j*(mx-2) + i] = -Prod;


				double rr=sqrt(pow(coor[k+1][j+1][i+1].x-xc_tb,2)+pow(coor[k+1][j+1][i+1].y-yc_tb,2));
			       	double area=sqrt( zet[k+1][j+1][i+1].x*zet[k+1][j+1][i+1].x + zet[k+1][j+1][i+1].y*zet[k+1][j+1][i+1].y + zet[k+1][j+1][i+1].z*zet[k+1][j+1][i+1].z );

				if (rr<r_tb) {
					A_sum += area;
					Val_sum += -Prod*area;
				}	
	
			}

			TurbProd[k+1]=Val_sum/(A_sum+1.e-19);

		}

		if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		I = TECDAT100(&III, x, &DIsDouble);
	

	
		/* Dissipation */

		for (k=zs; k<ze-2; k++) {
			double A_sum=0.0, Val_sum=0.0;
	
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) {
				double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
				double eta0= eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
				double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
				double ajc = aj[k+1][j+1][i+1];
				double dudc, dude, dudz;
				double dvdc, dvde, dvdz;
				double dwdc, dwde, dwdz;
				double du_dx, du_dy, du_dz;
				double dv_dx, dv_dy, dv_dz;
				double dw_dx, dw_dy, dw_dz;
			
				Compute_du_center (i+1, j+1, k+1, mx, my, mz, ucat, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
				Compute_du_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
								
				double Diss = 0;
				Diss += pow(du_dx, 2.) + pow(du_dy, 2.) + pow(du_dz, 2.);
				Diss += pow(dv_dx, 2.) + pow(dv_dy, 2.) + pow(dv_dz, 2.);
				Diss += pow(dw_dx, 2.) + pow(dw_dy, 2.) + pow(dw_dz, 2.);
					
				Diss *= 1./user[bi].ren;
				Diss *= -1.;
				x[k * (mx-2)*(my-2) + j*(mx-2) + i] = Diss;


				double rr=sqrt(pow(coor[k+1][j+1][i+1].x-xc_tb,2)+pow(coor[k+1][j+1][i+1].y-yc_tb,2));
			       	double area=sqrt( zet[k+1][j+1][i+1].x*zet[k+1][j+1][i+1].x + zet[k+1][j+1][i+1].y*zet[k+1][j+1][i+1].y + zet[k+1][j+1][i+1].z*zet[k+1][j+1][i+1].z );

				if (rr<r_tb) {
					A_sum += area;
					Val_sum += Diss*area;
				}	
			}

			ViscDisp[k+1]=Val_sum/(A_sum+1.e-19);	

		}
		if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		I = TECDAT100(&III, x, &DIsDouble);
				

		
		/* Viscous Diffusion */
		
		// Make dk/dxj vector
		Vec GradK;
		Cmpnts ***gradk;
		
		DMCreateGlobalVector(user[bi].fda, &GradK);
		DMDAVecGetArray(user[bi].fda, GradK, &gradk);
		
		for (k=zs; k<ze-2; k++)
		for (j=ys; j<ye-2; j++)
		for (i=xs; i<xe-2; i++) {
			double _U = ucat[k+1][j+1][i+1].x;
			double _V = ucat[k+1][j+1][i+1].y;
			double _W = ucat[k+1][j+1][i+1].z;

			double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
			double eta0= eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
			double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
			double ajc = aj[k+1][j+1][i+1];
			double dudc, dude, dudz;
			double dvdc, dvde, dvdz;
			double dwdc, dwde, dwdz;
			double du_dx, du_dy, du_dz;
			double dv_dx, dv_dy, dv_dz;
			double dw_dx, dw_dy, dw_dz;
			double dk_dx, dk_dy, dk_dz;
		

			Compute_du_center (i+1, j+1, k+1, mx, my, mz, ucat, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
			Compute_du_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
	
			dk_dx = 0;
			dk_dx += _U*du_dx;
			dk_dx += _V*dv_dx;
			dk_dx += _W*dw_dx;
			
			dk_dy = 0;
			dk_dy += _U*du_dy;
			dk_dy += _V*dv_dy;
			dk_dy += _W*dw_dy;
			
			dk_dz = 0;
			dk_dz += _U*du_dz;
			dk_dz += _V*dv_dz;
			dk_dz += _W*dw_dz;
			
			gradk[k+1][j+1][i+1].x = dk_dx;
			gradk[k+1][j+1][i+1].y = dk_dy;
			gradk[k+1][j+1][i+1].z = dk_dz;
		}
		
		// Viscous Diffusion
		for (k=zs; k<ze-2; k++) {
			double A_sum=0.0, Val_sum=0.0;
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) {
				double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
				double eta0= eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
				double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
				double ajc = aj[k+1][j+1][i+1];
				double dudc, dude, dudz;
				double dvdc, dvde, dvdz;
				double dwdc, dwde, dwdz;
				double du_dx, du_dy, du_dz;
				double dv_dx, dv_dy, dv_dz;
				double dw_dx, dw_dy, dw_dz;
			
				Compute_du_center (i+1, j+1, k+1, mx, my, mz, gradk, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
				Compute_du_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
			
				double ddk_dxdx = du_dx;
				double ddk_dydy = dv_dy;
				double ddk_dzdz = dw_dz;
				double nu = 1./user[bi].ren;
			
				double ViscousDiffusion = 0;
			
				ViscousDiffusion = nu * (ddk_dxdx);
				ViscousDiffusion += nu * (ddk_dydy);
				ViscousDiffusion += nu * (ddk_dzdz);
			
				x[k * (mx-2)*(my-2) + j*(mx-2) + i] = ViscousDiffusion;



				double rr=sqrt(pow(coor[k+1][j+1][i+1].x-xc_tb,2)+pow(coor[k+1][j+1][i+1].y-yc_tb,2));
			       	double area=sqrt( zet[k+1][j+1][i+1].x*zet[k+1][j+1][i+1].x + zet[k+1][j+1][i+1].y*zet[k+1][j+1][i+1].y + zet[k+1][j+1][i+1].z*zet[k+1][j+1][i+1].z );

				if (rr<r_tb) {
					A_sum += area;
					Val_sum += ViscousDiffusion*area;
				}

			}
			
			ViscTrans[k+1]=Val_sum/(A_sum+1.e-19);
		}
		if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		I = TECDAT100(&III, x, &DIsDouble);
		
		DMDAVecRestoreArray(user[bi].fda, GradK, &gradk);
		VecDestroy(&GradK);
	
	
		/* Pressure Transport */
		Vec P_sum;
				
		PetscReal ***psum;
				
		DMCreateGlobalVector(user[bi].da, &P_sum);
				
		sprintf(filen, "%ssp_%06d_%1d.dat", prefix, ti, user[bi]._this);
		PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
		VecLoad(P_sum,viewer);
		PetscViewerDestroy(&viewer);
				
		DMDAVecGetArray(user[bi].da, P_sum, &psum);
				
		for (k=zs; k<ze-2; k++) {
			double A_sum=0.0, Val_sum=0.0;
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) {
				double dpdc, dpde, dpdz;
				double dp_dx, dp_dy, dp_dz;
				double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
				double eta0= eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
				double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
				double ajc = aj[k+1][j+1][i+1];
					
				Compute_dscalar_center (i+1, j+1, k+1, mx, my, mz, psum, nvert, &dpdc, &dpde, &dpdz);
				Compute_dscalar_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dpdc, dpde, dpdz, &dp_dx, &dp_dy, &dp_dz);

				dp_dx /= N;
				dp_dy /= N;
				dp_dz /= N;
			
				double U = ucat[k+1][j+1][i+1].x;
				double V = ucat[k+1][j+1][i+1].y;
				double W = ucat[k+1][j+1][i+1].z;
			
				double PressureTrans = 0;
				PressureTrans -= U * dp_dx;
				PressureTrans -= V * dp_dy;
				PressureTrans -= W * dp_dz;
					
				x[k * (mx-2)*(my-2) + j*(mx-2) + i] = PressureTrans;


				double rr=sqrt(pow(coor[k+1][j+1][i+1].x-xc_tb,2)+pow(coor[k+1][j+1][i+1].y-yc_tb,2));
			       	double area=sqrt( zet[k+1][j+1][i+1].x*zet[k+1][j+1][i+1].x + zet[k+1][j+1][i+1].y*zet[k+1][j+1][i+1].y + zet[k+1][j+1][i+1].z*zet[k+1][j+1][i+1].z );

				if (rr<r_tb) {
					A_sum += area;
					Val_sum += PressureTrans*area;
				}
			}
			PresTrans[k+1]=Val_sum/(A_sum+1.e-19);
		}
		if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		I = TECDAT100(&III, x, &DIsDouble);

		DMDAVecRestoreArray(user[bi].da, P_sum, &psum);
				
		VecDestroy(&P_sum);
	
	
		/* Turbulent Convection */
		Vec UUU_sum, UUU;
		Cmpnts ***uuu;
		
		DMCreateGlobalVector(user[bi].fda, &UUU);
		
		DMDAVecGetArray(user[bi].fda, UUU, &uuu);
		
		// a new vector for <ui'ui'uk'>
		for (k=zs; k<ze-2; k++)
		for (j=ys; j<ye-2; j++)
		for (i=xs; i<xe-2; i++) {
			double UU = uu[k+1][j+1][i+1];
			double VV = vv[k+1][j+1][i+1];
			double WW = ww[k+1][j+1][i+1];
			double UV = uv[k+1][j+1][i+1];
			double VW = vw[k+1][j+1][i+1];
			double WU = wu[k+1][j+1][i+1];
			
			double _U = ucat[k+1][j+1][i+1].x;
			double _V = ucat[k+1][j+1][i+1].y;
			double _W = ucat[k+1][j+1][i+1].z;
			
			uuu[k+1][j+1][i+1].x = _U*UU+_V*UV+_W*WU;
			uuu[k+1][j+1][i+1].y = _U*UV+_V*VV+_W*VW;
			uuu[k+1][j+1][i+1].z = _U*WU+_V*VW+_W*WW;
		}
		
		// calculate the divergence of <uuu>
		for (k=zs; k<ze-2; k++) {
			double A_sum=0.0, Val_sum=0.0;
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) {
				double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
				double eta0= eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
				double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
				double ajc = aj[k+1][j+1][i+1];
				double dudc, dude, dudz;
				double dvdc, dvde, dvdz;
				double dwdc, dwde, dwdz;
				double du_dx, du_dy, du_dz;
				double dv_dx, dv_dy, dv_dz;
				double dw_dx, dw_dy, dw_dz;
			
				Compute_du_center (i+1, j+1, k+1, mx, my, mz, uuu, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
				Compute_du_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
			
				double TurbulentConvection = du_dx + dv_dy + dw_dz;
				
				x[k * (mx-2)*(my-2) + j*(mx-2) + i] = -TurbulentConvection;


				double rr=sqrt(pow(coor[k+1][j+1][i+1].x-xc_tb,2)+pow(coor[k+1][j+1][i+1].y-yc_tb,2));
			       	double area=sqrt( zet[k+1][j+1][i+1].x*zet[k+1][j+1][i+1].x + zet[k+1][j+1][i+1].y*zet[k+1][j+1][i+1].y + zet[k+1][j+1][i+1].z*zet[k+1][j+1][i+1].z );

				if (rr<r_tb) {
					A_sum += area;
					Val_sum += -TurbulentConvection*area;
				}
			}
			TurbTrans[k+1] = Val_sum/(A_sum+1.e-19);
		}

		if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		I = TECDAT100(&III, x, &DIsDouble);
		
		//DMDAVecRestoreArray(user[bi].fda, UUU_sum, &uuusum);
		DMDAVecRestoreArray(user[bi].fda, UUU, &uuu);
		//VecDestroy(&UUU_sum);
		VecDestroy(&UUU);

		
		/* Mean Convection */
		for (k=zs; k<ze-2; k++) {
			double A_sum=0.0, Val_sum=0.0;
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) {
				double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
				double eta0= eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
				double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
				double ajc = aj[k+1][j+1][i+1];
				double duudc, duude, duudz;
				double dvvdc, dvvde, dvvdz;
				double dwwdc, dwwde, dwwdz;
				double duu_dx, duu_dy, duu_dz;
				double dvv_dx, dvv_dy, dvv_dz;
				double dww_dx, dww_dy, dww_dz;
				double dk_dx, dk_dy, dk_dz;
			
				Compute_dscalar_center (i+1, j+1, k+1, mx, my, mz, uu_mean, nvert, &duudc, &duude, &duudz);
				Compute_dscalar_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, duudc, duude, duudz, &duu_dx, &duu_dy, &duu_dz);
			
				Compute_dscalar_center (i+1, j+1, k+1, mx, my, mz, vv_mean, nvert, &dvvdc, &dvvde, &dvvdz);
				Compute_dscalar_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dvvdc, dvvde, dvvdz, &dvv_dx, &dvv_dy, &dvv_dz);
			
				Compute_dscalar_center (i+1, j+1, k+1, mx, my, mz, ww_mean, nvert, &dwwdc, &dwwde, &dwwdz);
				Compute_dscalar_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dwwdc, dwwde, dwwdz, &dww_dx, &dww_dy, &dww_dz);
			
				dk_dx = 0;
				dk_dx += 0.5 * duu_dx;
				dk_dx += 0.5 * dvv_dx;
				dk_dx += 0.5 * dww_dx;
			
				dk_dy = 0;
				dk_dy += 0.5 * duu_dy;
				dk_dy += 0.5 * dvv_dy;
				dk_dy += 0.5 * dww_dy;
			
				dk_dz = 0;
				dk_dz += 0.5 * duu_dz;
				dk_dz += 0.5 * dvv_dz;
				dk_dz += 0.5 * dww_dz;
			
				double MeanConvection = 0;
			
				MeanConvection += ucat[k][j][i].x * dk_dx;
				MeanConvection += ucat[k][j][i].y * dk_dy;
				MeanConvection += ucat[k][j][i].z * dk_dz;
				MeanConvection *= -1.0;

				x[k * (mx-2)*(my-2) + j*(mx-2) + i] = MeanConvection;

				double rr=sqrt(pow(coor[k+1][j+1][i+1].x-xc_tb,2)+pow(coor[k+1][j+1][i+1].y-yc_tb,2));
			       	double area=sqrt( zet[k+1][j+1][i+1].x*zet[k+1][j+1][i+1].x + zet[k+1][j+1][i+1].y*zet[k+1][j+1][i+1].y + zet[k+1][j+1][i+1].z*zet[k+1][j+1][i+1].z );

				if (rr<r_tb) {
					A_sum += area;
					Val_sum += MeanConvection*area;
				}
			}
		
			MeanConv[k+1]=Val_sum/(A_sum+1.e-19);

		}
		if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		I = TECDAT100(&III, x, &DIsDouble);
	

	
		FILE *f;
		char fname[128];
	   	sprintf(fname, "DiskAveraged_mkeBudget.plt");
		f = fopen(fname, "w");
	      	PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=\"z\", \"mke\", \"TurbulenceProduction\", \"TurbulenceConvection\", \"ViscousDissipation\", \"PressureTransport\", \"MeanConvection\", \"ViscousDiffusion\" \n");

		i=1;j=1;
		for (k=zs+1; k<ze-1; k++) {
      			PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le %le %le %le %le %le\n", coor[k][j][i].z, mke[k],TurbProd[k],TurbTrans[k],ViscDisp[k],PresTrans[k],MeanConv[k],ViscTrans[k]);
		}
	      	fclose(f);



	
		DMDAVecRestoreArray(fda, Coor, &coor);
		DMDAVecRestoreArray(user[bi].fda, Ucat, &ucat);
		
		DMDAVecRestoreArray(user[bi].da, UU, &uu);
		DMDAVecRestoreArray(user[bi].da, VV, &vv);
		DMDAVecRestoreArray(user[bi].da, WW, &ww);
		DMDAVecRestoreArray(user[bi].da, UV, &uv);
		DMDAVecRestoreArray(user[bi].da, VW, &vw);
		DMDAVecRestoreArray(user[bi].da, WU, &wu);
		DMDAVecRestoreArray(user[bi].da, UU_Mean, &uu_mean);
		DMDAVecRestoreArray(user[bi].da, VV_Mean, &vv_mean);
		DMDAVecRestoreArray(user[bi].da, WW_Mean, &ww_mean);
				
		DMDAVecRestoreArray(user[bi].fda, user[bi].Csi, &csi);
		DMDAVecRestoreArray(user[bi].fda, user[bi].Eta, &eta);
		DMDAVecRestoreArray(user[bi].fda, user[bi].Zet, &zet);
		DMDAVecRestoreArray(user[bi].da, user[bi].Aj, &aj);
		DMDAVecRestoreArray(user[bi].da, user[bi].Nvert, &nvert);
		
		DMDAVecRestoreArray(user[bi].fda, user[bi].Ucat_cross_sum, &u1sum);
		DMDAVecRestoreArray(user[bi].fda, user[bi].Ucat_square_sum, &u2sum);
		
		delete []x;
		
		VecDestroy(&user[bi].Ucat_cross_sum);
		VecDestroy(&user[bi].Ucat_square_sum);
		VecDestroy(&Ucat);
		VecDestroy(&UU);
		VecDestroy(&VV);
		VecDestroy(&WW);
		VecDestroy(&UV);
		VecDestroy(&VW);
		VecDestroy(&WU);
		VecDestroy(&UU_Mean);
		VecDestroy(&VV_Mean);
		VecDestroy(&WW_Mean);
		printf("MKE Budget 13 !\n");
	}
	
	I = TECEND100();
	
	return 0;
}



PetscErrorCode TECIOOut_zMomentum_Budget(UserCtx *user)
{
	PetscInt bi;

	char filen[80];
	
	if(phase) sprintf(filen, "%szMomentumBudget_Phase%03d-avg.plt", prefix, (int)phase);
	else sprintf(filen, "%szMomentumBudget%06d-avg.plt", prefix, (int)ti);

	printf("zMomentum Budget !\n");
	
	INTEGER4 I, Debug, VIsDouble, DIsDouble, III, IMax, JMax, KMax;
	VIsDouble = 0;
	DIsDouble = 0;
	Debug = 0;

	I = TECINI100((char*)"Averaging", "X Y Z MKE_z Convection_x Convection_y Convection_z dP Reynolds_x Reynolds_y Reynolds_z Viscous_x Viscous_y Viscous_z",  filen, (char*)".",  &Debug,  &VIsDouble);
	
	for (bi=0; bi<block_number; bi++) {
		DM da = user[bi].da, fda = user[bi].fda;
		DMDALocalInfo info = user[bi].info;

		PetscInt	xs = info.xs, xe = info.xs + info.xm;
		PetscInt  	ys = info.ys, ye = info.ys + info.ym;
		PetscInt	zs = info.zs, ze = info.zs + info.zm;
		PetscInt	mx = info.mx, my = info.my, mz = info.mz;
    
		PetscInt	lxs, lys, lzs, lxe, lye, lze;
		PetscInt	i, j, k;
		Cmpnts	***coor, ***u2sum, ***u1sum,  ***usum;
		Cmpnts	***ucat, ***csi, ***eta, ***zet; 
		PetscReal	***nvert, ***uu, ***vv, ***ww, ***uv, ***vw, ***wu, ***aj, ***uu_mean, ***vv_mean, ***ww_mean;
		Vec Coor;
		Vec Ucat, UU, VV, WW, UV, VW, WU, UU_Mean, VV_Mean, WW_Mean;	// define vectors

		INTEGER4	ZoneType=0, ICellMax=0, JCellMax=0, KCellMax=0;
		INTEGER4	IsBlock=1, NumFaceConnections=0, FaceNeighborMode=0;
		INTEGER4    ShareConnectivityFromZone=0;
		INTEGER4	LOC[100] = {1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; 
		/* 1 is cell-centered   0 is node centered */

		IMax = mx-1; JMax = my-1; KMax = mz-1;

		I = TECZNE100((char*)"Block 1",
			&ZoneType, 	/* Ordered zone */
			&IMax,
			&JMax,
			&KMax,
			&ICellMax,
			&JCellMax,
			&KCellMax,
			&IsBlock,	/* ISBLOCK  BLOCK format */
			&NumFaceConnections,
			&FaceNeighborMode,
			LOC,
			NULL,
			&ShareConnectivityFromZone); /* No connectivity sharing */

		float *x;
	
		x = new float [(mx-1)*(my-1)*(mz-1)];
		III = (mx-1) * (my-1) * (mz-1);

		DMDAVecGetArray(user[bi].fda, user[bi].Csi, &csi);
		DMDAVecGetArray(user[bi].fda, user[bi].Eta, &eta);
		DMDAVecGetArray(user[bi].fda, user[bi].Zet, &zet);
		DMDAVecGetArray(user[bi].da, user[bi].Aj, &aj);
		DMDAVecGetArray(user[bi].da, user[bi].Nvert, &nvert);
				
		DMDAGetCoordinates(da, &Coor);
		DMDAVecGetArray(fda, Coor, &coor);


//
		// X
		for (k=zs; k<ze-1; k++)
		for (j=ys; j<ye-1; j++)
		for (i=xs; i<xe-1; i++) x[k * (mx-1)*(my-1) + j*(mx-1) + i] = coor[k][j][i].x;
		I = TECDAT100(&III, x, &DIsDouble);

		// Y
		for (k=zs; k<ze-1; k++)
		for (j=ys; j<ye-1; j++)
		for (i=xs; i<xe-1; i++)	x[k * (mx-1)*(my-1) + j*(mx-1) + i] = coor[k][j][i].y;
		I = TECDAT100(&III, x, &DIsDouble);

		// Z
		for (k=zs; k<ze-1; k++)
		for (j=ys; j<ye-1; j++)
		for (i=xs; i<xe-1; i++)	x[k * (mx-1)*(my-1) + j*(mx-1) + i] = coor[k][j][i].z;
		I = TECDAT100(&III, x, &DIsDouble);
	
    
		double N=(double)tis+1.0;
		if(phase) N = phase;
	
		III = (mx-2) * (my-2) * (mz-2);

		PetscViewer viewer;
		char filen[128];
		
		DMCreateGlobalVector(user[bi].fda, &user[bi].Ucat_sum);
		DMCreateGlobalVector(user[bi].fda, &user[bi].Ucat_cross_sum);
		DMCreateGlobalVector(user[bi].fda, &user[bi].Ucat_square_sum);
		
		DMCreateGlobalVector(user[bi].fda, &Ucat);
		DMCreateGlobalVector(user[bi].da, &UU);
		DMCreateGlobalVector(user[bi].da, &VV);
		DMCreateGlobalVector(user[bi].da, &WW);
		DMCreateGlobalVector(user[bi].da, &UV);
		DMCreateGlobalVector(user[bi].da, &VW);
		DMCreateGlobalVector(user[bi].da, &WU);
		DMCreateGlobalVector(user[bi].da, &UU_Mean);
		DMCreateGlobalVector(user[bi].da, &VV_Mean);
		DMCreateGlobalVector(user[bi].da, &WW_Mean);
		/*
		if(phase) {
			sprintf(filen, "phase%03d_su0_%06d_%1d.dat", phase, ti_phase, user[bi]._this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoad(viewer, (user[bi].Ucat_sum));
			PetscViewerDestroy(&viewer);
				
			PetscPrintf(PETSC_COMM_WORLD, "Read %s\n", filen);

			sprintf(filen, "phase%03d_su1_%06d_%1d.dat", phase, ti_phase, user[bi]._this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoad(viewer, (user[bi].Ucat_cross_sum));
			PetscViewerDestroy(&viewer);
				
			PetscPrintf(PETSC_COMM_WORLD, "Read %s\n", filen);
			
			sprintf(filen, "phase%03d_su2_%06d_%1d.dat", phase, ti_phase, user[bi]._this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoad(viewer, (user[bi].Ucat_square_sum));
			PetscViewerDestroy(&viewer);
				
			PetscPrintf(PETSC_COMM_WORLD, "Read %s\n", filen);
		}
		else*/ 
		{
			sprintf(filen, "%ssu0_%06d_%1d.dat", prefix, ti, user[bi]._this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoad(user[bi].Ucat_sum,viewer);
			PetscViewerDestroy(&viewer);

			sprintf(filen, "%ssu1_%06d_%1d.dat", prefix, ti, user[bi]._this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoad(user[bi].Ucat_cross_sum,viewer);
			PetscViewerDestroy(&viewer);
			
			sprintf(filen, "%ssu2_%06d_%1d.dat", prefix, ti, user[bi]._this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoad(user[bi].Ucat_square_sum,viewer);
			PetscViewerDestroy(&viewer);
		}
			
		
		DMDAVecGetArray(user[bi].fda, user[bi].Ucat_cross_sum, &u1sum);
		DMDAVecGetArray(user[bi].fda, user[bi].Ucat_square_sum, &u2sum);
		
		DMDAVecGetArray(user[bi].da, UU, &uu);
		DMDAVecGetArray(user[bi].da, VV, &vv);
		DMDAVecGetArray(user[bi].da, WW, &ww);
		DMDAVecGetArray(user[bi].da, UV, &uv);
		DMDAVecGetArray(user[bi].da, VW, &vw);
		DMDAVecGetArray(user[bi].da, WU, &wu);
		DMDAVecGetArray(user[bi].da, UU_Mean, &uu_mean);
		DMDAVecGetArray(user[bi].da, VV_Mean, &vv_mean);
		DMDAVecGetArray(user[bi].da, WW_Mean, &ww_mean);

		DMDAVecGetArray(user[bi].fda, user[bi].Ucat_sum, &usum);


		// Reynolds stresses
		for (k=zs; k<ze-2; k++)
		for (j=ys; j<ye-2; j++)
		for (i=xs; i<xe-2; i++) {
			double _U = usum[k+1][j+1][i+1].x/N;
			double _V = usum[k+1][j+1][i+1].y/N;
			double _W = usum[k+1][j+1][i+1].z/N;
			double _UV = usum[k+1][j+1][i+1].x*usum[k+1][j+1][i+1].y / (N*N);
			double _VW = usum[k+1][j+1][i+1].y*usum[k+1][j+1][i+1].z / (N*N);
			double _WU = usum[k+1][j+1][i+1].z*usum[k+1][j+1][i+1].x / (N*N);
			
			uu[k+1][j+1][i+1] = ( u2sum[k+1][j+1][i+1].x/N - _U*_U );
			vv[k+1][j+1][i+1] = ( u2sum[k+1][j+1][i+1].y/N - _V*_V );
			ww[k+1][j+1][i+1] = ( u2sum[k+1][j+1][i+1].z/N - _W*_W );
			uv[k+1][j+1][i+1] = u1sum[k+1][j+1][i+1].x/N - _UV;
			vw[k+1][j+1][i+1] = u1sum[k+1][j+1][i+1].y/N - _VW;
			wu[k+1][j+1][i+1] = u1sum[k+1][j+1][i+1].z/N - _WU;
			uu_mean[k+1][j+1][i+1] = _U*_U;
			vv_mean[k+1][j+1][i+1] = _V*_V;
			ww_mean[k+1][j+1][i+1] = _W*_W;
		}
		DMDAVecRestoreArray(user[bi].fda, user[bi].Ucat_sum, &usum);
		
					
		// Mean velocity
		VecAXPY(Ucat, 1./N, user[bi].Ucat_sum);
		
		// Destroy
		VecDestroy(&user[bi].Ucat_sum);
		
		
		DMDAVecGetArray(user[bi].fda, Ucat, &ucat);
	
		double mke_z[mz], Convectionx_z[mz], Convectiony_z[mz], Convectionz_z[mz], dP_z[mz], Reynoldsx_z[mz], Reynoldsy_z[mz], Reynoldsz_z[mz], Viscousx_z[mz], Viscousy_z[mz], Viscousz_z[mz]; 
		
	
		/* MKE_z */
		for (k=zs; k<ze-2; k++) {
			double A_sum=0.0,Val_sum=0.0;
	
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) {

				double rr=sqrt(pow(coor[k+1][j+1][i+1].x-xc_tb,2)+pow(coor[k+1][j+1][i+1].y-yc_tb,2));
			       	double area=sqrt( zet[k+1][j+1][i+1].x*zet[k+1][j+1][i+1].x + zet[k+1][j+1][i+1].y*zet[k+1][j+1][i+1].y + zet[k+1][j+1][i+1].z*zet[k+1][j+1][i+1].z );
				double _W = ucat[k+1][j+1][i+1].z;

				x[k * (mx-2)*(my-2) + j*(mx-2) + i] = 0.5 * _W*_W ;

				if (rr<r_tb) {
					A_sum += area;
					Val_sum += 0.5 * _W*_W *area;
				}	
		
			}

			mke_z[k+1]=Val_sum/(A_sum+1.e-19);
		}

		if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		I = TECDAT100(&III, x, &DIsDouble);
	
	
		/* Convection x*/
		for (k=zs; k<ze-2; k++) {
			double A_sum=0.0, Val_sum=0.0;
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) {

				double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
				double eta0= eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
				double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
				double ajc = aj[k+1][j+1][i+1];
				double dudc, dude, dudz;
				double dvdc, dvde, dvdz;
				double dwdc, dwde, dwdz;
				double du_dx, du_dy, du_dz;
				double dv_dx, dv_dy, dv_dz;
				double dw_dx, dw_dy, dw_dz;
			
				Compute_du_center (i+1, j+1, k+1, mx, my, mz, ucat, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
				Compute_du_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
	
				double MeanConvection_x = 0;
			
				MeanConvection_x = -ucat[k][j][i].x * dw_dx;

				x[k * (mx-2)*(my-2) + j*(mx-2) + i] = MeanConvection_x;

				double rr=sqrt(pow(coor[k+1][j+1][i+1].x-xc_tb,2)+pow(coor[k+1][j+1][i+1].y-yc_tb,2));
			       	double area=sqrt( zet[k+1][j+1][i+1].x*zet[k+1][j+1][i+1].x + zet[k+1][j+1][i+1].y*zet[k+1][j+1][i+1].y + zet[k+1][j+1][i+1].z*zet[k+1][j+1][i+1].z );

				if (rr<r_tb) {
					A_sum += area;
					Val_sum += MeanConvection_x*area;
				}
			}
		
			Convectionx_z[k+1]=Val_sum/(A_sum+1.e-19);

		}
		if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		I = TECDAT100(&III, x, &DIsDouble);
	
	
		/* Convection y*/
		for (k=zs; k<ze-2; k++) {
			double A_sum=0.0, Val_sum=0.0;
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) {

				double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
				double eta0= eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
				double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
				double ajc = aj[k+1][j+1][i+1];
				double dudc, dude, dudz;
				double dvdc, dvde, dvdz;
				double dwdc, dwde, dwdz;
				double du_dx, du_dy, du_dz;
				double dv_dx, dv_dy, dv_dz;
				double dw_dx, dw_dy, dw_dz;
			
				Compute_du_center (i+1, j+1, k+1, mx, my, mz, ucat, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
				Compute_du_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
	
				double MeanConvection_y = 0;
			
				MeanConvection_y = -ucat[k][j][i].y * dw_dy;

				x[k * (mx-2)*(my-2) + j*(mx-2) + i] = MeanConvection_y;

				double rr=sqrt(pow(coor[k+1][j+1][i+1].x-xc_tb,2)+pow(coor[k+1][j+1][i+1].y-yc_tb,2));
			       	double area=sqrt( zet[k+1][j+1][i+1].x*zet[k+1][j+1][i+1].x + zet[k+1][j+1][i+1].y*zet[k+1][j+1][i+1].y + zet[k+1][j+1][i+1].z*zet[k+1][j+1][i+1].z );

				if (rr<r_tb) {
					A_sum += area;
					Val_sum += MeanConvection_y*area;
				}
			}
		
			Convectiony_z[k+1]=Val_sum/(A_sum+1.e-19);

		}
		if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		I = TECDAT100(&III, x, &DIsDouble);
	
	
		/* Convection z*/
		for (k=zs; k<ze-2; k++) {
			double A_sum=0.0, Val_sum=0.0;
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) {

				double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
				double eta0= eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
				double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
				double ajc = aj[k+1][j+1][i+1];
				double dudc, dude, dudz;
				double dvdc, dvde, dvdz;
				double dwdc, dwde, dwdz;
				double du_dx, du_dy, du_dz;
				double dv_dx, dv_dy, dv_dz;
				double dw_dx, dw_dy, dw_dz;
			
				Compute_du_center (i+1, j+1, k+1, mx, my, mz, ucat, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
				Compute_du_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
	
				double MeanConvection_z = 0;
			
				MeanConvection_z = -ucat[k][j][i].z * dw_dz;

				x[k * (mx-2)*(my-2) + j*(mx-2) + i] = MeanConvection_z;

				double rr=sqrt(pow(coor[k+1][j+1][i+1].x-xc_tb,2)+pow(coor[k+1][j+1][i+1].y-yc_tb,2));
			       	double area=sqrt( zet[k+1][j+1][i+1].x*zet[k+1][j+1][i+1].x + zet[k+1][j+1][i+1].y*zet[k+1][j+1][i+1].y + zet[k+1][j+1][i+1].z*zet[k+1][j+1][i+1].z );

				if (rr<r_tb) {
					A_sum += area;
					Val_sum += MeanConvection_z*area;
				}
			}
		
			Convectionz_z[k+1]=Val_sum/(A_sum+1.e-19);

		}
		if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		I = TECDAT100(&III, x, &DIsDouble);
	
	
		/* Pressure Transport */
		Vec P_sum;
				
		PetscReal ***psum;
				
		DMCreateGlobalVector(user[bi].da, &P_sum);
				
		sprintf(filen, "%ssp_%06d_%1d.dat", prefix, ti, user[bi]._this);
		PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
		VecLoad(P_sum,viewer);
		PetscViewerDestroy(&viewer);
				
		DMDAVecGetArray(user[bi].da, P_sum, &psum);
				
		for (k=zs; k<ze-2; k++) {
			double A_sum=0.0, Val_sum=0.0;
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) {
				double dpdc, dpde, dpdz;
				double dp_dx, dp_dy, dp_dz;
				double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
				double eta0= eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
				double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
				double ajc = aj[k+1][j+1][i+1];
					
				Compute_dscalar_center (i+1, j+1, k+1, mx, my, mz, psum, nvert, &dpdc, &dpde, &dpdz);
				Compute_dscalar_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dpdc, dpde, dpdz, &dp_dx, &dp_dy, &dp_dz);

				dp_dx /= N;
				dp_dy /= N;
				dp_dz /= N;
			
			
				double PressureGradient_z = 0;
				PressureGradient_z = -dp_dz;
					
				x[k * (mx-2)*(my-2) + j*(mx-2) + i] = PressureGradient_z;

				double rr=sqrt(pow(coor[k+1][j+1][i+1].x-xc_tb,2)+pow(coor[k+1][j+1][i+1].y-yc_tb,2));
			       	double area=sqrt( zet[k+1][j+1][i+1].x*zet[k+1][j+1][i+1].x + zet[k+1][j+1][i+1].y*zet[k+1][j+1][i+1].y + zet[k+1][j+1][i+1].z*zet[k+1][j+1][i+1].z );

				if (rr<r_tb) {
					A_sum += area;
					Val_sum += PressureGradient_z*area;
				}
			}
			dP_z[k+1]=Val_sum/(A_sum+1.e-19);
		}
		if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		I = TECDAT100(&III, x, &DIsDouble);

		DMDAVecRestoreArray(user[bi].da, P_sum, &psum);
				
		VecDestroy(&P_sum);
	
		
		/* Reynolds Stress Gradient x  */
		
		
		// calculate the divergence of <wu>
		for (k=zs; k<ze-2; k++) {
			double A_sum=0.0, Val_sum=0.0;
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) {

				double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
				double eta0= eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
				double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
				double ajc = aj[k+1][j+1][i+1];
				double dwudc, dwude, dwudz;
				double dwu_dx, dwu_dy, dwu_dz;
			
				Compute_dscalar_center (i+1, j+1, k+1, mx, my, mz, wu, nvert, &dwudc, &dwude, &dwudz);
				Compute_dscalar_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dwudc, dwude, dwudz, &dwu_dx, &dwu_dy, &dwu_dz);
				
				x[k * (mx-2)*(my-2) + j*(mx-2) + i] = -dwu_dx;

				double rr=sqrt(pow(coor[k+1][j+1][i+1].x-xc_tb,2)+pow(coor[k+1][j+1][i+1].y-yc_tb,2));
			       	double area=sqrt( zet[k+1][j+1][i+1].x*zet[k+1][j+1][i+1].x + zet[k+1][j+1][i+1].y*zet[k+1][j+1][i+1].y + zet[k+1][j+1][i+1].z*zet[k+1][j+1][i+1].z );

				if (rr<r_tb) {
					A_sum += area;
					Val_sum += -dwu_dx*area;
				}
			}
			Reynoldsx_z[k+1] = Val_sum/(A_sum+1.e-19);
		}

		if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		I = TECDAT100(&III, x, &DIsDouble);
		

		
		/* Reynolds Stress Gradient y  */
		
		
		// calculate the divergence of <vw>
		for (k=zs; k<ze-2; k++) {
			double A_sum=0.0, Val_sum=0.0;
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) {

				double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
				double eta0= eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
				double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
				double ajc = aj[k+1][j+1][i+1];
				double dvwdc, dvwde, dvwdz;
				double dvw_dx, dvw_dy, dvw_dz;
			
				Compute_dscalar_center (i+1, j+1, k+1, mx, my, mz, vw, nvert, &dvwdc, &dvwde, &dvwdz);
				Compute_dscalar_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dvwdc, dvwde, dvwdz, &dvw_dx, &dvw_dy, &dvw_dz);
				
				x[k * (mx-2)*(my-2) + j*(mx-2) + i] = -dvw_dy;

				double rr=sqrt(pow(coor[k+1][j+1][i+1].x-xc_tb,2)+pow(coor[k+1][j+1][i+1].y-yc_tb,2));
			       	double area=sqrt( zet[k+1][j+1][i+1].x*zet[k+1][j+1][i+1].x + zet[k+1][j+1][i+1].y*zet[k+1][j+1][i+1].y + zet[k+1][j+1][i+1].z*zet[k+1][j+1][i+1].z );

				if (rr<r_tb) {
					A_sum += area;
					Val_sum += -dvw_dy*area;
				}
			}
			Reynoldsy_z[k+1] = Val_sum/(A_sum+1.e-19);
		}

		if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		I = TECDAT100(&III, x, &DIsDouble);
			
		/* Reynolds Stress Gradient z  */
		
		
		// calculate the divergence of <ww>
		for (k=zs; k<ze-2; k++) {
			double A_sum=0.0, Val_sum=0.0;
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) {

				double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
				double eta0= eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
				double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
				double ajc = aj[k+1][j+1][i+1];
				double dwwdc, dwwde, dwwdz;
				double dww_dx, dww_dy, dww_dz;
			
				Compute_dscalar_center (i+1, j+1, k+1, mx, my, mz, vw, nvert, &dwwdc, &dwwde, &dwwdz);
				Compute_dscalar_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dwwdc, dwwde, dwwdz, &dww_dx, &dww_dy, &dww_dz);
				
				x[k * (mx-2)*(my-2) + j*(mx-2) + i] = -dww_dz;

				double rr=sqrt(pow(coor[k+1][j+1][i+1].x-xc_tb,2)+pow(coor[k+1][j+1][i+1].y-yc_tb,2));
			       	double area=sqrt( zet[k+1][j+1][i+1].x*zet[k+1][j+1][i+1].x + zet[k+1][j+1][i+1].y*zet[k+1][j+1][i+1].y + zet[k+1][j+1][i+1].z*zet[k+1][j+1][i+1].z );

				if (rr<r_tb) {
					A_sum += area;
					Val_sum += -dww_dz*area;
				}
			}
			Reynoldsz_z[k+1] = Val_sum/(A_sum+1.e-19);
		}

		if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		I = TECDAT100(&III, x, &DIsDouble);
		

	
		/* Viscous Diffusion */
		
		// Make dk/dxj vector
		Vec GradK;
		Cmpnts ***gradk;
		
		DMCreateGlobalVector(user[bi].fda, &GradK);
		DMDAVecGetArray(user[bi].fda, GradK, &gradk);
		
		for (k=zs; k<ze-2; k++)
		for (j=ys; j<ye-2; j++)
		for (i=xs; i<xe-2; i++) {

			double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
			double eta0= eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
			double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
			double ajc = aj[k+1][j+1][i+1];
			double dudc, dude, dudz;
			double dvdc, dvde, dvdz;
			double dwdc, dwde, dwdz;
			double du_dx, du_dy, du_dz;
			double dv_dx, dv_dy, dv_dz;
			double dw_dx, dw_dy, dw_dz;
			double dk_dx, dk_dy, dk_dz;
		

			Compute_du_center (i+1, j+1, k+1, mx, my, mz, ucat, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
			Compute_du_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
	
			dk_dx = 0;
			dk_dx += dw_dx;
			
			dk_dy = 0;
			dk_dy += dw_dy;
			
			dk_dz = 0;
			dk_dz += dw_dz;
			
			gradk[k+1][j+1][i+1].x = dk_dx;
			gradk[k+1][j+1][i+1].y = dk_dy;
			gradk[k+1][j+1][i+1].z = dk_dz;
		}
		
		// Viscous dissipation x 
		for (k=zs; k<ze-2; k++) {
			double A_sum=0.0, Val_sum=0.0;
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) {
				double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
				double eta0= eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
				double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
				double ajc = aj[k+1][j+1][i+1];
				double dudc, dude, dudz;
				double dvdc, dvde, dvdz;
				double dwdc, dwde, dwdz;
				double du_dx, du_dy, du_dz;
				double dv_dx, dv_dy, dv_dz;
				double dw_dx, dw_dy, dw_dz;
			
				Compute_du_center (i+1, j+1, k+1, mx, my, mz, gradk, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
				Compute_du_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
			
				double nu = 1./user[bi].ren;
			
				x[k * (mx-2)*(my-2) + j*(mx-2) + i] = nu*du_dx;

				double rr=sqrt(pow(coor[k+1][j+1][i+1].x-xc_tb,2)+pow(coor[k+1][j+1][i+1].y-yc_tb,2));
			       	double area=sqrt( zet[k+1][j+1][i+1].x*zet[k+1][j+1][i+1].x + zet[k+1][j+1][i+1].y*zet[k+1][j+1][i+1].y + zet[k+1][j+1][i+1].z*zet[k+1][j+1][i+1].z );

				if (rr<r_tb) {
					A_sum += area;
					Val_sum += nu*du_dx*area;
				}

			}
			
			Viscousx_z[k+1]=Val_sum/(A_sum+1.e-19);
		}
		if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		I = TECDAT100(&III, x, &DIsDouble);
	
		// Viscous dissipation y 
		for (k=zs; k<ze-2; k++) {
			double A_sum=0.0, Val_sum=0.0;
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) {
				double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
				double eta0= eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
				double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
				double ajc = aj[k+1][j+1][i+1];
				double dudc, dude, dudz;
				double dvdc, dvde, dvdz;
				double dwdc, dwde, dwdz;
				double du_dx, du_dy, du_dz;
				double dv_dx, dv_dy, dv_dz;
				double dw_dx, dw_dy, dw_dz;
			
				Compute_du_center (i+1, j+1, k+1, mx, my, mz, gradk, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
				Compute_du_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
			
				double nu = 1./user[bi].ren;
			
				x[k * (mx-2)*(my-2) + j*(mx-2) + i] = nu*dv_dy;

				double rr=sqrt(pow(coor[k+1][j+1][i+1].x-xc_tb,2)+pow(coor[k+1][j+1][i+1].y-yc_tb,2));
			       	double area=sqrt( zet[k+1][j+1][i+1].x*zet[k+1][j+1][i+1].x + zet[k+1][j+1][i+1].y*zet[k+1][j+1][i+1].y + zet[k+1][j+1][i+1].z*zet[k+1][j+1][i+1].z );

				if (rr<r_tb) {
					A_sum += area;
					Val_sum += nu*dv_dy*area;
				}

			}
			
			Viscousy_z[k+1]=Val_sum/(A_sum+1.e-19);
		}
		if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		I = TECDAT100(&III, x, &DIsDouble);


		// Viscous dissipation z 
		for (k=zs; k<ze-2; k++) {
			double A_sum=0.0, Val_sum=0.0;
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) {
				double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
				double eta0= eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
				double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
				double ajc = aj[k+1][j+1][i+1];
				double dudc, dude, dudz;
				double dvdc, dvde, dvdz;
				double dwdc, dwde, dwdz;
				double du_dx, du_dy, du_dz;
				double dv_dx, dv_dy, dv_dz;
				double dw_dx, dw_dy, dw_dz;
			
				Compute_du_center (i+1, j+1, k+1, mx, my, mz, gradk, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
				Compute_du_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
			
				double nu = 1./user[bi].ren;
			
				x[k * (mx-2)*(my-2) + j*(mx-2) + i] = nu*dw_dz;

				double rr=sqrt(pow(coor[k+1][j+1][i+1].x-xc_tb,2)+pow(coor[k+1][j+1][i+1].y-yc_tb,2));
			       	double area=sqrt( zet[k+1][j+1][i+1].x*zet[k+1][j+1][i+1].x + zet[k+1][j+1][i+1].y*zet[k+1][j+1][i+1].y + zet[k+1][j+1][i+1].z*zet[k+1][j+1][i+1].z );

				if (rr<r_tb) {
					A_sum += area;
					Val_sum += nu*dw_dz*area;
				}

			}
			
			Viscousz_z[k+1]=Val_sum/(A_sum+1.e-19);
		}
		if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		I = TECDAT100(&III, x, &DIsDouble);


		DMDAVecRestoreArray(user[bi].fda, GradK, &gradk);
		VecDestroy(&GradK);
	

	
		FILE *f;
		char fname[128];
	   	sprintf(fname, "DiskAveraged_z-momentum.plt");
		f = fopen(fname, "w");
	      	PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=\"z\", \"z-mke\", \"Convectionx\", \"Convectiony\", \"Convectionz\", \"dP\", \"Reynoldsx\", \"Reynoldsy\", \"Reynoldsz\", \"Viscousx\", \"Viscousy\", \"Viscousz\" \n");

		i=1;j=1;
		for (k=zs+1; k<ze-1; k++) {
      			PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le %le %le %le %le %le %le %le %le %le \n", coor[k][j][i].z, mke_z[k], Convectionx_z[k], Convectiony_z[k], Convectionz_z[k], dP_z[k], Reynoldsx_z[k], Reynoldsy_z[k], Reynoldsz_z[k], Viscousx_z[k], Viscousy_z[k], Viscousz_z[k]);
		}
	      	fclose(f);



	
		DMDAVecRestoreArray(fda, Coor, &coor);
		DMDAVecRestoreArray(user[bi].fda, Ucat, &ucat);
		
		DMDAVecRestoreArray(user[bi].da, UU, &uu);
		DMDAVecRestoreArray(user[bi].da, VV, &vv);
		DMDAVecRestoreArray(user[bi].da, WW, &ww);
		DMDAVecRestoreArray(user[bi].da, UV, &uv);
		DMDAVecRestoreArray(user[bi].da, VW, &vw);
		DMDAVecRestoreArray(user[bi].da, WU, &wu);
		DMDAVecRestoreArray(user[bi].da, UU_Mean, &uu_mean);
		DMDAVecRestoreArray(user[bi].da, VV_Mean, &vv_mean);
		DMDAVecRestoreArray(user[bi].da, WW_Mean, &ww_mean);
				
		DMDAVecRestoreArray(user[bi].fda, user[bi].Csi, &csi);
		DMDAVecRestoreArray(user[bi].fda, user[bi].Eta, &eta);
		DMDAVecRestoreArray(user[bi].fda, user[bi].Zet, &zet);
		DMDAVecRestoreArray(user[bi].da, user[bi].Aj, &aj);
		DMDAVecRestoreArray(user[bi].da, user[bi].Nvert, &nvert);
		
		DMDAVecRestoreArray(user[bi].fda, user[bi].Ucat_cross_sum, &u1sum);
		DMDAVecRestoreArray(user[bi].fda, user[bi].Ucat_square_sum, &u2sum);
		
		delete []x;
		
		VecDestroy(&user[bi].Ucat_cross_sum);
		VecDestroy(&user[bi].Ucat_square_sum);
		VecDestroy(&Ucat);
		VecDestroy(&UU);
		VecDestroy(&VV);
		VecDestroy(&WW);
		VecDestroy(&UV);
		VecDestroy(&VW);
		VecDestroy(&WU);
		VecDestroy(&UU_Mean);
		VecDestroy(&VV_Mean);
		VecDestroy(&WW_Mean);
		printf("MKE Budget 13 !\n");
	}
	
	I = TECEND100();
	
	return 0;
}



PetscErrorCode TECIOOut_MKE_Budget1(UserCtx *user)
{
	PetscInt bi;

	char filen[80];
	
	if(phase) sprintf(filen, "%sMKEBudget_Phase%03d-avg.plt", prefix, (int)phase);
	else sprintf(filen, "%sMKEBudget%06d-avg.plt", prefix, (int)ti);

	printf("MKE Budget !\n");
	
	INTEGER4 I, Debug, VIsDouble, DIsDouble, III, IMax, JMax, KMax;
	VIsDouble = 0;
	DIsDouble = 0;
	Debug = 0;

	I = TECINI100((char*)"Averaging", "X Y Z MKE TurbulenceProduction Dissipation ViscousDiffusion PressureTransport TurbulentConvection TurbulentConvection_x TurbulentConvection_y TurbulentConvetion_z MeanConvection dW3dz/3",  filen, (char*)".",  &Debug,  &VIsDouble);
	
	for (bi=0; bi<block_number; bi++) {
		DM da = user[bi].da, fda = user[bi].fda;
		DMDALocalInfo info = user[bi].info;

		PetscInt	xs = info.xs, xe = info.xs + info.xm;
		PetscInt  	ys = info.ys, ye = info.ys + info.ym;
		PetscInt	zs = info.zs, ze = info.zs + info.zm;
		PetscInt	mx = info.mx, my = info.my, mz = info.mz;
    
		PetscInt	lxs, lys, lzs, lxe, lye, lze;
		PetscInt	i, j, k;
		Cmpnts	***coor, ***u2sum, ***u1sum,  ***usum;
		Cmpnts	***ucat, ***csi, ***eta, ***zet; 
		PetscReal	***nvert, ***uu, ***vv, ***ww, ***uv, ***vw, ***wu, ***aj, ***uu_mean, ***vv_mean, ***ww_mean;
		Vec Coor;
		Vec Ucat, UU, VV, WW, UV, VW, WU, UU_Mean, VV_Mean, WW_Mean;	// define vectors

		INTEGER4	ZoneType=0, ICellMax=0, JCellMax=0, KCellMax=0;
		INTEGER4	IsBlock=1, NumFaceConnections=0, FaceNeighborMode=0;
		INTEGER4    ShareConnectivityFromZone=0;
		INTEGER4	LOC[100] = {1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; 
		/* 1 is cell-centered   0 is node centered */

		IMax = mx-1; JMax = my-1; KMax = mz-1;

		I = TECZNE100((char*)"Block 1",
			&ZoneType, 	/* Ordered zone */
			&IMax,
			&JMax,
			&KMax,
			&ICellMax,
			&JCellMax,
			&KCellMax,
			&IsBlock,	/* ISBLOCK  BLOCK format */
			&NumFaceConnections,
			&FaceNeighborMode,
			LOC,
			NULL,
			&ShareConnectivityFromZone); /* No connectivity sharing */

		float *x;
	
		x = new float [(mx-1)*(my-1)*(mz-1)];
		III = (mx-1) * (my-1) * (mz-1);

		DMDAVecGetArray(user[bi].fda, user[bi].Csi, &csi);
		DMDAVecGetArray(user[bi].fda, user[bi].Eta, &eta);
		DMDAVecGetArray(user[bi].fda, user[bi].Zet, &zet);
		DMDAVecGetArray(user[bi].da, user[bi].Aj, &aj);
		DMDAVecGetArray(user[bi].da, user[bi].Nvert, &nvert);
				
		DMDAGetCoordinates(da, &Coor);
		DMDAVecGetArray(fda, Coor, &coor);


//
		// X
		for (k=zs; k<ze-1; k++)
		for (j=ys; j<ye-1; j++)
		for (i=xs; i<xe-1; i++) x[k * (mx-1)*(my-1) + j*(mx-1) + i] = coor[k][j][i].x;
		I = TECDAT100(&III, x, &DIsDouble);

		// Y
		for (k=zs; k<ze-1; k++)
		for (j=ys; j<ye-1; j++)
		for (i=xs; i<xe-1; i++)	x[k * (mx-1)*(my-1) + j*(mx-1) + i] = coor[k][j][i].y;
		I = TECDAT100(&III, x, &DIsDouble);

		// Z
		for (k=zs; k<ze-1; k++)
		for (j=ys; j<ye-1; j++)
		for (i=xs; i<xe-1; i++)	x[k * (mx-1)*(my-1) + j*(mx-1) + i] = coor[k][j][i].z;
		I = TECDAT100(&III, x, &DIsDouble);
	
    
		double N=(double)tis+1.0;
		if(phase) N = phase;
	
		III = (mx-2) * (my-2) * (mz-2);

		PetscViewer viewer;
		char filen[128];
		
		DMCreateGlobalVector(user[bi].fda, &user[bi].Ucat_sum);
		DMCreateGlobalVector(user[bi].fda, &user[bi].Ucat_cross_sum);
		DMCreateGlobalVector(user[bi].fda, &user[bi].Ucat_square_sum);
		
		DMCreateGlobalVector(user[bi].fda, &Ucat);
		DMCreateGlobalVector(user[bi].da, &UU);
		DMCreateGlobalVector(user[bi].da, &VV);
		DMCreateGlobalVector(user[bi].da, &WW);
		DMCreateGlobalVector(user[bi].da, &UV);
		DMCreateGlobalVector(user[bi].da, &VW);
		DMCreateGlobalVector(user[bi].da, &WU);
		DMCreateGlobalVector(user[bi].da, &UU_Mean);
		DMCreateGlobalVector(user[bi].da, &VV_Mean);
		DMCreateGlobalVector(user[bi].da, &WW_Mean);
		/*
		if(phase) {
			sprintf(filen, "phase%03d_su0_%06d_%1d.dat", phase, ti_phase, user[bi]._this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoad(viewer, (user[bi].Ucat_sum));
			PetscViewerDestroy(&viewer);
				
			PetscPrintf(PETSC_COMM_WORLD, "Read %s\n", filen);

			sprintf(filen, "phase%03d_su1_%06d_%1d.dat", phase, ti_phase, user[bi]._this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoad(viewer, (user[bi].Ucat_cross_sum));
			PetscViewerDestroy(&viewer);
				
			PetscPrintf(PETSC_COMM_WORLD, "Read %s\n", filen);
			
			sprintf(filen, "phase%03d_su2_%06d_%1d.dat", phase, ti_phase, user[bi]._this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoad(viewer, (user[bi].Ucat_square_sum));
			PetscViewerDestroy(&viewer);
				
			PetscPrintf(PETSC_COMM_WORLD, "Read %s\n", filen);
		}
		else*/ 
		{
			sprintf(filen, "%ssu0_%06d_%1d.dat", prefix, ti, user[bi]._this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoad(user[bi].Ucat_sum,viewer);
			PetscViewerDestroy(&viewer);

			sprintf(filen, "%ssu1_%06d_%1d.dat", prefix, ti, user[bi]._this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoad(user[bi].Ucat_cross_sum,viewer);
			PetscViewerDestroy(&viewer);
			
			sprintf(filen, "%ssu2_%06d_%1d.dat", prefix, ti, user[bi]._this);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
			VecLoad(user[bi].Ucat_square_sum,viewer);
			PetscViewerDestroy(&viewer);
		}
			
		
		DMDAVecGetArray(user[bi].fda, user[bi].Ucat_cross_sum, &u1sum);
		DMDAVecGetArray(user[bi].fda, user[bi].Ucat_square_sum, &u2sum);
		
		DMDAVecGetArray(user[bi].da, UU, &uu);
		DMDAVecGetArray(user[bi].da, VV, &vv);
		DMDAVecGetArray(user[bi].da, WW, &ww);
		DMDAVecGetArray(user[bi].da, UV, &uv);
		DMDAVecGetArray(user[bi].da, VW, &vw);
		DMDAVecGetArray(user[bi].da, WU, &wu);
		DMDAVecGetArray(user[bi].da, UU_Mean, &uu_mean);
		DMDAVecGetArray(user[bi].da, VV_Mean, &vv_mean);
		DMDAVecGetArray(user[bi].da, WW_Mean, &ww_mean);

		DMDAVecGetArray(user[bi].fda, user[bi].Ucat_sum, &usum);


		// Reynolds stresses
		for (k=zs; k<ze-2; k++)
		for (j=ys; j<ye-2; j++)
		for (i=xs; i<xe-2; i++) {
			double _U = usum[k+1][j+1][i+1].x/N;
			double _V = usum[k+1][j+1][i+1].y/N;
			double _W = usum[k+1][j+1][i+1].z/N;
			double _UV = usum[k+1][j+1][i+1].x*usum[k+1][j+1][i+1].y / (N*N);
			double _VW = usum[k+1][j+1][i+1].y*usum[k+1][j+1][i+1].z / (N*N);
			double _WU = usum[k+1][j+1][i+1].z*usum[k+1][j+1][i+1].x / (N*N);
			
			uu[k+1][j+1][i+1] = ( u2sum[k+1][j+1][i+1].x/N - _U*_U );
			vv[k+1][j+1][i+1] = ( u2sum[k+1][j+1][i+1].y/N - _V*_V );
			ww[k+1][j+1][i+1] = ( u2sum[k+1][j+1][i+1].z/N - _W*_W );
			uv[k+1][j+1][i+1] = u1sum[k+1][j+1][i+1].x/N - _UV;
			vw[k+1][j+1][i+1] = u1sum[k+1][j+1][i+1].y/N - _VW;
			wu[k+1][j+1][i+1] = u1sum[k+1][j+1][i+1].z/N - _WU;
			uu_mean[k+1][j+1][i+1] = _U*_U;
			vv_mean[k+1][j+1][i+1] = _V*_V;
			ww_mean[k+1][j+1][i+1] = _W*_W;
		}
		DMDAVecRestoreArray(user[bi].fda, user[bi].Ucat_sum, &usum);
		
					
		// Mean velocity
		VecAXPY(Ucat, 1./N, user[bi].Ucat_sum);
		
		// Destroy
		VecDestroy(&user[bi].Ucat_sum);
		
		
		DMDAVecGetArray(user[bi].fda, Ucat, &ucat);
	
		double mke[mz], MeanConv[mz], PresTrans[mz], ViscTrans[mz], ViscDisp[mz], TurbTrans[mz], TurbProd[mz], dW3dz[mz], TurbTrans_x[mz], TurbTrans_y[mz], TurbTrans_z[mz]; 
		
	
		/* MKE */
		for (k=zs; k<ze-2; k++) {
			double A_sum=0.0,Val_sum=0.0;
	
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) {

				double rr=sqrt(pow(coor[k+1][j+1][i+1].x-xc_tb,2)+pow(coor[k+1][j+1][i+1].y-yc_tb,2));
			       	double area=sqrt( zet[k+1][j+1][i+1].x*zet[k+1][j+1][i+1].x + zet[k+1][j+1][i+1].y*zet[k+1][j+1][i+1].y + zet[k+1][j+1][i+1].z*zet[k+1][j+1][i+1].z );
				double _U = ucat[k+1][j+1][i+1].x;
				double _V = ucat[k+1][j+1][i+1].y;
				double _W = ucat[k+1][j+1][i+1].z;

				x[k * (mx-2)*(my-2) + j*(mx-2) + i] = 0.5 * ( _U*_U+_V*_V+_W*_W );

				if (rr<r_tb) {
					A_sum += area;
					Val_sum += 0.5 * ( _U*_U+_V*_V+_W*_W )*area;
				}	
		
			}

			mke[k+1]=Val_sum/(A_sum+1.e-19);
		}

		if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		I = TECDAT100(&III, x, &DIsDouble);
	

	
		/* Turbulence Production */
		for (k=zs; k<ze-2; k++) {
			double A_sum=0.0,Val_sum=0.0;
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) {

				double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
				double eta0= eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
				double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
				double ajc = aj[k+1][j+1][i+1];
				double dudc, dude, dudz;
				double dvdc, dvde, dvdz;
				double dwdc, dwde, dwdz;
				double du_dx, du_dy, du_dz;
				double dv_dx, dv_dy, dv_dz;
				double dw_dx, dw_dy, dw_dz;
			
				Compute_du_center (i+1, j+1, k+1, mx, my, mz, ucat, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
				Compute_du_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
			
				double Prod = 0;
				Prod -= uu[k+1][j+1][i+1] * du_dx;
				Prod -= uv[k+1][j+1][i+1] * du_dy;
				Prod -= wu[k+1][j+1][i+1] * du_dz;
				Prod -= uv[k+1][j+1][i+1] * dv_dx;
				Prod -= vv[k+1][j+1][i+1] * dv_dy;
				Prod -= vw[k+1][j+1][i+1] * dv_dz;
				Prod -= wu[k+1][j+1][i+1] * dw_dx;
				Prod -= vw[k+1][j+1][i+1] * dw_dy;
				Prod -= ww[k+1][j+1][i+1] * dw_dz;

				x[k * (mx-2)*(my-2) + j*(mx-2) + i] = -Prod;


				double rr=sqrt(pow(coor[k+1][j+1][i+1].x-xc_tb,2)+pow(coor[k+1][j+1][i+1].y-yc_tb,2));
			       	double area=sqrt( zet[k+1][j+1][i+1].x*zet[k+1][j+1][i+1].x + zet[k+1][j+1][i+1].y*zet[k+1][j+1][i+1].y + zet[k+1][j+1][i+1].z*zet[k+1][j+1][i+1].z );

				if (rr<r_tb) {
					A_sum += area;
					Val_sum += -Prod*area;
				}	
	
			}

			TurbProd[k+1]=Val_sum/(A_sum+1.e-19);

		}

		if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		I = TECDAT100(&III, x, &DIsDouble);
	

	
		/* Dissipation */

		for (k=zs; k<ze-2; k++) {
			double A_sum=0.0, Val_sum=0.0;
	
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) {
				double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
				double eta0= eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
				double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
				double ajc = aj[k+1][j+1][i+1];
				double dudc, dude, dudz;
				double dvdc, dvde, dvdz;
				double dwdc, dwde, dwdz;
				double du_dx, du_dy, du_dz;
				double dv_dx, dv_dy, dv_dz;
				double dw_dx, dw_dy, dw_dz;
			
				Compute_du_center (i+1, j+1, k+1, mx, my, mz, ucat, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
				Compute_du_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
								
				double Diss = 0;
				Diss += pow(du_dx, 2.) + pow(du_dy, 2.) + pow(du_dz, 2.);
				Diss += pow(dv_dx, 2.) + pow(dv_dy, 2.) + pow(dv_dz, 2.);
				Diss += pow(dw_dx, 2.) + pow(dw_dy, 2.) + pow(dw_dz, 2.);
					
				Diss *= 1./user[bi].ren;
				Diss *= -1.;
				x[k * (mx-2)*(my-2) + j*(mx-2) + i] = Diss;


				double rr=sqrt(pow(coor[k+1][j+1][i+1].x-xc_tb,2)+pow(coor[k+1][j+1][i+1].y-yc_tb,2));
			       	double area=sqrt( zet[k+1][j+1][i+1].x*zet[k+1][j+1][i+1].x + zet[k+1][j+1][i+1].y*zet[k+1][j+1][i+1].y + zet[k+1][j+1][i+1].z*zet[k+1][j+1][i+1].z );

				if (rr<r_tb) {
					A_sum += area;
					Val_sum += Diss*area;
				}	
			}

			ViscDisp[k+1]=Val_sum/(A_sum+1.e-19);	

		}
		if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		I = TECDAT100(&III, x, &DIsDouble);
				

		
		/* Viscous Diffusion */
		
		// Make dk/dxj vector
		Vec GradK;
		Cmpnts ***gradk;
		
		DMCreateGlobalVector(user[bi].fda, &GradK);
		DMDAVecGetArray(user[bi].fda, GradK, &gradk);
		
		for (k=zs; k<ze-2; k++)
		for (j=ys; j<ye-2; j++)
		for (i=xs; i<xe-2; i++) {
			double _U = ucat[k+1][j+1][i+1].x;
			double _V = ucat[k+1][j+1][i+1].y;
			double _W = ucat[k+1][j+1][i+1].z;

			double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
			double eta0= eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
			double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
			double ajc = aj[k+1][j+1][i+1];
			double dudc, dude, dudz;
			double dvdc, dvde, dvdz;
			double dwdc, dwde, dwdz;
			double du_dx, du_dy, du_dz;
			double dv_dx, dv_dy, dv_dz;
			double dw_dx, dw_dy, dw_dz;
			double dk_dx, dk_dy, dk_dz;
		

			Compute_du_center (i+1, j+1, k+1, mx, my, mz, ucat, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
			Compute_du_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
	
			dk_dx = 0;
			dk_dx += _U*du_dx;
			dk_dx += _V*dv_dx;
			dk_dx += _W*dw_dx;
			
			dk_dy = 0;
			dk_dy += _U*du_dy;
			dk_dy += _V*dv_dy;
			dk_dy += _W*dw_dy;
			
			dk_dz = 0;
			dk_dz += _U*du_dz;
			dk_dz += _V*dv_dz;
			dk_dz += _W*dw_dz;
			
			gradk[k+1][j+1][i+1].x = dk_dx;
			gradk[k+1][j+1][i+1].y = dk_dy;
			gradk[k+1][j+1][i+1].z = dk_dz;
		}
		
		// Viscous Diffusion
		for (k=zs; k<ze-2; k++) {
			double A_sum=0.0, Val_sum=0.0;
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) {
				double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
				double eta0= eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
				double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
				double ajc = aj[k+1][j+1][i+1];
				double dudc, dude, dudz;
				double dvdc, dvde, dvdz;
				double dwdc, dwde, dwdz;
				double du_dx, du_dy, du_dz;
				double dv_dx, dv_dy, dv_dz;
				double dw_dx, dw_dy, dw_dz;
			
				Compute_du_center (i+1, j+1, k+1, mx, my, mz, gradk, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
				Compute_du_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
			
				double ddk_dxdx = du_dx;
				double ddk_dydy = dv_dy;
				double ddk_dzdz = dw_dz;
				double nu = 1./user[bi].ren;
			
				double ViscousDiffusion = 0;
			
				ViscousDiffusion = nu * (ddk_dxdx);
				ViscousDiffusion += nu * (ddk_dydy);
				ViscousDiffusion += nu * (ddk_dzdz);
			
				x[k * (mx-2)*(my-2) + j*(mx-2) + i] = ViscousDiffusion;



				double rr=sqrt(pow(coor[k+1][j+1][i+1].x-xc_tb,2)+pow(coor[k+1][j+1][i+1].y-yc_tb,2));
			       	double area=sqrt( zet[k+1][j+1][i+1].x*zet[k+1][j+1][i+1].x + zet[k+1][j+1][i+1].y*zet[k+1][j+1][i+1].y + zet[k+1][j+1][i+1].z*zet[k+1][j+1][i+1].z );

				if (rr<r_tb) {
					A_sum += area;
					Val_sum += ViscousDiffusion*area;
				}

			}
			
			ViscTrans[k+1]=Val_sum/(A_sum+1.e-19);
		}
		if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		I = TECDAT100(&III, x, &DIsDouble);
		
		DMDAVecRestoreArray(user[bi].fda, GradK, &gradk);
		VecDestroy(&GradK);
	
	
		/* Pressure Transport */
		Vec P_sum;
				
		PetscReal ***psum;
				
		DMCreateGlobalVector(user[bi].da, &P_sum);
				
		sprintf(filen, "%ssp_%06d_%1d.dat", prefix, ti, user[bi]._this);
		PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
		VecLoad(P_sum,viewer);
		PetscViewerDestroy(&viewer);
				
		DMDAVecGetArray(user[bi].da, P_sum, &psum);
				
		for (k=zs; k<ze-2; k++) {
			double A_sum=0.0, Val_sum=0.0;
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) {
				double dpdc, dpde, dpdz;
				double dp_dx, dp_dy, dp_dz;
				double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
				double eta0= eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
				double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
				double ajc = aj[k+1][j+1][i+1];
					
				Compute_dscalar_center (i+1, j+1, k+1, mx, my, mz, psum, nvert, &dpdc, &dpde, &dpdz);
				Compute_dscalar_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dpdc, dpde, dpdz, &dp_dx, &dp_dy, &dp_dz);

				dp_dx /= N;
				dp_dy /= N;
				dp_dz /= N;
			
				double U = ucat[k+1][j+1][i+1].x;
				double V = ucat[k+1][j+1][i+1].y;
				double W = ucat[k+1][j+1][i+1].z;
			
				double PressureTrans = 0;
				PressureTrans -= U * dp_dx;
				PressureTrans -= V * dp_dy;
				PressureTrans -= W * dp_dz;
					
				x[k * (mx-2)*(my-2) + j*(mx-2) + i] = PressureTrans;


				double rr=sqrt(pow(coor[k+1][j+1][i+1].x-xc_tb,2)+pow(coor[k+1][j+1][i+1].y-yc_tb,2));
			       	double area=sqrt( zet[k+1][j+1][i+1].x*zet[k+1][j+1][i+1].x + zet[k+1][j+1][i+1].y*zet[k+1][j+1][i+1].y + zet[k+1][j+1][i+1].z*zet[k+1][j+1][i+1].z );

				if (rr<r_tb) {
					A_sum += area;
					Val_sum += PressureTrans*area;
				}
			}
			PresTrans[k+1]=Val_sum/(A_sum+1.e-19);
		}
		if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		I = TECDAT100(&III, x, &DIsDouble);

		DMDAVecRestoreArray(user[bi].da, P_sum, &psum);
				
		VecDestroy(&P_sum);
	
	
		/* Turbulent Convection */
		Vec UUU_sum, UUU;
		Cmpnts ***uuu;
		
		DMCreateGlobalVector(user[bi].fda, &UUU);
		
		DMDAVecGetArray(user[bi].fda, UUU, &uuu);
		
		// a new vector for <ui'ui'uk'>
		for (k=zs; k<ze-2; k++)
		for (j=ys; j<ye-2; j++)
		for (i=xs; i<xe-2; i++) {
			double UU = uu[k+1][j+1][i+1];
			double VV = vv[k+1][j+1][i+1];
			double WW = ww[k+1][j+1][i+1];
			double UV = uv[k+1][j+1][i+1];
			double VW = vw[k+1][j+1][i+1];
			double WU = wu[k+1][j+1][i+1];
			
			double _U = ucat[k+1][j+1][i+1].x;
			double _V = ucat[k+1][j+1][i+1].y;
			double _W = ucat[k+1][j+1][i+1].z;
			
			uuu[k+1][j+1][i+1].x = _U*UU+_V*UV+_W*WU;
			uuu[k+1][j+1][i+1].y = _U*UV+_V*VV+_W*VW;
			uuu[k+1][j+1][i+1].z = _U*WU+_V*VW+_W*WW;
		}
		
		// calculate the divergence of <uuu>
		for (k=zs; k<ze-2; k++) {
			double A_sum=0.0, Val_sum=0.0;
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) {
				double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
				double eta0= eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
				double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
				double ajc = aj[k+1][j+1][i+1];
				double dudc, dude, dudz;
				double dvdc, dvde, dvdz;
				double dwdc, dwde, dwdz;
				double du_dx, du_dy, du_dz;
				double dv_dx, dv_dy, dv_dz;
				double dw_dx, dw_dy, dw_dz;
			
				Compute_du_center (i+1, j+1, k+1, mx, my, mz, uuu, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
				Compute_du_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
			
				double TurbulentConvection = du_dx + dv_dy + dw_dz;
				
				x[k * (mx-2)*(my-2) + j*(mx-2) + i] = -TurbulentConvection;


				double rr=sqrt(pow(coor[k+1][j+1][i+1].x-xc_tb,2)+pow(coor[k+1][j+1][i+1].y-yc_tb,2));
			       	double area=sqrt( zet[k+1][j+1][i+1].x*zet[k+1][j+1][i+1].x + zet[k+1][j+1][i+1].y*zet[k+1][j+1][i+1].y + zet[k+1][j+1][i+1].z*zet[k+1][j+1][i+1].z );

				if (rr<r_tb) {
					A_sum += area;
					Val_sum += -TurbulentConvection*area;
				}
			}
			TurbTrans[k+1] = Val_sum/(A_sum+1.e-19);
		}

		if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		I = TECDAT100(&III, x, &DIsDouble);

		
		// for only W related turbulence convection
		for (k=zs; k<ze-2; k++)
		for (j=ys; j<ye-2; j++)
		for (i=xs; i<xe-2; i++) {
			double UU = uu[k+1][j+1][i+1];
			double VV = vv[k+1][j+1][i+1];
			double WW = ww[k+1][j+1][i+1];
			double UV = uv[k+1][j+1][i+1];
			double VW = vw[k+1][j+1][i+1];
			double WU = wu[k+1][j+1][i+1];
			
			double _U = ucat[k+1][j+1][i+1].x;
			double _V = ucat[k+1][j+1][i+1].y;
			double _W = ucat[k+1][j+1][i+1].z;
			
			uuu[k+1][j+1][i+1].x = _W*WU;
			uuu[k+1][j+1][i+1].y = _W*VW;
			uuu[k+1][j+1][i+1].z = _W*WW;
		}
	

		// calculate the divergence of <uuu>.x
		for (k=zs; k<ze-2; k++) {
			double A_sum=0.0, Val_sum=0.0;
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) {
				double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
				double eta0= eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
				double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
				double ajc = aj[k+1][j+1][i+1];
				double dudc, dude, dudz;
				double dvdc, dvde, dvdz;
				double dwdc, dwde, dwdz;
				double du_dx, du_dy, du_dz;
				double dv_dx, dv_dy, dv_dz;
				double dw_dx, dw_dy, dw_dz;
			
				Compute_du_center (i+1, j+1, k+1, mx, my, mz, uuu, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
				Compute_du_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
			
				double TurbulentConvection = du_dx; // + dv_dy + dw_dz;
				
				x[k * (mx-2)*(my-2) + j*(mx-2) + i] = -TurbulentConvection;


				double rr=sqrt(pow(coor[k+1][j+1][i+1].x-xc_tb,2)+pow(coor[k+1][j+1][i+1].y-yc_tb,2));
			       	double area=sqrt( zet[k+1][j+1][i+1].x*zet[k+1][j+1][i+1].x + zet[k+1][j+1][i+1].y*zet[k+1][j+1][i+1].y + zet[k+1][j+1][i+1].z*zet[k+1][j+1][i+1].z );

				if (rr<r_tb) {
					A_sum += area;
					Val_sum += -TurbulentConvection*area;
				}
			}
			TurbTrans_x[k+1] = Val_sum/(A_sum+1.e-19);
		}

		if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		I = TECDAT100(&III, x, &DIsDouble);

			
		// calculate the divergence of <uuu>.y
		for (k=zs; k<ze-2; k++) {
			double A_sum=0.0, Val_sum=0.0;
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) {
				double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
				double eta0= eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
				double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
				double ajc = aj[k+1][j+1][i+1];
				double dudc, dude, dudz;
				double dvdc, dvde, dvdz;
				double dwdc, dwde, dwdz;
				double du_dx, du_dy, du_dz;
				double dv_dx, dv_dy, dv_dz;
				double dw_dx, dw_dy, dw_dz;
			
				Compute_du_center (i+1, j+1, k+1, mx, my, mz, uuu, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
				Compute_du_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
			
				double TurbulentConvection = dv_dy; // + dv_dy + dw_dz;
				
				x[k * (mx-2)*(my-2) + j*(mx-2) + i] = -TurbulentConvection;


				double rr=sqrt(pow(coor[k+1][j+1][i+1].x-xc_tb,2)+pow(coor[k+1][j+1][i+1].y-yc_tb,2));
			       	double area=sqrt( zet[k+1][j+1][i+1].x*zet[k+1][j+1][i+1].x + zet[k+1][j+1][i+1].y*zet[k+1][j+1][i+1].y + zet[k+1][j+1][i+1].z*zet[k+1][j+1][i+1].z );

				if (rr<r_tb) {
					A_sum += area;
					Val_sum += -TurbulentConvection*area;
				}
			}
			TurbTrans_y[k+1] = Val_sum/(A_sum+1.e-19);
		}

		if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		I = TECDAT100(&III, x, &DIsDouble);
	

			
		// calculate the divergence of <uuu>.z
		for (k=zs; k<ze-2; k++) {
			double A_sum=0.0, Val_sum=0.0;
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) {
				double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
				double eta0= eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
				double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
				double ajc = aj[k+1][j+1][i+1];
				double dudc, dude, dudz;
				double dvdc, dvde, dvdz;
				double dwdc, dwde, dwdz;
				double du_dx, du_dy, du_dz;
				double dv_dx, dv_dy, dv_dz;
				double dw_dx, dw_dy, dw_dz;
			
				Compute_du_center (i+1, j+1, k+1, mx, my, mz, uuu, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
				Compute_du_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
			
				double TurbulentConvection = dw_dz; // + dv_dy + dw_dz;
				
				x[k * (mx-2)*(my-2) + j*(mx-2) + i] = -TurbulentConvection;


				double rr=sqrt(pow(coor[k+1][j+1][i+1].x-xc_tb,2)+pow(coor[k+1][j+1][i+1].y-yc_tb,2));
			       	double area=sqrt( zet[k+1][j+1][i+1].x*zet[k+1][j+1][i+1].x + zet[k+1][j+1][i+1].y*zet[k+1][j+1][i+1].y + zet[k+1][j+1][i+1].z*zet[k+1][j+1][i+1].z );

				if (rr<r_tb) {
					A_sum += area;
					Val_sum += -TurbulentConvection*area;
				}
			}
			TurbTrans_z[k+1] = Val_sum/(A_sum+1.e-19);
		}

		if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		I = TECDAT100(&III, x, &DIsDouble);

		
		//DMDAVecRestoreArray(user[bi].fda, UUU_sum, &uuusum);
		DMDAVecRestoreArray(user[bi].fda, UUU, &uuu);
		//VecDestroy(&UUU_sum);
		VecDestroy(&UUU);

		
		/* Mean Convection */
		for (k=zs; k<ze-2; k++) {
			double A_sum=0.0, Val_sum=0.0;
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) {
				double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
				double eta0= eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
				double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
				double ajc = aj[k+1][j+1][i+1];
				double duudc, duude, duudz;
				double dvvdc, dvvde, dvvdz;
				double dwwdc, dwwde, dwwdz;
				double duu_dx, duu_dy, duu_dz;
				double dvv_dx, dvv_dy, dvv_dz;
				double dww_dx, dww_dy, dww_dz;
				double dk_dx, dk_dy, dk_dz;
			
				Compute_dscalar_center (i+1, j+1, k+1, mx, my, mz, uu_mean, nvert, &duudc, &duude, &duudz);
				Compute_dscalar_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, duudc, duude, duudz, &duu_dx, &duu_dy, &duu_dz);
			
				Compute_dscalar_center (i+1, j+1, k+1, mx, my, mz, vv_mean, nvert, &dvvdc, &dvvde, &dvvdz);
				Compute_dscalar_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dvvdc, dvvde, dvvdz, &dvv_dx, &dvv_dy, &dvv_dz);
			
				Compute_dscalar_center (i+1, j+1, k+1, mx, my, mz, ww_mean, nvert, &dwwdc, &dwwde, &dwwdz);
				Compute_dscalar_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dwwdc, dwwde, dwwdz, &dww_dx, &dww_dy, &dww_dz);
			
				dk_dx = 0;
				dk_dx += 0.5 * duu_dx;
				dk_dx += 0.5 * dvv_dx;
				dk_dx += 0.5 * dww_dx;
			
				dk_dy = 0;
				dk_dy += 0.5 * duu_dy;
				dk_dy += 0.5 * dvv_dy;
				dk_dy += 0.5 * dww_dy;
			
				dk_dz = 0;
				dk_dz += 0.5 * duu_dz;
				dk_dz += 0.5 * dvv_dz;
				dk_dz += 0.5 * dww_dz;
			
				double MeanConvection = 0;
			
				MeanConvection += ucat[k][j][i].x * dk_dx;
				MeanConvection += ucat[k][j][i].y * dk_dy;
				MeanConvection += ucat[k][j][i].z * dk_dz;
				MeanConvection *= -1.0;

				x[k * (mx-2)*(my-2) + j*(mx-2) + i] = MeanConvection;

				double rr=sqrt(pow(coor[k+1][j+1][i+1].x-xc_tb,2)+pow(coor[k+1][j+1][i+1].y-yc_tb,2));
			       	double area=sqrt( zet[k+1][j+1][i+1].x*zet[k+1][j+1][i+1].x + zet[k+1][j+1][i+1].y*zet[k+1][j+1][i+1].y + zet[k+1][j+1][i+1].z*zet[k+1][j+1][i+1].z );

				if (rr<r_tb) {
					A_sum += area;
					Val_sum += MeanConvection*area;
				}
			}
		
			MeanConv[k+1]=Val_sum/(A_sum+1.e-19);

		}
		if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		I = TECDAT100(&III, x, &DIsDouble);
	
	
		/* dW3dz */
		for (k=zs; k<ze-2; k++) {
			double A_sum=0.0, Val_sum=0.0;
			for (j=ys; j<ye-2; j++)
			for (i=xs; i<xe-2; i++) {
				double csi0 = csi[k+1][j+1][i+1].x, csi1 = csi[k+1][j+1][i+1].y, csi2 = csi[k+1][j+1][i+1].z;
				double eta0= eta[k+1][j+1][i+1].x, eta1 = eta[k+1][j+1][i+1].y, eta2 = eta[k+1][j+1][i+1].z;
				double zet0 = zet[k+1][j+1][i+1].x, zet1 = zet[k+1][j+1][i+1].y, zet2 = zet[k+1][j+1][i+1].z;
				double ajc = aj[k+1][j+1][i+1];
				double duudc, duude, duudz;
				double dvvdc, dvvde, dvvdz;
				double dwwdc, dwwde, dwwdz;
				double duu_dx, duu_dy, duu_dz;
				double dvv_dx, dvv_dy, dvv_dz;
				double dww_dx, dww_dy, dww_dz;
				double dk_dx, dk_dy, dk_dz;
			
				Compute_dscalar_center (i+1, j+1, k+1, mx, my, mz, uu_mean, nvert, &duudc, &duude, &duudz);
				Compute_dscalar_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, duudc, duude, duudz, &duu_dx, &duu_dy, &duu_dz);
			
				Compute_dscalar_center (i+1, j+1, k+1, mx, my, mz, vv_mean, nvert, &dvvdc, &dvvde, &dvvdz);
				Compute_dscalar_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dvvdc, dvvde, dvvdz, &dvv_dx, &dvv_dy, &dvv_dz);
			
				Compute_dscalar_center (i+1, j+1, k+1, mx, my, mz, ww_mean, nvert, &dwwdc, &dwwde, &dwwdz);
				Compute_dscalar_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dwwdc, dwwde, dwwdz, &dww_dx, &dww_dy, &dww_dz);
			
				dk_dz = 0;
				dk_dz += 0.5 * dww_dz;
			
				double MeanConvection = 0;
			
				MeanConvection += ucat[k][j][i].z * dk_dz;
				MeanConvection *= -1.0;

				x[k * (mx-2)*(my-2) + j*(mx-2) + i] = MeanConvection;

				double rr=sqrt(pow(coor[k+1][j+1][i+1].x-xc_tb,2)+pow(coor[k+1][j+1][i+1].y-yc_tb,2));
			       	double area=sqrt( zet[k+1][j+1][i+1].x*zet[k+1][j+1][i+1].x + zet[k+1][j+1][i+1].y*zet[k+1][j+1][i+1].y + zet[k+1][j+1][i+1].z*zet[k+1][j+1][i+1].z );

				if (rr<r_tb) {
					A_sum += area;
					Val_sum += MeanConvection*area;
				}
			}
		
			dW3dz[k+1]=Val_sum/(A_sum+1.e-19);

		}
		if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
		I = TECDAT100(&III, x, &DIsDouble);
	
	
		FILE *f;
		char fname[128];
	   	sprintf(fname, "DiskAveraged_mkeBudget.plt");
		f = fopen(fname, "w");
	      	PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=\"z\", \"mke\", \"TurbulenceProduction\", \"TurbulenceConvection\", \"ViscousDissipation\", \"PressureTransport\", \"MeanConvection\", \"ViscousDiffusion\", \"TurbulenceConvection_x\", \"TurbulenceConvection_y\", \"TurbulenceConvection_z\", \"dw3dz/3\" \n");

		i=1;j=1;
		for (k=zs+1; k<ze-1; k++) {
      			PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le %le %le %le %le %le %le %le %le %le\n", coor[k][j][i].z, mke[k],TurbProd[k],TurbTrans[k],ViscDisp[k],PresTrans[k],MeanConv[k],ViscTrans[k], TurbTrans_x[k], TurbTrans_y[k], TurbTrans_z[k], dW3dz[k]);
		}
	      	fclose(f);


	
		DMDAVecRestoreArray(fda, Coor, &coor);
		DMDAVecRestoreArray(user[bi].fda, Ucat, &ucat);
		
		DMDAVecRestoreArray(user[bi].da, UU, &uu);
		DMDAVecRestoreArray(user[bi].da, VV, &vv);
		DMDAVecRestoreArray(user[bi].da, WW, &ww);
		DMDAVecRestoreArray(user[bi].da, UV, &uv);
		DMDAVecRestoreArray(user[bi].da, VW, &vw);
		DMDAVecRestoreArray(user[bi].da, WU, &wu);
		DMDAVecRestoreArray(user[bi].da, UU_Mean, &uu_mean);
		DMDAVecRestoreArray(user[bi].da, VV_Mean, &vv_mean);
		DMDAVecRestoreArray(user[bi].da, WW_Mean, &ww_mean);
				
		DMDAVecRestoreArray(user[bi].fda, user[bi].Csi, &csi);
		DMDAVecRestoreArray(user[bi].fda, user[bi].Eta, &eta);
		DMDAVecRestoreArray(user[bi].fda, user[bi].Zet, &zet);
		DMDAVecRestoreArray(user[bi].da, user[bi].Aj, &aj);
		DMDAVecRestoreArray(user[bi].da, user[bi].Nvert, &nvert);
		
		DMDAVecRestoreArray(user[bi].fda, user[bi].Ucat_cross_sum, &u1sum);
		DMDAVecRestoreArray(user[bi].fda, user[bi].Ucat_square_sum, &u2sum);
		
		delete []x;
		
		VecDestroy(&user[bi].Ucat_cross_sum);
		VecDestroy(&user[bi].Ucat_square_sum);
		VecDestroy(&Ucat);
		VecDestroy(&UU);
		VecDestroy(&VV);
		VecDestroy(&WW);
		VecDestroy(&UV);
		VecDestroy(&VW);
		VecDestroy(&WU);
		VecDestroy(&UU_Mean);
		VecDestroy(&VV_Mean);
		VecDestroy(&WW_Mean);
		printf("MKE Budget 13 !\n");
	}
	
	I = TECEND100();
	
	return 0;
}






// end add 
