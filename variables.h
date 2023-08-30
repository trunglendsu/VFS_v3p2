#ifndef _VARIABLES_H_
#define _VARIABLES_H_

#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
#include "HYPRE_struct_ls.h"
#include "HYPRE_sstruct_ls.h"
#include "HYPRE_IJ_mv.h"
#include <unistd.h>
#include<string.h>
#include<stdio.h> 

#include <stdlib.h>
#include "petscvec.h"
#include "petscdmda.h"
#include "petscksp.h"
#include "petscsnes.h"
#include <assert.h>
#include <time.h>


// Seokkoo
#if defined (__CPLUSPLUS__) || defined (__cplusplus)

#include <vector>
#include <algorithm>
#include <assert.h>
#include <complex>

#endif

#ifndef M_PI
#define M_PI 3.141592653589793238462643
#endif

typedef struct {
	PetscScalar x, y, z;
} Cmpnts;

typedef struct {
	double x, y, z;
	int active, color;
} Particle;

typedef struct {
	PetscScalar x, y;
} Cmpnts2;

typedef struct {
	PetscReal	t, f;
} FlowWave;

typedef struct {
	PetscReal	x, y;
} Cpt2D;

typedef struct {
	Vec	Ubcs; // An ugly hack, waste of memory
} BCS;

// xiaolei add SEDI
typedef struct {
	PetscInt	c1,c2,c3 ;
} Cell2Cell;


typedef struct {
	PetscInt	c[100];
} Node2Cell;


typedef struct {
	PetscInt	c[100];
} Cell2AllCell;


// xiaolei end add

// xiaolei add
//

typedef struct {
	int i,j,k;
} ijkIndex;


// xiaolei add end


/*
typedef struct {
	Cmpnts csi, eta, zet, ucont, ucat;
	double aj, p;
	int process, nvert;
} Neighbor_Element;	// seokkoo

typedef struct {
	Cmpnts csi, eta, zet, ucont, ucat;
	double aj, p;
	int nvert;
	struct Neighbor_Element *neighbor[6];
	struct Element *next;
} Element;	// seokkoo
*/
typedef struct {
  int i1, j1, k1;
  int i2, j2, k2;
  int i3, j3, k3;
  PetscReal	cr1, cr2, cr3; // coefficients
  PetscReal	d_i; // distance to interception point on grid cells
  int	imode; // interception mode

  int	ni, nj, nk;	// fluid cell
  PetscReal	d_s; // shortest distance to solid surfaces
  Cmpnts	pmin;
  int	cell; // shortest distance surface cell
  PetscReal	cs1, cs2, cs3;

  int	i11, j11, k11;
  int	i22, j22, k22;
  int	i33, j33, k33;
  PetscReal	cr11, cr22, cr33; // coefficients
  PetscReal	d_ii; // distance to interception point on grid cells
  int	iimode; // interception mode
  PetscReal	cs11, cs22, cs33;

  int	ii1, jj1, kk1;
  int	ii2, jj2, kk2;
  int	ii3, jj3, kk3;
  PetscReal	ct1, ct2, ct3; // coefficients
  //PetscReal	d_s; // distance to interception point on grid cells
  int	smode; // interception mode

  int	ii11, jj11, kk11;
  int	ii22, jj22, kk22;
  int	ii33, jj33, kk33;
  PetscReal	ct11, ct22, ct33; // coefficients
  PetscReal	d_ss; // distance to interception point on grid cells
  int	ssmode; // interception mode

  
/*   int      bi1, bj1, bk1; */
/*   int      bi2, bj2, bk2; */
/*   int      bi3, bj3, bk3; */
/*   int      bi4, bj4, bk4; */

/*   PetscReal     bcr1,bcr2,bcr3,bcr4; */

// sediment 
//


  PetscReal     ucomp, vcomp, wcomp, utaub; // velocity components and bed shear stress for each IB node
  PetscReal     CC; // Concentration above the bed cell



// end sediment


} IBMInfo;

typedef struct {
	int	nbnumber;
	int	n_v, n_elmt;	// number of vertices and number of elements
	int	my_n_v, my_n_elmt;	// seokkoo, my proc
	int	*nv1, *nv2, *nv3;	// Node index of each triangle
	PetscReal	*nf_x, *nf_y, *nf_z;	// Normal direction
	//PetscReal	*nf_x0, *nf_y0, *nf_z0;	// Normal direction
	PetscReal	*x_bp, *y_bp, *z_bp;	// Coordinates of IBM surface nodes
	PetscReal	*x_bp0, *y_bp0, *z_bp0;
        PetscReal       *x_bp_i, *y_bp_i, *z_bp_i;
	PetscReal	*x_bp_o, *y_bp_o, *z_bp_o;
	//PetscReal	x_bp_in[101][3270], y_bp_in[101][3270], z_bp_in[101][3270];
	Cmpnts	*u, *uold, *urm1;
  
	// Added 4/1/06 iman
	PetscReal     *dA ;         // area of an element
	// Added 4/2/06 iman
	PetscReal     *nt_x, *nt_y, *nt_z; //tangent dir
	PetscReal     *ns_x, *ns_y, *ns_z; //azimuthal dir
	// Added 6/4/06 iman
	//Cmpnts        *cent; // center of the element 
	PetscReal     *cent_x,*cent_y,*cent_z;

	// for radius check
	Cmpnts *qvec;
	PetscReal *radvec; 
	
	//seokkoo
	
	//int *_nv1, *_nv2, *_nv3;	// for rank 0
	PetscReal *count;
	//int *local2global_elmt;
	//int total_n_elmt, total_n_v;
	
  //PetscReal *_x_bp, *_y_bp, *_z_bp;	// for rank 0
	PetscReal *shear;
	PetscReal *mean_shear;
	PetscReal *reynolds_stress1;
	PetscReal *reynolds_stress2;
	PetscReal *reynolds_stress3;
	PetscReal *pressure;
	Cmpnts *rel_velocity;	// flow velocity - body velocity


	// add begin (xiaolei)
	/* for calculating surface force */
	int	*ib_elmt, *jb_elmt, *kb_elmt; // the lowest interpolation point for element

	PetscReal	*xb_elmt, *yb_elmt, *zb_elmt; // normal extension from the surface element center

	PetscReal	*p_elmt;

	Cmpnts		*tau_elmt;


        // add being (xiaolei)
        PetscReal 	*F_lagr_x, *F_lagr_y, *F_lagr_z; // force at the IB surface points (lagrange points)
        PetscReal 	*U_lagr_x, *U_lagr_y, *U_lagr_z;
        int        *i_min, *i_max, *j_min, *j_max, *k_min, *k_max;

        /* ACL */
        PetscReal       *angle_attack, *angle_pitch, *chord_blade; // twist angle and chord length at each grid point
        PetscReal       *CD, *CL;

	PetscReal 	U_ref;

	PetscReal	*dhx, *dhy, *dhz;

	PetscReal 	CD_bluff;

	PetscReal	indf_axis, Tipspeedratio;

	Cell2Cell	*c2c;
	Cell2AllCell	*c2ac;
	Node2Cell	*n2c;


        // add end (xiaolei)
	
	// sediment begin 
    // Added 1/10/2010 ali
    Cmpnts        *Bvel;
    PetscReal     *Bvel_u, *Bvel_v, *Bvel_w;
    PetscReal     *vx, *vy, *vz;
    PetscReal     *deltz_p_us, *deltz_p_ds, *A_us, *A_ds;
    PetscReal     *BShS, *BCShS, *Shvel;
    PetscReal     *SCont;
    PetscReal     *C;
    PetscReal     *max_bed_angle, *SedFlux;
    int      *Rigidity, *Mobility;
    PetscReal     *Delz, *dzbp,*dVol,*cent_z_old, *cent_zl, *cent_z_AVE, *atke_old, *z_bp_l;
    PetscReal      Deltab;
    PetscReal     *elmt_depth, *scc, *sbb;
    PetscReal     *u_grad_x,*v_grad_x,*c_grad_x;//,*z_grad_x;
    PetscReal     *u_grad_y,*v_grad_y,*c_grad_y;//,*z_grad_y;
    PetscReal     *u_max,*v_max,*c_max;//,*z_max;
    PetscReal     *netflux_old;
    PetscReal     *F_flux_12,*F_flux_13,*F_flux_23;
    PetscReal     dtime,time_bedchange;
	
	// xiaolei add SEDI
	
	int  *color;
	// end add SEDI

	// sediment end

} IBMNodes;

typedef struct node {
	int Node;
	struct node *next;
} node;

typedef struct list{
	node *head;
} LIST;


typedef struct list_node {
	int	index;
	struct list_node *next;
} Node_List;

typedef struct IBMListNode {
	IBMInfo ibm_intp;
	struct IBMListNode* next;
} IBMListNode;

typedef struct IBMList {
	IBMListNode *head;
} IBMList;


typedef struct UserCtx {
	
	DM da;	/* Data structure for scalars (include the grid geometry informaion, to obtain the grid information, use DMDAGetCoordinates) */
	DM fda;	// Data Structure for vectors
	DM fda2;	// Data Structure for vectors with 2 variables
	DMDALocalInfo info;

	Vec	Cent;	// Coordinates of cell centers
	Vec	Csi, Eta, Zet, Aj;
	Vec	ICsi, IEta, IZet, IAj;
	Vec	JCsi, JEta, JZet, JAj;
	Vec	KCsi, KEta, KZet, KAj;
	//Vec	/*Area,lArea,*/ Volume, lVolume;
	Vec	Ucont;	// Contravariant velocity components
	Vec	Ucat;	// Cartesian velocity components
	Vec	Ucat_o;
	Vec	Ucont_o, Ucont_rm1, Rhs, dUcont, pUcont;
  
	//Vec	lUcont_c;
	Vec	RHS_o;//, RHS_rm1;		// for AM2 scheme
	Vec	dP;
//	Vec	RHS_IRK[6];	// for implicit RK
	int	current_IRK_stage; 
	//Vec	Ucont_rm2;//Conv_o, Visc_o;//, Conv_rm1;// 	// seokkoo

	Vec Fp;
	Vec	Div1, Div2, Div3;
	Vec	Visc1, Visc2, Visc3;
  
	Vec	K_Omega, lK_Omega, K_Omega_o, lK_Omega_o;//, K_Omega_rm1;
	Vec	Conc, lConc, Conc_o, lConc_o;
	Vec	Distance;
	Vec	lNu_t, lF1;		// eddy viscosity, f1 constant for SST model
	Vec	lUcat_old;
  
	Vec	P, P_o;
	Vec	Phi;
	Vec	GridSpace;
	Vec	Nvert;
	Vec	Nvert_o;
	Vec 	lNvert_o_fixed;
	Vec	*lNvert_tmp;	//Nvert_tmp[ibi]
	Vec   Itfc, ItfcP;
	BCS	Bcs;

	int	rhs_count;
	Vec	rhsD;
  //Vec Density;
// seokkoo, Dirichlet BC RHS
	Vec	Levelset, Levelset_o, lLevelset, lDensity, lMu;	// seokkoo, density of fluids
	Vec	LevelsetGrad;
	Vec	lST;	//surface tension
	Vec	lCs;
	Vec	lUstar;
	Vec	Gid, Gidm;				// seokkoo, Global ID for local array containg ghost points
	Vec	Phi2, B2, Ucont2;
	
	int	local_Phi2_size, p_global_begin;
	int	reduced_p_size;
	FILE *fp_inflow_u, *fp_inflow_ko, *fp_inflow_t;
  
	double shear_stress[6], shear_force_avg[6];
	double mean_k_flux, mean_k_area;
	double mean_k_flux_ibnode, mean_k_area_ibnode;
	double *local_bulk_velocity; 
	//std::vector<double> k_area;
	double *k_area, *k_area_ibnode;
	double ustar_avg, Re_tau_avg;
	double ustar_now[6];
	
	double lA_cyl, A_cyl;
	double lA_cyl_x, A_cyl_x;
	double lA_cyl_z, A_cyl_z;
	double lFvx_cyl, lFvz_cyl, Fvx_cyl, Fvz_cyl;
	double lFpx_cyl, lFpz_cyl, Fpx_cyl, Fpz_cyl;
  
	SNES snes;
	Vec rhs;
	Mat J;
	
	Vec	lUcont, lUcat, lP, lPhi;
	Vec	lUcont_o, lUcont_rm1;//, ldUcont;
	Vec	lCsi, lEta, lZet, lAj;
	Vec	lICsi, lIEta, lIZet, lIAj;
	Vec	lJCsi, lJEta, lJZet, lJAj;
	Vec	lKCsi, lKEta, lKZet, lKAj;
	Vec	lGridSpace;
	Vec	lNvert, lNvert_o, lNFace;
	Vec	lCent;
	Vec   lItfc, lItfcP;
	Vec	inletU;
	Vec	nhostU, nhostP;
	Vec	DUold;
	Vec	Forcing;
	Vec	Ucont_MG;

	Vec	Dt;

	AO	ao;

	PetscReal	ren;	// Reynolds number
	PetscReal	dt; 	// time step
	PetscReal	st;	// Strouhal number
  
	// added for new poisson solver
	PetscReal     FluxIntpSum;
	// Added Iman
	Vec           psuedot;
	PetscReal     cfl, vnn;

	PetscReal	r[101], tin[101], uinr[101][1001];

	PetscReal *itfchostx, *itfchosty, *itfchostz;
	PetscReal FluxInSum, FluxOutSum;

	PetscErrorCode aotopetsc;
	PetscBool assignedA;

	int _this;
	PetscInt *idx_from;

	int bctype[6];
	int itfcptsnumber;
	int *itfcI, *itfcJ, *itfcK;
	int *itfchostI, *itfchostJ, *itfchostK, *itfchostB;
	int	IM, JM, KM; // dimensions of grid

	int ibmnumber;
	IBMInfo  *ibm_intp;
	Mat	A, C;
	KSP   ksp;

	IBMNodes *ibm;

	DM	*da_f, *da_c;
	struct UserCtx *user_f, *user_c;
	Vec	*lNvert_c;

	Vec	B;
	Vec	Rhsp, X, R;
  
	Mat	MR, MP;
	MatNullSpace nullsp;

  
  /* Variables for multi-nullspace case */
	int *KSKE;
	PetscBool multinullspace;

	IBMList *ibmlist;

	int thislevel;
	PetscInt mglevels;

	PetscBool *isc, *jsc, *ksc;
  
	FlowWave *inflow;
	int number_flowwave;
  
	
  #if defined (__CPLUSPLUS__) || defined (__cplusplus)
	/*
	int map_set;
	int solver_set;
	bool has_singletone;
	Epetra_Map *Map;
	Epetra_CrsMatrix *Mat_P;
	Epetra_Vector *Vec_P, *Vec_RHS;
	
	Epetra_LinearProblem *Poisson_Problem;
	ML_Epetra::MultiLevelPreconditioner *MLPrec;
	Epetra_CrsSingletonFilter *SingletonFilter;
	Epetra_LinearProblem *Reduced_Problem;
	*/
  #endif
	Vec Ucat_sum;		// sum of u, v, w over time
	Vec Ucat_cross_sum;	// sum of uv, vw, wu
	Vec Ucat_square_sum;	// sum of uu, vv, ww
	Vec P_sum;			// sum of p
	Vec P_square_sum;	// sun of p*p
	
	Vec Ucat_sum_phase;		// sum of u, v, w over time
	Vec Ucat_cross_sum_phase;	// sum of uv, vw, wu
	Vec Ucat_square_sum_phase;	// sum of uu, vv, ww
	Vec P_sum_phase;			// sum of p
	Vec P_square_sum_phase;	// sun of p*p
        Vec Conc_sum;          //sum of Concentration	
	//Vec Nut_sum; // sum of eddy viscosity
	Vec K_sum; // sum of Kinetic energy in RANS
	Vec Levelset_sum; // sum of levelset
  Vec Levelset_square_sum; // sum of levelset^2
	/*Vec P_cross_sum;*/
	
	Vec Udp_sum; //size 1; u*dpdx + v*dpdy + w*dpdz
	Vec dU2_sum; //size 3; (dui_dx)^2 + (dui_dy)^2 + (dui_dz)^2;  ... (3*3=9 components)
	Vec UUU_sum; //size 3; (u^2+v^2+w^2)*ui
	//Vec tauS_sum; // size 1; sum of tau_ij Sij = 2 nu_t Sij Sij; tau_ij = +2 nu_t Sij
	
	Vec Vort_sum;
	Vec Vort_square_sum;
	
	Vec lSx, lSy, lSz, lS;
	Vec lLM, lMM, lNM;
	
	//int **free_surface_j;		// free_surf[i][k] = (top j cell)
	double **free_surface_location;	// y or z position[i][k]; only for rank=0
	//double **free_surface_p;	// Dirichlet pressure BC
	//double **vof;			// volume fraction at surface
	//int **bottom_surface_j;		// bottom_surf[i][k] = (bottom j cell)
	//double **bottom_surface_y;	// y position
	Cmpnts **ucat_plane;		// for pseudo periodic BC
	Cmpnts2 **komega_plane;
	std::vector<Particle> particles;

	// begin add (xiaolei)
        Vec     Visc1_tmprt, Visc2_tmprt, Visc3_tmprt;
  
        Vec     Visc1_wm, Visc2_wm, Visc3_wm;   
        Vec     Visc1_wmtmprt, Visc2_wmtmprt, Visc3_wmtmprt;  

        Vec     lVisc1, lVisc2, lVisc3;
        Vec     lVisc1_tmprt, lVisc2_tmprt, lVisc3_tmprt;

        Vec     lVisc1_wm, lVisc2_wm, lVisc3_wm;   
        Vec     lVisc1_wmtmprt, lVisc2_wmtmprt, lVisc3_wmtmprt;   

        PetscReal       Pr;     // Prantle Number
        PetscReal       Ri_x;   // Rechadson Number 
        PetscReal       Ri_y;   // Rechadson Number 
        PetscReal       Ri_z;   // Rechadson Number

        double **tmprt_plane;


	// levelset inflow
        double **level_plane;

	// rotor_model
        Vec F_eul;              // Force on the background grids
        Vec lF_eul;

        Vec Force_wm;
        Vec lForce_wm;

        PetscReal *p_jmin, *p_jmax, *friction_jmin, *friction_jmax;

        // varibles for temperature

        Vec Tmprt, lTmprt, Tmprt_o, lTmprt_o, Tmprt_rm1, lTmprt_rm1;//, K_Omega_rm1;

        Vec FTmprt;             //buoyancy force
        Vec lFTmprt;

        Vec T_sum;

        Vec TT_sum, TU_sum, TP_sum, Prt_sum;

        Vec Pr_t, lPr_t;

        Vec Force_wmtmprt;
        Vec lForce_wmtmprt;

        Vec     lTau;  // shear stress on the solid surface
        Vec     lQ_scalar;  // scalar flux on the solid surface

        int bctype_tmprt[6];

        PetscReal tmprts[6]; // the temperature value/flux at the six boundary

	// Sponge layer for levelset (Porous medium)

	Vec lFSponge;	
	Vec FSponge;	

	
	Vec Curv;
	Vec Grad_abs;
	Vec Heaviside;

	double *areak_air, *areak_water, *fluxk_air, *fluxk_water; 

	// add end (xiaolei)


	// sediment    
	//
	

	//Vec	Conc, lConc, Conc_o, lConc_o;

	// end sediment


} UserCtx;

typedef struct {
	UserCtx *user;
	int thislevel;
} MGCtx;

typedef struct {
	PetscInt mglevels;
	int thislevel;

	PetscBool  isc, jsc, ksc;
	MGCtx *mgctx;
} UserMG;

typedef struct {
	PetscReal	P;    //Press on the surface elmt
	int		n_P; //number of Press Pts on the elmt
	PetscReal	Tow_ws, Tow_wt; //wall shear stress of the elmt
	PetscReal	Tow_wn; // normal stress 

	int		Clsnbpt_i,Clsnbpt_j,Clsnbpt_k; //Closest Near Bndry Pt to each surf elmt 
	int		icell,jcell,kcell;
	int		FoundAroundcell;
	int		Need3rdPoint;
	//int      Aroundcellnum;
} SurfElmtInfo;

typedef struct {
	PetscReal    S_new[6],S_old[6],S_real[6],S_realm1[6];
	PetscReal    S_ang_n[6],S_ang_o[6],S_ang_r[6],S_ang_rm1[6];
	PetscReal    red_vel, damp, mu_s; // reduced vel, damping coeff, mass coeff
	PetscReal    F_x,F_y,F_z, A_tot; //Forces & Area
	PetscReal    F_x_old,F_y_old,F_z_old; //Forces & Area
	PetscReal    F_x_real,F_y_real,F_z_real; //Forces & Area
	PetscReal    M_x,M_y,M_z; // Moments
	PetscReal    M_x_old,M_y_old,M_z_old; //Forces & Area
	PetscReal    M_x_real,M_y_real,M_z_real; //Forces & Area
	PetscReal    M_x_rm2,M_y_rm2,M_z_rm2; //Forces & Area
	PetscReal    M_x_rm3,M_y_rm3,M_z_rm3; //Forces & Area
	PetscReal    x_c,y_c,z_c; // center of rotation(mass)
	PetscReal    Mdpdn_x, Mdpdn_y,Mdpdn_z;
	PetscReal    Mdpdn_x_old, Mdpdn_y_old,Mdpdn_z_old;

	// Aitkin's iteration
	PetscReal    dS[6],dS_o[6],atk,atk_o;

	// for force calculation
	SurfElmtInfo  *elmtinfo;
	IBMInfo       *fsi_intp;

	//Max & Min of ibm domain where forces are calculated
	PetscReal    Max_xbc,Min_xbc;
	PetscReal    Max_ybc,Min_ybc;
	PetscReal    Max_zbc,Min_zbc;

	// CV bndry
	int CV_ys,CV_ye,CV_zs,CV_ze;

	// add begin (xiaolei)
        PetscReal    omega_x, omega_y, omega_z; 
        
        PetscReal    nx_tb, ny_tb, nz_tb; // direction vector of rotor axis rotor_model

        PetscReal    angvel_z, angvel_x, angvel_y, angvel_axis;

        PetscReal    x_c0,y_c0,z_c0;


	// add end (xiaolei)
	
}	FSInfo;

// add begin (xiaolei)
typedef struct {

        /* for actuator line model */

        int        num_AC, num_CD, num_CL; // num_AC: number of input points for angle of pitch and chord length of blade
        PetscReal       *r_AC, *angle_pitchInp, *chord_bladeInp;  // pitch angel and chord length from input
        PetscReal       *ang_CD, *CDInp, *ang_CL, *CLInp;
        PetscReal       r_beg, r_end; // This type airfoil begins at r=r_beg and ends at r=r_end


} ACL;
// end add (xiaolei)


//#define COEF_TIME_ACCURACY 1.5
extern PetscReal COEF_TIME_ACCURACY;

#if defined (__CPLUSPLUS__) || defined (__cplusplus)
void Contra2Cart(UserCtx *user);
void Contra2Cart_2(UserCtx *user);
void Contra2Contra(UserCtx *user, Vec lUcont_center);
void Contra2Cart_single(Cmpnts &csi, Cmpnts &eta, Cmpnts &zet, Cmpnts &ucont, Cmpnts *ucat);

void DestroyIBMList(IBMList *ilist);
void AddIBMNode(IBMList *ilist, IBMInfo ibm_intp);
void InitIBMList(IBMList *ilist);
void insertnode(LIST *ilist, int Node);
void destroy(LIST *ilist);
PetscErrorCode Blank_Interface(UserCtx *user);
int intsect_triangle(PetscReal orig[3], PetscReal dir[3], PetscReal vert0[3], PetscReal vert1[3], PetscReal vert2[3], PetscReal *t, PetscReal *u, PetscReal *v);
PetscBool ISLineTriangleIntp(Cmpnts p1, Cmpnts p2, IBMNodes *ibm, int ln_v);
int ISPointInTriangle(Cmpnts p, Cmpnts p1, Cmpnts p2, Cmpnts p3, PetscReal nfx, PetscReal nfy, PetscReal nfz);
PetscErrorCode Dis_P_Line(Cmpnts p, Cmpnts p1, Cmpnts p2, Cmpnts *po, PetscReal *d);
PetscErrorCode triangle_intp2(Cpt2D p, Cpt2D p1, Cpt2D p2, Cpt2D p3, IBMInfo *ibminfo);
PetscErrorCode triangle_intp(Cpt2D p, Cpt2D p1, Cpt2D p2, Cpt2D p3, IBMInfo *ibminfo);
PetscErrorCode triangle_intpp(Cpt2D p, Cpt2D p1, Cpt2D p2, Cpt2D p3, IBMInfo *ibminfo);
void initlist(LIST *ilist);
PetscErrorCode InflowWaveFormRead(UserCtx *user);
PetscErrorCode FormMetrics(UserCtx *user);
PetscErrorCode GridRestriction(int i, int j, int k, int *ih, int *jh, int *kh, UserCtx *user);
int lidx(int i, int j, int k, UserCtx *user);
PetscErrorCode MG_Initial(UserMG *usermg, IBMNodes *ibm);
//PetscErrorCode ibm_read_ucd(IBMNodes *ibm);
PetscErrorCode ibm_read_ucd(IBMNodes *ibm, int ibi);	// ibm_io.c
PetscErrorCode FsiInitialize(int n_elmt, FSInfo *fsi,int ibi);
PetscErrorCode FSI_DATA_Input(FSInfo *FSinf, int ti);
PetscErrorCode Elmt_Move_FSI_ROT(FSInfo *FSinfo, IBMNodes *ibm, PetscReal dt, int ibi);
PetscErrorCode Elmt_Move_FSI_TRANS(FSInfo *FSinfo, IBMNodes *ibm);
PetscErrorCode ibm_surface_out(IBMNodes *ibm, int ti, int ibi);
PetscErrorCode ibm_search_advanced(UserCtx *user, IBMNodes *ibm, int ibi);
PetscErrorCode ibm_interpolation_advanced(UserCtx *user,int Add_dUndt);
PetscErrorCode fluxin(UserCtx *user);
PetscErrorCode Struc_Solver(UserMG *usermg,IBMNodes *ibm, FSInfo *fsi, int itr_sc, PetscBool *DoSCLoop);
PetscErrorCode Flow_Solver(UserMG *usermg,IBMNodes *ibm, FSInfo *fsi, IBMNodes *wtm, ACL *acl, FSInfo *fsi_wt, IBMNodes *ibm_ACD, FSInfo *fsi_IBDelta,
IBMNodes *ibm_IBDelta);
PetscErrorCode MG_Finalize(UserMG *usermg);
PetscErrorCode GhostNodeVelocity(UserCtx *user);
PetscErrorCode InflowFlux(UserCtx *user) ;
PetscErrorCode OutflowFlux(UserCtx *user);
PetscErrorCode FormBCS(UserCtx *user, FSInfo *fsi);
PetscErrorCode Block_Interface_U(UserCtx *user);
PetscErrorCode ComputeRHS(UserCtx *user, int istage);
//PetscErrorCode FormFunction1(Vec Ucont, Vec Rhs, UserCtx *user);
//PetscErrorCode FormFunction1_FVM(Vec Ucont, Vec Rhs, UserCtx *user);
PetscErrorCode Spectral(UserCtx *user);
PetscErrorCode Convection(UserCtx *user, Vec Ucont, Vec Ucat, Vec Conv);
PetscErrorCode Viscous(UserCtx *user, Vec Ucont, Vec Ucat, Vec Visc);
PetscErrorCode OutflowVelocity(UserCtx *user, Vec Ucont);
PetscErrorCode Find_fsi_2nd_interp_Coeff(int i, int j, int k, int elmt, Cmpnts p, IBMInfo *ibminfo,UserCtx *user, IBMNodes *ibm);
PetscErrorCode Find_fsi_2nd_interp_Coeff2(int i, int j, int k, int elmt, Cmpnts p, IBMInfo *ibminfo, UserCtx *user, IBMNodes *ibm);
PetscReal detmnt(PetscReal a[3][3]);
PetscErrorCode MyFieldRestriction(UserCtx *user);
PetscErrorCode Calc_forces_SI(FSInfo *FSinfo,UserCtx *user, IBMNodes *ibm,int ti, int ibi, int bi);
PetscErrorCode Calc_FSI_pos_SC(FSInfo *FSinfo,IBMNodes *ibm, PetscReal dt, PetscReal dtime, PetscReal Re);
PetscErrorCode SwingCylinder(FSInfo *fsi, IBMNodes *ibm);
PetscErrorCode Calc_FSI_Ang(FSInfo *FSinfo,IBMNodes *ibm, PetscReal dt, PetscReal dtime,int ibi, UserCtx *user);
PetscErrorCode Calc_FSI_Ang_intg(FSInfo *FSinfo,IBMNodes *ibm, PetscReal dt,int itrSC, int ibi, UserCtx *user);
PetscErrorCode FSI_DATA_output(FSInfo *FSinf, int ti);
PetscErrorCode FSI_DATA_Output(FSInfo *FSinfo, int ibi);
PetscErrorCode CollisionDetectionOfCylinders(FSInfo *fsi, int NumberOfBodies);
PetscErrorCode MyNvertRestriction(UserCtx *user_h, UserCtx *user_c);
PetscErrorCode ImplicitMomentumSolver(UserCtx *user, IBMNodes *ibm, FSInfo *fsi);
PetscErrorCode ImplicitMomentumSolver1(UserCtx *user, IBMNodes *ibm, FSInfo *fsi);
PetscErrorCode ImpRK(UserCtx *user, IBMNodes *ibm, FSInfo *fsi);
PetscErrorCode RungeKutta(UserCtx *user, IBMNodes *ibm, FSInfo *fsi);
PetscErrorCode PoissonSolver_MG(UserMG *usermg, IBMNodes *ibm, IBMInfo *ibminfo);
PetscErrorCode PoissonSolver_MG_original(UserMG *usermg, IBMNodes *ibm, IBMInfo *ibminfo);
PetscErrorCode UpdatePressure(UserCtx *user);
PetscErrorCode Projection(UserCtx *user);
PetscErrorCode Divergence(UserCtx *user);
PetscErrorCode Ucont_P_Binary_Output(UserCtx *user);
PetscErrorCode Ucat_Binary_Output(UserCtx *user);
PetscErrorCode SetInitialGuessToOne(UserCtx *user);
PetscErrorCode VolumeFlux(UserCtx *user, Vec lUcor, PetscReal *ibm_Flux, PetscReal *ibm_Area);
PetscErrorCode triangle_intp3D(double x, double y, double z, 
						double x1, double y1, double z1, 
						double x2, double y2, double z2,
						double x3, double y3, double z3, IBMInfo *ibminfo);
PetscErrorCode RungeKutta_Advection(UserCtx *user, IBMNodes *ibm, FSInfo *fsi);
void Compute_du_center (int i, int j, int k,  int mx, int my, int mz, Cmpnts ***ucat, PetscReal ***nvert, double *dudc, double *dvdc, double *dwdc, double *dude, double *dvde, double *dwde, double *dudz, double *dvdz, double *dwdz);
void Compute_dscalar_center (int i, int j, int k,  int mx, int my, int mz, PetscReal ***K, PetscReal ***nvert, double *dkdc, double *dkde, double *dkdz);
//void Compute_dscalar_center (int i, int j, int k,  int mx, int my, int mz, PetscReal ***K, PetscReal ***nvert, double *dkdc, double *dkde, double *dkdz, double solid);
void Compute_du_dxyz (	double csi0, double csi1, double csi2, double eta0, double eta1, double eta2, double zet0, double zet1, double zet2, double ajc,
					double dudc, double dvdc, double dwdc, double dude, double dvde, double dwde, double dudz, double dvdz, double dwdz,
					double *du_dx, double *dv_dx, double *dw_dx, double *du_dy, double *dv_dy, double *dw_dy, double *du_dz, double *dv_dz, double *dw_dz );
void Compute_du_dxyz ( int i, int j, int k, int mx, int my, int mz, Cmpnts ***ucat, PetscReal ***nvert, Cmpnts ***csi, Cmpnts ***eta, Cmpnts ***zet, PetscReal ***aj, 
		       double *du_dx, double *du_dy, double *du_dz, double *dv_dx, double *dv_dy, double *dv_dz, double *dw_dx, double *dw_dy, double *dw_dz);
void Compute_dscalar_dxyz ( double csi0, double csi1, double csi2, double eta0, double eta1, double eta2, double zet0, double zet1, double zet2, double ajc,
							double dkdc, double dkde, double dkdz, double *dk_dx, double *dk_dy, double *dk_dz);
void Conv_Diff_IC(UserCtx *user);
void Conv_Diff_BC(UserCtx *user);
void Conv_Diff_Function_RHS(UserCtx *user);

void Compute_Distance_Function(UserCtx *user);
void K_Omega_IC(UserCtx *user);
void Compute_Smagorinsky_Constant_1(UserCtx *user, Vec Ucont, Vec Ucat);
void Compute_eddy_viscosity_LES(UserCtx *user);
PetscErrorCode PoissonLHSNew(UserCtx *user, IBMNodes *ibm, IBMInfo *ibminfo);
PetscErrorCode PoissonLHSNew(UserCtx *user, HYPRE_IJMatrix &Ap);
void Convert_Phi2_Phi(UserCtx *user);
void PoissonSolver_Hypre(UserCtx *user, IBMNodes *ibm, IBMInfo *ibminfo);
PetscErrorCode PoissonRHS2(UserCtx *user, Vec B);

void wall_function (double nu, double sc, double sb, Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, double nx, double ny, double nz);
void wall_function_fsi (double nu, double sc, double sb, Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, double nx, double ny, double nz);
void wall_function_slip (double nu, double sc, double sb, Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, double nx, double ny, double nz);
void wall_function_roughness (double nu, double ks, double sc, double sb, Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, double nx, double ny, double nz);
void wall_function_roughness_loglaw (double nu, double ks, double sc, double sb, Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, double nx, double ny, double nz);
double find_utau_loglaw_roughness(double nu, double u, double y, double guess, double ks);

void noslip (UserCtx *user, double sc, double sb, Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, double nx, double ny, double nz);
void freeslip (UserCtx *user, double sc, double sb, Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub, double nx, double ny, double nz);

void Create_Hypre_Solver();
void Create_Hypre_P_Matrix_Vector(UserCtx *user);
void MatHYPRE_IJMatrixCopy(Mat v,HYPRE_IJMatrix &ij);
void Petsc_to_Hypre_Vector(Vec A, HYPRE_IJVector &B, PetscInt i_lower);
void Hypre_to_Petsc_Vector(HYPRE_IJVector &B, Vec A, PetscInt i_lower);
void Calculate_normal(Cmpnts csi, Cmpnts eta, Cmpnts zet, double ni[3], double nj[3], double nk[3]);
void Calculate_normal_and_area(Cmpnts csi, Cmpnts eta, Cmpnts zet, double ni[3], double nj[3], double nk[3], double *Ai, double *Aj, double *Ak);

void Compute_du_i (int i, int j, int k, int mx, int my, int mz, Cmpnts ***ucat, PetscReal ***nvert, 
				double *dudc, double *dvdc, double *dwdc, 
				double *dude, double *dvde, double *dwde,
				double *dudz, double *dvdz, double *dwdz);
				
void Compute_du_j (int i, int j, int k, int mx, int my, int mz, Cmpnts ***ucat, PetscReal ***nvert, 
				double *dudc, double *dvdc, double *dwdc, 
				double *dude, double *dvde, double *dwde,
				double *dudz, double *dvdz, double *dwdz);

void Compute_du_k (int i, int j, int k, int mx, int my, int mz, Cmpnts ***ucat, PetscReal ***nvert, 
				double *dudc, double *dvdc, double *dwdc, 
				double *dude, double *dvde, double *dwde,
				double *dudz, double *dvdz, double *dwdz);
				
void Compute_du_dxyz (	double csi0, double csi1, double csi2, double eta0, double eta1, double eta2, double zet0, double zet1, double zet2, double ajc,
					double dudc, double dvdc, double dwdc, double dude, double dvde, double dwde, double dudz, double dvdz, double dwdz,
					double *du_dx, double *dv_dx, double *dw_dx, double *du_dy, double *dv_dy, double *dw_dy, double *du_dz, double *dv_dz, double *dw_dz );

void Compute_dscalar_i (int i, int j, int k, int mx, int my, int mz, PetscReal ***K, PetscReal ***nvert, double *dkdc, double *dkde, double *dkdz );
void Compute_dscalar_j (int i, int j, int k, int mx, int my, int mz, PetscReal ***K, PetscReal ***nvert, double *dkdc, double *dkde, double *dkdz );
void Compute_dscalar_k (int i, int j, int k, int mx, int my, int mz, PetscReal ***K, PetscReal ***nvert, double *dkdc, double *dkde, double *dkdz );

void Calculate_dxdydz(PetscReal ajc, Cmpnts csi, Cmpnts eta, Cmpnts zet, double *dx, double *dy, double *dz);
void AxByC ( double a, Cmpnts &X, double b, Cmpnts &Y, Cmpnts *C);
void AxC ( double a, Cmpnts &X, Cmpnts *C);
double PPM(double WW, double W, double E, double EE, double a);
void Compute_Q(UserCtx *user, Vec Q);
void Levelset_BC(UserCtx *user);
void Distance_Function_RHS (UserCtx *user, Vec Levelset_RHS, int wall_distance);
void Distance_Function_RHS_new (UserCtx *user, Vec Levelset_RHS, int wall_distance);
void Init_Levelset_Vectors(UserCtx *user);
void Destroy_LevelSet_Vectors(UserCtx *user);
void Reinit_Levelset(UserCtx *user);
void Levelset_Function_IC(UserCtx *user);
void Advect_Levelset(UserCtx *user, double dt);

double time_coeff();
void IB_BC_Ucat(UserCtx *user);
void Calc_ShearStress(UserCtx *user);
void write_shear_stress_ibm();
void Calc_k_Flux(UserCtx *user);//, double *Flux, double *Area);

void Compute_Density(UserCtx *user);
double find_utau_Cabot(double nu,  double u, double y, double guess, double dpdn);
double find_utau_Cabot_roughness(double nu, double u, double y, double guess, double dpdn, double ks);
double u_Cabot(double nu, double y, double utau, double dpdn);
double u_Cabot_roughness(double nu, double y, double utau, double dpdn, double ks_plus);
PetscErrorCode Formfunction_2(UserCtx *user, Vec Rhs, double scale);
PetscErrorCode FormFunction_SNES(SNES snes, Vec Ucont, Vec Rhs, void *ptr);;
void initial_guess_for_snes(UserCtx *user, Vec U);
void Pressure_Gradient(UserCtx *user, Vec dP);
double Calc_Minimum_dt (UserCtx *user);
void write_data (UserCtx *user);
void write_IBdata(UserCtx *user, IBMNodes *ibm); // xiaolei
void Set ( Cmpnts *A, double a );
double H (double p, double dx);	// Heaviside function
double dH (double p, double dx);	// its derivative
double mean ( double A, double B );
//double mean0 ( double A, double B );
double sign(double a);
double sign1(double a, double dx);
double mod_sign(double d0, double grad, double e);
void get_weight ( int i, int j, int k, int mx, int my, int mz, PetscReal ***aj, PetscReal ***nvert, double nv, double weight[3][3][3]);
double integrate_testfilter(double val[3][3][3], double w[3][3][3]);
double integrate_testfilter_simpson(double val[3][3][3], double w[3][3][3]);
void Compute_Surface_Tension(UserCtx *user);
void Compute_dlevel_center_levelset (int i, int j, int k,  int mx, int my, int mz, double sgn, int wall_distance, PetscReal ***level, PetscReal ***nvert, double *dldc, double *dlde, double *dldz);
void Mass_Correction_Levelset (UserCtx *user, Vec D);
double weno3(double f0, double f1, double f2, double f3, double wavespeed);
double eno2(double f0, double f1, double f2, double f3, double a);
void Calc_Inlet_Area(UserCtx *user);
void pre_integrate ();
int file_exist(char *str);
double Upwind(double W, double E, double a);
void Keep_Constant_Flux(UserCtx *user);
//PetscErrorCode VolumeFlux2(UserCtx *user, PetscReal *ibm_Flux, PetscReal *ibm_Area, int flg);
PetscErrorCode VolumeFlux(UserCtx *user, PetscReal *ibm_Flux, PetscReal *ibm_Area, int flg);
void Add_IBMFlux_to_Outlet(UserCtx *user, PetscReal ibm_Flux);
PetscErrorCode ibm_surface_out_with_pressure(IBMNodes *ibm, int ibi);
void Initialize_free_surface_location_vector(UserCtx *user);
void Calc_free_surface_location(UserCtx *user);
void IB_BC(UserCtx *user);
void GlobalSum_Root(PetscReal* local,PetscReal* result,MPI_Comm comm);
void GlobalMax_Root(PetscReal* local,PetscReal* result,MPI_Comm comm);
void GlobalMin_All(PetscScalar* local,PetscScalar* result,MPI_Comm comm);
void GlobalSum_All(PetscScalar* local,PetscScalar* result,MPI_Comm comm);
void GlobalMax_All(PetscScalar* local,PetscScalar* result,MPI_Comm comm);
void n_phase(int *phase, int *previous_ti);
PetscErrorCode Phase_Averaging_Output(UserCtx *user);
void Compute_Water_Volume(UserCtx *user, double *vol);
void str_to_buffer(char *str, std::vector<char> &large_buffer);
double near_wall_eddy_viscosity(double yplus);
double near_wall_eddy_viscosity(double yplus, double yp_shift);
PetscErrorCode MySNESMonitor(SNES snes, PetscInt n,PetscReal rnorm,void *dummy);

extern Vec LevelSet, LevelSet0, LevelSet_o;
extern HYPRE_Solver pcg_solver_p, precon_p;
extern HYPRE_IJMatrix Ap;
extern HYPRE_ParCSRMatrix par_Ap;
extern HYPRE_IJVector Vec_p, Vec_p_rhs;
extern HYPRE_ParVector par_Vec_p, par_Vec_p_rhs;

extern double imp_free_tol;
extern IBMNodes	*ibm_ptr;
extern FSInfo   *fsi_ptr;
extern UserCtx	*user_ptr;
extern PetscInt user_m, user_n, user_p;
extern double mean_pressure_gradient;
extern PetscInt testfilter_ik, testfilter_1d;
extern PetscInt i_periodic, j_periodic, k_periodic;
extern PetscInt ii_periodic, jj_periodic, kk_periodic;
extern PetscInt periodic;
extern PetscInt poisson_it, amg_agg, amg_coarsentype;
extern double amg_thresh;
extern PetscInt i_homo_filter, j_homo_filter, k_homo_filter;
extern PetscInt clark;
extern PetscInt inletprofile;
extern double inlet_flux;
extern PetscInt tistart;
extern PetscInt les, wallfunction, rans, lowRe, slipbody, central, second_order;
extern PetscInt ti, tiout,tiout_bed, mixed, tistart, inflow_recycle_perioid;
extern PetscInt NumberOfBodies, NumberOfRotatingBodies;
extern PetscInt immersed, invicid;
extern PetscInt TwoD, mixed, averaging, phase_averaging;
extern char path[256], gridfile[256];
extern double dx_min, di_min, dj_min, dk_min;
extern double di_max, dj_max, dk_max;

extern int my_rank;
extern int block_number;
extern int outflow_scale;

extern PetscReal max_cs;//, Fr;
extern IBMNodes	*ibm_ptr;
extern FSInfo *fsi_ptr;
extern PetscInt implicit, initialzero;
extern PetscInt movefsi, rotatefsi;
extern PetscBool dpdz_set;
extern int ib_bctype[128];
extern PetscInt dynamic_freq;
extern PetscInt laplacian, rotational, skew, tiout_ufield, tiend_ufield;
extern PetscInt levelset;
extern double rho_water, rho_air;
extern double mu_water, mu_air;
extern double dthick, dtau_ratio, d_parameter;
extern PetscBool dthick_set;
extern PetscInt inviscid, surface_tension, poisson;
extern double gravity_x, gravity_y, gravity_z;
extern double inlet_y, outlet_y, inlet_z, outlet_z;
extern PetscInt fix_outlet, fix_inlet, fix_level;
extern int setup_lidx2(UserCtx *user);
extern PetscReal FluxInSum, FluxOutSum;
extern PetscReal FluxInSum_gas, FluxOutSum_gas;
extern PetscBool inlet_y_flag, inlet_z_flag;
extern PetscInt ibm_search, thin;
extern PetscInt inlet_buffer_k;
extern PetscBool rough_set;
extern double roughness_size;
extern PetscReal roughness_ice;
extern PetscInt i_proc, j_proc, k_proc, binary_input, xyz_input;
extern PetscReal L_dim;
extern PetscReal  max_angle;
extern PetscInt cop, regime, fish, MHV;
extern int ucat_plane_allocated;
extern double kappa;

extern PetscBool rstart_flg;
extern PetscInt levelset_it;
extern PetscInt cross_diffusion;
extern PetscInt rotdir, radi, fish;
extern double St_exp;
extern double angvel;
extern PetscInt ti_lastsave;
extern double x_r, y_r, z_r;
extern PetscInt sloshing, bubble;
extern double sloshing_a, sloshing_b, sloshing_d;
extern double bubble_d, bubble_z, bubble_ws;
extern double water_vol, water_vol_o;
extern double fluct_rms;
extern PetscInt subdt_levelset;
extern PetscReal Flux_in, angle;
extern PetscInt initial_perturbation;
extern double poisson_tol;
extern PetscInt STRONG_COUPLING, dgf_z,dgf_y,dgf_x, dgf_az,dgf_ay,dgf_ax;
extern PetscInt save_inflow, save_inflow_period, save_inflow_minus, pseudo_periodic, inletprofile, read_inflow_period;
extern double dthick_const, dt_inflow;
extern int save_ksection[1000];
extern int nsave_ksection;

double eps_levelset2();
double d_levelset2();
extern void save_inflow_section(UserCtx *user);
extern void read_inflow_section(UserCtx *user);
extern PetscErrorCode Calc_Moments(FSInfo *FSinfo, IBMNodes *ibm, SurfElmtInfo *elmtinfo, PetscReal Re, int ti);
extern PetscErrorCode ibm_Surf_stress(UserCtx *user, IBMNodes *ibm, SurfElmtInfo *elmtinfo, int ibi);
extern void Set_Nvert_for_Wallfunction(UserCtx *user);
std::complex<double> ERF(std::complex<double> z);
void Correct_Volume_Flux_Levelset(UserCtx *user);
extern void Calculate_Covariant_metrics(double g[3][3], double G[3][3]);
extern void Calculate_normal(Cmpnts csi, Cmpnts eta, Cmpnts zet, double ni[3], double nj[3], double nk[3]);
extern void Convert_Phi2_Phi(UserCtx *user);
extern void Convert_Ucont2_Ucont(UserCtx *user);
extern void Convert_RHS_RHS2(UserCtx *user, Vec RHS, Vec RHS2);
extern void Convert_Ucont_Ucont2(UserCtx *user);
extern void K_Omega_Set_Constant(UserCtx *user);
extern void Solve_K_Omega(UserCtx *user);
extern void K_Omega_IC(UserCtx *user);
extern void Compute_Distance_Function(UserCtx *user);
extern PetscErrorCode Implicit_MatrixFree(UserCtx *user, IBMNodes *ibm, FSInfo *fsi);	  
extern PetscErrorCode Do_averaging(UserCtx *user);
extern PetscErrorCode KE_Output(UserCtx *user);
extern void store_k_section(UserCtx *user, int kplane);
extern double M(double a, double b);
extern void write_torque_ibm();
extern double levelset_thickness(double aj, Cmpnts csi, Cmpnts eta, Cmpnts zet);
void Compute_LevelsetGrad0(UserCtx *user, Vec V);

// add begin (xiaolei)
extern PetscErrorCode Calc_F_lagr(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int NumberOfObjects);
extern PetscErrorCode Calc_F_eul(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int NumberOfObjects);
extern PetscErrorCode Calc_forces_rotor(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int bi);
extern PetscErrorCode Pre_process(UserCtx *user, IBMNodes *ibm, int NumberOfObjects); 

extern void Force_wallfunc (double nu, double nu_t, double sb, double sc,
                            Cmpnts Ub, Cmpnts Uc, Cmpnts Ua, PetscInt bctype, double ks,
                            double nx, double ny, double nz, double *f_i, double *f_j, double *f_k, PetscReal *ustar);

extern double utau_powerlaw(double nu, double ut_mag, double sc);

extern void wallmodel_2( double dpdx_b, double dpdy_b, double dpdz_b, double nu, double *nu_t_b, double *nu_t_c, double sb, double sc, Cmpnts *Ub, Cmpnts Uc, Cmpnts Ua, PetscInt bctype, double ks,double nx, double ny, double nz, double *tau_w);


extern void wallmodel_1( double dpdx, double dpdy, double dpdz, double nu, double *nu_t_b, double nu_t_c, double sb, double sc, Cmpnts *Ub, Cmpnts Uc, Cmpnts Ua, PetscInt bctype, double ks, double nx, double ny, double nz, double *f_i, double *f_j, double *f_k, double *f_i_c, double *f_j_c, double *f_k_c, PetscReal *ustar, Cmpnts *Ug);

extern void wallmodel_tmprt( double sb, double nx, double ny, double nz, Cmpnts Ub, Cmpnts Ua, double Tb, double Ta, PetscInt bctype, double ks, double *qs);
extern PetscErrorCode Calc_F_wm(UserCtx *user);

extern PetscErrorCode Calc_F_wmIB(UserCtx *user,IBMNodes *ibm );


extern PetscErrorCode TECIOOut_F_wm(UserCtx *user); 

extern PetscErrorCode TECIOOut_rhs(UserCtx *user, Vec Rhs, char fname[80]);

extern PetscErrorCode RungeKutta1(UserCtx *user, IBMNodes *ibm,
                          FSInfo *fsi);

extern PetscErrorCode force_jwall(UserCtx *user);

extern PetscErrorCode ibm_surface_out_with_pressure1(IBMNodes *ibm, int ibi);

extern PetscErrorCode surfaceforce_preprocess(UserCtx *user, IBMNodes *ibm);

extern PetscErrorCode cal_surfaceforce(UserCtx *user, IBMNodes *ibm);


extern PetscErrorCode ACD_read(IBMNodes *ibm, int ibi, FSInfo *fsi, int OneDmZ);

extern PetscErrorCode UrefACD_read(IBMNodes *ibm, int ibi, FSInfo *fsi, int OneDmZ);

extern PetscErrorCode ibm_read_ACL(IBMNodes *ibm, int ibi, FSInfo *fsi);

extern PetscErrorCode ACL_read_ucd(IBMNodes *ibm, int ibi, FSInfo *fsi);

extern PetscErrorCode ibmDelta_read_ucd(IBMNodes *ibm, int ibi, FSInfo *fsi, int iname);

extern PetscErrorCode Calc_F_lagr_ACL(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, ACL *acl, int NumberOfObjects);

extern PetscErrorCode airfoil_ACL(ACL *acl, IBMNodes *ibm,  FSInfo *fsi);

extern PetscErrorCode Calc_forces_ACL(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int bi);


extern double integrate_xztestfilter_simpson(double val[3][3][3], double w[3][3][3]);  //xyang


extern PetscErrorCode Calc_F_lagr_noslip(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int NumberOfTurbines); // xyang


extern PetscErrorCode Uref_ACL(UserCtx *user, IBMNodes *ibm_ACL, IBMNodes *ibm_ACD, FSInfo *fsi_wt, int NumberOfObjects);

void Compute_dt_center (int i, int j, int k,  int mx, int my, int mz, PetscReal ***tmprt, PetscReal ***nvert,
                                double *dtdc, double *dtde, double *dtdz );


void Compute_dt_dxyz (        double csi0, double csi1, double csi2, double eta0, double eta1, double eta2, double zet0, double zet1, double zet2, double ajc, double dtdc, double dtde, double dtdz, double *dt_dx, double *dt_dy, double *dt_dz);

void RHS_Tmprt(UserCtx *user, Vec Tmprt_RHS);

void Tmprt_IC(UserCtx *user);

void Tmprt_BC(UserCtx *user);

extern PetscErrorCode FormFunction_Tmprt(SNES snes, Vec Tmprt, Vec Rhs, void *ptr);

void Solve_Tmprt(UserCtx *user);

void Force_Tmprt(UserCtx *user);


extern void Compute_Prt(UserCtx *user);

extern PetscErrorCode TECIOOut_rhs1(UserCtx *user, Vec Rhs);


extern PetscErrorCode Pre_process_eg(UserCtx *user, IBMNodes *ibm);

extern PetscErrorCode Calc_Flagrnslip_g(UserCtx *user, IBMNodes *ibm, FSInfo *fsi);

extern PetscErrorCode Calc_Ulagr_g(UserCtx *user, IBMNodes *ibm);

extern PetscErrorCode Calc_Ulagr_e(UserCtx *user, IBMNodes *ibm);

extern PetscErrorCode Calc_Feul_g(UserCtx *user, IBMNodes *ibm, FSInfo *fsi);

extern PetscErrorCode Add_fluc(UserCtx *user);
extern PetscErrorCode Add_fluc_tmprt(UserCtx *user);

extern void Comput_du_wmlocal(double nx, double ny, double nz, double t1x, double t1y, double t1z, double t2x, double t2y, double t2z, double du_dx,double dv_dx,double dw_dx,double du_dy,double dv_dy,double dw_dy,double du_dz,double dv_dz,double dw_dz, double *dut1dn, double *dut2dn, double *dundn, double *dut1dt1, double *dut2dt1, double *dundt1, double *dut1dt2, double *dut2dt2, double *dundt2);

extern void Comput_du_Compgrid(double dxdc, double dxde, double dxdz, double dydc, double dyde, double dydz, double dzdc, double dzde, double dzdz, double nx, double ny, double nz, double t1x, double t1y, double t1z, double t2x, double t2y, double t2z, double dut1dn, double dut2dn, double dundn, double dut1dt1, double dut2dt1, double dundt1, double dut1dt2, double dut2dt2, double dundt2, double *dudc, double *dvdc, double *dwdc, double *dude, double *dvde, double *dwde, double *dudz, double *dvdz, double *dwdz);

extern void Comput_JacobTensor_i(int i, int j, int k, int mx, int my, int mz, Cmpnts ***coor, double *dxdc, double *dxde, double *dxdz, double *dydc, double *dyde, double *dydz, double *dzdc, double *dzde, double *dzdz);

extern void Comput_JacobTensor_j(int i, int j, int k, int mx, int my, int mz, Cmpnts ***coor, double *dxdc, double *dxde, double *dxdz, double *dydc, double *dyde, double *dydz, double *dzdc, double *dzde, double *dzdz);

extern void Comput_JacobTensor_k(int i, int j, int k, int mx, int my, int mz, Cmpnts ***coor, double *dxdc, double *dxde, double *dxdz, double *dydc, double *dyde, double *dydz, double *dzdc, double *dzde, double *dzdz);

extern PetscErrorCode Formfunction_wm(UserCtx *user, Vec Rhs, double scale);
extern PetscErrorCode Visc_wm(UserCtx *user);

extern PetscErrorCode Formfunction_wm_tmprt(UserCtx *user, Vec Rhs, double scale);
extern PetscErrorCode Visc_wm_tmprt(UserCtx *user, IBMNodes *ibm);


extern void Computewm_du_i (int i, int j, int k, int mx, int my, int mz, Cmpnts ***ucat, PetscReal ***nvert,
                                double *dudc, double *dvdc, double *dwdc,
                                double *dude, double *dvde, double *dwde,
                                double *dudz, double *dvdz, double *dwdz,
                                double dudc_wm, double dvdc_wm, double dwdc_wm,
                                double dude_wm, double dvde_wm, double dwde_wm,
                                double dudz_wm, double dvdz_wm, double dwdz_wm);
extern void Computewm_du_j (int i, int j, int k, int mx, int my, int mz, Cmpnts ***ucat, PetscReal ***nvert,
                                double *dudc, double *dvdc, double *dwdc,
                                double *dude, double *dvde, double *dwde,
                                double *dudz, double *dvdz, double *dwdz,
                                double dudc_wm, double dvdc_wm, double dwdc_wm,
                                double dude_wm, double dvde_wm, double dwde_wm,
                                double dudz_wm, double dvdz_wm, double dwdz_wm);
extern void Computewm_du_k (int i, int j, int k, int mx, int my, int mz, Cmpnts ***ucat, PetscReal ***nvert,
                                double *dudc, double *dvdc, double *dwdc,
                                double *dude, double *dvde, double *dwde,
                                double *dudz, double *dvdz, double *dwdz,
                                double dudc_wm, double dvdc_wm, double dwdc_wm,
                                double dude_wm, double dvde_wm, double dwde_wm,
                                double dudz_wm, double dvdz_wm, double dwdz_wm);
extern void Compute1_du_i (int i, int j, int k, int mx, int my, int mz, Cmpnts ***ucat, PetscReal ***nvert,
                                double *dudc, double *dvdc, double *dwdc,
                                double *dude, double *dvde, double *dwde,
                                double *dudz, double *dvdz, double *dwdz);
extern void Compute1_du_j (int i, int j, int k, int mx, int my, int mz, Cmpnts ***ucat, PetscReal ***nvert,
                                double *dudc, double *dvdc, double *dwdc,
                                double *dude, double *dvde, double *dwde,
                                double *dudz, double *dvdz, double *dwdz);
extern void Compute1_du_k (int i, int j, int k, int mx, int my, int mz, Cmpnts ***ucat, PetscReal ***nvert,
                                double *dudc, double *dvdc, double *dwdc,
                                double *dude, double *dvde, double *dwde,
                                double *dudz, double *dvdz, double *dwdz);
extern void tridag(int n, double *a, double *b, double *c, double *r, double *u);

extern PetscErrorCode Pre_process_wtm(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int NumberOfObjects);
extern PetscErrorCode UpdateXYZ_MoveFrame(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int NumberOfObjects);

extern PetscErrorCode Trans1DUp_XYZ(IBMNodes *ibm, int NumberOfObjects);


extern PetscErrorCode rotor_Rot(FSInfo *FSinfo, IBMNodes *ibm, PetscReal dt, int ibi);

extern void rotate_xyz (double ti, double dt, double angvel, double x_c, double y_c, double z_c, double x_bp0, double y_bp0, double z_bp0, double *x_bp, double *y_bp, double *z_bp, double *rot_angle);



extern PetscInt imin_wm, imax_wm, jmin_wm, jmax_wm, kmin_wm, kmax_wm;
extern PetscInt imin_wmtmprt, imax_wmtmprt, jmin_wmtmprt, jmax_wmtmprt, kmin_wmtmprt, kmax_wmtmprt;

extern PetscInt NumberOfTurbines; //xyang
extern PetscInt NumberOfIBDelta; //xyang
extern PetscInt rotatewt; //xyang
extern PetscInt IB_delta; // xyang
extern PetscInt NumIBPerLoc;

extern PetscInt IB_wm;
extern PetscInt IB_wmtmprt;
extern PetscInt i_periodicIB, j_periodicIB, k_periodicIB;
extern PetscInt ApproxBC_wm;
extern PetscInt Shear_wm;
extern PetscInt Force_wm;
extern PetscInt infRe;

extern double xmin,xmax,ymin,ymax,zmin,zmax; // the range of domain used for moving frame with wind turbines

extern PetscInt rotor_model;  // xyang 12-7-2010 1: actuator disk model 2, 3: actuator line model
extern PetscReal indf_ax; // xyang 12-16-2010
extern PetscInt surface_p_out; //xyang
extern PetscInt num_blade; // xyang 1-21-2011
extern PetscInt num_foiltype; // The number of airfoil types used in the blade xyang 03-18-2011
extern PetscReal dh1_wm; // the first off-wall grid spacing for wall model 4-11-2011
extern PetscReal dhratio_wm; // the ratio of grid spacing in wall model 
extern PetscReal reflength_wt;  
extern PetscReal reflength_IBDelta; 
extern PetscInt temperature; 
extern PetscInt deltafunc;
extern PetscInt add_fluctuations;
extern PetscInt add_fluctuations_tmprt;
extern PetscReal tipspeedratio;
extern PetscReal r_rotor;
extern PetscReal r_nacelle;

extern PetscInt les_prt;

extern PetscInt AL_Noslip;
extern PetscReal u_frame, v_frame, w_frame;
extern PetscInt MoveFrame;

extern PetscInt ii_periodicWT, jj_periodicWT, kk_periodicWT; // periodic WT, a row/column of ghost wind turbines needs to be added
extern PetscReal Sx_WT, Sy_WT, Sz_WT;
extern PetscInt Nx_WT, Ny_WT, Nz_WT;

// Idealized water wave
extern PetscReal a_iww, lamda_iww, C_iww;

extern char path_inflow[256];

extern PetscReal prt_eps;

extern PetscInt New_wallmodel;

extern PetscInt tisteps;

extern PetscReal les_eps;

extern PetscReal r_nacell;

extern PetscErrorCode TECIOOut_rhs_da(UserCtx *user, Vec Rhs, char fname[80]);	

extern void Compute_Force_SpongeLayer(UserCtx *user);


extern PetscInt SpongeLayer;
extern PetscInt SpongeDistance;

extern PetscInt MoveCylinderTest;

double DDelta(double p, double dx);

extern PetscErrorCode disk_read_ucd(IBMNodes *ibm, int ibi, FSInfo *fsi, int OneDmZ, char fname[80]);

extern PetscInt inflow_levelset;


extern PetscInt sublevel;
extern PetscReal smoothlevel;
extern void Compute_Surface_Curv(UserCtx *user);


extern void Compute_Force_Subgridlevel(UserCtx *user, Vec Force_sublevel);
double dH_1 (double p, double dx);	// its derivative


extern void wallmodel_s( double nu, double sb, double sc, Cmpnts Uc, Cmpnts *Ub,  Cmpnts Ua, PetscInt bctype, double ks,double nx, double ny, double nz, double *tau_w, PetscReal *ustar, double dpdx, double dpdy, double dpdz,double *nut_2sb, double nut_c);


extern PetscReal dfunc_wd;

extern PetscErrorCode Connectivity_ib( UserCtx * user, IBMNodes * ibm );

PetscErrorCode ibm_read_xpatch(IBMNodes *ibm, PetscInt ibi);

PetscErrorCode Scale_InitFlow(UserCtx *user);

extern PetscInt ExtractPoints_BG;
extern PetscInt ExtractPoints_IB;
extern int NumberOfEP_BG;
extern int NumberOfEP_IB;

extern PetscInt FixTipSpeedRatio;


extern void wallfunction_s( double nu, double sb, double sc, Cmpnts Uc, Cmpnts *Ub,  Cmpnts Ua,double ks,double nx, double ny, double nz, PetscReal *ustar, double dpdx, double dpdy, double dpdz);
extern void wall_function_s (double nu, double ks, double sc, double sb, Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, double nx, double ny, double nz);

extern PetscErrorCode preprocess_ib(UserCtx *user, IBMNodes *ibm);

extern PetscReal dt_bed;
extern PetscInt bed_avg;
extern PetscReal particle_D50;
extern PetscReal particle_dens;
extern PetscInt smooth_bed;
extern PetscInt particlevel_model;
extern PetscInt LS_bedConstr;
extern PetscReal Angle_repose;
extern PetscReal bed_porosity;
extern PetscReal sb_bed;
extern PetscInt number_smooth;
extern PetscReal deltab, fdepth, U_bulk_vel;
extern PetscReal rs_bed;
extern void Levelset_smooth(UserCtx *user); 
// add end (xiaolei)
//
// sediment begin
//ali added 2 June 2011
PetscErrorCode Scour(UserCtx * user, IBMNodes * ibm, PetscInt tistart, int ti, int itr_sc);
PetscErrorCode ibm_intp_pj_centroid(UserCtx * user, IBMNodes  * ibm);
PetscErrorCode Projecting_new(UserCtx * user, IBMNodes  * ibm);
PetscErrorCode flow_variables_ref_level_new(UserCtx * user, IBMNodes  * ibm);
PetscErrorCode flow_variables_ref_level(UserCtx * user, IBMNodes  * ibm);
PetscErrorCode Projecting(UserCtx * user, IBMNodes  * ibm);
PetscErrorCode Finalizing_Projecting(UserCtx * user, IBMNodes  * ibm);
PetscErrorCode sediment_variables_projection(UserCtx * user, IBMNodes  * ibm);
PetscErrorCode bed_concentration(UserCtx * user, IBMNodes  * ibm);
PetscErrorCode sediment_flux(UserCtx * user, IBMNodes  * ibm);
PetscErrorCode outlet_sediment_flux_bend(UserCtx * user, IBMNodes  * ibm);
PetscErrorCode outlet_sediment_flux_contra(UserCtx * user, IBMNodes  * ibm);
PetscErrorCode check_correct_new_elev(UserCtx * user, IBMNodes  * ibm);
PetscErrorCode read_bed_cell_depth(IBMNodes  * ibm, PetscInt tistart);

extern PetscInt sediment, input_ib_depth, conv_diff, SuspendedParticles, mobile_bed, projection_method, periodic_morpho, density_current;
extern PetscReal k_ss, w_s, Nseg, U_Bulk;
extern double inlet_sediment_flux;
extern PetscInt LiveBed;
extern PetscInt RigidBed;
extern PetscInt zero_grad;
extern PetscInt bed_roughness;

PetscErrorCode calc_ibm_volumeFlux(IBMNodes *ibm, PetscReal delti, PetscReal *VolumeFlux);

PetscErrorCode ibm_read_ucd_sedi(IBMNodes *ibm, PetscInt ibi);
void wall_function_roughness_a (UserCtx *user, double ks, double sc, double sb, Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, double nx, double ny, double nz);


PetscErrorCode	BC_elevation( UserCtx * user, IBMNodes * ibm );  // xiaolei SEDI
PetscErrorCode	BC_delz( UserCtx * user, IBMNodes * ibm );  // xiaolei SEDI
PetscErrorCode bedvel_model( UserCtx * user, IBMNodes * ibm );
// end of added
// sediment end 

 #include "inline.h"

/*
static double Butcher_ESDIRK2_A[4][4] = {
		{ 0., 0., 0., 0. },
		{ 1767732205903./4055673282236., 1767732205903./4055673282236., 0., 0. },
		{ 2746238789719./10658868560708., -640167445237./6845629431997., 1767732205903./4055673282236., 0. },
		{ 1471266399579./7840856788654., -4482444167858./7529755066697., 11266239266428./11593286722821., 1767732205903./4055673282236. },
};

static double Butcher_ESDIRK2_B[4] = { 1471266399579./7840856788654., -4482444167858./7529755066697., 11266239266428./11593286722821., 1767732205903./4055673282236. };
	
static double Butcher_ESDIRK4_A[6][6] = {
		{ 0., 0., 0., 0., 0., 0. },
		{ 1./4., 1./4., 0., 0., 0., 0. } , 
		{ 8611./62500., -1743./31250., 1./4., 0., 0., 0. }, 
		{ 5012029./34652500., -654441./2922500., 174375./388108., 1./4., 0., 0. },
		{ 15267082809./155376265600., -71443401./120774400., 730878875./902184768., 2285395./8070912., 1./4., 0.},
		{ 82889./524892., 0., 15625./83664., 69875./102672., -2260./8211., 1./4. },
};

static double Butcher_ARK4_A[6][6] = {
		{ 0., 0., 0., 0., 0., 0. },
		{ 1./2., 0., 0., 0., 0., 0. },
		{ 13861./62500., 6889./62500., 0., 0., 0., 0. },
		{ -116923316275./2393684061468., -2731218467317./15368042101831., 9408046702089./11113171139209., 0., 0., 0. },
		{ -451086348788./2902428689909., -2682348792572./7519795681897., 12662868775082./11960479115383., 3355817975965./11060851509271., 0., 0. },
		{ 647845179188./3216320057751., 73281519250./8382639484533., 552539513391./3454668386233., 3354512671639./8306763924573., 4040./17871., 0.},
};

static double Butcher_ESDIRK4_B[6] = { 82889./524892., 0., 15625./83664., 69875./102672., -2260./8211., 1./4. };
*/
#endif
#endif
