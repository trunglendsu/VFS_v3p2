/*Modeling the rotor blades */

#include "variables.h"

extern	double dfunc_2h(double r);
extern	double dfunc_4h(double r);
extern  double dfunc_s3h(double r);

extern	PetscErrorCode Calc_U_lagr(UserCtx *user, IBMNodes *ibm); 
extern	PetscInt ti;
extern 	PetscReal indf_ax;
double  coef_cr1 = 2.5, coef_cr2 = 1.0;


PetscErrorCode Pre_process(UserCtx *user, IBMNodes *ibm)
{
  DA              da = user->da, fda = user->fda;
  DALocalInfo     info;
  PetscInt        xs, xe, ys, ye, zs, ze; // Local grid information
  PetscInt        mx, my, mz; // Dimensions in three directions
  PetscInt        i, j, k, l, ibi;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;
  PetscReal	xb_min, xb_max, yb_min, yb_max, zb_min, zb_max;

  Cmpnts	***coor;

  Vec		Coor;

  DAGetLocalInfo(da, &info);
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


  DAGetGhostedCoordinates(da, &Coor);
  DAVecGetArray(fda, Coor, &coor);

  for (ibi=0; ibi<NumberOfTurbines; ibi++) {
    xb_min = 1.0e6; xb_max= -1.0e6;        
    yb_min = 1.0e6; yb_max= -1.0e6;
    zb_min = 1.0e6; zb_max= -1.0e6;

    for (l=0; l<ibm[ibi].n_elmt; l++) {
      xb_min = PetscMin(xb_min, ibm[ibi].cent_x[l]);
      xb_max = PetscMax(xb_max, ibm[ibi].cent_x[l]);

      yb_min = PetscMin(yb_min, ibm[ibi].cent_y[l]);
      yb_max = PetscMax(yb_max, ibm[ibi].cent_y[l]);

      zb_min = PetscMin(zb_min, ibm[ibi].cent_z[l]);
      zb_max = PetscMax(zb_max, ibm[ibi].cent_z[l]);
    }

/*    xb_min = xb_min - 2.0*fabs(coor[lzs][lys][lxs+1].x - coor[lzs][lys][lxs].x);
    xb_max = xb_max + 2.0*fabs(coor[lzs][lys][lxs+1].x - coor[lzs][lys][lxs].x);

    yb_min = yb_min - 2.0*fabs(coor[lzs][lys+1][lxs].y - coor[lzs][lys][lxs].y);
    yb_max = yb_max + 2.0*fabs(coor[lzs][lys+1][lxs].y - coor[lzs][lys][lxs].y);

    zb_min = zb_min - 2.0*fabs(coor[lzs+1][lys][lxs].z - coor[lzs][lys][lxs].z);
    zb_max = zb_max + 2.0*fabs(coor[lzs+1][lys][lxs].z - coor[lzs][lys][lxs].z);
*/
//    PetscPrintf(PETSC_COMM_SELF, "**** the xb_min, xb_max: %le %le \n", xb_min, xb_max);
//    PetscPrintf(PETSC_COMM_SELF, "**** the yb_min, yb_max: %le %le \n", yb_min, yb_max);
//    PetscPrintf(PETSC_COMM_SELF, "**** the zb_min, zb_max: %le %le \n", zb_min, zb_max);

    j = lys;
    k = lzs;
    ibm[ibi].i_min = lxe; ibm[ibi].i_max = lxs;
    for (i=lxs; i<lxe; i++) {
      if (xb_min >= coor[k][j][i-1].x && xb_min <= coor[k][j][i].x) ibm[ibi].i_min = i - 3;
      if (xb_max >= coor[k][j][i-1].x && xb_max <= coor[k][j][i].x) ibm[ibi].i_max = i + 3;
    }
    if (ibm[ibi].i_min != lxe && ibm[ibi].i_max == lxs) ibm[ibi].i_max = lxe;
    if (ibm[ibi].i_max != lxs && ibm[ibi].i_min == lxe) ibm[ibi].i_min = lxs;
    if (ibm[ibi].i_max == lxs && ibm[ibi].i_min == lxe) {
      ibm[ibi].i_min = lxs; ibm[ibi].i_max = lxs;
    }
    ibm[ibi].i_min = PetscMax(ibm[ibi].i_min, lxs);
    ibm[ibi].i_max = PetscMin(ibm[ibi].i_max, lxe);

    i = lxs;
    k = lzs;
    ibm[ibi].j_min = lye; ibm[ibi].j_max = lys;
    for (j=lys; j<lye; j++) {
      if (yb_min >= coor[k][j-1][i].y && yb_min <= coor[k][j][i].y) ibm[ibi].j_min = j - 3;
      if (yb_max >= coor[k][j-1][i].y && yb_max <= coor[k][j][i].y) ibm[ibi].j_max = j + 3;
    }
    if (ibm[ibi].j_min != lye && ibm[ibi].j_max == lys) ibm[ibi].j_max = lye;
    if (ibm[ibi].j_max != lys && ibm[ibi].j_min == lye) ibm[ibi].j_min = lys;
    if (ibm[ibi].j_max == lys && ibm[ibi].j_min == lye) {
      ibm[ibi].j_min = lys; ibm[ibi].j_max = lys;
    }

    ibm[ibi].j_min = PetscMax(ibm[ibi].j_min, lys);    
    ibm[ibi].j_max = PetscMin(ibm[ibi].j_max, lye);


    j = lys;
    i = lxs;
    ibm[ibi].k_min = lze; ibm[ibi].k_max = lzs;
    for (k=lzs; k<lze; k++) {
      if (zb_min >= coor[k-1][j][i].z && zb_min <= coor[k][j][i].z) ibm[ibi].k_min = k - 3;
      if (zb_max >= coor[k-1][j][i].z && zb_max <= coor[k][j][i].z) ibm[ibi].k_max = k + 3;
    }

    if (ibm[ibi].k_min != lze && ibm[ibi].k_max == lzs) ibm[ibi].k_max = lze;
    if (ibm[ibi].k_max != lzs && ibm[ibi].k_min == lze) ibm[ibi].k_min = lzs;
    if (ibm[ibi].k_max == lys && ibm[ibi].k_min == lye) {
      ibm[ibi].k_min = lys; ibm[ibi].k_max = lys;
    }

    ibm[ibi].k_min = PetscMax(ibm[ibi].k_min, lzs);   
    ibm[ibi].k_max = PetscMin(ibm[ibi].k_max, lze);

//    PetscPrintf(PETSC_COMM_SELF, "**** the i_min, i_max: %i %i \n", ibm[ibi].i_min, ibm[ibi].i_max);
//    PetscPrintf(PETSC_COMM_SELF, "**** the j_min, j_max: %i %i \n", ibm[ibi].j_min, ibm[ibi].j_max);
//    PetscPrintf(PETSC_COMM_SELF, "**** the k_min, k_max: %i %i \n", ibm[ibi].k_min, ibm[ibi].k_max);


  }


  DAVecRestoreArray(fda, Coor, &coor);

  return(0);

}


/* Read information of ACL, like chord, twist angle, lift and drag coefficients */
PetscErrorCode airfoil_ACL(ACL *acl, IBMNodes *ibm,  FSInfo *fsi)
{

	PetscInt	rank;
  	PetscInt	n_CD, n_CL;
	PetscInt	ifoil, i;

  	char string[128];

  	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  	if(!rank) { // root processor read in the data
    		FILE *fd;
    		PetscPrintf(PETSC_COMM_SELF, "READ airfoil data for ACL\n");
    		char filen[80]; 

		for(ifoil=0; ifoil<num_foiltype; ifoil++) {

    			sprintf(filen,"%s/FOIL%2.2d" , path, ifoil);
    			fd = fopen(filen, "r"); if (!fd) SETERRQ(1, "Cannot open airfoil data file");
			if (fd) {
      				fgets(string, 128, fd);
      				fgets(string, 128, fd);

      				fscanf(fd, "%i ",&(acl[ifoil].num_AC));

      				MPI_Bcast(&(acl[ifoil].num_AC), 1, MPI_INT, 0, PETSC_COMM_WORLD);

                        	PetscMalloc((acl[ifoil].num_AC)*sizeof(PetscReal), &(acl[ifoil].r_AC));  
                        	PetscMalloc((acl[ifoil].num_AC)*sizeof(PetscReal), &(acl[ifoil].angle_pitchInp));  
                        	PetscMalloc((acl[ifoil].num_AC)*sizeof(PetscReal), &(acl[ifoil].chord_bladeInp)); 

			 	for(i=0; i<acl[ifoil].num_AC; i++) {
					fscanf(fd, "%le %le %le", &(acl[ifoil].r_AC[i]), &(acl[ifoil].angle_pitchInp[i]), &(acl[ifoil].chord_bladeInp[i]));

					acl[ifoil].r_AC[i] = acl[ifoil].r_AC[i] / reflength_wt;
					acl[ifoil].chord_bladeInp[i] = acl[ifoil].chord_bladeInp[i] / reflength_wt;

				} 

				acl[ifoil].r_beg = acl[ifoil].r_AC[0];
				acl[ifoil].r_end = acl[ifoil].r_AC[acl[ifoil].num_AC-1];

				MPI_Bcast(acl[ifoil].r_AC, acl[ifoil].num_AC, MPIU_REAL, 0, PETSC_COMM_WORLD); 
				MPI_Bcast(acl[ifoil].angle_pitchInp, acl[ifoil].num_AC, MPIU_REAL, 0, PETSC_COMM_WORLD); 
				MPI_Bcast(acl[ifoil].chord_bladeInp, acl[ifoil].num_AC, MPIU_REAL, 0, PETSC_COMM_WORLD); 
				MPI_Bcast(&(acl[ifoil].r_beg), 1, MPIU_REAL, 0, PETSC_COMM_WORLD); 
				MPI_Bcast(&(acl[ifoil].r_end), 1, MPIU_REAL, 0, PETSC_COMM_WORLD); 

			}

                        sprintf(filen,"%s/CD%2.2d" , path, ifoil);
                        fd = fopen(filen, "r"); if (!fd) SETERRQ(1, "Cannot open CD file");
                        if (fd) {
                                fgets(string, 128, fd);
                                fgets(string, 128, fd);

                                fscanf(fd, "%i ",&(acl[ifoil].num_CD));

                                MPI_Bcast(&(acl[ifoil].num_CD), 1, MPI_INT, 0, PETSC_COMM_WORLD);

                                PetscMalloc((acl[ifoil].num_CD)*sizeof(PetscReal), &(acl[ifoil].ang_CD));
                                PetscMalloc((acl[ifoil].num_CD)*sizeof(PetscReal), &(acl[ifoil].CDInp));

                                for(i=0; i<acl[ifoil].num_CD; i++) {
                                        fscanf(fd, "%le %le", &(acl[ifoil].ang_CD[i]), &(acl[ifoil].CDInp[i]));
                                }

                                MPI_Bcast(acl[ifoil].ang_CD, acl[ifoil].num_CD, MPIU_REAL, 0, PETSC_COMM_WORLD);
                                MPI_Bcast(acl[ifoil].CDInp, acl[ifoil].num_CD, MPIU_REAL, 0, PETSC_COMM_WORLD);

                        }

                        sprintf(filen,"%s/CL%2.2d" , path, ifoil);
                        fd = fopen(filen, "r"); if (!fd) SETERRQ(1, "Cannot CL file");
                        if (fd) {
                                fgets(string, 128, fd);
                                fgets(string, 128, fd);

                                fscanf(fd, "%i ",&(acl[ifoil].num_CL));

                                MPI_Bcast(&(acl[ifoil].num_CL), 1, MPI_INT, 0, PETSC_COMM_WORLD);

                                PetscMalloc((acl[ifoil].num_CL)*sizeof(PetscReal), &(acl[ifoil].ang_CL));
                                PetscMalloc((acl[ifoil].num_CL)*sizeof(PetscReal), &(acl[ifoil].CLInp));

                                for(i=0; i<acl[ifoil].num_CL; i++) {
                                        fscanf(fd, "%le %le", &(acl[ifoil].ang_CL[i]), &(acl[ifoil].CLInp[i]));
                                }

                                MPI_Bcast(acl[ifoil].ang_CL, acl[ifoil].num_CL, MPIU_REAL, 0, PETSC_COMM_WORLD);
                                MPI_Bcast(acl[ifoil].CLInp, acl[ifoil].num_CL, MPIU_REAL, 0, PETSC_COMM_WORLD);

                        }



		}

 
 
	} else if(rank) {


		for(ifoil=0; ifoil<num_foiltype; ifoil++) {

			// angle of pitch , chord length of blade
      			MPI_Bcast(&(acl[ifoil].num_AC), 1, MPI_INT, 0, PETSC_COMM_WORLD);

                        PetscMalloc((acl[ifoil].num_AC)*sizeof(PetscReal), &(acl[ifoil].r_AC));  
                        PetscMalloc((acl[ifoil].num_AC)*sizeof(PetscReal), &(acl[ifoil].angle_pitchInp));  
                        PetscMalloc((acl[ifoil].num_AC)*sizeof(PetscReal), &(acl[ifoil].chord_bladeInp)); 

			
			MPI_Bcast(acl[ifoil].r_AC, acl[ifoil].num_AC, MPIU_REAL, 0, PETSC_COMM_WORLD); 
			MPI_Bcast(acl[ifoil].angle_pitchInp, acl[ifoil].num_AC, MPIU_REAL, 0, PETSC_COMM_WORLD); 
			MPI_Bcast(acl[ifoil].chord_bladeInp, acl[ifoil].num_AC, MPIU_REAL, 0, PETSC_COMM_WORLD); 
			MPI_Bcast(&(acl[ifoil].r_beg), 1, MPIU_REAL, 0, PETSC_COMM_WORLD); 
			MPI_Bcast(&(acl[ifoil].r_end), 1, MPIU_REAL, 0, PETSC_COMM_WORLD); 

			// drag coefficients
                        MPI_Bcast(&(acl[ifoil].num_CD), 1, MPI_INT, 0, PETSC_COMM_WORLD);

                        PetscMalloc((acl[ifoil].num_CD)*sizeof(PetscReal), &(acl[ifoil].ang_CD));
                        PetscMalloc((acl[ifoil].num_CD)*sizeof(PetscReal), &(acl[ifoil].CDInp));

                        MPI_Bcast(acl[ifoil].ang_CD, acl[ifoil].num_CD, MPIU_REAL, 0, PETSC_COMM_WORLD);
                        MPI_Bcast(acl[ifoil].CDInp, acl[ifoil].num_CD, MPIU_REAL, 0, PETSC_COMM_WORLD);

			// left coefficient
                        MPI_Bcast(&(acl[ifoil].num_CL), 1, MPI_INT, 0, PETSC_COMM_WORLD);

                        PetscMalloc((acl[ifoil].num_CL)*sizeof(PetscReal), &(acl[ifoil].ang_CL));
                        PetscMalloc((acl[ifoil].num_CL)*sizeof(PetscReal), &(acl[ifoil].CLInp));

                        MPI_Bcast(acl[ifoil].ang_CL, acl[ifoil].num_CL, MPIU_REAL, 0, PETSC_COMM_WORLD);
                        MPI_Bcast(acl[ifoil].CLInp, acl[ifoil].num_CL, MPIU_REAL, 0, PETSC_COMM_WORLD);

		}

	}

	int ibi, j;
	PetscReal r, fac1, fac2;
	double rr;
	for (ibi=0;ibi<NumberOfTurbines;ibi++) {
		for (i=0;i<ibm[ibi].n_elmt;i++) {
			r = sqrt( pow((ibm[ibi].cent_x[i]-fsi[ibi].x_c),2.0) + pow((ibm[ibi].cent_y[i]-fsi[ibi].y_c),2.0) + pow((ibm[ibi].cent_z[i]-fsi[ibi].z_c),2.0));	
			rr = sqrt( pow((ibm[ibi].cent_x[i]),2.0) + pow((ibm[ibi].cent_y[i]),2.0) + pow((ibm[ibi].cent_z[i]),2.0));	
//			if (!rank) printf("#### the center %d %le %le %le \n", i, fsi[ibi].x_c, fsi[ibi].y_c, fsi[ibi].z_c );

			for(ifoil=0; ifoil<num_foiltype; ifoil++) {

//			if (!rank) printf("#### relative r position %d %d %le %le %le %le \n", i, ifoil, acl[ifoil].r_beg, r, rr, acl[ifoil].r_end );
				if( r>=acl[ifoil].r_beg && r<=acl[ifoil].r_end ) {
					
					for(j=0; j<acl[ifoil].num_AC-1; j++) {
						if (r>=acl[ifoil].r_AC[j] && r<=acl[ifoil].r_AC[j+1]) {
							fac1 = (acl[ifoil].r_AC[j+1]-r)/(acl[ifoil].r_AC[j+1]-acl[ifoil].r_AC[j]);
							fac2 = (-acl[ifoil].r_AC[j]+r)/(acl[ifoil].r_AC[j+1]-acl[ifoil].r_AC[j]);
							ibm[ibi].angle_pitch[i]	= fac1*acl[ifoil].angle_pitchInp[j] + fac2*acl[ifoil].angle_pitchInp[j+1];
							ibm[ibi].chord_blade[i]	= fac1*acl[ifoil].chord_bladeInp[j] + fac2*acl[ifoil].chord_bladeInp[j+1];
						}
					}
				}
			}

//			if (!rank) printf("#### the chord and pitch angle %d %le %le \n", i,  ibm[ibi].chord_blade[i], ibm[ibi].angle_pitch[i]);
		}

	}


 

	return(0);






}





PetscErrorCode Calc_F_lagr(UserCtx *user, IBMNodes *ibm, FSInfo *fsi)
{

  PetscInt	l, ibi;
  double	pi = 3.141592653589793, a = 0.25;
  double 	C_T;  // = 4.0 / 3.0;
  double 	U_ref, A_sum;
  PetscReal	sign;
  
  C_T = 4.0 * indf_ax * (1-indf_ax);

  C_T = C_T / ( (1.0 - indf_ax)* (1.0 - indf_ax) );

  Calc_U_lagr(user, ibm);
  
  for (ibi=0; ibi<NumberOfTurbines; ibi++) {
    U_ref = 0.0;
    A_sum = 0.0;
    for (l=0; l<ibm[ibi].n_elmt; l++) {
      U_ref += ibm[ibi].U_lagr_z[l] * ibm[ibi].dA[l];
 //     PetscPrintf(PETSC_COMM_WORLD, "the U_lagr_z: %le \n", ibm[ibi].U_lagr_z[l] );
      A_sum += ibm[ibi].dA[l];
    }
    U_ref /= A_sum;

//    U_ref = 0.70;
    PetscPrintf(PETSC_COMM_WORLD, "**** The U_ref at %i th body: %le\n", ibi, U_ref);
//    U_ref = 1.0;

    sign = U_ref / fabs(U_ref);
    for (l=0; l<ibm[ibi].n_elmt; l++) {
      ibm[ibi].F_lagr_x[l] = 0.0;
      ibm[ibi].F_lagr_y[l] = 0.0;
      ibm[ibi].F_lagr_z[l] = -0.5 * C_T * (U_ref * U_ref) * sign; //* (1.0-exp(-(double)(ti-1)));    
    }
  }

 return(0);

}


PetscErrorCode Calc_F_lagr_noslip(UserCtx *user, IBMNodes *ibm, FSInfo *fsi)
{

  	PetscInt	l, ibi;
  	int 		nv1, nv2, nv3;
	double 		r1, r2, r3, rr, rx, ry, rz;


        PetscInt rank=0;
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);


  	Calc_U_lagr(user, ibm);
  
  	for (ibi=0; ibi<NumberOfBodies; ibi++) {
  	for (l=0; l<ibm[ibi].n_elmt; l++) {
             	nv1 = ibm[ibi].nv1[l];
            	nv2 = ibm[ibi].nv2[l];
              	nv3 = ibm[ibi].nv3[l];

               	rx = ibm[ibi].x_bp[nv2] - ibm[ibi].x_bp[nv1];
             	ry = ibm[ibi].y_bp[nv2] - ibm[ibi].y_bp[nv1];
             	rz = ibm[ibi].z_bp[nv2] - ibm[ibi].z_bp[nv1];

              	r1 = sqrt(rx*rx + ry*ry + rz*rz );
		
                rx = ibm[ibi].x_bp[nv3] - ibm[ibi].x_bp[nv2];
                ry = ibm[ibi].y_bp[nv3] - ibm[ibi].y_bp[nv2];
                rz = ibm[ibi].z_bp[nv3] - ibm[ibi].z_bp[nv2];

                r2 = sqrt(rx*rx + ry*ry + rz*rz );

                rx = ibm[ibi].x_bp[nv1] - ibm[ibi].x_bp[nv3];
                ry = ibm[ibi].y_bp[nv1] - ibm[ibi].y_bp[nv3];
                rz = ibm[ibi].z_bp[nv1] - ibm[ibi].z_bp[nv3];

                r3 = sqrt(rx*rx + ry*ry + rz*rz );

		rr = (r1 + r2 + r3)/3.0;

      		ibm[ibi].F_lagr_x[l] = rr * (ibm[ibi].u[l].x - ibm[ibi].U_lagr_x[l])/user->dt;
      		ibm[ibi].F_lagr_y[l] = rr * (ibm[ibi].u[l].y - ibm[ibi].U_lagr_y[l])/user->dt;
      		ibm[ibi].F_lagr_z[l] = rr * (ibm[ibi].u[l].z - ibm[ibi].U_lagr_z[l])/user->dt;    

//		if (l==0 && !rank) printf("The Lagrange Force %d %le %le %le \n", l, ibm[ibi].F_lagr_z[l], ibm[ibi].u[l].z,  ibm[ibi].U_lagr_z[l]);
    	}
  	}

 return(0);

}


// calculate force at largrangian points using actuator line model
PetscErrorCode Calc_F_lagr_ACL(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, ACL *acl)
{

  	PetscInt      	l, ibi, j;
  	double		pi = 3.141592653589793;
  	double		A_sum, U_ref, rr, rx, ry, rz, tmp;
	Cmpnts	      	U_b, n_relV, n_L, n_blade;	
	int		nv1, nv2, ifoil;
	PetscReal 	fac1, fac2, r;
	// n_relV: the direction of relative velocity
	// n_L: the direction of Lift

	int istall;

  	PetscInt rank=0;
  	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);


	Calc_U_lagr(user, ibm);



  	for (ibi=0; ibi<NumberOfTurbines; ibi++) {

		
		U_b.z = 0.0;
    		A_sum = 0.0;
    		for (l=0; l<ibm[ibi].n_elmt; l++) {
      			U_b.z += ibm[ibi].U_lagr_z[l] * ibm[ibi].dA[l];
      			A_sum += ibm[ibi].dA[l];

//			printf("dA %le \n", ibm[ibi].dA[l]);
    		}
   		U_b.z /= A_sum;

//		if (rotor_model == 3) U_b.z = ibm[ibi].U_ref;
		// for one turbine case, uniform inflow


    		PetscPrintf(PETSC_COMM_WORLD, "**** The U_ref in z-dir at %i th body: %le %le\n", ibi, U_b.z, A_sum);
		
    		
    		for (l=0; l<ibm[ibi].n_elmt; l++) {

			nv1 = ibm[ibi].nv1[l];
			nv2 = ibm[ibi].nv2[l];


			U_b.x = ibm[ibi].U_lagr_x[l] - 0.5*(ibm[ibi].u[nv1].x + ibm[ibi].u[nv2].x);
			U_b.y = ibm[ibi].U_lagr_y[l] - 0.5*(ibm[ibi].u[nv1].y + ibm[ibi].u[nv2].y);
			U_b.z = ibm[ibi].U_lagr_z[l];


			if (rotor_model == 3) U_b.z = ibm[ibi].U_ref;

//                        U_b.x = - 0.5*(ibm[ibi].u[nv1].x + ibm[ibi].u[nv2].x);
//                        U_b.y = - 0.5*(ibm[ibi].u[nv1].y + ibm[ibi].u[nv2].y);

			
			U_ref = sqrt(U_b.x*U_b.x + U_b.y*U_b.y + U_b.z*U_b.z);

			n_relV.x = U_b.x/U_ref; n_relV.y = U_b.y/U_ref; n_relV.z = U_b.z/U_ref; 

//			if (l==1) printf("the ref vel: %le \n", U_ref);

			tmp = n_relV.z;

			ibm[ibi].angle_attack[l] = 90.0 - acos(tmp) * 180.0 / pi - ibm[ibi].angle_pitch[l];
			
			if (!rank && ( ibm[ibi].angle_attack[l]<-45.0 || ibm[ibi].angle_attack[l]>45.0 ) ) PetscPrintf(PETSC_COMM_WORLD, "**** AOA > 45 or < -45: %d %le \n", l, ibm[ibi].angle_attack[l]);

 
			rx = ibm[ibi].x_bp[nv2] - ibm[ibi].x_bp[nv1]; 
			ry = ibm[ibi].y_bp[nv2] - ibm[ibi].y_bp[nv1]; 
			rz = ibm[ibi].z_bp[nv2] - ibm[ibi].z_bp[nv1];

			rr = sqrt(rx*rx + ry*ry + rz*rz );

			n_blade.x = rx/rr; n_blade.y = ry/rr; n_blade.z = rz/rr; 

        		n_L.x = n_relV.y*n_blade.z - n_relV.z*n_blade.y; 
			n_L.y = n_relV.z*n_blade.x - n_relV.x*n_blade.z; 
			n_L.z = n_relV.x*n_blade.y - n_relV.y*n_blade.x;
		
			tmp = U_b.z * n_L.z;
			
			if (tmp < 0.0) n_L.x = -n_L.x, n_L.y = -n_L.y, n_L.z = -n_L.z; 				 


                        r = sqrt( pow((ibm[ibi].cent_x[l]-fsi[ibi].x_c),2) + pow((ibm[ibi].cent_y[l]-fsi[ibi].y_c),2) + pow((ibm[ibi].cent_z[l]-fsi[ibi].z_c),2));

                        for(ifoil=0; ifoil<num_foiltype; ifoil++) {
				
//				if (!rank) printf("why zero %d %d %le %le %le \n",  l, ifoil, acl[ifoil].r_beg, r,  acl[ifoil].r_end  );				

                                if( r>=acl[ifoil].r_beg && r<=acl[ifoil].r_end ) {

                                        for(j=0; j<acl[ifoil].num_CD-1; j++) {
                                                if (ibm[ibi].angle_attack[l]>=acl[ifoil].ang_CD[j] && ibm[ibi].angle_attack[l]<=acl[ifoil].ang_CD[j+1]) {
                                                        fac1 = (acl[ifoil].ang_CD[j+1]-ibm[ibi].angle_attack[l])/(acl[ifoil].ang_CD[j+1]-acl[ifoil].ang_CD[j]);
                                                        fac2 = (-acl[ifoil].ang_CD[j]+ibm[ibi].angle_attack[l])/(acl[ifoil].ang_CD[j+1]-acl[ifoil].ang_CD[j]);
                                                        ibm[ibi].CD[l] = fac1*acl[ifoil].CDInp[j] + fac2*acl[ifoil].CDInp[j+1];
                                                }

						if (15.0>=acl[ifoil].ang_CD[j] && 15.0<=acl[ifoil].ang_CD[j+1]) {
							istall = j;						
						}
                                        }
					
//					if (ibm[ibi].angle_attack[l]<acl[ifoil].ang_CD[0]) ibm[ibi].CD[l] = 1.2;

                                        if (ibm[ibi].angle_attack[l]<acl[ifoil].ang_CD[0] && ibm[ibi].angle_attack[l]>-45.0) ibm[ibi].CD[l] = (-1.2 - acl[ifoil].CDInp[0]) * (ibm[ibi].angle_attack[l] - acl[ifoil].ang_CD[0]) / (-45.0-acl[ifoil].ang_CD[0]) + acl[ifoil].CDInp[0];
                                        if ( ibm[ibi].angle_attack[l]<=-45.0 && ibm[ibi].angle_attack[l]>=-90.0 ) ibm[ibi].CD[l] = -1.2;

					if ( ibm[ibi].angle_attack[l] < -90.0 ) ibm[ibi].CD[l] = -1.2;

 
					if (ibm[ibi].angle_attack[l]>acl[ifoil].ang_CD[acl[ifoil].num_CD-1] && ibm[ibi].angle_attack[l]<45.0) ibm[ibi].CD[l] = (1.2 - acl[ifoil].CDInp[acl[ifoil].num_CD-1]) * (ibm[ibi].angle_attack[l] - acl[ifoil].ang_CD[acl[ifoil].num_CD-1]) / (45.0-acl[ifoil].ang_CD[acl[ifoil].num_CD-1]) + acl[ifoil].CDInp[acl[ifoil].num_CD-1];
					if ( ibm[ibi].angle_attack[l]>=45.0 && ibm[ibi].angle_attack[l]<=90.0 ) ibm[ibi].CD[l] = 1.2;

//					if ( ibm[ibi].angle_attack[l]>acl[ifoil].ang_CD[acl[ifoil].num_CD-1] ) ibm[ibi].CD[l] = 1.2;


                                        for(j=0; j<acl[ifoil].num_CL-1; j++) {
                                                if (ibm[ibi].angle_attack[l]>=acl[ifoil].ang_CL[j] && ibm[ibi].angle_attack[l]<=acl[ifoil].ang_CL[j+1]) {
                                                        fac1 = (acl[ifoil].ang_CL[j+1]-ibm[ibi].angle_attack[l])/(acl[ifoil].ang_CL[j+1]-acl[ifoil].ang_CL[j]);
                                                        fac2 = (-acl[ifoil].ang_CL[j]+ibm[ibi].angle_attack[l])/(acl[ifoil].ang_CL[j+1]-acl[ifoil].ang_CL[j]);
                                                        ibm[ibi].CL[l] = fac1*acl[ifoil].CLInp[j] + fac2*acl[ifoil].CLInp[j+1];
                                                }

                                                if (15.0>=acl[ifoil].ang_CL[j] && 15.0<=acl[ifoil].ang_CL[j+1]) {
                                                        istall = j;
                                                }

                                        }

//                                        if (ibm[ibi].angle_attack[l]<acl[ifoil].ang_CL[0]) ibm[ibi].CL[l] = acl[ifoil].CLInp[0];  // not sure

                                        if (ibm[ibi].angle_attack[l]<acl[ifoil].ang_CL[0] && ibm[ibi].angle_attack[l]>-45.0) ibm[ibi].CL[l] = (-1.05 - acl[ifoil].CLInp[0] ) * (ibm[ibi].angle_attack[l] - acl[ifoil].ang_CL[0] ) / (-45.0-acl[ifoil].ang_CL[0] ) + acl[ifoil].CLInp[0] ;
                                        if ( ibm[ibi].angle_attack[l]<=-45.0 && ibm[ibi].angle_attack[l]>=-90.0 ) ibm[ibi].CL[l] = 1.05*sin(2.0*ibm[ibi].angle_attack[l]*pi/180.0);
					if ( ibm[ibi].angle_attack[l] < -90.0 ) ibm[ibi].CL[l] = 0.0;


                                        if (ibm[ibi].angle_attack[l]>acl[ifoil].ang_CL[acl[ifoil].num_CL-1] && ibm[ibi].angle_attack[l]<45.0) ibm[ibi].CL[l] = (1.05 - acl[ifoil].CLInp[acl[ifoil].num_CL-1] ) * (ibm[ibi].angle_attack[l] - acl[ifoil].ang_CL[acl[ifoil].num_CL-1] ) / (45.0-acl[ifoil].ang_CL[acl[ifoil].num_CL-1] ) + acl[ifoil].CLInp[acl[ifoil].num_CL-1] ;
					if ( ibm[ibi].angle_attack[l]>=45.0 && ibm[ibi].angle_attack[l]<=90.0 ) ibm[ibi].CL[l] = 1.05*sin(2.0*ibm[ibi].angle_attack[l]*pi/180.0);

// I think I read somewhere that a very approximate value for the lift coefficient of a flat plate is

// L=1.05*sin(2*a)


                        // add 3D and rotational effects
					if (ibm[ibi].angle_attack[l]>0.0 && ibm[ibi].angle_attack[l]<20.0)
					{
//                        			ibm[ibi].CD[l] += coef_cr1 * pow( ibm[ibi].chord_blade[l]/r, coef_cr2 ) * ( 0.0 - ibm[ibi].CD[l] );
//                        			ibm[ibi].CL[l] += coef_cr1 * pow( ibm[ibi].chord_blade[l]/r, coef_cr2 ) * ( 2.0*pi*pi*ibm[ibi].angle_attack[l]/180.0 - ibm[ibi].CL[l]);
					} 
					
//					if ( ibm[ibi].CD[l] < 0.0) ibm[ibi].CD[l] = 0.0;
					ibm[ibi].CD[l] = 0.0;
                                }



                        }

			// add 3D and rotational effects
//			ibm[ibi].CD[l] += coef_cr1 * pow( ibm[ibi].chord_blade[l]/r, coef_cr2 ) * ( ibm[ibi].CD[l] - acl[1].CDInp[0] );
//			ibm[ibi].CL[l] += coef_cr1 * pow( ibm[ibi].chord_blade[l]/r, coef_cr2 ) * ( 2.0*pi*pi*fabs(ibm[ibi].angle_attack[l])/180.0 - ibm[ibi].CL[l]);

//			ibm[ibi].CL[l] += coef_cr1 * pow( ibm[ibi].chord_blade[l]/r, coef_cr2 ) * ( 2.0*pi*pi*fabs(ibm[ibi].angle_attack[l])/180.0 - ibm[ibi].CL[l]);

			// force from C_D
      			ibm[ibi].F_lagr_x[l] = -0.5 * U_ref * U_ref * fabs(ibm[ibi].CD[l]) * ibm[ibi].chord_blade[l] * n_relV.x;
      			ibm[ibi].F_lagr_y[l] = -0.5 * U_ref * U_ref * fabs(ibm[ibi].CD[l]) * ibm[ibi].chord_blade[l] * n_relV.y;
      			ibm[ibi].F_lagr_z[l] = -0.5 * U_ref * U_ref * fabs(ibm[ibi].CD[l]) * ibm[ibi].chord_blade[l] * n_relV.z; 
			
//			if (!rank) printf("***** Check Lag CD %d %le %le %le %le %le %le \n", l, U_ref , ibm[ibi].CD[l] , ibm[ibi].chord_blade[l] , n_relV.x, n_relV.y, n_relV.z);
			// add force from C_L
                        ibm[ibi].F_lagr_x[l] += -0.5 * U_ref * U_ref * fabs(ibm[ibi].CL[l]) * ibm[ibi].chord_blade[l] * n_L.x;
                        ibm[ibi].F_lagr_y[l] += -0.5 * U_ref * U_ref * fabs(ibm[ibi].CL[l]) * ibm[ibi].chord_blade[l] * n_L.y;
                        ibm[ibi].F_lagr_z[l] += -0.5 * U_ref * U_ref * fabs(ibm[ibi].CL[l]) * ibm[ibi].chord_blade[l] * n_L.z;

//			if (!rank) printf("***** Check Lag CL %d %le %le %le %le %le %le \n", l, U_ref , ibm[ibi].CL[l] , ibm[ibi].chord_blade[l] , n_L.x, n_L.y, n_L.z);

/* noslip BCs*/
/*
			if ( r < 0.15*acl[num_foiltype-1].r_end ) {

                        	ibm[ibi].F_lagr_x[l] = - 0.0; //rr*rr*ibm[ibi].U_lagr_x[l] / user->dt;
                        	ibm[ibi].F_lagr_y[l] = - 0.0; //rr*rr*ibm[ibi].U_lagr_y[l] / user->dt;
                        	ibm[ibi].F_lagr_z[l] = - rr*rr*ibm[ibi].U_lagr_z[l] / user->dt;

			}
*/
//			printf("**** the Lag force: %le %le %le %le %le \n", ibm[ibi].F_lagr_x[l], ibm[ibi].F_lagr_y[l], ibm[ibi].F_lagr_z[l], rr, ibm[ibi].U_lagr_z[l]  );

    		}
  	}




}

/* Interpolate the velocity at the Lagrangian points */
PetscErrorCode Calc_U_lagr(UserCtx *user, IBMNodes *ibm)
{

  DA              da = user->da, fda = user->fda;
  DALocalInfo     info;
  PetscInt        xs, xe, ys, ye, zs, ze; // Local grid information
  PetscInt        mx, my, mz; // Dimensions in three directions
  PetscInt        i, j, k, l, ibi;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  Cmpnts	***ucat, ***coor, ***csi, ***eta, ***zet;

  Vec		Coor;

  PetscReal	dfunc;

  PetscReal	xc, yc, zc, hx, hy, hz, vol_eul, vol_lagr;

  double r1, r2, r3;

  DAGetLocalInfo(da, &info);
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


  DAGetGhostedCoordinates(da, &Coor);
  DAVecGetArray(fda, Coor, &coor);
  DAVecGetArray(fda, user->lUcat, &ucat);
  DAVecGetArray(fda, user->lCsi,  &csi);
  DAVecGetArray(fda, user->lEta,  &eta);
  DAVecGetArray(fda, user->lZet,  &zet);

  for (ibi=0; ibi<NumberOfTurbines; ibi++) {
    for (l=0; l<ibm[ibi].n_elmt; l++) {
      ibm[ibi].U_lagr_x[l] = 0.0;
      ibm[ibi].U_lagr_y[l] = 0.0;
      ibm[ibi].U_lagr_z[l] = 0.0;
    }
  }


  for (ibi=0; ibi<NumberOfTurbines; ibi++) {

    for (k = ibm[ibi].k_min; k<ibm[ibi].k_max; k++) {
      for (j = ibm[ibi].j_min; j<ibm[ibi].j_max; j++) {
        for (i= ibm[ibi].i_min; i<ibm[ibi].i_max; i++) {
          xc = (coor[k][j][i  ].x + coor[k][j-1][i  ].x + coor[k-1][j][i  ].x + coor[k-1][j-1][i  ].x +
                coor[k][j][i-1].x + coor[k][j-1][i-1].x + coor[k-1][j][i-1].x + coor[k-1][j-1][i-1].x ) * 0.125;
          yc = (coor[k][j  ][i].y + coor[k][j  ][i-1].y + coor[k-1][j  ][i].y + coor[k-1][j  ][i-1].y + 
                coor[k][j-1][i].y + coor[k][j-1][i-1].y + coor[k-1][j-1][i].y + coor[k-1][j-1][i-1].y ) * 0.125;
          zc = (coor[k  ][j][i].z + coor[k  ][j-1][i].z + coor[k  ][j][i-1].z + coor[k  ][j-1][i-1].z + 
                coor[k-1][j][i].z + coor[k-1][j-1][i].z + coor[k-1][j][i-1].z + coor[k-1][j-1][i-1].z ) * 0.125;

          hx = (coor[k][j][i  ].x + coor[k][j-1][i  ].x + coor[k-1][j][i  ].x + coor[k-1][j-1][i  ].x -
                coor[k][j][i-1].x - coor[k][j-1][i-1].x - coor[k-1][j][i-1].x - coor[k-1][j-1][i-1].x ) * 0.25;

          hy = (coor[k][j  ][i].y + coor[k][j  ][i-1].y + coor[k-1][j  ][i].y + coor[k-1][j  ][i-1].y -
                coor[k][j-1][i].y - coor[k][j-1][i-1].y - coor[k-1][j-1][i].y - coor[k-1][j-1][i-1].y ) * 0.25;

          hz = (coor[k  ][j][i].z + coor[k  ][j-1][i].z + coor[k  ][j][i-1].z + coor[k  ][j-1][i-1].z - 
                coor[k-1][j][i].z - coor[k-1][j-1][i].z - coor[k-1][j][i-1].z - coor[k-1][j-1][i-1].z ) * 0.25;

          for (l=0; l<ibm[ibi].n_elmt; l++) {
            r1 = (xc - ibm[ibi].cent_x[l])/hx; r2 = (yc - ibm[ibi].cent_y[l])/hy; r3 = (zc - ibm[ibi].cent_z[l])/hz;
            dfunc =  dfunc_s3h(r1) * dfunc_s3h(r2) * dfunc_s3h(r3);
            ibm[ibi].U_lagr_x[l] += ucat[k][j][i].x * dfunc;
            ibm[ibi].U_lagr_y[l] += ucat[k][j][i].y * dfunc;
            ibm[ibi].U_lagr_z[l] += ucat[k][j][i].z * dfunc;

          }
        }
      }
    }
  }

//  PetscPrintf(PETSC_COMM_SELF, "the U_lagr_z: %le \n", ibm[0].U_lagr_z[0] );
  PetscBarrier(PETSC_NULL);

  PetscReal	u_sum;
  for (ibi=0; ibi<NumberOfTurbines; ibi++) {
    for (l=0; l<ibm[ibi].n_elmt; l++) {
      PetscGlobalSum(&(ibm[ibi].U_lagr_x[l]), &u_sum, PETSC_COMM_WORLD);
      ibm[ibi].U_lagr_x[l] = u_sum;      

      PetscGlobalSum(&(ibm[ibi].U_lagr_y[l]), &u_sum, PETSC_COMM_WORLD);
      ibm[ibi].U_lagr_y[l] = u_sum;

      PetscGlobalSum(&(ibm[ibi].U_lagr_z[l]), &u_sum, PETSC_COMM_WORLD);
      ibm[ibi].U_lagr_z[l] = u_sum;

//      printf(" the ulagr %le \n", ibm[ibi].U_lagr_z[l]);
    }
  }

  DAVecRestoreArray(fda, Coor, &coor);
  DAVecRestoreArray(fda, user->lUcat, &ucat);
  DAVecRestoreArray(fda, user->lCsi,  &csi);
  DAVecRestoreArray(fda, user->lEta,  &eta);
  DAVecRestoreArray(fda, user->lZet,  &zet);

  return(0);

}


PetscErrorCode Calc_F_eul(UserCtx *user, IBMNodes *ibm, FSInfo *fsi)
{

  DA              da = user->da, fda = user->fda;
  DALocalInfo     info;
  PetscInt        xs, xe, ys, ye, zs, ze; // Local grid information
  PetscInt        mx, my, mz; // Dimensions in three directions
  PetscInt        i, j, k, l, ibi;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  Cmpnts	***lf_eul, ***coor, ***csi, ***eta, ***zet;

  Vec		Coor;

  PetscReal	xc, yc, zc, hx, hy, hz, vol_eul, vol_lagr;
  
  PetscReal	dfunc;

  double r1, r2, r3;

  DAGetLocalInfo(da, &info);
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

//  VecSet(user->lF_eul,0.0);

  DAGetGhostedCoordinates(da, &Coor);
  DAVecGetArray(fda, Coor, &coor);
  DAVecGetArray(fda, user->lF_eul, &lf_eul);
  DAVecGetArray(fda, user->lCsi,  &csi);
  DAVecGetArray(fda, user->lEta,  &eta);
  DAVecGetArray(fda, user->lZet,  &zet);
/*
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
 
        lf_eul[k][j][i].x = 0.0;
        lf_eul[k][j][i].y = 0.0;
        lf_eul[k][j][i].z = 0.0;
      }
    }
  }
*/
  for (ibi=0; ibi<NumberOfTurbines; ibi++) {

    for (k=ibm[ibi].k_min; k<ibm[ibi].k_max; k++) {
      for (j=ibm[ibi].j_min; j<ibm[ibi].j_max; j++) {
        for (i=ibm[ibi].i_min; i<ibm[ibi].i_max; i++) {
          xc = (coor[k][j][i].x + coor[k][j-1][i].x + coor[k-1][j][i].x + coor[k-1][j-1][i].x) * 0.25;
          yc = (coor[k][j][i].y + coor[k][j][i-1].y + coor[k-1][j][i].y + coor[k-1][j][i-1].y) * 0.25;
          zc = (coor[k][j][i].z + coor[k][j-1][i].z + coor[k][j][i-1].z + coor[k][j-1][i-1].z) * 0.25;

          hx = (coor[k][j][i+1].x + coor[k][j-1][i+1].x + coor[k-1][j][i+1].x + coor[k-1][j-1][i+1].x -
                coor[k][j][i-1].x - coor[k][j-1][i-1].x - coor[k-1][j][i-1].x - coor[k-1][j-1][i-1].x ) * 0.125;

          hy = (coor[k][j+1][i].y + coor[k][j+1][i-1].y + coor[k-1][j+1][i].y + coor[k-1][j+1][i-1].y -
                coor[k][j-1][i].y - coor[k][j-1][i-1].y - coor[k-1][j-1][i].y - coor[k-1][j-1][i-1].y ) * 0.125;

          hz = (coor[k+1][j][i].z + coor[k+1][j-1][i].z + coor[k+1][j][i-1].z + coor[k+1][j-1][i-1].z - 
                coor[k-1][j][i].z - coor[k-1][j-1][i].z - coor[k-1][j][i-1].z - coor[k-1][j-1][i-1].z ) * 0.125;

          vol_eul = 1.0/(hx * hy * hz);
          for (l=0; l<ibm[ibi].n_elmt; l++) {
            // vol_lagr = ibm[ibi].dA[l] * (hx + hy + hz) / 3.0;
            r1 = (xc - ibm[ibi].cent_x[l])/hx; r2 = (yc - ibm[ibi].cent_y[l])/hy; r3 = (zc - ibm[ibi].cent_z[l])/hz;
            dfunc = vol_eul * dfunc_s3h(r1) * dfunc_s3h(r2) * dfunc_s3h(r3);

//	    if (k==40 && j==40 && i==40) printf("** why zero %le %le %le \n", ibm[ibi].F_lagr_z[l] , dfunc , ibm[ibi].dA[l]);
            lf_eul[k][j][i].x += ibm[ibi].F_lagr_x[l] * dfunc * ibm[ibi].dA[l] * csi[k][j][i].x +
                                 ibm[ibi].F_lagr_y[l] * dfunc * ibm[ibi].dA[l] * csi[k][j][i].y +
                                 ibm[ibi].F_lagr_z[l] * dfunc * ibm[ibi].dA[l] * csi[k][j][i].z;
          }
        }
      }
    }
  }


  for (ibi=0; ibi<NumberOfTurbines; ibi++) {

    for (k=ibm[ibi].k_min; k<ibm[ibi].k_max; k++) {
      for (j=ibm[ibi].j_min; j<ibm[ibi].j_max; j++) {
        for (i=ibm[ibi].i_min; i<ibm[ibi].i_max; i++) {

          xc = (coor[k][j][i].x + coor[k][j-1][i].x + coor[k-1][j][i].x + coor[k-1][j-1][i].x) * 0.25;
          yc = (coor[k][j][i].y + coor[k][j][i-1].y + coor[k-1][j][i].y + coor[k-1][j][i-1].y) * 0.25;
          zc = (coor[k][j][i].z + coor[k][j-1][i].z + coor[k][j][i-1].z + coor[k][j-1][i-1].z) * 0.25;

          hx = (coor[k][j][i+1].x + coor[k][j-1][i+1].x + coor[k-1][j][i+1].x + coor[k-1][j-1][i+1].x -
                coor[k][j][i-1].x - coor[k][j-1][i-1].x - coor[k-1][j][i-1].x - coor[k-1][j-1][i-1].x ) * 0.125;

          hy = (coor[k][j+1][i].y + coor[k][j+1][i-1].y + coor[k-1][j+1][i].y + coor[k-1][j+1][i-1].y -
                coor[k][j-1][i].y - coor[k][j-1][i-1].y - coor[k-1][j-1][i].y - coor[k-1][j-1][i-1].y ) * 0.125;

          hz = (coor[k+1][j][i].z + coor[k+1][j-1][i].z + coor[k+1][j][i-1].z + coor[k+1][j-1][i-1].z -
                coor[k-1][j][i].z - coor[k-1][j-1][i].z - coor[k-1][j][i-1].z - coor[k-1][j-1][i-1].z ) * 0.125;

          vol_eul = 1.0/(hx * hy * hz);
          for (l=0; l<ibm[ibi].n_elmt; l++) {
            // vol_lagr = ibm[ibi].dA[l] * (hx + hy + hz) / 3.0;
            r1 = (xc - ibm[ibi].cent_x[l])/hx; r2 = (yc - ibm[ibi].cent_y[l])/hy; r3 = (zc - ibm[ibi].cent_z[l])/hz;
            dfunc = vol_eul * dfunc_s3h(r1) * dfunc_s3h(r2) * dfunc_s3h(r3);
            lf_eul[k][j][i].y += ibm[ibi].F_lagr_x[l] * dfunc * ibm[ibi].dA[l] *  eta[k][j][i].x +
                                 ibm[ibi].F_lagr_y[l] * dfunc * ibm[ibi].dA[l] *  eta[k][j][i].y +
                                 ibm[ibi].F_lagr_z[l] * dfunc * ibm[ibi].dA[l] *  eta[k][j][i].z;
          }
        }
      }
    }
  }

  for (ibi=0; ibi<NumberOfTurbines; ibi++) {

    for (k=ibm[ibi].k_min; k<ibm[ibi].k_max; k++) {
      for (j=ibm[ibi].j_min; j<ibm[ibi].j_max; j++) {
        for (i=ibm[ibi].i_min; i<ibm[ibi].i_max; i++) {

          xc = (coor[k][j][i].x + coor[k][j-1][i].x + coor[k-1][j][i].x + coor[k-1][j-1][i].x) * 0.25;
          yc = (coor[k][j][i].y + coor[k][j][i-1].y + coor[k-1][j][i].y + coor[k-1][j][i-1].y) * 0.25;
          zc = (coor[k][j][i].z + coor[k][j-1][i].z + coor[k][j][i-1].z + coor[k][j-1][i-1].z) * 0.25;

          hx = (coor[k][j][i+1].x + coor[k][j-1][i+1].x + coor[k-1][j][i+1].x + coor[k-1][j-1][i+1].x -
                coor[k][j][i-1].x - coor[k][j-1][i-1].x - coor[k-1][j][i-1].x - coor[k-1][j-1][i-1].x ) * 0.125;

          hy = (coor[k][j+1][i].y + coor[k][j+1][i-1].y + coor[k-1][j+1][i].y + coor[k-1][j+1][i-1].y -
                coor[k][j-1][i].y - coor[k][j-1][i-1].y - coor[k-1][j-1][i].y - coor[k-1][j-1][i-1].y ) * 0.125;

          hz = (coor[k+1][j][i].z + coor[k+1][j-1][i].z + coor[k+1][j][i-1].z + coor[k+1][j-1][i-1].z -
                coor[k-1][j][i].z - coor[k-1][j-1][i].z - coor[k-1][j][i-1].z - coor[k-1][j-1][i-1].z ) * 0.125;
  
          vol_eul = 1.0/(hx * hy * hz);
          for (l=0; l<ibm[ibi].n_elmt; l++) {
            // vol_lagr = ibm[ibi].dA[l] * (hx + hy + hz) / 3.0;
            r1 = (xc - ibm[ibi].cent_x[l])/hx; r2 = (yc - ibm[ibi].cent_y[l])/hy; r3 = (zc - ibm[ibi].cent_z[l])/hz;
            dfunc = vol_eul * dfunc_s3h(r1) * dfunc_s3h(r2) * dfunc_s3h(r3);
            lf_eul[k][j][i].z += ibm[ibi].F_lagr_x[l] * dfunc * ibm[ibi].dA[l] *  zet[k][j][i].x +
                                 ibm[ibi].F_lagr_y[l] * dfunc * ibm[ibi].dA[l] *  zet[k][j][i].y +
                                 ibm[ibi].F_lagr_z[l] * dfunc * ibm[ibi].dA[l] *  zet[k][j][i].z;

//		if (l==0) printf("the area %d %le \n", l, ibm[ibi].dA[l]);
          }
        }
      }
    }
  }

/*
  for (ibi=0; ibi<NumberOfTurbines; ibi++) {
    PetscReal f_sum = 0.0;
    PetscReal f_sum0;
    for (k=ibm[ibi].k_min; k<ibm[ibi].k_max; k++) {
      for (j=ibm[ibi].j_min; j<ibm[ibi].j_max; j++) {
        for (i=ibm[ibi].i_min; i<ibm[ibi].i_max; i++) {

          hx = (coor[k][j][i+1].x + coor[k][j-1][i+1].x + coor[k-1][j][i+1].x + coor[k-1][j-1][i+1].x -
                coor[k][j][i-1].x - coor[k][j-1][i-1].x - coor[k-1][j][i-1].x - coor[k-1][j-1][i-1].x ) * 0.125;

          hy = (coor[k][j+1][i].y + coor[k][j+1][i-1].y + coor[k-1][j+1][i].y + coor[k-1][j+1][i-1].y -
              coor[k][j-1][i].y - coor[k][j-1][i-1].y - coor[k-1][j-1][i].y - coor[k-1][j-1][i-1].y ) * 0.125;

          hz = (coor[k+1][j][i].z + coor[k+1][j-1][i].z + coor[k+1][j][i-1].z + coor[k+1][j-1][i-1].z -
              coor[k-1][j][i].z - coor[k-1][j-1][i].z - coor[k-1][j][i-1].z - coor[k-1][j-1][i-1].z ) * 0.125;

          vol_eul = hx * hy * hz;

          f_sum += lf_eul[k][j][i].z * vol_eul / zet[k][j][i].z;
        }
      }
    }
    PetscGlobalSum(&f_sum, &f_sum0, PETSC_COMM_WORLD);
    PetscPrintf(PETSC_COMM_WORLD, "*** Force summation on eul at %i th body: %le \n", ibi, f_sum0);
  }
*/

  
  DAVecRestoreArray(fda, Coor, &coor);
  DAVecRestoreArray(fda, user->lF_eul, &lf_eul);
  DAVecRestoreArray(fda, user->lCsi,  &csi);
  DAVecRestoreArray(fda, user->lEta,  &eta);
  DAVecRestoreArray(fda, user->lZet,  &zet);
  
  DALocalToGlobal(fda, user->lF_eul, INSERT_VALUES, user->F_eul);

  return(0);
}


PetscErrorCode Calc_forces_rotor(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, PetscInt bi)
{
  PetscInt      l, ibi;
  double        pi = 3.141592653589793;
  PetscReal	F_xSum, F_ySum, F_zSum, A_xSum, A_ySum, A_zSum, P_Sum, U_Sum;
  PetscReal F_z, P; 


  for (ibi=0; ibi<NumberOfTurbines; ibi++) {
    F_xSum = 0.0; F_ySum = 0.0; F_zSum = 0.0; A_xSum = 0.0; A_ySum = 0.0; A_zSum = 0.0; P_Sum = 0.0; U_Sum = 0.0;


    for (l=0; l<ibm[ibi].n_elmt; l++) {
      F_xSum += ibm[ibi].F_lagr_x[l] * ibm[ibi].dA[l] ;
      F_ySum += ibm[ibi].F_lagr_y[l] * ibm[ibi].dA[l] ;
      F_zSum += ibm[ibi].F_lagr_z[l] * ibm[ibi].dA[l] ;

      A_xSum += ibm[ibi].dA[l] ;
      A_ySum += ibm[ibi].dA[l] ;
      A_zSum += ibm[ibi].dA[l] ;

      P_Sum += ibm[ibi].F_lagr_z[l] * ibm[ibi].dA[l] * ibm[ibi].U_lagr_z[l];

      U_Sum += ibm[ibi].U_lagr_z[l] * ibm[ibi].dA[l] ;

//      printf("*** check %le %le %le \n", ibm[ibi].F_lagr_z[l] , ibm[ibi].dA[l] , ibm[ibi].nf_z[l]);
    }

    F_z = 4.0 * indf_ax * (1.0 - indf_ax);
    P = F_z * (1.0 - indf_ax);
   
    PetscInt rank=0;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if (!rank) {
      FILE *f;
      char filen[80];
      sprintf(filen, "Force_Coeff_WT_SI%2.2d_%2.2d",ibi,bi);
      f = fopen(filen, "a");
      PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le \n",ti, F_xSum, F_ySum, F_zSum, A_zSum, F_z);
      fclose(f);

      sprintf(filen, "Power_WT_SI%2.2d_%2.2d",ibi,bi);
      f = fopen(filen, "a");
      PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le \n", ti, P_Sum, F_zSum, U_Sum, A_zSum, P);
      fclose(f);

    }

  }
}


PetscErrorCode Calc_forces_ACL(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, PetscInt bi)
{
  PetscInt      l, ibi;
  double        pi = 3.141592653589793;
  PetscReal	F_xSum, F_ySum, F_zSum, A_xSum, A_ySum, A_zSum, P_Sum, U_Sum, M_xSum, M_ySum, M_zSum;
  PetscReal F_z, P, r, rx, ry, rz; 

  PetscInt rank=0;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);


  for (ibi=0; ibi<NumberOfTurbines; ibi++) {
    F_xSum = 0.0; F_ySum = 0.0; F_zSum = 0.0; A_xSum = 0.0; A_ySum = 0.0; A_zSum = 0.0; P_Sum = 0.0; U_Sum = 0.0;
    M_xSum = 0.0; M_ySum = 0.0; M_zSum = 0.0; 

    for (l=0; l<ibm[ibi].n_elmt; l++) {
      F_xSum += ibm[ibi].F_lagr_x[l] * ibm[ibi].dA[l] ;
      F_ySum += ibm[ibi].F_lagr_y[l] * ibm[ibi].dA[l] ;
      F_zSum += ibm[ibi].F_lagr_z[l] * ibm[ibi].dA[l] ;

      A_xSum += ibm[ibi].dA[l] ;
      A_ySum += ibm[ibi].dA[l] ;
      A_zSum += ibm[ibi].dA[l] ;

      P_Sum += ibm[ibi].F_lagr_z[l] * ibm[ibi].dA[l] * ibm[ibi].U_lagr_z[l];

      U_Sum += ibm[ibi].U_lagr_z[l] * ibm[ibi].dA[l] ;


      r = sqrt( pow((ibm[ibi].cent_x[l]-fsi[ibi].x_c),2) + pow((ibm[ibi].cent_y[l]-fsi[ibi].y_c),2) + pow((ibm[ibi].cent_z[l]-fsi[ibi].z_c),2));	
      rx = r * ibm[ibi].nf_x[l]; ry = r * ibm[ibi].nf_y[l]; rz = r * ibm[ibi].nf_z[l];

      M_xSum += ry*ibm[ibi].F_lagr_z[l] * ibm[ibi].dA[l] - rz*ibm[ibi].F_lagr_y[l] * ibm[ibi].dA[l];  
      M_ySum += rz*ibm[ibi].F_lagr_x[l] * ibm[ibi].dA[l] - rx*ibm[ibi].F_lagr_z[l] * ibm[ibi].dA[l];  
      M_zSum += rx*ibm[ibi].F_lagr_y[l] * ibm[ibi].dA[l] - ry*ibm[ibi].F_lagr_x[l] * ibm[ibi].dA[l];  

//      if (!rank) printf("*** check F_zSum  %d %le %le %le \n", l, ibm[ibi].F_lagr_z[l] , ibm[ibi].dA[l] , ibm[ibi].nf_z[l]);
    }

    F_z = 4.0 * indf_ax * (1.0 - indf_ax);
    P = F_z * (1.0 - indf_ax);
   
    if (!rank) {
      FILE *f;
      char filen[80];
      sprintf(filen, "Force_Coeff_WT_SI%2.2d_%2.2d",ibi,bi);
      f = fopen(filen, "a");
      PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le \n",ti, F_xSum, F_ySum, F_zSum, A_zSum, F_z);
      fclose(f);

      sprintf(filen, "Moment_Coeff_WT_SI%2.2d_%2.2d",ibi,bi);
      f = fopen(filen, "a");
      PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le \n",ti, M_xSum, M_ySum, M_zSum, A_zSum);
      fclose(f);


      sprintf(filen, "Power_WT_SI%2.2d_%2.2d",ibi,bi);
      f = fopen(filen, "a");
      PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le \n", ti, M_zSum*angvel, P_Sum);
      fclose(f);

    }

  }
}


PetscErrorCode Uref_ACL(UserCtx *user, IBMNodes *ibm_ACL, IBMNodes *ibm_ACD)
{
  	PetscInt      l, ibi;
  	double        pi = 3.141592653589793;
  	PetscReal     F_zSum, A_zSum, U_zSum, CT1, indf_a;
                PetscInt rank=0;
                MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	Calc_U_lagr(user, ibm_ACD);

  	for (ibi=0; ibi<NumberOfTurbines; ibi++) {

    		F_zSum = 0.0;  
    		for (l=0; l<ibm_ACL[ibi].n_elmt; l++) {
      			F_zSum += ibm_ACL[ibi].F_lagr_z[l] * ibm_ACL[ibi].dA[l] ;
    		}


    		U_zSum = 0.0; A_zSum = 0.0; 
                for (l=0; l<ibm_ACD[ibi].n_elmt; l++) {
                        U_zSum += ibm_ACD[ibi].U_lagr_z[l] * ibm_ACD[ibi].dA[l] ;
                        A_zSum += ibm_ACD[ibi].dA[l] ;

//		if (!rank) printf("CT1 and U_ref %d %le %le %le \n", l, F_zSum, A_zSum, U_zSum); 
                }

                U_zSum /= A_zSum;
		
		CT1 = -2.0 * F_zSum / ( A_zSum * U_zSum * U_zSum );
		

		indf_a = CT1 / ( 4.0 + CT1 );

		
		ibm_ACL[ibi].U_ref = U_zSum; // / (1.0 - indf_a);


		if (!rank) printf("the axial a, CT1 and U_ref %d %le %le %le \n", ibi, indf_a, CT1, ibm_ACL[ibi].U_ref); 

  	}


}










double dfunc_4h(double r)
{
  if (fabs(r) <= 1.0) {
    return (3.0-2.0*fabs(r)+sqrt(1.0+4.0*fabs(r)-4.0*r*r))/8.0;
  } else if (fabs(r) >= 1.0 && fabs(r) <= 2.0) {
    return (5.0-2.0*fabs(r)-sqrt(-7.0+12.0*fabs(r)-4.0*r*r))/8.0;
  } else {
    return 0.0;
  }
}

double dfunc_s3h(double r)
{
  if (fabs(r) <= 1.0) {
    return 17.0/48.0 + sqrt(3.0)*3.14159265/108.0 + fabs(r)/4.0 - r*r/4.0 + (1.0-2.0*fabs(r))*sqrt(-12.0*r*r+12.0*fabs(r)+1.0)/16.0 - sqrt(3.0)*asin(sqrt(3.0)*(2.0*fabs(r)-1.0)/2.0)/12.0;
  } else if (fabs(r) >= 1.0 && fabs(r) <= 2.0) {
    return 55.0/48.0 - sqrt(3.0)*3.14159265/108.0 - 13.0*fabs(r)/12.0 + r*r/4.0 + (2.0*fabs(r)-3.0)*sqrt(-12.0*r*r+36.0*fabs(r)-23.0)/48.0 + sqrt(3.0)*asin(sqrt(3.0)*(2.0*fabs(r)-3.0)/2.0)/36.0;
  } else { 
    return 0.0;
  }

}


double dfunc_2h(double r)
{
  if (fabs(r) < 1.0) {
    return 1.0-fabs(r);
  } else {
    return 0.0;
  }
}

