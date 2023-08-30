#include "variables.h"

PetscErrorCode FSI_DATA_Input(FSInfo *FSinf, int ibi)
{
  int  i;
  PetscReal t;

  FILE *f;
  char filen[80];  
  sprintf(filen, "DATA_FSI%5.5d_%2.2d.dat", (int)ti, ibi);
  f = fopen(filen, "r");
  if (!f) {
    //SETERRQ(1, "Cannot open FSI DATA file");
    PetscPrintf(PETSC_COMM_WORLD, "FSI_data cannot open file !!!!!!!!!!!!\n");
  }
  PetscPrintf(PETSC_COMM_WORLD, "FSI_data input begin %d %s\n",ti,filen);
  //  fscanf(f, "%le %le %le", &(FSinf->red_vel), &(FSinf->damp), &(FSinf->mu_s));	  
  fscanf(f, "%le %le %le", &t, &t, &t);	  
  PetscPrintf(PETSC_COMM_WORLD, "FSI_data input red vel damp mu %le %le %le \n",FSinf->red_vel,FSinf->damp,FSinf->mu_s);
  fscanf(f, "%le %le %le", &(FSinf->x_c), &(FSinf->y_c), &(FSinf->z_c));	  
  fscanf(f, "%le %le %le \n", &(FSinf->F_x),&(FSinf->F_y), &(FSinf->F_z));	 
  fscanf(f, "%le %le %le \n", &(FSinf->M_x),&(FSinf->M_y), &(FSinf->M_z));	  
  //PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le \n", FSinfo->x_c, FSinfo->y_c, FSinfo->z_c);	  

  for (i=0; i<6; i++) {
    fscanf(f, "%le %le %le %le", &(FSinf->S_new[i]),&(FSinf->S_old[i]), &(FSinf->S_real[i]), &(FSinf->S_realm1[i]));
    fscanf(f, "%le %le %le %le", &(FSinf->S_ang_n[i]),&(FSinf->S_ang_o[i]), &(FSinf->S_ang_r[i]), &(FSinf->S_ang_rm1[i]));
  }
  fclose(f);
  PetscPrintf(PETSC_COMM_WORLD, "FSI_data input z, dz/dt  %le %le %le %le\n",FSinf->S_new[4],FSinf->S_new[5],FSinf->red_vel,FSinf->damp);
  PetscPrintf(PETSC_COMM_WORLD, "FSI_data input y, dy/dt  %le %le %le %le\n",FSinf->S_new[2],FSinf->S_new[3],FSinf->red_vel,FSinf->damp);
/*   PetscPrintf(PETSC_COMM_WORLD, "FSI_data input angle_x, dang_x/dt  %le %le %le %le\n",FSinf->S_ang_n[0],FSinf->S_ang_n[1],FSinf->red_vel,FSinf->damp); */
  return(0);
}
/* ==================================================================================             */

/* ==================================================================================             */
PetscErrorCode FSI_DATA_Output(FSInfo *FSinfo, int ibi)
{
  int rank, i;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  PetscBarrier(PETSC_NULL);
  if (!rank) {
    FILE *f;
    char filen[80];
    sprintf(filen, "FSI_position%2.2d",ibi);
    f = fopen(filen, "a");
/*     PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le %le\n",ti, S_new[4],S_new[5],F_z, S_real[4], S_real[5], S_realm1[4], S_realm1[5]); */
/*     PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le %le\n",ti, FSinfo->S_new[2],FSinfo->S_new[3],FSinfo->F_y, FSinfo->S_real[2], FSinfo->S_real[3], FSinfo->S_realm1[2], FSinfo->S_realm1[3]); */
    PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le %le %le %le\n",ti, FSinfo->S_new[2],FSinfo->S_new[3],FSinfo->F_y, FSinfo->S_new[4],FSinfo->S_new[5],FSinfo->F_z,FSinfo->S_new[0],FSinfo->S_new[1],FSinfo->F_x);
    fclose(f);

/*     if (dgf_z) { */
/*     sprintf(filen, "FSI_position_x"); */
/*     f = fopen(filen, "a"); */
/*     PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le %le\n",ti, FSinfo->S_new[4],FSinfo->S_new[5],FSinfo->F_z, FSinfo->S_real[4], FSinfo->S_real[5], FSinfo->S_realm1[4], FSinfo->S_realm1[5]); */
/* /\*     PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le %le\n",ti, FSinfo->S_new[2],FSinfo->S_new[3],FSinfo->F_y, FSinfo->S_real[2], FSinfo->S_real[3], FSinfo->S_realm1[2], FSinfo->S_realm1[3]); *\/ */
/*     fclose(f); */
/*     } */

    sprintf(filen, "FSI_Agnle%2.2d",ibi);
    f = fopen(filen, "a");
    PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le %le\n",ti, FSinfo->S_ang_n[0],FSinfo->S_ang_n[1],FSinfo->M_x,FSinfo->S_ang_r[0],FSinfo->S_ang_r[1],FSinfo->S_ang_rm1[0],FSinfo->S_ang_rm1[1]);
    fclose(f);

    sprintf(filen, "DATA_FSI%5.5d_%2.2d.dat", (int)ti, ibi);
    f = fopen(filen, "w");
    PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le \n", FSinfo->red_vel, FSinfo->damp, FSinfo->mu_s);	  
    PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le \n", FSinfo->x_c, FSinfo->y_c, FSinfo->z_c);	  
    PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le \n", FSinfo->F_x, FSinfo->F_y, FSinfo->F_z);	  
    PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le \n", FSinfo->M_x, FSinfo->M_y, FSinfo->M_z);	  
    for (i=0; i<6; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le %le \n", FSinfo->S_new[i],FSinfo->S_old[i], FSinfo->S_real[i], FSinfo->S_realm1[i]);
      PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le %le \n", FSinfo->S_ang_n[i],FSinfo->S_ang_o[i], FSinfo->S_ang_r[i], FSinfo->S_ang_rm1[i]);
    }
    fclose(f);
  }
  return(0);
}

/* ==================================================================================             */
/*    Old subroutine use FSI_calc_pos instead! */
/* ==================================================================================             */
PetscErrorCode FSI_LinMom(FSInfo *FSinfo,IBMNodes *ibm, PetscReal dt, PetscReal dtime, PetscReal Re,PetscReal ti) 
{ 
  int     i;
  PetscReal    pi=3.141592654;
  PetscReal    S_new[6],S_old[6],S_real[6],S_realm1[6];  
  PetscReal    red_vel, damp, mu_s; // reduced vel, damping coeff, mass coeff
  PetscReal    F_x,F_y,F_z; //Forces andOA Area
  
  // Calc Forces
  //Calc_forces(&FSinfo,&ibm,Re,ti);
  
  // init values
  for (i=0;i<6;i++) {
    S_new[i]=FSinfo->S_new[i];
    S_real[i]=FSinfo->S_real[i];
    S_realm1[i]=FSinfo->S_realm1[i];
    FSinfo->S_old[i]=FSinfo->S_new[i];
    S_old[i]=S_new[i];
  }
  
  red_vel=FSinfo->red_vel;
  damp   =FSinfo->damp   ;
  mu_s   =FSinfo->mu_s   ;

  F_x = FSinfo->F_x;
  F_y = FSinfo->F_y;
  F_z = FSinfo->F_z;

  // solve lin mom equ
  S_new[0]=S_new[0]; // x
  S_new[1]=S_new[1]; // dx/dt
  S_new[2]=S_new[2]-dt/2./dtime*(3*S_new[3]-4*S_real[2]+S_realm1[2])+S_new[3]*dt; // y
  // dy/dt
  S_new[3]=S_new[3]-dt/2./dtime*(3*S_new[3]-4*S_real[3]+S_realm1[3])
    +dt*(-2.*damp*(red_vel)*S_new[3] 
	 -(red_vel*red_vel)*S_old[2] 
	 + mu_s*F_y);

  S_new[4]=S_new[4]; //z
  S_new[5]=S_new[5]; //dz/dt

  // store results
  for (i=0;i<6;i++){
    FSinfo->S_new[i]=S_new[i];
  }

  return(0);
}
/* ==================================================================================             */

/* ==================================================================================             */
/*    Old subroutine! use FSI_calc_ang instead! */
/* ==================================================================================             */
PetscErrorCode FSI_AngMom(FSInfo *FSinfo,IBMNodes *ibm, PetscReal dt, PetscReal dtime, PetscReal Re) 
{ 
  
  int     i;
  PetscReal    pi=3.141592654;
  PetscReal    S_ang_n[6],S_ang_o[6],S_ang_r[6],S_ang_rm1[6];  
  PetscReal    red_vel, damp, mu_s; // reduced vel, damping coeff, mass coeff
  PetscReal    M_x,M_y,M_z; //Forces and Area
  
  // Calc Forces
  //Calc_Moments(&FSinfo,&ibm,Re);
  
  // init values
  for (i=0;i<6;i++) {
    S_ang_n[i]=FSinfo->S_ang_n[i];
    S_ang_r[i]=FSinfo->S_ang_r[i];
    S_ang_rm1[i]=FSinfo->S_ang_rm1[i];
    FSinfo->S_ang_o[i]=FSinfo->S_ang_n[i];
    S_ang_o[i]=S_ang_n[i];
  }
  
  red_vel=FSinfo->red_vel;
  damp   =FSinfo->damp   ;
  mu_s   =FSinfo->mu_s   ;

  M_x = FSinfo->M_x;
  M_y = FSinfo->M_y;
  M_z = FSinfo->M_z;

  // solve Ang mom equ
  S_ang_n[0]=S_ang_n[0]-dt/2./dtime*(3*S_ang_n[0]-4*S_ang_r[0]+S_ang_rm1[0])+S_ang_n[1]*dt; // x
  S_ang_n[1]=S_ang_n[1]-dt/2./dtime*(3*S_ang_n[1]-4*S_ang_r[1]+S_ang_rm1[1])
    +dt*(-2.*damp*(red_vel)*S_ang_n[1] 
	 -(red_vel*red_vel)*S_ang_o[0] 
	 + mu_s*M_x); // dx/dt
  S_ang_n[2]=S_ang_n[2];//-dt/2./dtime*(3*S_ang_n[3]-4*S_ang_r[2]+S_ang_rm1[2])+S_ang_n[3]*dt; // y
  // dy/dt
  S_ang_n[3]=S_ang_n[3];/*-dt/2./dtime*(3*S_ang_n[3]-4*S_ang_r[3]+S_ang_rm1[3])
    +dt*(-2.*damp*(red_vel)*S_ang_n[3] 
	 -(red_vel*red_vel)*S_ang_o[2] 
	 + mu_s*M_y);*/

  S_ang_n[4]=S_ang_n[4]; //z
  S_ang_n[5]=S_ang_n[5]; //dz/dt

  // store results
  for (i=0;i<6;i++){
    FSinfo->S_ang_n[i]=S_ang_n[i];
  }
  
  return(0);
}

/* ==================================================================================             */
PetscErrorCode Calc_FSI_pos(FSInfo *FSinfo,IBMNodes *ibm, 
			    PetscReal dt, PetscReal dtime, 
			    PetscReal Re) 
{ 
  int     i,j;
  int     itr=13;
  PetscReal    pi=3.141592654;
  PetscReal    S_new[6],S_old[6],S_real[6],S_realm1[6];  
  PetscReal    red_vel, damp, mu_s; // reduced vel, damping coeff, mass coeff
  PetscReal    F_x,F_y,F_z; //Forces and Area
  
  // Calc Forces
  //PetscPrintf(PETSC_COMM_WORLD, "n_elmt in calc_pos  %d\n", ibm->n_elmt);
  //Calc_forces(&FSinfo,&ibm,Re,ti);
  
  // init values
  for (i=0;i<6;i++) {
    S_new[i]=FSinfo->S_new[i];
    S_real[i]=FSinfo->S_real[i];
    S_realm1[i]=FSinfo->S_realm1[i];
    FSinfo->S_old[i]=FSinfo->S_new[i];
    //S_old[i]=S_new[i];
  }
  
  red_vel=FSinfo->red_vel;
  damp   =FSinfo->damp   ;
  mu_s   =FSinfo->mu_s   ;

  F_x = 0.5*(FSinfo->F_x + FSinfo->F_x_old);
  F_y = 0.5*(FSinfo->F_y + FSinfo->F_y_old);
  F_z = 0.5*(FSinfo->F_z + FSinfo->F_z_old);

  // solve lin mom equ
  for (i=0; i< itr;i++) { 

    for (j=0;j<6;j++) {
      S_old[j]=S_new[j];
    }

    S_new[0]=S_new[0]; // x
    S_new[1]=S_new[1]; // dx/dt
    S_new[2]=S_new[2]-dt/2./dtime*
      (3.*S_new[2]-4.*S_real[2]+S_realm1[2])+S_new[3]*dt; // y
    // dy/dt
    S_new[3]=S_new[3]-dt/2./dtime*
         (3.*S_new[3]-4.*S_real[3]+S_realm1[3])+
      dt*(-2.*damp*(red_vel)*S_new[3]
	  -(red_vel*red_vel)*S_old[2]
	  + mu_s*F_y);
    
    S_new[4]=S_new[4];//-dt/2./dtime*(3*S_new[4]-4*S_real[4]+S_realm1[4])+S_new[5]*dt; //z
    S_new[5]=S_new[5];//-dt/2./dtime*(3*S_new[5]-4*S_real[5]+S_realm1[5])
    /*     +dt*(-2.*damp*(red_vel)*S_new[5] */
    /* 	 -(red_vel*red_vel)*(S_old[4]) */
    /* 	 + mu_s*F_z); //dz/dt */
  }
  
  // FSI convergence
  //PetscPrintf(PETSC_COMM_SELF, "FSI convergence z: %le  u_z:%le\n", S_new[4]-S_old[4],S_new[5]-S_old[5]);
  PetscPrintf(PETSC_COMM_WORLD, "FSI convergence z: %le  u_z:%le\n", S_new[2]-S_old[2],S_new[3]-S_old[3]);
  
  // store results
  for (i=0;i<6;i++){
    FSinfo->S_realm1[i]=FSinfo->S_real[i];
    FSinfo->S_real[i]=S_new[i];
    FSinfo->S_new[i]=S_new[i];
  }

  // output values
  PetscPrintf(PETSC_COMM_WORLD, "z, dz/dt  %le %le\n",S_new[4],S_new[5]);
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  PetscBarrier(PETSC_NULL);
  if (!rank) {
    FILE *f;
    char filen[80];
    sprintf(filen, "FSI_position");
    f = fopen(filen, "a");
/*     PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le %le\n",ti, S_new[4],S_new[5],F_z, S_real[4], S_real[5], S_realm1[4], S_realm1[5]); */
    PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le %le\n",ti, S_new[2],S_new[3],F_z, S_real[2], S_real[3], S_realm1[2], S_realm1[3]);
    fclose(f);

    sprintf(filen, "DATA_FSI%5.5d.dat", (int)ti);
    f = fopen(filen, "w");
    PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le \n", FSinfo->red_vel, FSinfo->damp, FSinfo->mu_s);	  
    PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le \n", FSinfo->x_c, FSinfo->y_c, FSinfo->z_c);	  
    PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le \n", FSinfo->F_x, FSinfo->F_y, FSinfo->F_z);	  
    PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le \n", FSinfo->M_x, FSinfo->M_y, FSinfo->M_z);	  
    for (i=0; i<6; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le %le \n", FSinfo->S_new[i],FSinfo->S_old[i], FSinfo->S_real[i], FSinfo->S_realm1[i]);
      PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le %le \n", FSinfo->S_ang_n[i],FSinfo->S_ang_o[i], FSinfo->S_ang_r[i], FSinfo->S_ang_rm1[i]);
    }
    fclose(f);
  }
  return(0);
}

/* ==================================================================================             */
/* For strong-coupling */
/* ==================================================================================             */
PetscErrorCode Calc_FSI_pos_SC(FSInfo *FSinfo,IBMNodes *ibm, 
			       PetscReal dt, PetscReal dtime, 
			       PetscReal Re) 
{ 
  int     i,j;
  int     itr=23;
  PetscReal    pi=3.141592654;
  PetscReal    S_new[6],S_old[6],S_real[6],S_realm1[6];  
  PetscReal    red_vel, damp, mu_s; // reduced vel, damping coeff, mass coeff
  PetscReal    F_x,F_y,F_z; //Forces and Area
  
  // Calc Forces
  //PetscPrintf(PETSC_COMM_WORLD, "n_elmt in calc_pos  %d\n", ibm->n_elmt);
  //Calc_forces(&FSinfo,&ibm,Re,ti);
  
  // init values
  for (i=0;i<6;i++) {
    S_new[i]=FSinfo->S_real[i];//FSinfo->S_new[i];
    S_real[i]=FSinfo->S_real[i];
    S_realm1[i]=FSinfo->S_realm1[i];
    //FSinfo->S_old[i]=FSinfo->S_new[i];
    //S_old[i]=S_new[i];
  }
  
  red_vel=FSinfo->red_vel;
  damp   =FSinfo->damp   ;
  mu_s   =FSinfo->mu_s   ;

  F_x = FSinfo->F_x;
  F_y = FSinfo->F_y;
  F_z = FSinfo->F_z;
/*   F_x = 0.5*(FSinfo->F_x + FSinfo->F_x_real); */
/*   F_y = 0.5*(FSinfo->F_y + FSinfo->F_y_real); */
/*   F_z = 0.5*(FSinfo->F_z + FSinfo->F_z_real); */

  PetscPrintf(PETSC_COMM_WORLD, "FSI  %le %le %le %le %le %le %le %le %le\n",red_vel,damp,mu_s,dt,dtime,F_y,F_z, S_real[2],S_realm1[2] );

  // solve lin mom equ
  for (i=0; i< itr;i++) { 

    for (j=0;j<6;j++) {
      S_old[j]=S_new[j];
    }

    if (dgf_x) {
    S_new[0]=S_new[0]-dt/2./dtime*
      (3.*S_new[0]-4.*S_real[0]+S_realm1[0])+S_new[1]*dt; // x
    S_new[1]=S_new[1]-dt/2./dtime*
         (3.*S_new[1]-4.*S_real[1]+S_realm1[1])+
      dt*(-2.*damp*(red_vel)*S_new[1]
	  -(red_vel*red_vel)*S_old[0]
	  + mu_s*F_x); // dx/dt
    }
    if (dgf_y) {
    S_new[2]=S_new[2]-dt/2./dtime*
      (3.*S_new[2]-4.*S_real[2]+S_realm1[2])+S_new[3]*dt; // y    
    // dy/dt
    S_new[3]=S_new[3]-dt/2./dtime*
         (3.*S_new[3]-4.*S_real[3]+S_realm1[3])+
      dt*(-2.*damp*(red_vel)*S_new[3]
	  -(red_vel*red_vel)*S_old[2]
	  + mu_s*F_y);
    }
    if (dgf_z) {
    S_new[4]=S_new[4]-dt/2./dtime*
      (3*S_new[4]-4*S_real[4]+S_realm1[4])+S_new[5]*dt; //z
    S_new[5]=S_new[5]-dt/2./dtime*(3*S_new[5]-4*S_real[5]+S_realm1[5])
              +dt*(-2.*damp*(red_vel)*S_new[5]
		   -(red_vel*red_vel)*(S_old[4])
		   + mu_s*F_z); //dz/dt
    }

    // FSI convergence
    //PetscPrintf(PETSC_COMM_SELF, "FSI convergence z: %le  u_z:%le\n", S_new[4]-S_old[4],S_new[5]-S_old[5]);
    PetscPrintf(PETSC_COMM_WORLD, "FSI convergence y: %le  u_y:%le\n", S_new[2]-S_old[2],S_new[3]-S_old[3]);

  }
    
  // store results
  for (i=0;i<6;i++){
    FSinfo->S_new[i]=S_new[i];
  }

  // output values
  PetscPrintf(PETSC_COMM_WORLD, "z, dz/dt %le %le %le\n",S_new[4],S_new[5], F_z);
  PetscPrintf(PETSC_COMM_WORLD, "y, dy/dt %le %le %le\n",S_new[2],S_new[3], F_y);
  PetscPrintf(PETSC_COMM_WORLD, "x, dx/dt %le %le %le\n",S_new[0],S_new[1], F_x);

  return(0);
}

/* ==================================================================================             */
/*  integral equation For strong/weak-coupling */
/* ==================================================================================             */
PetscErrorCode Calc_FSI_pos_intg(FSInfo *FSinfo,IBMNodes *ibm, 
			       PetscReal dt,
			       PetscReal Re) 
{ 
  int     i,j;
  int     itr=23;
  PetscReal    pi=3.141592654;
  PetscReal    S_new[6],S_old[6],S_real[6],S_realm1[6];  
  PetscReal    red_vel, damp, mu_s; // reduced vel, damping coeff, mass coeff
  PetscReal    F_x,F_y,F_z; //Forces and Area
  PetscReal    w=0.75;

  // Calc Forces
  //PetscPrintf(PETSC_COMM_WORLD, "n_elmt in calc_pos  %d\n", ibm->n_elmt);
  //Calc_forces(&FSinfo,&ibm,Re,ti);
  
  // init values
  for (i=0;i<6;i++) {
    //S_new[i]=FSinfo->S_real[i];//FSinfo->S_new[i];
    S_real[i]=FSinfo->S_real[i];
    S_realm1[i]=FSinfo->S_realm1[i];
    S_old[i]=FSinfo->S_old[i];    
  }
  
  red_vel=FSinfo->red_vel;
  damp   =FSinfo->damp   ;
  mu_s   =FSinfo->mu_s   ;

  F_x = 0.5*(FSinfo->F_x + FSinfo->F_x_real);
  F_y = 0.5*(FSinfo->F_y + FSinfo->F_y_real);
  F_z = 0.5*(FSinfo->F_z + FSinfo->F_z_real);

  PetscPrintf(PETSC_COMM_WORLD, "FSI  %le %le %le %le %le %le %le %le\n",red_vel,damp,mu_s,dt,F_y,F_z, S_real[2],S_realm1[2] );

  // solve lin mom equ
  if (dgf_x) {
    S_new[1]=S_real[1]+ dt*(mu_s*F_x); // u=u_r + int(F/mdt)

    S_new[0]=S_real[0]+0.5*(S_new[1]+S_real[1])*dt; //x=x_r+u_avedt
  }
  if (dgf_y) {
    S_new[3]=S_real[3]+ dt*(mu_s*F_y); // u=u_r + int(F/mdt)

    S_new[2]=S_real[2]+0.5*(S_new[3]+S_real[3])*dt; //y=y_r+u_avedt
  }
  if (dgf_z) {
    S_new[5]=S_real[5]+ dt*(mu_s*F_z); // u=u_r + int(F/mdt)

    S_new[4]=S_real[4]+0.5*(S_new[5]+S_real[5])*dt; //x=x_r+u_avedt
  }
    
  
  // store results
  for (i=0;i<6;i++){
    // Relaxation 
    S_new[i]=w*S_new[i]+(1-w)*S_old[i];
    FSinfo->S_new[i]=S_new[i];
  }

  // output values
  PetscPrintf(PETSC_COMM_WORLD, "z, dz/dt %le %le %le\n",S_new[4],S_new[5], F_z);
  PetscPrintf(PETSC_COMM_WORLD, "y, dy/dt %le %le %le\n",S_new[2],S_new[3], F_y);
  PetscPrintf(PETSC_COMM_WORLD, "x, dx/dt %le %le %le\n",S_new[0],S_new[1], F_x);

  return(0);
}
/* ==================================================================================             */
/* ==================================================================================             */

PetscErrorCode Forced_Motion(FSInfo *fsi,PetscReal A, 
			     PetscReal dt)
{
  //  PetscReal  pi = 3.141592653589793;
  PetscReal t,w;
  t = (ti)*dt;
  w=fsi->red_vel;
  fsi->S_new[2]= A*sin(w*t);
  fsi->S_new[3]= A*w*cos(w*t);
  PetscPrintf(PETSC_COMM_WORLD, "y, dy/dt %le %le %le\n",t,fsi->S_new[2],fsi->S_new[3]);

  return(0);
}


/* ==================================================================================             */
PetscErrorCode Calc_FSI_Ang(FSInfo *FSinfo,IBMNodes *ibm, 
			    PetscReal dt, PetscReal dtime,
			    int ibi, UserCtx *user)
{  
  int     i,itr=12,j,nv;
  PetscReal    pi=3.141592654;
  PetscReal    S_ang_n[6],S_ang_o[6],S_ang_r[6],S_ang_rm1[6];  
  PetscReal    red_vel, damp, mu_s; // reduced vel, damping coeff, mass coeff
  PetscReal    M_x,M_y,M_z; //Forces and Area
  PetscReal    rx,ry,rz;
  PetscReal    x_c=FSinfo->x_c, y_c=FSinfo->y_c, z_c=FSinfo->z_c;
  PetscReal    wx=0., wy=0., wz=0.;
  PetscReal    w=0.82;

  // Calc Forces
  //Calc_Moments(&FSinfo,&ibm,Re);
  
  // init values
  for (i=0;i<6;i++) {
    S_ang_n[i]=FSinfo->S_ang_r[i];
    S_ang_r[i]=FSinfo->S_ang_r[i];
    S_ang_rm1[i]=FSinfo->S_ang_rm1[i];
    //FSinfo->S_ang_o[i]=FSinfo->S_ang_n[i];
    //S_ang_o[i]=S_ang_n[i];
  }
  
  red_vel=FSinfo->red_vel;
  damp   =FSinfo->damp   ;
  mu_s   =FSinfo->mu_s   ;

  M_x = FSinfo->M_x;
  M_y = FSinfo->M_y;
  M_z = FSinfo->M_z;

  // solve Ang mom equw
  for (i=0; i< itr;i++) { 

  for (j=0;j<6;j++) {
    S_ang_o[j]=S_ang_n[j];
  }

  if (dgf_ax) {
  S_ang_n[0]=S_ang_n[0]-dt/2./dtime*
    (3*S_ang_n[0]-4*S_ang_r[0]+S_ang_rm1[0])+S_ang_n[1]*dt; // x
  S_ang_n[1]=S_ang_n[1]-dt/2./dtime*(3*S_ang_n[1]-4*S_ang_r[1]+S_ang_rm1[1])
    +dt*(-2.*damp*(red_vel)*S_ang_n[1] 
	 //-(red_vel*red_vel)*S_ang_o[0] 
	 + mu_s*M_x); // dx/dt
  }
  if (dgf_ay) {
  S_ang_n[2]=S_ang_n[2]-dt/2./dtime*(3*S_ang_n[2]-4*S_ang_r[2]+S_ang_rm1[2])+S_ang_n[3]*dt; // y
  // dy/dt
  S_ang_n[3]=S_ang_n[3]-dt/2./dtime*(3*S_ang_n[3]-4*S_ang_r[3]+S_ang_rm1[3])
    +dt*(-2.*damp*(red_vel)*S_ang_n[3] 
	 -(red_vel*red_vel)*S_ang_o[2] 
	 + mu_s*M_y);
  }
  if (dgf_az) {
  S_ang_n[4]=S_ang_n[4]-dt/2./dtime*(3*S_ang_n[4]-4*S_ang_r[4]+S_ang_rm1[4])+S_ang_n[5]*dt; //z
  S_ang_n[5]=S_ang_n[5]-dt/2./dtime*(3*S_ang_n[5]-4*S_ang_r[5]+S_ang_rm1[5])
    +dt*(-2.*damp*(red_vel)*S_ang_n[5] 
	 -(red_vel*red_vel)*S_ang_o[4] 
	 + mu_s*M_z); //dz/dt
  }

/* ==================================================================================             */
/*    making Moment from dUn/dt implicit !!! */
/* ==================================================================================             */

/*   // 1) Calc new vel of all surf nodes */
/*   wx= S_ang_n[1]; */
/*   for (nv=0; nv<ibm->n_v; i++) { */
/*     rx = ibm->x_bp[nv]-x_c; */
/*     ry = ibm->y_bp[nv]-y_c; */
/*     rz = ibm->z_bp[nv]-z_c;       */
/*     ibm->u[nv].x =   ry*wz-wy*rz  ; */
/*     ibm->u[nv].y =-( rx*wz-wx*rz ); */
/*     ibm->u[nv].z =   rx*wy-wx*ry  ;       */
/*   } */

/*   // 2) ibm interpolation to calc new P+dU/dt on all ibm nodes */
/*   PetscPrintf(PETSC_COMM_SELF, "IBM_INP\n"); */
/*   ibm_interpolation_advanced(user, ibm, FSinfo, ibi); */


/*   // 3) calc new moments */
/*   PetscPrintf(PETSC_COMM_SELF, "Calc_f\n"); */
/*   Calc_forces_SI(FSinfo,user,ibm, ti, ibi, 0); */

/*   M_z = FSinfo->M_z; */

/* ==================================================================================             */
/* ==================================================================================             */

  }

  // FSI convergence
  PetscPrintf(PETSC_COMM_WORLD, "FSI convergence ang_x: %le  w_x:%le\n", S_ang_n[0]-S_ang_o[0],S_ang_n[1]-S_ang_o[1]);

  // store results
  for (i=0;i<6;i++){
/*     FSinfo->S_ang_rm1[i]=FSinfo->S_ang_r[i]; */
/*     FSinfo->S_ang_r[i]=S_ang_n[i]; */

    // Relaxation
    //    S_ang_n[i]=w*S_ang_n[i]+(1-w)*FSinfo->S_ang_o[i];

    FSinfo->S_ang_n[i]=S_ang_n[i];
  }
  
  // output values
  PetscPrintf(PETSC_COMM_WORLD, "Ang_x, dAng_x/dt M_x  %le %le %le %le %le\n",S_ang_n[0],S_ang_n[1], M_x, S_ang_r[0],S_ang_rm1[0]);
/*   PetscPrintf(PETSC_COMM_WORLD, "Ang_y, dAng_y/dt M_y  %le %le %le %le %le\n",S_ang_n[2],S_ang_n[3], M_y, S_ang_r[2],S_ang_rm1[2]); */
/*   PetscPrintf(PETSC_COMM_WORLD, "Ang_z, dAng_z/dt M_z  %le %le %le %le %le\n",S_ang_n[4],S_ang_n[5], M_z, S_ang_r[4],S_ang_rm1[4]); */
/*   int rank; */
/*   MPI_Comm_rank(PETSC_COMM_WORLD, &rank); */
/*   PetscBarrier(PETSC_NULL); */
/*   if (!rank) { */
/*     FILE *f; */
/*     char filen[80]; */
/*     sprintf(filen, "FSI_Agnle"); */
/*     f = fopen(filen, "a"); */
/*     PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le %le\n",ti, S_ang_n[0],S_ang_n[1],M_x,S_ang_r[0],S_ang_r[1],S_ang_rm1[0],S_ang_rm1[1]); */
/*     fclose(f); */

/*     sprintf(filen, "DATA_FSI%5.5d.dat",ti); */
/*     f = fopen(filen, "w"); */
/*     PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le \n", FSinfo->red_vel, FSinfo->damp, FSinfo->mu_s);	   */
/*     PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le \n", FSinfo->x_c, FSinfo->y_c, FSinfo->z_c);	   */
/*     PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le \n", FSinfo->F_x, FSinfo->F_y, FSinfo->F_z);	   */
/*     PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le \n", FSinfo->M_x, FSinfo->M_y, FSinfo->M_z);	   */
/*     for (i=0; i<6; i++) { */
/*       PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le %le \n", FSinfo->S_new[i],FSinfo->S_old[i], FSinfo->S_real[i], FSinfo->S_realm1[i]); */
/*       PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le %le \n", FSinfo->S_ang_n[i],FSinfo->S_ang_o[i], FSinfo->S_ang_r[i], FSinfo->S_ang_rm1[i]); */
/*     } */
/*     fclose(f); */
/*   } */

  return(0);
}

/* ==================================================================================             */
PetscErrorCode Calc_FSI_Ang_intg(FSInfo *FSinfo,IBMNodes *ibm, 
				 PetscReal dt,int itrSC,
				 int ibi, UserCtx *user)
{  
  int     i,itr=12,j,nv;
  PetscReal    pi=3.141592654;
  PetscReal    S_ang_n[6],S_ang_r[6],S_ang_rm1[6],S_ang_o[6];  
  PetscReal    red_vel, damp, mu_s; // reduced vel, damping coeff, mass coeff
  PetscReal    M_x,M_y,M_z; //Forces and Area
  PetscReal    rx,ry,rz;
  PetscReal    x_c=FSinfo->x_c, y_c=FSinfo->y_c, z_c=FSinfo->z_c;
  PetscReal    wx=0., wy=0., wz=0.;
  PetscReal    w=.5, wf=1.;//0.82;
  PetscReal    Mdpdn_x;

  // Calc Forces
  //Calc_Moments(&FSinfo,&ibm,Re);
  
  // init values
  for (i=0;i<6;i++) {
    S_ang_o[i]=FSinfo->S_ang_o[i];
    S_ang_r[i]=FSinfo->S_ang_r[i];
    S_ang_rm1[i]=FSinfo->S_ang_rm1[i];
  }
  
  red_vel=FSinfo->red_vel;
  damp   =FSinfo->damp   ;
  mu_s   =FSinfo->mu_s   ;
  Mdpdn_x=FSinfo->Mdpdn_x;

  // Trapezoidal rule
/*   if (STRONG_COUPLING) */
/*   M_x = 0.5*(wf*FSinfo->M_x+(1-wf)*FSinfo->M_x_old */
/* 	     + FSinfo->M_x_real); */
/*     M_x = FSinfo->M_x; // 1st order */
/*   else */
/*     M_x = (55.*FSinfo->M_x - 59.*FSinfo->M_x_real + 37.*FSinfo->M_x_rm2 -9.*FSinfo->M_x_rm3)/24.; */

  M_x = 0.5*(FSinfo->M_x + FSinfo->M_x_real);
  M_y = 0.5*(FSinfo->M_y + FSinfo->M_y_real);
  M_z = 0.5*(FSinfo->M_z + FSinfo->M_z_real);

/*   M_x = (2.*FSinfo->M_x - FSinfo->M_x_real); */
/*   M_y = 0.5*(FSinfo->M_y + FSinfo->M_y_real); */
/*   M_z = 0.5*(FSinfo->M_z + FSinfo->M_z_real); */

/*   if (itrSC==1) { */
/*     M_x = FSinfo->M_x; */
/*     M_y = FSinfo->M_y; */
/*     M_z = FSinfo->M_z; */
/*   } else { */
/* /\*     M_x = 0.5*(FSinfo->M_x + FSinfo->M_x_old); *\/ */
/* /\*     M_y = 0.5*(FSinfo->M_y + FSinfo->M_y_old); *\/ */
/* /\*     M_z = 0.5*(FSinfo->M_z + FSinfo->M_z_old); *\/ */
/*   M_x = 0.5*(FSinfo->M_x + FSinfo->M_x_real); */
/*   M_y = 0.5*(FSinfo->M_y + FSinfo->M_y_real); */
/*   M_z = 0.5*(FSinfo->M_z + FSinfo->M_z_real); */
/*   } */

  // solve Ang mom equ
  if (dgf_ax) {
    S_ang_n[1]=(1-damp*dt)/(1.+damp*dt)*S_ang_r[1]
                      + dt/(1.+damp*dt)*(mu_s*M_x);///(1.+mu_s*Mdpdn_x); // w=w_r + int(M/Idt)
    S_ang_n[0]=S_ang_r[0]+0.5*(S_ang_n[1]+S_ang_r[1])*dt; //ang=ang_r+w_avedt
/*     S_ang_n[1]=S_ang_rm1[1]-damp*2.*dt*S_ang_r[1] */
/*                       + 2*dt*(mu_s*M_x);///(1.+mu_s*Mdpdn_x); // w=w_r + int(M/Idt) */
/*     S_ang_n[0]=S_ang_rm1[0]+2*S_ang_r[1]*dt; //ang=ang_r+w_avedt */
  }
  if (dgf_ay) {
    S_ang_n[3]=S_ang_r[3]+ dt*(mu_s*M_y); // w=w_r + int(M/Idt)

    S_ang_n[2]=S_ang_r[2]+0.5*(S_ang_n[3]+S_ang_r[3])*dt; //ang=ang_r+w_avedt
  }
  if (dgf_az) {
    S_ang_n[5]=S_ang_r[5]+ dt*(mu_s*M_z); // w=w_r + int(M/Idt)

    S_ang_n[4]=S_ang_r[4]+0.5*(S_ang_n[5]+S_ang_r[5])*dt; //ang=ang_r+w_avedt
  }

  
    // Relaxation
  if (STRONG_COUPLING) {
    FSinfo->atk=0.;
    for (i=1;i<6;i+=2){
      FSinfo->dS[i]=S_ang_o[i]-S_ang_n[i];
    
      if (fabs(FSinfo->dS[i]-FSinfo->dS_o[i])>1e-8  &&
     	  FSinfo->atk_o!=0.3) {
	FSinfo->atk+=(FSinfo->dS[i])/
	  (FSinfo->dS_o[i]-FSinfo->dS[i]);
      }
    }
    FSinfo->atk=FSinfo->atk_o+(FSinfo->atk_o-1)*FSinfo->atk;
    if (FSinfo->atk>.9) FSinfo->atk=.9;
    if (FSinfo->atk<-.2) FSinfo->atk=-0.2;
    
    w=1.-FSinfo->atk;
    for (i=1;i<6;i+=2){
      S_ang_n[i]=w*S_ang_n[i]+(1.-w)*S_ang_o[i];
      S_ang_n[i-1]=S_ang_r[i-1]+0.5*(S_ang_n[i]+S_ang_r[i])*dt;
    }
  }

  // store results
  for (i=0;i<6;i++){
    FSinfo->S_ang_n[i]=S_ang_n[i];
  }
  
  // output values
  //PetscPrintf(PETSC_COMM_WORLD, "Ang_x, dAng_x/dt M_x w %le %le %le %le %le %le\n",S_ang_n[0],S_ang_n[1], M_x, S_ang_r[0],S_ang_rm1[0],w);
  PetscPrintf(PETSC_COMM_WORLD, "Ang_x, dAng_x/dt M_x  %le %le %le %le %le %le %le\n",S_ang_n[0],S_ang_n[1], M_x, S_ang_r[0],S_ang_\
rm1[0],w,FSinfo->dS[1]);
  return(0);
}

/* ==================================================================================             */
PetscErrorCode Calc_FSI_Ang_staggered(FSInfo *FSinfo,IBMNodes *ibm, 
				 PetscReal dt, 
				 int ibi, UserCtx *user)
{  
  int     i,itr=12,j,nv;
  PetscReal    pi=3.141592654;
  PetscReal    S_ang_n[6],S_ang_r[6],S_ang_rm1[6],S_ang_o[6];  
  PetscReal    red_vel, damp, mu_s; // reduced vel, damping coeff, mass coeff
  PetscReal    M_x,M_y,M_z; //Forces and Area
  PetscReal    rx,ry,rz;
  PetscReal    x_c=FSinfo->x_c, y_c=FSinfo->y_c, z_c=FSinfo->z_c;
  PetscReal    wx=0., wy=0., wz=0.;
  PetscReal    w=0.72, wf=0.82;
  PetscReal    Mdpdn_x;

  // Calc Forces
  //Calc_Moments(&FSinfo,&ibm,Re);
  
  // init values
  for (i=0;i<6;i++) {
    S_ang_o[i]=FSinfo->S_ang_o[i];
    S_ang_r[i]=FSinfo->S_ang_r[i];
    S_ang_rm1[i]=FSinfo->S_ang_rm1[i];
  }
  
  red_vel=FSinfo->red_vel;
  damp   =FSinfo->damp   ;
  mu_s   =FSinfo->mu_s   ;
  Mdpdn_x=FSinfo->Mdpdn_x;

  M_x = FSinfo->M_x;
  M_y = FSinfo->M_y;
  M_z = FSinfo->M_z;

  // solve Ang mom equ
  if (dgf_ax) {
    S_ang_n[1]=(1-damp*dt)/(1.+damp*dt)*S_ang_r[1]
                      + dt/(1.+damp*dt)*(mu_s*M_x);///(1.+mu_s*Mdpdn_x); // w=w_r + int(M/Idt)

    S_ang_n[0]=S_ang_r[0]+0.5*(S_ang_n[1]+S_ang_r[1])*dt; //ang=ang_r+w_avedt
  }
  if (dgf_ay) {
    S_ang_n[3]=S_ang_r[3]+ dt*(mu_s*M_y); // w=w_r + int(M/Idt)

    S_ang_n[2]=S_ang_r[2]+0.5*(S_ang_n[3]+S_ang_r[3])*dt; //ang=ang_r+w_avedt
  }
  if (dgf_az) {
    S_ang_n[5]=S_ang_r[5]+ dt*(mu_s*M_z); // w=w_r + int(M/Idt)

    S_ang_n[4]=S_ang_r[4]+0.5*(S_ang_n[5]+S_ang_r[5])*dt; //ang=ang_r+w_avedt
  }

  
  // store results
  for (i=0;i<6;i++){
    FSinfo->S_ang_n[i]=S_ang_n[i];
  }
  
  // output values
  PetscPrintf(PETSC_COMM_WORLD, "Ang_x, dAng_x/dt M_x  %le %le %le %le %le\n",S_ang_n[0],S_ang_n[1], M_x, S_ang_r[0],S_ang_rm1[0]);

  return(0);
}

/* ==================================================================================             */

/* ==================================================================================             */
PetscErrorCode Forced_Rotation(FSInfo *fsi,PetscReal A, 
			       PetscReal w, PetscReal dt)
{
  PetscReal  pi = 3.141592653589793;
  PetscReal t;
  t = (ti-11050)*dt;

  fsi->S_ang_n[0]= A*sin(w*t);
  fsi->S_ang_n[1]= A*w*cos(w*t);

  return(0);
}


/* ==================================================================================             */
PetscErrorCode Elmt_Move_FSI_TRANS(FSInfo *FSinfo, IBMNodes *ibm)
{
  int n_v = ibm->n_v, n_elmt = ibm->n_elmt;
  int i;

/*   for (i=0; i<n_v; i++) { */
/*     ibm->x_bp_o[i] = ibm->x_bp[i]; */
/*     ibm->y_bp_o[i] = ibm->y_bp[i]; */
/*     ibm->z_bp_o[i] = ibm->z_bp[i]; */
/*   } */

  PetscPrintf(PETSC_COMM_WORLD, "MOVE BODY x: %le  y:%le z:%le\n", FSinfo->S_new[0],FSinfo->S_new[2],FSinfo->S_new[4]);

  for (i=0; i<n_v; i++) {
    // change for stat case 4/9/06 iman
    ibm->x_bp[i] = ibm->x_bp0[i]+(FSinfo->S_new[0]);//-FSinfo->S_real[0]);
    ibm->y_bp[i] = ibm->y_bp0[i]+(FSinfo->S_new[2]);//-FSinfo->S_real[2]);
    ibm->z_bp[i] = ibm->z_bp0[i]+(FSinfo->S_new[4]);//-FSinfo->S_real[4]);
  }
  
  for (i=0; i<n_v; i++) {
    ibm->u[i].x = FSinfo->S_new[1];
    ibm->u[i].y = FSinfo->S_new[3];
    ibm->u[i].z = FSinfo->S_new[5];
  }

  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank) {
    if (ti == (ti/tiout)*tiout) {
      FILE *f;
      char filen[80];
      sprintf(filen, "surface%3.3d.dat",(int)ti);
      f = fopen(filen, "w");
      PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=x,y,z\n");
      PetscFPrintf(PETSC_COMM_WORLD, f, "ZONE T='TRIANGLES', N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n", n_v, n_elmt);
      for (i=0; i<n_v; i++) {

	/*    ibm->x_bp[i] = ibm->x_bp0[i];
	      ibm->y_bp[i] = ibm->y_bp0[i];
	      ibm->z_bp[i] = ibm->z_bp0[i] + z0;*/
	PetscFPrintf(PETSC_COMM_WORLD, f, "%e %e %e\n", ibm->x_bp[i], ibm->y_bp[i], ibm->z_bp[i]);
      }
      for (i=0; i<n_elmt; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d %d\n", ibm->nv1[i]+1, ibm->nv2[i]+1, ibm->nv3[i]+1);
      }
      fclose(f);
    }
  }

  return(0);
}
/* ==================================================================================             */

void rotate_xyz (double ti, double dt, double angvel, double x_c, double y_c, double z_c, double x_bp0, double y_bp0, double z_bp0, double *x_bp, double *y_bp, double *z_bp, double *rot_angle)
{
	int i;
	if(rotdir==0) { // rotate around x-axis
		double rot_x = angvel * dt * ti;
		
		*x_bp = x_bp0;
		*y_bp = y_c + (y_bp0-y_c)*cos(rot_x) - (z_bp0-z_c)*sin(rot_x);
		*z_bp = z_c + (y_bp0-y_c)*sin(rot_x) + (z_bp0-z_c)*cos(rot_x);
		
		*rot_angle = rot_x;
	}
	else if(rotdir==1) { // rotate around y-axis
		double rot_y = angvel * dt * ti;

		*y_bp = y_bp0;
		*z_bp = z_c + (z_bp0-z_c)*cos(rot_y) - (x_bp0-x_c)*sin(rot_y);
		*x_bp = x_c + (z_bp0-z_c)*sin(rot_y) + (x_bp0-x_c)*cos(rot_y);
		
		*rot_angle = rot_y;
	}
	else { // rotate around z-axis
		double rot_z = angvel * dt * ti;
		
		*z_bp = z_bp0;
		*x_bp = x_c + (x_bp0-x_c)*cos(rot_z) - (y_bp0-y_c)*sin(rot_z);
		*y_bp = y_c + (x_bp0-x_c)*sin(rot_z) + (y_bp0-y_c)*cos(rot_z);
		
		
		*rot_angle = rot_z;
	}
}

/* ==================================================================================             */
PetscErrorCode Elmt_Move_FSI_ROT(FSInfo *FSinfo, IBMNodes *ibm, 
				 PetscReal dt, int ibi)
{
  //if(ibi!=0) return 0;

  int n_v = ibm->n_v, n_elmt = ibm->n_elmt;
  PetscReal  x_c=FSinfo->x_c, y_c=FSinfo->y_c, z_c=FSinfo->z_c;
  //PetscReal  rot_x=FSinfo->S_ang_n[0];//-FSinfo->S_ang_o[0];
  int i;
  int n1e, n2e, n3e;
  PetscReal dx12, dy12, dz12, dx13, dy13, dz13, dr;
  PetscReal rx,ry,rz;
  //PetscReal wx=FSinfo->S_ang_n[1], wy=FSinfo->S_ang_n[3], wz=FSinfo->S_ang_n[5];
  
	double rot_angle;
	
	for (i=0; i<n_v; i++) {
	  rotate_xyz (ti+ti_lastsave, dt, angvel, x_c, y_c, z_c, ibm->x_bp0[i], ibm->y_bp0[i], ibm->z_bp0[i], &ibm->x_bp[i], &ibm->y_bp[i], &ibm->z_bp[i], &rot_angle);
	  /*	
	  ibm->x_bp[i] = ibm->x_bp0[i];
	  ibm->y_bp[i] = ibm->y_bp0[i];
	  ibm->z_bp[i] = ibm->z_bp0[i];
	  */
		double x1, y1, z1;
		double x2, y2, z2;
		double tmp, eps=1.e-6;
		rotate_xyz (ti+ti_lastsave-eps, dt, angvel, x_c, y_c, z_c, ibm->x_bp0[i], ibm->y_bp0[i], ibm->z_bp0[i], &x1, &y1, &z1, &tmp);
		rotate_xyz (ti+ti_lastsave+eps, dt, angvel, x_c, y_c, z_c, ibm->x_bp0[i], ibm->y_bp0[i], ibm->z_bp0[i], &x2, &y2, &z2, &tmp);
		ibm->u[i].x = (x2 - x1) / dt / eps * 0.5;
		ibm->u[i].y = (y2 - y1) / dt / eps * 0.5;
		ibm->u[i].z = (z2 - z1) / dt / eps * 0.5;
		
		if(ti==tistart) {
		  rotate_xyz (ti+ti_lastsave-1-eps, dt, angvel, x_c, y_c, z_c, ibm->x_bp0[i], ibm->y_bp0[i], ibm->z_bp0[i], &x1, &y1, &z1, &tmp);
		  rotate_xyz (ti+ti_lastsave-1+eps, dt, angvel, x_c, y_c, z_c, ibm->x_bp0[i], ibm->y_bp0[i], ibm->z_bp0[i], &x2, &y2, &z2, &tmp);
		  ibm->uold[i].x = (x2 - x1) / dt / eps * 0.5;
		  ibm->uold[i].y = (y2 - y1) / dt / eps * 0.5;
		  ibm->uold[i].z = (z2 - z1) / dt / eps * 0.5;
		}
	}
	
  // only rot around x
  //wy=0.;wz=0.;
	if(rotdir==0) FSinfo->S_ang_n[0] = rot_angle;
	if(rotdir==1) FSinfo->S_ang_n[2] = rot_angle;
	if(rotdir==2) FSinfo->S_ang_n[4] = rot_angle;
	
	if(rotdir==0) FSinfo->S_ang_n[1] = angvel;
	if(rotdir==1) FSinfo->S_ang_n[3] = angvel;
	if(rotdir==2) FSinfo->S_ang_n[5] = angvel;
	
  //PetscPrintf(PETSC_COMM_WORLD, "Body Rotate: %le %le %le, Center %le %le %le\n", rot_x, rot_y, rot_z, x_c,y_c,z_c);
  
  
     
  for (i=0; i<n_elmt; i++) {

    n1e = ibm->nv1[i]; n2e =ibm->nv2[i]; n3e =ibm->nv3[i];
    dx12 = ibm->x_bp[n2e] - ibm->x_bp[n1e]; 
    dy12 = ibm->y_bp[n2e] - ibm->y_bp[n1e]; 
    dz12 = ibm->z_bp[n2e] - ibm->z_bp[n1e]; 
    
    dx13 = ibm->x_bp[n3e] - ibm->x_bp[n1e]; 
    dy13 = ibm->y_bp[n3e] - ibm->y_bp[n1e]; 
    dz13 = ibm->z_bp[n3e] - ibm->z_bp[n1e]; 
    
    ibm->nf_x[i] = dy12 * dz13 - dz12 * dy13;
    ibm->nf_y[i] = -dx12 * dz13 + dz12 * dx13;
    ibm->nf_z[i] = dx12 * dy13 - dy12 * dx13;

    dr = sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i] + 
	      ibm->nf_z[i]*ibm->nf_z[i]);
    
    ibm->nf_x[i] /=dr; ibm->nf_y[i]/=dr; ibm->nf_z[i]/=dr;
    
    // Addedd 5/7/06 iman
    // ns = nf x k
    if ((((1.-ibm->nf_z[i])<=1e-6 )&&((-1.+ibm->nf_z[i])<1e-6))||
	(((ibm->nf_z[i]+1.)<=1e-6 )&&((-1.-ibm->nf_z[i])<1e-6))) {
      ibm->ns_x[i] = 1.;
      ibm->ns_y[i] = 0.;
      ibm->ns_z[i] = 0 ;

      // nt = ns x nf
      ibm->nt_x[i] = 0.;
      ibm->nt_y[i] = 1.;
      ibm->nt_z[i] = 0.;
      } else {
      ibm->ns_x[i] =  ibm->nf_y[i]/ sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i]);
      ibm->ns_y[i] = -ibm->nf_x[i]/ sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i]);
      ibm->ns_z[i] = 0 ;

      // nt = ns x nf
      ibm->nt_x[i] = -ibm->nf_x[i]*ibm->nf_z[i]/ sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i]);
      ibm->nt_y[i] = -ibm->nf_y[i]*ibm->nf_z[i]/ sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i]);
      ibm->nt_z[i] = sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i]);
      }
    //Added 4/1/06 iman
    ibm->dA[i] = dr/2.;

    // Added 6/4/06 iman
    // Calc the center of the element
    ibm->cent_x[i]= (ibm->x_bp[n1e]+ibm->x_bp[n2e]+ibm->x_bp[n3e])/3.;
    ibm->cent_y[i]= (ibm->y_bp[n1e]+ibm->y_bp[n2e]+ibm->y_bp[n3e])/3.;
    ibm->cent_z[i]= (ibm->z_bp[n1e]+ibm->z_bp[n2e]+ibm->z_bp[n3e])/3.;

  }
  


  for (i=0; i<n_v; i++) {
    rx = ibm->x_bp[i]-x_c;
    ry = ibm->y_bp[i]-y_c;
    rz = ibm->z_bp[i]-z_c;  
    //if(ti==tistart) 
	{    /*
		ibm->u[i].x =   ry*wz-wy*rz  ;
		ibm->u[i].y =-( rx*wz-wx*rz );
		ibm->u[i].z =   rx*wy-wx*ry  ;
		*/
	}
  }
  
  
  return(0);
}

PetscErrorCode ibm_surface_out_with_pressure(IBMNodes *ibm, int ibi)
{
	int n_v = ibm->n_v, n_elmt = ibm->n_elmt;
	std::vector<double> pressure_buf(n_elmt);
	
	/*MPI_Allreduce(&ibm->pressure[0], &pressure_buf[0], n_elmt, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);*/

	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  
	if (!rank) {
		if (ti == (ti/tiout)*tiout) {
			FILE *f;
			char filen[80];
			sprintf(filen, "surface%06d_%1d.dat",(int)ti,ibi);
			f = fopen(filen, "w");
			
			int N_block=100*1024*1024; // 100Mb
			char str[256];
			char carriage_return = '\n';
			std::vector<char> large_buffer;
			
			//PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=x,y,z,u,v,w,p\n");
			//PetscFPrintf(PETSC_COMM_WORLD, f, "ZONE T='TRIANGLES', N=%d, E=%d, F=FEBLOCK, ET=TRIANGLE, VARLOCATION=([1-6]=NODAL,[7]=CELLCENTERED)\n", n_v, n_elmt);
			
			sprintf(str, "Variables=x,y,z,u,v,w,p\n");
			str_to_buffer(str, large_buffer);
			large_buffer.push_back(carriage_return);
			
			sprintf(str, "ZONE T='TRIANGLES', N=%d, E=%d, F=FEBLOCK, ET=TRIANGLE, VARLOCATION=([1-6]=NODAL,[7]=CELLCENTERED)\n", n_v, n_elmt);
			str_to_buffer(str, large_buffer);
			large_buffer.push_back(carriage_return);
			
			for (int i=0; i<n_v; i++) {
				//PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->x_bp[i]);
				sprintf(str, "%e ", ibm->x_bp[i]);
				str_to_buffer(str, large_buffer);
				if( (i+1)%10==0 || i==n_v-1) large_buffer.push_back(carriage_return);
				
				if(large_buffer.size()>N_block) {
					fwrite(&large_buffer[0], sizeof(char), large_buffer.size(), f);
					large_buffer.resize(0);
				}
			}
			for (int i=0; i<n_v; i++) {
				//PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->y_bp[i]);
				sprintf(str, "%e ", ibm->y_bp[i]);
				str_to_buffer(str, large_buffer);
				if( (i+1)%10==0 || i==n_v-1) large_buffer.push_back(carriage_return);
				
				if(large_buffer.size()>N_block) {
					fwrite(&large_buffer[0], sizeof(char), large_buffer.size(), f);
					large_buffer.resize(0);
				}
			}
			for (int i=0; i<n_v; i++) {
				//PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->z_bp[i]);
				sprintf(str, "%e ", ibm->z_bp[i]);
				str_to_buffer(str, large_buffer);
				if( (i+1)%10==0 || i==n_v-1) large_buffer.push_back(carriage_return);
				
				if(large_buffer.size()>N_block) {
					fwrite(&large_buffer[0], sizeof(char), large_buffer.size(), f);
					large_buffer.resize(0);
				}
			}
			for (int i=0; i<n_v; i++) {
				//PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->u[i].x);
				sprintf(str, "%e ", ibm->u[i].x);
				str_to_buffer(str, large_buffer);
				if( (i+1)%10==0 || i==n_v-1) large_buffer.push_back(carriage_return);
				
				if(large_buffer.size()>N_block) {
					fwrite(&large_buffer[0], sizeof(char), large_buffer.size(), f);
					large_buffer.resize(0);
				}
			}
			for (int i=0; i<n_v; i++) {
				//PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->u[i].y);
				sprintf(str, "%e ", ibm->u[i].y);
				str_to_buffer(str, large_buffer);
				if( (i+1)%10==0 || i==n_v-1) large_buffer.push_back(carriage_return);
				
				if(large_buffer.size()>N_block) {
					fwrite(&large_buffer[0], sizeof(char), large_buffer.size(), f);
					large_buffer.resize(0);
				}
			}
			for (int i=0; i<n_v; i++) {
				//PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->u[i].z);
				sprintf(str, "%e ", ibm->u[i].z);
				str_to_buffer(str, large_buffer);
				if( (i+1)%10==0 || i==n_v-1) large_buffer.push_back(carriage_return);
				
				if(large_buffer.size()>N_block) {
					fwrite(&large_buffer[0], sizeof(char), large_buffer.size(), f);
					large_buffer.resize(0);
				}
			}
			//for (i=0; i<n_elmt; i++) PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->nf_x[i]);
			//for (i=0; i<n_elmt; i++) PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->nf_y[i]);
			//for (i=0; i<n_elmt; i++) PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->nf_z[i]);
			for (int i=0; i<n_elmt; i++) {
				//PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", pressure_buf[i]);
				sprintf(str, "%e ", pressure_buf[i]);
				str_to_buffer(str, large_buffer);
				if( (i+1)%10==0 || i==n_elmt-1) large_buffer.push_back(carriage_return);
				
				if(large_buffer.size()>N_block) {
					fwrite(&large_buffer[0], sizeof(char), large_buffer.size(), f);
					large_buffer.resize(0);
				}
			}
			for (int i=0; i<n_elmt; i++) {
				//PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d %d\n", ibm->nv1[i]+1, ibm->nv2[i]+1, ibm->nv3[i]+1);
				sprintf(str, "%d %d %d\n", ibm->nv1[i]+1, ibm->nv2[i]+1, ibm->nv3[i]+1);
				str_to_buffer(str, large_buffer);
				large_buffer.push_back(carriage_return);
				
				if(large_buffer.size()>N_block) {
					fwrite(&large_buffer[0], sizeof(char), large_buffer.size(), f);
					large_buffer.resize(0);
				}
			}
			
			if(large_buffer.size()) fwrite(&large_buffer[0], sizeof(char), large_buffer.size(), f);
			fclose(f);
		}
	}
	return(0);
}

/* ==================================================================================             */
/* ==================================================================================             */

PetscErrorCode CollisionDetectionOfCylinders(FSInfo *fsi, 
					     int NumberOfBodies)

{
  PetscReal     x_c,y_c,z_c;
  PetscReal     x_c2,y_c2,z_c2;
  PetscReal     l_c;
  PetscReal     n_x,n_y,n_z; //collision direction
  PetscReal     v_x=0.,v_y=0.,v_z=0.;
  PetscReal     v_x2=0.,v_y2=0.,v_z2=0.;
  PetscReal     v_n1, v_t1; //vel in collision direction
  PetscReal     v_n2, v_t2;

  int	ibi,ibi2;

  for (ibi=0;ibi<NumberOfBodies;ibi++) {
    for (ibi2=ibi+1;ibi2<NumberOfBodies;ibi2++){
      
      x_c=fsi[ibi].x_c  ; y_c=fsi[ibi].y_c  ; z_c=fsi[ibi].z_c;
      x_c2=fsi[ibi2].x_c; y_c2=fsi[ibi2].y_c; z_c2=fsi[ibi2].z_c;

      l_c=sqrt((x_c-x_c2)*(x_c-x_c2) + 
	       (y_c-y_c2)*(y_c-y_c2) +
	       (z_c-z_c2)*(z_c-z_c2));

      if (l_c < 1.) { // collission !!!
	PetscPrintf(PETSC_COMM_WORLD, "Collision Detected!!!! cylinder %d with %d\n", ibi, ibi2);
	
	// Collision Direction
	n_x = 0.;//(x_c-x_c2)/l_c;
	n_y = (y_c-y_c2)/l_c;
	n_z = (z_c-z_c2)/l_c;
	
	/* Move the 2nd cyl to the 1D distance of 1st Cyl */
	fsi[ibi2].x_c= x_c + n_x;
	fsi[ibi2].y_c= y_c + n_y;
	fsi[ibi2].z_c= z_c + n_z;

	/* Change the Vel to the collsion Vel! */
	/*

	v'.t=v.t
	v'2.t=v2.t
        
	v'.n=v2.n
        v'2.n=v.n

        ' is the new vel (after collision)
        n is the collision direction vector
        t is the bi-normal of collision direction vector

	*/

	v_x=fsi[ibi].S_new[1];
	v_y=fsi[ibi].S_new[3];
	v_z=fsi[ibi].S_new[5];
	v_n1 = v_x*n_x + v_y*n_y + v_z*n_z;
	v_t1 = v_x*n_x + v_y*n_z - v_z*n_y;
	  
	v_x2=fsi[ibi2].S_new[1];
	v_y2=fsi[ibi2].S_new[3];
	v_z2=fsi[ibi2].S_new[5];
	v_n2 = v_x2*n_x + v_y2*n_y + v_z2*n_z;
	v_t2 = v_x2*n_x + v_y2*n_z - v_z2*n_y;

	PetscPrintf(PETSC_COMM_WORLD,"          Velocity!!!! cyl1 %le %le %le  cyl2 %le %le %le\n", v_x,v_y,v_z,v_x2,v_y2,v_z2);
/* 	PetscBarrier(PETSC_NULL); */

	v_x = v_n2*n_x + v_t1 *n_x;
	v_y = v_n2*n_y + v_t1 *n_z;
	v_z = v_n2*n_z - v_t1 *n_y;

	v_x2 = v_n1*n_x + v_t2 *n_x;
	v_y2 = v_n1*n_y + v_t2 *n_z;
	v_z2 = v_n1*n_z - v_t2 *n_y;

	fsi[ibi].S_new[1]=v_x;
	fsi[ibi].S_new[3]=v_y;
	fsi[ibi].S_new[5]=v_z;

	fsi[ibi2].S_new[1]=v_x2;
	fsi[ibi2].S_new[3]=v_y2;
	fsi[ibi2].S_new[5]=v_z2;

	PetscPrintf(PETSC_COMM_WORLD, "Collision Velocity!!!! cyl1 %le %le %le  cyl2 %le %le %le\n", v_x,v_y,v_z,v_x2,v_y2,v_z2);

	
      }
    }
  }
  return(0);
}

/* ==================================================================================             */
/* ==================================================================================             */

PetscErrorCode SwingCylinder(FSInfo *fsi, IBMNodes *ibm)
{
  int n_v = ibm->n_v, n_elmt = ibm->n_elmt;
  PetscReal  x_c=fsi->x_c, y_c=fsi->y_c, z_c=fsi->z_c;
  PetscReal  rot_x=fsi->S_ang_n[0], rot_y=0, rot_z=fsi->S_ang_n[4];//fsi->S_ang_n[2]
  int i;
  int n1e, n2e, n3e;
  PetscReal dx12, dy12, dz12, dx13, dy13, dz13, dr;
  PetscReal rx,ry,rz;
  PetscReal wx=fsi->S_ang_n[1], wy=fsi->S_ang_n[3], wz=fsi->S_ang_n[5];
  PetscReal lrot,denom;

  denom=sqrt(1.+tan(rot_z)*tan(rot_z)+tan(rot_y)*tan(rot_y));

  PetscPrintf(PETSC_COMM_WORLD, "Body Rotate y,z: %le %le Center %le %le %le %le max z y %le %le\n",rot_y,rot_z, x_c,y_c,z_c, denom, 31/denom*tan(rot_y),-31/denom*tan(rot_z));

  for (i=0; i<n_v; i++) {          
    lrot= ibm->x_bp0[i] / denom; 
    ibm->x_bp[i] = lrot;
    ibm->z_bp[i] = ibm->z_bp0[i]+ lrot*tan(rot_y);
    ibm->y_bp[i] = ibm->y_bp0[i]- lrot*tan(rot_z);
  }

  /*   calc the new normal and assign vel */
  for (i=0; i<n_elmt; i++) {
    n1e = ibm->nv1[i]; n2e =ibm->nv2[i]; n3e =ibm->nv3[i];
    dx12 = ibm->x_bp[n2e] - ibm->x_bp[n1e]; 
    dy12 = ibm->y_bp[n2e] - ibm->y_bp[n1e]; 
    dz12 = ibm->z_bp[n2e] - ibm->z_bp[n1e]; 
    
    dx13 = ibm->x_bp[n3e] - ibm->x_bp[n1e]; 
    dy13 = ibm->y_bp[n3e] - ibm->y_bp[n1e]; 
    dz13 = ibm->z_bp[n3e] - ibm->z_bp[n1e]; 
    
    ibm->nf_x[i] = dy12 * dz13 - dz12 * dy13;
    ibm->nf_y[i] = -dx12 * dz13 + dz12 * dx13;
    ibm->nf_z[i] = dx12 * dy13 - dy12 * dx13;

    dr = sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i] + 
	      ibm->nf_z[i]*ibm->nf_z[i]);
    
    ibm->nf_x[i] /=dr; ibm->nf_y[i]/=dr; ibm->nf_z[i]/=dr;
    
    // ns = nf x k
    if ((((1.-ibm->nf_z[i])<=1e-6 )&&((-1.+ibm->nf_z[i])<1e-6))||
	(((ibm->nf_z[i]+1.)<=1e-6 )&&((-1.-ibm->nf_z[i])<1e-6))) {
      ibm->ns_x[i] = 1.;     
      ibm->ns_y[i] = 0.;     
      ibm->ns_z[i] = 0 ;

      ibm->nt_x[i] = 0.;
      ibm->nt_y[i] = 1.;
      ibm->nt_z[i] = 0.;
      } else {
      ibm->ns_x[i] =  ibm->nf_y[i]/ sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i]);      
      ibm->ns_y[i] = -ibm->nf_x[i]/ sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i]);     
      ibm->ns_z[i] = 0 ;

      // nt = ns x nf
      ibm->nt_x[i] = -ibm->nf_x[i]*ibm->nf_z[i]/ sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i]);
      ibm->nt_y[i] = -ibm->nf_y[i]*ibm->nf_z[i]/ sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i]);
      ibm->nt_z[i] = sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i]);
      }

  }

  // 2nd order approx. 
  if (ti>0) {
    for (i=0; i<n_v; i++) {
      rx = ibm->x_bp[i]-x_c;
      ry = ibm->y_bp[i]-y_c;
      rz = ibm->z_bp[i]-z_c;      
      ibm->u[i].x =   ry*wz-wy*rz  ;
      ibm->u[i].y =-( rx*wz-wx*rz );
      ibm->u[i].z =   rx*wy-wx*ry  ;      
    }
  }
  else {
    for (i=0; i<n_v; i++) {
      ibm->u[i].x = 0.;
      ibm->u[i].y = 0.;
      ibm->u[i].z = 0.;
    }
  }

  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank) {
    if (ti == (ti/tiout)*tiout) {
      FILE *f;
      char filen[80];
      sprintf(filen, "surface%3.3d.dat",(int)ti);
      f = fopen(filen, "w");
      PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=x,y,z\n");
      PetscFPrintf(PETSC_COMM_WORLD, f, "ZONE T='TRIANGLES', N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n", n_v, n_elmt);
      for (i=0; i<n_v; i++) {

	/*    ibm->x_bp[i] = ibm->x_bp0[i];
	      ibm->y_bp[i] = ibm->y_bp0[i];
	      ibm->z_bp[i] = ibm->z_bp0[i] + z0;*/
	PetscFPrintf(PETSC_COMM_WORLD, f, "%e %e %e\n", ibm->x_bp[i], ibm->y_bp[i], ibm->z_bp[i]);
      }
      for (i=0; i<n_elmt; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d %d\n", ibm->nv1[i]+1, ibm->nv2[i]+1, ibm->nv3[i]+1);
      }
      fclose(f);
    }
  }
  return(0);
}


// xiaolei
void Vel_rotate (double angvel, double xc, double yc, double zc, double xb, double yb, double zb, double *ux, double *uy, double *uz)
{
	int i;
	if(rotdir==0) { // rotate around x-axis
		double rr=sqrt(pow(yb-yc,2)+pow(zb-zc,2))+1.e-19;
		double ny=(yb-yc)/rr;
		double nz=(zb-zc)/rr;
		double ut=angvel*rr;

		*ux = 0.0;
		*uy = -ut*nz;
		*uz = ut*ny;
	}
	else if(rotdir==1) { // rotate around y-axis
		double rr=sqrt(pow(xb-xc,2)+pow(zb-zc,2))+1.e-19;
		double nx=(xb-xc)/rr;
		double nz=(zb-zc)/rr;
		double ut=angvel*rr;

		*ux = ut*nz;
		*uy = 0.0;
		*uz = -ut*nx;
	}
	else { // rotate around z-axis
		double rr=sqrt(pow(xb-xc,2)+pow(yb-yc,2))+1.e-19;
		double nx=(xb-xc)/rr;
		double ny=(yb-yc)/rr;
		double ut=angvel*rr;

		*ux = -ut*ny;
		*uy = ut*nx;
		*uz = 0.0;

	}
}



/* ==================================================================================             */
/* ==================================================================================             */
