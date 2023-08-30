#include "variables.h"

double threshold=0.1;

#define INLET 5
#define OUTLET 4
#define SOLIDWALL 1	
#define SYMMETRIC 3
#define FARFIELD 6

PetscReal InletInterpolation(PetscReal r, UserCtx *user)
{
  PetscInt i;
  PetscReal temp;
  PetscInt tstep;
  tstep = ti - 1000 * (ti/1000);

/*   return(0); */

/*   return (2*(1-(r*r)*4)); */
/*   if (ti>0)  */
/*   return(1.); */
/*   else */
/*    return 1.; */
  return(Flux_in);

/*   if (r>1.) { */
/*     temp = 0.; */
/*     return (temp); */
/*   } */
/*   else { */
/*     temp = 2 * (1 - (r*2) * (r*2)); */
/*     return(temp); */
/*   } */
  
/*   if (ti==0) { */
/*     temp = 0.; */
/*     return(temp); */
/*   } */
  if (r>1.) {
    temp = user->uinr[99][tstep];
    //    PetscPrintf(PETSC_COMM_WORLD, "Inlet %e", temp);
    return(temp);
  }
  for (i=0; i<100; i++) {
    if (r>= (user->r[i]) && r< (user->r[i+1])) {
      temp = user->uinr[i][tstep] + (user->uinr[i+1][tstep] - user->uinr[i][tstep]) *
	(r-user->r[i]) / (user->r[i+1]-user->r[i]);
      //      PetscPrintf(PETSC_COMM_WORLD, "Inlet %e", temp);
      return (temp);
    }
  }   
}

double inletArea=0;
int k_area_allocated=0;

void Calc_Inlet_Area(UserCtx *user)
{
	inletArea=0;
		
	PetscInt i, j, k;

	Cmpnts	***csi, ***eta, ***zet;
	
	DM da = user->da, fda = user->fda;
	DMDALocalInfo	info = user->info;
	PetscInt	xs = info.xs, xe = info.xs + info.xm;
	PetscInt  	ys = info.ys, ye = info.ys + info.ym;
	PetscInt	zs = info.zs, ze = info.zs + info.zm;
	PetscInt	mx = info.mx, my = info.my, mz = info.mz;
	PetscInt	lxs, lxe, lys, lye, lzs, lze;
 
	PetscReal	***nvert, ***level, ***rho, ***aj;
	
	
	lxs = xs; lxe = xe;
	lys = ys; lye = ye;
	lzs = zs; lze = ze;
  
	if (xs==0) lxs = xs+1;
	if (ys==0) lys = ys+1;
	if (zs==0) lzs = zs+1;
  
	if (xe==mx) lxe = xe-1;
	if (ye==my) lye = ye-1;
	if (ze==mz) lze = ze-1;

	if(levelset)  {
		DMDAVecGetArray(da, user->lLevelset, &level);
		DMDAVecGetArray (da, user->lDensity, &rho);
	}
	DMDAVecGetArray(fda, user->lCsi, &csi);
	DMDAVecGetArray(fda, user->lEta, &eta);
	DMDAVecGetArray(fda, user->lZet, &zet);
	DMDAVecGetArray(da, user->lNvert, &nvert);	//seokkoo 
	DMDAVecGetArray(da, user->lAj, &aj);

	if(!k_area_allocated) {
		k_area_allocated=1;
		user->k_area = new double [mz];
		user->k_area_ibnode = new double [mz];
	}
        
	std::vector<double> lArea(mz), lArea_ibm(mz);
	
	std::fill ( lArea.begin(), lArea.end(), 0 );
	std::fill ( lArea_ibm.begin(), lArea_ibm.end(), 0 );

        for(k=lzs; k<lze; k++)
	for (j=ys; j<ye; j++)
	for (i=xs; i<xe; i++) {
		if(j>=1 && j<=my-2 && i>=1 && i<=mx-2) {
			double area = sqrt( zet[k][j][i].x*zet[k][j][i].x + zet[k][j][i].y*zet[k][j][i].y + zet[k][j][i].z*zet[k][j][i].z );
			if (nvert[k+1][j][i]+nvert[k][j][i] < 0.1) lArea[k] += area;
			else if (nvert[k+1][j][i]+nvert[k][j][i] < 2.1) lArea_ibm[k] += area;
		}
	}

        MPI_Allreduce( &lArea[0], &user->k_area[0], mz, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
        MPI_Allreduce( &lArea_ibm[0], &user->k_area_ibnode[0], mz, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
	
	user->k_area[0] = user->k_area[1];

	user->mean_k_area = 0;
	user->mean_k_area_ibnode = 0;
	
	for (k=1; k<=mz-2; k++) {
		user->mean_k_area += user->k_area[k];
		user->mean_k_area_ibnode += user->k_area_ibnode[k];
	}
	
        user->mean_k_area /= (double)(mz-2);
        user->mean_k_area_ibnode /= (double)(mz-2);
	
	if (user->bctype[4] == INLET || k_periodic || kk_periodic) {
		double lArea=0;
		if (zs==0) {
			k = 1;
			for (j=lys; j<lye; j++)
			for (i=lxs; i<lxe; i++) {
				double k_area = sqrt( zet[k][j][i].x*zet[k][j][i].x + zet[k][j][i].y*zet[k][j][i].y + zet[k][j][i].z*zet[k][j][i].z );
				
				if (nvert[k][j][i] < threshold) {				  
					if( levelset ) {
						double vf=1.0;
						   
						if(levelset==1) {
							double dx = levelset_thickness ( aj[k][j][i], csi[k][j][i], eta[k][j][i], zet[k][j][i] );
							vf = H(level[k][j][i], dx);
						}
						else if(levelset==2) {
							vf = level[k][j][i];
						}
						
						lArea += k_area * vf;
					}
					else {
						lArea += k_area;
					}
				}
			}
		}
		GlobalSum_All(&lArea, &inletArea, PETSC_COMM_WORLD);    
	}
	
	if (user->bctype[5] == INLET) {
		double lArea=0;
		if (ze==mz) {
			k = mz-2;
			for (j=lys; j<lye; j++)
			for (i=lxs; i<lxe; i++) {
				double k_area = sqrt( zet[k][j][i].x*zet[k][j][i].x + zet[k][j][i].y*zet[k][j][i].y + zet[k][j][i].z*zet[k][j][i].z );

				if (nvert[k][j][i] < threshold) {				  
					if( levelset ) {
						double vf=1.0;
						if(levelset==1) {
							double dx = levelset_thickness ( aj[k][j][i], csi[k][j][i], eta[k][j][i], zet[k][j][i] );
							vf = H(level[k][j][i], dx);
						}
						else if(levelset==2) {
							vf = level[k][j][i];
						}
						lArea += k_area * vf;
					}
					else {
						lArea += k_area;
					}
				}
			}
		}
		GlobalSum_All(&lArea, &inletArea, PETSC_COMM_WORLD);    
	}
       
	DMDAVecRestoreArray(fda, user->lCsi, &csi);
	DMDAVecRestoreArray(fda, user->lEta, &eta);
	DMDAVecRestoreArray(fda, user->lZet, &zet);
	DMDAVecRestoreArray(da, user->lNvert, &nvert);	//seokkoo 
	DMDAVecRestoreArray(da, user->lAj, &aj);
	if(levelset) {
		DMDAVecRestoreArray (da, user->lLevelset, &level);
		DMDAVecRestoreArray (da, user->lDensity, &rho);
	}
	
	PetscPrintf(PETSC_COMM_WORLD, "\n*** (Fluid) Inlet Area:%f, (Fluid) Inlet Flux: %f\n\n", inletArea, inlet_flux);
};

PetscErrorCode InflowFlux(UserCtx *user) 
{
	
  
	PetscInt i, j, k;
	PetscReal r, uin, xc, yc;
	Vec Coor;
	Cmpnts	***ucont, ***ubcs, ***ucat, ***coor, ***csi, ***eta, ***zet;
	Cmpnts ***icsi;
	
	PetscReal  /*H=4.1,*/ Umax=1.5;

	DM da = user->da, fda = user->fda;
	DMDALocalInfo	info = user->info;
	PetscInt	xs = info.xs, xe = info.xs + info.xm;
	PetscInt  	ys = info.ys, ye = info.ys + info.ym;
	PetscInt	zs = info.zs, ze = info.zs + info.zm;
	PetscInt	mx = info.mx, my = info.my, mz = info.mz;
	PetscInt	lxs, lxe, lys, lye, lzs, lze;
	
	PetscReal	***nvert, ***level, ***aj;	//seokkoo
	
	DMDAVecGetArray(da, user->lNvert, &nvert);	//seokkoo 
	if(levelset) DMDAVecGetArray(da, user->lLevelset, &level);
	
	lxs = xs; lxe = xe;
	lys = ys; lye = ye;
	lzs = zs; lze = ze;
  
	if (xs==0) lxs = xs+1;
	if (ys==0) lys = ys+1;
	if (zs==0) lzs = zs+1;
  
	if (xe==mx) lxe = xe-1;
	if (ye==my) lye = ye-1;
	if (ze==mz) lze = ze-1;
  
	DMDAGetGhostedCoordinates(da, &Coor);
	DMDAVecGetArray(fda, Coor, &coor);
	DMDAVecGetArray(fda, user->Ucont, &ucont);
	DMDAVecGetArray(fda, user->Bcs.Ubcs, &ubcs);
	DMDAVecGetArray(fda, user->Ucat,  &ucat);
  
	DMDAVecGetArray(fda, user->lCsi,  &csi);
	DMDAVecGetArray(fda, user->lEta,  &eta);
	DMDAVecGetArray(fda, user->lZet,  &zet);
	DMDAVecGetArray(da, user->lAj,  &aj);

	DMDAVecGetArray(fda, user->lICsi,  &icsi);
  
	PetscReal FluxIn=0., FluxIn_gas=0.;
	
	double lFluxIn0=0, sumFluxIn0=0;
	double lFluxIn1=0, sumFluxIn1=0;
  
	srand( time(NULL)) ;	// seokkoo
	int fluct=0;
      
	if (user->bctype[4] == INLET) {
		if (zs==0) {
			k = 0;
			for (j=lys; j<lye; j++)
			for (i=lxs; i<lxe; i++) {
				xc = (coor[k+1][j][i].x + coor[k+1][j-1][i].x + coor[k+1][j][i-1].x + coor[k+1][j-1][i-1].x) * 0.25;
				yc = (coor[k+1][j][i].y + coor[k+1][j-1][i].y + coor[k+1][j][i-1].y + coor[k+1][j-1][i-1].y) * 0.25;
				r = sqrt(xc * xc + yc * yc);
				double area = sqrt( zet[k][j][i].x*zet[k][j][i].x + zet[k][j][i].y*zet[k][j][i].y + zet[k][j][i].z*zet[k][j][i].z );
				
				if (inletprofile == 0) uin = InletInterpolation(r, user);
				else if (inletprofile==1) {
					if(inlet_flux<0) uin=1.;
					else uin = inlet_flux/inletArea;
				}
				else if (inletprofile==2) {	// uniform flow with noise
					if(inlet_flux<0) uin=1.;
					else uin = inlet_flux/inletArea;
					fluct=1;
				}
				else if (inletprofile == 3) { // parabolic with channel height 2, bulk = 1
					if(inlet_flux<0) uin=1.;
					else uin = inlet_flux/inletArea;
					
					uin *= 1.5 * (2 * yc - yc * yc);
				}
				else if (inletprofile == 4) { // parabolic with channel height 1, bulk = 1
				  if(inlet_flux<0) uin=1.;
				  else uin = inlet_flux/inletArea;

				  uin *= 6 * (yc - yc * yc);
                                }
				else if (inletprofile==10) {	// Power law for hemisphere case
					//fluct=1;
					double delta = 0.45263, a=5.99;
					if( yc>=delta ) uin=1.0;
					else if(yc<=0) uin=0.0;
					else uin = pow( yc/delta, 1./a );
				} 
				else if (inletprofile==11) {	// backward facing step
					fluct=1;
					yc -= 3;
					double delta=1.2;
					if( 2-fabs(yc)<delta ) uin = pow( (2-fabs(yc))/delta, 1./7 );
					else uin=1;
				} 
				else if (inletprofile==12) {	// pipe shear stress test
					//		fluct=1;
					uin = 2*(1 - pow(r/0.5,2.0) );
				}
				else if (inletprofile==13) {	// periodic channel flow
				}
				else if(inletprofile==14) {	// curved pipe flow (Anwer)
					double R[11] = {0.000,0.062,0.125,0.188,0.250,0.312,0.354,0.399,0.438,0.469,0.500};
					double W[11] = {1.120,1.115,1.110,1.093,1.050,1.010,0.963,0.915,0.825,0.730,0.000};
					
					int ii;
					for(ii=1; ii<11; ii++) {
						if( R[ii]>=r && R[ii-1]<r ) break;
					}
					uin = ( W[ii] - W[ii-1] ) / ( R[ii] - R[ii-1] ) * ( r - R[ii-1] ) + W[ii-1];
					
					if(r>0.5) uin = 0;
				}
				else if(inletprofile==15) {	// round jet (Longmire)
					if(r>0.5) uin = 0;
					else uin = 1;
				}
				else if (inletprofile==16) {	// periodic pipe flow
				}
				else if(inletprofile==20){}// enright test
				else if (inletprofile==21) {    // laminar poiseulle flow profile
				  double w_bulk=1.0;//inlet_flux/SumArea;
				  uin = 1.5 * w_bulk * yc * ( 2 - yc );
				}
				else if (inletprofile == 100) {}	// saved data for LES
				
				else {
					PetscPrintf(PETSC_COMM_SELF, "WRONG INLET PROFILE TYPE!!!! U_in = 0\n");
					uin = 0.;
				}
				
				if(pseudo_periodic || inletprofile==100) {	// pseudo-periodic BC in k-direction
					fluct=0;
				}
				
				if (nvert[k+1][j][i] < threshold) {	//seokkoo
					if(pseudo_periodic || inletprofile==100) {
						double u = user->ucat_plane[j][i].x;
						double v = user->ucat_plane[j][i].y;
						double w = user->ucat_plane[j][i].z;
						
						ucat[k][j][i].x = u;
						ucat[k][j][i].y = v;
						ucat[k][j][i].z = w;
						ucont[k][j][i].z = 0.5*(ucat[k][j][i].x+ucat[k+1][j][i].x) * zet[k][j][i].x + 0.5*(ucat[k][j][i].y+ucat[k+1][j][i].y) * zet[k][j][i].y + 0.5*(ucat[k][j][i].z+ucat[k+1][j][i].z) * zet[k][j][i].z;
					}
					else {
						ucont[k][j][i].z = uin * area;
						ucont[k][j][i].x = 0;
						ucont[k][j][i].y = 0;
						
						Cmpnts u;
						Contra2Cart_single(csi[k][j][i], eta[k][j][i], zet[k][j][i], ucont[k][j][i], &u);
						
						if(fluct) {
						  double f = fluct_rms;            // 1% noise
						  int n1 = rand() % 20000 - 10000;      // RAND_MAX = 65535
						  int n2 = rand() % 20000 - 10000;
						  int n3 = rand() % 20000 - 10000;
						  u.x += ((double)n1)/10000. * f * uin;
						  u.y += ((double)n2)/10000. * f * uin;
						  u.z += ((double)n3)/10000. * f * uin;
						}
						ucat[k][j][i].x = - ucat[k+1][j][i].x + 2*u.x; //important 110312
						ucat[k][j][i].y = - ucat[k+1][j][i].y + 2*u.y;
						ucat[k][j][i].z = - ucat[k+1][j][i].z + 2*u.z;
										
						//ucat[k][j][i] = ucat[k+1][j][i]; // removed 110312
						
					}
					ubcs[k][j][i] = ucat[k][j][i];
				}
				else {
					ucat[k][j][i].z = 0;	//seokkoo
					ubcs[k][j][i].z = 0;
					ucont[k][j][i].z = 0;
				}
				/*if(j>=0 && j<=my-2 && i>=0 && i<=mx-2) */
				{
					lFluxIn0 += ucont[k][j][i].z;

					if(fluct) {
					  double f = fluct_rms;		// 1% noise
					  int n1 = rand() % 20000 - 10000;	// RAND_MAX = 65535
					  ucont[k][j][i].z *= ( 1 + ((double)n1)/10000. * f );	// uin * (1+-0.xx)
					}
				
					if( levelset ) {
						double vf;
						
						if(levelset==1) {
							double dx = levelset_thickness ( aj[k][j][i], csi[k][j][i], eta[k][j][i], zet[k][j][i] );
							vf = H(level[k][j][i], dx); // xiaolei k+1-->k
						}
						else if(levelset==2) vf = level[k][j][i]; // xiaolei k+1-->k

						lFluxIn1 += ucont[k][j][i].z * vf;
					}
					else 
					lFluxIn1 += ucont[k][j][i].z;
				}
			}
		}
		
		GlobalSum_All(&lFluxIn0, &sumFluxIn0, PETSC_COMM_WORLD);
		GlobalSum_All(&lFluxIn1, &sumFluxIn1, PETSC_COMM_WORLD);
		
		if(inlet_flux<0) sumFluxIn0=1*inletArea;
		else sumFluxIn0=inlet_flux;

		if(pseudo_periodic || inletprofile==100 || inletprofile==14) {
		  //if(inlet_flux<0) inlet_flux=1.0*inletArea;	// bulk velocity is unity.
			PetscPrintf(PETSC_COMM_WORLD,  "\nConstant Flux is %f !\n\n", inlet_flux);
			//sumFluxIn0=inlet_flux;
		}
		
		PetscPrintf(PETSC_COMM_WORLD, "\n***** Fluxin0:%f, Fluxin1:%f, Area:%f\n\n", sumFluxIn0, sumFluxIn1, inletArea);
		
		// xiaolei
		if (zs==0 && fabs(sumFluxIn0 - sumFluxIn1)>1.e-9 && inlet_flux>0/* && (fluct || pseudo_periodic || inletprofile==100 || inletprofile==14) && inlet_flux>0*/ ) {
			
			PetscPrintf(PETSC_COMM_WORLD, "The inflow is corrected to Flux=%f \n\n", sumFluxIn0);
			k = 0;
			for (j=lys; j<lye; j++)
			for (i=lxs; i<lxe; i++) {
				//double A = zet[k][j][i].z;
				// xiaolei add vf
				double vf=1.0;

				if( levelset ) {
					if(levelset==1) {
						double dx = levelset_thickness ( aj[k][j][i], csi[k][j][i], eta[k][j][i], zet[k][j][i] );
						vf = H(level[k][j][i], dx);
					}
					else if(levelset==2) vf = level[k+1][j][i];
				}
				// end

			
				double A = sqrt( zet[k][j][i].x*zet[k][j][i].x + zet[k][j][i].y*zet[k][j][i].y + zet[k][j][i].z*zet[k][j][i].z );
				if (nvert[k+1][j][i] < threshold) ucont[k][j][i].z += (sumFluxIn0 - sumFluxIn1) * A * vf / inletArea;
			}
		}
		
		FluxIn = 0;
		if (zs==0) {
			k = 0;
			for (j=lys; j<lye; j++) 
			for (i=lxs; i<lxe; i++) {
				if (nvert[k+1][j][i] < threshold) {
					FluxIn += ucont[k][j][i].z;
				}
			}
		}
		
       }
	else if (user->bctype[5] == INLET) { // made just for subcritical levelset flow
		if (ze==mz) {
			k = mz-1;
			for (j=lys; j<lye; j++)
			for (i=lxs; i<lxe; i++) {
				double area = sqrt( zet[k-1][j][i].x*zet[k-1][j][i].x + zet[k-1][j][i].y*zet[k-1][j][i].y + zet[k-1][j][i].z*zet[k-1][j][i].z );
	  
				if (inletprofile == 1) {
					if(inlet_flux<0) uin=1.;
					else uin = inlet_flux/inletArea;
				}
				/*
				else if (inletprofile == -1) {
					uin = -1.;
				} 
				else if (inletprofile == 2) {
				  //uin = 4.*Umax*yc*(H-yc)/(H*H);//InletInterpolation(r, user);
				} 
				else if (inletprofile == 3) {
					uin = InletInterpolation(r, user);
				} 
				else {
					PetscPrintf(PETSC_COMM_SELF, "WRONG INLET PROFILE TYPE!!!! U_in = 0\n");
					uin = 0.;
				}
				*/
				if(ti==0) {
					ucat[k][j][i].z = uin;
					ubcs[k][j][i].z = uin;
					ucont[k-1][j][i].z = uin * area;
					ucont[k-1][j][i].x = 0;
					ucont[k-1][j][i].y = 0;
				}
				else {
					// treated in implicitsolver.c and correct with outflow scale
				}
				
				{
					lFluxIn0 += ucont[k-1][j][i].z;

					if( levelset ) {
						double vf;
						if(levelset==1) {
							double dx = levelset_thickness ( aj[k][j][i], csi[k][j][i], eta[k][j][i], zet[k][j][i] );
							vf = H(level[k-1][j][i], dx);
							lFluxIn1 += ucont[k-1][j][i].z * vf;
						}
						else if(levelset==2) vf = level[k-1][j][i];
					}
					else lFluxIn1 += ucont[k-1][j][i].z;
				}
			}
		}
		
		GlobalSum_All(&lFluxIn0, &sumFluxIn0, PETSC_COMM_WORLD);
		GlobalSum_All(&lFluxIn1, &sumFluxIn1, PETSC_COMM_WORLD);
		PetscPrintf(PETSC_COMM_WORLD, "\n***** Fluxin0:%f, Fluxin1:%f, Area:%f\n\n", sumFluxIn0, sumFluxIn1, inletArea);
		
		FluxIn = 0;
		if (ze==mz) {
			k = mz-2;
			for (j=lys; j<lye; j++)
			for (i=lxs; i<lxe; i++) {
				if (nvert[k][j][i] < threshold) FluxIn += ucont[k][j][i].z;
			}
		}
	}
        
	else if(user->bctype[0]==11) {
		FluxIn = 0;
		if (xs==0) {
			i = 0;
			for (j=lys; j<lye; j++) 
			for (k=lzs; k<lze; k++) {
				if (nvert[k][j][i+1] < threshold) {
					double zc = ( coor[k][j][i+1].z + coor[k-1][j][i+1].z + coor[k][j-1][i+1].z + coor[k-1][j-1][i+1].z ) * 0.25;
					
					if( zc <= 0 ) {
						double u=0, v=0, w=1.;
						ucont[k][j][i].x = u * icsi[k][j][i].x + v * icsi[k][j][i].y + w * icsi[k][j][i].z;
						
						FluxIn += ucont[k][j][i].x;
					}
					else { }	// outflow
				}
			}
		}
	}
	//exit(0);
	GlobalSum_All(&FluxIn, &FluxInSum, PETSC_COMM_WORLD);
	GlobalSum_All(&FluxIn_gas, &FluxInSum_gas, PETSC_COMM_WORLD);	// 100203
	
	user->FluxInSum = FluxInSum;
  
	DMDAVecRestoreArray(fda, Coor, &coor);
	DMDAVecRestoreArray(fda, user->Ucont, &ucont);
	DMDAVecRestoreArray(fda, user->Bcs.Ubcs, &ubcs);
	DMDAVecRestoreArray(fda, user->Ucat,  &ucat);
  
	DMDAVecRestoreArray(fda, user->lCsi,  &csi);
	DMDAVecRestoreArray(fda, user->lEta,  &eta);
	DMDAVecRestoreArray(fda, user->lZet,  &zet);
	DMDAVecRestoreArray(da, user->lAj,  &aj);

	DMDAVecRestoreArray(fda, user->lICsi,  &icsi);
	
	DMDAVecRestoreArray(da, user->lNvert, &nvert);	//seokkoo
	if(levelset) DMDAVecRestoreArray(da, user->lLevelset, &level);
	
	return 0;
}

PetscErrorCode OutflowFlux(UserCtx *user) {
  
  PetscInt i, j, k;
  PetscReal FluxOut;
  Vec Coor;
  Cmpnts	***ucont, ***ubcs, ***ucat, ***coor;
  


  DM da = user->da, fda = user->fda;
  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;
  
  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;
  
  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;
  
  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;
  
  DMDAGetGhostedCoordinates(da, &Coor);
  DMDAVecGetArray(fda, Coor, &coor);
  DMDAVecGetArray(fda, user->Ucont, &ucont);
  DMDAVecGetArray(fda, user->Bcs.Ubcs, &ubcs);
  DMDAVecGetArray(fda, user->Ucat,  &ucat);
 
	PetscReal	***nvert;	//seokkoo
	DMDAVecGetArray(da, user->lNvert, &nvert);	//seokkoo  
	
/*     DMDAVecGetArray(fda, user->Csi,  &csi); */
/*     DMDAVecGetArray(fda, user->Eta,  &eta); */
/*     DMDAVecGetArray(fda, user->Zet,  &zet); */

  FluxOut = 0;
  
  if (user->bctype[5] == 4) {    
    if (ze==mz) {
      k = mz-2;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
		if (nvert[k][j][i] < threshold) //seokkoo
			FluxOut += ucont[k][j][i].z;
	}
      }
    }
    else {
      FluxOut = 0;
    }

  } else if (user->bctype[4] == 4) {    
    if (zs==0) {
      k = 0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
		if (nvert[k+1][j][i] < threshold) //seokkoo
			FluxOut += ucont[k][j][i].z;
	}
      }
    }
    else {
      FluxOut = 0;
    }
  }

  GlobalSum_All(&FluxOut, &FluxOutSum, PETSC_COMM_WORLD);
  user->FluxOutSum = FluxOutSum;

  DMDAVecRestoreArray(fda, Coor, &coor);
  DMDAVecRestoreArray(fda, user->Ucont, &ucont);
  DMDAVecRestoreArray(fda, user->Bcs.Ubcs, &ubcs);
  DMDAVecRestoreArray(fda, user->Ucat,  &ucat);
  
  DMDAVecRestoreArray(da, user->lNvert, &nvert);	//seokkoo
  //  VecDestroy(&Coor);
/*     DMDAVecRestoreArray(fda, user->Csi,  &csi); */
/*     DMDAVecRestoreArray(fda, user->Eta,  &eta); */
/*     DMDAVecRestoreArray(fda, user->Zet,  &zet); */
	return 0;
}

PetscErrorCode Blank_Interface(UserCtx *user) {
  PetscInt counter=0;
  PetscInt ci, cj, ck;
  
  PetscInt i;

  PetscReal ***nvert;

/*   for (bi=0; bi<block_number; bi++) { */
    DM		da =user->da;
    DMDALocalInfo	info = user->info;

    PetscInt	xs = info.xs, xe = info.xs + info.xm;
    PetscInt  	ys = info.ys, ye = info.ys + info.ym;
    PetscInt	zs = info.zs, ze = info.zs + info.zm;
    
    PetscInt	mx = info.mx, my = info.my, mz = info.mz;
    PetscInt	lxs, lxe, lys, lye, lzs, lze;

    lxs = xs; lxe = xe;
    lys = ys; lye = ye;
    lzs = zs; lze = ze;
    
    if (xs==0) lxs = xs+2;
    if (ys==0) lys = ys+2;
    if (zs==0) lzs = zs+2;
    
    if (xe==mx) lxe = xe-1;
    if (ye==my) lye = ye-3;
    if (ze==mz) lze = ze-3;
    
    DMDAVecGetArray(user->da, user->Nvert, &nvert);

    for (i=0; i<user->itfcptsnumber; i++) {
      ci = user->itfcI[i];
      cj = user->itfcJ[i];
      ck = user->itfcK[i];

      if (ci>xs && ci<lxe &&
	  cj>lys && cj<lye &&
	  ck>lzs && ck<lze) {
	counter++;
	nvert[ck][cj][ci]=1.5;

      }

    }

    PetscPrintf(PETSC_COMM_WORLD, "Interface pts blanked!!!! %i,\n", counter);

    DMDAVecRestoreArray(user->da, user->Nvert, &nvert);

    DMGlobalToLocalBegin(da, user->Nvert, INSERT_VALUES, user->lNvert);
    DMGlobalToLocalEnd(da, user->Nvert, INSERT_VALUES, user->lNvert);

    return(0);
}

PetscErrorCode Block_Interface_U(UserCtx *user) {
  PetscInt bi;
  PetscInt ci, cj, ck;
  PetscInt hi, hj, hk, hb;
  PetscInt i, j, k;
  PetscReal x, y, z;

  Vec	hostU, nhostU;
  Cmpnts ***itfc, ***ucat;
  PetscReal *hostu;

/*   for (bi=0; bi<block_number; bi++) { */
/*     VecCreateSeq(PETSC_COMM_SELF, */
/* 		 3*user[bi].info.mx*user[bi].info.my*user[bi].info.mz, */
/* 		 &(user[bi].nhostU)); */
/*   } */
  // First calculate Phi components at grid nodes
  for (bi=0; bi<block_number; bi++) {
    DMDALocalInfo	info = user[bi].info;
    PetscInt	xs = info.xs, xe = info.xs + info.xm;
    PetscInt  	ys = info.ys, ye = info.ys + info.ym;
    PetscInt	zs = info.zs, ze = info.zs + info.zm;
    PetscInt	mx = info.mx, my = info.my, mz = info.mz;
    DMDAVecGetArray(user[bi].fda, user[bi].Itfc, &itfc);
    DMDAVecGetArray(user[bi].fda, user[bi].lUcat, &ucat);

    DMDACreateNaturalVector(user[bi].fda, &nhostU);
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	for (i=xs; i<xe; i++) {
	  if (k<mz-1 && j<my-1 && i<mx-1) {
	    itfc[k][j][i].x = 0.125 * (ucat[k  ][j  ][i  ].x +
				       ucat[k  ][j  ][i+1].x +
				       ucat[k  ][j+1][i  ].x +
				       ucat[k  ][j+1][i+1].x +
				       ucat[k+1][j  ][i  ].x +
				       ucat[k+1][j  ][i+1].x +
				       ucat[k+1][j+1][i  ].x +
				       ucat[k+1][j+1][i+1].x);

	    itfc[k][j][i].y = 0.125 * (ucat[k  ][j  ][i  ].y +
				       ucat[k  ][j  ][i+1].y +
				       ucat[k  ][j+1][i  ].y +
				       ucat[k  ][j+1][i+1].y +
				       ucat[k+1][j  ][i  ].y +
				       ucat[k+1][j  ][i+1].y +
				       ucat[k+1][j+1][i  ].y +
				       ucat[k+1][j+1][i+1].y);
	    itfc[k][j][i].z = 0.125 * (ucat[k  ][j  ][i  ].z +
				       ucat[k  ][j  ][i+1].z +
				       ucat[k  ][j+1][i  ].z +
				       ucat[k  ][j+1][i+1].z +
				       ucat[k+1][j  ][i  ].z +
				       ucat[k+1][j  ][i+1].z +
				       ucat[k+1][j+1][i  ].z +
				       ucat[k+1][j+1][i+1].z);
	  }
	}
      }
    }
    DMDAVecRestoreArray(user[bi].fda, user[bi].Itfc, &itfc);
    DMDAVecRestoreArray(user[bi].fda, user[bi].lUcat, &ucat);

    DMDAGlobalToNaturalBegin(user[bi].fda, user[bi].Itfc, INSERT_VALUES, nhostU);
    DMDAGlobalToNaturalEnd(user[bi].fda, user[bi].Itfc, INSERT_VALUES, nhostU);
    //    DMDALocalToGlobal(user[bi].fda, user[bi].lItfc, INSERT_VALUES, user[bi].Itfc);
    VecScatter tolocalall;
    VecScatterCreateToAll(nhostU, &tolocalall, &(user[bi].nhostU));
    VecScatterBegin(tolocalall, nhostU, user[bi].nhostU, INSERT_VALUES, SCATTER_FORWARD);

    VecScatterEnd(tolocalall, nhostU, user[bi].nhostU, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterDestroy(&tolocalall);
    VecDestroy(&nhostU);
  }

  VecCreateSeq(PETSC_COMM_SELF, 24, &hostU);
  for (bi=0; bi<block_number; bi++) {
    DMDALocalInfo	info = user[bi].info;

    PetscInt	xs = info.xs, xe = info.xs + info.xm;
    PetscInt  	ys = info.ys, ye = info.ys + info.ym;
    PetscInt	zs = info.zs, ze = info.zs + info.zm;
    

    PetscInt    lmx, lmy, lmz;
    DMDAVecGetArray(user[bi].fda, user[bi].lItfc, &itfc);
    for (i=0; i<user[bi].itfcptsnumber; i++) {
      ci = user[bi].itfcI[i];
      cj = user[bi].itfcJ[i];
      ck = user[bi].itfcK[i];

      hi = user[bi].itfchostI[i];
      hj = user[bi].itfchostJ[i];
      hk = user[bi].itfchostK[i];
      hb = user[bi].itfchostB[i];

      x = user[bi].itfchostx[i];
      y = user[bi].itfchosty[i];
      z = user[bi].itfchostz[i];


      if (ci>=xs && ci<xe &&
	  cj>=ys && cj<ye &&
	  ck>=zs && ck<ze) {

	VecGetArray(user[hb].nhostU, &hostu);
	lmx = user[hb].info.mx; lmy = user[hb].info.my; lmz = user[hb].info.mz;
	itfc[ck][cj][ci].x = (hostu[((hk  )*lmx*lmy + (hj  )*lmx + (hi  )) * 3]
			      * (1-x) * (1-y) * (1-z) + //i,j,k
			      hostu[((hk  )*lmx*lmy + (hj  )*lmx + (hi+1)) * 3]
			      * x     * (1-y) * (1-z) + //i+1,j,k
 			      hostu[((hk  )*lmx*lmy + (hj+1)*lmx + (hi  )) * 3]
			      * (1-x) * y     * (1-z) + //i,j+1,k
 			      hostu[((hk  )*lmx*lmy + (hj+1)*lmx + (hi+1)) * 3]
			      * x     * y     * (1-z) + //i+1,j+1,k
 			      hostu[((hk+1)*lmx*lmy + (hj  )*lmx + (hi  )) * 3]
			      * (1-x) * (1-y) * z     + //i,j,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj  )*lmx + (hi+1)) * 3]
			      * x     * (1-y) * z     + //i+1,j,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj+1)*lmx + (hi  )) * 3]
			      * (1-x) * y     * z     + //i,j+1,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj+1)*lmx + (hi+1)) * 3]
			      * x     * y     * z); //i+1,j+1,k+1

	itfc[ck][cj][ci].y = (hostu[((hk  )*lmx*lmy + (hj  )*lmx + (hi  ))*3+1]
			      * (1-x) * (1-y) * (1-z) + //i,j,k
			      hostu[((hk  )*lmx*lmy + (hj  )*lmx + (hi+1))*3+1]
			      * x     * (1-y) * (1-z) + //i+1,j,k
 			      hostu[((hk  )*lmx*lmy + (hj+1)*lmx + (hi  ))*3+1]
			      * (1-x) * y     * (1-z) + //i,j+1,k
 			      hostu[((hk  )*lmx*lmy + (hj+1)*lmx + (hi+1))*3+1]
			      * x     * y     * (1-z) + //i+1,j+1,k
 			      hostu[((hk+1)*lmx*lmy + (hj  )*lmx + (hi  ))*3+1]
			      * (1-x) * (1-y) * z     + //i,j,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj  )*lmx + (hi+1))*3+1]
			      * x     * (1-y) * z     + //i+1,j,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj+1)*lmx + (hi  ))*3+1]
			      * (1-x) * y     * z     + //i,j+1,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj+1)*lmx + (hi+1))*3+1]
			      * x     * y     * z); //i+1,j+1,k+1

	itfc[ck][cj][ci].z = (hostu[((hk  )*lmx*lmy + (hj  )*lmx + (hi  ))*3+2]
			      * (1-x) * (1-y) * (1-z) + //i,j,k
			      hostu[((hk  )*lmx*lmy + (hj  )*lmx + (hi+1))*3+2]
			      * x     * (1-y) * (1-z) + //i+1,j,k
 			      hostu[((hk  )*lmx*lmy + (hj+1)*lmx + (hi  ))*3+2]
			      * (1-x) * y     * (1-z) + //i,j+1,k
 			      hostu[((hk  )*lmx*lmy + (hj+1)*lmx + (hi+1))*3+2]
			      * x     * y     * (1-z) + //i+1,j+1,k
 			      hostu[((hk+1)*lmx*lmy + (hj  )*lmx + (hi  ))*3+2]
			      * (1-x) * (1-y) * z     + //i,j,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj  )*lmx + (hi+1))*3+2]
			      * x     * (1-y) * z     + //i+1,j,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj+1)*lmx + (hi  ))*3+2]
			      * (1-x) * y     * z     + //i,j+1,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj+1)*lmx + (hi+1))*3+2]
			      * x     * y     * z); //i+1,j+1,k+1

	VecRestoreArray(user[hb].nhostU, &hostu);
      }
      // Is the point a local point?
      // Get the host cell information from the host CPU
      // Update
    }
    PetscBarrier(PETSC_NULL);
    DMDAVecRestoreArray(user[bi].fda, user[bi].lItfc, &itfc);
    DMDALocalToLocalBegin(user[bi].fda, user[bi].lItfc, INSERT_VALUES,
			user[bi].lItfc);
    DMDALocalToLocalEnd(user[bi].fda, user[bi].lItfc, INSERT_VALUES,
		      user[bi].lItfc);
  }
  VecDestroy(&hostU);
  
  for (bi=0; bi<block_number; bi++) {

    DMDALocalInfo	info = user[bi].info;
    PetscInt	xs = info.xs, xe = info.xs + info.xm;
    PetscInt  	ys = info.ys, ye = info.ys + info.ym;
    PetscInt	zs = info.zs, ze = info.zs + info.zm;
    PetscInt	mx = info.mx, my = info.my, mz = info.mz;

    Cmpnts ***ucont, ***kzet, ***jeta;

    DMDAVecGetArray(user[bi].fda, user[bi].Ucat, &ucat);
    DMDAVecGetArray(user[bi].fda, user[bi].lItfc, &itfc);
    DMDAVecGetArray(user[bi].fda, user[bi].Ucont, &ucont);
    DMDAVecGetArray(user[bi].fda, user[bi].lKZet, &kzet);
    DMDAVecGetArray(user[bi].fda, user[bi].lJEta, &jeta);
    if (user[bi].bctype[0] == 0 && xs==0) {
      i=1;
      for (k=zs; k<ze; k++) {
	for (j=ys; j<ye; j++) {
	  if (k>0 && j>0) {
	    ucat[k][j][i].x = 0.125 * (itfc[k  ][j  ][i  ].x +
				       itfc[k  ][j  ][i-1].x +
				       itfc[k  ][j-1][i  ].x +
				       itfc[k-1][j  ][i  ].x +
				       itfc[k  ][j-1][i-1].x +
				       itfc[k-1][j  ][i-1].x +
				       itfc[k-1][j-1][i  ].x +
				       itfc[k-1][j-1][i-1].x);
	
	    ucat[k][j][i].y = 0.125 * (itfc[k  ][j  ][i  ].y +
				       itfc[k  ][j  ][i-1].y +
				       itfc[k  ][j-1][i  ].y +
				       itfc[k-1][j  ][i  ].y +
				       itfc[k  ][j-1][i-1].y +
				       itfc[k-1][j  ][i-1].y +
				       itfc[k-1][j-1][i  ].y +
				       itfc[k-1][j-1][i-1].y);
	
	    ucat[k][j][i].z = 0.125 * (itfc[k  ][j  ][i  ].z +
				       itfc[k  ][j  ][i-1].z +
				       itfc[k  ][j-1][i  ].z +
				       itfc[k-1][j  ][i  ].z +
				       itfc[k  ][j-1][i-1].z +
				       itfc[k-1][j  ][i-1].z +
				       itfc[k-1][j-1][i  ].z +
				       itfc[k-1][j-1][i-1].z);
	  }
	}
      }
    }

    if (user[bi].bctype[1] == 0 && xe==mx) {
      i=mx-2;
      for (k=zs; k<ze; k++) {
	for (j=ys; j<ye; j++) {
	  if (k>0 && j>0) {
	    ucat[k][j][i].x = 0.125 * (itfc[k  ][j  ][i  ].x +
				       itfc[k  ][j  ][i-1].x +
				       itfc[k  ][j-1][i  ].x +
				       itfc[k-1][j  ][i  ].x +
				       itfc[k  ][j-1][i-1].x +
				       itfc[k-1][j  ][i-1].x +
				       itfc[k-1][j-1][i  ].x +
				       itfc[k-1][j-1][i-1].x);
	
	    ucat[k][j][i].y = 0.125 * (itfc[k  ][j  ][i  ].y +
				       itfc[k  ][j  ][i-1].y +
				       itfc[k  ][j-1][i  ].y +
				       itfc[k-1][j  ][i  ].y +
				       itfc[k  ][j-1][i-1].y +
				       itfc[k-1][j  ][i-1].y +
				       itfc[k-1][j-1][i  ].y +
				       itfc[k-1][j-1][i-1].y);
	
	    ucat[k][j][i].z = 0.125 * (itfc[k  ][j  ][i  ].z +
				       itfc[k  ][j  ][i-1].z +
				       itfc[k  ][j-1][i  ].z +
				       itfc[k-1][j  ][i  ].z +
				       itfc[k  ][j-1][i-1].z +
				       itfc[k-1][j  ][i-1].z +
				       itfc[k-1][j-1][i  ].z +
				       itfc[k-1][j-1][i-1].z);
	  }
	}
      }
    }

    if (user[bi].bctype[2] == 0 && ys==0) {
      j=1;
      for (k=zs; k<ze; k++) {
	for (i=xs; i<xe; i++) {
	  if (k>0 && i>0) {
	    ucat[k][j][i].x = 0.125 * (itfc[k  ][j  ][i  ].x +
				       itfc[k  ][j  ][i-1].x +
				       itfc[k  ][j-1][i  ].x +
				       itfc[k-1][j  ][i  ].x +
				       itfc[k  ][j-1][i-1].x +
				       itfc[k-1][j  ][i-1].x +
				       itfc[k-1][j-1][i  ].x +
				       itfc[k-1][j-1][i-1].x);
	
	    ucat[k][j][i].y = 0.125 * (itfc[k  ][j  ][i  ].y +
				       itfc[k  ][j  ][i-1].y +
				       itfc[k  ][j-1][i  ].y +
				       itfc[k-1][j  ][i  ].y +
				       itfc[k  ][j-1][i-1].y +
				       itfc[k-1][j  ][i-1].y +
				       itfc[k-1][j-1][i  ].y +
				       itfc[k-1][j-1][i-1].y);
	
	    ucat[k][j][i].z = 0.125 * (itfc[k  ][j  ][i  ].z +
				       itfc[k  ][j  ][i-1].z +
				       itfc[k  ][j-1][i  ].z +
				       itfc[k-1][j  ][i  ].z +
				       itfc[k  ][j-1][i-1].z +
				       itfc[k-1][j  ][i-1].z +
				       itfc[k-1][j-1][i  ].z +
				       itfc[k-1][j-1][i-1].z);
	  }
	}
      }
    }

    if (user[bi].bctype[3] == 0 && ye==my) {
      j=my-2;
      for (k=zs; k<ze; k++) {
	for (i=xs; i<xe; i++) {
	  if (k>0 && i>0) {
	    ucat[k][j][i].x = 0.125 * (itfc[k  ][j  ][i  ].x +
				       itfc[k  ][j  ][i-1].x +
				       itfc[k  ][j-1][i  ].x +
				       itfc[k-1][j  ][i  ].x +
				       itfc[k  ][j-1][i-1].x +
				       itfc[k-1][j  ][i-1].x +
				       itfc[k-1][j-1][i  ].x +
				       itfc[k-1][j-1][i-1].x);
	
	    ucat[k][j][i].y = 0.125 * (itfc[k  ][j  ][i  ].y +
				       itfc[k  ][j  ][i-1].y +
				       itfc[k  ][j-1][i  ].y +
				       itfc[k-1][j  ][i  ].y +
				       itfc[k  ][j-1][i-1].y +
				       itfc[k-1][j  ][i-1].y +
				       itfc[k-1][j-1][i  ].y +
				       itfc[k-1][j-1][i-1].y);
	
	    ucat[k][j][i].z = 0.125 * (itfc[k  ][j  ][i  ].z +
				       itfc[k  ][j  ][i-1].z +
				       itfc[k  ][j-1][i  ].z +
				       itfc[k-1][j  ][i  ].z +
				       itfc[k  ][j-1][i-1].z +
				       itfc[k-1][j  ][i-1].z +
				       itfc[k-1][j-1][i  ].z +
				       itfc[k-1][j-1][i-1].z);
	  }
	}
      }
    }




    if (user[bi].bctype[4] == 0 && zs==0) {
      k=0;
      for (j=ys; j<ye; j++) {
	for (i=xs; i<xe; i++) {
	  if (j>0 && i>0) {
	    ucont[k][j][i].z = (0.25 * (itfc[k  ][j  ][i  ].x +
					itfc[k  ][j  ][i-1].x +
					itfc[k  ][j-1][i  ].x +
					itfc[k  ][j-1][i-1].x) *
				kzet[k][j][i].x +
				0.25 * (itfc[k  ][j  ][i  ].y +
					itfc[k  ][j  ][i-1].y +
					itfc[k  ][j-1][i  ].y +
					itfc[k  ][j-1][i-1].y) *
				kzet[k][j][i].y +
				0.25 * (itfc[k  ][j  ][i  ].z +
					itfc[k  ][j  ][i-1].z +
					itfc[k  ][j-1][i  ].z +
					itfc[k  ][j-1][i-1].z) *
				kzet[k][j][i].z);
				    
/* 	    ucat[k][j][i].x = 0.125 * (itfc[k  ][j  ][i  ].x + */
/* 				       itfc[k  ][j  ][i-1].x + */
/* 				       itfc[k  ][j-1][i  ].x + */
/* 				       itfc[k-1][j  ][i  ].x + */
/* 				       itfc[k  ][j-1][i-1].x + */
/* 				       itfc[k-1][j  ][i-1].x + */
/* 				       itfc[k-1][j-1][i  ].x + */
/* 				       itfc[k-1][j-1][i-1].x); */
	
/* 	    ucat[k][j][i].y = 0.125 * (itfc[k  ][j  ][i  ].y + */
/* 				       itfc[k  ][j  ][i-1].y + */
/* 				       itfc[k  ][j-1][i  ].y + */
/* 				       itfc[k-1][j  ][i  ].y + */
/* 				       itfc[k  ][j-1][i-1].y + */
/* 				       itfc[k-1][j  ][i-1].y + */
/* 				       itfc[k-1][j-1][i  ].y + */
/* 				       itfc[k-1][j-1][i-1].y); */
	
/* 	    ucat[k][j][i].z = 0.125 * (itfc[k  ][j  ][i  ].z + */
/* 				       itfc[k  ][j  ][i-1].z + */
/* 				       itfc[k  ][j-1][i  ].z + */
/* 				       itfc[k-1][j  ][i  ].z + */
/* 				       itfc[k  ][j-1][i-1].z + */
/* 				       itfc[k-1][j  ][i-1].z + */
/* 				       itfc[k-1][j-1][i  ].z + */
/* 				       itfc[k-1][j-1][i-1].z); */

/* 	    if (i==1 || j==1 || i==mx-2 || j==my-2) { */
/* 	      ucat[k][j][i].x *= 2.; */
/* 	      ucat[k][j][i].y *= 2.; */
/* 	      ucat[k][j][i].z *= 2.; */
/* 	    } */

/* 	    ucx = (ucat[k][j][i].x + ucat[k+1][j][i].x) * 0.5; */
/* 	    ucy = (ucat[k][j][i].y + ucat[k+1][j][i].y) * 0.5; */
/* 	    ucz = (ucat[k][j][i].z + ucat[k+1][j][i].z) * 0.5; */

/* 	    ucont[k][j][i].z = ucx * kzet[k][j][i].x + */
/* 	      ucy * kzet[k][j][i].y + ucz * kzet[k][j][i].z; */
	  }
	}
      }

    }

    if (user[bi].bctype[5] == 0 && ze==mz) {
      k=mz-2;
      for (j=ys; j<ye; j++) {
	for (i=xs; i<xe; i++) {
	  if (j>0 && i>0) {
	    ucont[k][j][i].z = (0.25 * (itfc[k  ][j  ][i  ].x +
					itfc[k  ][j  ][i-1].x +
					itfc[k  ][j-1][i  ].x +
					itfc[k  ][j-1][i-1].x) *
				kzet[k][j][i].x +
				0.25 * (itfc[k  ][j  ][i  ].y +
					itfc[k  ][j  ][i-1].y +
					itfc[k  ][j-1][i  ].y +
					itfc[k  ][j-1][i-1].y) *
				kzet[k][j][i].y +
				0.25 * (itfc[k  ][j  ][i  ].z +
					itfc[k  ][j  ][i-1].z +
					itfc[k  ][j-1][i  ].z +
					itfc[k  ][j-1][i-1].z) *
				kzet[k][j][i].z);
/* 	    ucat[k][j][i].x = 0.125 * (itfc[k  ][j  ][i  ].x + */
/* 				       itfc[k  ][j  ][i-1].x + */
/* 				       itfc[k  ][j-1][i  ].x + */
/* 				       itfc[k-1][j  ][i  ].x + */
/* 				       itfc[k  ][j-1][i-1].x + */
/* 				       itfc[k-1][j  ][i-1].x + */
/* 				       itfc[k-1][j-1][i  ].x + */
/* 				       itfc[k-1][j-1][i-1].x); */
	
/* 	    ucat[k][j][i].y = 0.125 * (itfc[k  ][j  ][i  ].y + */
/* 				       itfc[k  ][j  ][i-1].y + */
/* 				       itfc[k  ][j-1][i  ].y + */
/* 				       itfc[k-1][j  ][i  ].y + */
/* 				       itfc[k  ][j-1][i-1].y + */
/* 				       itfc[k-1][j  ][i-1].y + */
/* 				       itfc[k-1][j-1][i  ].y + */
/* 				       itfc[k-1][j-1][i-1].y); */
	
/* 	    ucat[k][j][i].z = 0.125 * (itfc[k  ][j  ][i  ].z + */
/* 				       itfc[k  ][j  ][i-1].z + */
/* 				       itfc[k  ][j-1][i  ].z + */
/* 				       itfc[k-1][j  ][i  ].z + */
/* 				       itfc[k  ][j-1][i-1].z + */
/* 				       itfc[k-1][j  ][i-1].z + */
/* 				       itfc[k-1][j-1][i  ].z + */
/* 				       itfc[k-1][j-1][i-1].z); */

/* 	    if (i==1 || j==1 || i==mx-2 || j==my-2) { */
/* 	      ucat[k][j][i].x *= 2.; */
/* 	      ucat[k][j][i].y *= 2.; */
/* 	      ucat[k][j][i].z *= 2.; */
/* 	    } */


/* 	    ucx = (ucat[k][j][i].x + ucat[k-1][j][i].x) * 0.5; */
/* 	    ucy = (ucat[k][j][i].y + ucat[k-1][j][i].y) * 0.5; */
/* 	    ucz = (ucat[k][j][i].z + ucat[k-1][j][i].z) * 0.5; */

/* 	    ucont[k-1][j][i].z = ucx * kzet[k-1][j][i].x + */
/* 	      ucy * kzet[k-1][j][i].y + ucz * kzet[k-1][j][i].z; */
	  }
	}
      }
    }

    if (user[bi].bctype[2] == 0 && ys==0) {
      j=0;
      for (k=zs; k<ze; k++) {
	for (i=xs; i<xe; i++) {
	  if (k>0 && i>0) {
	
	    ucont[k][j][i].y = (0.25 * (itfc[k-1][j  ][i  ].x +
					itfc[k  ][j  ][i-1].x +
					itfc[k  ][j  ][i  ].x +
					itfc[k-1][j  ][i-1].x) *
				jeta[k][j][i].x +
				0.25 * (itfc[k-1][j  ][i  ].y +
					itfc[k  ][j  ][i-1].y +
					itfc[k  ][j  ][i  ].y +
					itfc[k-1][j  ][i-1].y) *
				jeta[k][j][i].y +
				0.25 * (itfc[k-1][j  ][i  ].z +
					itfc[k  ][j  ][i-1].z +
					itfc[k  ][j  ][i  ].z +
					itfc[k-1][j  ][i-1].z) *
				jeta[k][j][i].z);
				
	
	  }
	}
      }
    }

    if (user[bi].bctype[3] == 0 && ye==my) {
      j=my-2;
      for (k=zs; k<ze; k++) {
	for (i=xs; i<xe; i++) {
	  if (k>0 && i>0) {
	
	    ucont[k][j][i].y = (0.25 * (itfc[k-1][j  ][i  ].x +
					itfc[k  ][j  ][i-1].x +
					itfc[k  ][j  ][i  ].x +
					itfc[k-1][j  ][i-1].x) *
				jeta[k][j][i].x +
				0.25 * (itfc[k-1][j  ][i  ].y +
					itfc[k  ][j  ][i-1].y +
					itfc[k  ][j  ][i  ].y +
					itfc[k-1][j  ][i-1].y) *
				jeta[k][j][i].y +
				0.25 * (itfc[k-1][j  ][i  ].z +
					itfc[k  ][j  ][i-1].z +
					itfc[k  ][j  ][i  ].z +
					itfc[k-1][j  ][i-1].z) *
				jeta[k][j][i].z);

	
	  }
	}
      }
    }



    DMDAVecRestoreArray(user[bi].fda, user[bi].lItfc, &itfc);
    DMDAVecRestoreArray(user[bi].fda, user[bi].Ucat, &ucat);
    DMDAVecRestoreArray(user[bi].fda, user[bi].Ucont, &ucont);
    DMDAVecRestoreArray(user[bi].fda, user[bi].lKZet, &kzet);
    DMDAVecRestoreArray(user[bi].fda, user[bi].lJEta, &jeta);
    
  }

  for (bi=0; bi<block_number; bi++) {
    VecDestroy(&user[bi].nhostU);
    Contra2Cart(&(user[bi]));
  }

  
	return 0;
}

PetscErrorCode Block_Interface_P(UserCtx *user) {
  PetscInt bi;
  PetscInt ci, cj, ck;
  PetscInt hi, hj, hk, hb;
  PetscInt i, j, k;
  PetscReal x, y, z;
  //  VecScatter tolocal;
  Vec	hostP;
  Vec	nhostP;
  PetscReal ***itfcp, *hostp, ***phi;

  VecCreateSeq(PETSC_COMM_SELF, 8, &hostP);

  
  for (bi=0; bi<block_number; bi++) {
    VecCreateSeq(PETSC_COMM_SELF,
		 user[bi].info.mx*user[bi].info.my*user[bi].info.mz,
		 &user[bi].nhostP);
/*     PetscInt N; */
/*     VecGetSize(user[bi].nhostP, &N); */
/*     PetscPrintf(PETSC_COMM_SELF, "Number %i\n", N); */
  }

  // First calculate Phi components at grid nodes
  for (bi=0; bi<block_number; bi++) {
    DMDALocalInfo	info = user[bi].info;
    PetscInt	xs = info.xs, xe = info.xs + info.xm;
    PetscInt  	ys = info.ys, ye = info.ys + info.ym;
    PetscInt	zs = info.zs, ze = info.zs + info.zm;
    PetscInt	mx = info.mx, my = info.my, mz = info.mz;
    DMDAVecGetArray(user[bi].da, user[bi].ItfcP, &itfcp);
    DMDAVecGetArray(user[bi].da, user[bi].lPhi, &phi);

    DMDACreateNaturalVector(user[bi].da, &nhostP);

    PetscBarrier(PETSC_NULL);
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	for (i=xs; i<xe; i++) {
	  if (k<mz-1 && j<my-1 && i<mx-1) {
	    itfcp[k][j][i] = 0.125 * (phi[k  ][j  ][i  ] +
				      phi[k  ][j  ][i+1] +
				      phi[k  ][j+1][i  ] +
				      phi[k  ][j+1][i+1] +
				      phi[k+1][j  ][i  ] +
				      phi[k+1][j  ][i+1] +
				      phi[k+1][j+1][i  ] +
				      phi[k+1][j+1][i+1]);
	  }
	}
      }
    }
    if (bi==0 && zs!=0) {
      PetscPrintf(PETSC_COMM_SELF, "%le PPP\n", itfcp[28][21][21]);
    }
    DMDAVecRestoreArray(user[bi].da, user[bi].ItfcP, &itfcp);
    DMDAVecRestoreArray(user[bi].da, user[bi].lPhi, &phi);

    
    PetscBarrier(PETSC_NULL);

    DMDAGlobalToNaturalBegin(user[bi].da, user[bi].ItfcP, INSERT_VALUES, nhostP);
    DMDAGlobalToNaturalEnd(user[bi].da, user[bi].ItfcP, INSERT_VALUES, nhostP);
    
/*     DMDALocalToGlobal(user[bi].da, user[bi].lItfcP, INSERT_VALUES, user[bi].ItfcP); */
    VecScatter tolocalall;
    VecScatterCreateToAll(nhostP, &tolocalall, &(user[bi].nhostP));
/*     DMDAGlobalToNaturalAllCreate(user[bi].da, &tolocalall); */
    VecScatterBegin(tolocalall, nhostP, user[bi].nhostP, INSERT_VALUES, SCATTER_FORWARD);
    
    VecScatterEnd(tolocalall, nhostP, user[bi].nhostP, INSERT_VALUES,  SCATTER_FORWARD);
    VecScatterDestroy(&tolocalall);

    VecGetArray(user[bi].nhostP, &hostp);
    if (bi==0) {
      PetscPrintf(PETSC_COMM_SELF, "%i %le PPP\n", zs, hostp[28*42*42+21*42+21]);
    }

    VecRestoreArray(user[bi].nhostP, &hostp);
    VecDestroy(&nhostP);
  }

  for (bi=0; bi<block_number; bi++) {
    DMDALocalInfo	info = user[bi].info;

    PetscInt	xs = info.xs, xe = info.xs + info.xm;
    PetscInt  	ys = info.ys, ye = info.ys + info.ym;
    PetscInt	zs = info.zs, ze = info.zs + info.zm;
    

    PetscInt    lmx, lmy, lmz;

    DMDAVecGetArray(user[bi].da, user[bi].lItfcP, &itfcp);
    for (i=0; i<user[bi].itfcptsnumber; i++) {
      ci = user[bi].itfcI[i];
      cj = user[bi].itfcJ[i];
      ck = user[bi].itfcK[i];

      hi = user[bi].itfchostI[i];
      hj = user[bi].itfchostJ[i];
      hk = user[bi].itfchostK[i];
      hb = user[bi].itfchostB[i];

      x = user[bi].itfchostx[i];
      y = user[bi].itfchosty[i];
      z = user[bi].itfchostz[i];
      
      if (ci>=xs && ci<xe &&
	  cj>=ys && cj<ye &&
	  ck>=zs && ck<ze) {
	
	VecGetArray(user[hb].nhostP, &hostp);
	lmx = user[hb].info.mx; lmy = user[hb].info.my; lmz = user[hb].info.mz;
	itfcp[ck][cj][ci] = (hostp[((hk  )*lmx*lmy + (hj  )*lmx + (hi  ))]
			     * (1-x) * (1-y) * (1-z) + //i,j,k
			     hostp[((hk  )*lmx*lmy + (hj  )*lmx + (hi+1))]
			     * x     * (1-y) * (1-z) + //i+1,j,k
			     hostp[((hk  )*lmx*lmy + (hj+1)*lmx + (hi  ))]
			     * (1-x) * y     * (1-z) + //i,j+1,k
			     hostp[((hk  )*lmx*lmy + (hj+1)*lmx + (hi+1))]
			     * x     * y     * (1-z) + //i+1,j+1,k
			     hostp[((hk+1)*lmx*lmy + (hj  )*lmx + (hi  ))]
			     * (1-x) * (1-y) * z     + //i,j,k+1
			     hostp[((hk+1)*lmx*lmy + (hj  )*lmx + (hi+1))]
			     * x     * (1-y) * z     + //i+1,j,k+1
			     hostp[((hk+1)*lmx*lmy + (hj+1)*lmx + (hi  ))]
			     * (1-x) * y     * z     + //i,j+1,k+1
			     hostp[((hk+1)*lmx*lmy + (hj+1)*lmx + (hi+1))]
			     * x     * y     * z); //i+1,j+1,k+1


	VecRestoreArray(user[hb].nhostP, &hostp);
      }
    }
    DMDAVecRestoreArray(user[bi].da, user[bi].lItfcP, &itfcp);
    DMDALocalToLocalBegin(user[bi].da, user[bi].lItfcP, INSERT_VALUES,
			user[bi].lItfcP);
    DMDALocalToLocalEnd(user[bi].da, user[bi].lItfcP, INSERT_VALUES,
		      user[bi].lItfcP);

  }
  VecDestroy(&hostP);
  
  PetscReal ***p;

  for (bi=0; bi<block_number; bi++) {
    DMDALocalInfo	info = user[bi].info;
    PetscInt	xs = info.xs, xe = info.xs + info.xm;
    PetscInt  	ys = info.ys, ye = info.ys + info.ym;
    PetscInt	zs = info.zs, ze = info.zs + info.zm;
    PetscInt	mx = info.mx, my = info.my, mz = info.mz;
    Cmpnts ***kcsi, ***keta, ***kzet, ***ucont;
    PetscReal ***kaj;
    DMDAVecGetArray(user[bi].da, user[bi].Phi, &p);
    DMDAVecGetArray(user[bi].da, user[bi].lItfcP, &itfcp);
    
    DMDAVecGetArray(user[bi].fda, user[bi].lKCsi, &kcsi);
    DMDAVecGetArray(user[bi].fda, user[bi].lKEta, &keta);
    DMDAVecGetArray(user[bi].fda, user[bi].lKZet, &kzet);
    DMDAVecGetArray(user[bi].da, user[bi].lKAj, &kaj);
    DMDAVecGetArray(user[bi].fda, user[bi].Ucont, &ucont);
    if (user[bi].bctype[0] == 0 && xs==0) {
      i=1;
      for (k=zs; k<ze; k++) {
	for (j=ys; j<ye; j++) {
	  p[k][j][i] = 0.125 * (itfcp[k  ][j  ][i  ] +
				itfcp[k  ][j  ][i-1] +
				itfcp[k  ][j-1][i  ] +
				itfcp[k-1][j  ][i  ] +
				itfcp[k  ][j-1][i-1] +
				itfcp[k-1][j  ][i-1] +
				itfcp[k-1][j-1][i  ] +
				itfcp[k-1][j-1][i-1]);
	
	}
      }
    }

    if (user[bi].bctype[1] == 0 && xe==mx) {
      i=mx-2;
      for (k=zs; k<ze; k++) {
	for (j=ys; j<ye; j++) {
	  p[k][j][i] = 0.125 * (itfcp[k  ][j  ][i  ] +
				itfcp[k  ][j  ][i-1] +
				itfcp[k  ][j-1][i  ] +
				itfcp[k-1][j  ][i  ] +
				itfcp[k  ][j-1][i-1] +
				itfcp[k-1][j  ][i-1] +
				itfcp[k-1][j-1][i  ] +
				itfcp[k-1][j-1][i-1]);
	}
      }
    }

    if (user[bi].bctype[2] == 0 && ys==0) {
      j=1;
      for (k=zs; k<ze; k++) {
	for (i=xs; i<xe; i++) {
	  p[k][j][i] = 0.125 * (itfcp[k  ][j  ][i  ] +
				itfcp[k  ][j  ][i-1] +
				itfcp[k  ][j-1][i  ] +
				itfcp[k-1][j  ][i  ] +
				itfcp[k  ][j-1][i-1] +
				itfcp[k-1][j  ][i-1] +
				itfcp[k-1][j-1][i  ] +
				itfcp[k-1][j-1][i-1]);
	}
      }
    }

    if (user[bi].bctype[3] == 0 && ye==my) {
      j=my-2;
      for (k=zs; k<ze; k++) {
	for (i=xs; i<xe; i++) {
	  p[k][j][i] = 0.125 * (itfcp[k  ][j  ][i  ] +
				itfcp[k  ][j  ][i-1] +
				itfcp[k  ][j-1][i  ] +
				itfcp[k-1][j  ][i  ] +
				itfcp[k  ][j-1][i-1] +
				itfcp[k-1][j  ][i-1] +
				itfcp[k-1][j-1][i  ] +
				itfcp[k-1][j-1][i-1]);
	}
      }
    }

    if (user[bi].bctype[4] == 0 && zs==0) {
      k=1;
      for (j=ys; j<ye; j++) {
	for (i=xs; i<xe; i++) {
	  p[k][j][i] = 0.125 * (itfcp[k  ][j  ][i  ] +
				itfcp[k  ][j  ][i-1] +
				itfcp[k  ][j-1][i  ] +
				itfcp[k-1][j  ][i  ] +
				itfcp[k  ][j-1][i-1] +
				itfcp[k-1][j  ][i-1] +
				itfcp[k-1][j-1][i  ] +
				itfcp[k-1][j-1][i-1]);
	}
      }
      if (ti>0) {
	PetscPrintf(PETSC_COMM_WORLD, "PP %le\n", p[k][21][21]);
      }
/*       k=1; */
/*       for (j=ys+1; j<ye-1; j++) { */
/* 	for (i=xs+1; i<xe-1; i++) { */
/* 	  ucont[k][j][i].z -= ((p[k+1][j][i]-p[k][j][i]) * */
/* 	    		 (kzet[k][j][i].x * kzet[k][j][i].x + */
/* 			  kzet[k][j][i].y * kzet[k][j][i].y + */
/* 			  kzet[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i]) */
/* 	    * user->dt * user->st; */

/* 	  ucont[k][j][i].z -= */
/* 	    (0.25 * (p[k][j][i+1] + p[k+1][j][i+1] - */
/* 		     p[k][j][i-1] - p[k+1][j][i-1]) * */
/* 	     (kcsi[k][j][i].x * kzet[k][j][i].x + */
/* 	      kcsi[k][j][i].y * kzet[k][j][i].y + */
/* 	      kcsi[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i]) */
/* 	    * user->dt * user->st; */
/* 	  ucont[k][j][i].z -=	       */
/* 	    (0.25 * (p[k][j+1][i] + p[k+1][j+1][i] - */
/* 		     p[k][j-1][i] - p[k+1][j-1][i]) * */
/* 	     (keta[k][j][i].x * kzet[k][j][i].x + */
/* 	      keta[k][j][i].y * kzet[k][j][i].y + */
/* 	      keta[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i]) */
/* 	    * user->dt * user->st; */
/* 	} */
/*       } */

    }

    if (user[bi].bctype[5] == 0 && ze==mz) {
      k=mz-2;
      for (j=ys; j<ye; j++) {
	for (i=xs; i<xe; i++) {
	  p[k][j][i] = 0.125 * (itfcp[k  ][j  ][i  ] +
				itfcp[k  ][j  ][i-1] +
				itfcp[k  ][j-1][i  ] +
				itfcp[k-1][j  ][i  ] +
				itfcp[k  ][j-1][i-1] +
				itfcp[k-1][j  ][i-1] +
				itfcp[k-1][j-1][i  ] +
				itfcp[k-1][j-1][i-1]);
	}
      }
      if (ti>0) {
	PetscPrintf(PETSC_COMM_SELF, "PP0 %le\n", p[k][21][21]);
      }
/*       k=mz-3; */
/*       for (j=ys+1; j<ye-1; j++) { */
/* 	for (i=xs+1; i<xe-1; i++) { */
/* 	  ucont[k][j][i].z -= ((p[k+1][j][i]-p[k][j][i]) * */
/* 	    		 (kzet[k][j][i].x * kzet[k][j][i].x + */
/* 			  kzet[k][j][i].y * kzet[k][j][i].y + */
/* 			  kzet[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i]) */
/* 	    * user->dt * user->st; */

/* 	  ucont[k][j][i].z -= */
/* 	    (0.25 * (p[k][j][i+1] + p[k+1][j][i+1] - */
/* 		     p[k][j][i-1] - p[k+1][j][i-1]) * */
/* 	     (kcsi[k][j][i].x * kzet[k][j][i].x + */
/* 	      kcsi[k][j][i].y * kzet[k][j][i].y + */
/* 	      kcsi[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i]) */
/* 	    * user->dt * user->st; */
/* 	  ucont[k][j][i].z -=	       */
/* 	    (0.25 * (p[k][j+1][i] + p[k+1][j+1][i] - */
/* 		     p[k][j-1][i] - p[k+1][j-1][i]) * */
/* 	     (keta[k][j][i].x * kzet[k][j][i].x + */
/* 	      keta[k][j][i].y * kzet[k][j][i].y + */
/* 	      keta[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i]) */
/* 	    * user->dt * user->st; */
/* 	} */
/*       } */
    }
    DMDAVecRestoreArray(user[bi].da, user[bi].lItfcP, &itfcp);
    DMDAVecRestoreArray(user[bi].da, user[bi].Phi, &p);
    DMDAVecRestoreArray(user[bi].fda, user[bi].lKCsi, &kcsi);
    DMDAVecRestoreArray(user[bi].fda, user[bi].lKEta, &keta);
    DMDAVecRestoreArray(user[bi].fda, user[bi].lKZet, &kzet);
    DMDAVecRestoreArray(user[bi].da, user[bi].lKAj, &kaj);
    DMDAVecRestoreArray(user[bi].fda, user[bi].Ucont, &ucont);
  }

  for (bi=0; bi<block_number; bi++) {
    VecDestroy(&user[bi].nhostP);
  }
	return 0;
}


/* Boundary condition defination (array user->bctype[0-5]):
   0:	interpolation
   1:	solid wall (not moving)
   2:	moving solid wall (U=1)
   5:	Inlet
   4:	Outlet
   8:   Characteristic BC
*/
int outflow_scale=1;

PetscErrorCode FormBCS(UserCtx *user, FSInfo *fsi)
{
  DM da = user->da, fda = user->fda;
  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;
  PetscInt	i, j, k;

  Vec		Coor;
  Cmpnts	***ucont, ***ubcs, ***ucat, ***coor, ***csi, ***eta, ***zet;
  PetscScalar	FluxIn, FluxOut, ratio;
  PetscScalar   lArea, AreaSum, ***level;
  PetscScalar   FarFluxIn=0., FarFluxOut=0., FarFluxInSum, FarFluxOutSum;
  PetscScalar   FarAreaIn=0., FarAreaOut=0., FarAreaInSum, FarAreaOutSum;
  PetscScalar   FluxDiff, VelDiffIn, VelDiffOut;
  Cmpnts        V_frame;
  PetscInt      moveframe=1;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

	double ibm_Flux=0, ibm_Area=0;
	//CalcVolumeFlux(user, user->Ucont, &ibm_Flux, &ibm_Area);
	ibm_Flux=0, ibm_Area=0;
    
  DMDAGetGhostedCoordinates(da, &Coor);
  DMDAVecGetArray(fda, Coor, &coor);
/*   DMDAVecGetArray(fda, user->Ucont, &ucont); */
  DMDAVecGetArray(fda, user->Bcs.Ubcs, &ubcs);
	if(levelset) DMDAVecGetArray(da, user->lLevelset,  &level);

  DMDAVecGetArray(fda, user->lCsi,  &csi);
  DMDAVecGetArray(fda, user->lEta,  &eta);
  DMDAVecGetArray(fda, user->lZet,  &zet);


  //PetscInt ttemp;
  /*for (ttemp=0; ttemp<5; ttemp++) */{	// 5 times? why. seokkoo
  Contra2Cart(user);
  DMDAVecGetArray(fda, user->Ucat,  &ucat);

/* ==================================================================================             */
/*   FAR-FIELD BC */
/* ==================================================================================             */
  DMDAVecGetArray(fda, user->lUcont, &ucont);


  if (user->bctype[5]==6) {
    if (moveframe) {
      V_frame.x= -(fsi->S_new[1]-fsi->S_old[1]);
      V_frame.y= -(fsi->S_new[3]-fsi->S_old[3]);
      V_frame.z= -(fsi->S_new[5]-fsi->S_old[5]);
    } else {
      V_frame.x=0.;
      V_frame.y=0.;
      V_frame.z=0.;
    }

    if (moveframe) {
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
	  for (i=lxs; i<lxe; i++) {
	    if (k>1 && k<ze-1)
	    ucont[k-1][j][i].z += V_frame.z * 0.5*(zet[k-1][j  ][i  ].z + zet[k][j][i].z) +
	                          V_frame.y * 0.5*(zet[k-1][j  ][i  ].y + zet[k][j][i].y) +
	                          V_frame.x * 0.5*(zet[k-1][j  ][i  ].x + zet[k][j][i].x)  ;
	    if (j>1 && j<ye-1)
	    ucont[k][j-1][i].y += V_frame.z * eta[k  ][j-1][i  ].z +
	                          V_frame.y * eta[k  ][j-1][i  ].y +
	                          V_frame.x * eta[k  ][j-1][i  ].x;
	    if (i>1 && i<xe-1)
	    ucont[k][j][i-1].x += V_frame.z * csi[k  ][j  ][i-1].z +
	                          V_frame.y * csi[k  ][j  ][i-1].y +
	                          V_frame.x * csi[k  ][j  ][i-1].x;
	  }
	}
      }
      
    }
        
  }

  DMDAVecRestoreArray(fda, user->lUcont, &ucont);

  DMDAVecGetArray(fda, user->Ucont, &ucont);

  if (user->bctype[0]==6) {
    if (xs == 0) {
      i= xs;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
	  ubcs[k][j][i].x = ucat[k][j][i+1].x;
	  ubcs[k][j][i].y = ucat[k][j][i+1].y;
	  ubcs[k][j][i].z = ucat[k][j][i+1].z;	
	  ucont[k][j][i].x = ubcs[k][j][i].x * csi[k][j][i+1].x;
	  FarFluxIn += ucont[k][j][i].x;
	  FarAreaIn += csi[k][j][i].x;
	}
      }
    }
  }
  
  if (user->bctype[1]==6) {
    if (xe==mx) {
      i= xe-1;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
	  ubcs[k][j][i].x = ucat[k][j][i-1].x;
	  ubcs[k][j][i].y = ucat[k][j][i-1].y;
	  ubcs[k][j][i].z = ucat[k][j][i-1].z;
	  ucont[k][j][i-1].x = ubcs[k][j][i].x * csi[k][j][i-1].x;
/* 	  FarFluxIn -= ucont[k][j][i-1].x; */
	  FarFluxOut += ucont[k][j][i-1].x;
	  FarAreaOut += csi[k][j][i-1].x;
	}
      }
    }
  }

  if (user->bctype[2]==6) {
    if (ys==0) {
      j= ys;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = ucat[k][j+1][i].x;
	  ubcs[k][j][i].y = ucat[k][j+1][i].y;
	  ubcs[k][j][i].z = ucat[k][j+1][i].z;
	  ucont[k][j][i].y = ubcs[k][j][i].y * eta[k][j+1][i].y;
	  FarFluxIn += ucont[k][j][i].y;
	  FarAreaIn += eta[k][j][i].y;
	}
      }
    }
  }
  
  if (user->bctype[3]==6) {
    if (ye==my) {
      j=ye-1;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = ucat[k][j-1][i].x;
	  ubcs[k][j][i].y = ucat[k][j-1][i].y;
	  ubcs[k][j][i].z = ucat[k][j-1][i].z;
	  ucont[k][j-1][i].y = ubcs[k][j][i].y * eta[k][j-1][i].y;
/* 	  FarFluxIn -= ucont[k][j-1][i].y; */
	  FarFluxOut += ucont[k][j-1][i].y;
	  FarAreaOut += eta[k][j-1][i].y;
	}
      }
    }
  }

  if (user->bctype[4]==6) {
    if (zs==0) {
      k = 0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {  
	  ubcs[k][j][i].x = ucat[k+1][j][i].x;
	  ubcs[k][j][i].y = ucat[k+1][j][i].y;
	  ubcs[k][j][i].z = ucat[k+1][j][i].z;
	  ucont[k][j][i].z = ubcs[k][j][i].z * zet[k+1][j][i].z;
	  FarFluxIn += ucont[k][j][i].z;
	  FarAreaIn += zet[k][j][i].z;
	}
      }
    }
  }

  if (user->bctype[5]==6) {
    if (ze==mz) {
      k = ze-1;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {  
	  ubcs[k][j][i].x = ucat[k-1][j][i].x;
	  ubcs[k][j][i].y = ucat[k-1][j][i].y;
	  ubcs[k][j][i].z = ucat[k-1][j][i].z;
	  ucont[k-1][j][i].z = ubcs[k][j][i].z * zet[k-1][j][i].z;
/* 	  FarFluxIn -= ucont[k-1][j][i].z; */
	  FarFluxOut += ucont[k-1][j][i].z;
	  FarAreaOut += zet[k-1][j][i].z;
	}
      }
    }
  }

  GlobalSum_All(&FarFluxIn, &FarFluxInSum, PETSC_COMM_WORLD);
  GlobalSum_All(&FarFluxOut, &FarFluxOutSum, PETSC_COMM_WORLD);

  GlobalSum_All(&FarAreaIn, &FarAreaInSum, PETSC_COMM_WORLD);
  GlobalSum_All(&FarAreaOut, &FarAreaOutSum, PETSC_COMM_WORLD);

  if (user->bctype[5]==6) {
    FluxDiff = 0.5*(FarFluxInSum - FarFluxOutSum) ;
    VelDiffIn  = FluxDiff / FarAreaInSum ;
    if (fabs(FluxDiff) < 1.e-6) VelDiffIn = 0.;
    if (fabs(FarAreaInSum) <1.e-6) VelDiffIn = 0.;

    VelDiffOut  = FluxDiff / FarAreaOutSum ;
    if (fabs(FluxDiff) < 1.e-6) VelDiffOut = 0.;
    if (fabs(FarAreaOutSum) <1.e-6) VelDiffOut = 0.;
    if (moveframe) {
      V_frame.x= -(fsi->S_new[1]-fsi->S_old[1]);
      V_frame.y= -(fsi->S_new[3]-fsi->S_old[3]);
      V_frame.z= -(fsi->S_new[5]-fsi->S_old[5]);
    } else {
      V_frame.x=0.;
      V_frame.y=0.;
      V_frame.z=0.;
    }

    PetscPrintf(PETSC_COMM_WORLD, "Far Flux Diff %d %le %le %le %le %le %le %le\n", ti, FarFluxInSum, FarFluxOutSum, FluxDiff, FarAreaInSum, FarAreaOutSum, VelDiffIn, VelDiffOut);
    PetscPrintf(PETSC_COMM_WORLD, "Cop Vel  Diff %d %le %le %le \n", ti, V_frame.x,V_frame.y,V_frame.z);    
        
  }


  // scale global mass conservation

  if (user->bctype[5]==6) {
    if (ze==mz) {
      k = ze-1;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].z = ucat[k-1][j][i].z + VelDiffOut + V_frame.z;
	  ucont[k-1][j][i].z = ubcs[k][j][i].z * zet[k-1][j][i].z;
	}
      }
    }
  }

  if (user->bctype[3]==6) {
    if (ye==my) {
      j=ye-1;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].y = ucat[k][j-1][i].y + VelDiffOut + V_frame.y;
	  ucont[k][j-1][i].y = ubcs[k][j][i].y * eta[k][j-1][i].y;
	}
      }
    }
  }
    
  if (user->bctype[1]==6) {
    if (xe==mx) {
      i= xe-1;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
	  ubcs[k][j][i].x = ucat[k][j][i-1].x + VelDiffOut + V_frame.x;
	  ucont[k][j][i-1].x = ubcs[k][j][i].x * csi[k][j][i-1].x;
	}
      }
    }
  }


  if (user->bctype[0]==6) {
    if (xs == 0) {
      i= xs;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
	  ubcs[k][j][i].x = ucat[k][j][i+1].x - VelDiffIn + V_frame.x;
	  ucont[k][j][i].x = ubcs[k][j][i].x * csi[k][j][i+1].x;
	}
      }
    }
  }
  

  if (user->bctype[2]==6) {
    if (ys==0) {
      j= ys;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].y = ucat[k][j+1][i].y - VelDiffIn + V_frame.y;
	  ucont[k][j][i].y = ubcs[k][j][i].y * eta[k][j+1][i].y;
	}
      }
    }
  }
  

  if (user->bctype[4]==6) {
    if (zs==0) {
      k = 0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].z = ucat[k+1][j][i].z - VelDiffIn + V_frame.z;
	  ucont[k][j][i].z = ubcs[k][j][i].z * zet[k+1][j][i].z;
	}
      }
    }
  }


/* ==================================================================================             */
/*     CHARACTERISTIC OUTLET BC :8 */
/* ==================================================================================             */

  if (user->bctype[5]==8) {
    if (ze == mz) {
      k = ze-2;
      FluxOut = 0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  FluxOut += ucont[k][j][i].z;
	}
      }
    }
    else {
      FluxOut = 0.;
    }
    
    FluxIn = FluxInSum + FarFluxInSum;
    GlobalSum_All(&FluxOut, &FluxOutSum, PETSC_COMM_WORLD);

    //ratio = FluxInSum / FluxOutSum;
    ratio = FluxIn / FluxOutSum;
    if (fabs(FluxOutSum) < 1.e-6) ratio = 1.;
    //if (fabs(FluxInSum) <1.e-6) ratio = 0.;
    if (fabs(FluxIn) <1.e-6) ratio = 0.;
    PetscPrintf(PETSC_COMM_WORLD, "Char Ratio %d %le %le %le %le %d %d\n", ti, ratio, FluxIn, FluxOutSum, FarFluxInSum,zs, ze);

    if (ze==mz) {
      k = ze-1;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {  
	  ubcs[k][j][i].x = ucat[k-1][j][i].x;
	  ubcs[k][j][i].y = ucat[k-1][j][i].y;
	  if (ti==0 || ti==1) 
	    if (inletprofile<0) 
	      ubcs[k][j][i].z = -1.;
	    else if (user->bctype[4]==6) 
	      ubcs[k][j][i].z = 0.;
	    else
	      ubcs[k][j][i].z = 1.;//ubcs[0][j][i].z;//-1.;//1.;
	  
	  else 
	    ucont[k-1][j][i].z = ucont[k-1][j][i].z*ratio;
	  ubcs[k][j][i].z = ucont[k-1][j][i].z / zet[k-1][j][i].z;
	}
      }
    }
  }


/* ==================================================================================             */
/*     OUTLET BC :4 */
/* ==================================================================================             */
/*
// temp for plate !!
   for (k=lzs; k<lze; k++) { 
     for (j=lys; j<lye; j++) { 
       for (i=lxs; i<lxe; i++) { 
 	if (k==121 && j>96 && j<145)
 	  ucont[k][j][i].z = 0.;
       } 
     } 
   } 
*/
	double FluxOut_gas=0, lArea_gas=0, AreaSum_gas=0, ratio_gas;
  
	if(user->bctype[0]==11) {
		PetscReal	***nvert;
		DMDAVecGetArray(da, user->lNvert, &nvert);
		
		lArea=0.;
		if (xs==0) {
			i = 0;
			/*
			for (j=lys; j<lye; j++) 
			for (k=lzs; k<lze; k++) {
				double zc = ( coor[k][j][i+1].z + coor[k-1][j][i+1].z + coor[k][j-1][i+1].z + coor[k-1][j-1][i+1].z ) * 0.25;
				if( zc > 0 ) ubcs[k][j][i] = ucat[k][j][i+1];
			}
			*/
			
			FluxOut = 0;
			for (j=lys; j<lye; j++) 
			for (k=lzs; k<lze; k++) {
				double zc = ( coor[k][j][i+1].z + coor[k-1][j][i+1].z + coor[k][j-1][i+1].z + coor[k-1][j-1][i+1].z ) * 0.25;
				if( zc > 0 && nvert[k][j][i+1] < threshold) {
					
					double u=ucat[k][j][i+1].x, v=ucat[k][j][i+1].y, w=ucat[k][j][i+1].z;
					ucat[k][j][i].x=u;
					ucat[k][j][i].y=v;
					ucat[k][j][i].z=w;
					
					ucont[k][j][i].x = u*csi[k][j][i].x + v*csi[k][j][i].y + w*csi[k][j][i].z;
					
					FluxOut +=  ucont[k][j][i].x;
					//lArea += sqrt( csi[k][j][i].x*csi[k][j][i].x + csi[k][j][i].y*csi[k][j][i].y + csi[k][j][i].z*csi[k][j][i].z );
					lArea += fabs(csi[k][j][i].z);
				}
			}
		}
		else FluxOut = 0.;
		
		FluxIn = FluxInSum + FarFluxInSum;
		GlobalSum_All(&FluxOut, &FluxOutSum, PETSC_COMM_WORLD);
		GlobalSum_All(&lArea, &AreaSum, PETSC_COMM_WORLD);
		 
		//here
		FluxOutSum *= -1;
		ratio = (FluxInSum - FluxOutSum) / AreaSum;
		
		double FluxOut_new=0, FluxOut_new_sum;
		
		if(outflow_scale) {
			PetscPrintf(PETSC_COMM_WORLD, "Time %d, Vel correction=%e, FluxIn=%e, FluxOut=%e, Area=%f\n", ti, ratio, FluxInSum, FluxOutSum, AreaSum);
		
			if (xs==0) {
				i=0;
				for (j=lys; j<lye; j++) 
				for (k=lzs; k<lze; k++) {
					double zc = ( coor[k][j][i+1].z + coor[k-1][j][i+1].z + coor[k][j-1][i+1].z + coor[k-1][j-1][i+1].z ) * 0.25;
					if( zc > 0 && nvert[k][j][i+1] < threshold) {
						double Area = sqrt( csi[k][j][i+1].x*csi[k][j][i+1].x + csi[k][j][i+1].y*csi[k][j][i+1].y + csi[k][j][i+1].z*csi[k][j][i+1].z );
						Area = csi[k][j][i+1].z;
						ucont[k][j][i].x += (FluxInSum - FluxOutSum) * Area / AreaSum;
						FluxOut_new += ucont[k][j][i].x;
					}
				}
			}
		}		
		GlobalSum_All(&FluxOut_new, &FluxOut_new_sum, PETSC_COMM_WORLD);
		PetscPrintf(PETSC_COMM_WORLD, "Corrected FluxOut=%e\n", FluxOut_new_sum);
		DMDAVecRestoreArray(da, user->lNvert, &nvert);	//seokkoo 
	}
		
	if (user->bctype[5]==4 || user->bctype[5]==5) {
	  PetscReal	***nvert, ***level;	//seokkoo

		if(levelset) DMDAVecGetArray(da, user->lLevelset, &level);
		DMDAVecGetArray(da, user->lNvert, &nvert);	//seokkoo 
		
		lArea=0.;
		if (ze==mz) {
			k = ze-1;
			for (j=lys; j<lye; j++) 
			for (i=lxs; i<lxe; i++) {  
				ubcs[k][j][i].x = ucat[k-1][j][i].x;
				ubcs[k][j][i].y = ucat[k-1][j][i].y;
				if (nvert[k-1][j][i] < threshold) ubcs[k][j][i].z = ucat[k-1][j][i].z;
				else ubcs[k][j][i].z = ucat[k-1][j][i].z * 0;
			}
						
			FluxOut = 0;
			for (j=lys; j<lye; j++) 
			for (i=lxs; i<lxe; i++) {
				if (nvert[k-1][j][i] < threshold) //seokkoo
				{
					FluxOut +=  ucont[k-1][j][i].z;
					lArea += sqrt( zet[k-1][j][i].x*zet[k-1][j][i].x + zet[k-1][j][i].y*zet[k-1][j][i].y + zet[k-1][j][i].z*zet[k-1][j][i].z );
				}
			}
		}
		else	FluxOut = 0.;
		
		FluxIn = FluxInSum + FarFluxInSum;
		GlobalSum_All(&FluxOut, &FluxOutSum, PETSC_COMM_WORLD);
		GlobalSum_All(&FluxOut_gas, &FluxOutSum_gas, PETSC_COMM_WORLD);
		
		GlobalSum_All(&lArea, &AreaSum, PETSC_COMM_WORLD);
		GlobalSum_All(&lArea_gas, &AreaSum_gas, PETSC_COMM_WORLD);
		 
		ratio = (FluxInSum - FluxOutSum) / AreaSum;
		ratio_gas = (FluxInSum_gas - FluxOutSum_gas) / AreaSum_gas;
		
		if(outflow_scale) {
			PetscPrintf(PETSC_COMM_WORLD, "Time %d, Vel correction=%e, FluxIn=%e, FluxOut=%e, Area=%f\n", ti, ratio, FluxInSum, FluxOutSum, AreaSum);
		
			if (ze==mz) {
				k = ze-1;
				for (j=lys; j<lye; j++) 
				for (i=lxs; i<lxe; i++) {
					double Area = sqrt( zet[k-1][j][i].x*zet[k-1][j][i].x + zet[k-1][j][i].y*zet[k-1][j][i].y + zet[k-1][j][i].z*zet[k-1][j][i].z );
					
					if (nvert[k-1][j][i] < threshold) {
						// xiaolei deactivate levelset if 
						/*if (levelset && user->bctype[5]==4) {  
						}*/
						/*else if(levelset && user->bctype[5]==5) {
						  if(level[k-1][j][i]>0) ucont[k-1][j][i].z += (FluxInSum - FluxOutSum) * Area / AreaSum;
						  }*/
						//else 
						ucont[k-1][j][i].z += (FluxInSum - FluxOutSum) * Area / AreaSum;
					}
				}
			}
		}		
		if(levelset) DMDAVecRestoreArray(da, user->lLevelset, &level);
		DMDAVecRestoreArray(da, user->lNvert, &nvert);	//seokkoo 
	}
	
  else if (user->bctype[5]==0) {
    if (ze==mz) {
      k = ze-1;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {  
	  ubcs[k][j][i].x = ucat[k-1][j][i].x;
	  ubcs[k][j][i].y = ucat[k-1][j][i].y;
	  ubcs[k][j][i].z = ucat[k-1][j][i].z;
	}
      }
    }
  } 
  else if (user->bctype[5]==2) {
  /* Designed for driven cavity problem (top(k=kmax) wall moving)
   u_x = 1 at k==kmax */
    if (ze==mz) {
      k = ze-1;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = 1.;// - ucat[k-1][j][i].x;
	  ubcs[k][j][i].y = 0.;//- ucat[k-1][j][i].y;
	  ubcs[k][j][i].z = 0.;//- ucat[k-1][j][i].z;
	}
      }
    }
  }
	
  // slip
	Cmpnts ***lucont;	// for use of ucont[k-1] etc..
	DMDAVecGetArray(fda, user->lUcont, &lucont);
  
  
	if ( (user->bctype[0]==1 || user->bctype[0]==10) && xs==0) {
		i= 0;
		for (k=lzs; k<lze; k++) 	/* use lzs */
		for (j=lys; j<lye; j++) {
			ucont[k][j][i].x=0;
		}
	}
      
	if ( (user->bctype[1]==1 || user->bctype[1]==10) && xe==mx) {
		i= xe-1;
		for (k=lzs; k<lze; k++)
		for (j=lys; j<lye; j++) {
			ucont[k][j][i-1].x=0;
		}
	}
	
	if ( (user->bctype[2]==1 || user->bctype[2]==10) && ys==0) {
		j= 0;
		for (k=lzs; k<lze; k++) 	/* use lzs */
		for (i=lxs; i<lxe; i++) {
			ucont[k][j][i].y=0;
		}
	}
	
	if ( (user->bctype[3]==1 || user->bctype[3]==10) && ye==my) {
		j= ye-1;
		for (k=lzs; k<lze; k++) 	/* use lzs */
		for (i=lxs; i<lxe; i++) {
			ucont[k][j-1][i].y=0;
		}
	}
	/*
	if ( (user->bctype[4]==1 || user->bctype[4]==10) && zs==0) {
		k= 0;
		for (i=lxs; i<lxe; i++)
		for (j=lys; j<lye; j++) {
			ucont[k][j][i].z=0;
		}
	}
	
	if ( (user->bctype[5]==1 || user->bctype[5]==10) && ze==mz) {
		k= ze-1;
		for (i=lxs; i<lxe; i++)
		for (j=lys; j<lye; j++) {
			ucont[k-1][j][i].z=0;
		}
	}*/
    
  DMDAVecRestoreArray(fda, user->lUcont, &lucont);
//  end slip
  
  DMDAVecRestoreArray(fda, user->Ucont, &ucont);
  DMGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont);
  DMGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont);
    
  DMDAVecRestoreArray(fda, user->Ucat, &ucat);
  
  Contra2Cart(user);
  DMDAVecGetArray(fda, user->Ucat, &ucat);
  



/* ==================================================================================             */
/*   SYMMETRY BC */
/* ==================================================================================             */
  if (user->bctype[0]==3) {
	  
    if (xs==0) {
    i= xs;

    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	ubcs[k][j][i].x = 0.;
	ubcs[k][j][i].y = ucat[k][j][i+1].y;
	ubcs[k][j][i].z = ucat[k][j][i+1].z;
      }
    }
    }
  }

  if (user->bctype[1]==3) {
    if (xe==mx) {
    i= xe-1;

    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	ubcs[k][j][i].x = 0.;
	ubcs[k][j][i].y = ucat[k][j][i-1].y;
	ubcs[k][j][i].z = ucat[k][j][i-1].z;
      }
    }
    }
  }

  if (user->bctype[2]==3) {
    if (ys==0) {
    j= ys;

    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	ubcs[k][j][i].x = ucat[k][j+1][i].x;
	ubcs[k][j][i].y = 0.;
	ubcs[k][j][i].z = ucat[k][j+1][i].z;
      }
    }
    }
  }

  if (user->bctype[3]==3) {
    if (ye==my) {
    j=ye-1;

    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	ubcs[k][j][i].x = ucat[k][j-1][i].x;
	ubcs[k][j][i].y = 0.;
	ubcs[k][j][i].z = ucat[k][j-1][i].z;
      }
    }
    }
  }

/* ==================================================================================             */
/*   INTERFACE BC */
/* ==================================================================================             */
  if (user->bctype[0]==0) {
    if (xs==0) {
    i= xs;

    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	ubcs[k][j][i].x = ucat[k][j][i+1].x;
	ubcs[k][j][i].y = ucat[k][j][i+1].y;
	ubcs[k][j][i].z = ucat[k][j][i+1].z;
      }
    }
    }
  }

  if (user->bctype[1]==0) {
    if (xe==mx) {
    i= xe-1;

    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	ubcs[k][j][i].x = ucat[k][j][i-1].y;
	ubcs[k][j][i].y = ucat[k][j][i-1].y;
	ubcs[k][j][i].z = ucat[k][j][i-1].z;
      }
    }
    }
  }

  if (user->bctype[2]==0) {
    if (ys==0) {
    j= ys;

    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	ubcs[k][j][i].x = ucat[k][j+1][i].x;
	ubcs[k][j][i].y = ucat[k][j+1][i].y;
	ubcs[k][j][i].z = ucat[k][j+1][i].z;
      }
    }
    }
  }

  if (user->bctype[3]==0) {
    if (ye==my) {
    j=ye-1;

    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	ubcs[k][j][i].x = ucat[k][j-1][i].x;
	ubcs[k][j][i].y = ucat[k][j-1][i].y;
	ubcs[k][j][i].z = ucat[k][j-1][i].z;
      }
    }
    }
  }

  if (user->bctype[4]==0) {
    if (zs==0) {
      k = 0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {  
	  ubcs[k][j][i].x = ucat[k+1][j][i].x;
	  ubcs[k][j][i].y = ucat[k+1][j][i].y;
	  ubcs[k][j][i].z = ucat[k+1][j][i].z;
	}
      }
    }
  }

  if (user->bctype[5]==0) {
    if (ze==mz) {
      k = ze-1;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = ucat[k-1][j][i].x;
	  ubcs[k][j][i].y = ucat[k-1][j][i].y;
	  ubcs[k][j][i].z = ucat[k-1][j][i].z;
	}
      }
    }
  }


/* ==================================================================================             */
  // boundary conditions on ghost nodes
  
  //removed by seokkoo
	  /*
  if (xs==0) {
    i = xs;
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k][j][i+1].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k][j][i+1].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k][j][i+1].z;
      }
    }
  }

  if (xe==mx) {
    i = xe-1;
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k][j][i-1].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k][j][i-1].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k][j][i-1].z;
      }
    }
  }


  if (ys==0) {
    j = ys;
    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k][j+1][i].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k][j+1][i].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k][j+1][i].z;
      }
    }
  }

  if (ye==my) {
    j = ye-1;
    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k][j-1][i].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k][j-1][i].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k][j-1][i].z;
      }
    }
  }

  if (zs==0) {
    k = zs;
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k+1][j][i].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k+1][j][i].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k+1][j][i].z;
      }
    }
  }

  if (ze==mz) {
    k = ze-1;
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k-1][j][i].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k-1][j][i].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k-1][j][i].z;
      }
    }
  }
*/

  
  // 0 velocity on the corner point
  if (zs==0) {
    k=0;
    if (xs==0) {
      i=0;
      for (j=ys; j<ye; j++) {
	ucat[k][j][i].x = 0.;
	ucat[k][j][i].y = 0.;
	ucat[k][j][i].z = 0.;
      }
    }
    if (xe == mx) {
      i=mx-1;
      for (j=ys; j<ye; j++) {
	ucat[k][j][i].x = 0.;
	ucat[k][j][i].y = 0.;
	ucat[k][j][i].z = 0.;
      }
    }

    if (ys==0) {
      j=0;
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 0.;
	ucat[k][j][i].y = 0.;
	ucat[k][j][i].z = 0.;
      }
    }

    if (ye==my) {
      j=my-1;
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 0.;
	ucat[k][j][i].y = 0.;
	ucat[k][j][i].z = 0.;
      }
    }

  }

  if (ze==mz) {
    k=mz-1;
    if (xs==0) {
      i=0;
      for (j=ys; j<ye; j++) {
	ucat[k][j][i].x = 0.;
	ucat[k][j][i].y = 0.;
	ucat[k][j][i].z = 0.;
      }
    }
    if (xe == mx) {
      i=mx-1;
      for (j=ys; j<ye; j++) {
	ucat[k][j][i].x = 0.;
	ucat[k][j][i].y = 0.;
	ucat[k][j][i].z = 0.;
      }
    }

    if (ys==0) {
      j=0;
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 0.;
	ucat[k][j][i].y = 0.;
	ucat[k][j][i].z = 0.;
      }
    }

    if (ye==my) {
      j=my-1;
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 0.;
	ucat[k][j][i].y = 0.;
	ucat[k][j][i].z = 0.;
      }
    }

  }

  if (ys==0) {
    j=0;
    if (xs==0) {
      i=0;
      for (k=zs; k<ze; k++) {
	ucat[k][j][i].x = 0.;
	ucat[k][j][i].y = 0.;
	ucat[k][j][i].z = 0.;
      }
    }

    if (xe==mx) {
      i=mx-1;
      for (k=zs; k<ze; k++) {
	ucat[k][j][i].x = 0.;
	ucat[k][j][i].y = 0.;
	ucat[k][j][i].z = 0.;
      }
    }
  }

  if (ye==my) {
    j=my-1;
    if (xs==0) {
      i=0;
      for (k=zs; k<ze; k++) {
	ucat[k][j][i].x = 0.;
	ucat[k][j][i].y = 0.;
	ucat[k][j][i].z = 0.;
      }
    }

    if (xe==mx) {
      i=mx-1;
      for (k=zs; k<ze; k++) {
	ucat[k][j][i].x = 0.;
	ucat[k][j][i].y = 0.;
	ucat[k][j][i].z = 0.;
      }
    }
  }
  DMDAVecRestoreArray(fda, user->Ucat,  &ucat);
  }

  DMDAVecRestoreArray(fda, user->Bcs.Ubcs, &ubcs);
  DMDAVecRestoreArray(fda, Coor, &coor);

  DMDAVecRestoreArray(fda, user->lCsi,  &csi);
  DMDAVecRestoreArray(fda, user->lEta,  &eta);
  DMDAVecRestoreArray(fda, user->lZet,  &zet);
  
  if(levelset) DMDAVecRestoreArray(da, user->lLevelset,  &level);

  //  DMDAVecRestoreArray(fda, user->Ucont_o, &ucont_o);

  DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
  DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);
  return(0);
}

PetscErrorCode fluxin(UserCtx *user)
{
  PetscInt  iRotate;

  PetscInt ts_p_cycle;
  PetscInt opening, closing;
  PetscInt open_steps, close_steps;

/*   ts_p_cycle = 500; */
/*   opening = 0; */
/*   open_steps = 50; */

/*   closing = 225; */
/*   close_steps = 50; */

  //ts_p_cycle = 1290;
  ts_p_cycle = 2500;//10000;//5000;

  opening = 10;
  open_steps = 100;

  closing = 580;
  close_steps = 80;

  PetscReal t_rel;

  iRotate = ti - ((ti / ts_p_cycle) * ts_p_cycle);

  // if (angle>.0 && iRotate>1058) iRotate-=angle;

  t_rel = iRotate * (1. / ts_p_cycle) * 860 + 6.8  ; //+7.15;
/*   t_rel = (iRotate-940) * (1. / ts_p_cycle) * 860/2 + */
/*                    940. * (1. / 2500.     ) * 860 + 6.8;     */

  PetscInt i;
  PetscBool interpolated = PETSC_FALSE;
  //PetscPrintf(PETSC_COMM_WORLD, "Inflow00 Rate %d %e %e %i\n",ti, Flux_in, t_rel, user->number_flowwave);
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank) {
    for (i=0; i<user->number_flowwave-1; i++) {
/*       PetscPrintf(PETSC_COMM_WORLD, "Inflow Rate %e %e\n", Flux_in, user->inflow[i].t); */
      if (t_rel >= user->inflow[i].t && t_rel <= user->inflow[i+1].t) {
	Flux_in = user->inflow[i].f + (user->inflow[i+1].f - user->inflow[i].f) /
	  (user->inflow[i+1].t - user->inflow[i].t) *
	  (t_rel - user->inflow[i].t);
/* 	PetscPrintf(PETSC_COMM_SELF, "Inflow Rate %i %e %e %e %e %e %e\n", i, Flux_in, t_rel, user->inflow[i].f, user->inflow[i].t, user->inflow[i+1].f, user->inflow[i+1].t); */
	interpolated = PETSC_TRUE;
      }
      if (interpolated) break;
    }  
    //if (t_rel > 350) Flux_in = 0.;

    MPI_Bcast(&Flux_in, 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
  }
  else {
    MPI_Bcast(&Flux_in, 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
  }

/*   if (Flux_in<0.) { */
/*     user->bctype[5]= 5; */
/*     user->bctype[4]= 4; */
/*     PetscPrintf(PETSC_COMM_WORLD, "IINNFLOW change inlet!!!!!!!%e\n", Flux_in); */
/*   } else { */
/*     user->bctype[5]= 4; */
/*     user->bctype[4]= 5; */
/*   } */
    

/*   PetscPrintf(PETSC_COMM_SELF, "IINNFLOW %e\n", Flux_in); */
/*   Flux_in = PetscMax(0.000, sin( (iRotate) * (1./ts_p_cycle) * 2 * 3.1415926)); */

  /*
  if (iRotate >= opening && iRotate<=opening + open_steps) {
    angle = -rg + rg/(PetscReal)open_steps * (iRotate);
  }
  else if (iRotate>(closing - close_steps) && iRotate <=closing) {
    angle = -rg/(PetscReal)close_steps  * (iRotate-(closing-close_steps));
  }
  else if (iRotate>closing) {
    angle = -rg;
  }
  else {
    angle = 0;
  }
  */

/*   angle = 0.; */
  //Flux_in = 1.;
  PetscPrintf(PETSC_COMM_WORLD, "Angle %d %le %le flux-in %le intp%d\n",ti, t_rel, angle, Flux_in, interpolated);

  return 0;
}

PetscErrorCode OutflowVelocity(UserCtx *user, Vec Ucont)
{
  DM fda = user->fda;
  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  PetscReal	lFluxOut = 0., ratio;
  Cmpnts ***ucont, ***zet, ***ucat, ***ubcs;
  
  PetscInt i, j, k;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  
  if (user->bctype[5] == 4) {

    //    OutflowFlux(user);

    Contra2Cart(user);
    DMDAVecGetArray(fda, Ucont, &ucont);
    DMDAVecGetArray(fda, user->lZet, &zet);
    DMDAVecGetArray(fda, user->Ucat, &ucat);
    DMDAVecGetArray(fda, user->Bcs.Ubcs, &ubcs);
    /* Inflow flux at previous time step is 0, under this condition, it's assumed
       the flux difference at two time steps is uniformly distributed 
       to all outflow boundary nodes.*/
    if (ti==1) {//Flux_in_old < 1.e-6) {
      PetscReal lArea = 0., AreaSum;
      
      if (ze == mz) {
	k = ze-2;
	for (j=lys; j<lye; j++) {
	  for (i=lxs; i<lxe; i++) {
	    lArea += zet[k][j][i].z;
	  }
	}
      }
      GlobalSum_All(&lArea, &AreaSum, PETSC_COMM_WORLD);

      PetscReal vd;
      vd = (user->FluxInSum) / AreaSum;
      PetscPrintf(PETSC_COMM_SELF, "FluxOOOO %e %e %e\n", vd, Flux_in, Flux_in);
      
      if (ze==mz) {
	k = ze-1;
	for (j=lys; j<lye; j++) {
	  for (i=lxs; i<lxe; i++) {
	    ubcs[k][j][i].z = vd;
	    ucont[k-1][j][i].z = vd * zet[k-1][j][i].z;
	  }
	}
      }
    }
    /* Scale the outflow flux to ensure global flux conservation */
    else {
      lFluxOut = 0.;
      if (ze==mz) {
	k = ze-2;
	for (j=lys; j<lye; j++) {
	  for (i=lxs; i<lxe; i++) {
	    lFluxOut += (ucat[k][j][i].x * (zet[k][j][i].x + zet[k-1][j][i].x) +
			 ucat[k][j][i].y * (zet[k][j][i].y + zet[k-1][j][i].y) +
			 ucat[k][j][i].z * (zet[k][j][i].z + zet[k-1][j][i].z))
	      * 0.5;
	  }
	}
      }
      GlobalSum_All(&lFluxOut, &FluxOutSum, PETSC_COMM_WORLD);
      PetscBarrier(PETSC_NULL);
      ratio = user->FluxInSum / FluxOutSum;

      PetscPrintf(PETSC_COMM_WORLD, "Ratio %e %e\n", ratio, FluxOutSum);
      if (ze==mz) {
	k = ze-1;
	for (j=lys; j<lye; j++) {
	  for (i=lxs; i<lxe; i++) {
	    ubcs[k][j][i].x = ucat[k-1][j][i].x;
	    ubcs[k][j][i].y = ucat[k-1][j][i].y;
	    ubcs[k][j][i].z = ucat[k-1][j][i].z * ratio;
	    ucont[k-1][j][i].z = ubcs[k][j][i].z * zet[k-1][j][i].z;
	  }
	}
      }
      
    }
    DMDAVecRestoreArray(fda, user->Bcs.Ubcs, &ubcs);
    DMDAVecRestoreArray(fda, Ucont, &ucont);
    DMDAVecRestoreArray(fda, user->lZet, &zet);
    DMDAVecRestoreArray(fda, user->Ucat, &ucat);

/*     DMGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont); */
/*     DMGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont); */

/*     Contra2Cart(user, user->lUcont, user->Ucat); */
  }
  return 0;
}

PetscErrorCode SetInitialGuessToOne(UserCtx *user)
{
	DM da = user->da, fda = user->fda;
	DMDALocalInfo	info = user->info;
	PetscInt	xs = info.xs, xe = info.xs + info.xm;
	PetscInt  	ys = info.ys, ye = info.ys + info.ym;
	PetscInt	zs = info.zs, ze = info.zs + info.zm;
	PetscInt	mx = info.mx, my = info.my, mz = info.mz;
	PetscInt	lxs, lxe, lys, lye, lzs, lze;

	Cmpnts ***ucont, ***cent;
	Cmpnts ***icsi, ***jeta, ***kzet, ***zet;
  
	PetscInt i, j, k;

	PetscReal	***nvert, ***p, ***level;	//seokkoo
	
	
	lxs = xs; lxe = xe;
	lys = ys; lye = ye;
	lzs = zs; lze = ze;

	if (xs==0) lxs = xs+1;
	if (ys==0) lys = ys+1;
	if (zs==0) lzs = zs+1;

	if (xe==mx) lxe = xe-1;
	if (ye==my) lye = ye-1;
	if (ze==mz) lze = ze-1;
  
	Vec Coor;
	Cmpnts	***coor;
	DMDAGetGhostedCoordinates(da, &Coor);
	
	if(levelset) DMDAVecGetArray(da, user->lLevelset, &level);
	DMDAVecGetArray(da, user->lNvert, &nvert);
	DMDAVecGetArray(da, user->P, &p);
	DMDAVecGetArray(fda, Coor, &coor);
	DMDAVecGetArray(fda, user->Ucont, &ucont);
	DMDAVecGetArray(fda, user->lICsi,  &icsi);
	DMDAVecGetArray(fda, user->lJEta,  &jeta);
	DMDAVecGetArray(fda, user->lKZet,  &kzet);
	DMDAVecGetArray(fda, user->lZet,  &zet);
	DMDAVecGetArray(fda, user->lCent,  &cent);
  
	double lArea=0, SumArea;
	if(zs==0) {
		k=0;
		for (j=lys; j<lye; j++)
		for (i=lxs; i<lxe; i++) if(nvert[k][j][i]+nvert[k+1][j][i]<0.1) lArea+=sqrt( kzet[k][j][i].x*kzet[k][j][i].x + kzet[k][j][i].y*kzet[k][j][i].y + kzet[k][j][i].z*kzet[k][j][i].z );
	}
	GlobalSum_All(&lArea, &SumArea, PETSC_COMM_WORLD);
		      
	
	  
	
	double a = 2.*M_PI;
	double lambda = user->ren/2. - sqrt ( pow(user->ren/2., 2.0) + pow(a, 2.0) );
	
	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=xs; i<lxe; i++) {
		double xi = (coor[k  ][j  ][i].x + coor[k-1][j  ][i].x + coor[k  ][j-1][i].x + coor[k-1][j-1][i].x) * 0.25;	// centx[].x
		double yi = (coor[k  ][j  ][i].y + coor[k-1][j  ][i].y + coor[k  ][j-1][i].y + coor[k-1][j-1][i].y) * 0.25;	// centx[].y
		double zi = (coor[k  ][j  ][i].z + coor[k-1][j  ][i].z + coor[k  ][j-1][i].z + coor[k-1][j-1][i].z) * 0.25; 

		if(inletprofile==17) {	// 2D Taylor-Green vortex
			ucont[k][j][i].x = - cos (xi*a) * sin (yi*a) * icsi[k][j][i].x;
		}
		else if(inletprofile==18) {	// 2D Kovasznay flow
			ucont[k][j][i].x = ( 1.0 - exp ( lambda*xi )*cos( a*yi ) ) * icsi[k][j][i].x;
		}
		else if(inletprofile==20) {// Enright test
		  ucont[k][j][i].x = 2 * pow( sin (M_PI*xi), 2.) * sin (2.*M_PI*yi) * sin (2.*M_PI*zi) * icsi[k][j][i].x;
		}
		else if( inletprofile==0 ) ucont[k][j][i].x = 0;
	}
	
	for (k=lzs; k<lze; k++)
	for (j=ys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {	
		double xj = (coor[k  ][j][i  ].x + coor[k-1][j][i  ].x + coor[k  ][j][i-1].x + coor[k-1][j][i-1].x) * 0.25;
		double yj = (coor[k  ][j][i  ].y + coor[k-1][j][i  ].y + coor[k  ][j][i-1].y + coor[k-1][j][i-1].y) * 0.25;
		double zj = (coor[k  ][j][i  ].z + coor[k-1][j][i  ].z + coor[k  ][j][i-1].z + coor[k-1][j][i-1].z) * 0.25;

		if(inletprofile==17) {	// 2D Taylor-Green vortex
			ucont[k][j][i].y = sin (xj*a) * cos (yj*a) * jeta[k][j][i].y;
		}
		else if(inletprofile==18) {	// 2D Kovasznay flow
			ucont[k][j][i].y = ( lambda/a * exp ( lambda*xj )*sin( a*yj ) ) * jeta[k][j][i].y;
		}
		else if(inletprofile==20) {// Enright test
		  ucont[k][j][i].y = - sin (2.*M_PI*xj) * pow(sin (M_PI*yj), 2.) * sin (2.*M_PI*zj) * jeta[k][j][i].y;
                }
		else if( inletprofile==0 ) ucont[k][j][i].y = 0;

	}
	
	for (k=zs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
		double xk = (coor[k][j][i].x + coor[k][j-1][i].x + coor[k][j][i-1].x + coor[k][j-1][i-1].x) * 0.25;
		double yk = (coor[k][j][i].y + coor[k][j-1][i].y + coor[k][j][i-1].y + coor[k][j-1][i-1].y) * 0.25;
		double zk = (coor[k][j][i].z + coor[k][j-1][i].z + coor[k][j][i-1].z + coor[k][j-1][i-1].z) * 0.25;
		if(inletprofile==20) {// Enright test
			ucont[k][j][i].z = - sin (2.*M_PI*xk) * sin (2.*M_PI*yk) * pow( sin (M_PI*zk), 2.) * kzet[k][j][i].z;
		}
		else if( inletprofile==0 ) ucont[k][j][i].z = 0;
		else if(inletprofile==17) {	// 2D Taylor-Green vortex
			ucont[k][j][i].z = 0;
			if(k) {
				p[k][j][i] = - 0.25 * ( cos(2.*a*cent[k][j][i].x) + cos(2.*a*cent[k][j][i].y) );
			}
		}
		else if(inletprofile==18) {	// 18: 2D Kovasznay flow
			ucont[k][j][i].z = 0;
		}
		else if(inletprofile==20) {}
		else if(nvert[k][j][i]+nvert[k+1][j][i]<0.1) {
			double w;
		      
			//double area = sqrt( kzet[k][j][i].x*kzet[k][j][i].x + kzet[k][j][i].y*kzet[k][j][i].y + kzet[k][j][i].z*kzet[k][j][i].z );
			double area = sqrt( zet[k][j][i].x*zet[k][j][i].x + zet[k][j][i].y*zet[k][j][i].y + zet[k][j][i].z*zet[k][j][i].z );
			//double xc = (coor[k+1][j][i].x + coor[k+1][j-1][i].x + coor[k+1][j][i-1].x + coor[k+1][j-1][i-1].x) * 0.25;
			//double yc = (coor[k+1][j][i].y + coor[k+1][j-1][i].y + coor[k+1][j][i-1].y + coor[k+1][j-1][i-1].y) * 0.25;
		      	double xc = (coor[k][j][i].x + coor[k][j-1][i].x + coor[k][j][i-1].x + coor[k][j-1][i-1].x) * 0.25;
			double yc = (coor[k][j][i].y + coor[k][j-1][i].y + coor[k][j][i-1].y + coor[k][j-1][i-1].y) * 0.25;
		 
			if(inletprofile==10) {
				double delta = 0.45263, a=5.99;
				if( yc>=delta ) w=1.0;
				else if(yc<=0) w=0.0;
				else w = pow( yc/delta, 1./a );
			}
			else if(inletprofile==12) {	// pipe shear stress test
				double r = sqrt(xc * xc + yc * yc);
				w = 2*(1 - pow(r/0.5,2.0) );
			}
			else if(inletprofile==13) {	// periodic channel flow
				double w_bulk=inlet_flux/SumArea;
				w = 1.5 * w_bulk * yc * ( 2 - yc );
			}
			else if(inletprofile==15) w = 0;	// jet
			else if(inletprofile==16) {		// periodic pipe flow
				double r = sqrt(xc * xc + yc * yc);
				double w_bulk=inlet_flux/SumArea;
				w = w_bulk*2.*(1. - pow(r/0.5,2.0) );
			}
			else if(inletprofile==-9) w=0;
			else if(inlet_flux>0) {
				/*
				w=inlet_flux/user->k_area[k];*/
				w=inlet_flux/inletArea;
			}
			// else if (inletprofile==100) w=0;
			else w=1;
						
			ucont[k][j][i].z = w * area;
			//	printf("%f %f %f\n", w, inlet_flux,inletArea);
		       
			if(inletprofile==19) {
				double u=0, v=0, w=1;
				ucont[k][j][i].x = u * icsi[k][j][i].x  + v  * icsi[k][j][i].y + w * icsi[k][j][i].z;
				ucont[k][j][i].y = u * jeta[k][j][i].x  + v  * jeta[k][j][i].y + w * jeta[k][j][i].z;
				ucont[k][j][i].z = u * kzet[k][j][i].x  + v  * kzet[k][j][i].y + w * kzet[k][j][i].z;
			}
			
			
		}
	}
    
	srand( time(NULL)) ;	// seokkoo
	for (i = 0; i < (rand() % 3000); i++) (rand() % 3000);	//seokkoo
  
	if( initial_perturbation) {
		PetscPrintf(PETSC_COMM_WORLD, "\nGenerating initial perturbation\n");
		for(k=lzs; k<lze; k++) {	// for all CPU
			for (j=lys; j<lye; j++)
			for (i=lxs; i<lxe; i++) {
				if (nvert[k][j][i]+nvert[k+1][j][i] < threshold) {
				  int n1, n2, n3;
					double F;
					
					F  = 1.00; // 100%
					n1 = rand() % 20000 - 10000;
					n2 = rand() % 20000 - 10000;
					n3 = rand() % 20000 - 10000;
					/*
					if(inletprofile==13) {
						ucont[k][j][i].z += 10*((double)n)/10000.* F * kzet[k][j][i].z;
						
						F  = 0.1;
						n = rand() % 20000; n -= 10000;
						ucont[k][j][i].x = ((double)n)/10000. * F * icsi[k][j][i].x;
					
						F  = 0.1;
						n = rand() % 20000; n -= 10000;
						ucont[k][j][i].y = ((double)n)/10000. * F * jeta[k][j][i].y;
					}
					else */
					//ucont[k][j][i].z *= ( 1 + ((double)n1)/10000.*F );
					ucont[k][j][i].x = ((double)n3)/10000. * 0.1 * ucont[k][j][i].z;
					ucont[k][j][i].y = ((double)n2)/10000. * 0.1 * ucont[k][j][i].z;
					ucont[k][j][i].z *= ( 1 + ((double)n1)/10000.*F );
				}
			}
		}
	}
	
	if(levelset) DMDAVecRestoreArray(da, user->lLevelset, &level);
	DMDAVecRestoreArray(da, user->lNvert, &nvert);
	DMDAVecRestoreArray(da, user->P, &p);
	DMDAVecRestoreArray(fda, Coor, &coor);
	DMDAVecRestoreArray(fda, user->Ucont, &ucont);
	DMDAVecRestoreArray(fda, user->lICsi,  &icsi);
	DMDAVecRestoreArray(fda, user->lJEta,  &jeta);
	DMDAVecRestoreArray(fda, user->lKZet,  &kzet);
	DMDAVecRestoreArray(fda, user->lZet,  &zet);
	DMDAVecRestoreArray(fda, user->lCent,  &cent);
	
	DMGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont);
	DMGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont);
	
	DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
	DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);
	
	Contra2Cart(user);
		

	VecCopy(user->Ucont, user->Ucont_o);
        DMGlobalToLocalBegin(user->fda, user->Ucont_o, INSERT_VALUES, user->lUcont_o);
        DMGlobalToLocalEnd(user->fda, user->Ucont_o, INSERT_VALUES, user->lUcont_o);

        DMGlobalToLocalBegin(user->fda, user->Ucat, INSERT_VALUES, user->lUcat);
        DMGlobalToLocalEnd(user->fda, user->Ucat, INSERT_VALUES, user->lUcat);

        DMGlobalToLocalBegin(user->fda, user->Ucat, INSERT_VALUES, user->lUcat_old);
        DMGlobalToLocalEnd(user->fda, user->Ucat, INSERT_VALUES, user->lUcat_old);


	
	//char fname[80];
	//sprintf(fname,"Ucat");
	//TECIOOut_rhs(user, user->Ucat, fname);


/*
	if(initialzero) {
		VecSet(user->Ucont, 0);
		VecSet(user->Ucont_o, 0);
		VecSet(user->lUcont, 0);
		VecSet(user->Ucat, 0);
		VecSet(user->lUcat, 0);
	}
*/	
	return 0;
}



// xiaolei add  
// Inflow ajustment for subcritical flows
//

PetscErrorCode AjustInflow(UserCtx *user) 
{
	
  
	int i, j, k;
	Vec Coor;
	Cmpnts	***ucont, ***cent, ***ucat, ***coor, ***csi, ***eta, ***zet;
	Cmpnts ***icsi;
	
	DM da = user->da, fda = user->fda;
	DMDALocalInfo	info = user->info;
	PetscInt	xs = info.xs, xe = info.xs + info.xm;
	PetscInt  	ys = info.ys, ye = info.ys + info.ym;
	PetscInt	zs = info.zs, ze = info.zs + info.zm;
	PetscInt	mx = info.mx, my = info.my, mz = info.mz;
	PetscInt	lxs, lxe, lys, lye, lzs, lze;
	
	PetscReal	***nvert, ***level, ***aj;	
	
	lxs = xs; lxe = xe;
	lys = ys; lye = ye;
	lzs = zs; lze = ze;
  
	if (xs==0) lxs = xs+1;
	if (ys==0) lys = ys+1;
	if (zs==0) lzs = zs+1;
  
	if (xe==mx) lxe = xe-1;
	if (ye==my) lye = ye-1;
	if (ze==mz) lze = ze-1;
  
	DMDAGetGhostedCoordinates(da, &Coor);
	DMDAVecGetArray(fda, Coor, &coor);
	DMDAVecGetArray(fda, user->Ucont, &ucont);
	DMDAVecGetArray(fda, user->Ucat,  &ucat);
  
	DMDAVecGetArray(fda, user->lCsi,  &csi);
	DMDAVecGetArray(fda, user->lEta,  &eta);
	DMDAVecGetArray(fda, user->lZet,  &zet);
	DMDAVecGetArray(da, user->lAj,  &aj);

	DMDAVecGetArray(fda, user->lICsi,  &icsi);
 
	DMDAVecGetArray(fda, user->lCent,  &cent);	

	DMDAVecGetArray(da, user->lNvert, &nvert);
	if(levelset) DMDAVecGetArray(da, user->lLevelset, &level);










	DMDAVecRestoreArray(fda, Coor, &coor);
	DMDAVecRestoreArray(fda, user->Ucont, &ucont);
	DMDAVecRestoreArray(fda, user->Ucat,  &ucat);
  
	DMDAVecRestoreArray(fda, user->lCsi,  &csi);
	DMDAVecRestoreArray(fda, user->lEta,  &eta);
	DMDAVecRestoreArray(fda, user->lZet,  &zet);
	DMDAVecRestoreArray(da, user->lAj,  &aj);

	DMDAVecRestoreArray(fda, user->lICsi,  &icsi);
 
	DMDAVecRestoreArray(fda, user->lCent,  &cent);	

	DMDAVecRestoreArray(da, user->lNvert, &nvert);
	if(levelset) DMDAVecRestoreArray(da, user->lLevelset, &level);



}	




// add xiaolei
PetscErrorCode Scale_InitFlow(UserCtx *user)
{
	PetscInt i, j, k;

	Cmpnts	***csi, ***eta, ***zet;
	
	DM da = user->da, fda = user->fda;
	DMDALocalInfo	info = user->info;
	PetscInt	xs = info.xs, xe = info.xs + info.xm;
	PetscInt  	ys = info.ys, ye = info.ys + info.ym;
	PetscInt	zs = info.zs, ze = info.zs + info.zm;
	PetscInt	mx = info.mx, my = info.my, mz = info.mz;
	PetscInt	lxs, lxe, lys, lye, lzs, lze;
 
	PetscReal	***nvert, ***level, ***rho, ***aj;
	
	
	lxs = xs; lxe = xe;
	lys = ys; lye = ye;
	lzs = zs; lze = ze;
  
	if (xs==0) lxs = xs+1;
	if (ys==0) lys = ys+1;
	if (zs==0) lzs = zs+1;
  
	if (xe==mx) lxe = xe-1;
	if (ye==my) lye = ye-1;
	if (ze==mz) lze = ze-1;

	if(levelset)  {
		DMDAVecGetArray(da, user->lLevelset, &level);
		DMDAVecGetArray (da, user->lDensity, &rho);
	}
	DMDAVecGetArray(fda, user->lCsi, &csi);
	DMDAVecGetArray(fda, user->lEta, &eta);
	DMDAVecGetArray(fda, user->lZet, &zet);
	DMDAVecGetArray(da, user->lNvert, &nvert);	//seokkoo 
	DMDAVecGetArray(da, user->lAj, &aj);

	Cmpnts ***ucont;
	DMDAVecGetArray(fda, user->Ucont, &ucont);

	user->areak_air = new double [mz];
	user->areak_water = new double [mz];
	user->fluxk_air = new double [mz];
	user->fluxk_water = new double [mz];


	// Area         
	std::vector<double> lArea(mz), lArea_water(mz);
	std::vector<double> Sum_lArea(mz), Sum_lArea_water(mz);
	
	std::fill ( lArea.begin(), lArea.end(), 0 );
	std::fill ( lArea_water.begin(), lArea_water.end(), 0 );
	std::fill ( Sum_lArea.begin(), Sum_lArea.end(), 0 );
	std::fill ( Sum_lArea_water.begin(), Sum_lArea_water.end(), 0 );

      PetscPrintf(PETSC_COMM_WORLD,"Complete zeros out the Area\n"); 
       	for(k=lzs; k<lze; k++)
	for (j=ys; j<ye; j++)
	for (i=xs; i<xe; i++) {
		if(j>=1 && j<=my-2 && i>=1 && i<=mx-2) {
			double vf=1.0;
			double k_area = sqrt( zet[k][j][i].x*zet[k][j][i].x + zet[k][j][i].y*zet[k][j][i].y + zet[k][j][i].z*zet[k][j][i].z );
			if(levelset==1) {
				double dx = levelset_thickness ( aj[k][j][i], csi[k][j][i], eta[k][j][i], zet[k][j][i] );
				vf = H(level[k][j][i], dx);
			}
			else if(levelset==2) {
				vf = level[k][j][i];
			}
			if (nvert[k+1][j][i]+nvert[k][j][i] < 0.1) lArea[k] += k_area;
			if (nvert[k+1][j][i]+nvert[k][j][i] < 0.1) lArea_water[k] += k_area * vf;
		}
	}

        MPI_Allreduce( &lArea[0], &Sum_lArea[0], mz, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
        MPI_Allreduce( &lArea_water[0], &Sum_lArea_water[0], mz, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
	
	Sum_lArea[0] = Sum_lArea[1];
	Sum_lArea_water[0] = Sum_lArea_water[1];

	for(k=0;k<mz-1;k++) {
		user->areak_water[k]=Sum_lArea_water[k];
		user->areak_air[k]=Sum_lArea[k]-Sum_lArea_water[k];
	}

	// Flux
        
	std::vector<double> lFlux(mz), lFlux_water(mz);
	std::vector<double> Sum_lFlux(mz), Sum_lFlux_water(mz);
	
	std::fill ( lFlux.begin(), lFlux.end(), 0 );
	std::fill ( lFlux_water.begin(), lFlux_water.end(), 0 );
	std::fill ( Sum_lFlux.begin(), Sum_lFlux.end(), 0 );
	std::fill ( Sum_lFlux_water.begin(), Sum_lFlux_water.end(), 0 );

        for(k=lzs; k<lze; k++)
	for (j=ys; j<ye; j++)
	for (i=xs; i<xe; i++) {
		if(j>=1 && j<=my-2 && i>=1 && i<=mx-2) {
			double vf=1.0;
			if(levelset==1) {
				double dx = levelset_thickness ( aj[k][j][i], csi[k][j][i], eta[k][j][i], zet[k][j][i] );
				vf = H(level[k][j][i], dx);
			}
			else if(levelset==2) {
				vf = level[k][j][i];
			}
			if (nvert[k+1][j][i]+nvert[k][j][i] < 0.1) lFlux[k] += ucont[k][j][i].z;
			if (nvert[k+1][j][i]+nvert[k][j][i] < 0.1) lFlux_water[k] += ucont[k][j][i].z * vf;
		}
	}

        MPI_Allreduce( &lFlux[0], &Sum_lFlux[0], mz, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
        MPI_Allreduce( &lFlux_water[0], &Sum_lFlux_water[0], mz, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
	
	Sum_lFlux[0] = Sum_lFlux[1];
	Sum_lFlux_water[0] = Sum_lFlux_water[1];

	for(k=0;k<mz-1;k++) {
		user->fluxk_water[k]=Sum_lFlux_water[k];
		user->fluxk_air[k]=Sum_lFlux[k]-Sum_lFlux_water[k];
	}

	
	for (k=lzs; k<lze; k++) 
	for (j=lys; j<lye; j++) 
	for (i=lxs; i<lxe; i++) {			
		double vf=1.0;
		if(levelset==1) {
			double dx = levelset_thickness ( aj[k][j][i], csi[k][j][i], eta[k][j][i], zet[k][j][i] );
			vf = H(level[k][j][i], dx);
		}
		else if(levelset==2) {
			vf = level[k][j][i];
		}

		double Area = sqrt( zet[k-1][j][i].x*zet[k-1][j][i].x + zet[k-1][j][i].y*zet[k-1][j][i].y + zet[k-1][j][i].z*zet[k-1][j][i].z );
		double AreaSum=user->areak_water[k];
		double FluxInlet=user->fluxk_water[0];
		double Flux=user->fluxk_water[k];
		if (nvert[k+1][j][i]+nvert[k][j][i] < 0.1) ucont[k][j][i].z += (FluxInlet - Flux) * Area * vf / AreaSum;

		AreaSum=user->areak_air[k];
		FluxInlet=user->fluxk_air[0];
		Flux=user->fluxk_air[k];
		if (levelset && nvert[k+1][j][i]+nvert[k][j][i] < 0.1) ucont[k][j][i].z += (FluxInlet - Flux) * Area * (1.0-vf) / AreaSum;

	}		
		

      PetscPrintf(PETSC_COMM_WORLD,"Complete scaling Flux\n"); 

       	if(levelset) {
		DMDAVecRestoreArray (da, user->lLevelset, &level);
		DMDAVecRestoreArray (da, user->lDensity, &rho);
	}

      PetscPrintf(PETSC_COMM_WORLD,"End Levelset\n"); 



	DMDAVecRestoreArray(fda, user->lCsi, &csi);
	DMDAVecRestoreArray(fda, user->lEta, &eta);
	DMDAVecRestoreArray(fda, user->lZet, &zet);
	DMDAVecRestoreArray(da, user->lNvert, &nvert);	//seokkoo 
	DMDAVecRestoreArray(da, user->lAj, &aj);
	DMDAVecRestoreArray(fda, user->Ucont, &ucont);


      PetscPrintf(PETSC_COMM_WORLD,"Before ContratoCart\n"); 


  	Contra2Cart(user);


      PetscPrintf(PETSC_COMM_WORLD,"End ContraToCart\n"); 


};


