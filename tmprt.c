#include "variables.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
using namespace std;
extern PetscInt immersed, NumberOfBodies, ti, tistart, wallfunction;
extern PetscErrorCode MyKSPMonitor1(KSP ksp,PetscInt n,PetscReal rnorm,void *dummy);
extern PetscInt inletprofile;

//double solid=0.1;

void RHS_Tmprt(UserCtx *user, Vec Tmprt_RHS)
{
	DM		da = user->da, fda = user->fda, fda2 = user->fda2;
	DMDALocalInfo	info;
	PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
	PetscInt	mx, my, mz; // Dimensions in three directions
	PetscInt	i, j, k;
	
	PetscReal	***aj;
	
	PetscInt	lxs, lxe, lys, lye, lzs, lze;
	
	PetscReal	***tmprt, ***tmprt_o, ***tmprt_rhs;
	Cmpnts	***ucont, ***ucat;
	Cmpnts	***csi, ***eta, ***zet;
	PetscReal	***nvert, ***distance, ***lf1, ***lnu_t, ***pr_t;
 
	Vec Fp1, Fp2, Fp3;
	PetscReal ***fp1, ***fp2, ***fp3;
	PetscReal ***visc1, ***visc2, ***visc3;
	
	Cmpnts	***icsi, ***ieta, ***izet;
	Cmpnts	***jcsi, ***jeta, ***jzet;
	Cmpnts	***kcsi, ***keta, ***kzet;
	PetscReal	***iaj, ***jaj, ***kaj, ***rho, ***mu;

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
	
	VecDuplicate(user->lTmprt, &Fp1);
	VecDuplicate(user->lTmprt, &Fp2);
	VecDuplicate(user->lTmprt, &Fp3);
	
	VecSet(Fp1,0);
	VecSet(Fp2,0);
	VecSet(Fp3,0);
	VecSet(user->lVisc1_tmprt,0);
	VecSet(user->lVisc2_tmprt,0);
	VecSet(user->lVisc3_tmprt,0);
	
	
	if(levelset) {
		DMDAVecGetArray(da, user->lDensity, &rho);
		DMDAVecGetArray(da, user->lMu, &mu);
	}
	
	//DMDAVecGetArray(da, user->lSrans, &lS);
	if(les_prt) DMDAVecGetArray(da, user->lNu_t, &lnu_t);

	if(les_prt) DMDAVecGetArray(da, user->lPr_t, &pr_t);

	
	DMDAVecGetArray(fda, user->lUcont, &ucont);
	DMDAVecGetArray(fda, user->lUcat,  &ucat);
	
	DMDAVecGetArray(fda, user->lCsi, &csi);
	DMDAVecGetArray(fda, user->lEta, &eta);
	DMDAVecGetArray(fda, user->lZet, &zet);
	
	DMDAVecGetArray(fda, user->lICsi, &icsi);
	DMDAVecGetArray(fda, user->lIEta, &ieta);
	DMDAVecGetArray(fda, user->lIZet, &izet);
	DMDAVecGetArray(fda, user->lJCsi, &jcsi);
	DMDAVecGetArray(fda, user->lJEta, &jeta);
	DMDAVecGetArray(fda, user->lJZet, &jzet);
	DMDAVecGetArray(fda, user->lKCsi, &kcsi);
	DMDAVecGetArray(fda, user->lKEta, &keta);
	DMDAVecGetArray(fda, user->lKZet, &kzet);
	
	DMDAVecGetArray(da, user->lNvert, &nvert);
	
	DMDAVecGetArray(da, user->lAj, &aj);
	DMDAVecGetArray(da, user->lIAj, &iaj);
	DMDAVecGetArray(da, user->lJAj, &jaj);
	DMDAVecGetArray(da, user->lKAj, &kaj);
	
	DMDAVecGetArray(da, user->lTmprt, &tmprt);
	DMDAVecGetArray(da, user->lTmprt_o, &tmprt_o);
	DMDAVecGetArray(da, Tmprt_RHS, &tmprt_rhs);
		
	DMDAVecGetArray(da, Fp1, &fp1);
	DMDAVecGetArray(da, Fp2, &fp2);
	DMDAVecGetArray(da, Fp3, &fp3);
	DMDAVecGetArray(da, user->lVisc1_tmprt, &visc1);
	DMDAVecGetArray(da, user->lVisc2_tmprt, &visc2);
	DMDAVecGetArray(da, user->lVisc3_tmprt, &visc3);		
	

	double nu_t, D_m, D_t;

	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs-1; i<lxe; i++) {
		double csi0 = icsi[k][j][i].x, csi1 = icsi[k][j][i].y, csi2 = icsi[k][j][i].z;
		double eta0 = ieta[k][j][i].x, eta1 = ieta[k][j][i].y, eta2 = ieta[k][j][i].z;
		double zet0 = izet[k][j][i].x, zet1 = izet[k][j][i].y, zet2 = izet[k][j][i].z;
		double ajc = iaj[k][j][i];
		
		double g11 = csi0 * csi0 + csi1 * csi1 + csi2 * csi2;
		double g21 = eta0 * csi0 + eta1 * csi1 + eta2 * csi2;
		double g31 = zet0 * csi0 + zet1 * csi1 + zet2 * csi2;
		
		double dtdc, dtde, dtdz;
		
		
		dtdc = tmprt[k][j][i+1] - tmprt[k][j][i];

		if ((nvert[k][j+1][i])> 1.1 || (nvert[k][j+1][i+1])> 1.1) {
			dtde = (tmprt[k][j  ][i+1] + tmprt[k][j  ][i] - tmprt[k][j-1][i+1] - tmprt[k][j-1][i]) * 0.5;
		}
		else if  ((nvert[k][j-1][i])> 1.1 || (nvert[k][j-1][i+1])> 1.1) {
			dtde = (tmprt[k][j+1][i+1] + tmprt[k][j+1][i] - tmprt[k][j  ][i+1] - tmprt[k][j  ][i]) * 0.5;
		}
		else {
			dtde = (tmprt[k][j+1][i+1] + tmprt[k][j+1][i] - tmprt[k][j-1][i+1] - tmprt[k][j-1][i]) * 0.25;
		}	  

                if ((nvert[k+1][j][i])> 1.1 || (nvert[k+1][j][i+1])> 1.1) {
                        dtdz = (tmprt[k  ][j][i+1] + tmprt[k  ][j][i] - tmprt[k-1][j][i+1] - tmprt[k-1][j][i]) * 0.5;
			
		}
		else if ((nvert[k-1][j][i])> 1.1 || (nvert[k-1][j][i+1])> 1.1) {
			dtdz = (tmprt[k+1][j][i+1] + tmprt[k+1][j][i] - tmprt[k  ][j][i+1] - tmprt[k  ][j][i]) * 0.5;
		}
		else {
			dtdz = (tmprt[k+1][j][i+1] + tmprt[k+1][j][i] - tmprt[k-1][j][i+1] - tmprt[k-1][j][i]) * 0.25;
		}
		
//		double wL=1./aj[k][j][i], wR=1./aj[k][j][i+1];
//		double Tm_o = ( tmprt_o[k][j][i]*wL + tmprt_o[k][j][i+1]*wR ) / (wL+wR); //	Tm_o = PetscMax ( Tm_o, 0 );
//		double sigma_k=0, sigma_o=0;
		
		if (les_prt) {
                        if ( i==0 || nvert[k][j][i]>0.1 ) D_t = pr_t[k][j][i+1];
                        else if( i==mx-2 || nvert[k][j][i+1]>0.1 ) D_t = pr_t[k][j][i];
                        else D_t = 0.5 * (pr_t[k][j][i] + pr_t[k][j][i+1]);
		}

//		if( nvert[k][j][i]+nvert[k][j][i+1]>0.1 || i==0 || i==mx-2 || periodic ) {
//			fp1[k][j][i] = -ucont[k][j][i].x * Upwind ( tmprt[k][j][i], tmprt[k][j][i+1], ucont[k][j][i].x);
//		}
//		else {
                        fp1[k][j][i] = - ucont[k][j][i].x * 0.5 * ( tmprt[k][j][i] + tmprt[k][j][i+1]);
//			fp1[k][j][i] = -ucont[k][j][i].x * weno3 ( tmprt[k][j][i-1], tmprt[k][j][i], tmprt[k][j][i+1], tmprt[k][j][i+2], ucont[k][j][i].x );
//		}
		
		D_m = 1./ (user->ren * user->Pr);
//		D_t = D_m;
		if(levelset) {
//			if(nvert[k][j][i]>0.1 || i==0) nu=mu[k][j][i+1];
//			else if(nvert[k][j][i+1]>0.1 || i==mx-2) nu=mu[k][j][i];
//			else nu = 0.5 * ( mu[k][j][i] + mu[k][j][i+1] );
		}
		
		visc1[k][j][i] = (g11 * dtdc + g21 * dtde + g31 * dtdz) * ajc * (D_m + D_t);
	}
	
	for (k=lzs; k<lze; k++)
	for (j=lys-1; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
		double csi0 = jcsi[k][j][i].x, csi1 = jcsi[k][j][i].y, csi2 = jcsi[k][j][i].z;
		double eta0= jeta[k][j][i].x, eta1 = jeta[k][j][i].y, eta2 = jeta[k][j][i].z;
		double zet0 = jzet[k][j][i].x, zet1 = jzet[k][j][i].y, zet2 = jzet[k][j][i].z;
		double ajc = jaj[k][j][i];
		
		double g11 = csi0 * eta0 + csi1 * eta1 + csi2 * eta2;
		double g21 = eta0 * eta0 + eta1 * eta1 + eta2 * eta2;
		double g31 = zet0 * eta0 + zet1 * eta1 + zet2 * eta2;
		
		double dtdc, dtde, dtdz;
		
		if ((nvert[k][j][i+1])> 1.1 || (nvert[k][j+1][i+1])> 1.1) {
			dtdc = (tmprt[k][j+1][i  ] + tmprt[k][j][i  ] - tmprt[k][j+1][i-1] - tmprt[k][j][i-1]) * 0.5;
		}
		else if ((nvert[k][j][i-1])> 1.1 || (nvert[k][j+1][i-1])> 1.1) {
			dtdc = (tmprt[k][j+1][i+1] + tmprt[k][j][i+1] - tmprt[k][j+1][i  ] - tmprt[k][j][i  ]) * 0.5;
		}
		else {
			dtdc = (tmprt[k][j+1][i+1] + tmprt[k][j][i+1] - tmprt[k][j+1][i-1] - tmprt[k][j][i-1]) * 0.25;
		}

		//if (j==0) printf("tmprt bc at j==0 %le \n", tmprt[k][j][i]);
		//if (j==1) printf("tmprt bc at j==1 %le \n", tmprt[k][j][i]);
		dtde = tmprt[k][j+1][i] - tmprt[k][j][i];

		if ((nvert[k+1][j][i])> 1.1 || (nvert[k+1][j+1][i])> 1.1) {
			dtdz = (tmprt[k  ][j+1][i] + tmprt[k  ][j][i] - tmprt[k-1][j+1][i] - tmprt[k-1][j][i]) * 0.5;
		}
		else if ((nvert[k-1][j][i])> 1.1 || (nvert[k-1][j+1][i])> 1.1) {
			dtdz = (tmprt[k+1][j+1][i] + tmprt[k+1][j][i] - tmprt[k  ][j+1][i] - tmprt[k  ][j][i]) * 0.5;
		}
		else {
			dtdz = (tmprt[k+1][j+1][i] + tmprt[k+1][j][i] - tmprt[k-1][j+1][i] - tmprt[k-1][j][i]) * 0.25;
		}
		
//		double wL=1./aj[k][j][i], wR=1./aj[k][j+1][i];
//		double Tm_o = ( tmprt_o[k][j][i]*wL + tmprt_o[k][j+1][i]*wR ) / (wL+wR); //	Tm_o = PetscMax ( Tm_o, 0 );
//		double sigma_k=0, sigma_o=0;


		if (les_prt) {
                        if ( j==0 || nvert[k][j][i]>0.1 ) D_t = pr_t[k][j+1][i];
                        else if( j==my-2 || nvert[k][j+1][i]>0.1 ) D_t = pr_t[k][j][i];
                        else D_t = 0.5 * (pr_t[k][j][i] + pr_t[k][j+1][i]);
		}

//		if( nvert[k][j][i]+nvert[k][j+1][i]>0.1 || j==0 || j==my-2 || periodic ) {
//			fp2[k][j][i] = -ucont[k][j][i].y * Upwind ( tmprt[k][j][i], tmprt[k][j+1][i], ucont[k][j][i].y);
//		}
//		else {
                        fp2[k][j][i] = -ucont[k][j][i].y * 0.5 * ( tmprt[k][j][i] + tmprt[k][j+1][i]);
//			fp2[k][j][i] = -ucont[k][j][i].y * weno3 ( tmprt[k][j-1][i], tmprt[k][j][i], tmprt[k][j+1][i], tmprt[k][j+2][i], ucont[k][j][i].y );
//		}
		
		D_m = 1./(user->ren * user->Pr);
//		D_t = D_m;
		if(levelset) {
//			if(nvert[k][j][i]>0.1 || j==0) nu=mu[k][j+1][i];
//			else if(nvert[k][j+1][i]>0.1 || j==my-2) nu=mu[k][j][i];
//			else nu = 0.5 * ( mu[k][j][i] + mu[k][j+1][i] );
		}
		
		visc2[k][j][i] = (g11 * dtdc + g21 * dtde + g31 * dtdz) * ajc * (D_m + D_t);
		
	}
	
	for (k=lzs-1; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
		double csi0 = kcsi[k][j][i].x, csi1 = kcsi[k][j][i].y, csi2 = kcsi[k][j][i].z;
		double eta0 = keta[k][j][i].x, eta1 = keta[k][j][i].y, eta2 = keta[k][j][i].z;
		double zet0 = kzet[k][j][i].x, zet1 = kzet[k][j][i].y, zet2 = kzet[k][j][i].z;
		double ajc = kaj[k][j][i];
		
		double g11 = csi0 * zet0 + csi1 * zet1 + csi2 * zet2;
		double g21 = eta0 * zet0 + eta1 * zet1 + eta2 * zet2;
		double g31 = zet0 * zet0 + zet1 * zet1 + zet2 * zet2;
		
		double dtdc, dtde, dtdz;
		double nu_t, prt;
		
		if ((nvert[k][j][i+1])> 1.1 || (nvert[k+1][j][i+1])> 1.1) {
			dtdc = (tmprt[k+1][j][i  ] + tmprt[k][j][i  ] - tmprt[k+1][j][i-1] - tmprt[k][j][i-1]) * 0.5;
		}
		else if ((nvert[k][j][i-1])> 1.1 || (nvert[k+1][j][i-1])> 1.1) {
			dtdc = (tmprt[k+1][j][i+1] + tmprt[k][j][i+1] - tmprt[k+1][j][i  ] - tmprt[k][j][i  ]) * 0.5;
		}
		else {
			dtdc = (tmprt[k+1][j][i+1] + tmprt[k][j][i+1] - tmprt[k+1][j][i-1] - tmprt[k][j][i-1]) * 0.25;
		}

		if ((nvert[k][j+1][i])> 1.1 || (nvert[k+1][j+1][i])> 1.1) {
			dtde = (tmprt[k+1][j  ][i] + tmprt[k][j  ][i] - tmprt[k+1][j-1][i] - tmprt[k][j-1][i]) * 0.5;
		}
		else if ((nvert[k][j-1][i])> 1.1 || (nvert[k+1][j-1][i])> 1.1) {
			dtde = (tmprt[k+1][j+1][i] + tmprt[k][j+1][i] - tmprt[k+1][j  ][i] - tmprt[k][j  ][i]) * 0.5;
		}
		else {
			dtde = (tmprt[k+1][j+1][i] + tmprt[k][j+1][i] - tmprt[k+1][j-1][i] - tmprt[k][j-1][i]) * 0.25;
		}

		dtdz = tmprt[k+1][j][i] - tmprt[k][j][i];
		
//		double wL=1./aj[k][j][i], wR=1./aj[k+1][j][i];
//		double Tm_o = ( tmprt_o[k][j][i]*wL + tmprt_o[k+1][j][i]*wR ) / (wL+wR); //	Tm_o = PetscMax ( Tm_o, 0 );
//		double sigma_k=0, sigma_o=0;
		
		if (les_prt) {
                        if ( k==0 || nvert[k][j][i]>0.1 ) D_t = pr_t[k+1][j][i];
                        else if( k==mz-2 || nvert[k+1][j][i]>0.1 ) D_t = pr_t[k][j][i];
                        else D_t = 0.5 * (pr_t[k][j][i] + pr_t[k+1][j][i]);
		}

		
//		if( nvert[k][j][i]+nvert[k+1][j][i]>0.1 || k==0 || k==mz-2 || periodic ) {
//			fp3[k][j][i] = -ucont[k][j][i].z * Upwind ( tmprt[k][j][i], tmprt[k+1][j][i], ucont[k][j][i].z);
//		}
//		else {
                        fp3[k][j][i] = -ucont[k][j][i].z * 0.5 * ( tmprt[k][j][i] + tmprt[k+1][j][i]);
//			fp3[k][j][i] = -ucont[k][j][i].z * weno3 ( tmprt[k-1][j][i], tmprt[k][j][i], tmprt[k+1][j][i], tmprt[k+2][j][i], ucont[k][j][i].z );
//:		}
		
		D_m = 1./ (user->ren * user->Pr);
//		D_t = D_m;
		if(levelset) {
//			if(nvert[k][j][i]>0.1 || k==0) nu=mu[k+1][j][i];
//			else if(nvert[k+1][j][i]>0.1 || k==mz-2) nu=mu[k][j][i];
//			else nu = 0.5 * ( mu[k][j][i] + mu[k+1][j][i] );
		}
		
		visc3[k][j][i] = (g11 * dtdc + g21 * dtde + g31 * dtdz) * ajc * (D_m + D_t);
	}
	
	DMDAVecRestoreArray(da, Fp1, &fp1);
	DMDAVecRestoreArray(da, Fp2, &fp2);
	DMDAVecRestoreArray(da, Fp3, &fp3);
	DMDAVecRestoreArray(da, user->lVisc1_tmprt, &visc1);
	DMDAVecRestoreArray(da, user->lVisc2_tmprt, &visc2);
	DMDAVecRestoreArray(da, user->lVisc3_tmprt, &visc3);

	DMDALocalToLocalBegin(da, Fp1, INSERT_VALUES, Fp1);
	DMDALocalToLocalEnd(da, Fp1, INSERT_VALUES, Fp1);
	
	DMDALocalToLocalBegin(da, Fp2, INSERT_VALUES, Fp2);
	DMDALocalToLocalEnd(da, Fp2, INSERT_VALUES, Fp2);
	
	DMDALocalToLocalBegin(da, Fp3, INSERT_VALUES, Fp3);
	DMDALocalToLocalEnd(da, Fp3, INSERT_VALUES, Fp3);
	
	DMDALocalToLocalBegin(da, user->lVisc1_tmprt, INSERT_VALUES, user->lVisc1_tmprt);
	DMDALocalToLocalEnd(da, user->lVisc1_tmprt, INSERT_VALUES, user->lVisc1_tmprt);
	
	DMDALocalToLocalBegin(da, user->lVisc2_tmprt, INSERT_VALUES, user->lVisc2_tmprt);
	DMDALocalToLocalEnd(da, user->lVisc2_tmprt, INSERT_VALUES, user->lVisc2_tmprt);
	
	DMDALocalToLocalBegin(da, user->lVisc3_tmprt, INSERT_VALUES, user->lVisc3_tmprt);
	DMDALocalToLocalEnd(da, user->lVisc3_tmprt, INSERT_VALUES, user->lVisc3_tmprt);
	
	DMDAVecGetArray(da, Fp1, &fp1);
	DMDAVecGetArray(da, Fp2, &fp2);
	DMDAVecGetArray(da, Fp3, &fp3);
	DMDAVecGetArray(da, user->lVisc1_tmprt, &visc1);
	DMDAVecGetArray(da, user->lVisc2_tmprt, &visc2);
	DMDAVecGetArray(da, user->lVisc3_tmprt, &visc3);
	
	if(periodic)
	for (k=zs; k<ze; k++)
	for (j=ys; j<ye; j++)
	for (i=xs; i<xe; i++) {
		int a=i, b=j, c=k;

		int flag=0;
		
		if(i_periodic && i==0) a=mx-2, flag=1;
		else if(i_periodic && i==mx-1) a=1, flag=1;
		
		if(j_periodic && j==0) b=my-2, flag=1;
		else if(j_periodic && j==my-1) b=1, flag=1;
		
		if(k_periodic && k==0) c=mz-2, flag=1;
		else if(k_periodic && k==mz-1) c=1, flag=1;
		
		if(ii_periodic && i==0) a=-2, flag=1;
		else if(ii_periodic && i==mx-1) a=mx+1, flag=1;
		
		if(jj_periodic && j==0) b=-2, flag=1;
		else if(jj_periodic && j==my-1) b=my+1, flag=1;
		
		if(kk_periodic && k==0) c=-2, flag=1;
		else if(kk_periodic && k==mz-1) c=mz+1, flag=1;
				
		if(flag) {
			fp1[k][j][i] = fp1[c][b][a];
			fp2[k][j][i] = fp2[c][b][a];
			fp3[k][j][i] = fp3[c][b][a];
			visc1[k][j][i] = visc1[c][b][a];
			visc2[k][j][i] = visc2[c][b][a];
			visc3[k][j][i] = visc3[c][b][a];
		}
	}
	
	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
		if ( i==0 || i==mx-1 || j==0 || j==my-1 || k==0 || k==mz-1 || nvert[k][j][i]>0.1 ) {
			tmprt_rhs[k][j][i] = 0.0 ;
			continue;
		}
			
		double dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz;
		double du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz;
		
//		double Tm = tmprt[k][j][i];
//		double Tm_o = tmprt_o[k][j][i];
		
//		Km_o = PetscMax ( Km_o, 0 );
//		Om_o = PetscMax ( Om_o, 1.e-5 );

		double csi0 = csi[k][j][i].x, csi1 = csi[k][j][i].y, csi2 = csi[k][j][i].z;
		double eta0= eta[k][j][i].x, eta1 = eta[k][j][i].y, eta2 = eta[k][j][i].z;
		double zet0 = zet[k][j][i].x, zet1 = zet[k][j][i].y, zet2 = zet[k][j][i].z;
		double ajc = aj[k][j][i];//, d = distance[k][j][i];
		

		double dtdc, dtde, dtdz;
		double dt_dx, dt_dy, dt_dz;
		
		Compute_dscalar_center(i, j, k, mx, my, mz, tmprt, nvert,  &dtdc, &dtde, &dtdz );
		Compute_dscalar_dxyz (csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc,dtdc, dtde, dtdz, &dt_dx, &dt_dy, &dt_dz );
		
		Compute_du_center (i, j, k, mx, my, mz, ucat, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz); 

		Compute_du_dxyz (	csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz,
							&du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
		
		double Sxx = 0.5*( du_dx + du_dx ), Sxy = 0.5*(du_dy + dv_dx), Sxz = 0.5*(du_dz + dw_dx);
		double Syx = Sxy, Syy = 0.5*(dv_dy + dv_dy),	Syz = 0.5*(dv_dz + dw_dy);
		double Szx = Sxz, Szy=Syz, Szz = 0.5*(dw_dz + dw_dz);
		double S = sqrt( 2.0*( Sxx*Sxx + Sxy*Sxy + Sxz*Sxz + Syx*Syx + Syy*Syy + Syz*Syz + Szx*Szx + Szy*Szy + Szz*Szz) );
		
		double nu_t, inv_nu_t, f1=1;
		
		
							
		if ( nvert[k][j][i] < 0.1 ) {
			double r = 1.;
			
			if(levelset) r = rho[k][j][i];
			
			//if (j==1) printf("Visc2 bc at j==0 %le \n", visc2[k][j-1][i]);
			//if (j==1) printf("Visc2 bc at j==1 %le \n", visc2[k][j][i]);
			tmprt_rhs[k][j][i] = ( fp1[k][j][i] - fp1[k][j][i-1] + fp2[k][j][i] - fp2[k][j-1][i] + fp3[k][j][i] - fp3[k-1][j][i] ) * ajc;	// advection
			tmprt_rhs[k][j][i] += ( visc1[k][j][i] - visc1[k][j][i-1] + visc2[k][j][i] - visc2[k][j-1][i] + visc3[k][j][i] - visc3[k-1][j][i]
) * ajc / r;	// diffusion
			

		}

		//if (j==0) printf("rhs at j==0 %le \n", tmprt_rhs[k][j][i]);
		//if (i==10 && j==1 && k==10) printf("in RHS Calc rhs at i=10 j==1 k=10 %le \n", tmprt_rhs[k][j][i]);
		//if (i==10 && j==50 && k==10) printf("in RHS Calc rhs at i=10 j==50 k=10 %le \n", tmprt_rhs[k][j][i]);
		//if (j==2) printf("rhs at j==2 %le \n", tmprt_rhs[k][j][i]);


	}

	/*
	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {

		if (i==10 && j==1 && k==10) printf("in RHS Calc rhs at i=10 j==1 k=10 %le \n", tmprt_rhs[k][j][i]);
		if (i==10 && j==50 && k==10) printf("in RHS Calc rhs at i=10 j==50 k=10 %le \n", tmprt_rhs[k][j][i]);
	}	
	*/	


	if (les_prt) DMDAVecRestoreArray(da, user->lNu_t, &lnu_t);
	if(les_prt) DMDAVecRestoreArray(da, user->lPr_t, &pr_t);

	DMDAVecRestoreArray(da, Fp1, &fp1);
	DMDAVecRestoreArray(da, Fp2, &fp2);
	DMDAVecRestoreArray(da, Fp3, &fp3);
	DMDAVecRestoreArray(da, user->lVisc1_tmprt, &visc1);
	DMDAVecRestoreArray(da, user->lVisc2_tmprt, &visc2);
	DMDAVecRestoreArray(da, user->lVisc3_tmprt, &visc3);
	
	DMDAVecRestoreArray(fda, user->lUcont, &ucont);
	DMDAVecRestoreArray(fda, user->lUcat,  &ucat);
	
	DMDAVecRestoreArray(fda, user->lCsi, &csi);
	DMDAVecRestoreArray(fda, user->lEta, &eta);
	DMDAVecRestoreArray(fda, user->lZet, &zet);
	
	DMDAVecRestoreArray(fda, user->lICsi, &icsi);
	DMDAVecRestoreArray(fda, user->lIEta, &ieta);
	DMDAVecRestoreArray(fda, user->lIZet, &izet);
	DMDAVecRestoreArray(fda, user->lJCsi, &jcsi);
	DMDAVecRestoreArray(fda, user->lJEta, &jeta);
	DMDAVecRestoreArray(fda, user->lJZet, &jzet);
	DMDAVecRestoreArray(fda, user->lKCsi, &kcsi);
	DMDAVecRestoreArray(fda, user->lKEta, &keta);
	DMDAVecRestoreArray(fda, user->lKZet, &kzet);
	
	DMDAVecRestoreArray(da, user->lNvert, &nvert);
	
	DMDAVecRestoreArray(da, user->lAj, &aj);
	DMDAVecRestoreArray(da, user->lIAj, &iaj);
	DMDAVecRestoreArray(da, user->lJAj, &jaj);
	DMDAVecRestoreArray(da, user->lKAj, &kaj);
	
	DMDAVecRestoreArray(da, user->lTmprt, &tmprt);
	DMDAVecRestoreArray(da, user->lTmprt_o, &tmprt_o);
	DMDAVecRestoreArray(da, Tmprt_RHS, &tmprt_rhs);
	
	if(levelset) {
		DMDAVecRestoreArray(da, user->lDensity, &rho);
		DMDAVecRestoreArray(da, user->lMu, &mu);
	}


//	TECIOOut_rhs1(user, Tmprt_RHS);
//        int aaa;
//        cout << "here \n";
//        cin >> aaa;


	VecDestroy(&Fp1);
	VecDestroy(&Fp2);
	VecDestroy(&Fp3);
};


void Tmprt_IC(UserCtx *user)
{
	DM		da = user->da;
	DMDALocalInfo	info;
	PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
	PetscInt	mx, my, mz; // Dimensions in three directions
	PetscInt	i, j, k;
	double nu = 1./(user->ren * user->Pr);
	
	PetscInt	lxs, lxe, lys, lye, lzs, lze;

	
	PetscReal 	***tmprt;
	PetscReal	***nvert;

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

	DMDAVecGetArray(da, user->lNvert, &nvert);
	DMDAVecGetArray(user->da, user->lTmprt, &tmprt);

        Cmpnts          ***coor;
        Vec             Coor;

        DMDAGetCoordinates(da, &Coor);

	DMDAVecGetArray(user->fda, Coor, &coor);

	
	for (k=zs; k<ze; k++)
	for (j=ys; j<ye; j++)
	for (i=xs; i<xe; i++) {	// pressure node; cell center

		//tmprt[k][j][i] = coor[k][j][i].y*user->tmprts[3];
		tmprt[k][j][i] = 0.0;

                if(initial_perturbation) {
                  	double F = 1.0;                        // 3% noise
                   	int n = rand() % 20000;         // RAND_MAX = 65535
                      	n -= 10000;
//                 	tmprt[k][j][i] *= 0.0; //( 1.0 + ((double)n)/10000. * F );     // uin * (1+-0.xx)
          	}


		
		if( nvert[k][j][i]>0.1) {
//			tmprt[k][j][i] = 0.0;
		}
		else if(i==0 || i==mx-1 || j==0 || j==my-1 || k==0 || k==mz-1 ) {
//			tmprt[k][j][i] = 0.0;
		}


	}
	
	DMDAVecRestoreArray(user->fda, Coor, &coor);
	DMDAVecRestoreArray(da, user->lNvert, &nvert);
	DMDAVecRestoreArray(user->da, user->lTmprt, &tmprt);

//	Tmprt_BC(user);	

	DMLocalToGlobalBegin(user->da, user->lTmprt, INSERT_VALUES, user->Tmprt);
	DMLocalToGlobalEnd(user->da, user->lTmprt, INSERT_VALUES, user->Tmprt);
	VecCopy(user->Tmprt, user->Tmprt_o);


	
};


void Force_Tmprt(UserCtx *user)
{
	DM		da = user->da;
	DMDALocalInfo	info;
	PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
	PetscInt	mx, my, mz; // Dimensions in three directions
	PetscInt	i, j, k;
	double nu = 1./(user->ren * user->Pr);
	double Ri_x = user->Ri_x;
	double Ri_y = user->Ri_y;
	double Ri_z = user->Ri_z;

	
	PetscInt	lxs, lxe, lys, lye, lzs, lze;

	
	PetscReal 	***tmprt;
	PetscReal	***nvert;

	Cmpnts		***ftmprt, ***csi, ***eta, ***zet;

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

	DMDAVecGetArray(da, user->lNvert, &nvert);
	DMDAVecGetArray(user->da, user->lTmprt, &tmprt);
	DMDAVecGetArray(user->fda, user->lFTmprt, &ftmprt);
  	DMDAVecGetArray(user->fda, user->lCsi,  &csi);
  	DMDAVecGetArray(user->fda, user->lEta,  &eta);
  	DMDAVecGetArray(user->fda, user->lZet,  &zet);

//	PetscReal tmprt_xzAvg[my];

//	for (j=0; j<my; j++) {
//		tmprt_xzAvg[j] = 0.0;
//	}

	
//        for (k=lzs; k<lze; k++)
//        for (j=lys; j<lye; j++)
//        for (i=lxs; i<lxe; i++) {       // pressure node; cell center

//		tmprt_xzAvg[j]+=tmprt[k][j][i];

//	}


//  	PetscBarrier(PETSC_NULL);

//  	PetscReal     u_sum;
//    	for (j=0; j<my; j++) {
      	
//		PetscGlobalSum(&(tmprt_xzAvg[j]), &u_sum, PETSC_COMM_WORLD);
//      		tmprt_xzAvg[j] = u_sum/((double)(mx*mz));

//		printf("the avg temperature %d %le \n", j, tmprt_xzAvg[j]);
//    	}



	double force_x, force_y, force_z;
	
	for (k=zs; k<ze; k++)
	for (j=ys; j<ye; j++)
	for (i=xs; i<xe; i++) {	// pressure node; cell center

		force_x = 0.0;
		force_y = 0.0;
		force_z = 0.0;

		force_x = Ri_x * tmprt[k][j][i]; 
		force_y = Ri_y * tmprt[k][j][i]; 
		force_z = Ri_z * tmprt[k][j][i]; 

	      	ftmprt[k][j][i].x = force_x *  csi[k][j][i].x + force_y *  csi[k][j][i].y + force_z *  csi[k][j][i].z;
	      	ftmprt[k][j][i].y = force_x *  eta[k][j][i].x + force_y *  eta[k][j][i].y + force_z *  eta[k][j][i].z;
	      	ftmprt[k][j][i].z = force_x *  zet[k][j][i].x + force_y *  zet[k][j][i].y + force_z *  zet[k][j][i].z;

//		if (k==3 && i==3) printf("the force  %d %d %le %le %le %le\n", my-2, j, ftmprt[k][j][i].y, force_x, force_y, force_z);

                if( nvert[k][j][i]>0.1) {
                  	ftmprt[k][j][i].x = 0.0;
                	ftmprt[k][j][i].y = 0.0;
                	ftmprt[k][j][i].z = 0.0;
                }
//                else if(i==0 || i==mx-1 || j==0 || j==my-1 || k==0 || k==mz-1 || j==my-2) {
//                        ftmprt[k][j][i].x = 0.0;
//                        ftmprt[k][j][i].y = 0.0;
//                        ftmprt[k][j][i].z = 0.0;
//                }

//		if (user->bctype_tmprt[1] && i==mx-2) ftmprt[k][j][i].x = 0.0;
//		if (user->bctype_tmprt[3] && j==my-2) ftmprt[k][j][i].y = 0.0;
//		if (user->bctype_tmprt[5] && k==mz-2) ftmprt[k][j][i].z = 0.0;


	}


	if (xs ==0) {
		i = 0;
		for (k=zs; k<ze; k++)
		for (j=ys; j<ye; j++) {
			ftmprt[k][j][i].x = 0;
			ftmprt[k][j][i].y = 0;
			ftmprt[k][j][i].z = 0;
		}
	}

	if (xe == mx) {
		for (k=zs; k<ze; k++)
		for (j=ys; j<ye; j++) {
			if(!i_periodic && !ii_periodic) {
				i = mx-2;
				ftmprt[k][j][i].x = 0;
			}
			i = mx-1;
			ftmprt[k][j][i].x = 0;
			ftmprt[k][j][i].y = 0;
			ftmprt[k][j][i].z = 0;
		}
	}


	if (ys == 0) {
		for (k=zs; k<ze; k++)
		for (i=xs; i<xe; i++) {
			j=0;
			ftmprt[k][j][i].x = 0;
			ftmprt[k][j][i].y = 0;
			ftmprt[k][j][i].z = 0;
		}
	}
  
	if (ye == my) {
		for (k=zs; k<ze; k++) 
		for (i=xs; i<xe; i++) {
			if(!j_periodic && !jj_periodic) {
				j=my-2;
				ftmprt[k][j][i].y = 0;
			}
			j=my-1;
			ftmprt[k][j][i].x = 0;
			ftmprt[k][j][i].y = 0;
			ftmprt[k][j][i].z = 0;
		}
	}
	
	
	if (zs == 0) {
		k=0;
		for (j=ys; j<ye; j++)
		for (i=xs; i<xe; i++) {
			ftmprt[k][j][i].x = 0;
			ftmprt[k][j][i].y = 0;
			ftmprt[k][j][i].z = 0;
		}
	}
  
	if (ze == mz) {
		for (j=ys; j<ye; j++)
		for (i=xs; i<xe; i++) {
			if(!k_periodic && !kk_periodic) {
				k=mz-2;
				ftmprt[k][j][i].z = 0;
			}
			k=mz-1;
			ftmprt[k][j][i].x = 0;
			ftmprt[k][j][i].y = 0;
			ftmprt[k][j][i].z = 0;
		}
	}


	DMDAVecRestoreArray(da, user->lNvert, &nvert);
	DMDAVecRestoreArray(user->da, user->lTmprt, &tmprt);
	DMDAVecRestoreArray(user->fda, user->lFTmprt, &ftmprt);

	DMLocalToGlobalBegin(user->fda, user->lFTmprt, INSERT_VALUES, user->FTmprt);
	DMLocalToGlobalEnd(user->fda, user->lFTmprt, INSERT_VALUES, user->FTmprt);

//	TECIOOut_rhs(user,  user->FTmprt);

        DMDAVecRestoreArray(user->fda, user->lCsi,  &csi);
        DMDAVecRestoreArray(user->fda, user->lEta,  &eta);
        DMDAVecRestoreArray(user->fda, user->lZet,  &zet);
	
};


void Tmprt_BC(UserCtx *user)
{
	DM		da = user->da, fda = user->fda;
	DMDALocalInfo	info;
	PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
	PetscInt	mx, my, mz; // Dimensions in three directions
	PetscInt	i, j, k, ibi;
	//double nu = 1./user->ren;
	PetscReal	***aj, ***distance, ***rho, ***mu;
	
	PetscInt	lxs, lxe, lys, lye, lzs, lze;

	
	PetscReal ***tmprt;
	Cmpnts	***csi, ***eta, ***zet, ***ucat;
	PetscReal	***nvert, ***ustar;


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

	DMDAVecGetArray(fda, user->lCsi, &csi);
	DMDAVecGetArray(fda, user->lEta, &eta);
	DMDAVecGetArray(fda, user->lZet, &zet);
	DMDAVecGetArray(fda, user->lUcat, &ucat);
	
	DMDAVecGetArray(da, user->lNvert, &nvert);
	DMDAVecGetArray(da, user->lAj, &aj);
	DMDAVecGetArray(da, user->lUstar, &ustar);
	
	DMDAVecGetArray(user->da, user->lTmprt, &tmprt);
	
	if(levelset) {
		DMDAVecGetArray(da, user->lDensity, &rho);
		DMDAVecGetArray(da, user->lMu, &mu);
	}


	
	// BC for temperature
	for (k=zs; k<ze; k++)
	for (j=ys; j<ye; j++)
	for (i=xs; i<xe; i++) {	// pressure node; cell center
		double ren = user->ren;
		//if(levelset) ren = rho[k][j][i]/ mu[k][j][i];
		
		// from saved inflow file
		if(inletprofile==100 && k==1) {
			tmprt[k-1][j][i] = user->tmprt_plane[j][i];
		}
		else if ( user->bctype_tmprt[4]==5 && k==1 && nvert[k][j][i]<0.1) {	// inflow
			tmprt[k-1][j][i] =2.0*user->tmprts[4]-tmprt[k][j][i];
		}
			
		if ( user->bctype_tmprt[5]==4 && k==mz-2) {	// outflow
			tmprt[k+1][j][i] =tmprt[k][j][i];
		}	
		// NO energy flux
		if ( user->bctype_tmprt[0] == 10 && i==1 ) tmprt[k][j][i-1] = tmprt[k][j][i];
		if ( user->bctype_tmprt[1] == 10 && i==mx-2 ) tmprt[k][j][i+1] = tmprt[k][j][i];
		if ( user->bctype_tmprt[2] == 10 && j==1 ) tmprt[k][j-1][i] = tmprt[k][j][i];
		if ( user->bctype_tmprt[3] == 10 && j==my-2 ) tmprt[k][j+1][i] = tmprt[k][j][i];
		if ( user->bctype_tmprt[4] == 10 && k==1 ) tmprt[k-1][j][i] = tmprt[k][j][i];
		if ( user->bctype_tmprt[5] == 10 && k==mz-2 ) tmprt[k+1][j][i] = tmprt[k][j][i];

		// fixed temperature 

                if (user->bctype_tmprt[0] == 1 && i==1    ) {
			tmprt[k][j][i-1] = 2.0*user->tmprts[0]-tmprt[k][j][i];
		}

                if (user->bctype_tmprt[1] == 1 && i==mx-2 ) {
			tmprt[k][j][i+1] = 2.0*user->tmprts[1]-tmprt[k][j][i];
		}

                if (user->bctype_tmprt[2] == 1 && j==1    ) {
			tmprt[k][j-1][i] = 2.0*user->tmprts[2]-tmprt[k][j][i];
			//tmprt[k][j][i] = user->tmprts[2];
			//printf("tmprt BC jmin %le Tmprt(j=0)=%le Tmprt(j=1)=%le  \n", user->tmprts[2],  tmprt[k][j-1][i],  tmprt[k][j][i]);
		}

                if (user->bctype_tmprt[3] == 1 && j==my-2 ) {
			tmprt[k][j+1][i] = 2.0*user->tmprts[3]-tmprt[k][j][i];
			//tmprt[k][j][i] = user->tmprts[3];
			//printf("tmprt BC jmax %le Tmprt(j=my-1)=%le Tmprt(j=my-2)=%le  \n", user->tmprts[3],  tmprt[k][j+1][i],  tmprt[k][j][i]);
		}

                if (user->bctype_tmprt[4] == 1 && k==1    ) {
			tmprt[k-1][j][i] = 2.0*user->tmprts[4]-tmprt[k][j][i];
		}

                if (user->bctype_tmprt[5] == 1 && k==mz-2 ) {
			tmprt[k+1][j][i] = 2.0*user->tmprts[5]-tmprt[k][j][i];
		}


		// Uniform heat flux

		double D_m = 1./ (user->ren * user->Pr);
		if ( user->bctype_tmprt[0] == 2 && i==1 ) {
	                double area = sqrt( csi[k][j][i].x*csi[k][j][i].x + csi[k][j][i].y*csi[k][j][i].y + csi[k][j][i].z*csi[k][j][i].z );
        	        double hx = 1.0/aj[k][j][i]/area;

			tmprt[k][j][i-1] = -user->tmprts[0]/D_m + tmprt[k][j][i];
		}

                if ( user->bctype_tmprt[1] == 2 && i==mx-2 ) {
                        double area = sqrt( csi[k][j][i].x*csi[k][j][i].x + csi[k][j][i].y*csi[k][j][i].y + csi[k][j][i].z*csi[k][j][i].z );
                        double hx = 1.0/aj[k][j][i]/area;
                        tmprt[k][j][i+1] =  user->tmprts[1]/D_m + tmprt[k][j][i];
                }


                if ( user->bctype_tmprt[2] == 2 && j==1 ) {
                        double area = sqrt( eta[k][j][i].x*eta[k][j][i].x + eta[k][j][i].y*eta[k][j][i].y + eta[k][j][i].z*eta[k][j][i].z );
                        double hx = 1.0/aj[k][j][i]/area;
                        tmprt[k][j-1][i] = -user->tmprts[2]/D_m + tmprt[k][j][i];
                }

                if ( abs(user->bctype_tmprt[3]) == 2 && j==my-2 ) {
                        double area = sqrt( eta[k][j][i].x*eta[k][j][i].x + eta[k][j][i].y*eta[k][j][i].y + eta[k][j][i].z*eta[k][j][i].z );
                        double hx = 1.0/aj[k][j][i]/area;
                        tmprt[k][j+1][i] = user->tmprts[3]/D_m + tmprt[k][j][i];
                }


                if ( abs(user->bctype_tmprt[4]) == 2 && k==1 ) {
                        double area = sqrt( zet[k][j][i].x*zet[k][j][i].x + zet[k][j][i].y*zet[k][j][i].y + zet[k][j][i].z*zet[k][j][i].z );
                        double hx = 1.0/aj[k][j][i]/area;
                        tmprt[k-1][j][i] = -user->tmprts[4]/D_m + tmprt[k][j][i];
                }

                if ( abs(user->bctype_tmprt[5]) == 2 && k==mz-2 ) {
                        double area = sqrt( zet[k][j][i].x*zet[k][j][i].x + zet[k][j][i].y*zet[k][j][i].y + zet[k][j][i].z*zet[k][j][i].z );
                        double hx = 1.0/aj[k][j][i]/area;
                        tmprt[k+1][j][i] = user->tmprts[5]/D_m + tmprt[k][j][i];
                }




/*			
		
		// wall function
		if( nvert[k][j][i]<1.1 && ( ( (user->bctype[0]==-1 || user->bctype[0]==-2) && i==1) || ( (user->bctype[1]==-1 || user->bctype[1]==-2) &&  i==mx-2) ) && (j!=0 && j!=my-1 && k!=0 && k!=mz-1) ) {
		  
		}
		
		if( nvert[k][j][i]<1.1 && ( ( (user->bctype[2]==-1 || user->bctype[2]==-2) && j==1) || ( (user->bctype[3]==-1 || user->bctype[3]==-2) &&  j==my-2) ) && (i!=0 && i!=mx-1 && k!=0 && k!=mz-1)) {
		  
		
		if ( nvert[k][j][i] > 1.1 ) K_Omega[k][j][i].x = K_Omega[k][j][i].y = 0.;
*/

	}
	DMDAVecRestoreArray(user->da, user->lTmprt, &tmprt);

        DMDALocalToLocalBegin(user->da, user->lTmprt, INSERT_VALUES, user->lTmprt);
        DMDALocalToLocalEnd(user->da, user->lTmprt, INSERT_VALUES, user->lTmprt);

        DMDAVecGetArray(user->da, user->lTmprt, &tmprt);
	
	for (k=zs; k<ze; k++)
	for (j=ys; j<ye; j++)
	for (i=xs; i<xe; i++) {
		int flag=0, a=i, b=j, c=k;
		
		if(i_periodic && i==0) a=mx-2, flag=1;
		else if(i_periodic && i==mx-1) a=1, flag=1;
		
		if(j_periodic && j==0) b=my-2, flag=1;
		else if(j_periodic && j==my-1) b=1, flag=1;
		
		if(k_periodic && k==0) c=mz-2, flag=1;
		else if(k_periodic && k==mz-1) c=1, flag=1;
		
		
		if(ii_periodic && i==0) a=-2, flag=1;
		else if(ii_periodic && i==mx-1) a=mx+1, flag=1;
		
		if(jj_periodic && j==0) b=-2, flag=1;
		else if(jj_periodic && j==my-1) b=my+1, flag=1;
		
		if(kk_periodic && k==0) c=-2, flag=1;
		else if(kk_periodic && k==mz-1) c=mz+1, flag=1;
	
		if(flag) tmprt[k][j][i] = tmprt[c][b][a];
	}

/*	
	if(immersed)
	for(ibi=0; ibi<NumberOfBodies; ibi++)
	{
		extern IBMNodes *ibm_ptr;
		IBMNodes *ibm = &ibm_ptr[ibi];
		
		IBMListNode *current;
		current = user->ibmlist[ibi].head;
		while (current) {
			IBMInfo *ibminfo = &current->ibm_intp;
			//int ni = ibminfo->cell;
			current = current->next;
			double sb = ibminfo->d_s, sc = sb + ibminfo->d_i;
			int ip1 = ibminfo->i1, jp1 = ibminfo->j1, kp1 = ibminfo->k1;
			int ip2 = ibminfo->i2, jp2 = ibminfo->j2, kp2 = ibminfo->k2;
			int ip3 = ibminfo->i3, jp3 = ibminfo->j3, kp3 = ibminfo->k3;
			i = ibminfo->ni, j= ibminfo->nj, k = ibminfo->nk;
			double sk1  = ibminfo->cr1, sk2 = ibminfo->cr2, sk3 = ibminfo->cr3;
			double Kc = (K_Omega[kp1][jp1][ip1].x * sk1 + K_Omega[kp2][jp2][ip2].x * sk2 + K_Omega[kp3][jp3][ip3].x * sk3);
			
			double ren = user->ren;
			if(levelset) ren = rho[k][j][i]/ mu[k][j][i];
				
			if(wallfunction && ti>0) {
				double utau = ustar[k][j][i];
				const double yplus_min = 0.25;
				sb = PetscMax( sb, 1./ren * yplus_min / utau );	// prevent the case sb=0
				double K, Ob;//, Om;
				
				// Wilcox pp.108-109
				K = utau*utau/sqrt(0.09);
				Ob = utau/sqrt(0.09)/(kappa*sb);
				
				K_Omega[k][j][i].x = K;
				K_Omega[k][j][i].y = Ob;
			}
			else {
				const double yplus_min = 0.25;
				
				double utau = ustar[k][j][i];
				
				K_Omega[k][j][i].x = Kc * sb / sc;
				sb = PetscMax( sb, 1./ren * yplus_min / utau );	// prevent the case sb=0
				K_Omega[k][j][i].y = wall_omega(ren, sb);	
				
				if ( K_Omega[k][j][i].x < 0 ) K_Omega[k][j][i].x = utau*utau/sqrt(0.09);
			}
			if(user->bctype[4]==5 && k==1) K_Omega[k][j][i]=K_Omega[k-1][j][i];
		};
	}
*/

	
	DMDAVecRestoreArray(fda, user->lCsi, &csi);
	DMDAVecRestoreArray(fda, user->lEta, &eta);
	DMDAVecRestoreArray(fda, user->lZet, &zet);
	DMDAVecRestoreArray(fda, user->lUcat, &ucat);
	
	DMDAVecRestoreArray(da, user->lNvert, &nvert);
	DMDAVecRestoreArray(da, user->lAj, &aj);
	DMDAVecRestoreArray(da, user->lUstar, &ustar);
	
	DMDAVecRestoreArray(user->da, user->lTmprt, &tmprt);
	
	if(levelset) {
		DMDAVecRestoreArray(da, user->lDensity, &rho);
		DMDAVecRestoreArray(da, user->lMu, &mu);
	}
	
	DMDALocalToLocalBegin(user->da, user->lTmprt, INSERT_VALUES, user->lTmprt);
	DMDALocalToLocalEnd(user->da, user->lTmprt, INSERT_VALUES, user->lTmprt);
};

PetscErrorCode FormFunction_Tmprt(SNES snes, Vec Tmprt, Vec Rhs, void *ptr)
{
	UserCtx *user = (UserCtx*)ptr;

	DMDALocalInfo	info;
	PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
	PetscInt	mx, my, mz; // Dimensions in three directions
	PetscInt	i, j, k;
	PetscInt	lxs, lxe, lys, lye, lzs, lze;
	
	PetscReal	***nvert;
	PetscReal ***rhs;

	Vec Rhs_wm;	

	DMDAGetLocalInfo(user->da, &info);
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
	

	DMGlobalToLocalBegin(user->da, Tmprt, INSERT_VALUES, user->lTmprt);
	DMGlobalToLocalEnd(user->da, Tmprt, INSERT_VALUES, user->lTmprt);

	
	Tmprt_BC(user);
	RHS_Tmprt(user, Rhs);
	
	if (temperature && Shear_wm && (imin_wmtmprt != 0 || imax_wmtmprt != 0 || jmin_wmtmprt != 0 || jmax_wmtmprt != 0 || IB_wmtmprt != 0)) {
        	VecDuplicate(Rhs, &Rhs_wm);
		VecSet(Rhs_wm,0.0);
		Formfunction_wm_tmprt(user, Rhs_wm, 1.0);
		VecAXPY(Rhs, 1, Rhs_wm);
	}
	
//	VecSet(Rhs, 0.0);
//        TECIOOut_rhs1(user, Tmprt);
//        int aaa;
//        cout << "here \n";
//        cin >> aaa;

        double coeff = time_coeff();

	//PetscPrintf(PETSC_COMM_WORLD, "coeff %le \n", coeff);	
		

        if( coeff>0.9 && coeff<1.1 ) {
         	VecAXPY(Rhs, -1/user->dt, Tmprt);
	        VecAXPY(Rhs, 1/user->dt, user->Tmprt_o);

        }
        else/* if( coeff > 1.4 && coeff < 1.6 )*/ {
	        VecAXPY(Rhs, -1.5/user->dt, Tmprt);
        	VecAXPY(Rhs, 2./user->dt, user->Tmprt_o);
	        VecAXPY(Rhs, -0.5/user->dt, user->Tmprt_rm1);

        }


	/*	
	PetscReal ***tmprt_o, ***tmprt;

	DMDAVecGetArray(user->da, user->Tmprt, &tmprt);
	DMDAVecGetArray(user->da, user->Tmprt_o, &tmprt_o);
		
	for (k=zs; k<ze; k++)
	for (j=ys; j<ye; j++)
	for (i=xs; i<xe; i++) {

		if (i==10 && k==10 && j==1) printf("at i=10 j==1 k=10 tmprt=%le, tmprt_o=%le \n", tmprt[k][j][i], tmprt_o[k][j][i]);
		if (i==10 && k==10 && j==50) printf("at i=10 j==50 k=10 tmprt=%le, tmprt_o=%le \n", tmprt[k][j][i], tmprt_o[k][j][i]);

	}


	DMDAVecRestoreArray(user->da, user->Tmprt, &tmprt);
	DMDAVecRestoreArray(user->da, user->Tmprt_o, &tmprt_o);
	*/


	DMDAVecGetArray(user->da, user->lNvert, &nvert);
	DMDAVecGetArray(user->da, Rhs, &rhs);
	
	for (k=zs; k<ze; k++)
	for (j=ys; j<ye; j++)
	for (i=xs; i<xe; i++) {

		//if (j==0) printf("rhs at j==0 %le \n", rhs[k][j][i]);
		//if (i==10 && k==10 && j==1) printf("FOrm Func rhs at i=10 j==1 k=10 %le \n", rhs[k][j][i]);
		//if (i==10 && k==10 && j==50) printf("FOrm Func rhs at i=10 j==50 k=10 %le \n", rhs[k][j][i]);
		//if (j==2) printf("rhs at j==2 %le \n", rhs[k][j][i]);

		if(i==0 || i==mx-1 || j==0 || j==my-1 || k==0 || k==mz-1 || nvert[k][j][i]>0.1) {
			rhs[k][j][i] = 0;
		}

/*		
		// couette
		if ( user->bctype[3] == 12 && j==my-2 ) rhs[k][j][i].y = 0;
	
		//wall_omega
		if( i<=1 && user->bctype[0]==1 ) rhs[k][j][i].y = 0;
		if( i>=mx-2 && user->bctype[1]==1 ) rhs[k][j][i].y = 0;
		
		if( j<=1 && user->bctype[2]==1 ) rhs[k][j][i].y = 0;
		if( j>=my-2 && user->bctype[3]==1 ) rhs[k][j][i].y = 0;
		
		if( k==1 && user->bctype[4]==1 ) rhs[k][j][i].y = 0;
		if( k==mz-2 && user->bctype[5]==1 ) rhs[k][j][i].y = 0;
		
		// wall function k, omega
		if ( 	( i==1 && user->bctype[0] == -1 ) || ( i==mx-2 && user->bctype[1] == -1 ) ||
			( j==1 && user->bctype[2] == -1 ) || ( j==my-2 && user->bctype[3] == -1 ) ||
			( k==1 && user->bctype[4] == -1 ) || ( k==mz-2 && user->bctype[5] == -1 ) ||
			( i==1 && user->bctype[0] == -2 ) || ( i==mx-2 && user->bctype[1] == -2 ) ||
                        ( j==1 && user->bctype[2] == -2 ) || ( j==my-2 && user->bctype[3] == -2 ) ||
                        ( k==1 && user->bctype[4] == -2 ) || ( k==mz-2 && user->bctype[5] == -2 )
			) {
			rhs[k][j][i].x = 0;
			rhs[k][j][i].y = 0;
		}
*/

	}
	DMDAVecRestoreArray(user->da, user->lNvert, &nvert);
	DMDAVecRestoreArray(user->da, Rhs, &rhs);

//        TECIOOut_rhs1(user, Rhs);
//        cout << "here 2 \n";
//        cin >> aaa;


	if (temperature && Shear_wm && (imin_wmtmprt != 0 || imax_wmtmprt != 0 || jmin_wmtmprt != 0 || jmax_wmtmprt != 0 || IB_wmtmprt != 0)) VecDestroy(&Rhs_wm);

	return(0);
}









int snes_Tmprt_created=0;
Vec r_Tmprt;
Mat J_Tmprt;
SNES snes_Tmprt_eq;

void Solve_Tmprt(UserCtx *user)
{
	
	KSP ksp;
	PC pc;
	
	
	double norm;
	
	int bi=0;
	double tol=1.e-6;//1.e-6
	
	if(!snes_Tmprt_created) {
		snes_Tmprt_created=1;
		
		VecDuplicate(user[bi].Tmprt, &r_Tmprt);
		SNESCreate(PETSC_COMM_WORLD,&snes_Tmprt_eq);
		SNESSetFunction(snes_Tmprt_eq,r_Tmprt,FormFunction_Tmprt,(void *)&user[bi]);
		MatCreateSNESMF(snes_Tmprt_eq, &J_Tmprt);
		SNESSetJacobian(snes_Tmprt_eq,J_Tmprt,J_Tmprt,MatMFFDComputeJacobian,(void *)&user[bi]);
		SNESSetType(snes_Tmprt_eq, SNESTR);			//SNESTR,SNESLS	
		SNESSetMaxLinearSolveFailures(snes_Tmprt_eq,10000);
		SNESSetMaxNonlinearStepFailures(snes_Tmprt_eq,10000);		
		SNESKSPSetUseEW(snes_Tmprt_eq, PETSC_TRUE);
		SNESKSPSetParametersEW(snes_Tmprt_eq,3,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
		SNESSetTolerances(snes_Tmprt_eq,PETSC_DEFAULT,tol,PETSC_DEFAULT,50,50000);
			
		SNESGetKSP(snes_Tmprt_eq, &ksp);
		KSPSetType(ksp, KSPGMRES);
		//KSPGMRESSetPreAllocateVectors(ksp);
		
		KSPGetPC(ksp,&pc);
		PCSetType(pc,PCNONE);
		
		int maxits=20/*10000*/;	double rtol=tol, atol=PETSC_DEFAULT, dtol=PETSC_DEFAULT;
		KSPSetTolerances(ksp,rtol,atol,dtol,maxits);
	}
	
	extern PetscErrorCode MySNESMonitor(SNES snes,PetscInt n,PetscReal rnorm,void *dummy);
	SNESMonitorSet(snes_Tmprt_eq,MySNESMonitor,PETSC_NULL,PETSC_NULL);
	
	PetscPrintf(PETSC_COMM_WORLD, "\nSolving Temperature...\n");
	
	SNESSolve(snes_Tmprt_eq, PETSC_NULL, user[bi].Tmprt);
	
	SNESGetFunctionNorm(snes_Tmprt_eq, &norm);
	PetscPrintf(PETSC_COMM_WORLD, "\nTemperature SNES residual norm=%.5e\n\n", norm);

	DMDALocalInfo	info;
	PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
	PetscInt	mx, my, mz; // Dimensions in three directions
	PetscInt	i, j, k;
	PetscInt	lxs, lxe, lys, lye, lzs, lze;

	DMDAGetLocalInfo(user->da, &info);
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
	
	DMGlobalToLocalBegin(user->da, user->Tmprt, INSERT_VALUES, user->lTmprt);
	DMGlobalToLocalEnd(user->da, user->Tmprt, INSERT_VALUES, user->lTmprt);

	PetscReal ***tmprt, ***ltmprt;
	
	DMDAVecGetArray(user->da, user->Tmprt, &tmprt);
	DMDAVecGetArray(user->da, user->lTmprt, &ltmprt);

	if(periodic)
	for (k=zs; k<ze; k++)
	for (j=ys; j<ye; j++)
	for (i=xs; i<xe; i++) {
		int flag=0, a=i, b=j, c=k;

		if(i_periodic && i==0) a=mx-2, flag=1;
		else if(i_periodic && i==mx-1) a=1, flag=1;
		
		if(j_periodic && j==0) b=my-2, flag=1;
		else if(j_periodic && j==my-1) b=1, flag=1;
		
		if(k_periodic && k==0) c=mz-2, flag=1;
		else if(k_periodic && k==mz-1) c=1, flag=1;
		
		
		if(ii_periodic && i==0) a=-2, flag=1;
		else if(ii_periodic && i==mx-1) a=mx+1, flag=1;
		
		if(jj_periodic && j==0) b=-2, flag=1;
		else if(jj_periodic && j==my-1) b=my+1, flag=1;
		
		if(kk_periodic && k==0) c=-2, flag=1;
		else if(kk_periodic && k==mz-1) c=mz+1, flag=1;
		
		/*
                if(k_periodic && k==0) c=mz-2, flag=1;
                else if(k_periodic && k==mz-1) c=1, flag=1;
                if(j_periodic && j==0) b=my-2, flag=1;
                else if(j_periodic && j==my-1) b=1, flag=1;
                if(i_periodic && i==0) a=mx-2, flag=1;
                else if(i_periodic && i==mx-1) a=1, flag=1;
		
		if(kk_periodic && k==0) c=-2, flag=1;
                else if(kk_periodic && k==mz-1) c=mz+1, flag=1;
		*/
		if(flag) tmprt[k][j][i] = ltmprt[c][b][a];
	}
	DMDAVecRestoreArray(user->da, user->Tmprt, &tmprt);
	DMDAVecRestoreArray(user->da, user->lTmprt, &ltmprt);

	DMGlobalToLocalBegin(user->da, user->Tmprt, INSERT_VALUES, user->lTmprt);
	DMGlobalToLocalEnd(user->da, user->Tmprt, INSERT_VALUES, user->lTmprt);
	Tmprt_BC(user);
	DMLocalToGlobalBegin(user->da, user->lTmprt, INSERT_VALUES, user->Tmprt);
	DMLocalToGlobalEnd(user->da, user->lTmprt, INSERT_VALUES, user->Tmprt);


	/*
	char fname[80];
	sprintf(fname,"Tmprt");
	TECIOOut_rhs_da(user, user->Tmprt, fname);
	*/


//		TECIOOut_rhs1(user,  user->Tmprt_o);
//	TECIOOut_rhs1(user,  user->Tmprt);

};


PetscErrorCode TECIOOut_rhs1(UserCtx *user, Vec Rhs)	
{


	PetscInt IMax, JMax, KMax;
        PetscReal x, y, z, F_wm_x, F_wm_y, F_wm_z;

        int rank;
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	char filen[80];
	sprintf(filen, "rhs%03d.plt", rank);


        FILE * pFile;
        pFile = fopen(filen, "w");


        fprintf(pFile, " TITLE = \"rhs\" \n");
        fprintf(pFile, " VARIABLES = \"X\", \"Y\", \"Z\", \"rhs\" \n" );

        DM              da = user->da, fda = user->fda;
        DMDALocalInfo     info;

        DMDAGetLocalInfo(da, &info);

	PetscInt	xs = info.xs, xe = info.xs + info.xm;
	PetscInt  	ys = info.ys, ye = info.ys + info.ym;
	PetscInt	zs = info.zs, ze = info.zs + info.zm;
	PetscInt	mx = info.mx, my = info.my, mz = info.mz;
    
	PetscInt	lxs, lys, lzs, lxe, lye, lze;
	PetscInt	i, j, k;
	Cmpnts		***coor;
	PetscReal	***rhs;
	Vec		Coor;
		
	IMax = mx-2;
	JMax = my-2;
	KMax = mz-2;

	fprintf(pFile, "ZONE I=%d, J=%d, K=%d, DATAPACKING=POINT \n", xe-xs-1, ye-ys-1, ze-zs-1 );
		
	DMDAGetCoordinates(da, &Coor);
	DMDAVecGetArray(fda, Coor, &coor);

	DMDAVecGetArray(da, Rhs,  &rhs);


  	for (k=zs; k<ze-1; k++){
        for (j=ys; j<ye-1; j++){
        for (i=xs; i<xe-1; i++){
		x = coor[k][j][i].x; y = coor[k][j][i].y; z = coor[k][j][i].z;
 		fprintf(pFile, "%le %le %le %le  \n", x, y, z, rhs[k][j][i] );

	}
	}
	}

	fclose(pFile);

        DMDAVecRestoreArray(fda, Coor, &coor);		

	DMDAVecRestoreArray(da, Rhs,  &rhs);

	return 0;
}


PetscErrorCode TE_Output(UserCtx *user)
{
	DMDALocalInfo info = user->info;
	PetscInt	xs = info.xs, xe = info.xs + info.xm;
	PetscInt  	ys = info.ys, ye = info.ys + info.ym;
	PetscInt	zs = info.zs, ze = info.zs + info.zm;
	PetscInt	mx = info.mx, my = info.my, mz = info.mz;

	PetscInt i, j, k;
	PetscInt	lxs, lys, lzs, lxe, lye, lze;
  
	PetscReal	***tmprt;
	PetscReal ***aj;
	
	lxs = xs; lxe = xe;
	lys = ys; lye = ye;
	lzs = zs; lze = ze;

	if (xs==0) lxs = xs+1;
	if (ys==0) lys = ys+1;
	if (zs==0) lzs = zs+1;

	if (xe==mx) lxe = xe-1;
	if (ye==my) lye = ye-1;
	if (ze==mz) lze = ze-1;
  
	DMDAVecGetArray(user->da, user->lTmprt, &tmprt);
	DMDAVecGetArray(user->da, user->lAj, &aj);
		
	double local_sum=0, sum=0;
	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
		local_sum += 0.5 * tmprt[k][j][i] * tmprt[k][j][i] / aj[k][j][i];
	}
//	PetscGlobalSum(&local_sum, &sum, PETSC_COMM_WORLD);
        MPI_Allreduce (&local_sum, &sum, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
	
	DMDAVecRestoreArray(user->da, user->lTmprt, &tmprt);
	DMDAVecRestoreArray(user->da, user->lAj, &aj);
	
	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	if (!rank) {	
		char filen[80];
		sprintf(filen, "%s/Tmperature_Square.dat", path);
		FILE *f = fopen(filen, "a");
		PetscFPrintf(PETSC_COMM_WORLD, f, "%d\t%.7e\n", ti, sum);
		fclose(f);
	}
	
	return 0;
}

// xyang add fluctuations 
PetscErrorCode Add_fluc_tmprt(UserCtx *user)
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
	DMDAVecGetArray(da, user->Tmprt, &p);
  
    
	srand( time(NULL)) ;	// seokkoo
	for (i = 0; i < (rand() % 3000); i++) (rand() % 3000);	//seokkoo
  
	if( initial_perturbation) {
		PetscPrintf(PETSC_COMM_WORLD, "\nGenerating initial perturbation\n");
		for(k=lzs; k<lze; k++) {	// for all CPU
			for (j=lys; j<lye; j++)
			for (i=lxs; i<lxe; i++) {
				if (nvert[k][j][i]+nvert[k+1][j][i] < 0.1) {
				  int n1, n2, n3;
					double F;
					
					F  = 0.1; // 100% 
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
                 			p[k][j][i] *= ( 1.0 + ((double)n3)/10000. * F );     // uin * (1+-0.xx)
//					p[k][j][i] += p[k][j][i]*F*((double)n3)/10000.;
				}
			}
		}
	}
	
	if(levelset) DMDAVecRestoreArray(da, user->lLevelset, &level);
	DMDAVecRestoreArray(da, user->lNvert, &nvert);
	DMDAVecRestoreArray(da, user->Tmprt, &p);
	
	DMGlobalToLocalBegin(da, user->Tmprt, INSERT_VALUES, user->lTmprt);
	DMGlobalToLocalEnd(da, user->Tmprt, INSERT_VALUES, user->lTmprt);
	
		
	return 0;
}
