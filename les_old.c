#include "variables.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

extern double integrate_hat(int np, double *val, double *w);
extern double integrate_testfilter_i(double val[3][3][3], double vol[3][3][3]);
extern double integrate_testfilter_j(double val[3][3][3], double vol[3][3][3]);
extern double integrate_testfilter_k(double val[3][3][3], double vol[3][3][3]);


double integrate_gridfilter(double val[3][3][3], double vol[3][3][3]);
extern void Calculate_Covariant_tensor(double g[3][3], double G[3][3]);
extern void Calculate_Covariant_metrics(double g[3][3], double G[3][3]);

extern int mixed;
extern int les, ti, tistart;
extern PetscBool rstart_flg;

const double les_eps=1.e-4;// osl 1.e-7
const double wall_cs=0.001;

void get_weight ( int i, int j, int k, int mx, int my, int mz, PetscReal ***aj, PetscReal ***nvert, PetscReal nv, double weight[3][3][3])
{
		for(int p=-1; p<=1; p++)
		for(int q=-1; q<=1; q++)
		for(int r=-1; r<=1; r++) {
			int R=r+1, Q=q+1, P=p+1;
			int K=k+r, J=j+q, I=i+p;

			if( nvert[K][J][I]>nv ) weight[R][Q][P]=0;
			else weight[R][Q][P] = 1./aj[K][J][I];
			/*
			if( i_periodic ) {
				if( I==0 ) I=mx-2;
				else if( I==mx-1 ) I=1;
			}
			else if( ii_periodic ) {
				if( I==0 ) I=-2;
				else if( I==mx-1 ) I=mx+1;
			}
			else if( I==0 || I==mx-1) weight[R][Q][P]=0;
			
			if( j_periodic ) {
				if( J==0 ) J=my-2;
				else if( J==my-1 ) J=1;
			}
			else if( jj_periodic ) {
				if( J==0 ) J=-2;
				else if( J==my-1 ) J=my+1;
			}
			else if(J==0 || J==my-1) weight[R][Q][P]=0;
			
			if( k_periodic ) {
				if( K==0 ) K=mz-2;
				else if( K==mz-1 ) K=1;
			}
			else if( kk_periodic ) {
				if( K==0 ) K=-2;
				else if( K==mz-1 ) K=mz+1;
			}
			else if( K==0 || K==mz-1 ) weight[R][Q][P]=0;
			*/
		}
}

void Compute_Smagorinsky_Constant_1(UserCtx *user, Vec Ucont, Vec Ucat)
{
	if(ti<2 && tistart==0 && !rstart_flg){
		VecSet(user->lCs, 0.0);
		return;
	}
	
	if(les==1) {
		VecSet(user->lCs, 0.01);
		return;
	}
	
	Vec Csi = user->lCsi, Eta = user->lEta, Zet = user->lZet;

	Cmpnts	***ucont, ***ucat;
	Cmpnts	***csi, ***eta, ***zet;
	
	PetscReal	***nvert, ***Cs;//, ***lnu_t;

	DM		da = user->da, fda = user->fda;
	DMDALocalInfo	info;
	int	xs, xe, ys, ye, zs, ze; // Local grid information
	int	mx, my, mz; // Dimensions in three directions
	int	i, j, k;
  
	PetscReal	***aj;
	Cmpnts ***Ax, ***Ay, ***Az, ***ucat_f, ***cent;//, ;
	double ***LM, ***MM;//, ***NM;
	PetscReal ***Sabs;//, ***Cs, ***lCs_o;//, ***nu_t;

	int	lxs, lxe, lys, lye, lzs, lze;

	PetscReal	ajc;

	PetscReal	dudc, dude, dudz, dvdc, dvde, dvdz, dwdc, dwde, dwdz;
	
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

	Vec lSx;		// Sx :	 Sxx, Sxy, Sxz
	Vec lSy;		// Sy :	 Syx, Syy, Syz
	Vec lSz;		// Sz :	 Szx, Szy, Szz
	Vec lS;
	Vec lUcat_f;
	
	DMDAVecGetArray(fda, Ucont, &ucont);
	DMDAVecGetArray(fda, Ucat,  &ucat);
	DMDAVecGetArray(fda, Csi, &csi);
	DMDAVecGetArray(fda, Eta, &eta);
	DMDAVecGetArray(fda, Zet, &zet);
	DMDAVecGetArray(da, user->lNvert, &nvert);
	DMDAVecGetArray(da, user->lAj, &aj); 
	
	DMDAVecGetArray(da, user->lCs, &Cs);
	//DMDAVecGetArray(da, user->lNu_t, &lnu_t);
	
	DMDAVecGetArray(fda, user->lCent, &cent);

	Cmpnts	***icsi, ***ieta, ***izet;
	Cmpnts	***jcsi, ***jeta, ***jzet;
	Cmpnts	***kcsi, ***keta, ***kzet;
	PetscReal	***iaj, ***jaj, ***kaj;

	DMDAVecGetArray(da, user->lIAj, &iaj);  
	DMDAVecGetArray(da, user->lJAj, &jaj);  
	DMDAVecGetArray(da, user->lKAj, &kaj);  
	
	DMDAVecGetArray(fda, user->lICsi, &icsi);
	DMDAVecGetArray(fda, user->lIEta, &ieta);
	DMDAVecGetArray(fda, user->lIZet, &izet);

	DMDAVecGetArray(fda, user->lJCsi, &jcsi);
	DMDAVecGetArray(fda, user->lJEta, &jeta);
	DMDAVecGetArray(fda, user->lJZet, &jzet);

	DMDAVecGetArray(fda, user->lKCsi, &kcsi);
	DMDAVecGetArray(fda, user->lKEta, &keta);
	DMDAVecGetArray(fda, user->lKZet, &kzet);
	//
  
	VecDuplicate(user->lP, &user->lLM);
	VecDuplicate(user->lP, &user->lMM);
	//VecDuplicate(user->lP, &user->lNM);
	
	VecSet(user->lLM, 0);
	VecSet(user->lMM, 0);
	//VecSet(user->lNM, 0);
	
	VecDuplicate(user->lUcont, &lSx);
	VecDuplicate(user->lUcont, &lSy);
	VecDuplicate(user->lUcont, &lSz);
	VecDuplicate(user->lNvert, &lS);
	VecDuplicate(user->lUcont, &lUcat_f);

	VecSet(lSx, 0);  
	VecSet(lSy, 0);  
	VecSet(lSz, 0);	
	VecSet(lS, 0);
	
	DMDAVecGetArray(da, user->lLM, &LM);//
	DMDAVecGetArray(da, user->lMM, &MM);//
	//DMDAVecGetArray(da, user->lNM, &NM);//
  	
	DMDAVecGetArray(fda, lSx, &Ax);//
	DMDAVecGetArray(fda, lSy, &Ay);//
	DMDAVecGetArray(fda, lSz, &Az);//
	DMDAVecGetArray(da, lS, &Sabs);//
	DMDAVecGetArray(fda, lUcat_f, &ucat_f);//
  	
	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
		if( nvert[k][j][i]>1.1) continue;
		
		ajc = aj[k][j][i];
		double csi0 = csi[k][j][i].x, csi1 = csi[k][j][i].y, csi2 = csi[k][j][i].z;
		double eta0 = eta[k][j][i].x, eta1 = eta[k][j][i].y, eta2 = eta[k][j][i].z;
		double zet0 = zet[k][j][i].x, zet1 = zet[k][j][i].y, zet2 = zet[k][j][i].z;

				
		double dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz;
		double du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz;
		
		Compute_du_center (i, j, k, mx, my, mz, ucat, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
		Compute_du_dxyz (	csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz,
							&du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
		
		double Sxx = 0.5*( du_dx + du_dx ), Sxy = 0.5*(du_dy + dv_dx), Sxz = 0.5*(du_dz + dw_dx);
		double Syx = Sxy, Syy = 0.5*(dv_dy + dv_dy),	Syz = 0.5*(dv_dz + dw_dy);
		double Szx = Sxz, Szy=Syz, Szz = 0.5*(dw_dz + dw_dz);
	
		Sabs[k][j][i] = sqrt( 2.0*( Sxx*Sxx + Sxy*Sxy + Sxz*Sxz + Syx*Syx + Syy*Syy + Syz*Syz + Szx*Szx + Szy*Szy + Szz*Szz ) );
		
		Ax[k][j][i].x = du_dx;	Ax[k][j][i].y = du_dy;	Ax[k][j][i].z = du_dz;
		Ay[k][j][i].x = dv_dx;	Ay[k][j][i].y = dv_dy;	Ay[k][j][i].z = dv_dz;
		Az[k][j][i].x = dw_dx;	Az[k][j][i].y = dw_dy;	Az[k][j][i].z = dw_dz;
		
		double weight[3][3][3];
		double u[3][3][3];
		double v[3][3][3];
		double w[3][3][3];
		get_weight (i, j, k, mx, my, mz, aj, nvert, 1.1, weight);
		
		for(int p=-1; p<=1; p++)
		for(int q=-1; q<=1; q++)
		for(int r=-1; r<=1; r++) {
			int R=r+1, Q=q+1, P=p+1;
			int K=k+r, J=j+q, I=i+p;
			u[R][Q][P] = ucat[K][J][I].x;
			v[R][Q][P] = ucat[K][J][I].y;
			w[R][Q][P] = ucat[K][J][I].z;
		}
		
		ucat_f[k][j][i].x = integrate_testfilter(u, weight);
		ucat_f[k][j][i].y = integrate_testfilter(v, weight);
		ucat_f[k][j][i].z = integrate_testfilter(w, weight);
	}
	
	DMDAVecRestoreArray(fda, lSx, &Ax);//
	DMDAVecRestoreArray(fda, lSy, &Ay);//
	DMDAVecRestoreArray(fda, lSz, &Az);//
	DMDAVecRestoreArray(da, lS, &Sabs);//
	DMDAVecRestoreArray(fda, lUcat_f, &ucat_f);//
		
	DMDALocalToLocalBegin(fda, lSx, INSERT_VALUES, lSx);
	DMDALocalToLocalEnd(fda, lSx, INSERT_VALUES, lSx);
 
	DMDALocalToLocalBegin(fda, lSy, INSERT_VALUES, lSy);
	DMDALocalToLocalEnd(fda, lSy, INSERT_VALUES, lSy);
 
	DMDALocalToLocalBegin(fda, lSz, INSERT_VALUES, lSz);
	DMDALocalToLocalEnd(fda, lSz, INSERT_VALUES, lSz);
 
	DMDALocalToLocalBegin(da, lS, INSERT_VALUES, lS);
	DMDALocalToLocalEnd(da, lS, INSERT_VALUES, lS);
	
	DMDALocalToLocalBegin(fda, lUcat_f, INSERT_VALUES, lUcat_f);
	DMDALocalToLocalEnd(fda, lUcat_f, INSERT_VALUES, lUcat_f);
	
	DMDAVecGetArray(fda, lSx, &Ax);//
	DMDAVecGetArray(fda, lSy, &Ay);//
	DMDAVecGetArray(fda, lSz, &Az);//
	DMDAVecGetArray(da, lS, &Sabs);//
	DMDAVecGetArray(fda, lUcat_f, &ucat_f);//
	
	for (k=zs; k<ze; k++)
	for (j=ys; j<ye; j++)
	for (i=xs; i<xe; i++) {
		int c=k, b=j, a=i, flag=0;
		
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
			Sabs[k][j][i] = Sabs[c][b][a];
			Ax[k][j][i] = Ax[c][b][a];
			Ay[k][j][i] = Ay[c][b][a];
			Az[k][j][i] = Az[c][b][a];
			ucat_f[k][j][i] = ucat_f[c][b][a];
		}
	}
	
 	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
		if(nvert[k][j][i]>1.1) {
			LM[k][j][i]=MM[k][j][i]=0;//NM[k][j][i]=0;
			continue;
		}
		ajc = aj[k][j][i];
	
		double csi0 = csi[k][j][i].x, csi1 = csi[k][j][i].y, csi2 = csi[k][j][i].z;
		double eta0 = eta[k][j][i].x, eta1 = eta[k][j][i].y, eta2 = eta[k][j][i].z;
		double zet0 = zet[k][j][i].x, zet1 = zet[k][j][i].y, zet2 = zet[k][j][i].z;
	
		int a, b;
		double Lij[3][3], Sij_hat[3][3], SSij_hat[3][3], Mij[3][3], Nij[3][3], Nij_cat[3][3];
		double Lij_cat[3][3], Mij_cat[3][3];
		
		double filter, test_filter;
		double S[3][3][3], S_hat;
		double u[3][3][3], v[3][3][3], w[3][3][3];
		
		double U[3][3][3], V[3][3][3], W[3][3][3];
		double Uu[3][3][3], Uv[3][3][3], Uw[3][3][3];
		double Vu[3][3][3], Vv[3][3][3], Vw[3][3][3];
		double Wu[3][3][3], Wv[3][3][3], Ww[3][3][3];
		
		double uu[3][3][3], uv[3][3][3], uw[3][3][3];
		double vv[3][3][3], vw[3][3][3], ww[3][3][3];
		double S11[3][3][3], S12[3][3][3], S13[3][3][3], S21[3][3][3], S22[3][3][3], S23[3][3][3], S31[3][3][3], S32[3][3][3], S33[3][3][3];
		double SS11[3][3][3], SS12[3][3][3], SS13[3][3][3], SS21[3][3][3], SS22[3][3][3], SS23[3][3][3], SS31[3][3][3], SS32[3][3][3], SS33[3][3][3];
		
		double dU_dx[3][3][3], dU_dy[3][3][3], dU_dz[3][3][3];
		double dV_dx[3][3][3], dV_dy[3][3][3], dV_dz[3][3][3];
		double dW_dx[3][3][3], dW_dy[3][3][3], dW_dz[3][3][3];
		double T11[3][3][3], T12[3][3][3], T13[3][3][3], T21[3][3][3], T22[3][3][3], T23[3][3][3], T31[3][3][3], T32[3][3][3], T33[3][3][3];	// Clark model term
		
		double weight[3][3][3];
		
		double dx,dy,dz;
		Calculate_dxdydz(aj[k][j][i], csi[k][j][i], eta[k][j][i], zet[k][j][i], &dx, &dy, &dz);
		double dx2=dx*dx, dy2=dy*dy, dz2=dz*dz;
		

		get_weight (i, j, k, mx, my, mz, aj, nvert, 1.1, weight);
		
		int p,q,r;
		for(p=-1; p<=1; p++)
		for(q=-1; q<=1; q++)
		for(r=-1; r<=1; r++) {
			int R=r+1, Q=q+1, P=p+1;
			int K=k+r, J=j+q, I=i+p;

			u[R][Q][P] = v[R][Q][P] = w[R][Q][P] = 0;
			uu[R][Q][P] = uv[R][Q][P] = uw[R][Q][P] = 0;
			vv[R][Q][P] = vw[R][Q][P] = ww[R][Q][P] = 0;
			
			S11[R][Q][P] = S12[R][Q][P] = S13[R][Q][P] = 0;
			S21[R][Q][P] = S22[R][Q][P] = S23[R][Q][P] = 0;
			S31[R][Q][P] = S32[R][Q][P] = S33[R][Q][P] = 0;
						
			S[R][Q][P] = 0;
				
			SS11[R][Q][P] = SS12[R][Q][P] = SS13[R][Q][P] = 0;
			SS21[R][Q][P] = SS22[R][Q][P] = SS23[R][Q][P] = 0;
			SS31[R][Q][P] = SS32[R][Q][P] = SS33[R][Q][P] = 0;
						
			u[R][Q][P] = ucat[K][J][I].x;
			v[R][Q][P] = ucat[K][J][I].y;
			w[R][Q][P] = ucat[K][J][I].z;
			
			// metric tensors are also test-filtered : Big difference
			
			/*if(I==0)*/ U[R][Q][P] = u[R][Q][P]*csi[K][J][I].x + v[R][Q][P]*csi[K][J][I].y + w[R][Q][P]*csi[K][J][I].z;
			//else U[R][Q][P] = 0.5 * ( ucont[K][J][I].x + ucont[K][J][I-1].x );
			
			/*if(J==0)*/ V[R][Q][P] = u[R][Q][P]*eta[K][J][I].x + v[R][Q][P]*eta[K][J][I].y + w[R][Q][P]*eta[K][J][I].z;
			//else V[R][Q][P] = 0.5 * ( ucont[K][J][I].y + ucont[K][J-1][I].y );
			
			/*if(K==0)*/ W[R][Q][P] = u[R][Q][P]*zet[K][J][I].x + v[R][Q][P]*zet[K][J][I].y + w[R][Q][P]*zet[K][J][I].z;
			//else W[R][Q][P] = 0.5 * ( ucont[K][J][I].z + ucont[K-1][J][I].z );
			
			Uu[R][Q][P] = U[R][Q][P]*u[R][Q][P];
			Uv[R][Q][P] = U[R][Q][P]*v[R][Q][P];
			Uw[R][Q][P] = U[R][Q][P]*w[R][Q][P];
			
			Vu[R][Q][P] = V[R][Q][P]*u[R][Q][P];
			Vv[R][Q][P] = V[R][Q][P]*v[R][Q][P];
			Vw[R][Q][P] = V[R][Q][P]*w[R][Q][P];
			
			Wu[R][Q][P] = W[R][Q][P]*u[R][Q][P];
			Wv[R][Q][P] = W[R][Q][P]*v[R][Q][P];
			Ww[R][Q][P] = W[R][Q][P]*w[R][Q][P];
			
			uu[R][Q][P] = u[R][Q][P]*u[R][Q][P];	
			uv[R][Q][P] = u[R][Q][P]*v[R][Q][P];	
			uw[R][Q][P] = u[R][Q][P]*w[R][Q][P];
			vv[R][Q][P] = v[R][Q][P]*v[R][Q][P];	
			vw[R][Q][P] = v[R][Q][P]*w[R][Q][P];	
			ww[R][Q][P] = w[R][Q][P]*w[R][Q][P];
			
			const double du_dx = Ax[K][J][I].x, du_dy = Ax[K][J][I].y, du_dz = Ax[K][J][I].z;
			const double dv_dx = Ay[K][J][I].x, dv_dy = Ay[K][J][I].y, dv_dz = Ay[K][J][I].z;
			const double dw_dx = Az[K][J][I].x, dw_dy = Az[K][J][I].y, dw_dz = Az[K][J][I].z;
		
			const double Sxx = 0.5*( du_dx + du_dx ), Sxy = 0.5*(du_dy + dv_dx), Sxz = 0.5*(du_dz + dw_dx);
			const double Syx = Sxy, Syy = 0.5*(dv_dy + dv_dy),	Syz = 0.5*(dv_dz + dw_dy);
			const double Szx = Sxz, Szy=Syz, Szz = 0.5*(dw_dz + dw_dz);
			
			dU_dx[R][Q][P] = du_dx, dU_dy[R][Q][P] = du_dy, dU_dz[R][Q][P] = du_dz;
			dV_dx[R][Q][P] = dv_dx, dV_dy[R][Q][P] = dv_dy, dV_dz[R][Q][P] = dv_dz;
			dW_dx[R][Q][P] = dw_dx, dW_dy[R][Q][P] = dw_dy, dW_dz[R][Q][P] = dw_dz;
			
			T11[R][Q][P] = ( du_dx * du_dx * dx2 + du_dy * du_dy * dy2 + du_dz * du_dz * dz2 ) / 12.;
			T12[R][Q][P] = ( du_dx * dv_dx * dx2 + du_dy * dv_dy * dy2 + du_dz * dv_dz * dz2 ) / 12.;
			T13[R][Q][P] = ( du_dx * dw_dx * dx2 + du_dy * dw_dy * dy2 + du_dz * dw_dz * dz2 ) / 12.;
			T21[R][Q][P] = T12[R][Q][P];
			T22[R][Q][P] = ( dv_dx * dv_dx * dx2 + dv_dy * dv_dy * dy2 + dv_dz * dv_dz * dz2 ) / 12.;
			T23[R][Q][P] = ( dv_dx * dw_dx * dx2 + dv_dy * dw_dy * dy2 + dv_dz * dw_dz * dz2 ) / 12.;
			T31[R][Q][P] = T13[R][Q][P];
			T32[R][Q][P] = T23[R][Q][P];
			T33[R][Q][P] = ( dw_dx * dw_dx * dx2 + dw_dy * dw_dy * dy2 + dw_dz * dw_dz * dz2 ) / 12.;
		
			S11[R][Q][P] = Sxx, S12[R][Q][P] = Sxy, S13[R][Q][P] = Sxz;
			S21[R][Q][P] = Syx, S22[R][Q][P] = Syy, S23[R][Q][P] = Syz;
			S31[R][Q][P] = Szx, S32[R][Q][P] = Szy, S33[R][Q][P] = Szz;
			
			S[R][Q][P] = Sabs[K][J][I];
			
			SS11[R][Q][P] = S11[R][Q][P]*S[R][Q][P], SS12[R][Q][P] = S12[R][Q][P]*S[R][Q][P], SS13[R][Q][P] = S13[R][Q][P]*S[R][Q][P];
			SS21[R][Q][P] = S21[R][Q][P]*S[R][Q][P], SS22[R][Q][P] = S22[R][Q][P]*S[R][Q][P], SS23[R][Q][P] = S23[R][Q][P]*S[R][Q][P];
			SS31[R][Q][P] = S31[R][Q][P]*S[R][Q][P], SS32[R][Q][P] = S32[R][Q][P]*S[R][Q][P], SS33[R][Q][P] = S33[R][Q][P]*S[R][Q][P];
		}
		
		double sum_weight=0;
		double coef[3][3][3]={
			0.125, 0.250, 0.125, 
			0.250, 0.500, 0.250, 
			0.125, 0.250, 0.125, 
				
			0.250, 0.500, 0.250,
			0.500, 1.000, 0.500,
			0.250, 0.500, 0.250,
				
			0.125, 0.250, 0.125, 
			0.250, 0.500, 0.250,
			0.125, 0.250, 0.125
		};
		
		for(p=-1; p<=1; p++)
		for(q=-1; q<=1; q++)
		for(r=-1; r<=1; r++) {
			sum_weight += weight[r+1][q+1][p+1] * coef[r+1][q+1][p+1];
		}
		
		filter = pow( 1./aj[k][j][i], 1./3. );
		
		if(testfilter_ik) test_filter = pow(5.0, 1./3.) * filter;
		else {
			//test_filter = 2.0 * filter;
			test_filter = pow( sum_weight, 1./3. );
		}
		/*		
		double _U=integrate_testfilter(U, weight);
		double _V=integrate_testfilter(V, weight);
		double _W=integrate_testfilter(W, weight);
		
		double _u=integrate_testfilter(u, weight);
		double _v=integrate_testfilter(v, weight);
		double _w=integrate_testfilter(w, weight);
		*/
		double _u=ucat_f[k][j][i].x;
		double _v=ucat_f[k][j][i].y;
		double _w=ucat_f[k][j][i].z;
		double _U = _u * csi[k][j][i].x + _v * csi[k][j][i].y + _w * csi[k][j][i].z;
		double _V = _u * eta[k][j][i].x + _v * eta[k][j][i].y + _w * eta[k][j][i].z;
		double _W = _u * zet[k][j][i].x + _v * zet[k][j][i].y + _w * zet[k][j][i].z;
		
		double _dudc, _dvdc, _dwdc, _dude, _dvde, _dwde, _dudz, _dvdz, _dwdz;
		double _du_dx, _du_dy, _du_dz, _dv_dx, _dv_dy, _dv_dz, _dw_dx, _dw_dy, _dw_dz;
		Compute_du_center ( i, j, k, mx, my, mz, ucat_f, nvert, &_dudc, &_dvdc, &_dwdc, &_dude, &_dvde, &_dwde, &_dudz, &_dvdz, &_dwdz);        // derivative of ucat_f
		Compute_du_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, _dudc, _dvdc, _dwdc, _dude, _dvde, _dwde, _dudz, _dvdz, _dwdz, &_du_dx, &_dv_dx, &_dw_dx, &_du_dy, &_dv_dy, &_dw_dy, &_du_dz, &_dv_dz, &_dw_dz );
		double _Sxx = 0.5*( _du_dx + _du_dx ), _Sxy = 0.5*(_du_dy + _dv_dx), _Sxz = 0.5*(_du_dz + _dw_dx);
		double _Syx = _Sxy, _Syy = 0.5*(_dv_dy + _dv_dy), _Syz = 0.5*(_dv_dz + _dw_dy);
		double _Szx = _Sxz, _Szy=_Syz, _Szz = 0.5*(_dw_dz + _dw_dz);
		Sij_hat[0][0] = _Sxx, Sij_hat[0][1] = _Sxy, Sij_hat[0][2] = _Sxz;
		Sij_hat[1][0] = _Syx, Sij_hat[1][1] = _Syy, Sij_hat[1][2] = _Syz;
		Sij_hat[2][0] = _Szx, Sij_hat[2][1] = _Szy, Sij_hat[2][2] = _Szz;
		S_hat = sqrt( 2.0 *( _Sxx*_Sxx + _Sxy*_Sxy + _Sxz*_Sxz + _Syx*_Syx + _Syy*_Syy + _Syz*_Syz + _Szx*_Szx + _Szy*_Szy + _Szz*_Szz ) );

		double _T11=integrate_testfilter(T11, weight);
		double _T12=integrate_testfilter(T12, weight);
		double _T13=integrate_testfilter(T13, weight);
		double _T21=integrate_testfilter(T21, weight);
		double _T22=integrate_testfilter(T22, weight);
		double _T23=integrate_testfilter(T23, weight);
		double _T31=integrate_testfilter(T31, weight);
		double _T32=integrate_testfilter(T32, weight);
		double _T33=integrate_testfilter(T33, weight);
		/*
		double _dU_dx=integrate_testfilter(dU_dx, weight);
		double _dV_dx=integrate_testfilter(dV_dx, weight);
		double _dW_dx=integrate_testfilter(dW_dx, weight);
		double _dU_dy=integrate_testfilter(dU_dy, weight);
		double _dV_dy=integrate_testfilter(dV_dy, weight);
		double _dW_dy=integrate_testfilter(dW_dy, weight);
		double _dU_dz=integrate_testfilter(dU_dz, weight);
		double _dV_dz=integrate_testfilter(dV_dz, weight);
		double _dW_dz=integrate_testfilter(dW_dz, weight);
		*/
		double _R11 = ( _du_dx * _du_dx * dx2 + _du_dy * _du_dy * dy2 + _du_dz * _du_dz * dz2 ) / 12.;
		double _R12 = ( _du_dx * _dv_dx * dx2 + _du_dy * _dv_dy * dy2 + _du_dz * _dv_dz * dz2 ) / 12.;
		double _R13 = ( _du_dx * _dw_dx * dx2 + _du_dy * _dw_dy * dy2 + _du_dz * _dw_dz * dz2 ) / 12.;
		double _R21 = _R12;
		double _R22 = ( _dv_dx * _dv_dx * dx2 + _dv_dy * _dv_dy * dy2 + _dv_dz * _dv_dz * dz2 ) / 12.;
		double _R23 = ( _dv_dx * _dw_dx * dx2 + _dv_dy * _dw_dy * dy2 + _dv_dz * _dw_dz * dz2 ) / 12.;
		double _R31 = _R13;
		double _R32 = _R23;
		double _R33 = ( _dw_dx * _dw_dx * dx2 + _dw_dy * _dw_dy * dy2 + _dw_dz * _dw_dz * dz2 ) / 12.;
		
		Lij[0][0] = integrate_testfilter(Uu, weight) - _U*_u;
		Lij[0][1] = integrate_testfilter(Uv, weight) - _U*_v;
		Lij[0][2] = integrate_testfilter(Uw, weight) - _U*_w;
		Lij[1][0] = integrate_testfilter(Vu, weight) - _V*_u;
		Lij[1][1] = integrate_testfilter(Vv, weight) - _V*_v;
		Lij[1][2] = integrate_testfilter(Vw, weight) - _V*_w;
		Lij[2][0] = integrate_testfilter(Wu, weight) - _W*_u;
		Lij[2][1] = integrate_testfilter(Wv, weight) - _W*_v;
		Lij[2][2] = integrate_testfilter(Ww, weight) - _W*_w;
		
		
		Nij_cat[0][0] = 4.0 * _R11  - _T11;
		Nij_cat[0][1] = 4.0 * _R12  - _T12;
		Nij_cat[0][2] = 4.0 * _R13  - _T13;
		Nij_cat[1][0] = 4.0 * _R21  - _T21;
		Nij_cat[1][1] = 4.0 * _R22  - _T22;
		Nij_cat[1][2] = 4.0 * _R23  - _T23;
		Nij_cat[2][0] = 4.0 * _R31  - _T31;
		Nij_cat[2][1] = 4.0 * _R32  - _T32;
		Nij_cat[2][2] = 4.0 * _R33  - _T33;
		
		Nij[0][0] = Nij_cat[0][0] * csi0 + Nij_cat[0][1] * csi1 + Nij_cat[0][2] * csi2;
		Nij[0][1] = Nij_cat[0][0] * eta0 + Nij_cat[0][1] * eta1 + Nij_cat[0][2] * eta2;
		Nij[0][2] = Nij_cat[0][0] * zet0 + Nij_cat[0][1] * zet1 + Nij_cat[0][2] * zet2;
		Nij[1][0] = Nij_cat[1][0] * csi0 + Nij_cat[1][1] * csi1 + Nij_cat[1][2] * csi2;
		Nij[1][1] = Nij_cat[1][0] * eta0 + Nij_cat[1][1] * eta1 + Nij_cat[1][2] * eta2;
		Nij[1][2] = Nij_cat[1][0] * zet0 + Nij_cat[1][1] * zet1 + Nij_cat[1][2] * zet2;
		Nij[2][0] = Nij_cat[2][0] * csi0 + Nij_cat[2][1] * csi1 + Nij_cat[2][2] * csi2;
		Nij[2][1] = Nij_cat[2][0] * eta0 + Nij_cat[2][1] * eta1 + Nij_cat[2][2] * eta2;
		Nij[2][2] = Nij_cat[2][0] * zet0 + Nij_cat[2][1] * zet1 + Nij_cat[2][2] * zet2;
		
		/*
		double _dudc, _dvdc, _dwdc, _dude, _dvde, _dwde, _dudz, _dvdz, _dwdz;
		double _du_dx, _du_dy, _du_dz, _dv_dx, _dv_dy, _dv_dz, _dw_dx, _dw_dy, _dw_dz;
		Compute_du_center ( i, j, k, mx, my, mz, ucat_f, nvert, &_dudc, &_dvdc, &_dwdc, &_dude, &_dvde, &_dwde, &_dudz, &_dvdz, &_dwdz);	// derivative of ucat_f
		Compute_du_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, _dudc, _dvdc, _dwdc, _dude, _dvde, _dwde, _dudz, _dvdz, _dwdz, &_du_dx, &_dv_dx, &_dw_dx, &_du_dy, &_dv_dy, &_dw_dy, &_du_dz, &_dv_dz, &_dw_dz );
		double _Sxx = 0.5*( _du_dx + _du_dx ), _Sxy = 0.5*(_du_dy + _dv_dx), _Sxz = 0.5*(_du_dz + _dw_dx);
		double _Syx = _Sxy, _Syy = 0.5*(_dv_dy + _dv_dy), _Syz = 0.5*(_dv_dz + _dw_dy);
		double _Szx = _Sxz, _Szy=_Syz, _Szz = 0.5*(_dw_dz + _dw_dz);
		Sij_hat[0][0] = _Sxx, Sij_hat[0][1] = _Sxy, Sij_hat[0][2] = _Sxz;
		Sij_hat[1][0] = _Syx, Sij_hat[1][1] = _Syy, Sij_hat[1][2] = _Syz;
		Sij_hat[2][0] = _Szx, Sij_hat[2][1] = _Szy, Sij_hat[2][2] = _Szz;
		S_hat = sqrt( 2.0 *( _Sxx*_Sxx + _Sxy*_Sxy + _Sxz*_Sxz + _Syx*_Syx + _Syy*_Syy + _Syz*_Syz + _Szx*_Szx + _Szy*_Szy + _Szz*_Szz ) );
		*/
		//
		
		/*
		Sij_hat[0][0] = integrate_testfilter(S11, weight);	
		Sij_hat[0][1] = integrate_testfilter(S12, weight);	
		Sij_hat[0][2] = integrate_testfilter(S13, weight);
		Sij_hat[1][0] = integrate_testfilter(S21, weight);	
		Sij_hat[1][1] = integrate_testfilter(S22, weight);	
		Sij_hat[1][2] = integrate_testfilter(S23, weight);
		Sij_hat[2][0] = integrate_testfilter(S31, weight);	
		Sij_hat[2][1] = integrate_testfilter(S32, weight);	
		Sij_hat[2][2] = integrate_testfilter(S33, weight);
		
		S_hat=0;
		for(a=0; a<3; a++)
		for(b=0; b<3; b++) {
			S_hat += pow( Sij_hat[a][b], 2. );
		}
		S_hat = sqrt ( 2 * S_hat );
		*/
		
		SSij_hat[0][0] = integrate_testfilter(SS11, weight);	
		SSij_hat[0][1] = integrate_testfilter(SS12, weight);	
		SSij_hat[0][2] = integrate_testfilter(SS13, weight);
		SSij_hat[1][0] = integrate_testfilter(SS21, weight);
		SSij_hat[1][1] = integrate_testfilter(SS22, weight);	
		SSij_hat[1][2] = integrate_testfilter(SS23, weight);
		SSij_hat[2][0] = integrate_testfilter(SS31, weight);	
		SSij_hat[2][1] = integrate_testfilter(SS32, weight);	
		SSij_hat[2][2] = integrate_testfilter(SS33, weight);
		
		
		
		//S_hat = integrate_testfilter(S, weight);
		
		
		double gg[3][3], ggc[3][3], G[3][3];
		double xcsi, xeta, xzet, ycsi, yeta, yzet, zcsi, zeta, zzet;

		gg[0][0]=csi0, gg[0][1]=csi1, gg[0][2]=csi2;
		gg[1][0]=eta0, gg[1][1]=eta1, gg[1][2]=eta2;
		gg[2][0]=zet0, gg[2][1]=zet1, gg[2][2]=zet2;
		Calculate_Covariant_metrics(gg, ggc);
		xcsi=ggc[0][0], xeta=ggc[0][1], xzet=ggc[0][2];
		ycsi=ggc[1][0], yeta=ggc[1][1], yzet=ggc[1][2];
		zcsi=ggc[2][0], zeta=ggc[2][1], zzet=ggc[2][2];
		G[0][0] = xcsi * xcsi + ycsi * ycsi + zcsi * zcsi;
		G[1][1] = xeta * xeta + yeta * yeta + zeta * zeta;
		G[2][2] = xzet * xzet + yzet * yzet + zzet * zzet;
		G[0][1] = G[1][0] = xeta * xcsi + yeta * ycsi + zeta * zcsi;
		G[0][2] = G[2][0] = xzet * xcsi + yzet * ycsi + zzet * zcsi;
		G[1][2] = G[2][1] = xeta * xzet + yeta * yzet + zeta * zzet;
		
		for(a=0; a<3; a++)
		for(b=0; b<3; b++) {
		  	Mij_cat[a][b] = -pow( test_filter, 2. ) * S_hat * Sij_hat[a][b] + pow( filter, 2. ) * SSij_hat[a][b];
			/*
			Mij_cat_i[a][b] = -4. * S_hat_i * Sij_hat_i[a][b] + SSij_hat_i[a][b];
			Mij_cat_j[a][b] = -4. * S_hat_j * Sij_hat_j[a][b] + SSij_hat_j[a][b];
			Mij_cat_k[a][b] = -4. * S_hat_k * Sij_hat_k[a][b] + SSij_hat_k[a][b];*/
		}
		
		Mij[0][0] = Mij_cat[0][0] * csi0 + Mij_cat[0][1] * csi1 + Mij_cat[0][2] * csi2;
		Mij[0][1] = Mij_cat[0][0] * eta0 + Mij_cat[0][1] * eta1 + Mij_cat[0][2] * eta2;
		Mij[0][2] = Mij_cat[0][0] * zet0 + Mij_cat[0][1] * zet1 + Mij_cat[0][2] * zet2;
		Mij[1][0] = Mij_cat[1][0] * csi0 + Mij_cat[1][1] * csi1 + Mij_cat[1][2] * csi2;
		Mij[1][1] = Mij_cat[1][0] * eta0 + Mij_cat[1][1] * eta1 + Mij_cat[1][2] * eta2;
		Mij[1][2] = Mij_cat[1][0] * zet0 + Mij_cat[1][1] * zet1 + Mij_cat[1][2] * zet2;
		Mij[2][0] = Mij_cat[2][0] * csi0 + Mij_cat[2][1] * csi1 + Mij_cat[2][2] * csi2;
		Mij[2][1] = Mij_cat[2][0] * eta0 + Mij_cat[2][1] * eta1 + Mij_cat[2][2] * eta2;
		Mij[2][2] = Mij_cat[2][0] * zet0 + Mij_cat[2][1] * zet1 + Mij_cat[2][2] * zet2;
					
		double num=0, num1=0, denom=0;
		int m, n, l;
		
		/*
			g11 ~ csi*csi ~ dx^4
			G11 ~ dx^-4
		*/
	
	
		for(q=0; q<3; q++)
		for(a=0; a<3; a++)
		for(b=0; b<3; b++) {
			num += Lij[a][b] * Mij[a][q] * G[b][q];
			if(clark) num -= Nij[a][b] * Mij[a][q] * G[b][q];
		}
		
		for(m=0; m<3; m++)
		for(n=0; n<3; n++)
		for(l=0; l<3; l++) {
			denom += Mij[n][m] * Mij[n][l] * G[m][l];
		}
	
		//printf("%f %f\n", num, denom);
		
		LM[k][j][i] = num;
		MM[k][j][i] = denom;
    	}
	
	DMDAVecRestoreArray(da, user->lLM, &LM);//
	DMDAVecRestoreArray(da, user->lMM, &MM);//
	//DMDAVecRestoreArray(da, user->lNM, &NM);//
	
	DMDALocalToLocalBegin(da, user->lLM, INSERT_VALUES, user->lLM);
	DMDALocalToLocalEnd(da, user->lLM, INSERT_VALUES, user->lLM);
	DMDALocalToLocalBegin(da, user->lMM, INSERT_VALUES, user->lMM);
	DMDALocalToLocalEnd(da, user->lMM, INSERT_VALUES, user->lMM);
	//DMDALocalToLocalBegin(da, user->lNM, INSERT_VALUES, user->lNM);
	//DMDALocalToLocalEnd(da, user->lNM, INSERT_VALUES, user->lNM);
	  
	DMDAVecGetArray(da, user->lLM, &LM);//
	DMDAVecGetArray(da, user->lMM, &MM);//
	//DMDAVecGetArray(da, user->lNM, &NM);//
	
	for (k=zs; k<ze; k++)
	for (j=ys; j<ye; j++)
	for (i=xs; i<xe; i++) {
		int c=k, b=j, a=i, flag=0;
		
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
			LM[k][j][i] = LM[c][b][a];
			MM[k][j][i] = MM[c][b][a];
		}
	}
	
	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
		if(nvert[k][j][i]>1.1) {
			Cs[k][j][i] = 0;
			continue;
		}
				
		double weight[3][3][3];
		double LM0[3][3][3], MM0[3][3][3];
		int a, b, c;
			
		for(a=-1; a<=1; a++)
		for(b=-1; b<=1; b++)
		for(c=-1; c<=1; c++) {
			int R=c+1, Q=b+1, P=a+1;
			int K=k+c, J=j+b, I=i+a;
			
			weight[R][Q][P] = 1./aj[K][J][I];
			
			if( nvert[K][J][I]>1.1 ) weight[R][Q][P]=0;
			
			if( i_periodic ) {
				if( I==0 ) I=mx-2;
				else if( I==mx-1 ) I=1;
			}
			else if( ii_periodic ) {
				if( I==0 ) I=-2;
				else if( I==mx-1 ) I=mx+1;
			}
			else if( I==0 || I==mx-1) weight[R][Q][P]=0;
			
			if( j_periodic ) {
				if( J==0 ) J=my-2;
				else if( J==my-1 ) J=1;
			}
			else if( jj_periodic ) {
				if( J==0 ) J=-2;
				else if( J==my-1 ) J=my+1;
			}
			else if( J==0 || j==my-1) weight[R][Q][P]=0;
			
			if( k_periodic ) {
				if( K==0 ) K=mz-2;
				else if( K==mz-1 ) K=1;
			}
			else if( kk_periodic ) {
				if( K==0 ) K=-2;
				else if( K==mz-1 ) K=mz+1;
			}
			else if( K==0 || K==mz-1) weight[R][Q][P]=0;
			
			LM0[R][Q][P] = LM[K][J][I];
			MM0[R][Q][P] = MM[K][J][I];
		}
			
		double C=0;
			
		double LM_avg, MM_avg;//, NM_avg;

		if ( i_homo_filter || j_homo_filter || k_homo_filter || les==3 ) {
			LM_avg = LM[k][j][i];
			MM_avg = MM[k][j][i];
		}
		else {
			LM_avg = integrate_testfilter(LM0,weight);
			MM_avg = integrate_testfilter(MM0,weight);
		}
		
		C = 0.5 * LM_avg / (MM_avg + les_eps );
		
		if ( les==3 ) {
			if(ti<100 && tistart==0 && !rstart_flg) { }
			else {
				C = (1.0 - 0.001) * Cs[k][j][i] + 0.001 * C;
			}
		}
		
		if(les==1) Cs[k][j][i] = 0.01;
		else Cs[k][j][i] = PetscMax(C, 0);
	}
	
	if( les==3 ) {}
	else if ( i_homo_filter && k_homo_filter ) {
		std::vector<int> count, total_count;
		std::vector<double> J_LM(my), J_MM(my), LM_tmp(my), MM_tmp(my);
		
		count.resize(my);
		total_count.resize(my);
		
		for(j=0; j<my; j++) {
			LM_tmp[j] = 0;
			MM_tmp[j] = 0;
			count[j] = total_count[j] = 0;
		}
		
		for (k=lzs; k<lze; k++)
		for (j=lys; j<lye; j++)
		for (i=lxs; i<lxe; i++) {
			if( nvert[k][j][i]<0.1 ) {
				LM_tmp[j] += LM[k][j][i];
				MM_tmp[j] += MM[k][j][i];
				count[j] ++;
			}
		}
		
		MPI_Allreduce( &LM_tmp[0], &J_LM[0], my, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
		MPI_Allreduce( &MM_tmp[0], &J_MM[0], my, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
		MPI_Allreduce( &count[0], &total_count[0], my, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);	
		
		for(j=0; j<my; j++) {
			if( total_count[j]>0) {
				J_LM[j] /= (double) (total_count[j]);
				J_MM[j] /= (double) (total_count[j]);
			}
		}
		
		for (j=lys; j<lye; j++)
		for (k=lzs; k<lze; k++)
		for (i=lxs; i<lxe; i++) {
			Cs[k][j][i] = 0.5 * J_LM[j] / ( J_MM[j]+les_eps);
		}
	}
	else if (i_homo_filter || j_homo_filter || k_homo_filter) {
		int plane_size;
		
		if(i_homo_filter) plane_size = my*mz;
		else if(j_homo_filter) plane_size = mx*mz;
		else if(k_homo_filter) plane_size = mx*my;
		
		std::vector<int> count(plane_size), total_count(plane_size);
		std::vector<double> J_LM(plane_size), J_MM(plane_size), LM_tmp(plane_size), MM_tmp(plane_size);
		int pos;
		
		pos=0;
		
		std::fill( LM_tmp.begin(), LM_tmp.end(), 0.);
		std::fill( MM_tmp.begin(), MM_tmp.end(), 0.);
		std::fill( count.begin(), count.end(), 0);
		std::fill( total_count.begin(), total_count.end(), 0);
		
		for(pos=0; pos<plane_size; pos++) {
			LM_tmp[pos] = 0;
			MM_tmp[pos] = 0;
			count[pos] = total_count[pos] = 0;
		}
		
		pos=0;
		
		if(i_homo_filter)  {
			for(k=0; k<mz; k++)
			for(j=0; j<my; j++) {
				if( k>=lzs && k<lze && j>=lys && j<lye) {
					for (i=lxs; i<lxe; i++) {
						if( nvert[k][j][i]<0.1 ) {
							LM_tmp[pos] += LM[k][j][i];
							MM_tmp[pos] += MM[k][j][i];
							count[pos] ++;
						}
					}
				}
				pos++;
			}
		}
		else if(j_homo_filter)  {
			for(k=0; k<mz; k++)
			for(i=0; i<mx; i++) {
				if( i>=lxs && i<lxe && k>=lzs && k<lze) {
					for (j=lys; j<lye; j++) {
						if( nvert[k][j][i]<0.1 ) {
							LM_tmp[pos] += LM[k][j][i];
							MM_tmp[pos] += MM[k][j][i];
							count[pos] ++;
						}
					}
				}
				pos++;
			}
		}
		else if(k_homo_filter)  {
			for(j=0; j<my; j++)
			for(i=0; i<mx; i++) {
				if( i>=lxs && i<lxe && j>=lys && j<lye) {
					for (k=lzs; k<lze; k++) {
						if( nvert[k][j][i]<0.1 ) {
							LM_tmp[pos] += LM[k][j][i];
							MM_tmp[pos] += MM[k][j][i];
							count[pos] ++;
						}
					}
				}
				pos++;
			}
		}
		
		MPI_Allreduce( &LM_tmp[0], &J_LM[0], plane_size, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
		MPI_Allreduce( &MM_tmp[0], &J_MM[0], plane_size, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
		MPI_Allreduce( &count[0], &total_count[0], plane_size, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);	
		
		pos=0;
		
		for(pos=0; pos<plane_size; pos++) {
			if( total_count[pos]>0) {
				double N = (double) (total_count[pos]);
				J_LM[pos] /= N;
				J_MM[pos] /= N;
			}
		}
		
		pos=0;
		if(i_homo_filter)  {
			for(k=0; k<mz; k++)
			for(j=0; j<my; j++) {
				if( k>=lzs && k<lze && j>=lys && j<lye) {
					for (i=lxs; i<lxe; i++) {
						if( nvert[k][j][i]<0.1 ) {
							Cs[k][j][i] = 0.5 * J_LM[pos] / ( J_MM[pos] + les_eps );
						}
					}
				}
				pos++;
			}
		}
		else if(j_homo_filter)  {
			for(k=0; k<mz; k++)
			for(i=0; i<mx; i++) {
				if( i>=lxs && i<lxe && k>=lzs && k<lze) {
					for (j=lys; j<lye; j++) {
						if( nvert[k][j][i]<0.1 ) {
							Cs[k][j][i] = 0.5 * J_LM[pos] / ( J_MM[pos] + les_eps );
						}
					}
				}
				pos++;
			}
		}
		else if(k_homo_filter)  {
			for(j=0; j<my; j++)
			for(i=0; i<mx; i++) {
				if( i>=lxs && i<lxe && j>=lys && j<lye) {
					for (k=lzs; k<lze; k++) {
						if( nvert[k][j][i]<0.1 ) {
							Cs[k][j][i] = 0.5 * J_LM[pos] / ( J_MM[pos] + les_eps );
						}
					}
				}
				pos++;
			}
		}
	}
	
	for (k=zs; k<ze; k++)
	for (j=ys; j<ye; j++)
	for (i=xs; i<xe; i++) {
		
		if(nvert[k][j][i]>1.1 || k==0 || k==mz-1 || j==0 || j==my-1 || i==0 || i==mx-1) {
			Cs[k][j][i] = 0;
		}
		else {
			if(nvert[k][j][i]>0.1 && nvert[k][j][i]<1.1) {
				Cs[k][j][i] = PetscMax(wall_cs, Cs[k][j][i]);	// stabilize at high Re, osl 0.005
			}
			Cs[k][j][i] = PetscMin( PetscMax ( Cs[k][j][i], 0 ), max_cs);
		}
	}
	
	DMDAVecRestoreArray(fda, lSx, &Ax);//
	DMDAVecRestoreArray(fda, lSy, &Ay);//
	DMDAVecRestoreArray(fda, lSz, &Az);//
	DMDAVecRestoreArray(da, lS, &Sabs);//
	DMDAVecRestoreArray(fda, lUcat_f, &ucat_f);//
	
	DMDAVecRestoreArray(da, user->lLM, &LM);//
	DMDAVecRestoreArray(da, user->lMM, &MM);//

	DMDAVecRestoreArray(fda, Ucont, &ucont);
	DMDAVecRestoreArray(fda, Ucat,  &ucat);
	DMDAVecRestoreArray(fda, Csi, &csi);
	DMDAVecRestoreArray(fda, Eta, &eta);
	DMDAVecRestoreArray(fda, Zet, &zet);
	DMDAVecRestoreArray(da, user->lNvert, &nvert);
	DMDAVecRestoreArray(da, user->lAj, &aj); 
	DMDAVecRestoreArray(da, user->lCs, &Cs);
	
	DMDAVecRestoreArray(fda, user->lCent, &cent);
	VecDestroy(&user->lLM);//
	VecDestroy(&user->lMM);//
  
	VecDestroy(&lSx);//
	VecDestroy(&lSy);//
	VecDestroy(&lSz);//
	VecDestroy(&lS);//
	VecDestroy(&lUcat_f);//

	DMDAVecRestoreArray(da, user->lIAj, &iaj);  
	DMDAVecRestoreArray(da, user->lJAj, &jaj);  
	DMDAVecRestoreArray(da, user->lKAj, &kaj);  
	
	DMDAVecRestoreArray(fda, user->lICsi, &icsi);
	DMDAVecRestoreArray(fda, user->lIEta, &ieta);
	DMDAVecRestoreArray(fda, user->lIZet, &izet);

	DMDAVecRestoreArray(fda, user->lJCsi, &jcsi);
	DMDAVecRestoreArray(fda, user->lJEta, &jeta);
	DMDAVecRestoreArray(fda, user->lJZet, &jzet);

	DMDAVecRestoreArray(fda, user->lKCsi, &kcsi);
	DMDAVecRestoreArray(fda, user->lKEta, &keta);
	DMDAVecRestoreArray(fda, user->lKZet, &kzet);
	
	DMDALocalToLocalBegin(da, user->lCs, INSERT_VALUES, user->lCs);
	DMDALocalToLocalEnd(da, user->lCs, INSERT_VALUES, user->lCs);
	
	DMDAVecGetArray(da, user->lCs, &Cs);
	for(k=zs; k<ze; k++)
	for(j=ys; j<ye; j++)
	for(i=xs; i<xe; i++) {
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


		if(flag) Cs[k][j][i] = Cs[c][b][a];
	}
	DMDAVecRestoreArray(da, user->lCs, &Cs);
	/*
	if ( ! i_homo_filter && ! j_homo_filter && ! k_homo_filter && les!=3 ) {
		Vec lCs_tmp;
		PetscReal ***Cs_tmp;
		
		VecDuplicate(user->lCs, &lCs_tmp);
		VecCopy(user->lCs, lCs_tmp);
		
		DMDAVecGetArray(da, user->lAj, &aj);
		DMDAVecGetArray(da, user->lNvert, &nvert);
		DMDAVecGetArray(da, user->lCs, &Cs);
		DMDAVecGetArray(da, lCs_tmp, &Cs_tmp);
		for(k=lzs; k<lze; k++)
		for(j=lys; j<lye; j++)
		for(i=lxs; i<lxe; i++) {
			if(nvert[k][j][i]<0.1) {
				double weight[3][3][3], C[3][3][3];
				int a, b, c;

				for(a=-1; a<=1; a++)
				for(b=-1; b<=1; b++)
				for(c=-1; c<=1; c++) {
					int R=c+1, Q=b+1, P=a+1;
					int K=k+c, J=j+b, I=i+a;
					
					weight[R][Q][P] = 1./aj[K][J][I];
					
					if( nvert[K][J][I]>1.1 ) weight[R][Q][P]=0;
					
					if( i_periodic ) {
						if( I==0 ) I=mx-2;
						else if( I==mx-1 ) I=1;
					}
					else if( ii_periodic ) {
						if( I==0 ) I=-2;
						else if( I==mx-1 ) I=mx+1;
					}
					else if( I==0 || I==mx-1) weight[R][Q][P]=0;
					
					if( j_periodic ) {
						if( J==0 ) J=my-2;
						else if( J==my-1 ) J=1;
					}
					else if( jj_periodic ) {
						if( J==0 ) J=-2;
						else if( J==my-1 ) J=my+1;
					}
					else if( J==0 || j==my-1) weight[R][Q][P]=0;
					
					if( k_periodic ) {
						if( K==0 ) K=mz-2;
						else if( K==mz-1 ) K=1;
					}
					else if( kk_periodic ) {
						if( K==0 ) K=-2;
						else if( K==mz-1 ) K=mz+1;
					}
					else if( K==0 || K==mz-1) weight[R][Q][P]=0;
					
					C[R][Q][P] = Cs_tmp[K][J][I];
				}
				Cs[k][j][i] = integrate_testfilter(C,weight);
			}
		}
		DMDAVecRestoreArray(da, user->lAj, &aj);
		DMDAVecRestoreArray(da, user->lNvert, &nvert);
		DMDAVecRestoreArray(da, user->lCs, &Cs);
		DMDAVecRestoreArray(da, lCs_tmp, &Cs_tmp);
		
		VecDestroy(&lCs_tmp);
	}
	*/
	double lmax_norm=0, max_norm;
	int p;
	
	if(testfilter_ik) PetscPrintf(PETSC_COMM_WORLD, "\nFilter type : Box filter homogeneous\n");
	else PetscPrintf(PETSC_COMM_WORLD, "Filter type : Box filter 3D\n");
	if(clark) PetscPrintf(PETSC_COMM_WORLD, "Clark model\n");
	
	VecMax(user->lCs, &p, &lmax_norm);
	GlobalMax_All(&lmax_norm, &max_norm, PETSC_COMM_WORLD);
	PetscPrintf(PETSC_COMM_WORLD, "Max Cs  = %e \n", max_norm);
	
};

void Compute_eddy_viscosity_LES(UserCtx *user)
{
	DM		da = user->da, fda = user->fda;
	DMDALocalInfo	info;
	int	xs, xe, ys, ye, zs, ze; // Local grid information
	int	mx, my, mz; // Dimensions in three directions
	int	i, j, k;
	
	PetscReal ***Cs, ***lnu_t, ***nvert, ***aj;
	Cmpnts ***csi, ***eta, ***zet, ***ucat;
	
	DMDAGetLocalInfo(da, &info);
	mx = info.mx, my = info.my, mz = info.mz;
	xs = info.xs, xe = xs + info.xm;
	ys = info.ys, ye = ys + info.ym;
	zs = info.zs, ze = zs + info.zm;

	int lxs = xs, lxe = xe;
	int lys = ys, lye = ye;
	int lzs = zs, lze = ze;

	if (xs==0) lxs = xs+1;
	if (ys==0) lys = ys+1;
	if (zs==0) lzs = zs+1;

	if (xe==mx) lxe = xe-1;
	if (ye==my) lye = ye-1;
	if (ze==mz) lze = ze-1;

	DMDAVecGetArray(fda, user->lUcat,  &ucat);
	DMDAVecGetArray(fda, user->lCsi, &csi);
	DMDAVecGetArray(fda, user->lEta, &eta);
	DMDAVecGetArray(fda, user->lZet, &zet);
	
	DMDAVecGetArray(da, user->lNvert, &nvert);
	DMDAVecGetArray(da, user->lAj, &aj);
	DMDAVecGetArray(da, user->lNu_t, &lnu_t);
	DMDAVecGetArray(da, user->lCs, &Cs);
	
	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
		if(nvert[k][j][i]>1.1) {
			lnu_t[k][j][i]=0;
			continue;
		}
		double ajc = aj[k][j][i];
		double csi0 = csi[k][j][i].x, csi1 = csi[k][j][i].y, csi2 = csi[k][j][i].z;
		double eta0 = eta[k][j][i].x, eta1 = eta[k][j][i].y, eta2 = eta[k][j][i].z;
		double zet0 = zet[k][j][i].x, zet1 = zet[k][j][i].y, zet2 = zet[k][j][i].z;
		double dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz;
		double du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz;
		Compute_du_center (i, j, k, mx, my, mz, ucat, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
		Compute_du_dxyz (	csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz,
							&du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
		
		double Sxx = 0.5*( du_dx + du_dx ), Sxy = 0.5*(du_dy + dv_dx), Sxz = 0.5*(du_dz + dw_dx);
		double Syx = Sxy, Syy = 0.5*(dv_dy + dv_dy),	Syz = 0.5*(dv_dz + dw_dy);
		double Szx = Sxz, Szy=Syz, Szz = 0.5*(dw_dz + dw_dz);
	
		double Sabs = sqrt( 2.0*( Sxx*Sxx + Sxy*Sxy + Sxz*Sxz + Syx*Syx + Syy*Syy + Syz*Syz + Szx*Szx + Szy*Szy + Szz*Szz ) );
		
		double filter  = pow( 1./aj[k][j][i],1./3.);
		lnu_t[k][j][i] = Cs[k][j][i] * pow ( filter, 2.0 ) * Sabs;
	}
	
	DMDAVecRestoreArray(fda, user->lUcat,  &ucat);
	DMDAVecRestoreArray(fda, user->lCsi, &csi);
	DMDAVecRestoreArray(fda, user->lEta, &eta);
	DMDAVecRestoreArray(fda, user->lZet, &zet);
	
	DMDAVecRestoreArray(da, user->lNvert, &nvert);
	DMDAVecRestoreArray(da, user->lAj, &aj);
	DMDAVecRestoreArray(da, user->lNu_t, &lnu_t);
	DMDAVecRestoreArray(da, user->lCs, &Cs);
	
	DMDALocalToLocalBegin(da, user->lNu_t, INSERT_VALUES, user->lNu_t);
	DMDALocalToLocalEnd(da, user->lNu_t, INSERT_VALUES, user->lNu_t);
	
	DMDAVecGetArray(da, user->lNu_t, &lnu_t);
	for(k=zs; k<ze; k++)
	for(j=ys; j<ye; j++)
	for(i=xs; i<xe; i++) {
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


		if(flag) {
			lnu_t[k][j][i] = lnu_t[c][b][a];
		}
	}
	DMDAVecRestoreArray(da, user->lNu_t, &lnu_t);
	
	double lmax_norm=0, max_norm;
	int p;
	
	VecMax(user->lNu_t, &p, &lmax_norm);
	GlobalMax_All(&lmax_norm, &max_norm, PETSC_COMM_WORLD);
	PetscPrintf(PETSC_COMM_WORLD, "Max Nu_t = %e\n", max_norm);
};
