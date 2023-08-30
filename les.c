#include "variables.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

// const double les_eps=1.e-4;// osl 1.e-7
// This is the limit for wall by SK
//const double wall_cs=0.001;
const double wall_cs=0.001;

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
  
	PetscReal	***aj, ***rho, ***level;
	Cmpnts ***Ax, ***Ay, ***Az, ***cent;
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
	
	if(levelset) {
		DMDAVecGetArray(da, user->lDensity, &rho);
		DMDAVecGetArray(da, user->lLevelset, &level);
	}
	
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
	}
	
	DMDAVecRestoreArray(fda, lSx, &Ax);//
	DMDAVecRestoreArray(fda, lSy, &Ay);//
	DMDAVecRestoreArray(fda, lSz, &Az);//
	DMDAVecRestoreArray(da, lS, &Sabs);//
		
	DMDALocalToLocalBegin(fda, lSx, INSERT_VALUES, lSx);
	DMDALocalToLocalEnd(fda, lSx, INSERT_VALUES, lSx);
 
	DMDALocalToLocalBegin(fda, lSy, INSERT_VALUES, lSy);
	DMDALocalToLocalEnd(fda, lSy, INSERT_VALUES, lSy);
 
	DMDALocalToLocalBegin(fda, lSz, INSERT_VALUES, lSz);
	DMDALocalToLocalEnd(fda, lSz, INSERT_VALUES, lSz);
 
	DMDALocalToLocalBegin(da, lS, INSERT_VALUES, lS);
	DMDALocalToLocalEnd(da, lS, INSERT_VALUES, lS);

	DMDAVecGetArray(fda, lSx, &Ax);//
	DMDAVecGetArray(fda, lSy, &Ay);//
	DMDAVecGetArray(fda, lSz, &Az);//
	DMDAVecGetArray(da, lS, &Sabs);//
	
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
		double Lij[3][3], Sij_hat[3][3], SSij_hat[3][3], rhoSSij_hat[3][3], Mij[3][3], Nij[3][3], Nij_cat[3][3];
		double Lij_cat[3][3], Mij_cat[3][3];
		
		double filter, test_filter;
		double S[3][3][3], S_hat;
		
		double u[3][3][3], v[3][3][3], w[3][3][3], d[3][3][3];
		double U[3][3][3], V[3][3][3], W[3][3][3];
		double Uu[3][3][3], Uv[3][3][3], Uw[3][3][3];
		double Vu[3][3][3], Vv[3][3][3], Vw[3][3][3];
		double Wu[3][3][3], Wv[3][3][3], Ww[3][3][3];
		
		double S11[3][3][3], S12[3][3][3], S13[3][3][3], S21[3][3][3], S22[3][3][3], S23[3][3][3], S31[3][3][3], S32[3][3][3], S33[3][3][3];
		double SS11[3][3][3], SS12[3][3][3], SS13[3][3][3], SS21[3][3][3], SS22[3][3][3], SS23[3][3][3], SS31[3][3][3], SS32[3][3][3], SS33[3][3][3];
		
		double rhoUu[3][3][3], rhoUv[3][3][3], rhoUw[3][3][3];
		double rhoVu[3][3][3], rhoVv[3][3][3], rhoVw[3][3][3];
		double rhoWu[3][3][3], rhoWv[3][3][3], rhoWw[3][3][3];
		double rhoSS11[3][3][3], rhoSS12[3][3][3], rhoSS13[3][3][3], rhoSS21[3][3][3], rhoSS22[3][3][3], rhoSS23[3][3][3], rhoSS31[3][3][3], rhoSS32[3][3][3], rhoSS33[3][3][3];
		
		
		
		int p,q,r;
		for(p=-1; p<=1; p++)
		for(q=-1; q<=1; q++)
		for(r=-1; r<=1; r++) {
			int R=r+1, Q=q+1, P=p+1;
			int K=k+r, J=j+q, I=i+p;
			
			u[R][Q][P] = ucat[K][J][I].x;
			v[R][Q][P] = ucat[K][J][I].y;
			w[R][Q][P] = ucat[K][J][I].z;
			
			// metric tensors are also test-filtered
			
			U[R][Q][P] = u[R][Q][P]*csi[K][J][I].x + v[R][Q][P]*csi[K][J][I].y + w[R][Q][P]*csi[K][J][I].z;
			V[R][Q][P] = u[R][Q][P]*eta[K][J][I].x + v[R][Q][P]*eta[K][J][I].y + w[R][Q][P]*eta[K][J][I].z;
			W[R][Q][P] = u[R][Q][P]*zet[K][J][I].x + v[R][Q][P]*zet[K][J][I].y + w[R][Q][P]*zet[K][J][I].z;
			
			Uu[R][Q][P] = U[R][Q][P] * u[R][Q][P];
			Uv[R][Q][P] = U[R][Q][P] * v[R][Q][P];
			Uw[R][Q][P] = U[R][Q][P] * w[R][Q][P];
			
			Vu[R][Q][P] = V[R][Q][P] * u[R][Q][P];
			Vv[R][Q][P] = V[R][Q][P] * v[R][Q][P];
			Vw[R][Q][P] = V[R][Q][P] * w[R][Q][P];
			
			Wu[R][Q][P] = W[R][Q][P] * u[R][Q][P];
			Wv[R][Q][P] = W[R][Q][P] * v[R][Q][P];
			Ww[R][Q][P] = W[R][Q][P] * w[R][Q][P];
			
			const double du_dx = Ax[K][J][I].x, du_dy = Ax[K][J][I].y, du_dz = Ax[K][J][I].z;
			const double dv_dx = Ay[K][J][I].x, dv_dy = Ay[K][J][I].y, dv_dz = Ay[K][J][I].z;
			const double dw_dx = Az[K][J][I].x, dw_dy = Az[K][J][I].y, dw_dz = Az[K][J][I].z;
		
			const double Sxx = 0.5*( du_dx + du_dx ), Sxy = 0.5*(du_dy + dv_dx), Sxz = 0.5*(du_dz + dw_dx);
			const double Syx = Sxy, Syy = 0.5*(dv_dy + dv_dy),	Syz = 0.5*(dv_dz + dw_dy);
			const double Szx = Sxz, Szy=Syz, Szz = 0.5*(dw_dz + dw_dz);
					
			S11[R][Q][P] = Sxx, S12[R][Q][P] = Sxy, S13[R][Q][P] = Sxz;
			S21[R][Q][P] = Syx, S22[R][Q][P] = Syy, S23[R][Q][P] = Syz;
			S31[R][Q][P] = Szx, S32[R][Q][P] = Szy, S33[R][Q][P] = Szz;
			
			S[R][Q][P] = Sabs[K][J][I];
			
			SS11[R][Q][P] = S11[R][Q][P]*S[R][Q][P], SS12[R][Q][P] = S12[R][Q][P]*S[R][Q][P], SS13[R][Q][P] = S13[R][Q][P]*S[R][Q][P];
			SS21[R][Q][P] = S21[R][Q][P]*S[R][Q][P], SS22[R][Q][P] = S22[R][Q][P]*S[R][Q][P], SS23[R][Q][P] = S23[R][Q][P]*S[R][Q][P];
			SS31[R][Q][P] = S31[R][Q][P]*S[R][Q][P], SS32[R][Q][P] = S32[R][Q][P]*S[R][Q][P], SS33[R][Q][P] = S33[R][Q][P]*S[R][Q][P];
			
			if(levelset) {
				d[R][Q][P] = rho[K][J][I];
				
				rhoUu[R][Q][P] = U[R][Q][P] * u[R][Q][P] * d[R][Q][P];	// (rho * U) * (rho * u) / rho = rho * U * u
				rhoUv[R][Q][P] = U[R][Q][P] * v[R][Q][P] * d[R][Q][P];
				rhoUw[R][Q][P] = U[R][Q][P] * w[R][Q][P] * d[R][Q][P];
				
				rhoVu[R][Q][P] = V[R][Q][P] * u[R][Q][P] * d[R][Q][P];
				rhoVv[R][Q][P] = V[R][Q][P] * v[R][Q][P] * d[R][Q][P];
				rhoVw[R][Q][P] = V[R][Q][P] * w[R][Q][P] * d[R][Q][P];
				
				rhoWu[R][Q][P] = W[R][Q][P] * u[R][Q][P] * d[R][Q][P];
				rhoWv[R][Q][P] = W[R][Q][P] * v[R][Q][P] * d[R][Q][P];
				rhoWw[R][Q][P] = W[R][Q][P] * w[R][Q][P] * d[R][Q][P];
			
				rhoSS11[R][Q][P] = S11[R][Q][P]*S[R][Q][P]*d[R][Q][P], rhoSS12[R][Q][P] = S12[R][Q][P]*S[R][Q][P]*d[R][Q][P], rhoSS13[R][Q][P] = S13[R][Q][P]*S[R][Q][P]*d[R][Q][P];
				rhoSS21[R][Q][P] = S21[R][Q][P]*S[R][Q][P]*d[R][Q][P], rhoSS22[R][Q][P] = S22[R][Q][P]*S[R][Q][P]*d[R][Q][P], rhoSS23[R][Q][P] = S23[R][Q][P]*S[R][Q][P]*d[R][Q][P];
				rhoSS31[R][Q][P] = S31[R][Q][P]*S[R][Q][P]*d[R][Q][P], rhoSS32[R][Q][P] = S32[R][Q][P]*S[R][Q][P]*d[R][Q][P], rhoSS33[R][Q][P] = S33[R][Q][P]*S[R][Q][P]*d[R][Q][P];
			}
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
		
		double weight[3][3][3];
		double sum_vol=0;
		//get_weight (i, j, k, mx, my, mz, aj, nvert, 0.1, weight);
		
		
		for(p=-1; p<=1; p++)
		for(q=-1; q<=1; q++)
		for(r=-1; r<=1; r++) {
			int R=r+1, Q=q+1, P=p+1;
			int K=k+r, J=j+q, I=i+p;
			
			sum_weight += weight[R][Q][P] * coef[R][Q][P];
			if(nvert[K][J][I]<0.1) {
				sum_vol += 1./aj[K][J][I] * coef[R][Q][P];
				weight[R][Q][P] = 1;
			}
			else weight[R][Q][P] = 0;
		}
		
		filter = pow( 1./aj[k][j][i], 1./3. );
		
		if(testfilter_ik) test_filter = pow(5.0, 1./3.) * filter;
		else {
			//test_filter = 2.0 * filter;
			//test_filter = pow( sum_weight, 1./3. );
			test_filter = pow( sum_vol, 1./3. );
		}
				
		double _U=integrate_testfilter_simpson(U, weight);
		double _V=integrate_testfilter_simpson(V, weight);
		double _W=integrate_testfilter_simpson(W, weight);
		
		double _u=integrate_testfilter_simpson(u, weight);
		double _v=integrate_testfilter_simpson(v, weight);
		double _w=integrate_testfilter_simpson(w, weight);
		double _d=1;
		if(levelset) _d=integrate_testfilter_simpson(d, weight);
		/*
		double _U = _u * csi[k][j][i].x + _v * csi[k][j][i].y + _w * csi[k][j][i].z;
		double _V = _u * eta[k][j][i].x + _v * eta[k][j][i].y + _w * eta[k][j][i].z;
		double _W = _u * zet[k][j][i].x + _v * zet[k][j][i].y + _w * zet[k][j][i].z;
		*/
		if(levelset) {
			Lij[0][0] = integrate_testfilter_simpson(rhoUu, weight) - _U*_u*_d;
			Lij[0][1] = integrate_testfilter_simpson(rhoUv, weight) - _U*_v*_d;
			Lij[0][2] = integrate_testfilter_simpson(rhoUw, weight) - _U*_w*_d;
			Lij[1][0] = integrate_testfilter_simpson(rhoVu, weight) - _V*_u*_d;
			Lij[1][1] = integrate_testfilter_simpson(rhoVv, weight) - _V*_v*_d;
			Lij[1][2] = integrate_testfilter_simpson(rhoVw, weight) - _V*_w*_d;
			Lij[2][0] = integrate_testfilter_simpson(rhoWu, weight) - _W*_u*_d;
			Lij[2][1] = integrate_testfilter_simpson(rhoWv, weight) - _W*_v*_d;
			Lij[2][2] = integrate_testfilter_simpson(rhoWw, weight) - _W*_w*_d;
		} 
		else {
			Lij[0][0] = integrate_testfilter_simpson(Uu, weight) - _U*_u;
			Lij[0][1] = integrate_testfilter_simpson(Uv, weight) - _U*_v;
			Lij[0][2] = integrate_testfilter_simpson(Uw, weight) - _U*_w;
			Lij[1][0] = integrate_testfilter_simpson(Vu, weight) - _V*_u;
			Lij[1][1] = integrate_testfilter_simpson(Vv, weight) - _V*_v;
			Lij[1][2] = integrate_testfilter_simpson(Vw, weight) - _V*_w;
			Lij[2][0] = integrate_testfilter_simpson(Wu, weight) - _W*_u;
			Lij[2][1] = integrate_testfilter_simpson(Wv, weight) - _W*_v;
			Lij[2][2] = integrate_testfilter_simpson(Ww, weight) - _W*_w;
		}		
				
		/***************/
		Sij_hat[0][0] = integrate_testfilter_simpson(S11, weight);	
		Sij_hat[0][1] = integrate_testfilter_simpson(S12, weight);	
		Sij_hat[0][2] = integrate_testfilter_simpson(S13, weight);
		Sij_hat[1][0] = integrate_testfilter_simpson(S21, weight);	
		Sij_hat[1][1] = integrate_testfilter_simpson(S22, weight);	
		Sij_hat[1][2] = integrate_testfilter_simpson(S23, weight);
		Sij_hat[2][0] = integrate_testfilter_simpson(S31, weight);	
		Sij_hat[2][1] = integrate_testfilter_simpson(S32, weight);	
		Sij_hat[2][2] = integrate_testfilter_simpson(S33, weight);
		
		S_hat=0;
		for(a=0; a<3; a++)
		for(b=0; b<3; b++) {
			S_hat += pow( Sij_hat[a][b], 2. );
		}
		S_hat = sqrt ( 2 * S_hat );
		/***************/
		
		if(levelset) {
			rhoSSij_hat[0][0] = integrate_testfilter_simpson(rhoSS11, weight);	
			rhoSSij_hat[0][1] = integrate_testfilter_simpson(rhoSS12, weight);	
			rhoSSij_hat[0][2] = integrate_testfilter_simpson(rhoSS13, weight);
			rhoSSij_hat[1][0] = integrate_testfilter_simpson(rhoSS21, weight);
			rhoSSij_hat[1][1] = integrate_testfilter_simpson(rhoSS22, weight);	
			rhoSSij_hat[1][2] = integrate_testfilter_simpson(rhoSS23, weight);
			rhoSSij_hat[2][0] = integrate_testfilter_simpson(rhoSS31, weight);	
			rhoSSij_hat[2][1] = integrate_testfilter_simpson(rhoSS32, weight);	
			rhoSSij_hat[2][2] = integrate_testfilter_simpson(rhoSS33, weight);
		}
		else {
			SSij_hat[0][0] = integrate_testfilter_simpson(SS11, weight);	
			SSij_hat[0][1] = integrate_testfilter_simpson(SS12, weight);	
			SSij_hat[0][2] = integrate_testfilter_simpson(SS13, weight);
			SSij_hat[1][0] = integrate_testfilter_simpson(SS21, weight);
			SSij_hat[1][1] = integrate_testfilter_simpson(SS22, weight);	
			SSij_hat[1][2] = integrate_testfilter_simpson(SS23, weight);
			SSij_hat[2][0] = integrate_testfilter_simpson(SS31, weight);	
			SSij_hat[2][1] = integrate_testfilter_simpson(SS32, weight);	
			SSij_hat[2][2] = integrate_testfilter_simpson(SS33, weight);
		}
		
		//S_hat = integrate_testfilter_simpson(S, weight);
		
		
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
			if(levelset) {
				Mij_cat[a][b] = - _d * pow( test_filter, 2. ) * S_hat * Sij_hat[a][b] + pow( filter, 2. ) * rhoSSij_hat[a][b];
			}
			else {
				Mij_cat[a][b] = - pow( test_filter, 2. ) * S_hat * Sij_hat[a][b] + pow( filter, 2. ) * SSij_hat[a][b];
			}
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
			num += Lij[b][a] * Mij[a][q] * G[b][q];
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
			if(user->bctype[5]==4) {
				LM_avg = (1.0 * LM0[0][1][1] + 4.0 * LM0[1][1][1] + 1.0 * LM0[2][1][1]) / 6.;
				MM_avg = (1.0 * MM0[0][1][1] + 4.0 * MM0[1][1][1] + 1.0 * MM0[2][1][1]) / 6.;
			}
			else {
				LM_avg = integrate_testfilter_simpson(LM0,weight);
				MM_avg = integrate_testfilter_simpson(MM0,weight);
			}
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
	
	DMDAVecRestoreArray(da, user->lLM, &LM);//
	DMDAVecRestoreArray(da, user->lMM, &MM);//

	if(levelset) {
		DMDAVecRestoreArray(da, user->lDensity, &rho);
		DMDAVecRestoreArray(da, user->lLevelset, &level);
	}

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

	double lmax_norm=0, max_norm;
	PetscInt p;
	
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
	
	PetscReal ***Cs, ***lnu_t, ***nvert, ***aj, ***ustar, ***rho;
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

	VecSet(user->lNu_t, 0);
	
	if(levelset) {
		DMDAVecGetArray(da, user->lDensity, &rho);
	}
	
	DMDAVecGetArray(fda, user->lUcat,  &ucat);
	DMDAVecGetArray(fda, user->lCsi, &csi);
	DMDAVecGetArray(fda, user->lEta, &eta);
	DMDAVecGetArray(fda, user->lZet, &zet);
	
	DMDAVecGetArray(da, user->lNvert, &nvert);
	DMDAVecGetArray(da, user->lAj, &aj);
	DMDAVecGetArray(da, user->lNu_t, &lnu_t);
	DMDAVecGetArray(da, user->lCs, &Cs);
	DMDAVecGetArray(da, user->lUstar, &ustar);
	
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
		Compute_du_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz,
							&du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
		
		double Sxx = 0.5*( du_dx + du_dx ), Sxy = 0.5*(du_dy + dv_dx), Sxz = 0.5*(du_dz + dw_dx);
		double Syx = Sxy, Syy = 0.5*(dv_dy + dv_dy),	Syz = 0.5*(dv_dz + dw_dy);
		double Szx = Sxz, Szy=Syz, Szz = 0.5*(dw_dz + dw_dz);
	
		double Sabs = sqrt( 2.0*( Sxx*Sxx + Sxy*Sxy + Sxz*Sxz + Syx*Syx + Syy*Syy + Syz*Syz + Szx*Szx + Szy*Szy + Szz*Szz ) );
		
		double filter  = pow( 1./aj[k][j][i],1./3.);
	
		lnu_t[k][j][i] = Cs[k][j][i] * pow ( filter, 2.0 ) * Sabs;
		

		if(levelset) lnu_t[k][j][i] *= rho[k][j][i];

		if(les && wallfunction==2 && nvert[k][j][i]+nvert[k][j][i+1]+nvert[k][j][i-1]+nvert[k][j+1][i]+nvert[k][j-1][i]+nvert[k+1][j][i]+nvert[k-1][j][i]>0.1) lnu_t[k][j][i]=0;
		
	}
	
	if(immersed && wallfunction==2) {
		DMDAVecRestoreArray(da, user->lNu_t, &lnu_t);
		
		DMDALocalToLocalBegin(da, user->lNu_t, INSERT_VALUES, user->lNu_t);
		DMDALocalToLocalEnd(da, user->lNu_t, INSERT_VALUES, user->lNu_t);
		
		DMDAVecGetArray(da, user->lNu_t, &lnu_t);
		
		for(int ibi=0; ibi<NumberOfBodies; ibi++)
		{
			IBMNodes *ibm = &ibm_ptr[ibi];
			
			IBMListNode *current;
			current = user->ibmlist[ibi].head;
			while (current) {
				IBMInfo *ibminfo = &current->ibm_intp;
				
				current = current->next;
				double sb = ibminfo->d_s, sc = sb + ibminfo->d_i;
				int ni = ibminfo->cell;
				int ip1 = ibminfo->i1, jp1 = ibminfo->j1, kp1 = ibminfo->k1;
				int ip2 = ibminfo->i2, jp2 = ibminfo->j2, kp2 = ibminfo->k2;
				int ip3 = ibminfo->i3, jp3 = ibminfo->j3, kp3 = ibminfo->k3;
				i = ibminfo->ni, j= ibminfo->nj, k = ibminfo->nk;
				double sk1  = ibminfo->cr1, sk2 = ibminfo->cr2, sk3 = ibminfo->cr3;
				double cs1 = ibminfo->cs1, cs2 = ibminfo->cs2, cs3 = ibminfo->cs3;
				double nx = ibm->nf_x[ni], ny = ibm->nf_y[ni], nz = ibm->nf_z[ni];
				
				Cmpnts Ua, Uc;
				if (ni>=0) {
					Ua.x = ibm->u[ibm->nv1[ni]].x * cs1 + ibm->u[ibm->nv2[ni]].x * cs2 + ibm->u[ibm->nv3[ni]].x * cs3;
					Ua.y = ibm->u[ibm->nv1[ni]].y * cs1 + ibm->u[ibm->nv2[ni]].y * cs2 + ibm->u[ibm->nv3[ni]].y * cs3;
					Ua.z = ibm->u[ibm->nv1[ni]].z * cs1 + ibm->u[ibm->nv2[ni]].z * cs2 + ibm->u[ibm->nv3[ni]].z * cs3;
				}
				else {
					Ua.x = Ua.y = Ua.z = 0;
				}
				
				Uc.x = (ucat[kp1][jp1][ip1].x * sk1 + ucat[kp2][jp2][ip2].x * sk2 + ucat[kp3][jp3][ip3].x * sk3);
				Uc.y = (ucat[kp1][jp1][ip1].y * sk1 + ucat[kp2][jp2][ip2].y * sk2 + ucat[kp3][jp3][ip3].y * sk3);
				Uc.z = (ucat[kp1][jp1][ip1].z * sk1 + ucat[kp2][jp2][ip2].z * sk2 + ucat[kp3][jp3][ip3].z * sk3);
				
				//double yplus = ustar[k][j][i] * sb * user->ren;
				double nu = 1./user->ren;
				double nu_t_c = (lnu_t[kp1][jp1][ip1] * sk1 + lnu_t[kp2][jp2][ip2] * sk2 + lnu_t[kp3][jp3][ip3] * sk3);
				double eps=1.e-5;
				double dUc_ds = u_Cabot(nu, sc+eps, ustar[k][j][i], 0) - u_Cabot(nu, sc-eps, ustar[k][j][i], 0);
				dUc_ds /= (2.0 * eps);
				double f1 = ustar[k][j][i] * ustar[k][j][i];
				double f2 = (nu+nu_t_c) *  dUc_ds;
								
				double u_c = Uc.x - Ua.x, v_c = Uc.y - Ua.y, w_c = Uc.z - Ua.z;
				double un = u_c * nx + v_c * ny + w_c * nz;
				double ut = u_c - un * nx, vt = v_c - un * ny, wt = w_c - un * nz;
				double ut_mag_c = sqrt( ut*ut + vt*vt + wt*wt );
				double ut_mag_b  = u_Cabot(nu, sb, ustar[k][j][i], 0);
				
				//lnu_t[k][j][i] = nu_t_c * sb / sc;
				//lnu_t[k][j][i] = - nu + (f1*(sc-sb) + f2*sb) * (sc-sb) / sc / (ut_mag_c - ut_mag_b);
				//printf("nu_t = %f %f %f\n", ut_mag_c, ut_mag_b, ut_mag_c - (f1*(sc-sb) + f2*sb) * (sc-sb) / sc / (0 + nu));
				
				
				//ut_mag_b =  ut_mag_c - (f1*(sc-sb) + f2*sb) * (sc-sb) / sc / (lnu_t[k][j][i] + nu);
				//lnu_t[k][j][i] = near_wall_eddy_viscosity(yplus) / user->ren;
				/*
					Uc_t = Uc - (Uc . n) n
					
					
					f1 = ustar*ustar;
					f2 = shear_stress at C;
					(nu+nu_t)*(Uc - Ub)/(sc-sb) = (f1*(sc-sb) + f2*sb) / sc
					nu + nu_t = (f1*(sc-sb) + f2*sb) * (sc-sb) / sc / (Uc - Ub)
					nu_t = - nu + (f1*(sc-sb) + f2*sb) * (sc-sb) / sc / (Uc - Ub)
				*/
			};
		}
	}
	
	if(levelset) {
		DMDAVecRestoreArray(da, user->lDensity, &rho);
	}
	
	DMDAVecRestoreArray(fda, user->lUcat,  &ucat);
	DMDAVecRestoreArray(fda, user->lCsi, &csi);
	DMDAVecRestoreArray(fda, user->lEta, &eta);
	DMDAVecRestoreArray(fda, user->lZet, &zet);
	
	DMDAVecRestoreArray(da, user->lNvert, &nvert);
	DMDAVecRestoreArray(da, user->lAj, &aj);
	DMDAVecRestoreArray(da, user->lNu_t, &lnu_t);
	DMDAVecRestoreArray(da, user->lCs, &Cs);
	DMDAVecRestoreArray(da, user->lUstar, &ustar);
	
	DMDALocalToLocalBegin(da, user->lNu_t, INSERT_VALUES, user->lNu_t);
	DMDALocalToLocalEnd(da, user->lNu_t, INSERT_VALUES, user->lNu_t);
	
	if(periodic) {
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
	}
	
	double lmax_norm=0, max_norm;
	PetscInt p;
	
	VecMax(user->lNu_t, &p, &lmax_norm);
	GlobalMax_All(&lmax_norm, &max_norm, PETSC_COMM_WORLD);
	PetscPrintf(PETSC_COMM_WORLD, "Max Nu_t = %e\n", max_norm);
};



void Compute_Prt(UserCtx *user)
{
	
	Vec Csi = user->lCsi, Eta = user->lEta, Zet = user->lZet;
	Vec Ucont = user->lUcont, Ucat = user->lUcat, Tmprt = user->lTmprt;
	Cmpnts	***ucont, ***ucat, ***dtmprt;
	PetscReal ***tmprt, ***tmprt_f;	
	Cmpnts	***csi, ***eta, ***zet;
	PetscReal	***nvert, ***pr_t, ***Cs, ***lnu_t;

	DM		da = user->da, fda = user->fda;
	DMDALocalInfo	info;
	PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
	PetscInt	mx, my, mz; // Dimensions in three directions
	PetscInt	i, j, k;
  
	PetscReal	***aj;
	Cmpnts ***ucat_f, ***cent;//, ;
	double ***LM, ***MM;//, ***NM;
	PetscReal ***Sabs;//, ***Cs, ***lCs_o;//, ***nu_t;

	PetscInt	lxs, lxe, lys, lye, lzs, lze;
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

	Vec lS;
	Vec lUcat_f;
	Vec Tmprt_f;
	Vec dTmprt;	

	DMDAVecGetArray(fda, Ucont, &ucont);
	DMDAVecGetArray(fda, Ucat,  &ucat);
	DMDAVecGetArray(fda, Csi, &csi);
	DMDAVecGetArray(fda, Eta, &eta);
	DMDAVecGetArray(fda, Zet, &zet);
	DMDAVecGetArray(da, user->lNvert, &nvert);
	DMDAVecGetArray(da, user->lAj, &aj); 
	DMDAVecGetArray(da, user->lCs, &Cs);
	DMDAVecGetArray(da, user->lNu_t, &lnu_t);
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

	DMDAVecGetArray(da, Tmprt, &tmprt);
	VecSet(user->lPr_t, 0);
	DMDAVecGetArray(da, user->lPr_t, &pr_t);
  
	VecDuplicate(user->lP, &user->lLM);
	VecDuplicate(user->lP, &user->lMM);
	
	VecSet(user->lLM, 0);
	VecSet(user->lMM, 0);
	
	VecDuplicate(user->lNvert, &lS);
	VecDuplicate(user->lUcont, &lUcat_f);
	VecDuplicate(user->lUcont, &dTmprt);
	VecDuplicate(Tmprt, &Tmprt_f);

	VecSet(lS, 0);
	
	DMDAVecGetArray(da, user->lLM, &LM);//
	DMDAVecGetArray(da, user->lMM, &MM);//
	//DAVecGetArray(da, user->lNM, &NM);//
  	
	DMDAVecGetArray(da, lS, &Sabs);//
	DMDAVecGetArray(fda, lUcat_f, &ucat_f);//
	DMDAVecGetArray(fda, dTmprt, &dtmprt);//
	DMDAVecGetArray(da, Tmprt_f, &tmprt_f);

  	
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
		
		double weight[3][3][3];
		double u[3][3][3];
		double v[3][3][3];
		double w[3][3][3];
		double t[3][3][3];

		get_weight (i, j, k, mx, my, mz, aj, nvert, 1.1, weight);
		
		for(int p=-1; p<=1; p++)
		for(int q=-1; q<=1; q++)
		for(int r=-1; r<=1; r++) {
			int R=r+1, Q=q+1, P=p+1;
			int K=k+r, J=j+q, I=i+p;
			u[R][Q][P] = ucat[K][J][I].x;
			v[R][Q][P] = ucat[K][J][I].y;
			w[R][Q][P] = ucat[K][J][I].z;

			t[R][Q][P] = tmprt[K][J][I];
		}
		
//		ucat_f[k][j][i].x = integrate_testfilter_simpson(u, weight);
//		ucat_f[k][j][i].y = integrate_testfilter_simpson(v, weight);
//		ucat_f[k][j][i].z = integrate_testfilter_simpson(w, weight);

//		tmprt_f[k][j][i] = integrate_testfilter_simpson(t, weight);


                ucat_f[k][j][i].x = integrate_xztestfilter_simpson(u, weight);
                ucat_f[k][j][i].y = integrate_xztestfilter_simpson(v, weight);
                ucat_f[k][j][i].z = integrate_xztestfilter_simpson(w, weight);

                tmprt_f[k][j][i] = integrate_xztestfilter_simpson(t, weight);


                double dtdc, dtde, dtdz;
                double dt_dx, dt_dy, dt_dz;

                Compute_dscalar_center (i, j, k, mx, my, mz, tmprt, nvert,  &dtdc, &dtde, &dtdz );
                Compute_dscalar_dxyz (csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dtdc, dtde, dtdz, &dt_dx, &dt_dy, &dt_dz );

		dtmprt[k][j][i].x = dt_dx;
		dtmprt[k][j][i].y = dt_dy;
		dtmprt[k][j][i].z = dt_dz;


	}
	
	DMDAVecRestoreArray(da, lS, &Sabs);//
	DMDAVecRestoreArray(fda, lUcat_f, &ucat_f);//
	DMDAVecRestoreArray(fda, dTmprt, &dtmprt);//
	DMDAVecRestoreArray(da, Tmprt_f, &tmprt_f);
		
	DMDALocalToLocalBegin(da, lS, INSERT_VALUES, lS);
	DMDALocalToLocalEnd(da, lS, INSERT_VALUES, lS);
	
	DMDALocalToLocalBegin(fda, lUcat_f, INSERT_VALUES, lUcat_f);
	DMDALocalToLocalEnd(fda, lUcat_f, INSERT_VALUES, lUcat_f);

        DMDALocalToLocalBegin(da, Tmprt_f, INSERT_VALUES, Tmprt_f);
        DMDALocalToLocalEnd(da, Tmprt_f, INSERT_VALUES, Tmprt_f);

        DMDALocalToLocalBegin(fda, dTmprt, INSERT_VALUES, dTmprt);
        DMDALocalToLocalEnd(fda, dTmprt, INSERT_VALUES, dTmprt);


	DMDAVecGetArray(da, lS, &Sabs);//
	DMDAVecGetArray(fda, lUcat_f, &ucat_f);//
	DMDAVecGetArray(fda, dTmprt, &dtmprt);//
	DMDAVecGetArray(da, Tmprt_f, &tmprt_f);


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
			ucat_f[k][j][i] = ucat_f[c][b][a];
			tmprt_f[k][j][i] = tmprt_f[c][b][a];
			dtmprt[k][j][i] = dtmprt[c][b][a];
		}
	}
///




	
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
		
		double filter, test_filter;
		double S[3][3][3], S_hat;
		double U[3][3][3], V[3][3][3], W[3][3][3];
		double u[3][3][3], v[3][3][3], w[3][3][3], t[3][3][3];		

		double Ut[3][3][3], Vt[3][3][3], Wt[3][3][3], dt_dx, dt_dy, dt_dz;
		double dhSdtdx[3][3][3], dhSdtdy[3][3][3], dhSdtdz[3][3][3];
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
						
			S[R][Q][P] = 0;

                        u[R][Q][P] = ucat[K][J][I].x;
                        v[R][Q][P] = ucat[K][J][I].y;
                        w[R][Q][P] = ucat[K][J][I].z;
                        t[R][Q][P] = tmprt[K][J][I];
			
			
			U[R][Q][P] = u[R][Q][P]*csi[K][J][I].x + v[R][Q][P]*csi[K][J][I].y + w[R][Q][P]*csi[K][J][I].z;			
			V[R][Q][P] = u[R][Q][P]*eta[K][J][I].x + v[R][Q][P]*eta[K][J][I].y + w[R][Q][P]*eta[K][J][I].z;			
			W[R][Q][P] = u[R][Q][P]*zet[K][J][I].x + v[R][Q][P]*zet[K][J][I].y + w[R][Q][P]*zet[K][J][I].z;
			
			Ut[R][Q][P] = U[R][Q][P]*t[R][Q][P];
			Vt[R][Q][P] = V[R][Q][P]*t[R][Q][P];
			Wt[R][Q][P] = W[R][Q][P]*t[R][Q][P];

                        //Ut[R][Q][P] = u[R][Q][P]*t[R][Q][P];
                        //Vt[R][Q][P] = v[R][Q][P]*t[R][Q][P];
                        //Wt[R][Q][P] = w[R][Q][P]*t[R][Q][P];

			
			S[R][Q][P] = Sabs[K][J][I];

			filter = pow( 1./aj[K][J][I], 1./3. );

			dhSdtdx[R][Q][P] = pow(filter,2)*S[R][Q][P]*dtmprt[K][J][I].x;
			dhSdtdy[R][Q][P] = pow(filter,2)*S[R][Q][P]*dtmprt[K][J][I].y;
			dhSdtdz[R][Q][P] = pow(filter,2)*S[R][Q][P]*dtmprt[K][J][I].z;
			
		}

		double _dhSdtdx[3];
		
//		_dhSdtdx[0] = integrate_testfilter_simpson(dhSdtdx, weight);
//		_dhSdtdx[1] = integrate_testfilter_simpson(dhSdtdy, weight);
//		_dhSdtdx[2] = integrate_testfilter_simpson(dhSdtdz, weight);

                _dhSdtdx[0] = integrate_xztestfilter_simpson(dhSdtdx, weight);
                _dhSdtdx[1] = integrate_xztestfilter_simpson(dhSdtdy, weight);
                _dhSdtdx[2] = integrate_xztestfilter_simpson(dhSdtdz, weight);




		double _Ut[3];
//                _Ut[0] = integrate_testfilter_simpson(Ut, weight);
//                _Ut[1] = integrate_testfilter_simpson(Vt, weight);
//                _Ut[2] = integrate_testfilter_simpson(Wt, weight);

                _Ut[0] = integrate_xztestfilter_simpson(Ut, weight);
                _Ut[1] = integrate_xztestfilter_simpson(Vt, weight);
                _Ut[2] = integrate_xztestfilter_simpson(Wt, weight);


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
		

		double area = sqrt( eta[k][j][i].x*eta[k][j][i].x + eta[k][j][i].y*eta[k][j][i].y + eta[k][j][i].z*eta[k][j][i].z );
                /*if (!testfilter_ik)*/ filter = pow( 1./aj[k][j][i], 1./3. );
//                if (testfilter_ik) filter = pow(area,1.0/2.0);

		
		test_filter = 2.0*filter;

                if(testfilter_ik) test_filter = pow(5.0, 1./3.) * filter; 
		else {
			//test_filter = 2.0 * filter;
			test_filter = pow( sum_weight, 1./3. );
		}
		
		double _u=ucat_f[k][j][i].x;
		double _v=ucat_f[k][j][i].y;
		double _w=ucat_f[k][j][i].z;
		double _t=tmprt_f[k][j][i];
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
		S_hat = sqrt( 2.0 *( _Sxx*_Sxx + _Sxy*_Sxy + _Sxz*_Sxz + _Syx*_Syx + _Syy*_Syy + _Syz*_Syz + _Szx*_Szx + _Szy*_Szy + _Szz*_Szz ) );

		double _dtdc, _dtde, _dtdz;
		double _dt_dx, _dt_dy, _dt_dz;

                Compute_dscalar_center (i, j, k, mx, my, mz, tmprt_f, nvert,  &_dtdc, &_dtde, &_dtdz );
                Compute_dscalar_dxyz (csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc,_dtdc, _dtde, _dtdz, &_dt_dx, &_dt_dy, &_dt_dz );

		double _U_t[3];

                _U_t[0] = _U*_t;
                _U_t[1] = _V*_t;
                _U_t[2] = _W*_t;

                //_U_t[0] = _u*_t;
                //_U_t[1] = _v*_t;
                //_U_t[2] = _w*_t;

	
                double dtdx1[3];
	
		dtdx1[0]=dtmprt[k][j][i].x; 
		dtdx1[1]=dtmprt[k][j][i].y; 
		dtdx1[2]=dtmprt[k][j][i].z; 
			

		double _dt_dx1[3];	
		
		_dt_dx1[0]=_dt_dx;
		_dt_dx1[1]=_dt_dy;
		_dt_dx1[2]=_dt_dz;

		
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
	

		double NN[3];

		NN[0] = (pow(test_filter,2)*S_hat*_dt_dx1[0]-_dhSdtdx[0])*csi0;	
		NN[1] = (pow(test_filter,2)*S_hat*_dt_dx1[1]-_dhSdtdx[1])*eta1;	
		NN[2] = (pow(test_filter,2)*S_hat*_dt_dx1[2]-_dhSdtdx[2])*zet2;	

		for(p==0;p<3;p++){	
                        MM[k][j][i] += pow (NN[p],2);
                        LM[k][j][i] += (_U_t[p]-_Ut[p])*NN[p];

//                        MM[k][j][i] += NN[p]*dtdx1[p];
//                        LM[k][j][i] += (_U_t[p]-_Ut[p])*dtdx1[p];


		}
    	}
	
	DMDAVecRestoreArray(da, user->lLM, &LM);//
	DMDAVecRestoreArray(da, user->lMM, &MM);//
	//DAVecRestoreArray(da, user->lNM, &NM);//
	
	DMDALocalToLocalBegin(da, user->lLM, INSERT_VALUES, user->lLM);
	DMDALocalToLocalEnd(da, user->lLM, INSERT_VALUES, user->lLM);
	DMDALocalToLocalBegin(da, user->lMM, INSERT_VALUES, user->lMM);
	DMDALocalToLocalEnd(da, user->lMM, INSERT_VALUES, user->lMM);
	//DALocalToLocalBegin(da, user->lNM, INSERT_VALUES, user->lNM);
	//DALocalToLocalEnd(da, user->lNM, INSERT_VALUES, user->lNM);
	  
	DMDAVecGetArray(da, user->lLM, &LM);//
	DMDAVecGetArray(da, user->lMM, &MM);//
	//DAVecGetArray(da, user->lNM, &NM);//
	
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
			pr_t[k][j][i] = 0.0;
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

		if ( i_homo_filter || j_homo_filter || k_homo_filter || les_prt==3 ) {
			LM_avg = LM[k][j][i];
			MM_avg = MM[k][j][i];
		}
		else {
//                        LM_avg = integrate_testfilter_simpson(LM0,weight);
//                        MM_avg = integrate_testfilter_simpson(MM0,weight);

			LM_avg = integrate_xztestfilter_simpson(LM0,weight);
			MM_avg = integrate_xztestfilter_simpson(MM0,weight);
		}
		
		C = 0.5*LM_avg / (MM_avg + prt_eps );
	
//		printf("LM_avg, MM_avg, Cs %le %le %le \n", LM_avg, MM_avg, Cs[k][j][i]);


//		if ( les==3 ) {
//			if(ti<100 && tistart==0 && !rstart_flg) { }
//			else {
//				C = (1.0 - 0.001) * [k][j][i] + 0.001 * C;
//			}
//		}
		
		double filter = pow( 1./aj[k][j][i], 1./3. );
		pr_t[k][j][i] = C*pow(filter,2)*Sabs[k][j][i];
	}
	
	if( les_prt==3 ) {}
	else if ( i_homo_filter && k_homo_filter ) {


	PetscPrintf(PETSC_COMM_WORLD, "ikfilter \n");	

		double C=0.0;

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
//			PetscPrintf(PETSC_COMM_WORLD, "%i %i %le %le \n ", j, total_count[j], J_LM[j], J_MM[j] );	
			if( total_count[j]>0) {
				J_LM[j] /= (double) (total_count[j]);
				J_MM[j] /= (double) (total_count[j]);
			}
		}
		
		for (j=lys; j<lye; j++)
		for (k=lzs; k<lze; k++)
		for (i=lxs; i<lxe; i++) {
			double filter = pow( 1./aj[k][j][i], 1./3. );
			C =  0.5*J_LM[j] / ( J_MM[j]+prt_eps);
			// pr_t[k][j][i] = C*Cs[k][j][i]; //C*pow(filter,2)*Sabs[k][j][i]; //C*Cs[k][j][i];
			pr_t[k][j][i] = C*pow(filter,2)*Sabs[k][j][i]; // RDT
		}
	}


	else if (i_homo_filter || j_homo_filter || k_homo_filter) {


	PetscPrintf(PETSC_COMM_WORLD, "i or kfiler \n");	
		int plane_size;

                double C=0.0;
		
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

							C =  0.5*J_LM[pos] / ( J_MM[pos] + prt_eps );

                					double filter = pow( 1./aj[k][j][i], 1./3. );
				                        //pr_t[k][j][i] = C*Cs[k][j][i]; //C*pow(filter,2)*Sabs[k][j][i]; //C*Cs[k][j][i];
				                        pr_t[k][j][i] = C*pow(filter,2)*Sabs[k][j][i];

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
                                                        C =  0.5*J_LM[pos] / ( J_MM[pos] + prt_eps );
                					double filter = pow( 1./aj[k][j][i], 1./3. );
				                        //pr_t[k][j][i] = C*Cs[k][j][i]; //C*pow(filter,2)*Sabs[k][j][i]; //C*Cs[k][j][i];
                                                        //pr_t[k][j][i] = C*Cs[k][j][i]; //C*pow(filter,2)*Sabs[k][j][i]; // C*Cs[k][j][i];
                                                        pr_t[k][j][i] = C*pow(filter,2)*Sabs[k][j][i]; 

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
                                                        C =  0.5*J_LM[pos] / ( J_MM[pos] + prt_eps );
                					double filter = pow( 1./aj[k][j][i], 1./3. );
				                        //pr_t[k][j][i] = C*Cs[k][j][i]; //C*pow(filter,2)*Sabs[k][j][i]; //C*Cs[k][j][i];
                                                        //pr_t[k][j][i] = C*Cs[k][j][i]; //C*pow(filter,2)*Sabs[k][j][i]; //C*Cs[k][j][i];
                                                        pr_t[k][j][i] = C*pow(filter,2)*Sabs[k][j][i]; 

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
			pr_t[k][j][i] = 0.0;
		}
		else {
			if(nvert[k][j][i]>0.1 && nvert[k][j][i]<1.1) {
//				pr_t[k][j][i] = Cs[k][j][i] * PetscMax(wall_cs, Cs[k][j][i]);	// stabilize at high Re, osl 0.005
			}
			pr_t[k][j][i] = PetscMax ( pr_t[k][j][i], 0);

		}


		// pr_t[k][j][i] = lnu_t[k][j][i]/(pr_t[k][j][i]+prt_eps);
	}
	
	DMDAVecRestoreArray(da, lS, &Sabs);//
	DMDAVecRestoreArray(fda, lUcat_f, &ucat_f);//
	

	DMDAVecRestoreArray(da, Tmprt_f, &tmprt_f);//
	DMDAVecRestoreArray(fda, dTmprt, &dtmprt);//

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
	DMDAVecRestoreArray(da, user->lNu_t, &lnu_t);
	
	DMDAVecRestoreArray(fda, user->lCent, &cent);
	VecDestroy(&user->lLM);//
	VecDestroy(&user->lMM);//
  
	VecDestroy(&lS);//
	VecDestroy(&lUcat_f);//
	VecDestroy(&Tmprt_f);//
	VecDestroy(&dTmprt);


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
	
	DMDAVecRestoreArray(da, Tmprt, &tmprt);
	DMDAVecRestoreArray(da, user->lPr_t, &pr_t);

	DMDALocalToLocalBegin(da, user->lPr_t, INSERT_VALUES, user->lPr_t);
	DMDALocalToLocalEnd(da, user->lPr_t, INSERT_VALUES, user->lPr_t);
	
	DMDAVecGetArray(da, user->lPr_t, &pr_t);
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


		if(flag) pr_t[k][j][i] = pr_t[c][b][a];
	}
	DMDAVecRestoreArray(da, user->lPr_t, &pr_t);

	double lmax_norm=0, max_norm;

	PetscInt p;
	
	if(testfilter_ik) PetscPrintf(PETSC_COMM_WORLD, "\nFilter type : Box filter homogeneous\n");
	else PetscPrintf(PETSC_COMM_WORLD, "Filter type : Box filter 3D\n");
	
	VecMax(user->lPr_t, &p, &lmax_norm);
//	PetscGlobalMax(&lmax_norm, &max_norm, PETSC_COMM_WORLD);

	GlobalMax_All(&lmax_norm, &max_norm, PETSC_COMM_WORLD);

	PetscPrintf(PETSC_COMM_WORLD, "Max Pr_t  = %e \n", max_norm);

//        TECIOOut_rhs1(user,  user->lPr_t);

//        PetscFinalize();
//        exit(0);

	
};

