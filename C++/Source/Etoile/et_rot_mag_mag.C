/*
 * Computes magnetic fields and derived quantities for rotating equilibrium
 *
 * (see file et_rot_mag.h for documentation)
 *
 */

/*
 *   Copyright (c) 2002 Emmanuel Marcq
 *   Copyright (c) 2002 Jerome Novak
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   LORENE is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with LORENE; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

char et_rot_mag_mag_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.7  2002/05/20 10:31:59  j_novak
 * *** empty log message ***
 *
 * Revision 1.6  2002/05/17 15:08:01  e_marcq
 *
 * Rotation progressive plug-in, units corrected, Q and a_j new member data
 *
 * Revision 1.5  2002/05/16 10:02:09  j_novak
 * Errors in stress energy tensor corrected
 *
 * Revision 1.4  2002/05/15 09:54:00  j_novak
 * First operational version
 *
 * Revision 1.3  2002/05/14 13:38:36  e_marcq
 *
 *
 * Unit update, new outputs
 *
 * Revision 1.1  2002/05/10 09:26:52  j_novak
 * Added new class Et_rot_mag for magnetized rotating neutron stars (under development)
 *
 *
 * $Header$
 *
 */

// Headers C
#include <stdlib.h>
#include <math.h>

// Headers Lorene
#include "et_rot_mag.h"
#include "utilitaires.h"
#include "param.h"
#include "graphique.h"

extern "C" {
  void dgesv_(int *, int *, double [], int *, int [], double[], int *, int *) ;

}

// Algo du papier de 1995

void Et_rot_mag::magnet_comput(Cmp (*f_j)(const Cmp& x,const double a_j), 
			       Param& par_poisson_At, 
			       Param& par_poisson_Avect){
  double relax_mag = 0.5 ;

  int Z = mp.get_mg()->get_nzone();

  // Calcul de A_0t dans l'etoile (conducteur parfait)

  Cmp A_0t(- omega * A_phi) ;
  A_0t.annule(nzet,Z-1) ;
 
  Tenseur ATTENS(A_t) ; 
  Tenseur APTENS(A_phi) ;
  Tenseur BMN(-logn) ;
  BMN =  BMN + log(bbb) ;
  BMN.set_std_base() ;


  Cmp grad1(flat_scalar_prod_desal(ATTENS.gradient_spher(),nphi.gradient_spher())());
  Cmp grad2(flat_scalar_prod_desal(APTENS.gradient_spher(),nphi.gradient_spher())()) ;
  Cmp grad3(flat_scalar_prod_desal(ATTENS.gradient_spher(),BMN.gradient_spher())()+2*nphi()*flat_scalar_prod_desal(APTENS.gradient_spher(),BMN.gradient_spher())()) ;

  Cmp ATANT(A_phi.srdsdt()); // Constrction par copie pour mapping

  ATANT.va = ATANT.va.mult_ct().ssint() ;

  Cmp ttnphi(tnphi()) ;
  ttnphi.mult_rsint() ;
  Cmp BLAH(- b_car()/(nnn()*nnn())*ttnphi*grad1)  ;
  BLAH -= (1+b_car()/(nnn()*nnn())*tnphi()*tnphi())*grad2  ;
  Cmp nphisr(nphi()) ;
  nphisr.div_r() ;
  Cmp npgrada(2*nphisr*(A_phi.dsdr()+ATANT )) ;
  npgrada.inc2_dzpuis() ;
  BLAH -=  grad3 + npgrada ;
  Cmp gtt(-nnn()*nnn()+b_car()*tnphi()*tnphi()) ;
  Cmp gtphi( - b_car()*ttnphi) ;

  // Calcul de j_t grace a Maxwell-Gauss

  j_t = ((BLAH - A_0t.laplacien())/a_car() - gtphi*j_phi)
    / gtt;
  j_t.annule(nzet, Z-1) ;

  // Calcul du courant j_phi

  j_phi = omega * j_t + (ener() + press())*f_j(A_phi, a_j) ;
  j_phi.std_base_scal() ;

  // Resolution de Maxwell Ampere (-> A_phi)
  // Calcul des termes sources avec A-t du pas precedent.

  Cmp grad4(flat_scalar_prod_desal(APTENS.gradient_spher(),
				   BMN.gradient_spher())());

  Tenseur source_tAphi(mp, 1, CON, mp.get_bvect_spher()) ;

  source_tAphi.set_etat_qcq() ;
  Cmp tjphi(j_phi) ;
  tjphi.mult_rsint() ;
  Cmp tgrad1(grad1) ;
  tgrad1.mult_rsint() ;
  Cmp d_grad4(grad4) ;
  d_grad4.div_rsint() ;
  source_tAphi.set(0)=0 ;
  source_tAphi.set(1)=0 ;
  source_tAphi.set(2)= -b_car()*a_car()*(tjphi-tnphi()*j_t)
    + b_car()/(nnn()*nnn())*(tgrad1+tnphi()*grad2)+d_grad4 ;

  source_tAphi.change_triad(mp.get_bvect_cart());

  Tenseur WORK_VECT(mp, 1, CON, mp.get_bvect_cart()) ;
    WORK_VECT.set_etat_qcq() ;
    for (int i=0; i<3; i++) {
      WORK_VECT.set(i) = 0 ; 
    }
  Tenseur WORK_SCAL(mp) ;
  WORK_SCAL.set_etat_qcq() ;
  WORK_SCAL.set() = 0 ;
  
  double lambda_mag = 0. ; // No 3D version !

  Tenseur AVECT(source_tAphi) ;
  if (source_tAphi.get_etat() != ETATZERO) {
    
    for (int i=0; i<3; i++) {
      if(source_tAphi(i).dz_nonzero()) {
	assert( source_tAphi(i).get_dzpuis() == 4 ) ; 
      }
      else{
	(source_tAphi.set(i)).set_dzpuis(4) ; 
      }
    }
    
  }
  source_tAphi.poisson_vect(lambda_mag, par_poisson_Avect, AVECT, WORK_VECT,
			    WORK_SCAL) ;
  AVECT.change_triad(mp.get_bvect_spher());
  Cmp A_phi_n(AVECT(2));
  A_phi_n.mult_rsint() ;

  // Resolution de Maxwell-Ampere : A_1

  Cmp source_A_1t(-a_car()*(j_t*gtt + j_phi*gtphi) + BLAH);

  Cmp A_1t(mp);
  A_1t = 0 ;
  source_A_1t.poisson(par_poisson_At, A_1t) ;

  int L = mp.get_mg()->get_nt(0);

  Tbl MAT(L,L) ;
  Tbl VEC(L) ;

  MAT.set_etat_qcq() ;
  VEC.set_etat_qcq() ;

  mp.tet.fait() ;
  Mtbl* theta = mp.tet.c ;
  const Map_radial* mpr = dynamic_cast<const Map_radial*>(&mp) ;
  assert (mpr != 0x0) ;
  Tbl leg(L,2*L) ;
  leg.set_etat_qcq() ;
  for(int k=0;k<L;k++){

    for(int l=0;l<2*L;l++){

      // leg[k,l] : legendre_l(cos(theta_k))
      // Construction par recurrence de degre 2

      if(l==0) leg.set(k,l)=1. ;
      if(l==1) leg.set(k,l)=cos((*theta)(l_surf()(0,k),0,k,0)) ;
      if(l>=2)
	leg.set(k,l)= double(2*l-1)/double(l)*cos((*theta)(l_surf()(0,k),0,k,0))
	  * leg(k,l-1)-double(l-1)/double(l)*leg(k,l-2) ;
    }
  }
  
  for(int k=0;k<L;k++){
    // Rsurf retourne la valeur du rayon en theta_k. via xi et ??
    double Rsurf = mpr->val_r_jk(l_surf()(0,k), xi_surf()(0,k), k, 0) ;

    // Valeurs a la surface trouvees via va.val_point_jk(l,xisurf,k,0)

    VEC.set(k) = A_0t.va.val_point_jk(l_surf()(0,k), xi_surf()(0,k), k, 0)
      -A_1t.va.val_point_jk(l_surf()(0,k), xi_surf()(0,k), k, 0);
    for(int l=0;l<L;l++){

      MAT.set(l,k) = leg(k,2*l)/pow(Rsurf,2*l+1);

    }
  }
  // appel fortran : 

  int* IPIV=new int[L] ;
  int INFO ;

  Tbl MAT_SAVE(MAT) ;
  Tbl VEC2(L) ;
  VEC2.set_etat_qcq() ;
  int un = 1 ;

  dgesv_(&L, &un, MAT.t, &L, IPIV, VEC.t, &L, &INFO) ; 

  // coeffs a_l dans VEC

  for(int k=0;k<L;k++) {VEC2.set(k)=1. ; }

  dgesv_(&L, &un, MAT_SAVE.t, &L, IPIV, VEC2.t, &L, &INFO) ;

  delete [] IPIV ;

  Cmp psi(mp);
  Cmp psi2(mp);
  psi.allocate_all() ;
  psi2.allocate_all() ;

  mp.r.fait() ;

  for(int nz=0;nz < Z; nz++){
    for(int i=0;i< mp.get_mg()->get_nr(nz);i++){
      for(int k=0;k<L;k++){
	psi.set(nz,0,k,i) = 0. ;
	psi2.set(nz,0,k,i) = 0. ;
	for(int l=0;l<L;l++){
	  psi.set(nz,0,k,i) += VEC(l)*leg(k,2*l) / 
	    pow((*mp.r.c)(nz,0,k,i),2*l+1);
	  psi2.set(nz,0,k,i) += VEC2(l)*leg(k,2*l)/
	    pow((*mp.r.c)(nz, 0, k,i),2*l+1);
	}
      }
    }
  }
  psi.annule(0) ; 
  psi2.annule(0) ;
  psi.std_base_scal() ;
  psi2.std_base_scal() ;
  
  assert(psi.get_dzpuis() == 0) ;
  int dif = A_1t.get_dzpuis() ;
  if (dif > 0) {
    for (int d=0; d<dif; d++) A_1t.dec_dzpuis() ;
  }

  Cmp A_t_ext(A_1t + psi) ;
  A_t_ext.annule(0,nzet-1) ;
  A_0t += A_t_ext ;
  A_0t.std_base_scal() ;

  Valeur** asymp = A_0t.asymptot(1) ;

  double Q_0 = -4*M_PI*(*asymp[1])(Z-1,0,0,0) ; // utilise A_0t plutot que E
  delete asymp[0] ;
  delete asymp[1] ;

  delete [] asymp ;

  asymp = psi2.asymptot(1) ;

  double Q_2 = -4*M_PI*(*asymp[1])(Z-1,0,0,0)  ; // A_2t = psi2 a l'infini
  delete asymp[0] ;
  delete asymp[1] ;

  delete [] asymp ;

  // solution definitive :

  double C = (Q-Q_0)/Q_2 ;

  assert(psi2.get_dzpuis() == 0) ;
  dif = A_0t.get_dzpuis() ;
  if (dif > 0) {
    for (int d=0; d<dif; d++) A_0t.dec_dzpuis() ;
  }
  Cmp A_t_n(A_0t + C) ;
  A_t_ext = A_0t + C*psi2 ;
  A_t_ext.annule(0,nzet-1) ;
  A_t_n.annule(nzet,Z-1) ;
  A_t_n += A_t_ext ;

  A_t_n.std_base_scal() ;

  asymp = A_t_n.asymptot(1) ;

  delete asymp[0] ;
  delete asymp[1] ;

  delete [] asymp ;
  A_t = relax_mag*A_t_n + (1.-relax_mag)*A_t ;
  A_phi = relax_mag*A_phi_n + (1. - relax_mag)*A_phi ;

}



