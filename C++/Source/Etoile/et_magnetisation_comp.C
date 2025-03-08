/*
 * Computational functions for magnetized rotating equilibrium
 *
 * (see file et_rot_mag.h for documentation)
 *
 */

/*
 *   Copyright (c) 2013 Debarati Chatterjee, Jerome Novak
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

 

/*
 * $Id$
 * $Log$
 * Revision 1.15  2016/12/05 16:17:53  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.14  2016/11/01 09:12:59  j_novak
 * Correction of a missing '-' in mom_quad_old().
 *
 * Revision 1.13  2015/06/12 12:38:25  j_novak
 * Implementation of the corrected formula for the quadrupole momentum.
 *
 * Revision 1.12  2014/10/21 09:23:54  j_novak
 * Addition of global functions mass_g(), angu_mom(), grv2/3() and mom_quad().
 *
 * Revision 1.11  2014/10/13 08:52:57  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.10  2014/07/04 12:08:02  j_novak
 * Added some filtering.
 *
 * Revision 1.9  2014/05/14 15:19:05  j_novak
 * The magnetisation field is now filtered.
 *
 * Revision 1.8  2014/05/13 15:37:12  j_novak
 * Updated to new magnetic units.
 *
 * Revision 1.7  2014/05/01 13:07:16  j_novak
 * Fixed two bugs: in the computation of F31,F32 and the triad of U_up.
 *
 * Revision 1.6  2014/04/29 13:46:07  j_novak
 * Addition of switches 'use_B_in_eos' and 'include_magnetisation' to control the model.
 *
 * Revision 1.5  2014/04/28 14:53:29  j_novak
 * Minor modif.
 *
 * Revision 1.4  2014/04/28 12:48:13  j_novak
 * Minor modifications.
 *
 * Revision 1.2  2013/12/19 17:05:40  j_novak
 * Corrected a dzpuis problem.
 *
 * Revision 1.1  2013/12/13 16:36:51  j_novak
 * Addition and computation of magnetisation terms in the Einstein equations.
 *
 *
 *
 * $Header$
 *
 */

// Headers C
#include <cstdlib>
#include <cmath>

// Headers Lorene
#include "et_rot_mag.h"
#include "metric.h"
#include "utilitaires.h"
#include "param.h"
#include "proto_f77.h"
#include "unites.h"

namespace Lorene {

  using namespace Unites_mag ;

// Algo du papier de 1995

void Et_magnetisation::magnet_comput(const int adapt_flag,
			       Cmp (*f_j)(const Cmp&, const double), 
			       Param& par_poisson_At, 
			       Param& par_poisson_Avect){
  double relax_mag = 0.5 ;

  int Z = mp.get_mg()->get_nzone();

  bool adapt(adapt_flag) ;
  /****************************************************************
   *  Assertion that all zones have same number of points in theta
   ****************************************************************/
  int nt = mp.get_mg()->get_nt(nzet-1) ; 
  for (int l=0; l<Z; l++) assert(mp.get_mg()->get_nt(l) == nt) ;
  
  Tbl Rsurf(nt) ;
  Rsurf.set_etat_qcq() ;
  mp.r.fait() ;
  mp.tet.fait() ;
  Mtbl* theta = mp.tet.c ;
  const Map_radial* mpr = dynamic_cast<const Map_radial*>(&mp) ;
  assert (mpr != 0x0) ;
  for (int j=0; j<nt; j++) 
    Rsurf.set(j) = mpr->val_r_jk(l_surf()(0,j), xi_surf()(0,j), j, 0) ;
  
  
  // Calcul de A_0t dans l'etoile (conducteur parfait)
  
  Cmp A_0t(- omega * A_phi) ;
  A_0t.annule(nzet,Z-1) ;
  
  Tenseur ATTENS(A_t) ; 
  Tenseur APTENS(A_phi) ;
  Tenseur BMN(-logn) ;
  BMN =  BMN + log(bbb) ;
  BMN.set_std_base() ;
  
  
  Cmp grad1(flat_scalar_prod_desal(ATTENS.gradient_spher(),
				   nphi.gradient_spher())());
  Cmp grad2(flat_scalar_prod_desal(APTENS.gradient_spher(),
				   nphi.gradient_spher())()) ;
  Cmp grad3(flat_scalar_prod_desal(ATTENS.gradient_spher(),
				   BMN.gradient_spher())()
	    + 2*nphi()*flat_scalar_prod_desal(APTENS.gradient_spher(),
					      BMN.gradient_spher())()) ;
  
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
  
  // Computation of j_t thanks to Maxwell-Gauss
  // modified to include Magnetisation
  // components of F
  Cmp F01 = 1/(a_car()*nnn()*nnn())*A_0t.dsdr() 
        +   1/(a_car()*nnn()*nnn())*nphi()*A_phi.dsdr() ;

  Cmp F02 =  1/(a_car()*nnn()*nnn())*A_0t.srdsdt()
        +   1/(a_car()*nnn()*nnn())*nphi()*A_phi.srdsdt() ;

  Cmp tmp = A_phi.dsdr() / (bbb() * bbb() * a_car() );
  tmp.div_rsint() ;
  tmp.div_rsint() ;
  Cmp F31 = 1/(a_car()*nnn()*nnn())*nphi()*nphi()*A_phi.dsdr()
    +    1/(a_car()*nnn()*nnn())*nphi()*A_0t.dsdr()  
    +    tmp ;

  tmp = A_phi.srdsdt() / (bbb() * bbb() * a_car() );
  tmp.div_rsint() ;
  tmp.div_rsint() ;
  Cmp F32 = 1/(a_car()*nnn()*nnn())*nphi()*nphi()*A_phi.srdsdt()
    +    1/(a_car()*nnn()*nnn())*nphi()*A_0t.srdsdt()
    +    tmp ;

  Cmp x = get_magnetisation();
  Cmp one_minus_x = 1 - x ;
  one_minus_x.std_base_scal() ;

  tmp = ((BLAH - A_0t.laplacien())*one_minus_x/a_car()
	   - gtphi*j_phi 
	   - gtt*(F01*x.dsdr()+F02*x.srdsdt())
	   - gtphi*(F31*x.dsdr()+F32*x.srdsdt()) ) / gtt ;
 
  tmp.annule(nzet, Z-1) ;
  if (adapt) {
    j_t = tmp ;
  }
  else {
    j_t.allocate_all() ;
    for (int j=0; j<nt; j++) 
      for (int l=0; l<nzet; l++) 
	for (int i=0; i<mp.get_mg()->get_nr(l); i++) 
	  j_t.set(l,0,j,i) = ( (*mp.r.c)(l,0,j,i) > Rsurf(j) ?
			       0. : tmp(l,0,j,i) ) ;
    j_t.annule(nzet,Z-1) ;
  }
  j_t.std_base_scal() ;
  
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

// modified to include Magnetisation
  Cmp phifac = (F31-nphi()*F01)*x.dsdr()
             + (F32-nphi()*F02)*x.srdsdt() ;
  phifac.mult_rsint();
  source_tAphi.set(2)= -b_car()*a_car()/one_minus_x
      *(tjphi-tnphi()*j_t + phifac)
    + b_car()/(nnn()*nnn())*(tgrad1+tnphi()*grad2)
    + d_grad4 ;
  
  source_tAphi.change_triad(mp.get_bvect_cart());

  // Filtering
  for (int i=0; i<3; i++) {
    Scalar tmp_filter = source_tAphi(i) ;
    tmp_filter.exponential_filter_r(0, 2, 1) ;
    tmp_filter.exponential_filter_ylm(0, 2, 1) ;
    source_tAphi.set(i) = tmp_filter ;
  }

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
  
  // Solution to Maxwell-Ampere : A_1
  // modified to include Magnetisation
  Cmp source_A_1t(-a_car()*( j_t*gtt + j_phi*gtphi 
   + gtt*(F01*x.dsdr()+F02*x.srdsdt())
   + gtphi*(F31*x.dsdr()+F32*x.srdsdt()) )/one_minus_x 
   + BLAH);
  Scalar tmp_filter = source_A_1t ;
  tmp_filter.exponential_filter_r(0, 2, 1) ;
  tmp_filter.exponential_filter_ylm(0, 2, 1) ;
  source_A_1t = tmp_filter ;
  
  Cmp A_1t(mp);
  A_1t = 0 ;
  source_A_1t.poisson(par_poisson_At, A_1t) ;

  int L = mp.get_mg()->get_nt(0);
  
  Tbl MAT(L,L) ;
  Tbl MAT_PHI(L,L);
  Tbl VEC(L) ;
  
  MAT.set_etat_qcq() ;
  VEC.set_etat_qcq() ;
  MAT_PHI.set_etat_qcq() ;
  
  Tbl leg(L,2*L) ;
  leg.set_etat_qcq() ;
  
  Cmp psi(mp);
  Cmp psi2(mp);
  psi.allocate_all() ;
  psi2.allocate_all() ;
  
  for (int p=0; p<mp.get_mg()->get_np(0); p++) {
    // leg[k,l] : legendre_l(cos(theta_k))
    // Construction par recurrence de degre 2
    for(int k=0;k<L;k++){
      for(int l=0;l<2*L;l++){
	
	if(l==0) leg.set(k,l)=1. ;
	if(l==1) leg.set(k,l)=cos((*theta)(l_surf()(p,k),p,k,0)) ;
	if(l>=2) leg.set(k,l) = double(2*l-1)/double(l)  
		   * cos((*theta)(l_surf()(p,k),p,k,0))
		   * leg(k,l-1)-double(l-1)/double(l)*leg(k,l-2) ;
      }
    }
    
    for(int k=0;k<L;k++){
      
      // Valeurs a la surface trouvees via va.val_point_jk(l,xisurf,k,p)
      
      VEC.set(k) = A_0t.va.val_point_jk(l_surf()(p,k), xi_surf()(p,k), k, p)
	-A_1t.va.val_point_jk(l_surf()(p,k), xi_surf()(p,k), k, p);
      
      for(int l=0;l<L;l++) MAT.set(l,k) = leg(k,2*l)/pow(Rsurf(k),2*l+1);
      
    }
    // appel fortran : 
    
    int* IPIV=new int[L] ;
    int INFO ;
    
    Tbl MAT_SAVE(MAT) ;
    Tbl VEC2(L) ;
    VEC2.set_etat_qcq() ;
    int un = 1 ;
    
    F77_dgesv(&L, &un, MAT.t, &L, IPIV, VEC.t, &L, &INFO) ; 
    
    // coeffs a_l dans VEC
    
    for(int k=0;k<L;k++) {VEC2.set(k)=1. ; }
    
    F77_dgesv(&L, &un, MAT_SAVE.t, &L, IPIV, VEC2.t, &L, &INFO) ;
    
    delete [] IPIV ;
    
    for(int nz=0;nz < Z; nz++){
      for(int i=0;i< mp.get_mg()->get_nr(nz);i++){
	for(int k=0;k<L;k++){
	  psi.set(nz,p,k,i) = 0. ;
	  psi2.set(nz,p,k,i) = 0. ;
	  for(int l=0;l<L;l++){
	    psi.set(nz,p,k,i) += VEC(l)*leg(k,2*l) / 
	      pow((*mp.r.c)(nz,p,k,i),2*l+1);
	    psi2.set(nz,p,k,i) += VEC2(l)*leg(k,2*l)/
	      pow((*mp.r.c)(nz, p, k,i),2*l+1);
	  }
	}
      }
    }
  }
  psi.std_base_scal() ;
  psi2.std_base_scal() ;
  
  assert(psi.get_dzpuis() == 0) ;
  int dif = A_1t.get_dzpuis() ;
  if (dif > 0) {
    for (int d=0; d<dif; d++) A_1t.dec_dzpuis() ;
  }
  
  if (adapt) {
    Cmp A_t_ext(A_1t + psi) ;
    A_t_ext.annule(0,nzet-1) ;
    A_0t += A_t_ext ;
  }
  else {
    tmp = A_0t ;
    A_0t.allocate_all() ;
    for (int j=0; j<nt; j++) 
      for (int l=0; l<Z; l++) 
	for (int i=0; i<mp.get_mg()->get_nr(l); i++) 
	  A_0t.set(l,0,j,i) = ( (*mp.r.c)(l,0,j,i) > Rsurf(j) ?
				A_1t(l,0,j,i) + psi(l,0,j,i) : tmp(l,0,j,i) ) ;
  } 
  A_0t.std_base_scal() ;
  
  tmp_filter = A_0t ;
  tmp_filter.exponential_filter_r(0, 2, 1) ;
  tmp_filter.exponential_filter_ylm(0, 2, 1) ;
  A_0t = tmp_filter ;

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
  
  // solution definitive de A_t:
  
  double C = (Q-Q_0)/Q_2 ;
  
  assert(psi2.get_dzpuis() == 0) ;
  dif = A_0t.get_dzpuis() ;
  if (dif > 0) {
    for (int d=0; d<dif; d++) A_0t.dec_dzpuis() ;
  }
  Cmp A_t_n(mp) ;
  if (adapt) {
    A_t_n = A_0t + C ;
    Cmp A_t_ext(A_0t + C*psi2) ;
    A_t_ext.annule(0,nzet-1) ;
    A_t_n.annule(nzet,Z-1) ;
    A_t_n += A_t_ext ;
  }
  else {
    A_t_n.allocate_all() ;
    for (int j=0; j<nt; j++) 
      for (int l=0; l<Z; l++) 
	for (int i=0; i<mp.get_mg()->get_nr(l); i++) 
	  A_t_n.set(l,0,j,i) = ( (*mp.r.c)(l,0,j,i) > Rsurf(j) ?
				 A_0t(l,0,j,i) + C*psi2(l,0,j,i) : 
				 A_0t(l,0,j,i) + C ) ;    
  }
  A_t_n.std_base_scal() ;
  tmp_filter = A_t_n ;
  tmp_filter.exponential_filter_r(0, 2, 1) ;
  tmp_filter.exponential_filter_ylm(0, 2, 1) ;
  A_t_n = tmp_filter ;
  
  asymp = A_t_n.asymptot(1) ;
  
  delete asymp[0] ;
  delete asymp[1] ;
  
  delete [] asymp ;
  A_t = relax_mag*A_t_n + (1.-relax_mag)*A_t ;
  A_phi = relax_mag*A_phi_n + (1. - relax_mag)*A_phi ;

}


void Et_magnetisation::MHD_comput() {
  // Computes the E-M terms of the stress-energy tensor...
  
  Tenseur ATTENS(A_t) ;

  Tenseur APTENS(A_phi) ;
  
  Tenseur ApAp ( flat_scalar_prod_desal(APTENS.gradient_spher(),
					APTENS.gradient_spher())() );
  Tenseur ApAt ( flat_scalar_prod_desal(APTENS.gradient_spher(),
					ATTENS.gradient_spher())() );
  Tenseur AtAt ( flat_scalar_prod_desal(ATTENS.gradient_spher(),
					ATTENS.gradient_spher())() );

  if (ApAp.get_etat() != ETATZERO) {
    ApAp.set().div_rsint() ;
    ApAp.set().div_rsint() ;
  }
  if (ApAt.get_etat() != ETATZERO) 
    ApAt.set().div_rsint() ;
  
  E_em = 0.5*mu0 * ( 1/(a_car*nnn*nnn) * (AtAt + 2*tnphi*ApAt)
	      + ( (tnphi*tnphi/(a_car*nnn*nnn)) + 1/(a_car*b_car) )*ApAp );
  Jp_em = -mu0 * (ApAt + tnphi*ApAp) /(a_car*nnn) ;
  if (Jp_em.get_etat() != ETATZERO) Jp_em.set().mult_rsint() ;
  Srr_em = 0 ;
  // Stt_em = -Srr_em
  Spp_em = E_em ;

  // ... and those corresponding to the magnetization.
  Tenseur Efield = Elec() ;
  Tenseur Bfield = Magn() ;

  Scalar EiEi ( flat_scalar_prod(Efield, Efield)() ) ;
  Scalar BiBi ( flat_scalar_prod(Bfield, Bfield)() ) ;

  Vector U_up(mp, CON, mp.get_bvect_cart()) ;
  for (int i=1; i<=3; i++) 
    U_up.set(i) = u_euler(i-1) ;
  U_up.change_triad(mp.get_bvect_spher()) ;

  Sym_tensor gamij(mp, COV, mp.get_bvect_spher()) ;
  for (int i=1; i<=3; i++)
    for (int j=1; j<i; j++) {
      gamij.set(i,j) = 0 ;
    }
  gamij.set(1,1) = a_car() ;
  gamij.set(2,2) = a_car() ;
  gamij.set(3,3) = b_car() ;
  Metric met(gamij) ;
  Vector Ui = U_up.down(0, met) ;

  Scalar fac = sqrt(a_car()) ;
  Vector B_up(mp, CON, mp.get_bvect_spher()) ;
  B_up.set(1) = Scalar(Bfield(0)) / fac ;
  B_up.set(2) = Scalar(Bfield(1)) / fac ;
  B_up.set(3) = 0 ;
  Vector Bi = B_up.down(0, met) ;

  fac = Scalar(gam_euler()*gam_euler()) ;

  E_I = get_magnetisation() * EiEi / mu0 ;

  J_I = get_magnetisation() * BiBi * Ui / mu0 ;
  Sij_I = get_magnetisation() 
    * ( (BiBi / fac) * gamij  + BiBi*Ui*Ui - Bi*Bi / fac ) / mu0 ;

  for (int i=1; i<=3; i++)
    for (int j=i; j<=3; j++)
      Sij_I.set(i,j).set_dzpuis(0) ;

}

			//----------------------------//
			//	Gravitational mass    //
			//----------------------------//

double Et_magnetisation::mass_g() const {

  if (p_mass_g == 0x0) {    // a new computation is required
    
    if (relativistic) {
      
      // Magnetisation: S_{rr} + S_{\theta\theta}
      Tenseur SrrplusStt( Cmp(Sij_I(1, 1) + Sij_I(2, 2)) ) ; 
      SrrplusStt = SrrplusStt / a_car ; // S^r_r + S^\theta_\theta
      
      Tenseur Spp (Cmp(Sij_I(3, 3))) ; // Magnetisation: S_{\phi\phi}
      Spp = Spp / b_car ; // S^\phi_\phi

      Cmp temp(E_I) ;
      Tenseur E_i (temp) ;
      Tenseur J_i (Cmp(J_I(3))) ;
      
      Tenseur source = nnn * (ener_euler + E_em + E_i 
			      + s_euler + Spp_em + SrrplusStt + Spp) + 
	nphi * (Jp_em + J_i) 
	+ 2 * bbb * (ener_euler + press) * tnphi * uuu ;
      
      source = a_car * bbb * source ;
      
      source.set_std_base() ;
      
      p_mass_g = new double( source().integrale() ) ;


    }
    else{  // Newtonian case 
      p_mass_g = new double( mass_b() ) ;   // in the Newtonian case
      //  M_g = M_b
    }
  }
    
  return *p_mass_g ; 

} 

			//----------------------------//
			//	Angular momentum      //
			//----------------------------//

double Et_magnetisation::angu_mom() const {

  if (p_angu_mom == 0x0) {    // a new computation is required
	
    Cmp dens = uuu() ; 

    dens.mult_r() ;			//  Multiplication by
    dens.va = (dens.va).mult_st() ;	//    r sin(theta)

    if (relativistic) {
      dens = a_car() * (b_car() * (ener_euler() + press()) 
			* dens + bbb() * (Jp_em() + Cmp(J_I(3)) ) ) ; 
    }
    else {    // Newtonian case 
      dens = nbar() * dens ; 
    }

    dens.std_base_scal() ; 

    p_angu_mom = new double( dens.integrale() ) ;
    
  }
    
  return *p_angu_mom ; 

}

			//----------------------------//
			//	     GRV2	      //
			//----------------------------//

double Et_magnetisation::grv2() const {

    if (p_grv2 == 0x0) {    // a new computation is required
	
      // To get qpig:	
      using namespace Unites ;

      Tenseur Spp (Cmp(Sij_I(3, 3))) ; //S_{\phi\phi}
      Spp = Spp / b_car ; // S^\phi_\phi

      Tenseur sou_m =  2 * qpig * a_car * (press + (ener_euler+press)
					   * uuu*uuu + Spp) ;
      
      Tenseur sou_q =   2 * qpig * a_car * Spp_em + 1.5 * ak_car
	- flat_scalar_prod(logn.gradient_spher(), logn.gradient_spher() ) ;

      p_grv2 = new double( double(1) - lambda_grv2(sou_m(), sou_q()) ) ; 

    }
    
    return *p_grv2 ; 

}


			//----------------------------//
			//	     GRV3	      //
			//----------------------------//

double Et_magnetisation::grv3(ostream* ost) const {

    if (p_grv3 == 0x0) {    // a new computation is required

	// To get qpig:	
      using namespace Unites ;
      
      Tenseur source(mp) ; 
	
	// Gravitational term [cf. Eq. (43) of Gourgoulhon & Bonazzola
	// ------------------	    Class. Quantum Grav. 11, 443 (1994)]

	if (relativistic) {
	    Tenseur alpha = dzeta - logn ; 
	    Tenseur beta = log( bbb ) ; 
	    beta.set_std_base() ; 
	    
	    source = 0.75 * ak_car 
		     - flat_scalar_prod(logn.gradient_spher(),
					logn.gradient_spher() )
		     + 0.5 * flat_scalar_prod(alpha.gradient_spher(),
					      beta.gradient_spher() ) ; 
	    
	    Cmp aa = alpha() - 0.5 * beta() ; 
	    Cmp daadt = aa.srdsdt() ;	// 1/r d/dth
	    
	    // What follows is valid only for a mapping of class Map_radial : 
	    const Map_radial* mpr = dynamic_cast<const Map_radial*>(&mp) ; 
	    if (mpr == 0x0) {
		cout << "Etoile_rot::grv3: the mapping does not belong"
		     << " to the class Map_radial !" << endl ; 
		abort() ; 
	    }
		
	    // Computation of 1/tan(theta) * 1/r daa/dtheta
	    if (daadt.get_etat() == ETATQCQ) {
		Valeur& vdaadt = daadt.va ; 
		vdaadt = vdaadt.ssint() ;	// division by sin(theta)
		vdaadt = vdaadt.mult_ct() ;	// multiplication by cos(theta)
	    }
	    
	    Cmp temp = aa.dsdr() + daadt ; 
	    temp = ( bbb() - a_car()/bbb() ) * temp ; 
	    temp.std_base_scal() ; 
	    
	    // Division by r 
	    Valeur& vtemp = temp.va ; 
	    vtemp = vtemp.sx() ;    // division by xi in the nucleus
				    // Id in the shells
				    // division by xi-1 in the ZEC
	    vtemp = (mpr->xsr) * vtemp ; // multiplication by xi/r in the nucleus
					 //		  by 1/r in the shells
					 //		  by r(xi-1) in the ZEC

	    // In the ZEC, a multiplication by r has been performed instead
	    //   of the division: 			
	    temp.set_dzpuis( temp.get_dzpuis() + 2 ) ;  
	    
	    source = bbb() * source() + 0.5 * temp ; 

	}
	else{
	    source = - 0.5 * flat_scalar_prod(logn.gradient_spher(),
					      logn.gradient_spher() ) ; 
	}
	
	source.set_std_base() ; 

	double int_grav = source().integrale() ; 

	// Matter term
	// -----------
	
	if (relativistic) {   
 
	  // S_{rr} + S_{\theta\theta}
	  Tenseur SrrplusStt( Cmp(Sij_I(1, 1) + Sij_I(2, 2)) ) ; 
	  SrrplusStt = SrrplusStt / a_car ; // S^r_r + S^\theta_\theta
	  
	  Tenseur Spp (Cmp(Sij_I(3, 3))) ; //S_{\phi\phi}
	  Spp = Spp / b_car ; // S^\phi_\phi
  
	    source  = qpig * a_car * bbb * ( s_euler + Spp_em + SrrplusStt + Spp ) ;
	}
	else{
	    source = qpig * ( 3 * press + nbar * uuu * uuu ) ; 
	}

	source.set_std_base() ; 

	double int_mat = source().integrale() ; 

	// Virial error
	// ------------
	if (ost != 0x0) {
	    *ost << "Et_magnetisation::grv3 : gravitational term : " << int_grav 
		 << endl ;
	    *ost << "Et_magnetisation::grv3 : matter term :        " << int_mat 
		 << endl ;
	}

	p_grv3 = new double( (int_grav + int_mat) / int_mat ) ; 	 

    }
    
    return *p_grv3 ; 

}

			//----------------------------//
			//     Quadrupole moment      //
			//----------------------------//

  double Et_magnetisation::mom_quad_old() const {
    
    if (p_mom_quad_old == 0x0) {    // a new computation is required
      
      // To get qpig:	
      using namespace Unites ;
      
      // Source for of the Poisson equation for nu
      // -----------------------------------------
      
      Tenseur source(mp) ; 
      
      if (relativistic) {
	// S_{rr} + S_{\theta\theta}
	Tenseur SrrplusStt( Cmp(Sij_I(1, 1) + Sij_I(2, 2)) ) ; 
	SrrplusStt = SrrplusStt / a_car ; // S^r_r + S^\theta_\theta
	
	Tenseur Spp (Cmp(Sij_I(3, 3))) ; //S_{\phi\phi}
	Spp = Spp / b_car ; // S^\phi_\phi
	
	Cmp temp(E_I) ;
	Tenseur E_i(temp) ;
	
	Tenseur beta = log(bbb) ; 
	beta.set_std_base() ; 
	source =  qpig * a_car *( ener_euler + E_em + E_i
				  + s_euler + Spp_em + SrrplusStt + Spp) 
	  + ak_car - flat_scalar_prod(logn.gradient_spher(), 
		     logn.gradient_spher() + beta.gradient_spher()) ; 
      }
      else {
	source = qpig * nbar ; 
      }
      source.set_std_base() ; 	
      
      // Multiplication by -r^2 P_2(cos(theta))
      //  [cf Eq.(7) of Salgado et al. Astron. Astrophys. 291, 155 (1994) ]
      // ------------------------------------------------------------------
      
      // Multiplication by r^2 : 
      // ----------------------
      Cmp& csource = source.set() ; 
      csource.mult_r() ; 
      csource.mult_r() ; 
      if (csource.check_dzpuis(2)) {
	csource.inc2_dzpuis() ; 
      }
      
      // Muliplication by cos^2(theta) :
      // -----------------------------
      Cmp temp = csource ; 
      
      // What follows is valid only for a mapping of class Map_radial : 
      assert( dynamic_cast<const Map_radial*>(&mp) != 0x0 ) ; 
      
      if (temp.get_etat() == ETATQCQ) {
	Valeur& vtemp = temp.va ; 
	vtemp = vtemp.mult_ct() ;	// multiplication by cos(theta)
	vtemp = vtemp.mult_ct() ;	// multiplication by cos(theta)
      }
      
      // Muliplication by -P_2(cos(theta)) :
      // ----------------------------------
      source = 0.5 * source() - 1.5 * temp ; 
      
      // Final result
      // ------------
      p_mom_quad_old = new double( - source().integrale() / qpig ) ; 	  
    }
    return *p_mom_quad_old ; 
  }


  double Et_magnetisation::mom_quad_Bo() const {
    
    using namespace Unites ;
    
    if (p_mom_quad_Bo == 0x0) {    // a new computation is required

      // S_{rr} + S_{\theta\theta} =  A^2*(S^r_r + S^\theta_\theta)
      Tenseur SrrplusStt( Cmp(Sij_I(1, 1) + Sij_I(2, 2)) ) ; 
      
      Cmp dens = a_car() * press() ;
      dens = bbb() * nnn() * (SrrplusStt() + 2*dens) ; 
      dens.mult_rsint() ;
      dens.std_base_scal() ; 
      
      p_mom_quad_Bo = new double( - 16. * dens.integrale() / qpig  ) ;  
    }
    return *p_mom_quad_Bo ;  
  }



}
