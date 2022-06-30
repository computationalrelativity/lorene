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

 

/*
 * $Id$
 * $Log$
 * Revision 1.19  2022/06/30 16:24:16  j_novak
 * Corrected a bug in the matching of A_t potential in the case np>1 (loops were not correctly set).
 *
 * Revision 1.18  2016/12/05 16:17:54  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.17  2014/10/13 08:52:58  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.16  2014/09/03 15:33:42  j_novak
 * Filtering of Maxwell sources is now optional.
 *
 * Revision 1.15  2014/07/04 12:15:12  j_novak
 * Added filtering.
 *
 * Revision 1.14  2005/06/03 15:31:56  j_novak
 * Better computation when more than one point in phi.
 *
 * Revision 1.13  2003/10/03 15:58:47  j_novak
 * Cleaning of some headers
 *
 * Revision 1.12  2002/09/09 13:00:39  e_gourgoulhon
 * Modification of declaration of Fortran 77 prototypes for
 * a better portability (in particular on IBM AIX systems):
 * All Fortran subroutine names are now written F77_* and are
 * defined in the new file C++/Include/proto_f77.h.
 *
 * Revision 1.11  2002/06/05 15:15:59  j_novak
 * The case of non-adapted mapping is treated.
 * parmag.d and parrot.d have been merged.
 *
 * Revision 1.10  2002/06/03 13:23:16  j_novak
 * The case when the mapping is not adapted is now treated
 *
 * Revision 1.9  2002/06/03 13:00:45  e_marcq
 *
 * conduc parameter read in parmag.d
 *
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
#include <cstdlib>
#include <cmath>

// Headers Lorene
#include "et_rot_mag.h"
#include "utilitaires.h"
#include "param.h"
#include "proto_f77.h"
#include "graphique.h"
#include "tensor.h"

namespace Lorene {
// Local prototype (for drawings only)
Cmp raccord_c1(const Cmp& uu, int l1) ; 

// Algo du papier de 1995

void Et_rot_mag::magnet_comput(const int adapt_flag,
			       Cmp (*f_j)(const Cmp&, const double), 
			       Param& par_poisson_At, 
			       Param& par_poisson_Avect){
  double relax_mag = 0.5 ;

  int mag_filter = 0 ;
  if (par_poisson_At.get_n_int() > 1)
    mag_filter = par_poisson_At.get_int(1) ;

  int Z = mp.get_mg()->get_nzone();

  if(is_conduct()) {
    bool adapt(adapt_flag) ;
    /****************************************************************
     *  Assertion that all zones have same number of points in theta
     ****************************************************************/
    int nt = mp.get_mg()->get_nt(nzet-1) ; 
    int np = mp.get_mg()->get_np(nzet-1) ; 
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
    //A_0t.annule(nzet,Z-1) ;
 
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

    // Calcul de j_t grace a Maxwell-Gauss
    Cmp tmp(((BLAH - A_0t.laplacien())/a_car() - gtphi*j_phi)
	    / gtt);
    tmp.annule(nzet, Z-1) ;
    if (adapt) {
      j_t = tmp ;
    }
    else {
      j_t.allocate_all() ;
      for (int k=0; k<np; k++) 
	for (int j=0; j<nt; j++) 
	  for (int l=0; l<nzet; l++) 
	    for (int i=0; i<mp.get_mg()->get_nr(l); i++) 
	      j_t.set(l,k,j,i) = ( (*mp.r.c)(l,k,j,i) > Rsurf(j) ?
				   0. : tmp(l,k,j,i) ) ;
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
    source_tAphi.set(2)= -b_car()*a_car()*(tjphi-tnphi()*j_t)
      + b_car()/(nnn()*nnn())*(tgrad1+tnphi()*grad2)+d_grad4 ;

    source_tAphi.change_triad(mp.get_bvect_cart());
    if (mag_filter == 1) {
      for (int i=0; i<3; i++) {
	Scalar tmp_filter = source_tAphi(i) ;
	tmp_filter.exponential_filter_r(0, 2, 1) ;
	tmp_filter.exponential_filter_ylm(0, 2, 1) ;
	source_tAphi.set(i) = tmp_filter ;
      }
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

    // Resolution de Maxwell-Ampere : A_1

    Cmp source_A_1t(-a_car()*(j_t*gtt + j_phi*gtphi) + BLAH);
    if (mag_filter == 1) {
      Scalar tmp_filter = source_A_1t ;
      tmp_filter.exponential_filter_r(0, 2, 1) ;
      tmp_filter.exponential_filter_ylm(0, 2, 1) ;
      source_A_1t = tmp_filter ;
    }

    Cmp A_1t(mp);
    A_1t = 0 ;

    source_A_1t.poisson(par_poisson_At, A_1t) ;

    int L = mp.get_mg()->get_nt(0) ;

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
    Cmp psi3(mp);
    psi.allocate_all() ;
    psi2.allocate_all() ;
    psi3.allocate_all() ;

    Tbl VEC3(L) ;
    VEC3.set_etat_qcq() ;
    for (int i=0; i<L; i++)
      VEC3.set(i) = 1. / double(i+1) ;

    for (int p=0; p<np; p++) {
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
		    psi3.set(nz,p,k,i) = 0. ;
		    for(int l=0;l<L;l++){
			psi.set(nz,p,k,i) += VEC(l)*leg(k,2*l) / 
			    pow((*mp.r.c)(nz,p,k,i),2*l+1);
			psi2.set(nz,p,k,i) += VEC2(l)*leg(k,2*l)/
			    pow((*mp.r.c)(nz, p, k,i),2*l+1);
			psi3.set(nz,p,k,i) += VEC3(l)*leg(k,2*l)/
			  (pow((*mp.r.c)(nz, p, k,i),2*l+1)) ;
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
      for (int k=0; k<np; k++)
	for (int j=0; j<nt; j++) 
	  for (int l=0; l<Z; l++) 
	    for (int i=0; i<mp.get_mg()->get_nr(l); i++) 
	      A_0t.set(l,k,j,i) = ( (*mp.r.c)(l,k,j,i) > Rsurf(j) ?
				  A_1t(l,k,j,i) + psi(l,k,j,i) : tmp(l,k,j,i) ) ;
    } 
    A_0t.std_base_scal() ;

   if (mag_filter == 1) {
      Scalar tmp_filter = A_0t ;
      tmp_filter.exponential_filter_r(0, 2, 1) ;
      tmp_filter.exponential_filter_ylm(0, 2, 1) ;
      A_0t = tmp_filter ;
    }

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
      for (int k=0; k<np; k++)
	for (int j=0; j<nt; j++) 
	  for (int l=0; l<Z; l++) 
	    for (int i=0; i<mp.get_mg()->get_nr(l); i++) {
	      A_t_n.set(l,k,j,i) = ( (*mp.r.c)(l,k,j,i) > Rsurf(j) ?
				   A_0t(l,k,j,i) + C*psi2(l,k,j,i) : 
				   A_0t(l,k,j,i) + C ) ;    
	  }
    }
    A_t_n.std_base_scal() ;
    if (mag_filter == 1) {
      Scalar tmp_filter = A_t_n ;
      tmp_filter.exponential_filter_r(0, 2, 1) ;
      tmp_filter.exponential_filter_ylm(0, 2, 1) ;
      A_t_n = tmp_filter ;
    }

    asymp = A_t_n.asymptot(1) ;

    delete asymp[0] ;
    delete asymp[1] ;

    delete [] asymp ;
    A_t = relax_mag*A_t_n + (1.-relax_mag)*A_t ;
    A_phi = relax_mag*A_phi_n + (1. - relax_mag)*A_phi ;

  } // End of perfect conductor case
  
  else
    {
    
    /***************
     * CAS ISOLANT *
     ***************/								    

    // Calcul de j_t
    j_t = Q*nbar() + (ener()+press())*f_j(omega* A_phi - A_t,a_j) ;
    j_t.annule(nzet,Z-1) ;
    j_t.std_base_scal() ;
    
    // Calcul de j_phi
    j_phi = omega * j_t ;
    j_phi.std_base_scal() ;
    
    // Resolution de A_t
    
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

    Cmp ATANT(A_phi.srdsdt());

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

    Cmp source_A_t_n(mp);
    if (relativistic) {
      source_A_t_n = (-a_car()*(j_t*gtt + j_phi*gtphi) + BLAH);
      source_A_t_n.std_base_scal();}
    else{
      source_A_t_n = j_t;}

    Cmp A_t_n(A_t) ;
    A_t_n = 0 ;
    A_t_n.std_base_scal() ;

    source_A_t_n.poisson(par_poisson_At, A_t_n) ;

    // Resolution de A_phi

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

    if (relativistic) {
      source_tAphi.set(2)= -b_car()*a_car()*(tjphi-tnphi()*j_t)
	+ b_car()/(nnn()*nnn())*(tgrad1+tnphi()*grad2)+d_grad4 ;}
    else{
      source_tAphi.set(2)= - tjphi ;}

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

  // Relaxation

    A_t = relax_mag*A_t_n + (1.-relax_mag)*A_t ;
    A_phi = relax_mag*A_phi_n + (1. - relax_mag)*A_phi ;

  }


}













}
