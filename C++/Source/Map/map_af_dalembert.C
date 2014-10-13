/*
 *   Copyright (c) 2000-2004 Jerome Novak
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


char map_af_dalembert_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.18  2014/10/13 08:53:02  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.17  2014/10/06 15:13:11  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.16  2008/08/27 08:55:31  jl_cornou
 * Added R_JACO02 case
 *
 * Revision 1.15  2007/11/06 14:42:20  j_novak
 * Copy of field at previous time-steps to local variables to deal with the
 * dzpuis.
 *
 * Revision 1.14  2006/08/31 08:56:37  j_novak
 * Added the possibility to have a shift in the quantum number l in the operator.
 *
 * Revision 1.13  2004/06/08 14:01:27  j_novak
 * *** empty log message ***
 *
 * Revision 1.11  2004/03/04 15:15:48  e_gourgoulhon
 * Treatment of case fj in state ETATZERO at the end.
 *
 * Revision 1.10  2004/03/01 13:30:28  j_novak
 * Corrected some errors
 *
 * Revision 1.9  2004/03/01 09:57:03  j_novak
 * the wave equation is solved with Scalars. It now accepts a grid with a
 * compactified external domain, which the solver ignores and where it copies
 * the values of the field from one time-step to the next.
 *
 * Revision 1.8  2003/07/22 13:24:48  j_novak
 * *** empty log message ***
 *
 * Revision 1.7  2003/06/20 10:08:12  j_novak
 * *** empty log message ***
 *
 * Revision 1.6  2003/06/20 09:27:10  j_novak
 * Modif commentaires.
 *
 * Revision 1.5  2003/06/19 16:16:38  j_novak
 * Parabolic approximation for a non flat dalembert operator
 *
 * Revision 1.4  2003/06/18 08:45:27  j_novak
 * In class Mg3d: added the member get_radial, returning only a radial grid
 * For dAlembert solver: the way the coefficients of the operator are defined has been changed.
 *
 * Revision 1.3  2002/01/03 15:30:28  j_novak
 * Some comments modified.
 *
 * Revision 1.2  2002/01/02 14:07:57  j_novak
 * Dalembert equation is now solved in the shells. However, the number of
 * points in theta and phi must be the same in each domain. The solver is not
 * completely tested (beta version!).
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 1.7  2001/10/16  10:04:22  novak
 * cleaning (no more source terms for enhanced BC)
 *
 * Revision 1.6  2001/07/19 14:07:15  novak
 * tentative for new outgoing boundary condition
 *
 * Revision 1.5  2000/12/04 15:01:34  novak
 * *** empty log message ***
 *
 * Revision 1.4  2000/12/04 14:20:36  novak
 * odd case enabled
 *
 * Revision 1.3  2000/11/27 14:54:51  novak
 * 3D boundary conditions operational
 *
 * Revision 1.2  2000/10/24 16:18:34  novak
 * Outgoing wave boundary conditions and addition of the Tbl coeff
 *
 * Revision 1.1  2000/10/19 14:17:39  novak
 * Initial revision
 *
 *
 * $Header$
 *
 */
//Header C++
#include <cmath>

// Header Lorene:
#include "tensor.h"
#include "param.h"
#include "proto.h"

//**************************************************************************

namespace Lorene {

void Map_af::dalembert(Param& par, Scalar& fjp1, const Scalar& fj, const Scalar& fjm1,
		       const Scalar& source) const {


    assert(source.get_etat() != ETATNONDEF) ; 
    assert(source.get_mp().get_mg() == mg) ; 
    assert(fj.get_etat() != ETATNONDEF) ; 
    assert(fj.get_mp().get_mg() == mg) ; 
    assert(fjm1.get_etat() != ETATNONDEF) ; 
    assert(fjm1.get_mp().get_mg() == mg) ; 
    assert(fjp1.get_mp().get_mg() == mg) ;

    assert(par.get_n_double() == 1) ;
    assert(par.get_n_int() >= 1) ;
    assert(par.get_n_int_mod() == 1) ;
    int& nap = par.get_int_mod() ;
    assert ((nap == 0) || (par.get_n_tbl_mod() > 1)) ;

    int nz = mg->get_nzone() ;
    bool ced = (mg->get_type_r(nz-1) == UNSURR) ;
    int nz0 = (ced ? nz - 1 : nz) ;
    double dt = par.get_double() ;

    Scalar fj_local = fj ;
    Scalar fjm1_local = fjm1 ;
    if (ced) {
	fj_local.annule_domain(nz-1) ;
	fjm1_local.annule_domain(nz-1) ;
    }
    Scalar sigma = 2*fj_local - fjm1_local ; // The source (first part)     

    // Coefficients
    //-------------
    
    Tbl* coeff ;
    if (nap == 0) {
      coeff = new Tbl(12,nz);
      coeff->set_etat_qcq() ;
      par.add_tbl_mod(*coeff) ;
    }
    else 
      coeff= &par.get_tbl_mod() ;
    Tbl a1(nz) ; a1 = 1 ; //Flat dalembertian
    Tbl a2(nz) ; a2 = 0 ;
    Tbl a3(nz) ; a3 = 0 ;

    if (par.get_n_tensor_mod() > 0) { // Metric in front of the dalembertian
      assert(par.get_n_tensor_mod() == 1) ;
      Scalar* metri = dynamic_cast<Scalar*>(&par.get_tensor_mod()) ;
      assert(metri != 0x0) ;
      assert (metri->get_etat() == ETATQCQ) ;

      const Map_af* tmap ; //Spherically symmetric grid and mapping
      if (nap == 0) {
	double* bornes = new double[nz+1] ;
	bornes[0] = beta[0] ;
	for (int i=0; i<nz; i++) bornes[i+1] = alpha[i] + beta[i] ;
	tmap = new Map_af(*mg->get_radial() , bornes) ;
	par.add_map(*tmap) ;
	delete [] bornes ;
      }
      else {
	tmap = dynamic_cast<const Map_af*>(&par.get_map()) ;
	assert (tmap != 0x0) ;
      }
      metri->set_spectral_va().ylm() ;

      Scalar xmetr(*tmap) ; // l=0 part of the potential in front of the Laplacian
      xmetr.set_etat_qcq() ;
      xmetr.std_spectral_base() ;
      xmetr.set_spectral_va().set_base_t(T_LEG_PP) ; // Only l=0 matters in any case...
      xmetr.set_spectral_va().set_etat_cf_qcq() ;
      Mtbl_cf* mt = xmetr.set_spectral_va().c_cf ; 
      mt->annule_hard() ;
      for (int lz=0; lz<nz0; lz++) 
	for (int ir=0; ir<mg->get_nr(lz); ir++) 
	  mt->set(lz,0,0,ir) = (*metri->get_spectral_va().c_cf)(lz, 0, 0, ir) ; //only l=0
	
      if (mg->get_nt(0) != 1) xmetr = xmetr / sqrt(double(2)) ; //!!!
      xmetr.set_spectral_va().ylm_i() ;
      xmetr.set_spectral_va().coef_i() ;
      const Mtbl& erre = this->r ;

      a1.set_etat_qcq() ;
      a2.set_etat_qcq() ;
      a3.set_etat_qcq() ;
      Scalar mime(*this) ;
      mime.annule_hard() ;
      for (int lz=0; lz<nz0; lz++) {
	int nr = mg->get_nr(lz);
	double r1 = erre(lz, 0, 0, nr-1) ;
	double rm1 = erre(lz, 0, 0, 0) ;
	double x1 = xmetr.val_grid_point(lz, 0, 0, nr-1) ;
	double xm1 = xmetr.val_grid_point(lz, 0, 0, 0) ;
	
	if (mg->get_type_r(lz) == RARE) { //In the nucleus, no a2*r
      	  a1.set(lz) = xm1 ;
    	  a2.set(lz) = 0 ;
      	  a3.set(lz) = (x1 - a1(lz)) / (r1 * r1);
	}
	else { // In the shells, general case
	  int i0 = (nr-1)/2 ;
	  double r0 = erre(lz, 0, 0, i0) ;
	  double x0 = xmetr.val_grid_point(lz, 0, 0, i0) ;
	  double p1 = (r1 - rm1)*(r1 - r0) ;
	  double pm1 = (r0 - rm1)*(r1 - rm1) ;
	  double p0 = (r0 - rm1)*(r1 - r0) ;
	  a1.set(lz) = xm1*r1*r0/pm1 + x1*rm1*r0/p1 - x0*rm1*r1/p0 ;
	  a2.set(lz) = x0*(rm1 + r1)/p0 - xm1*(r1 + r0)/pm1 
	    - x1*(rm1 + r0)/p1 ;
	  a3.set(lz) = xm1/pm1+x1/p1-x0/p0 ;
	}

	for (int k=0; k<mg->get_np(lz)+2; k++) 
	  for (int j=0; j<mg->get_nt(lz); j++) 
	    for (int i=0; i<nr; i++)
	      mime.set_grid_point(lz, k, j, i) = a1(lz) + erre(lz, 0, 0, i)*
		(a2(lz) + erre(lz, 0, 0, i)*a3(lz)) ;

	Tbl diff = metri->domain(lz) - mime.domain(lz) ;
	double offset = max(diff) ; // Not sure that this is really 
  	a1.set(lz) += offset ;  // necessary (supposed to ensure stability).
  	mime.set_domain(lz) += offset ;
      }

      Scalar reste = (*metri - mime)*fj_local.laplacian() ;
      if (ced) reste.annule_domain(nz-1) ;
      sigma += (dt*dt)*(source + reste) ;
      if (ced) sigma.annule_domain(nz-1) ;
      sigma +=  (0.5*dt*dt)*mime*fjm1_local.laplacian() ; //Source (2nd part)
    }
    else {
      sigma += (dt*dt) * source ;
      if (ced) sigma.annule_domain(nz-1) ;
      sigma += (0.5*dt*dt)*fjm1_local.laplacian() ;
      if (par.get_n_int() > 1) { //there is a shift in the quantum number l
	  int dl = -1 ;
	  int l_min = par.get_int(1) ;
	  sigma.set_spectral_va().ylm() ;      
	  Scalar tmp = fjm1_local ;
	  tmp.div_r() ; tmp.div_r() ; // f^(J-1) / r^2
	  tmp.set_spectral_va().ylm() ;
	  const Base_val& base = tmp.get_spectral_base() ;
	  int l_q, m_q, baser ;
	
	  for (int lz=0; lz<nz-1; lz++) {
	      int nt = mg->get_nt(lz) ;
	      int np = mg->get_np(lz) ;
	      for (int k=0; k<np+2; k++) 
		  for (int j=0; j<nt; j++) {
		      base.give_quant_numbers(lz, k, j, m_q, l_q, baser) ;
		      if ((nullite_plm(j, nt, k, np, base) == 1) && (l_q+dl >= l_min) ) {
			  for (int i=0; i<mg->get_nr(lz); i++) {
			      sigma.set_spectral_va().c_cf->set(lz, k, j, i) -=
			      0.5*dt*dt*dl*(2*l_q + dl +1)
				  *(*tmp.get_spectral_va().c_cf)(lz, k, j, i) ;
			  }
		      }
		  }
	  }
	  if (sigma.get_spectral_va().c != 0x0) {
	      delete sigma.set_spectral_va().c ;
	      sigma.set_spectral_va().c = 0x0 ;
	  }
      }
    }
    if (ced) sigma.annule_domain(nz-1) ;

    //--------------------------------------------
    // The operator reads
    // Id - 0.5dt^2*(a1 + a2 r + a3 r^2)Laplacian
    //--------------------------------------------
    for (int i=0; i<nz; i++) {
	coeff->set(1,i) = a1(i) ; 
	coeff->set(2,i) = a2(i) ;
	coeff->set(3,i) = a3(i) ;
	coeff->set(4,i) = 0. ;
	coeff->set(5,i) = 0. ;
	coeff->set(6,i) = 0. ;
	coeff->set(7,i) = 0. ;
	coeff->set(8,i) = 0. ;
	coeff->set(9,i) = 0. ; 
	coeff->set(10,i) = beta[i] ;
	coeff->set(11,i) = alpha[i] ;
    }

    // Defining the boundary conditions
    // --------------------------------
    double R = this->val_r(nz0-1, 1., 0., 0.) ;
    int nr = mg->get_nr(nz0-1) ;
    int nt = mg->get_nt(nz0-1) ;
    int np2 = mg->get_np(nz0-1) + 2;

    // For each pair of quantic numbers l, m one the result must satisfy
    // bc1 * f_{l,m} (R) + bc2 * f_{l,m}'(R) = tbc3_{l,m}
    // Memory is allocated for the parameter (par) at first call
    double* bc1 ;
    double* bc2 ;
    Tbl* tbc3 ;
    Tbl* phijm1 = 0x0 ;
    Tbl* phij = 0x0 ;
    if (nap == 0) { 
      bc1 = new double ;
      bc2 = new double ;
      tbc3 = new Tbl(np2,nt) ;
      par.add_double_mod(*bc1,1) ;
      par.add_double_mod(*bc2,2) ;
      par.add_tbl_mod(*tbc3,1) ;
      // Hereafter the enhanced outgoing-wave condition needs 2 auxiliary
      // functions phij and phijm1 to define the evolution on the boundary
      // surface (outer sphere).
      if (par.get_int(0) == 2) {
	phijm1 = new Tbl(np2,nt) ;
	phij = new Tbl(np2,nt) ;
	par.add_tbl_mod(*phijm1,2) ;
	par.add_tbl_mod(*phij,3) ;
	phij->annule_hard() ;
	phijm1->annule_hard() ;
      }
      nap = 1 ;
    }
    else {
      bc1 = &par.get_double_mod(1) ;
      bc2 = &par.get_double_mod(2) ;
      tbc3 = &par.get_tbl_mod(1) ;
      if (par.get_int(0) == 2) {
	phijm1 = &par.get_tbl_mod(2) ;
	phij = &par.get_tbl_mod(3) ;
      }
    }
    switch (par.get_int(0)) {
    case 0:   // Homogeneous boundary conditions (f(t,r=R) =0)
      *bc1 = 1 ;
      *bc2 = 0 ;
      
      *tbc3 = 0 ;
      
      break ;
    case 1:  { // Outgoing wave condition (f(t,r) = 1/r S(t-r/c))
      Valeur bound3(mg) ;
      bound3 = R*(4*fj_local.get_spectral_va() - fjm1_local.get_spectral_va()) ;
      if (bound3.get_etat() == ETATZERO) {
	*bc1 = 3*R + 2*dt ;
	*bc2 = 2*R*dt ;
	*tbc3 = 0 ;
      }
      else {
	if (nz0>1) bound3.annule(0,nz0-2) ;

	bound3.coef() ;
	bound3.ylm() ;
	
	*bc1 = 3*R + 2*dt ;
	*bc2 = 2*R*dt ;
	
	tbc3->set_etat_qcq() ;
	double val ;
	for (int k=0; k<np2; k++)
	  for (int j=0; j<nt; j++) {
	    val = 0. ;
	    for (int i=0; i<nr; i++) 
	      val += (*bound3.c_cf)(nz0-1,k,j,i) ;
	    tbc3->set(k,j) = val ;
	  }
      }
      break ;
    }
    /******************************************************************
     * Enhanced outgoing wave condition. 
     * Time integration of the wave equation on the sphere for the 
     * auxiliary function phij.
     *****************************************************************/
     case 2: { 
      Valeur souphi(mg) ;
      souphi = fj_local.get_spectral_va()/R - fj_local.dsdr().get_spectral_va() ;
      if (nz0>1) souphi.annule(0,nz0-2) ;
      souphi.coef() ;
      souphi.ylm() ;

      bool zero = (souphi.get_etat() == ETATZERO) ;
      if (zero) {
	Base_val base_ref(mg->std_base_scal()) ; //## Maybe not good...
	base_ref.dsdx() ;
	base_ref.ylm() ;
	souphi.set_base(base_ref) ;
      }

      int l_s, m_s, base_r ;
      double val ;
      int dl = (par.get_n_int() > 1) ? -1 : 0 ;
      for (int k=0; k<np2; k++) {
	for (int j=0; j<nt; j++) {
	  donne_lm(nz, nz0-1, j, k, souphi.base, m_s, l_s, base_r) ;
	  l_s += dl ;
	  val = 0. ;
	  if (!zero)
	      val = -4*dt*dt*l_s*(l_s+1)*souphi.c_cf->val_out_bound_jk(nz0-1, j, k) ;
	  double multi = 8*R*R + dt*dt*(6+3*l_s*(l_s+1)) + 12*R*dt ;
	  val = ( 16*R*R*(*phij)(k,j) -
		  (multi-24*R*dt)*(*phijm1)(k,j) 
		  + val)/multi ;
	  phijm1->set(k,j) = (*phij)(k,j) ; 
	  phij->set(k,j) = val ;
  	}
      }
      Valeur bound3(mg) ;
      *bc1 = 3*R + 2*dt ;
      *bc2 = 2*R*dt ;
      bound3 = R*(4*fj_local.get_spectral_va() - fjm1_local.get_spectral_va()) ;
      if (bound3.get_etat() == ETATZERO) *tbc3 = 0 ;
      else {
	if (nz0 > 1) bound3.annule(0,nz0-2) ;
	bound3.coef() ;
	bound3.ylm() ;
	tbc3->set_etat_qcq() ;
	for (int k=0; k<np2; k++)
	  for (int j=0; j<nt; j++) {
	    val = 0. ;
	    for (int i=0; i<nr; i++) 
	      val += (*bound3.c_cf)(nz0-1,k,j,i) ;
	    tbc3->set(k,j) = val + 2*R*dt*(*phij)(k,j);
	  }
      }
      break ;
    }
    default:
      cout << "ERROR: Map_af::dalembert" << endl ;
      cout << "The boundary condition par.get_int(0) = "<< par.get_int(0) 
	   << " is unknown!" << endl ;
      abort() ;
    }

    if (sigma.get_etat() == ETATZERO) {
      fjp1.set_etat_zero() ;
      return ;  
    }

    // Spherical harmonic expansion of the source
    // ------------------------------------------    
    Valeur& sourva = sigma.set_spectral_va() ; 

    // Spectral coefficients of the source
    assert(sourva.get_etat() == ETATQCQ) ; 
    sourva.ylm() ;			// spherical harmonic transforms 

    // Final result returned as a Scalar
    // ------------------------------
    fjp1.set_etat_zero() ;  // to call Scalar::del_t().
    fjp1.set_etat_qcq() ; 
    
    // Call to the Mtbl_cf version
    // ---------------------------
    fjp1.set_spectral_va() = sol_dalembert(par, *this, *(sourva.c_cf) ) ;
    fjp1.set_spectral_va().ylm_i() ; // Back to standard basis.	 

    if (ced) {
        if (fj.get_etat() == ETATZERO) {
            fjp1.annule_domain(nz-1) ; 
        }
        else {
            fjp1.set_domain(nz-1) = fj.domain(nz-1) ;
        }
	fjp1.set_dzpuis(fj.get_dzpuis()) ;
    }
}

}
