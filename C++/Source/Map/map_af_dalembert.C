/*
 *   Copyright (c) 2000-2001 Jerome Novak
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
 * Revision 1.1  2001/11/20 15:19:27  e_gourgoulhon
 * Initial revision
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
#include <math.h>

// Header Lorene:
#include "cmp.h"
#include "param.h"
#include "proto.h"

//**************************************************************************



void Map_af::dalembert(Param& par, Cmp& fjp1, const Cmp& fj, const Cmp& fjm1,
		       const Cmp& source) const {
    
    assert(source.get_etat() != ETATNONDEF) ; 
    assert(source.get_mp()->get_mg() == mg) ; 
    assert(fj.get_etat() != ETATNONDEF) ; 
    assert(fj.get_mp()->get_mg() == mg) ; 
    assert(fjm1.get_etat() != ETATNONDEF) ; 
    assert(fjm1.get_mp()->get_mg() == mg) ; 
    assert(fjp1.get_mp()->get_mg() == mg) ;

    assert(par.get_n_double() == 1) ;
    assert(par.get_n_int() == 1) ;
    static int nap = 0;
    assert ((nap == 0) || (par.get_n_tbl_mod() > 1)) ;

    // The source and explicit parts
    // -----------------------------

    Cmp sigma(2*fj) ;
    int nz = mg->get_nzone() ;
    double dt = par.get_double() ;
    sigma += dt*dt * (source + 0.5*fjm1.laplacien(0)) - fjm1 ;

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
    if (par.get_n_cmp_mod() > 0) { // Metric in front of the dalembertian
      cout << "Map_af_dalembert: ERROR:" <<endl ;
      cout << "only a flat dalembertian is implemented!" << endl ;
      exit(-1) ;
    }
    else  // Flat dalembertian
      for (int i=0; i<nz; i++) {
	coeff->set(1,i) = 1. ;
	coeff->set(2,i) = 0. ;
	coeff->set(3,i) = 0. ;
	coeff->set(4,i) = 0. ;
	coeff->set(5,i) = 0. ;
	coeff->set(6,i) = 2. ;
	coeff->set(7,i) = 0. ;
	coeff->set(8,i) = 0. ;
	coeff->set(9,i) = 0. ; // The -l(l+1) term will come later...
	coeff->set(10,i) = beta[i] ;
	coeff->set(11,i) = alpha[i] ;
      }

    // Defining the boundary conditions
    // --------------------------------

    double R = this->val_r(nz-1, 1., 0., 0.) ;
    int nr = mg->get_nr(nz-1) ;
    int nt = mg->get_nt(nz-1) ;
    int np2 = mg->get_np(nz-1) + 2;


    // For each pair of quantic numbers l, m one the result must satisfy
    // bc1 * f_{l,m} (R) + bc2 * f_{l,m}'(R) = tbc3_{l,m}
    // Memory is allocated for the parameter (par) at first call

    double* bc1 ;
    double* bc2 ;
    Tbl* tbc3 ;
    Tbl* phijm1 ;
    Tbl* phij ;
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
      bound3 = R*(4*fj.va - fjm1.va) ;
      if (bound3.get_etat() == ETATZERO) {
	*bc1 = 3*R + 2*dt ;
	*bc2 = 2*R*dt ;
	*tbc3 = 0 ;
      }
      else {
	if (nz>1) bound3.annule(0,nz-2) ;
	
	bound3.coef() ;
	bound3.ylm() ;
	
	*bc1 = 3*R + 2*dt ;
	*bc2 = 2*R*dt ;
	
	tbc3->set_etat_qcq() ;
	double val ;
	double* tmp = new double[nr] ;
	int b1,b2, base_r ;
	for (int k=0; k<np2; k++)
	  for (int j=0; j<nt; j++) {
	    donne_lm(nz, nz-1, j, k, bound3.base, b1, b2, base_r) ;
	    for (int i=0; i<nr; i++) 
	      tmp[i] = (*bound3.c_cf)(nz-1,k,j,i) ;
	    if (base_r == R_CHEBP) som_r_chebp(tmp, nr, 1, 1, 1., &val) ;
	    else som_r_chebi(tmp, nr, 1, 1, 1., &val) ;
	    tbc3->set(k,j) = val ;
	  }
	delete[] tmp ;
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
      souphi = fj.va/R - fj.dsdr().va ;
      souphi.coef() ;
      souphi.ylm() ;
      souphi = 4*dt*dt*souphi.lapang() ; 

      bool zero = (souphi.get_etat() == ETATZERO) ;
      if (zero) {
	Base_val base_ref(mg->std_base_scal()) ;
	base_ref.dsdx() ;
	base_ref.ylm() ;
	souphi.set_base(base_ref) ;
      }

      int l_s, m_s, base_r ;
      double val ;
      double* tmp = new double[nr] ;
      for (int k=0; k<np2; k++) {
	for (int j=0; j<nt; j++) {
	  donne_lm(nz, nz-1, j, k, souphi.base, m_s, l_s, base_r) ;
	  if (zero) {
	    val = 0 ;
	  }
	  else {
	    for (int i=0; i<nr; i++) 
	      tmp[i] = (*souphi.c_cf)(nz-1,k,j,i) ;
	    if (base_r == R_CHEBP) som_r_chebp(tmp, nr, 1, 1, 1., &val) ;
	    else som_r_chebi(tmp, nr, 1, 1, 1., &val) ;
	  }
	  double multi = 8*R*R + dt*dt*(6+3*l_s*(l_s+1)) + 12*R*dt ;
	  val = ( 16*R*R*(*phij)(k,j) -
		  (multi-24*R*dt)*(*phijm1)(k,j) 
		  + val)/multi ;
	  phijm1->set(k,j) = (*phij)(k,j) ; 
	  phij->set(k,j) = val ;
  	}
      }
      delete[] tmp ;
      Valeur bound3(mg) ;
      *bc1 = 3*R + 2*dt ;
      *bc2 = 2*R*dt ;
      bound3 = R*(4*fj.va - fjm1.va) ;
      if (bound3.get_etat() == ETATZERO) *tbc3 = 0 ;
      else {
	if (nz>1) bound3.annule(0,nz-2) ;

	bound3.coef() ;
	bound3.ylm() ;
	tbc3->set_etat_qcq() ;
	double* tmp = new double[nr] ;
	int b1,b2 ;
	for (int k=0; k<np2; k++)
	  for (int j=0; j<nt; j++) {
	    donne_lm(nz, nz-1, j, k, bound3.base, b1, b2, base_r) ;
	    for (int i=0; i<nr; i++) 
	      tmp[i] = (*bound3.c_cf)(nz-1,k,j,i) ;
	    if (base_r == R_CHEBP) som_r_chebp(tmp, nr, 1, 1, 1., &val) ;
	    else som_r_chebi(tmp, nr, 1, 1, 1., &val) ;
	    tbc3->set(k,j) = val + 2*R*dt*(*phij)(k,j);
	  }
	delete[] tmp ;
      }
      break ;
    }
    default:
      cout << "ERROR: Map_af::dalembert" << endl ;
      cout << "The boundary condition par.get_int(0) = "<< par.get_int(0) 
	   << " is unknown!" << endl ;
      abort() ;
    }

    // Spherical harmonic expansion of the source
    // ------------------------------------------
    
    Valeur& sourva = sigma.va ; 

    if (sourva.get_etat() == ETATZERO) {
	fjp1.set_etat_zero() ;
	return ;  
    }

    // Spectral coefficients of the source
    assert(sourva.get_etat() == ETATQCQ) ; 
    
    sourva.coef() ; 
    sourva.ylm() ;			// spherical harmonic transforms 


    // Call to the Mtbl_cf version
    // ---------------------------
    Mtbl_cf resu = sol_dalembert(par, *this, *(sourva.c_cf) ) ;
    
    // Final result returned as a Cmp
    // ------------------------------
    
    fjp1.set_etat_zero() ;  // to call Cmp::del_t().

    fjp1.set_etat_qcq() ; 
    
    fjp1.va = resu ;
    (fjp1.va).ylm_i() ; // Back to standard basis.	 

}


