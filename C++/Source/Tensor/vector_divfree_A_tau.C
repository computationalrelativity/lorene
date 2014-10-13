/*
 *   Copyright (c) 2005 Philippe Grandclement
 
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


char vector_divfree_A_tau[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.4  2014/10/13 08:53:45  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.3  2014/10/06 15:13:20  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2013/06/05 15:10:43  j_novak
 * Suppression of FINJAC sampling in r. This Jacobi(0,2) base is now
 * available by setting colloc_r to BASE_JAC02 in the Mg3d constructor.
 *
 * Revision 1.1  2008/08/27 09:01:27  jl_cornou
 * Methods for solving Dirac systems for divergence free vectors
 *
 * Revision 1.6  2007/12/14 10:19:34  jl_cornou
 * *** empty log message ***
 *
 * Revision 1.4  2005/11/24 14:07:54  j_novak
 * Use of Matrice::annule_hard()
 *
 * Revision 1.3  2005/08/26 14:02:41  p_grandclement
 * Modification of the elliptic solver that matches with an oscillatory exterior solution
 * small correction in Poisson tau also...
 *
 * Revision 1.2  2005/08/25 12:16:01  p_grandclement
 * *** empty log message ***
 *
 * 
 * $Header$
 *
 */

// C headers
#include <cassert>
#include <cmath>

// Lorene headers
#include "tensor.h"
#include "diff.h"
#include "proto.h"
#include "param.h"
#include "matrice.h"
#include "type_parite.h"


namespace Lorene {
void Vector_divfree::sol_Dirac_A_tau(const Scalar& aaa, Scalar& eta_tilde, Scalar& vr, const Param* par_bc) const {

	const Map_af* mp_aff = dynamic_cast<const Map_af*>(mp) ;
	assert(mp_aff != 0x0) ; // Only affine mapping for the moment

	const Mg3d& mgrid = *mp_aff->get_mg();
	assert((mgrid.get_type_r(0) == RARE) || (mgrid.get_type_r(0) == FIN));
	if (aaa.get_etat() == ETATZERO) {
		eta_tilde = 0;
		vr = 0 ;
		return ;
	}
	assert(aaa.get_etat() != ETATNONDEF) ;
	int nz = mgrid.get_nzone();
	int nzm1 = nz - 1 ;
	bool ced = (mgrid.get_type_r(nzm1) == UNSURR) ;
	int n_shell = ced ? nz-2 : nzm1 ;
	int nz_bc = nzm1 ;
	if (par_bc != 0x0)
		if (par_bc->get_n_int() > 0) nz_bc = par_bc->get_int();
	n_shell = (nz_bc < n_shell ? nz_bc : n_shell) ;
	//bool cedbc (ced && (nz_bc == nzm1));
	
//#ifndef NDEBUG
//    	if (!cedbc) {
//		assert(par_bc != 0x0) ;
//		assert(par_bc->get_n_tbl_mod() >= 3) ;
//    	}
//#endif
	int nt = mgrid.get_nt(0) ;
    	int np = mgrid.get_np(0) ;
	int nr ;	

	Scalar source = aaa ;
	Scalar source_coq = aaa ;
	source_coq.annule_domain(0);
	
	int dzpuis = source.get_dzpuis();

	if (ced) source_coq.annule_domain(nzm1) ;
    	source_coq.mult_r() ;
    	source.set_spectral_va().ylm() ;
    	source_coq.set_spectral_va().ylm() ;
    	Base_val base = source.get_spectral_base() ;
    	base.mult_x() ;

	eta_tilde.annule_hard() ;
    	eta_tilde.set_spectral_base(base) ;
    	eta_tilde.set_spectral_va().set_etat_cf_qcq() ;
    	eta_tilde.set_spectral_va().c_cf->annule_hard() ;
    	vr.annule_hard() ;
    	vr.set_spectral_base(base) ;
    	vr.set_spectral_va().set_etat_cf_qcq() ;
    	vr.set_spectral_va().c_cf->annule_hard() ;   
	Mtbl_cf& meta = *eta_tilde.set_spectral_va().c_cf ;
	Mtbl_cf& mvr = *vr.set_spectral_va().c_cf ; 


	int base_r ;


    	// Determination of the size of the systeme :
    	int size = 0 ;
    	int max_nr = 0 ;
    	for (int l=0 ; l<nz ; l++) { 
    		nr = mgrid.get_nr(l) ;
        	size += 2*nr ;
		if (nr > max_nr)
	    	max_nr = nr ;
    	}

	Matrice systeme (size, size) ;
    	systeme.set_etat_qcq() ;
    	Tbl sec_membre (size) ;
   
    	np = mgrid.get_np(0) ;
    	nt = mgrid.get_nt(0) ;
    	
	
	double alpha, beta ; 
	int l_quant, m_quant ;    

    	// On bosse pour chaque l, m :
    	for (int k=0 ; k<np+1 ; k++)
      		for (int j=0 ; j<nt ; j++)
      			if (nullite_plm(j, nt, k, np, base) == 1) {

	systeme.annule_hard() ;
	sec_membre.annule_hard() ;
	     
	int column_courant = 0 ;
	int ligne_courant = 0 ;
		
        	//--------------------------
		//       NUCLEUS
		//--------------------------
		
	nr = mgrid.get_nr(0) ; 
	alpha = (*mp_aff).get_alpha()[0] ;
        base.give_quant_numbers (0, k, j, m_quant, l_quant, base_r) ;
	Diff_dsdx odn(base_r, nr) ; const Matrice& mdn = odn.get_matrice() ;
	Diff_sx osn(base_r, nr) ; const Matrice& msn = osn.get_matrice() ;  
        
	

	int nbr_cl = 0 ;
	// RARE case	
	if (source.get_mp().get_mg()->get_type_r(0) == RARE) {
	// regularity conditions for eta :
	if (l_quant > 1) {
	     nbr_cl = 1 ;
	     if (l_quant%2==0) {
	        //Even case
		for (int col=0 ; col<nr ; col++) {
		    if (col%2==0) {
		        systeme.set(ligne_courant, col+column_courant) = 1 ;
			systeme.set(ligne_courant, col+column_courant+nr) = 0 ; }
		    else  {
		        systeme.set(ligne_courant, col+column_courant) = -1 ;
			systeme.set(ligne_courant, col+column_courant+nr) = 0 ; } 
		}
	    }
	    else {
	     //Odd case
	         for (int col=0 ; col<nr ; col++) {
		    if (col%2==0) {
		        systeme.set(ligne_courant, col+column_courant) = 2*col+1 ;
			systeme.set(ligne_courant, col+column_courant+nr) = 0 ; }
		    else {
		        systeme.set(ligne_courant, col+column_courant) = -(2*col+1) ;
			systeme.set(ligne_courant, col+column_courant+nr) = 0 ; }
		 }
	    }
	}
	}

	// R_JACO02 case
	else {
	  assert (base_r == R_JACO02) ;
	// regularity conditions for eta :
	if (l_quant == 0) {
	    nbr_cl = 1 ;
		for (int col=0 ; col<nr ; col++) {
		    systeme.set(ligne_courant, col+column_courant) =  col*(col+1)*(col+2)*(col+3)/double(12)*(2*(col%2)-1);
		    systeme.set(ligne_courant, col+column_courant+nr) = 0 ;
		}
	    }
	else if (l_quant == 1) {
	    nbr_cl = 1 ;
		for (int col=0 ; col<nr ; col++) {
		    systeme.set(ligne_courant, col+column_courant) = (col+1)*(col+2)/double(2)*(1-2*(col%2)) ;
		    systeme.set(ligne_courant, col+column_courant+nr) = 0 ;
		}
	    }
	else {
	    nbr_cl = 2 ;
		for (int col=0 ; col<nr ; col++) {
		    systeme.set(ligne_courant, col+column_courant) = (col+1)*(col+2)/double(2)*(1-2*(col%2)) ;
		    systeme.set(ligne_courant, col+column_courant+nr) = 0 ;
		    systeme.set(ligne_courant+1, col+column_courant) = col*(col+1)*(col+2)*(col+3)/double(12)*(2*(col%2)-1) ;
		    systeme.set(ligne_courant+1, col+column_courant+nr) = 0 ;
		}
	    }
	}	
	ligne_courant += nbr_cl ;

	// Premiere partie de l'operateur :
	for (int lig=0 ; lig<nr-1-nbr_cl ; lig++) {
	    for (int col=0 ; col<nr ; col++) {
	        systeme.set(lig+ligne_courant,col+column_courant) = mdn(lig,col) + msn(lig,col) ;
	        sec_membre.set(lig+ligne_courant) = alpha*alpha*(*source.get_spectral_va().c_cf)(0, k, j, lig) ;
		systeme.set(lig+ligne_courant,col+column_courant+nr)= - msn(lig,col) ;
	//	systeme.set(lig+ligne_courant+nr) = 0 ;
	    }
	}


	ligne_courant += nr-1-nbr_cl ;
	

	// RARE case	
	if (source.get_mp().get_mg()->get_type_r(0) == RARE) {
	// regularity conditions for vr :
	if (l_quant > 1) {
	     nbr_cl = 1 ;
	     if (l_quant%2==0) {
	        //Even case
		for (int col=0 ; col<nr ; col++) {
		    if (col%2==0) {
		        systeme.set(ligne_courant, col+column_courant) = 0 ;
			systeme.set(ligne_courant, col+column_courant+nr) = 1 ; }
		    else {
		        systeme.set(ligne_courant, col+column_courant) = 0 ;
			systeme.set(ligne_courant, col+column_courant+nr) = -1 ; }
		}
	    }
	
	     else {
	     //Odd case
	         for (int col=0 ; col<nr ; col++) {
		    if (col%2==0) {
		        systeme.set(ligne_courant, col+column_courant) = 0 ;
			systeme.set(ligne_courant, col+column_courant+nr) = 2*col+1 ; }
		    else {
		        systeme.set(ligne_courant, col+column_courant) = 0 ;
			systeme.set(ligne_courant, col+column_courant+nr) = -(2*col+1) ; }
		}
	    }
	  }
	}

	// R_JACO02 case
	else {
	  assert (base_r == R_JACO02) ;
	// regularity conditions for vr :
	if (l_quant == 0) {
	    nbr_cl = 1 ;
		for (int col=0 ; col<nr ; col++) {
		    systeme.set(ligne_courant, col+column_courant) = 0 ;
		    systeme.set(ligne_courant, col+column_courant+nr) = col*(col+1)*(col+2)*(col+3)/double(12)*(2*(col%2)-1) ;
		}
	    }
	else if (l_quant == 1) {
	    nbr_cl = 1 ;
		for (int col=0 ; col<nr ; col++) {
		    systeme.set(ligne_courant, col+column_courant) = 0 ;
		    systeme.set(ligne_courant, col+column_courant+nr) = (col+1)*(col+2)/double(2)*(1-2*(col%2)) ;
		}
	    }
	else {
	    nbr_cl = 2 ;
		for (int col=0 ; col<nr ; col++) {
		    systeme.set(ligne_courant, col+column_courant) = 0 ;
		    systeme.set(ligne_courant, col+column_courant+nr) = (col+1)*(col+2)/double(2)*(1-2*(col%2)) ;
		    systeme.set(ligne_courant+1, col+column_courant) = 0 ;
		    systeme.set(ligne_courant+1, col+column_courant+nr) = col*(col+1)*(col+2)*(col+3)/double(12)*(2*(col%2)-1) ;
		}
	    }
	}	
	ligne_courant += nbr_cl ;

	// Deuxieme partie de l'operateur
	for (int lig=0 ; lig<nr-1-nbr_cl ; lig++) {
	    for (int col=0 ; col<nr ; col++) {
	        systeme.set(lig+ligne_courant,col+column_courant) = - l_quant*(l_quant+1)*msn(lig,col) ;
	        sec_membre.set(lig+ligne_courant) = 0 ;
		systeme.set(lig+ligne_courant,col+column_courant+nr)= mdn(lig,col) + 2*msn(lig,col) ;
	//	systeme.set(lig+ligne_courant+nr) = 0 ;
	    }
	}
	
	
	ligne_courant += nr-1-nbr_cl ;
	
	  
	// Le raccord pour eta
	for (int col=0 ; col<nr ; col++) {
	     if (source.get_mp().get_mg()->get_type_r(0) == RARE) {
	     // La fonction eta
	     systeme.set(ligne_courant, col+column_courant) = 1 ;
	     systeme.set(ligne_courant, col+column_courant+nr) = 0 ;
	     // Sa derivee :
	     if (l_quant%2==0) {
	         systeme.set(ligne_courant+1, col+column_courant) = 4*col*col/alpha ;
		 systeme.set(ligne_courant+1, col+column_courant+nr) = 0 ;
	     }
	     else { 
	         systeme.set(ligne_courant+1, col+column_courant) = (2*col+1)*(2*col+1)/alpha ;
		 systeme.set(ligne_courant+1, col+column_courant+nr) = 0 ;
 		}
	      }
	     else { // Cas R_JACO02
	       assert( base_r == R_JACO02 ) ;
	     // La fonction eta
	     systeme.set(ligne_courant, col+column_courant) = 1 ;
	     systeme.set(ligne_courant, col+column_courant+nr) = 0 ; 
	     // Sa derivee :
	     systeme.set(ligne_courant+1, col+column_courant) = col*(col+3)/double(2)/alpha ;
	     systeme.set(ligne_courant+1, col+column_courant+nr) = 0 ;
	      }
	     }

	ligne_courant += 2 ;
	// Le raccord pour vr
	for (int col=0 ; col<nr ; col++) {
	     if (source.get_mp().get_mg()->get_type_r(0) == RARE) {
	     // La fonction vr
	     systeme.set(ligne_courant, col+column_courant) = 0 ;
	     systeme.set(ligne_courant, col+column_courant+nr) = 1 ;
	     // Sa derivee :
	     if (l_quant%2==0) {
	         systeme.set(ligne_courant+1, col+column_courant) = 0 ;
		 systeme.set(ligne_courant+1, col+column_courant+nr) = 4*col*col/alpha ;
	     }
	     else { 
	         systeme.set(ligne_courant+1, col+column_courant) = 0 ;
		 systeme.set(ligne_courant+1, col+column_courant+nr) = (2*col+1)*(2*col+1)/alpha ;
 		}
	      }
	     else { // Cas R_JACO02
	       assert( base_r == R_JACO02 ) ;
	     // La fonction vr
	     systeme.set(ligne_courant, col+column_courant) = 0 ;
	     systeme.set(ligne_courant, col+column_courant+nr) = 1 ; 
	     // Sa derivee :
	     systeme.set(ligne_courant+1, col+column_courant) = 0 ;
	     systeme.set(ligne_courant+1, col+column_courant+nr) = col*(col+3)/double(2)/alpha ;
	      }
	     }
	
	ligne_courant -= 2 ; // On retourne pour le raccord dans la prochaine zone

	column_courant += 2*nr ; // On va changer de zone

		//--------------------------
		//       SHELLS
		//--------------------------
	for (int l=1 ; l<nz-1 ; l++) {
	
		nr = mgrid.get_nr(l) ; 
		alpha = (*mp_aff).get_alpha()[l] ;
		beta = (*mp_aff).get_beta()[l] ;
		
         	base.give_quant_numbers (l, k, j, m_quant, l_quant, base_r) ;  
		Diff_id ois(base_r, nr) ; const Matrice& mis = ois.get_matrice() ;
		Diff_xdsdx oxds(base_r, nr) ; const Matrice& mxds = oxds.get_matrice() ;
	
	   // matching with previous domain :
	   for (int col=0 ; col<nr ; col++) {
	     // La fonction eta
	     if (col%2==0) {
	          systeme.set(ligne_courant, col+column_courant) = -1 ;	
		  systeme.set(ligne_courant, col+column_courant+nr) = 0 ; }
	     else  {
	          systeme.set(ligne_courant, col+column_courant) = 1 ;
		  systeme.set(ligne_courant, col+column_courant+nr) = 0 ; }
	     // Sa derivee :
	     if (col%2==0) {
	         systeme.set(ligne_courant+1, col+column_courant) = col*col/alpha ;
		 systeme.set(ligne_courant+1, col+column_courant+nr) = 0 ; }
	     else  {
	         systeme.set(ligne_courant+1, col+column_courant) = -col*col/alpha ;
		 systeme.set(ligne_courant+1, col+column_courant+nr) = 0 ; }
	  }
	  ligne_courant += 2 ;

	  // matching with previous domain :
	   for (int col=0 ; col<nr ; col++) {
	     // La fonction vr
	     if (col%2==0) {
	          systeme.set(ligne_courant, col+column_courant) = 0 ;	
		  systeme.set(ligne_courant, col+column_courant+nr) = -1 ; }
	     else  {
	          systeme.set(ligne_courant, col+column_courant) = 0 ;
		  systeme.set(ligne_courant, col+column_courant+nr) = 1 ; }
	     // Sa derivee :
	     if (col%2==0) {
	         systeme.set(ligne_courant+1, col+column_courant) = 0 ;
		 systeme.set(ligne_courant+1, col+column_courant+nr) = col*col/alpha ; }
	     else  {
	         systeme.set(ligne_courant+1, col+column_courant) = 0 ;
		 systeme.set(ligne_courant+1, col+column_courant+nr) = -col*col/alpha ; }
	  }
	  ligne_courant += 2 ;
	  

	  // L'operateur (premiere et deuxieme parties) :
	  
	  // source must be multiplied by (x+echelle)^2
	  Tbl source_aux(nr) ;
	  source_aux.set_etat_qcq() ;
	  for (int i=0 ; i<nr ; i++)
	      source_aux.set(i) = (*source.get_spectral_va().c_cf)(l,k,j,i)*alpha*alpha ;
    	Tbl xso(source_aux) ;
    	Tbl xxso(source_aux) ;
    	multx_1d(nr, &xso.t, R_CHEB) ;
    	multx_1d(nr, &xxso.t, R_CHEB) ;
    	multx_1d(nr, &xxso.t, R_CHEB) ;
    	source_aux = beta*beta/alpha/alpha*source_aux+2*beta/alpha*xso+xxso ;
	  
	for (int lig=0 ; lig<nr-2 ; lig++) {
	     for (int col=0 ; col<nr ; col++) {
	         systeme.set(lig+ligne_courant,col+column_courant) = mxds(lig,col) + mis(lig,col) ;
		 systeme.set(lig+ligne_courant,col+column_courant+nr) = -mis(lig,col) ; 
		 sec_membre.set(lig+ligne_courant) = source_aux(lig) ;
		 systeme.set(lig+ligne_courant+nr-2, col+column_courant) = -l_quant*(l_quant+1) ;
		 systeme.set(lig+ligne_courant+nr-2, col+column_courant+nr) = mxds(lig,col) + 2*mis(lig,col) ;
		 sec_membre.set(lig+ligne_courant+nr-2) = 0 ; }
	  }
	  
	  
	  ligne_courant += 2*nr-4 ;
	  // Matching with the next domain :
	  for (int col=0 ; col<nr ; col++) {
	     // La fonction eta
	     systeme.set(ligne_courant, col+column_courant) = 1 ;
	     systeme.set(ligne_courant, col+column_courant+nr) = 0 ;
	     // Sa derivee :
	     systeme.set(ligne_courant+1, col+column_courant) = col*col/alpha ;
	     systeme.set(ligne_courant+1, col+column_courant+nr) = 0 ;
	     // La fonction vr
	     systeme.set(ligne_courant+2, col+column_courant) = 0 ;
	     systeme.set(ligne_courant+2, col+column_courant+nr) = 1 ;
	     // Sa derivee
	     systeme.set(ligne_courant+3, col+column_courant) = 0 ;	
	     systeme.set(ligne_courant+3, col+column_courant+nr) = col*col/alpha ; 
	     }
	     
	  column_courant += 2*nr ;   
	  }
	  

		//--------------------------
		//       ZEC
		//--------------------------
	nr = mgrid.get_nr(nz-1) ; 
	alpha = (*mp_aff).get_alpha()[nz-1] ;
	beta = (*mp_aff).get_beta()[nz-1] ;
		
	base.give_quant_numbers (nz-1, k, j, m_quant, l_quant, base_r) ;
	//work = new Matrice(laplacien_mat(nr, l_quant, 0., dzpuis, base_r)) ;
	Diff_dsdx odzec(base_r, nr) ; const Matrice& mdzec = odzec.get_matrice() ;
	Diff_sx oszec(base_r, nr) ; const Matrice& mszec = oszec.get_matrice() ;

	
	// Matching with the previous domain :
	 for (int col=0 ; col<nr ; col++) {
	     // La fonction eta
	     if (col%2==0) {
	          systeme.set(ligne_courant, col+column_courant) = -1 ;
		  systeme.set(ligne_courant, col+column_courant+nr) = 0 ; }
	     else {
	          systeme.set(ligne_courant, col+column_courant) = 1 ;
		  systeme.set(ligne_courant, col+column_courant+nr) = 0 ; }
	     // Sa derivee :
	     if (col%2==0) {
	         systeme.set(ligne_courant+1, col+column_courant) = -4*alpha*col*col ;
		 systeme.set(ligne_courant+1, col+column_courant+nr) = 0 ; }
	     else {
	         systeme.set(ligne_courant+1, col+column_courant) = 4*alpha*col*col ;
		 systeme.set(ligne_courant+1, col+column_courant+nr) = 0 ; }
	    // La fonction vr
	    if (col%2==0) {
	          systeme.set(ligne_courant+2, col+column_courant) = -1 ;
		  systeme.set(ligne_courant+2, col+column_courant+nr) = 0 ; }
	     else {
	          systeme.set(ligne_courant+2, col+column_courant) = 1 ;
		  systeme.set(ligne_courant+2, col+column_courant+nr) = 0 ; }
	     // Sa derivee :
	     if (col%2==0) {
	         systeme.set(ligne_courant+3, col+column_courant) = -4*alpha*col*col ;
		 systeme.set(ligne_courant+3, col+column_courant+nr) = 0 ; }
	     else {
	         systeme.set(ligne_courant+3, col+column_courant) = 4*alpha*col*col ;
		 systeme.set(ligne_courant+3, col+column_courant+nr) = 0 ; }
	  }
	  ligne_courant += 4 ;	
	  
	  // Regularity and BC at infinity ?
	  nbr_cl =0 ;
	  switch (dzpuis) {
	       case 4 : 
	           if (l_quant==0) {
		       nbr_cl = 1 ;
		       // Only BC at infinity :
		       for (int col=0 ; col<nr ; col++) {
		           systeme.set(ligne_courant, col+column_courant) = 1 ;
			   systeme.set(ligne_courant, col+column_courant+nr) = 0 ;
			   systeme.set(ligne_courant+1, col+column_courant) = 0 ;
			   systeme.set(ligne_courant+1, col+column_courant+nr) = 1 ; }
		       }
		   else { 
		       nbr_cl = 2 ;
		       // BC at infinity :
		       for (int col=0 ; col<nr ; col++) {
		           systeme.set(ligne_courant, col+column_courant) = 1 ;
			   systeme.set(ligne_courant, col+column_courant+nr) = 0 ;
			   systeme.set(ligne_courant+1, col+column_courant) = 0 ;
			   systeme.set(ligne_courant+1, col+column_courant+nr) = 1; } 
		       // Regularity :
		       for (int col=0 ; col<nr ; col++) {
		           systeme.set(ligne_courant+2, col+column_courant) = -4*alpha*col*col ;
			   systeme.set(ligne_courant+2, col+column_courant+nr) = 0 ;
			   systeme.set(ligne_courant+3, col+column_courant) = 0 ;
			   systeme.set(ligne_courant+3, col+column_courant+nr) = -4*alpha*col*col ; }
		    }
		    break ;
	       
	        case 3 :
		    nbr_cl = 1 ;
		    // Only BC at infinity :
		    for (int col=0 ; col<nr ; col++) {
		       systeme.set(ligne_courant, col+column_courant) = 1 ;
		       systeme.set(ligne_courant, col+column_courant+nr) = 0 ;
		       systeme.set(ligne_courant+1, col+column_courant) = 0 ;
		       systeme.set(ligne_courant+1, col+column_courant+nr) = 1 ; }
		    break ;
		
		case 2 :
		     if (l_quant==0) {
		         nbr_cl = 1 ;
		        // Only BC at infinity :
		        for (int col=0 ; col<nr ; col++) {
		           systeme.set(ligne_courant, col+column_courant) = 1 ;
			   systeme.set(ligne_courant, col+column_courant+nr) = 0 ;
			   systeme.set(ligne_courant+1, col+column_courant) = 0 ;
			   systeme.set(ligne_courant+1, col+column_courant+1) = 1 ; }
			}
		     break ;
		default : 
		    cout << "Unknown dzpuis in vector_divfree_A ..." << endl ;
		    abort() ;
	}
	
	ligne_courant += nbr_cl ;
	
	// Multiplication of the source :
	double indic = 1 ;
	switch (dzpuis) {
	    case 4 : 
	        indic = alpha*alpha ;
		break ;
	    case 3 :
	        indic = alpha ;
		break ;
	default : 
	     break ;
	}
	
	// L'operateur :
	for (int lig=0 ; lig<nr-1-nbr_cl ; lig++) {
	    for (int col=0 ; col<nr ; col++) {
	       systeme.set(lig+ligne_courant,col+column_courant) = mdzec(lig,col)+mszec(lig,col) ;
	       systeme.set(lig+ligne_courant,col+column_courant+nr) = -mszec(lig,col) ;
	       sec_membre.set(lig+ligne_courant) = indic*(*source.get_spectral_va().c_cf)(nz-1, k, j, lig) ;
	       systeme.set(lig+ligne_courant+nr-1-nbr_cl,col+column_courant)= - l_quant*(l_quant+1)*mszec(lig,col) ;
	       systeme.set(lig+ligne_courant+nr-1-nbr_cl,col+column_courant+nr) = mdzec(lig,col)+2*mszec(lig,col) ;
	       sec_membre.set(lig+ligne_courant+nr-1-nbr_cl) = 0 ; }
	}
	

	// Solving the system:
	systeme.set_band (max_nr, max_nr) ;
	systeme.set_lu() ;
	Tbl soluce (systeme.inverse(sec_membre)) ;
	
	// On range :
	int conte = 0 ;
	for (int l=0 ; l<nz ; l++) {
	     nr = mgrid.get_nr(l) ;
	     for (int i=0 ; i<nr ; i++) {
	         meta.set(l,k,j,i) = soluce(conte) ;
		 mvr.set(l,k,i,j) = soluce(conte+nr) ;
		 conte ++ ;
		}
	}
}
if (eta_tilde.set_spectral_va().c != 0x0) 
	delete eta_tilde.set_spectral_va().c ;
    eta_tilde.set_spectral_va().c = 0x0 ;
    eta_tilde.set_spectral_va().ylm_i() ;

    if (vr.set_spectral_va().c != 0x0) 
	delete vr.set_spectral_va().c ;
    vr.set_spectral_va().c = 0x0 ;
    vr.set_spectral_va().ylm_i() ;
}
}
