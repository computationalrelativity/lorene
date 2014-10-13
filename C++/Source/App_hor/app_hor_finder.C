/*
 *  Function ah_finder
 * 
 * (see file app_hor.h for documentation)
 */

/*
 *   Copyright (c) 2005 Lap-Ming Lin & Jerome Novak
 *
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License version 2
 *   as published by the Free Software Foundation.
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

char app_hor_finder_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.11  2014/10/13 08:52:38  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.10  2014/10/06 15:12:56  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.9  2012/01/02 13:52:57  j_novak
 * New parameter 'verbose' to get less output if needed.
 *
 * Revision 1.8  2008/01/08 13:56:54  j_novak
 * Minor modif.
 *
 * Revision 1.7  2007/10/23 12:26:08  j_novak
 * Added a test for the case where there is no AH, h(theta,phi) is then going out of the grid
 *
 * Revision 1.6  2005/12/09 09:35:06  lm_lin
 *
 * Add more information to screen output if no convergence.
 *
 * Revision 1.5  2005/12/07 14:16:36  lm_lin
 *
 * Add option to turn off screen output if no horizon is found
 * (for performance reason in hydrodynamics simulation).
 *
 * Revision 1.4  2005/12/07 11:11:09  lm_lin
 *
 * Add option to turn off screen output during iterations.
 *
 * Revision 1.3  2005/11/17 15:53:28  lm_lin
 *
 * A tiny fix.
 *
 * Revision 1.2  2005/11/17 14:20:43  lm_lin
 *
 * Check the expansion function evaluated on the apparent horizon after the
 * iteration of the 2-surface converges.
 *
 * Revision 1.1  2005/10/13 08:51:15  j_novak
 * New stuff for apparent horizon finder. For the moment, there is only an
 * external function. A class should come soon...
 *
 *
 *
 * $Header$
 *
 */


// C headers
#include <cmath>
#include <cassert>

// Lorene headers
#include "app_hor.h"
#include "graphique.h"


namespace Lorene {
bool ah_finder(const Metric& gamma, const Sym_tensor& k_dd, Valeur& h, Scalar& ex_fcn,
	       double a_axis, double b_axis, double c_axis, bool verbose, bool print, 
	       double precis, double precis_exp, int step_max, int step_relax, 
	       double relax)
{

    bool ah_flag = false ;

 //Get the mapping, grid, base vector...etc
  const Map& map = gamma.get_mp() ;
  const Mg3d* mg = map.get_mg() ;
  const Mg3d* g_angu = mg->get_angu() ;
  int nz = mg->get_nzone() ;   

  const Base_vect_spher& bspher = map.get_bvect_spher() ;

  const Coord& rr = map.r ;
  const Coord& theta = map.tet ;
  const Coord& phi = map.phi ;
  const Coord& costh = map.cost ; 
  const Coord& cosph = map.cosp ;
  const Coord& sinth = map.sint ;
  const Coord& sinph = map.sinp ;

  double r_min = min(+rr)(0) ;
  double r_max = max(+rr)(nz-1) ;

  // Set up a triaxial ellipsoidal surface as the initial guess for h
  //------------------------------------------------------------------

  double aa = a_axis ; 
  double bb = b_axis ;
  double cc = c_axis ;

  Scalar ct(map) ;  
  ct = costh ;
  ct.std_spectral_base() ;

  Scalar st(map) ;
  st = sinth ;
  st.std_spectral_base() ;

  Scalar cp(map) ;
  cp = cosph ;
  cp.std_spectral_base() ;

  Scalar sp(map) ;
  sp = sinph ;
  sp.std_spectral_base() ;

  Scalar h_tmp(map) ;

  h_tmp = st*st*( cp*cp/aa/aa + sp*sp/bb/bb ) + ct*ct/cc/cc ;
  h_tmp = 1./sqrt( h_tmp ) ;
  h_tmp.std_spectral_base() ;

  Valeur h_guess(g_angu) ;
  h_guess.annule_hard() ;

  for (int l=0; l<nz; l++) {
    
    int jmax = mg->get_nt(l) ;
    int kmax = mg->get_np(l) ;

    for (int k=0; k<kmax; k++) {
      for (int j=0; j<jmax; j++) {		      
	  h_guess.set(l,k,j,0)  = h_tmp.val_grid_point(l,k,j,0) ; 		      
      }
    }
  }

  h_guess.std_base_scal() ;
 
  h = h_guess ;    // initialize h to h_guess 
  h.std_base_scal() ;

  //End setting initial guess for h
  //------------------------------------

  const Metric_flat& fmets = map.flat_met_spher() ; 


  // Define the conformal factor                           
  double  one_third = double(1) / double(3) ;

  Scalar psi4 = gamma.determinant() / fmets.determinant() ;
  psi4 = pow( psi4, one_third ) ;      
  psi4.std_spectral_base() ;

  // The expansion function at the n-1 iteration step
  Scalar ex_fcn_old(map) ;
  ex_fcn_old.set_etat_zero() ;

  // The expansion function evaulated on the 2 surface h
  Valeur ex_AH(g_angu) ;
  ex_AH.annule_hard() ;

  // Normal unit vector of the level set surface F = r - h(\theta,phi)
  Vector s_unit(map, CON, bspher) ; 

  double relax_prev = double(1) - relax ; 
  double diff_exfcn = 1. ;
  Tbl diff_h(nz) ;
  diff_h = 1. ;

  bool no_AH_in_grid = false ;

  //--------------------------------------------------------
  //         Start of iteration
  //--------------------------------------------------------

  for (int step=0 ; 
       (max(diff_h) > precis) && (step < step_max) && (!no_AH_in_grid);
       step++) {
      
      
      // ***To be fixed: the function "set_grid_point" does not delete the derived 
      //                 quantities of F.
      // Temporary fix: Define the level set F inside the iteration loop...
      //----------------------------------------------------------------
      
      Scalar F(map) ;  // level set function: F = r - h(theta,phi) 
      F.allocate_all() ;

      for (int l=0; l<nz; l++) {
	  
	  int imax = mg->get_nr(l) ;
	  int jmax = mg->get_nt(l) ;
	  int kmax = mg->get_np(l) ;
	  
	  for (int k=0; k<kmax; k++) {
	      for (int j=0; j<jmax; j++) {
		  for (int i=0; i<imax; i++) {
		      
		      // (+rr) converts rr to Mtbl
		      F.set_grid_point(l,k,j,i) = (+rr)(l,k,j,i) - h(l,k,j,0) ; 
		      
		  }
	      }
	  }
      }
      
      F.std_spectral_base() ;
      
      // Construct the unit normal vector s^i of the surface F
      Scalar dF_norm(map) ;   
      dF_norm = contract( contract(gamma.con(), 0, F.derive_cov(gamma), 0), 
			0, F.derive_cov(gamma), 0) ;
      dF_norm = sqrt( dF_norm ) ;
      dF_norm.std_spectral_base() ;
      
      s_unit = F.derive_con(gamma) / dF_norm ; 
      
      // The expansion function
      ex_fcn = s_unit.divergence(gamma) - k_dd.trace(gamma) + 
      contract( s_unit, 0, contract(s_unit, 0, k_dd, 1), 0) ; 
      
      // Construct the source term for the angular Laplace equation
      //---------------------------------------------------------

      Sym_tensor sou_1(map, CON, bspher) ;
      sou_1 = gamma.con() - fmets.con()/psi4  - s_unit*s_unit  ;
      
      Scalar source_tmp(map) ;
      source_tmp = contract( sou_1, 0, 1, F.derive_cov(fmets).derive_cov(fmets), 0, 1 ) ;
      source_tmp = source_tmp / dF_norm ;
      
      Sym_tensor d_gam(map, COV, bspher) ;
      d_gam = contract( gamma.connect().get_delta(), 0, F.derive_cov(fmets), 0) ; 
      
      source_tmp -= contract( gamma.con() - s_unit*s_unit, 0, 1,
			      d_gam, 0, 1) / dF_norm ; 
    
      source_tmp = psi4*dF_norm*( source_tmp - k_dd.trace(gamma) + 
				contract(s_unit, 0, contract(s_unit, 0, k_dd, 1), 0) ) ;
      
      source_tmp.std_spectral_base() ;
      
      
      Valeur sou_angu(g_angu) ; // source defined on the angular grid 
      // S(theta, phi) = S(h(theta,phi),theta,phi) 
      sou_angu.annule_hard() ;
      
      double h_min = min(h)(0) ;
      double h_max = max(h)(0) ;
      if ( (r_min < h_min) && (h_max < r_max) ) { 
	  
	  for (int l=0; l<nz; l++) {
	      
	      int jmax = mg->get_nt(l) ;
	      int kmax = mg->get_np(l) ;
	      for (int k=0; k<kmax; k++) {
		  for (int j=0; j<jmax; j++) {
		    sou_angu.set(l,k,j,0) = source_tmp.val_point(h(l,k,j,0)
					  ,(+theta)(l,k,j,0) ,(+phi)(l,k,j,0)) ;  
		}
	      }
	  }
	  sou_angu = h*h*sou_angu  ; // Final source term: psi4*dF_norm*h^2*(source_tmp)  
      } 
      else {
	  no_AH_in_grid = true ;
	  break ;
      }
      sou_angu.std_base_scal() ; 
      
      
      // Done with the source term
      //-----------------------------------------
      
      
      // Start solving the equation L^2h - 2h = source
      //-----------------------------------------------
      
      sou_angu.ylm() ;
      
      Valeur h_new = sou_angu ;
      
      h_new.c_cf->poisson_angu(-2.) ; 
      
      h_new.ylm_i() ;      
      
      if (h_new.c != 0x0)
	  delete h_new.c ;
      h_new.c = 0x0 ;
      h_new.coef_i() ;
      
      // Convergence condition:  
      diff_h = max(abs(h - h_new)) ;
      
      
      // Relaxations
      if (step >= step_relax) {
	  h_new = relax * h_new + relax_prev * h ;
      }
      
      // Recycling for the next step
      h = h_new ;

      
      if (print) 
      {
	  
	  cout << "-------------------------------------" << endl ;
	  cout << "App_hor iteration step: " << step << endl ;
	  cout << "    " << endl ;
	  
	  cout << "Difference in h : " << diff_h << endl ;
	  
	  // Check: calculate the difference between ex_fcn and ex_fcn_old
	  Tbl diff_exfcn_tbl = diffrel( ex_fcn, ex_fcn_old ) ;
	  diff_exfcn = diff_exfcn_tbl(0) ;
	  for (int l=1; l<nz; l++) {
	      diff_exfcn += diff_exfcn_tbl(l) ;
	  }
	  diff_exfcn /= nz ;
	  cout << "diff_exfcn : " << diff_exfcn << endl ;
	  
	  ex_fcn_old = ex_fcn ; // recycling 
	  // End check
    
      }
      
      if ( (step == step_max-1) && (max(diff_h) > precis) ) {
	  
	  
	  //Check: Evaluate the expansion function on the 2-surface
	  
      for (int l=0; l<nz; l++) {
	  
	  int jmax = mg->get_nt(l) ;
	  int kmax = mg->get_np(l) ;
	  
	  for (int k=0; k<kmax; k++) {
	      for (int j=0; j<jmax; j++) {
		  
		  ex_AH.set(l,k,j,0) = ex_fcn.val_point(h(l,k,j,0),(+theta)(l,k,j,0)
							,(+phi)(l,k,j,0))  ;  
	      }
	  }
      }
      
      if (verbose) {
	cout << " " << endl ;
	cout << "###############################################" << endl ;
	cout << "AH finder: maximum number of iteration reached!" << endl ;
	cout << "      No convergence in the 2-surface h!      " << endl ;
	cout << " max( difference in h ) > prescribed tolerance " << endl ;
	cout << " " << endl ;
	cout << " prescribed tolerance = " << precis << endl ;
	cout << " max( difference in h ) = " << max(diff_h) << endl ;
	cout << " max( expansion function on h ) = " << max(abs(ex_AH(0))) << endl ; 
	cout << "###############################################" << endl ;
	cout << " " << endl ;
	
      }
      }
      
      
  } // End of iteration
  
  //Done with the AH finder
  

  
  //Check: Evaluate the expansion function on the 2-surface
  
  if (no_AH_in_grid) {
      if (print) {
      cout << " " << endl ;
      cout << "###############################################" << endl ;
      cout << " AH finder: no horizon found inside the grid!" << endl ;
      cout << " Grid: rmin= " << r_min << ", rmax= " << r_max << endl ;
      cout << "###############################################" << endl ;
      cout << " " << endl ;
      }
  }
  else {
  for (int l=0; l<nz; l++) {
      
      int jmax = mg->get_nt(l) ;
      int kmax = mg->get_np(l) ;
      
      for (int k=0; k<kmax; k++) {
	for (int j=0; j<jmax; j++) {
	    
	    ex_AH.set(l,k,j,0) = ex_fcn.val_point(h(l,k,j,0),(+theta)(l,k,j,0)
						  ,(+phi)(l,k,j,0))  ;  
	}
      }
  }
  
  
  
  if ( (max(diff_h) < precis) && (max(abs(ex_AH(0))) < precis_exp) ) {
      
      ah_flag = true ; 

      if (verbose) {
	cout << "  " << endl ;
	cout << "################################################" << endl ;
	cout << " AH finder: Apparent horizon found!!!      " << endl ;
	cout << " Max error of the expansion function on h: " << endl ;
	cout << " max( expansion function on AH ) = " << max(abs(ex_AH(0))) << endl ; 
	cout << "################################################" << endl ;
	cout << " " << endl ;
      }
      
  }
  
  if ( (max(diff_h) < precis) && (max(abs(ex_AH(0))) > precis_exp) ) {

      if (print) {
	  cout << " " << endl ;
	  cout << "#############################################" << endl ;
	  cout << " AH finder: convergence in the 2 surface h.   " << endl ;
	  cout << " But max error of the expansion function evaulated on h > precis_exp"  << endl ;
	  cout << "   max( expansion function on AH ) =  " << max(abs(ex_AH(0))) << endl ;
	  cout << " Probably not an apparent horizon! " << endl ;
	  cout << "#############################################" << endl ;
	  cout << "   " << endl ;
      }

  }
  }
  return ah_flag ;
  
} // End ah_finder

}
