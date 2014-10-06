/*
 *  Methods of class Bin_hor
 *
 *   (see file bin_hor.h for documentation)
 *
 */

/*
 *   Copyright (c) 2004-2005 Francois Limousin
 *                           Jose Luis Jaramillo
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


char bin_hor_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.11  2014/10/06 15:13:00  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.10  2007/04/13 15:28:55  f_limousin
 * Lots of improvements, generalisation to an arbitrary state of
 * rotation, implementation of the spatial metric given by Samaya.
 *
 * Revision 1.9  2006/08/01 14:37:19  f_limousin
 * New version
 *
 * Revision 1.8  2006/06/29 08:51:00  f_limousin
 * *** empty log message ***
 *
 * Revision 1.7  2006/06/28 13:36:09  f_limousin
 * Convergence to a given irreductible mass
 *
 * Revision 1.6  2006/05/24 16:56:37  f_limousin
 * Many small modifs.
 *
 * Revision 1.5  2005/06/13 15:47:29  jl_jaramillo
 * Add some quatities in write_global()
 *
 * Revision 1.4  2005/06/09 16:12:04  f_limousin
 * Implementation of amg_mom_adm().
 *
 * Revision 1.3  2005/04/29 14:02:44  f_limousin
 * Important changes : manage the dependances between quantities (for
 * instance psi and psi4). New function write_global(ost).
 *
 * Revision 1.2  2005/03/04 09:38:41  f_limousin
 * Implement the constructor from a file, operator>>, operator<<
 * and function sauve.
 *
 * Revision 1.1  2004/12/29 16:11:02  f_limousin
 * First version
 *
 *
 * $Header$
 *
 */

//standard
#include <cstdlib>
#include <cmath>

// Lorene
#include "nbr_spx.h"
#include "tenseur.h"
#include "tensor.h"
#include "isol_hor.h"
#include "proto.h"
#include "utilitaires.h"
//#include "graphique.h"

// Standard constructor
// --------------------

Bin_hor::Bin_hor (Map_af& mp1, Map_af& mp2) :
	hole1(mp1), hole2(mp2), omega(0){

    holes[0] = &hole1 ;
    holes[1] = &hole2 ;
}

// Copy constructor
// ----------------

Bin_hor::Bin_hor (const Bin_hor& source) :
	    hole1(source.hole1), hole2(source.hole2), omega(source.omega) {
    
    holes[0] = &hole1 ;
    holes[1] = &hole2 ;
    }

// Constructor from a file
// -----------------------
    
Bin_hor::Bin_hor(Map_af& mp1, Map_af& mp2, FILE* fich)
    : hole1(mp1, fich),
      hole2(mp2, fich),
      omega(0) {

    fread_be(&omega, sizeof(double), 1, fich) ;
    holes[0] = &hole1 ;
    holes[1] = &hole2 ;

}

			    //--------------//
			    //  Destructor  //
			    //--------------//

Bin_hor::~Bin_hor () {
}

                    //-----------------------//
                    // Mutators / assignment //
                    //-----------------------//

void Bin_hor::operator= (const Bin_hor& source) {    
    hole1 = source.hole1 ;
    hole2 = source.hole2 ;
    
    omega = source.omega ;
}


                //--------------------------//
                //      Save in a file      //
                //--------------------------//

void Bin_hor::sauve(FILE* fich) const {

    hole1.sauve(fich) ;
    hole2.sauve(fich) ;
    fwrite_be (&omega, sizeof(double), 1, fich) ;
   
}


//Initialisation : Sum of two static BH
void Bin_hor::init_bin_hor() {
    set_omega (0) ;
    hole1.init_bhole() ;
    hole2.init_bhole() ;
    
    hole1.psi_comp_import(hole2) ;
    hole2.psi_comp_import(hole1) ;
    
    hole1.n_comp_import(hole2) ;
    hole2.n_comp_import(hole1) ;
    
    decouple() ;
    extrinsic_curvature() ;

}


void Bin_hor::write_global(ostream& ost, double lim_nn, int bound_nn,
			   int bound_psi, int bound_beta, double alpha) const {

  double distance = hole1.get_mp().get_ori_x() - hole2.get_mp().get_ori_x() ;
  double mass_adm = adm_mass() ;
  double mass_komar = komar_mass() ;
  double mass_area = sqrt(hole1.area_hor()/16/M_PI) + 
      sqrt(hole2.area_hor()/16/M_PI) ;
  double J_adm = ang_mom_adm() ;
  double J_hor = ang_mom_hor() ; //hole1.ang_mom_hor() + hole2.ang_mom_hor() ;
  double j1 = hole1.ang_mom_hor() ;
  double j2 = hole2.ang_mom_hor() ;
  double mass_ih1 = hole1.mass_hor() ;
  double mass_ih2 = hole2.mass_hor() ;
  double mass_ih = mass_ih1 + mass_ih2 ;
  double omega1 = hole1.omega_hor() ;
  double omega2 = hole2.omega_hor() ;

  // Verification of Smarr :
  // -----------------------

    // Les alignemenents pour le signe des CL.
  double orientation1 = hole1.mp.get_rot_phi() ;
  assert ((orientation1 == 0) || (orientation1 == M_PI)) ;
  int aligne1 = (orientation1 == 0) ? 1 : -1 ;
  
  Vector angular1 (hole1.mp, CON, hole1.mp.get_bvect_cart()) ;
  Scalar yya1 (hole1.mp) ;
  yya1 = hole1.mp.ya ;
  Scalar xxa1 (hole1.mp) ;
  xxa1 = hole1.mp.xa ;
  
  angular1.set(1) = aligne1 * omega * yya1 ;
  angular1.set(2) = - aligne1 * omega * xxa1 ;
  angular1.set(3).annule_hard() ;
  
  angular1.set(1).set_spectral_va()
    .set_base(*(hole1.mp.get_mg()->std_base_vect_cart()[0])) ;
  angular1.set(2).set_spectral_va()
    .set_base(*(hole1.mp.get_mg()->std_base_vect_cart()[1])) ;
  angular1.set(3).set_spectral_va()
    .set_base(*(hole1.mp.get_mg()->std_base_vect_cart()[2])) ;
  
  angular1.change_triad(hole1.mp.get_bvect_spher()) ;


  double orientation2 = hole2.mp.get_rot_phi() ;
  assert ((orientation2 == 0) || (orientation2 == M_PI)) ;
  int aligne2 = (orientation2 == 0) ? 1 : -1 ;
  
  Vector angular2 (hole2.mp, CON, hole2.mp.get_bvect_cart()) ;
  Scalar yya2 (hole2.mp) ;
  yya2 = hole2.mp.ya ;
  Scalar xxa2 (hole2.mp) ;
  xxa2 = hole2.mp.xa ;
  
  angular2.set(1) = aligne2 * omega * yya2 ;
  angular2.set(2) = - aligne2 * omega * xxa2 ;
  angular2.set(3).annule_hard() ;
  
  angular2.set(1).set_spectral_va()
    .set_base(*(hole2.mp.get_mg()->std_base_vect_cart()[0])) ;
  angular2.set(2).set_spectral_va()
    .set_base(*(hole2.mp.get_mg()->std_base_vect_cart()[1])) ;
  angular2.set(3).set_spectral_va()
    .set_base(*(hole2.mp.get_mg()->std_base_vect_cart()[2])) ;
  
  angular2.change_triad(hole2.mp.get_bvect_spher()) ;


  Scalar btilde1 (hole1.b_tilde() - 
 contract(angular1, 0, hole1.tgam.radial_vect().up_down(hole1.tgam), 0)) ;
  Scalar btilde2 (hole2.b_tilde() - 
 contract(angular2, 0, hole2.tgam.radial_vect().up_down(hole2.tgam), 0)) ;
  



  Vector integrand_un (hole1.mp, COV, hole1.mp.get_bvect_spher()) ;
  integrand_un = hole1.nn.derive_cov(hole1.ff)*pow(hole1.psi, 2)
    - btilde1*contract(hole1.get_k_dd(), 1,
			   hole1.tgam.radial_vect(), 0)*pow(hole1.psi, 2) ;
  integrand_un.std_spectral_base() ;
 
  Vector integrand_deux (hole2.mp, COV, hole2.mp.get_bvect_spher()) ;
  integrand_deux = hole2.nn.derive_cov(hole2.ff)*pow(hole2.psi, 2)
    - btilde2*contract(hole2.get_k_dd(), 1,
			   hole2.tgam.radial_vect(), 0)*pow(hole2.psi, 2) ;
  integrand_deux.std_spectral_base() ;
 
  double horizon = hole1.mp.integrale_surface(integrand_un(1),
					    hole1.get_radius())+
    hole2.mp.integrale_surface(integrand_deux(1), hole2.get_radius()) ;

  horizon /= 4*M_PI ;

  double J_smarr = (mass_komar - horizon) / 2. / omega ;

  ost.precision(8) ;
  ost << "# Grid : " << hole1.mp.get_mg()->get_nr(1) << "x" 
      << hole1.mp.get_mg()->get_nt(1) << "x" 
      << hole1.mp.get_mg()->get_np(1) << "    R_out(l) : " ;
      
  for (int ll=0; ll<hole1.mp.get_mg()->get_nzone(); ll++) {
    ost << " " << hole1.mp.val_r(ll, 1., M_PI/2, 0) ; 
  }
  ost << endl ; 
  ost << "# bound N, lim N : " << bound_nn << " " << lim_nn 
      << " - bound Psi : " << bound_psi << " - bound shift : " << bound_beta
      << " alpha = " << alpha << endl ;

  ost << "# distance  omega  Mass_ADM  Mass_K  M_area  J_ADM  J_hor" << endl ;
  ost << distance << " " ;
  ost << omega << " " ;
  ost << mass_adm << " " ;
  ost << mass_komar << " " ;
  ost << mass_area << " " ;
  ost << J_adm << " " ;
  ost << J_hor << endl ;
  ost << "# mass_ih1  mass_ih2  mass_ih  j1  J2  omega1  omega2" << endl ;
  ost << mass_ih1 << " " ;
  ost << mass_ih2 << " " ;
  ost << mass_ih << " " ;
  ost << j1 << " " ;
  ost << j2 << " " ;
  ost << omega1 << " " ;
  ost << omega2 << endl ;
  ost << "# ADM_mass/M_area  J/M_area2  omega*M_area" << endl ;
  ost << mass_adm / mass_area << " " ;
  ost << J_adm /mass_area / mass_area << " " ;
  ost << omega * mass_area << endl ;
  ost << "# Diff J_hor/J_ADM    Diff J_ADM/J_Smarr   Diff J_hor/J_smarr" 
      << endl ;
  ost << fabs(J_adm - J_hor) / J_adm << " " <<  fabs(J_adm - J_smarr) / J_adm 
      << " " << fabs(J_hor - J_smarr) / J_hor << endl ;


}
      
