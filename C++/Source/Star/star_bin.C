/*
 * Methods for the class Star_bin
 *
 * (see file star.h for documentation)
 */

/*
 *   Copyright (c) 2004 Francois Limousin
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


char star_bin_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.11  2004/07/21 11:49:03  f_limousin
 * Remove function sprod.
 *
 * Revision 1.10  2004/06/22 12:48:52  f_limousin
 * Change qq, qq_auto and qq_comp to beta, beta_auto and beta_comp.
 *
 * Revision 1.9  2004/04/08 16:32:28  f_limousin
 * The new variable is ln(Q) instead of Q=psi^2*N. It improves the
 * convergence of the code.
 *
 * Revision 1.8  2004/03/25 10:29:26  j_novak
 * All LORENE's units are now defined in the namespace Unites (in file unites.h).
 *
 * Revision 1.7  2004/03/23 09:54:54  f_limousin
 * Add comments
 *
 * Revision 1.6  2004/02/27 09:50:57  f_limousin
 * Scalars ssjm1_logn, ssjm1_qq ... have been added for all metric
 * quantities for the resolution of Poisson equations.
 * Members bsn and d_psi are now constructed on a cartesian triad.
 *
 * Revision 1.5  2004/02/21 17:05:13  e_gourgoulhon
 * Method Scalar::point renamed Scalar::val_grid_point.
 * Method Scalar::set_point renamed Scalar::set_grid_point.
 *
 * Revision 1.4  2004/01/22 10:07:18  f_limousin
 * Add methods set_logn_comp() and set_shift_auto().
 *
 * Revision 1.3  2004/01/20 15:17:34  f_limousin
 * First version
 *
 *
 * $Header$
 *
 */

// Headers C
#include "math.h"

// Headers Lorene
#include "etoile.h" 
#include "star.h"
#include "eos.h"
#include "unites.h"	    

// Local prototype
Cmp raccord_c1(const Cmp& uu, int l1) ; 

			    //--------------//
			    // Constructors //
			    //--------------//

// Standard constructor
// --------------------
Star_bin::Star_bin(Map& mpi, int nzet_i, const Eos& eos_i, 
		     bool irrot, bool conf_flat0)
    : Star(mpi, nzet_i, eos_i), 
      irrotational(irrot), 
      psi0(mpi), 
      d_psi(mpi, COV, mpi.get_bvect_cart()), 
      wit_w(mpi, CON, mpi.get_bvect_cart()), 
      loggam(mpi), 
      bsn(mpi, CON, mpi.get_bvect_cart()), 
      pot_centri(mpi), 
      logn_auto(mpi),
      logn_comp(mpi), 
      dcov_logn(mpi, COV, mpi.get_bvect_spher()),
      dcon_logn(mpi, CON, mpi.get_bvect_spher()),
      beta_auto(mpi),
      beta_comp(mpi),
      psi4(mpi),
      dcov_lnpsi(mpi, COV, mpi.get_bvect_spher()),
      dcon_lnpsi(mpi, CON, mpi.get_bvect_spher()),
      flat(mpi, mpi.get_bvect_spher()),
      gtilde(flat),
      shift_auto(mpi, CON, mpi.get_bvect_spher()), 
      shift_comp(mpi, CON, mpi.get_bvect_spher()), 
      hij(mpi, CON, mpi.get_bvect_spher()),
      hij_auto(mpi, CON, mpi.get_bvect_spher()),
      hij_comp(mpi, CON, mpi.get_bvect_spher()), 
      tkij_auto(mpi, CON, mpi.get_bvect_spher()), 
      tkij_comp(mpi, CON, mpi.get_bvect_spher()), 
      kcar_auto(mpi), 
      kcar_comp(mpi), 
      ssjm1_logn(mpi),
      ssjm1_beta(mpi),
      ssjm1_phi(mpi),
      ssjm1_khi(mpi),
      ssjm1_mu(mpi),
      ssjm1_wshift(mpi, CON, mpi.get_bvect_cart()),
      ssjm1_h11(mpi),
      ssjm1_h21(mpi),
      ssjm1_h31(mpi),
      ssjm1_h22(mpi),
      ssjm1_h32(mpi),
      ssjm1_h33(mpi),
      decouple(mpi),
      conf_flat(conf_flat0){
    
    // Pointers of derived quantities initialized to zero : 
    set_der_0x0() ;
    
    // All quantities are initialized to zero : 
    psi0 = 0 ; 
    d_psi.set_etat_zero() ; 
    wit_w.set_etat_zero() ; 
    loggam = 0 ; 
    bsn.set_etat_zero() ; 
    pot_centri = 0 ;
  
    logn_auto = 0 ;
    logn_comp = 0 ; 
    dcov_logn.set_etat_zero() ;
    dcon_logn.set_etat_zero() ;
    shift_auto.set_etat_zero() ; 
    shift_comp.set_etat_zero() ; 
    beta_auto = 0 ;
    beta_comp = 0 ;
    psi4 = 1 ;
    dcov_lnpsi.set_etat_zero() ;
    dcon_lnpsi.set_etat_zero() ;
    hij.set_etat_zero() ;
    hij_auto.set_etat_zero() ;
    hij_comp.set_etat_zero() ;

    tkij_auto.set_etat_zero() ; 
    tkij_comp.set_etat_zero() ; 
    kcar_auto = 0 ;
    kcar_comp = 0 ; 
    ssjm1_logn = 0 ;
    ssjm1_beta = 0 ;
    ssjm1_phi = 0 ;
    ssjm1_khi = 0 ;
    ssjm1_mu = 0 ;
    ssjm1_wshift.set_etat_zero() ;
    ssjm1_h11 = 0 ;
    ssjm1_h21 = 0 ;
    ssjm1_h31 = 0 ;
    ssjm1_h22 = 0 ;
    ssjm1_h32 = 0 ;
    ssjm1_h33 = 0 ;
}

// Copy constructor
// ----------------
Star_bin::Star_bin(const Star_bin& star)
		       : Star(star), 
			 irrotational(star.irrotational), 
			 psi0(star.psi0), 
			 d_psi(star.d_psi), 
			 wit_w(star.wit_w), 
			 loggam(star.loggam), 
			 bsn(star.bsn), 
			 pot_centri(star.pot_centri), 
			 logn_auto(star.logn_auto),
			 logn_comp(star.logn_comp), 
			 dcov_logn(star.dcov_logn),
			 dcon_logn(star.dcon_logn),
			 beta_auto(star.beta_auto),
			 beta_comp(star.beta_comp),
			 psi4(star.psi4),
			 dcov_lnpsi(star.dcov_lnpsi),
			 dcon_lnpsi(star.dcon_lnpsi),
			 flat(star.flat),
			 gtilde(star.gtilde),
			 shift_auto(star.shift_auto), 
			 shift_comp(star.shift_comp), 
			 hij(star.hij),
			 hij_auto(star.hij_auto),
			 hij_comp(star.hij_comp),
			 tkij_auto(star.tkij_auto), 
			 tkij_comp(star.tkij_comp), 
			 kcar_auto(star.kcar_auto), 
			 kcar_comp(star.kcar_comp), 
			 ssjm1_logn(star.ssjm1_logn),
			 ssjm1_beta(star.ssjm1_beta),
			 ssjm1_phi(star.ssjm1_phi),
			 ssjm1_khi(star.ssjm1_khi),
			 ssjm1_mu(star.ssjm1_mu),
			 ssjm1_wshift(star.ssjm1_wshift),
			 ssjm1_h11(star.ssjm1_h11),
			 ssjm1_h21(star.ssjm1_h21),
			 ssjm1_h31(star.ssjm1_h31),
			 ssjm1_h22(star.ssjm1_h22),
			 ssjm1_h32(star.ssjm1_h32),
			 ssjm1_h33(star.ssjm1_h33),
			 decouple(star.decouple),
			 conf_flat(star.conf_flat)
{
    set_der_0x0() ;    

}    

// Constructor from a file
// -----------------------
Star_bin::Star_bin(Map& mpi, const Eos& eos_i, FILE* fich)
		       : Star(mpi, eos_i, fich), 
			 psi0(mpi), 
			 d_psi(mpi, COV, mpi.get_bvect_cart()), 
			 wit_w(mpi, CON, mpi.get_bvect_cart()), 
			 loggam(mpi), 
			 bsn(mpi, CON, mpi.get_bvect_cart()), 
			 pot_centri(mpi), 
			 logn_auto(mpi, *(mpi.get_mg()), fich),
			 logn_comp(mpi), 
			 dcov_logn(mpi, COV, mpi.get_bvect_spher()),
			 dcon_logn(mpi, CON, mpi.get_bvect_spher()),
			 beta_auto(mpi, *(mpi.get_mg()), fich),
			 beta_comp(mpi),
			 psi4(mpi),
			 dcov_lnpsi(mpi, COV, mpi.get_bvect_spher()),
			 dcon_lnpsi(mpi, CON, mpi.get_bvect_spher()),
			 flat(mpi, mpi.get_bvect_spher()),
			 gtilde(flat),
			 shift_auto(mpi, mpi.get_bvect_spher(), fich), 
			 shift_comp(mpi, CON, mpi.get_bvect_spher()), 
			 hij(mpi, CON, mpi.get_bvect_spher()),
			 hij_auto(mpi, mpi.get_bvect_spher(), fich),
			 hij_comp(mpi, CON, mpi.get_bvect_spher()),
     			 tkij_auto(mpi, CON, mpi.get_bvect_spher()), 
			 tkij_comp(mpi, CON, mpi.get_bvect_spher()), 
			 kcar_auto(mpi), 
			 kcar_comp(mpi), 
			 ssjm1_logn(mpi, *(mpi.get_mg()), fich),
			 ssjm1_beta(mpi, *(mpi.get_mg()), fich),
			 ssjm1_phi(mpi, *(mpi.get_mg()), fich),
			 ssjm1_khi(mpi, *(mpi.get_mg()), fich),
			 ssjm1_mu(mpi, *(mpi.get_mg()), fich),
			 ssjm1_wshift(mpi, mpi.get_bvect_cart(), fich),
			 ssjm1_h11(mpi, *(mpi.get_mg()), fich),
			 ssjm1_h21(mpi, *(mpi.get_mg()), fich),
			 ssjm1_h31(mpi, *(mpi.get_mg()), fich),
			 ssjm1_h22(mpi, *(mpi.get_mg()), fich),
			 ssjm1_h32(mpi, *(mpi.get_mg()), fich),
			 ssjm1_h33(mpi, *(mpi.get_mg()), fich),
			 decouple(mpi){

    // Etoile parameters
    // -----------------

    // irrotational and conf_flat is read in the file:     
    fread(&irrotational, sizeof(bool), 1, fich) ;
    fread(&conf_flat, sizeof(bool), 1, fich) ;
    	  
   
    // Read of the saved fields:
    // ------------------------

    if (irrotational) {
	Scalar gam_euler_file(mp, *(mp.get_mg()), fich) ; 
	gam_euler = gam_euler_file ; 	

	Scalar psi0_file(mp, *(mp.get_mg()), fich) ; 
	psi0 = psi0_file ; 
    }

    // All other fields are initialized to zero : 
    // ----------------------------------------

    d_psi.set_etat_zero() ; 
    wit_w.set_etat_zero() ; 
    loggam = 0 ; 
    bsn.set_etat_zero() ; 
    pot_centri = 0 ;
    logn_comp = 0 ; 
    dcov_logn.set_etat_zero() ;
    dcon_logn.set_etat_zero() ;
    shift_comp.set_etat_zero() ; 
    beta_comp = 0 ;
    psi4 = 1 ;
    dcov_lnpsi.set_etat_zero() ;
    dcon_lnpsi.set_etat_zero() ;
    hij.set_etat_zero() ;
    hij_comp.set_etat_zero() ;

    tkij_auto.set_etat_zero() ; 
    tkij_comp.set_etat_zero() ; 
    kcar_auto = 0 ;
    kcar_comp = 0 ; 
 
    // Pointers of derived quantities initialized to zero 
    // --------------------------------------------------
    set_der_0x0() ;
    
}

			    //------------//
			    // Destructor //
			    //------------//

Star_bin::~Star_bin(){

    del_deriv() ; 

}

			//----------------------------------//
			// Management of derived quantities //
			//----------------------------------//

void Star_bin::del_deriv() const {

    Star::del_deriv() ; 

    if (p_xa_barycenter != 0x0) delete p_xa_barycenter ; 
    
    set_der_0x0() ; 
}			    




void Star_bin::set_der_0x0() const {

    Star::set_der_0x0() ;

    p_xa_barycenter = 0x0 ; 

}			    

void Star_bin::del_hydro_euler() {

    Star::del_hydro_euler() ; 

    del_deriv() ; 

}			    


			    //--------------//
			    //  Assignment  //
			    //--------------//

// Assignment to another Star_bin
// --------------------------------
void Star_bin::operator=(const Star_bin& star) {

    // Assignment of quantities common to the derived classes of Star
    Star::operator=(star) ;	    

    // Assignement of proper quantities of class Star_bin
    irrotational = star.irrotational ; 
    psi0 = star.psi0 ; 
    d_psi = star.d_psi ;
    wit_w = star.wit_w ; 
    loggam = star.loggam ;
    bsn = star.bsn ;
    pot_centri = star.pot_centri ;
    logn_auto = star.logn_auto ;    
    logn_comp = star.logn_comp ;
    dcov_logn = star.dcov_logn ;
    dcon_logn = star.dcon_logn ;
    beta_auto = star.beta_auto ;
    beta_comp = star.beta_comp ;
    psi4 = star.psi4 ;
    dcov_lnpsi = star.dcov_lnpsi ;
    dcon_lnpsi = star.dcon_lnpsi ;
    flat = star.flat ;
    gtilde = star.gtilde ;
    shift_auto = star.shift_auto ;
    shift_comp = star.shift_comp ; 
    hij = star.hij ;
    hij_auto = star.hij_auto ;
    hij_comp = star.hij_comp ; 
    tkij_auto = star.tkij_auto ;
    tkij_comp = star.tkij_comp ;
    kcar_auto = star.kcar_auto ;
    kcar_comp = star.kcar_comp ;
    ssjm1_logn = star.ssjm1_logn ;
    ssjm1_beta = star.ssjm1_beta ;
    ssjm1_phi = star.ssjm1_phi ;
    ssjm1_khi = star.ssjm1_khi ;
    ssjm1_mu = star.ssjm1_mu ;
    ssjm1_wshift = star.ssjm1_wshift ;
    ssjm1_h11 = star.ssjm1_h11 ;
    ssjm1_h21 = star.ssjm1_h21 ;
    ssjm1_h31 = star.ssjm1_h31 ;
    ssjm1_h22 = star.ssjm1_h22 ;
    ssjm1_h32 = star.ssjm1_h32 ;
    ssjm1_h33 = star.ssjm1_h33 ;
    decouple = star.decouple ;
    conf_flat = star.conf_flat ;
    
    del_deriv() ;  // Deletes all derived quantities

}	

Scalar& Star_bin::set_pot_centri() {

    del_deriv() ;	// sets to 0x0 all the derived quantities
    return pot_centri ;
    
} 

Scalar& Star_bin::set_logn_comp() {

    del_deriv() ;	// sets to 0x0 all the derived quantities
    return logn_comp ;
    
} 

Vector& Star_bin::set_shift_auto() {
    
    del_deriv() ;	// sets to 0x0 all the derived quantities
    return shift_auto ;
    
} 

Vector& Star_bin::set_shift() {

    del_deriv() ;	// sets to 0x0 all the derived quantities
   return shift ;
    
} 


			    //--------------//
			    //	  Outputs   //
			    //--------------//

// Save in a file
// --------------
void Star_bin::sauve(FILE* fich) const {
    
    Star::sauve(fich) ; 
    
    logn_auto.sauve(fich) ;
    beta_auto.sauve(fich) ;
    shift_auto.sauve(fich) ;
    hij_auto.sauve(fich) ;

    ssjm1_logn.sauve(fich) ;
    ssjm1_beta.sauve(fich) ;
    ssjm1_phi.sauve(fich) ;
    ssjm1_khi.sauve(fich) ;
    ssjm1_mu.sauve(fich) ;
    ssjm1_wshift.sauve(fich) ;
    ssjm1_h11.sauve(fich) ;
    ssjm1_h21.sauve(fich) ;
    ssjm1_h31.sauve(fich) ;
    ssjm1_h22.sauve(fich) ;
    ssjm1_h32.sauve(fich) ;
    ssjm1_h33.sauve(fich) ;

    fwrite(&irrotational, sizeof(bool), 1, fich) ;		
    fwrite(&conf_flat, sizeof(bool), 1, fich) ;		
    
    if (irrotational) {
	gam_euler.sauve(fich) ; // required to construct d_psi from psi0
	psi0.sauve(fich) ; 
    }
 
}

// Printing
// --------

ostream& Star_bin::operator>>(ostream& ost) const {
    
  using namespace Unites ;

    Star::operator>>(ost) ; 
    
    ost << endl ; 
    ost << "Star in a binary system" << endl ; 
    ost << "-----------------------" << endl ; 
    
    if (irrotational) {
	ost << "irrotational configuration" << endl ; 
    }
    else {
	ost << "corotating configuration" << endl ; 
    }
       
    ost << "Absolute abscidia of the stellar center: " <<
	mp.get_ori_x() / km << " km" << endl ; 
    
    ost << "Absolute abscidia of the barycenter of the baryon density : " <<
	xa_barycenter() / km << " km" << endl ; 
    
    double r_0 = 0.5 * ( ray_eq() + ray_eq_pi() ) ; 
    double d_ns = fabs( mp.get_ori_x() ) + ray_eq_pi() - r_0 ;
    double d_tilde = 2 * d_ns / r_0 ;  
    
    ost << "d_tilde : " << d_tilde << endl ; 

    ost << "Central value of gam_euler : " 
        << gam_euler.val_grid_point(0, 0, 0, 0)  << endl ; 

    ost << "Central u_euler (U^r, U^t, U^p) [c] : " 
	<< u_euler(1).val_grid_point(0, 0, 0, 0) << "  " 
	<< u_euler(2).val_grid_point(0, 0, 0, 0) << "  " 
	<< u_euler(3).val_grid_point(0, 0, 0, 0) << endl ; 

    if (irrotational) {
    ost << "Central d_psi (r, t, p) [c] :         " 
	    << d_psi(1).val_grid_point(0, 0, 0, 0) << "  " 
	    << d_psi(2).val_grid_point(0, 0, 0, 0) << "  " 
	    << d_psi(3).val_grid_point(0, 0, 0, 0) << endl ; 

	ost << "Central vel. / co-orb. (W^r, W^t, W^p) [c] : " 
	    << wit_w(1).val_grid_point(0, 0, 0, 0) << "  " 
	    << wit_w(2).val_grid_point(0, 0, 0, 0) << "  " 
	    << wit_w(3).val_grid_point(0, 0, 0, 0) << endl ; 

	ost << "Max vel. / co-orb. (W^r, W^t, W^p) [c] : " 
	    << max(max(wit_w(1))) << "  " 
	    << max(max(wit_w(2))) << "  " 
	    << max(max(wit_w(3))) << endl ; 

	ost << "Min vel. / co-orb. (W^r, W^t, W^p) [c] : " 
	    << min(min(wit_w(1))) << "  " 
	    << min(min(wit_w(2))) << "  " 
	    << min(min(wit_w(3))) << endl ; 

	double r_surf = mp.val_r(0,1.,M_PI/4,M_PI/4) ;

	ost << "Velocity at (r_surf,pi/4,pi/4) / co-orb. [c] : "
	    << wit_w(1).val_point(r_surf,M_PI/4,M_PI/4) << "  "
	    << wit_w(2).val_point(r_surf,M_PI/4,M_PI/4) << "  "
	    << wit_w(3).val_point(r_surf,M_PI/4,M_PI/4) << endl ;

	ost << "Central value of loggam : " 
	    << loggam.val_grid_point(0, 0, 0, 0)  << endl ; 	
    }


    ost << "Central value of log(N) auto, comp :         " 
	<< logn_auto.val_grid_point(0, 0, 0, 0) << "  " 
	<< logn_comp.val_grid_point(0, 0, 0, 0) << endl ; 

    ost << "Central value of shift (N^r, N^t, N^p) [c] : " 
	<< shift(1).val_grid_point(0, 0, 0, 0) << "  " 
	<< shift(2).val_grid_point(0, 0, 0, 0) << "  " 
	<< shift(3).val_grid_point(0, 0, 0, 0) << endl ; 

    ost << "  ... shift_auto part of it [c] :            " 
	<< shift_auto(1).val_grid_point(0, 0, 0, 0) << "  " 
	<< shift_auto(2).val_grid_point(0, 0, 0, 0) << "  " 
	<< shift_auto(3).val_grid_point(0, 0, 0, 0) << endl ; 

    ost << endl << "Central value of (B^r, B^t, B^p)/N [c] : " 
	<< bsn(1).val_grid_point(0, 0, 0, 0) << "  " 
	<< bsn(2).val_grid_point(0, 0, 0, 0) << "  " 
	<< bsn(3).val_grid_point(0, 0, 0, 0) << endl ; 


    ost << endl << "Central A^2 K^{ij} [c/km] : " << endl ; 
    ost << "  A^2 K^{rr} auto, comp : " 
	<< tkij_auto(1, 1).val_grid_point(0, 0, 0, 0) * km  << "  "
	<< tkij_comp(1, 1).val_grid_point(0, 0, 0, 0) * km << endl ; 
    ost << "  A^2 K^{rt} auto, comp : " 
	<< tkij_auto(1, 2).val_grid_point(0, 0, 0, 0) * km  << "  "
	<< tkij_comp(1, 2).val_grid_point(0, 0, 0, 0) * km << endl ; 
    ost << "  A^2 K^{rp} auto, comp : " 
	<< tkij_auto(1, 3).val_grid_point(0, 0, 0, 0) * km  << "  "
	<< tkij_comp(1, 3).val_grid_point(0, 0, 0, 0) * km << endl ; 
    ost << "  A^2 K^{tt} auto, comp : " 
	<< tkij_auto(2, 2).val_grid_point(0, 0, 0, 0) * km  << "  "
	<< tkij_comp(2, 2).val_grid_point(0, 0, 0, 0) * km << endl ; 
    ost << "  A^2 K^{tp} auto, comp : " 
	<< tkij_auto(2, 3).val_grid_point(0, 0, 0, 0) * km  << "  "
	<< tkij_comp(2, 3).val_grid_point(0, 0, 0, 0) * km << endl ; 
    ost << "  A^2 K^{pp} auto, comp : " 
	<< tkij_auto(3, 3).val_grid_point(0, 0, 0, 0) * km  << "  "
	<< tkij_comp(3, 3).val_grid_point(0, 0, 0, 0) * km << endl ; 

    ost << endl << "Central A^2 K_{ij} K^{ij} [c^2/km^2] : " << endl ; 
    ost << "   A^2 K_{ij} K^{ij}  auto, comp : " 
	<< kcar_auto.val_grid_point(0, 0, 0, 0) * km*km  << "  "
	<< kcar_comp.val_grid_point(0, 0, 0, 0) * km*km << endl ; 

    
    return ost ; 
}

			    //-------------------------//
			    //	Computational routines //
			    //-------------------------//

void Star_bin::fait_d_psi() {

    if (!irrotational) {
	d_psi.set_etat_nondef() ; 
	return ; 
    }

    // Specific relativistic enthalpy		    ---> hhh
    //----------------------------------
    
    Scalar hhh = exp(ent) ;  // = 1 at the Newtonian limit
 
    //  Computation of W^i = - h Gamma_n B^i/N
    //----------------------------------------------

    Vector www = - hhh * gam_euler * bsn * psi4 ; 
      
    // Constant value of W^i at the center of the star
    //-------------------------------------------------
    
    Vector v_orb(mp, COV, mp.get_bvect_cart()) ; 
    
    for (int i=1; i<=3; i++) {
	v_orb.set(i) = www(i).val_grid_point(0, 0, 0, 0) ; 
    }
    
    // Gradient of psi 
    //----------------

    Vector d_psi0 = psi0.derive_cov(flat) ; 
    
    d_psi0.change_triad( mp.get_bvect_cart() ) ; 
    d_psi0.std_spectral_base() ;

    d_psi = d_psi0 + v_orb ; 
    for (int i=1; i<=3; i++) {
	if (d_psi(i).get_etat() == ETATZERO)
	    d_psi.set(i).annule_hard() ;
    }

    // C^1 continuation of d_psi outside the star
    //  (to ensure a smooth enthalpy field accross the stellar surface)
    // ----------------------------------------------------------------
    
    int nzm1 = mp.get_mg()->get_nzone() - 1 ;    
    d_psi.annule(nzet, nzm1) ;	 
    for (int i=1; i<=3; i++) {
	Cmp d_psi_i (d_psi(i)) ;
	d_psi_i.va.set_base( d_psi0(i).get_spectral_va().base ) ; 
	d_psi_i = raccord_c1(d_psi_i, nzet) ; 
	d_psi.set(i) = d_psi_i ; 
    }

} 


void Star_bin::relaxation(const Star_bin& star_jm1, double relax_ent, 
			    double relax_met, int mer, int fmer_met) {
				
    double relax_ent_jm1 = 1. - relax_ent ; 
    double relax_met_jm1 = 1. - relax_met ; 

    ent = relax_ent * ent + relax_ent_jm1 * star_jm1.ent ; 

    if ( (mer != 0) && (mer % fmer_met == 0)) {

	logn_auto = relax_met * logn_auto + relax_met_jm1 * star_jm1.logn_auto ;
	beta_auto = relax_met * beta_auto + relax_met_jm1 * star_jm1.beta_auto ;
	
	shift_auto = relax_met * shift_auto 
				   + relax_met_jm1 * star_jm1.shift_auto ;
	

	hij_auto = relax_met * hij_auto + relax_met_jm1 * star_jm1.hij_auto ;
	
    }

    del_deriv() ; 
    
    equation_of_state() ; 

}

void Star_bin::test_K_Hi() const {

    int nr = mp.get_mg()->get_nr(0) ;
    int nt = mp.get_mg()->get_nt(0) ;
    int np = mp.get_mg()->get_np(0) ;
    
    Sym_tensor tkij = tkij_auto + tkij_comp ;
    
    cout << "Le maximal slicing est il bien satisfait ??" 
	 << endl ;
    cout << "Tensor Kij" << endl ;
    for (int i=1; i<=3; i++)
	cout << "  Comp. 1 1 : " << norme(tkij.down(1, gamma)(i,i)
					  /(nr*nt*np)) << endl ;

    cout << endl << "Trace of Kij : "<< norme(tkij.down(1, gamma).trace()
					      /(nr*nt*np))<< endl ;
     

    cout << "La jauge de Dirac est elle bien satisfaite ??" << endl ;
    cout << "Vector Hi" << endl ;
    for (int i=1; i<=3; i++)
	cout << "  Comp. " << i << " : " << norme((gtilde.con()
				   .divergence(flat))(i)/(nr*nt*np)) << endl ;
    
    
    cout << "Pour comparaison valeur de D_i(g^1i)" << endl ;
    for (int i=1; i<=3; i++)
	cout << "  i = " << i << " : " << norme((gtilde.con().derive_cov(flat))
					      (1, i, i)/(nr*nt*np)) << endl ;
    

    
}
