/*
 * Methods of class Star
 *
 * (see file star.h for documentation)
 */

/*
 *   Copyright (c) 2004 Francois Limousin
 *  
 *   Copyright (c) 2000-2001 Eric Gourgoulhon (for preceding class Etoile)
 *   Copyright (c) 2000-2001 Keisuke Taniguchi (for preceding class Etoile)
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


char star_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.10  2004/06/22 12:48:08  f_limousin
 * Change qq, qq_auto and qq_comp to beta, beta_auto and beta_comp.
 *
 * Revision 1.9  2004/06/07 16:21:35  f_limousin
 * Add outputs
 *
 * Revision 1.8  2004/04/08 16:32:10  f_limousin
 * The new variable is ln(Q) instead of Q=psi^2*N. It improves the
 * convergence of the code.
 *
 * Revision 1.7  2004/03/25 10:29:26  j_novak
 * All LORENE's units are now defined in the namespace Unites (in file unites.h).
 *
 * Revision 1.6  2004/03/08 11:48:00  f_limousin
 * Error in del_deriv() and set_der_0x0() : p_mass_b and p_mass_g were
 * missing. And so they were never recomputed.
 *
 * Revision 1.5  2004/02/27 09:36:46  f_limousin
 * u_euler is now constructed on a cartesian basis instead
 * of a spherical one.
 *
 * Revision 1.4  2004/02/21 17:05:13  e_gourgoulhon
 * Method Scalar::point renamed Scalar::val_grid_point.
 * Method Scalar::set_point renamed Scalar::set_grid_point.
 *
 * Revision 1.3  2004/01/20 15:16:58  f_limousin
 * First version
 *
 *
 * $Header$
 *
 */

// Headers C
#include "math.h"

// Headers Lorene
#include "star.h"
#include "eos.h"
#include "utilitaires.h"
#include "param.h"
#include "unites.h"


			    //--------------//
			    // Constructors //
			    //--------------//

// Standard constructor
// --------------------
Star::Star(Map& mpi, int nzet_i, const Eos& eos_i)
		 : mp(mpi), 
		   nzet(nzet_i), 
		   eos(eos_i), 
		   ent(mpi), 
		   nbar(mpi), 
		   ener(mpi), 
		   press(mpi),  
		   ener_euler(mpi), 
		   s_euler(mpi), 
		   gam_euler(mpi), 
		   u_euler(mpi, CON, mp.get_bvect_cart()), 
		   stress_euler(mpi, 2, CON, mp.get_bvect_spher()),
		   logn(mpi), 
		   nnn(mpi), 
		   shift(mpi, CON, mp.get_bvect_spher()),
		   beta(mpi),
		   gamma(mp.flat_met_spher()){
    
 
    // Check of the EOS
    const Eos_poly_newt* p_eos_poly_newt = 
			    dynamic_cast<const Eos_poly_newt*>( &eos ) ; 
	  
    const Eos_incomp_newt* p_eos_incomp_newt = 
			    dynamic_cast<const Eos_incomp_newt*>( &eos ) ; 
	  

    
    if (p_eos_poly_newt != 0x0) {
	cout << 
	    "Star::Star : the EOS Eos_poly_newt must not be employed"
	     << " for a relativistic star ! " << endl ; 
	cout << "(Use Eos_poly instead)" << endl ; 
	abort() ; 
	}
    if (p_eos_incomp_newt != 0x0) {
	cout << 
	    "Star::Star : the EOS Eos_incomp_newt must not be employed"
	     << " for a relativistic star ! " << endl ; 
	cout << "(Use Eos_incomp instead)" << endl ; 
	abort() ; 
    }
 
    // Pointers of derived quantities initialized to zero : 
    set_der_0x0() ;

    // All the matter quantities are initialized to zero :
    nbar = 0 ; 
    ener = 0 ; 
    press = 0 ; 
    ent = 0 ; 
    ener_euler = 0 ; 
    s_euler = 0 ; 
    gam_euler = 1. ; 
    gam_euler.std_spectral_base() ; 
    u_euler.set_etat_zero() ; 
    stress_euler.set_etat_zero() ;

    // The metric is initialized to the flat one : 
    Metric_flat flat(mp.flat_met_spher()) ;
    flat.cov() ;
    gamma = flat ;

    logn = 0 ; 
    nnn = 1. ; 
    nnn.std_spectral_base() ; 
    shift.set_etat_zero() ; 
    beta = 0 ;

}

// Copy constructor
// ----------------
Star::Star(const Star& et) 
		 : mp(et.mp), 
		   nzet(et.nzet), 
		   eos(et.eos), 
		   ent(et.ent), 
		   nbar(et.nbar), 
		   ener(et.ener), 
		   press(et.press),  
		   ener_euler(et.ener_euler), 
		   s_euler(et.s_euler), 
		   gam_euler(et.gam_euler), 
		   u_euler(et.u_euler), 
		   stress_euler(et.stress_euler),
		   logn(et.logn), 
		   nnn(et.nnn), 
		   shift(et.shift),
		   beta(et.beta),
		   gamma(et.gamma){
	       
    set_der_0x0() ;

}

// Constructor from a file
// -----------------------
Star::Star(Map& mpi, const Eos& eos_i, FILE* fich)
		 : mp(mpi), 
		   eos(eos_i), 
		   ent(mpi), 
		   nbar(mpi), 
		   ener(mpi), 
		   press(mpi),  
		   ener_euler(mpi), 
		   s_euler(mpi), 
		   gam_euler(mpi), 
		   u_euler(mpi, CON, mp.get_bvect_cart()), 
		   stress_euler(mpi, 2, CON, mp.get_bvect_spher()), 
		   logn(mpi, *(mpi.get_mg()), fich), 
		   nnn(mpi), 
		   shift(mpi, CON, mp.get_bvect_spher()),
		   beta(mpi, *(mpi.get_mg()), fich),
		   gamma(mpi.flat_met_spher()){

    // Star parameters
    // -----------------

    // nzet is read in the file:     
    int xx ; 
    fread_be(&xx, sizeof(int), 1, fich) ;	
    nzet = xx ;
    		
    // Equation of state
    // -----------------
    
    // Read of the saved EOS
    Eos* p_eos_file = Eos::eos_from_file(fich) ; 
    
    // Comparison with the assigned EOS:
    if (eos != *p_eos_file) {
	cout << 
	"Star::Star(const Map&, const Eos&, FILE*) : the EOS given in "
	<< endl << 
	" argument and that read in the file are different !" << endl ; 
	abort() ;  
    }
    
    // p_eos_file is no longer required (it was used only for checking the
    //  EOS compatibility)
    delete p_eos_file ;
    
    // Read of the saved fields:
    // ------------------------
    Scalar ent_file(mp, *(mp.get_mg()), fich) ; 
    ent = ent_file ; 
    u_euler.set_etat_zero() ; 
    stress_euler.set_etat_zero() ;
    nnn = 1 ;
    shift.set_etat_zero() ;

    // Pointers of derived quantities initialized to zero 
    // --------------------------------------------------
    set_der_0x0() ;
    
}

			    //------------//
			    // Destructor //
			    //------------//

Star::~Star(){

    del_deriv() ; 

}


			//----------------------------------//
			// Management of derived quantities //
			//----------------------------------//

void Star::del_deriv() const {

    if (p_mass_b != 0x0) delete p_mass_b ; 
    if (p_mass_g != 0x0) delete p_mass_g ; 
    if (p_ray_eq != 0x0) delete p_ray_eq ; 
    if (p_ray_eq_pis2 != 0x0) delete p_ray_eq_pis2 ; 
    if (p_ray_eq_pi != 0x0) delete p_ray_eq_pi ; 
    if (p_ray_pole != 0x0) delete p_ray_pole ; 
    if (p_l_surf != 0x0) delete p_l_surf ; 
    if (p_xi_surf != 0x0) delete p_xi_surf ; 

    Star::set_der_0x0() ; 
}			    




void Star::set_der_0x0() const {

    p_mass_b = 0x0 ; 
    p_mass_g = 0x0 ; 
    p_ray_eq = 0x0 ; 
    p_ray_eq_pis2 = 0x0 ; 
    p_ray_eq_pi = 0x0 ; 
    p_ray_pole = 0x0 ; 
    p_l_surf = 0x0 ; 
    p_xi_surf = 0x0 ; 

}			    

void Star::del_hydro_euler() {

    ener_euler.set_etat_nondef() ; 
    s_euler.set_etat_nondef() ; 
    gam_euler.set_etat_nondef() ; 
    u_euler.set_etat_nondef() ; 
    stress_euler.set_etat_nondef() ;

    del_deriv() ; 

}			    




			    //--------------//
			    //  Assignment  //
			    //--------------//

// Assignment to another Star
// ----------------------------
void Star::operator=(const Star& et) {

    assert( &(et.mp) == &mp ) ;		    // Same mapping
    assert( &(et.eos) == &eos ) ;	    // Same EOS
    
    nzet = et.nzet ; 
    ent = et.ent ;
    nbar = et.nbar ; 
    ener = et.ener ;
    press = et.press ;
    ener_euler = et.ener_euler ;
    s_euler = et.s_euler ;
    gam_euler = et.gam_euler ;
    u_euler = et.u_euler ;
    stress_euler = et.stress_euler ;
    logn = et.logn ;
    nnn = et.nnn ;
    shift = et.shift ;
    beta = et.beta ;
    gamma = et.gamma ;

    del_deriv() ;  // Deletes all derived quantities

}	

// Assignment of the enthalpy field
// --------------------------------

void Star::set_enthalpy(const Scalar& ent_i) {
    
    ent = ent_i ; 
    
    // Update of (nbar, ener, press) :
    equation_of_state() ; 
    
    // The derived quantities are obsolete:
    del_deriv() ; 
    
}

			    //--------------//
			    //	  Outputs   //
			    //--------------//

// Save in a file
// --------------
void Star::sauve(FILE* fich) const {
    
    logn.sauve(fich) ;
    beta.sauve(fich) ;

    int xx = nzet ;     
    fwrite_be(&xx, sizeof(int), 1, fich) ;			

    eos.sauve(fich) ; 
    ent.sauve(fich) ;     
}

// Printing
// --------

ostream& operator<<(ostream& ost, const Star& et)  {
    et >> ost ;
    return ost ;
}
    
ostream& Star::operator>>(ostream& ost) const {
    
  using namespace Unites ;

    ost << endl ; 
      
    ost << "Number of domains occupied by the star : " << nzet << endl ; 
    
    ost << "Equation of state : " << endl ; 
    ost << eos << endl ; 
    
    ost << endl << "Central enthalpy : " << ent.val_grid_point(0,0,0,0) << " c^2" << endl ; 
    ost << "Central proper baryon density : " << nbar.val_grid_point(0,0,0,0) 
	<< " x 0.1 fm^-3" << endl ; 
    ost << "Central proper energy density : " << ener.val_grid_point(0,0,0,0) 
	<< " rho_nuc c^2" << endl ; 
    ost << "Central pressure : " << press.val_grid_point(0,0,0,0) 
	<< " rho_nuc c^2" << endl ; 
    
    ost << endl ;
    ost << "Central lapse N :      " << nnn.val_grid_point(0,0,0,0) <<  endl ; 
    ost << "Central value of beta : " << beta.val_grid_point(0,0,0,0) <<  endl ; 
  
    ost << endl 
	<< "Coordinate equatorial radius (phi=0) a1 =    " 
	<< ray_eq()/km << " km" << endl ;  
    ost << "Coordinate equatorial radius (phi=pi/2) a2 = " 
	<< ray_eq_pis2()/km << " km" << endl ;  
    ost << "Coordinate equatorial radius (phi=pi):       " 
	<< ray_eq_pi()/km << " km" << endl ;  
    ost << "Coordinate polar radius a3 =                 " 
	<< ray_pole()/km << " km" << endl ;  
    ost << "Axis ratio a2/a1 = " << ray_eq_pis2() / ray_eq() 
	<< "  a3/a1 = " << ray_pole() / ray_eq() << endl ; 	
    ost << endl << "Baryon mass :        " << mass_b() / msol << " M_sol" << endl ; 
    ost << "Gravitational mass : " << mass_g() / msol << " M_sol" << endl ; 
    

    return ost ; 
}

		//-----------------------------------------//
		//	Computation of hydro quantities	   //
		//-----------------------------------------//

void Star::equation_of_state() {

	Scalar ent_eos = ent ;


    // Slight rescale of the enthalpy field in case of 2 domains inside the
    //  star


        double epsilon = 1.e-12 ;

	const Mg3d* mg = mp.get_mg() ;
        int nz = mg->get_nzone() ;

        Mtbl xi(mg) ;
        xi.set_etat_qcq() ;
        for (int l=0; l<nz; l++) {
        	xi.t[l]->set_etat_qcq() ;
        	for (int k=0; k<mg->get_np(l); k++) {
        		for (int j=0; j<mg->get_nt(l); j++) {
        			for (int i=0; i<mg->get_nr(l); i++) {
        				xi.set(l,k,j,i) =
        					mg->get_grille3d(l)->x[i] ;
        			}
        		}
        	}

        }

     	Scalar fact_ent(mp) ;
     	fact_ent.allocate_all() ;
     	
     	fact_ent.set_domain(0) = 1 + epsilon * xi(0) * xi(0) ;
     	fact_ent.set_domain(1) = 1 - 0.25 * epsilon * (xi(1) - 1) * (xi(1) - 1) ;
     	
     	for (int l=nzet; l<nz; l++) {
     		fact_ent.set_domain(l) = 1 ;
     	}

    if (nzet > 1) {

    	if (nzet > 2) {
    	
    		cout << "Star::equation_of_state: not ready yet for nzet > 2 !"
    		     << endl ;    	
    	}

    	ent_eos = fact_ent * ent_eos ;
    	ent_eos.std_spectral_base() ;
    }





    // Call to the EOS (the EOS is called domain by domain in order to
    //          allow for the use of MEos)

    Scalar tempo(mp) ;

    nbar.set_etat_qcq() ;
    nbar = 0 ;
    for (int l=0; l<nzet; l++) {

        Param par ;       // Paramater for multi-domain equation of state
        par.add_int(l) ;

        tempo =  eos.nbar_ent(ent_eos, 1, l, &par) ;

        nbar = nbar + tempo ;

    }

    ener.set_etat_qcq() ;
    ener = 0 ;
    for (int l=0; l<nzet; l++) {

        Param par ;    // Paramater for multi-domain equation of state
        par.add_int(l) ;

        tempo =  eos.ener_ent(ent_eos, 1, l, &par) ;

        ener = ener + tempo ;

    }

    press.set_etat_qcq() ;
    press = 0 ;
    for (int l=0; l<nzet; l++) {

        Param par ;     // Paramater for multi-domain equation of state
        par.add_int(l) ;

        tempo =  eos.press_ent(ent_eos, 1, l, &par) ;

        press = press + tempo ;

    }


    // Set the bases for spectral expansion
    nbar.std_spectral_base() ; 
    ener.std_spectral_base() ; 
    press.std_spectral_base() ; 

    // The derived quantities are obsolete
    del_deriv() ; 
    
}

void Star::hydro_euler() {
    
    cout << 
    "Star::hydro_euler : hydro_euler must be called via a derived class"
    << endl << " of Star !" << endl ; 
    
    abort() ;        
    
}
