/*
 * Methods for class Et_rot_diff.
 *
 * (see file et_rot_diff.h for documentation).
 *
 */

/*
 *   Copyright (c) 2001 Eric Gourgoulhon
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


char et_rot_diff_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:28  e_gourgoulhon
 * Initial revision
 *
 * Revision 1.3  2001/10/25  09:20:54  eric
 * Ajout de la fonction virtuelle display_poly.
 * Affichage de Max nbar, ener et press.
 *
 * Revision 1.2  2001/10/24  16:23:01  eric
 * *** empty log message ***
 *
 * Revision 1.1  2001/10/19  08:18:10  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */

// Headers C
#include "math.h"

// Headers Lorene
#include "et_rot_diff.h"
#include "eos.h"
#include "nbr_spx.h"

			    //--------------//
			    // Constructors //
			    //--------------//

// Standard constructor
// --------------------
Et_rot_diff::Et_rot_diff(Map& mp_i, int nzet_i, bool relat, const Eos& eos_i, 
			 double (*frot_i)(double, const Tbl&), 
			 double (*primfrot_i)(double, const Tbl&), 
			 const Tbl& par_frot_i)
		       : Etoile_rot(mp_i, nzet_i, relat, eos_i),
			 frot(frot_i), 
			 primfrot(primfrot_i), 
			 par_frot(par_frot_i), 
			 omega_field(mp_i), 
			 prim_field(mp_i)
			 {

    // To make sure that omega is not used
    omega = __infinity ; 
    
    // Initialization to a static state : 
    omega_field = 0 ; 
    prim_field = 0 ; 
    omega_min = 0 ;
    omega_max = 0 ; 
    
} 

// Copy constructor
// ----------------
Et_rot_diff::Et_rot_diff(const Et_rot_diff& et)
		       : Etoile_rot(et), 
		         frot(et.frot), 
		         primfrot(et.primfrot), 
			 par_frot(et.par_frot), 
			 omega_field(et.omega_field), 
			 omega_min(et.omega_min), 
			 omega_max(et.omega_max),  
			 prim_field(et.prim_field)
			 {}
			 

// Constructor from a file
// -----------------------
Et_rot_diff::Et_rot_diff(Map& mp_i, const Eos& eos_i, FILE* fich,
			 double (*frot_i)(double, const Tbl&), 
			 double (*primfrot_i)(double, const Tbl&) )
			: Etoile_rot(mp_i, eos_i, fich),
			  frot(frot_i), 
			  primfrot(primfrot_i), 
			  par_frot(fich), 
			  omega_field(mp_i), 
			  prim_field(mp_i)
			 {

    Tenseur omega_field_file(mp, fich) ; 
    omega_field = omega_field_file ; 
    fait_prim_field() ; 
			 
    // omega_min and omega_max are read in the file:     
    fread(&omega_min, sizeof(double), 1, fich) ;		
    fread(&omega_max, sizeof(double), 1, fich) ;		

}


			    //------------//
			    // Destructor //
			    //------------//

Et_rot_diff::~Et_rot_diff(){} 


			    //--------------//
			    //  Assignment  //
			    //--------------//

// Assignment to another Et_rot_diff
// ---------------------------------
void Et_rot_diff::operator=(const Et_rot_diff& et) {

    // Assignment of quantities common to all the derived classes of Etoile_rot
    Etoile_rot::operator=(et) ;	    
    
    // Assignment of proper quantities of class Etoile_rot
    frot = et.frot ; 
    primfrot = et.primfrot ; 
    par_frot = et.par_frot ; 
    omega_field = et.omega_field ; 
    prim_field = et.prim_field ; 
    omega_min = et.omega_min ; 
    omega_max = et.omega_max ; 

}

			    //--------------//
			    //	  Outputs   //
			    //--------------//

// Save in a file
// --------------

void Et_rot_diff::sauve(FILE* fich) const {
    
    Etoile_rot::sauve(fich) ; 
    
    par_frot.sauve(fich) ; 
    
    omega_field.sauve(fich) ; 
    
    fwrite(&omega_min, sizeof(double), 1, fich) ;		
    fwrite(&omega_max, sizeof(double), 1, fich) ;		

}


// Printing
// --------

ostream& Et_rot_diff::operator>>(ostream& ost) const {
    
    #include "unites.h"	    
    // To avoid some compilation warnings
    if (&ost == 0x0) {
	cout << qpig << msol << mevpfm3 << km << endl ; 
    }    

    Etoile_rot::operator>>(ost) ; 
    
    ost << endl ; 
    ost << "Differentially rotating  star" << endl ; 
    ost << "-----------------------------" << endl ; 
    
    ost << endl << "Parameters of F(Omega) : " << endl ; 
    ost << par_frot << endl ; 
    
    ost << "Min, Max of Omega/(2pi) : " << omega_min / (2*M_PI) * f_unit 
	    << " Hz,  " << omega_max / (2*M_PI) * f_unit << " Hz" << endl ; 
    int lsurf = nzet - 1; 
    int nt = mp.get_mg()->get_nt(lsurf) ; 
    int nr = mp.get_mg()->get_nr(lsurf) ; 
    ost << "Central, equatorial value of Omega:        "
	<< omega_field()(0, 0, 0, 0) * f_unit << " rad/s,   " 
	<< omega_field()(nzet-1, 0, nt-1, nr-1) * f_unit << " rad/s" << endl ; 
    
    ost << "Central, equatorial value of Omega/(2 Pi): "
	<< omega_field()(0, 0, 0, 0) * f_unit / (2*M_PI) << " Hz,      " 
	<< omega_field()(nzet-1, 0, nt-1, nr-1) * f_unit / (2*M_PI) 
	<< " Hz" << endl ; 
    
    double nbar_max = max( max( nbar() ) ) ; 
    double ener_max = max( max( ener() ) ) ; 
    double press_max = max( max( press() ) ) ; 
    ost << "Max prop. bar. dens. :          " << nbar_max 
 	<< " x 0.1 fm^-3 = " << nbar_max / nbar()(0, 0, 0, 0) << " central"
	<< endl ; 
    ost << "Max prop. ener. dens. (e_max) : " << ener_max
	<< " rho_nuc c^2 = " << ener_max / ener()(0, 0, 0, 0) << " central"
	<< endl ; 
    ost << "Max pressure :                  " << press_max
	<< " rho_nuc c^2 = " << press_max / press()(0, 0, 0, 0) << " central"
	<< endl ; 
   
    // Length scale set by the maximum energy density:
    double len_rho = 1. / sqrt( ggrav * ener_max ) ;
    ost << endl << "Value of A = par_frot(1) in units of e_max^{1/2} : " << 
	    par_frot(1) / len_rho << endl ; 
    
    ost << "Value of A / r_eq : " << 
	    par_frot(1) / ray_eq() << endl ; 
    
    ost << "r_p/r_eq : " << aplat() << endl ; 
    ost << "KEH l^2 = (c/G^2)^{2/3} J^2 e_max^{1/3} M_B^{-10/3} : " <<
	angu_mom() * angu_mom() / pow(len_rho, 0.6666666666666666)
	/ pow(mass_b(), 3.3333333333333333) 
	/ pow(ggrav, 1.3333333333333333) << endl ;
	 
    ost << "M e_max^{1/2} : " << ggrav * mass_g() / len_rho << endl ; 
    
    ost << "r_eq e_max^{1/2} : " << ray_eq() / len_rho << endl ; 
    
    ost << "T/W : " << tsw() << endl ; 
    
    ost << "Omega_c / e_max^{1/2} : " << get_omega_c() * len_rho << endl ; 
     
    display_poly(ost) ; 

    return ost ;
    
    
}

// display_poly
// ------------

void Et_rot_diff::display_poly(ostream& ost) const {

    #include "unites.h"	    
    // To avoid some compilation warnings
    if (&ost == 0x0) {
	cout << f_unit << qpig << msol << km << mevpfm3 << endl ; 
    }    

    Etoile_rot::display_poly( ost ) ; 
    
    const Eos_poly* p_eos_poly = dynamic_cast<const Eos_poly*>( &eos ) ; 	  

    if (p_eos_poly != 0x0) {

	double kappa = p_eos_poly->get_kap() ; 
	double gamma = p_eos_poly->get_gam() ;  ; 

	// kappa^{n/2}
	double kap_ns2 = pow( kappa,  0.5 /(gamma-1) ) ; 
    
	// Polytropic unit of length in terms of r_unit : 
	double r_poly = kap_ns2 / sqrt(ggrav) ; 
    
	// Polytropic unit of time in terms of t_unit :
	double t_poly = r_poly ; 

	// Polytropic unit of density in terms of rho_unit :
	double rho_poly = 1. / (ggrav * r_poly * r_poly) ;  

	ost.precision(10) ; 
	ost << "  n_max     : " << max( max( nbar() ) ) / rho_poly << endl ; 
	ost << "  e_max     : " << max( max( ener() ) ) / rho_poly << endl ; 
	ost << "  P_min     : " << 2.*M_PI / omega_max / t_poly << endl ; 
	ost << "  P_max     : " << 2.*M_PI / omega_min / t_poly << endl ; 
    
	int lsurf = nzet - 1; 
	int nt = mp.get_mg()->get_nt(lsurf) ; 
	int nr = mp.get_mg()->get_nr(lsurf) ; 
	ost << "  P_eq      : " << 2.*M_PI / 
		omega_field()(nzet-1, 0, nt-1, nr-1) / t_poly << endl ; 
	    
    }

}



		//-----------------------//
		// 	Miscellaneous	 //
		//-----------------------//

double Et_rot_diff::funct_omega(double omeg) const {

 	return frot(omeg, par_frot) ;
 	
}

double Et_rot_diff::prim_funct_omega(double omeg) const {

 	return primfrot(omeg, par_frot) ;
 	
}

double Et_rot_diff::get_omega_c() const {
    
    return omega_field()(0, 0, 0, 0) ; 
    
}


		
