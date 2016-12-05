/*
 *  Methods of class Star_rot_Dirac_diff
 *
 *    (see file star_rot_dirac_diff.h for documentation).
 *
 */

/*
 *   Copyright (c) 2005 Motoyuki Saijo
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

 

/*
 * $Header$
 *
 */


// C headers
#include <cmath>
#include <cassert>

// Lorene headers
#include "star_rot_dirac_diff.h"
#include "unites.h" 
#include "utilitaires.h"


                   //--------------//
                   // Constructors //
                   //--------------//

// Standard constructor
//-------------------------
namespace Lorene {
Star_rot_Dirac_diff::Star_rot_Dirac_diff(Map& mpi, int nzet_i, const Eos& eos_i,
			 double (*frot_i)(double, const Tbl&), 
			 double (*primfrot_i)(double, const Tbl&), 
			 const Tbl& par_frot_i)
                        :Star_rot_Dirac(mpi, nzet_i, eos_i),
			 frot(frot_i), 
			 primfrot(primfrot_i), 
			 par_frot(par_frot_i), 
			 omega_field(mpi), 
			 prim_field(mpi)
			 {

    // Initialization to a static state : 
    omega_field = 0 ; 
    prim_field = 0 ; 
    omega_min = 0 ;
    omega_max = 0 ; 
    
} 


// Copy constructor
//-----------------
Star_rot_Dirac_diff::Star_rot_Dirac_diff(const Star_rot_Dirac_diff& star)
                   : Star_rot_Dirac(star), 
		         frot(star.frot), 
		         primfrot(star.primfrot), 
			 par_frot(star.par_frot), 
			 omega_field(star.omega_field), 
			 omega_min(star.omega_min), 
			 omega_max(star.omega_max),  
			 prim_field(star.prim_field)
			 {}


//Constructor from a file //## to be more general...
//------------------------
Star_rot_Dirac_diff::Star_rot_Dirac_diff(Map& mpi, const Eos& eos_i, FILE* fich,
			 double (*frot_i)(double, const Tbl&), 
			 double (*primfrot_i)(double, const Tbl&) )
			: Star_rot_Dirac(mpi, eos_i, fich),
			  frot(frot_i), 
			  primfrot(primfrot_i), 
			  par_frot(fich), 
			  omega_field(mpi), 
			  prim_field(mpi)
			 {

    Scalar omega_field_file(mp, *(mp.get_mg()), fich) ; 
    omega_field = omega_field_file ; 
    fait_prim_field() ; 
			 
    // omega_min and omega_max are read in the file:     
    fread_be(&omega_min, sizeof(double), 1, fich) ;		
    fread_be(&omega_max, sizeof(double), 1, fich) ;		

}


                      //------------// 
                      // Destructor //
                      //------------//

Star_rot_Dirac_diff::~Star_rot_Dirac_diff(){}


                    //---------------//
                    //  Assignment   //
                    //---------------//   

// Assignment to another Star_rot_Dirac_diff
// ------------------------------------

void Star_rot_Dirac_diff::operator=(const Star_rot_Dirac_diff& star) {

     // Assignment of quantities common to all the derived classes of 
     // Star_rot_Dirac
     Star_rot_Dirac::operator=(star) ;

     // Assignment of proper quantities of class Star_rot_Dirac
    frot = star.frot ; 
    primfrot = star.primfrot ; 
    par_frot = star.par_frot ; 
    omega_field = star.omega_field ; 
    prim_field = star.prim_field ; 
    omega_min = star.omega_min ; 
    omega_max = star.omega_max ; 

}
     

                      //-----------//
                      //  Outputs  //
                      //-----------//

// Save in a file
// --------------

void Star_rot_Dirac_diff::sauve(FILE* fich) const {

      Star::sauve(fich) ;

      par_frot.sauve(fich) ; 
    
      omega_field.sauve(fich) ; 
    
      fwrite_be(&omega_min, sizeof(double), 1, fich) ;		
      fwrite_be(&omega_max, sizeof(double), 1, fich) ;		

      // What else to save? //## to be more general ...

}


// Printing
// ---------

ostream& Star_rot_Dirac_diff::operator>>(ostream& ost) const {

  using namespace Unites ;

     Star::operator>>(ost) ;

     ost << "Differentially rotating star in Dirac gauge" << '\n';

     // Only differentially rotating star for the moment....
     ost << '\n';
     ost << "Differentially rotating star" << '\n';
     ost << "-----------------------" << '\n';

    ost << '\n'<< "Parameters of F(Omega) : " << '\n'; 
    ost << par_frot << '\n'; 

    ost << "Min, Max of Omega/(2pi) : " << omega_min / (2*M_PI) * f_unit 
	    << " Hz,  " << omega_max / (2*M_PI) * f_unit << " Hz" << '\n'; 
    int lsurf = nzet - 1; 
    int nt = mp.get_mg()->get_nt(lsurf) ;
    int nr = mp.get_mg()->get_nr(lsurf) ;
    ost << "Central, equatorial value of Omega:        "
 	<< omega_field.val_grid_point(0, 0, 0, 0) * f_unit 
        << " rad/s,   " 
 	<< omega_field.val_grid_point(nzet-1, 0, nt-1, nr-1) * f_unit 
        << " rad/s" << '\n'; 
  
    ost << "Central, equatorial value of Omega/(2 Pi): "
 	<< omega_field.val_grid_point(0, 0, 0, 0) * f_unit / (2*M_PI) 
        << " Hz,      " 
 	<< omega_field.val_grid_point(nzet-1, 0, nt-1, nr-1) * 
             f_unit / (2*M_PI) 
 	<< " Hz" << '\n'; 

    ost << "Error on the virial identity GRV2 : " << '\n';
    ost << "GRV2 = " << grv2() << '\n';
    ost << "Error on the virial identity GRV3 : " << '\n';
    ost << "GRV3 = " << grv3() << '\n';
    
    ost << "Angular momentum J :    "
        << angu_mom()/( qpig / (4*M_PI) *msol*msol) << " G M_sol^2 / c"
        << '\n';
    ost << "c J / (G M^2) :         "
        << angu_mom()/( qpig / (4*M_PI) * pow(mass_g(), 2.) ) << '\n';

        if (omega != 0.) {
        double mom_iner = angu_mom() / omega ; 
        double mom_iner_38si = mom_iner * rho_unit * (pow(r_unit, double(5.)) 
          / double(1.e38) ) ; 
        ost << "Moment of inertia:       " << mom_iner_38si << " 10^38 kg m^2"
 	    << '\n'; 
        }

     ost << "Ratio T/W :              " << tsw() << '\n';

//     ost << "Omega_c / e_max^{1/2} : " << get_omega_c() * len_rho << '\n'; 

     ost << "Circumferential equatorial radius R_circ :     "
	 << r_circ()/km << " km" << '\n';
     ost << "Circumferential polar radius Rp_circ :     "
	 << rp_circ()/km << " km" << endl ;
     ost << "Coordinate equatorial radius r_eq : " << ray_eq()/km << " km"
	 << endl ;
     ost << "Flattening r_pole/r_eq :  " << aplat() << endl ;
     ost << "Ellipticity sqrt(1-(Rp_circ/R_circ)^2) :  " << ellipt() << endl ;
     double compact = qpig/(4.*M_PI) * mass_g() / r_circ() ;
     ost << "Compaction parameter M_g / R_circ : " << compact << '\n'; 
     

     // More to come here.....

     return ost ;

}


		//-----------------------//
		// 	Miscellaneous	 //
		//-----------------------//

double Star_rot_Dirac_diff::funct_omega(double omeg) const {

 	return frot(omeg, par_frot) ;
 	
}

double Star_rot_Dirac_diff::prim_funct_omega(double omeg) const {

 	return primfrot(omeg, par_frot) ;
 	
}

double Star_rot_Dirac_diff::get_omega_c() const {
    
    return omega_field.val_grid_point(0, 0, 0, 0) ; 
    
}
}
