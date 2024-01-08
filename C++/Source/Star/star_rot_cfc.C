/*
 *  Methods of class Star_rot_CFC
 *
 *    (see file star_rot_dirac.h for documentation).
 *
 */

/*
 *   Copyright (c) 2024  Jerome Novak
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

 


// Lorene headers
#include "star_rot_cfc.h"
#include "unites.h" 
#include "utilitaires.h"


                   //--------------//
                   // Constructors //
                   //--------------//

// Standard constructor
//-------------------------
namespace Lorene {
  Star_rot_CFC::Star_rot_CFC(Map& mpi, int nzet_i, const Eos& eos_i, int relat_i, int filter)
                   : Star(mpi, nzet_i, eos_i),
		     relat_type(relat_i),
		     spectral_filter(filter),
		     psi(mpi),
		     j_euler(mpi, CON, mpi.get_bvect_spher()),
		     v2(mpi),
		     psi4(mpi),
		     psi2(mpi),
		     flat(mpi.flat_met_spher()),
		     hatA(mpi, CON, mpi.get_bvect_spher()),
		     hatA_quad(mpi)
{
  assert (relat_type == 3) ;
  assert ((spectral_filter>=0) && (spectral_filter<100)) ;

  // Initialization to a static state
  omega = 0 ;
  v2 = 0 ;

  // All the matter quantities are initialized to zero
  j_euler.set_etat_zero() ;

  // Initialization to a flat case
  psi = 1 ;
  psi.std_spectral_base( );
  psi4 = 1 ;
  psi2 = 1 ;
  hatA.set_etat_zero() ;
  hatA_quad = 0 ;

  // Pointers of derived quantities initialized to zero : 
  set_der_0x0() ;

} 
		     

// Copy constructor
//-----------------
Star_rot_CFC::Star_rot_CFC(const Star_rot_CFC& star)
                   : Star(star),
		     relat_type(star.relat_type),
		     spectral_filter(star.spectral_filter),
		     psi(star.psi),
		     j_euler(star.j_euler),
		     v2(star.v2),
		     psi4(star.psi4),
		     psi2(star.psi2),
		     flat(star.flat),
		     hatA(star.hatA),
		     hatA_quad(star.hatA_quad)
{

  omega = star.omega ;

  // Pointers of derived quantities initialized to zero : 
  set_der_0x0() ;

}


//Constructor from a file 
//------------------------
Star_rot_CFC::Star_rot_CFC(Map& mpi, const Eos& eos_i, FILE* fich)
                  : Star(mpi, eos_i, fich),
		    psi(mpi, *(mpi.get_mg()), fich),
		    j_euler(mpi, CON, mpi.get_bvect_spher()),
		    v2(mpi),
		    psi4(mpi),
		    psi2(mpi),
		    flat(mpi.flat_met_spher()),
		    hatA(mpi, CON, mpi.get_bvect_spher()),
		    hatA_quad(mpi)
{

  // Pointers of derived quantities initialized to zero 
  //----------------------------------------------------
  set_der_0x0() ;

  fread_be(&relat_type, sizeof(int), 1, fich) ;
  fread_be(&spectral_filter, sizeof(int), 1, fich) ;

  // Metric fields are read in the file:
  fread_be(&omega, sizeof(double), 1, fich) ;
  Vector shift_tmp(mpi, mpi.get_bvect_spher(), fich) ;
  beta = shift_tmp ;

  update_metric() ;

  equation_of_state() ;

  hydro_euler() ;

  // Computation of \hat{A}^{ij} and its norm
  //------------------------------------------
  Vector sou_Xi = 2*Unites::qpig*psi4*psi4*psi2*j_euler ;
  double lambda = 1./3. ;
  Vector Xi_new = sou_Xi.poisson(lambda) ;
  
  hatA = Xi_new.ope_killing_conf(flat) ;
  hatA_quad = contract(hatA, 0, 1, hatA.up_down(flat), 0, 1) ;
    
}


                      //------------// 
                      // Destructor //
                      //------------//

Star_rot_CFC::~Star_rot_CFC(){

  Star_rot_CFC::del_deriv() ;
  
}


               //----------------------------------//
               // Management of derived quantities //
               //----------------------------------//

void Star_rot_CFC::del_deriv() const {

       if (p_angu_mom != 0x0) delete p_angu_mom ;
       if (p_grv2 != 0x0) delete p_grv2 ;
       if (p_grv3 != 0x0) delete p_grv3 ;
       if (p_tsw != 0x0) delete p_tsw ;
       if (p_r_circ != 0x0) delete p_r_circ ;
       if (p_rp_circ != 0x0) delete p_rp_circ ;

       set_der_0x0() ;

       Star::del_deriv() ;

}


void Star_rot_CFC::set_der_0x0() const {

       p_angu_mom = 0x0 ;
       p_grv2 = 0x0 ;
       p_grv3 = 0x0 ;
       p_tsw = 0x0 ;
       p_r_circ = 0x0 ;
       p_rp_circ = 0x0 ;

}


void Star_rot_CFC::del_hydro_euler() {

    j_euler.set_etat_nondef() ;
    v2.set_etat_nondef() ;

    del_deriv() ;
    
    Star::del_hydro_euler() ;

}



                    //---------------//
                    //  Assignment   //
                    //---------------//   

// Assignment to another Star_rot_CFC
// ------------------------------------

void Star_rot_CFC::operator=(const Star_rot_CFC& star) {

     // Assignment of quantities common to all the derived classes of Star
     Star::operator=(star) ;

     // Assignment of proper quantities of class Star_rot_CFC
     relat_type = star.relat_type ;
     spectral_filter = star.spectral_filter ;
     psi = star.psi ;
     omega = star.omega ;
     j_euler = star.j_euler ;
     v2 = star.v2 ;
     psi4 = star.psi4 ;
     psi2 = star.psi2 ;
     hatA = star.hatA ;
     hatA_quad = star.hatA_quad ;

     assert(&flat == &star.flat) ;

     del_deriv() ;    // Deletes all derived quantities

}
     

                      //-----------//
                      //  Outputs  //
                      //-----------//

// Save in a file
// --------------

void Star_rot_CFC::sauve(FILE* fich) const {

      Star::sauve(fich) ;

      psi.sauve(fich) ;
      
      fwrite_be(&relat_type, sizeof(int), 1, fich) ;
      fwrite_be(&spectral_filter, sizeof(int), 1, fich) ;
      fwrite_be(&omega, sizeof(double), 1, fich) ;
      beta.sauve(fich) ;

}


// Printing
// ---------

ostream& Star_rot_CFC::operator>>(ostream& ost) const {

  using namespace Unites ;

     Star::operator>>(ost) ;

     ost << "Rotating star " ;
     switch (relat_type) {
     case 0 : {
       ost << "in Newtonian theory" << endl ;
       break;
     }
     case 1 : {
       ost << "in Dirac gauge" << endl ;
       break ;
     }
     case 2 : {
       ost << "in Conformal Flatness Condition (CFC)" << endl ;
       break ;
     }
     case 3: {
       ost << "in Conformal Flatness Condition (CFC), with hat{A}^{ij}_{TT} = 0\n"
	   << "(see Cordero-Carrion et al. (2009) for details" << endl ;
       break ;
     }
     default: {
       cerr << "Star_rot_CFC::operator>> : unknown value for 'relat_type'" << endl ;
       cerr << "relat_type = " << relat_type << endl ;
       abort() ;
     }
     }

     // Only uniformly rotating star for the moment....
     ost << endl ;
     ost << "Uniformly rotating star" << endl ;
     ost << "-----------------------" << endl ;
     if (spectral_filter > 0)
	 ost << "hydro sources of equations are filtered\n"
	     << "with " << spectral_filter << "-order exponential filter" << endl ;

     double freq = omega/ (2.*M_PI) ;
     ost << "Omega : " << omega * f_unit
         << " rad/s    f : " << freq * f_unit << " Hz" << endl ;
     ost << "Rotation period : " << 1000. / (freq * f_unit) << " ms"
         << endl ;

     ost << "Error on the virial identity GRV2 : " << endl ;
     ost << "GRV2 = " << grv2() << endl ;
     ost << "Error on the virial identity GRV3 : " << endl ;
     ost << "GRV3 = " << grv3() << endl ;

     ost << "Angular momentum J :    "
         << angu_mom()/( qpig / (4*M_PI) *msol*msol) << " G M_sol^2 / c"
         << endl ;
     ost << "c J / (G M^2) :         "
         << angu_mom()/( qpig / (4*M_PI) * pow(mass_g(), 2.) ) << endl ;

     if (omega != 0.) {
       double mom_iner = angu_mom() / omega ; 
       double mom_iner_38si = mom_iner * rho_unit * (pow(r_unit, double(5.)) 
         / double(1.e38) ) ; 
       ost << "Moment of inertia:       " << mom_iner_38si << " 10^38 kg m^2"
	   << endl ; 
     }

     ost << "Ratio T/W :              " << tsw() << endl ;
     ost << "Circumferential equatorial radius R_circ :     "
      	 << r_circ()/km << " km" << endl ;
     if (mp.get_mg()->get_np(0) == 1) 
       ost << "Circumferential polar radius Rp_circ :     "
	   << rp_circ()/km << " km" << endl ;
     ost << "Coordinate equatorial radius r_eq : " << ray_eq()/km << " km"
      	 << endl ;
     ost << "Flattening r_pole/r_eq :  " << aplat() << endl ;
     if (mp.get_mg()->get_np(0) == 1) 
       ost << "Ellipticity sqrt(1-(Rp_circ/R_circ)^2) :  " << ellipt() << endl ;

     double compact = qpig/(4.*M_PI) * mass_g() / r_circ() ;
     ost << "Compaction parameter M_g / R_circ : " << compact << endl ; 
     

     // More to come here.....

     return ost ;

}
}
