/*
 *  Methods of class Bin_bhns_extr
 *
 *    (see file bin_bhns_extr.h for documentation).
 *
 */

/*
 *   Copyright (c) 2004-2005 Keisuke Taniguchi
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
 * $Id$
 * $Log$
 * Revision 1.11  2016/12/05 16:17:46  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.10  2014/10/13 08:52:41  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.9  2014/10/06 15:13:00  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.8  2005/02/28 23:06:16  k_taniguchi
 * Modification of some arguments to include the case of the conformally
 * flat background metric.
 *
 * Revision 1.7  2004/12/13 21:08:59  k_taniguchi
 * Addition of some outputs in display_poly.
 *
 * Revision 1.6  2004/12/03 20:17:52  k_taniguchi
 * Correction of "display_poly">
 *
 * Revision 1.5  2004/12/02 22:43:35  k_taniguchi
 * Modification of "display_poly".
 *
 * Revision 1.4  2004/12/02 17:36:42  k_taniguchi
 * Modification of "display_poly()".
 *
 * Revision 1.3  2004/12/01 16:37:11  k_taniguchi
 * Modification of the output again.
 *
 * Revision 1.2  2004/12/01 16:27:09  k_taniguchi
 * Modification of the output.
 *
 * Revision 1.1  2004/11/30 20:45:49  k_taniguchi
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// C headers
#include <cmath>

// Lorene headers
#include "bin_bhns_extr.h"
#include "eos.h"
#include "utilitaires.h"
#include "unites.h"

          //--------------------------------//
          //          Constructors          //
          //--------------------------------//

// Standard constructor
// --------------------
namespace Lorene {
Bin_bhns_extr::Bin_bhns_extr(Map& mp, int nzet, const Eos& eos, bool irrot,
			     bool relat, bool kerrs, bool multi)
    : ref_triad(0., "Absolute frame Cartesian basis"),
      star(mp, nzet, relat, eos, irrot, ref_triad, kerrs, multi)
{

    omega = 0. ;
    separ = 0. ;
    mass_bh = 0. ;

    // Pointers of derived quantities initialized to zero :
    set_der_0x0() ;

}

// Copy constructor
// ----------------
Bin_bhns_extr::Bin_bhns_extr(const Bin_bhns_extr& bibi)
    : ref_triad(0., "Absolute frame Cartesian basis"),
      star(bibi.star),
      omega(bibi.omega),
      separ(bibi.separ),
      mass_bh(bibi.mass_bh)
{

    // Pointers of derived quantities initialized to zero :
    set_der_0x0() ;

}

// Constructor from a file
// -----------------------
Bin_bhns_extr::Bin_bhns_extr(Map& mp, const Eos& eos, FILE* fich)
    : ref_triad(0., "Absolute frame Cartesian basis"),
      star(mp, eos, ref_triad, fich)
{

    // omega and x_axe are read in the file :
    fread_be(&omega, sizeof(double), 1, fich) ;
    fread_be(&separ, sizeof(double), 1, fich) ;
    fread_be(&mass_bh, sizeof(double), 1, fich) ;

    // Pointers of derived quantities initialized to zero :
    set_der_0x0() ;

}

          //------------------------------//
          //          Destructor          //
          //------------------------------//

Bin_bhns_extr::~Bin_bhns_extr()
{
    del_deriv() ;
}

          //----------------------------------------------------//
          //          Management of derived quantities          //
          //----------------------------------------------------//

void Bin_bhns_extr::del_deriv() const {

    if (p_xa_barycenter_extr != 0x0) delete p_xa_barycenter_extr ;
    if (p_ya_barycenter_extr != 0x0) delete p_ya_barycenter_extr ;
    if (p_mass_b_extr != 0x0) delete p_mass_b_extr ;

    set_der_0x0() ;

}

void Bin_bhns_extr::set_der_0x0() const {

    p_xa_barycenter_extr = 0x0 ;
    p_ya_barycenter_extr = 0x0 ;
    p_mass_b_extr = 0x0 ;

}

          //------------------------------//
          //          Assignment          //
          //------------------------------//

// Assignment to another Bin_bhns_extr
// -----------------------------------
void Bin_bhns_extr::operator=(const Bin_bhns_extr& bibi) {

    assert( bibi.ref_triad == ref_triad ) ;

    star = bibi.star ;

    omega = bibi.omega ;
    separ = bibi.separ ;
    mass_bh = bibi.mass_bh ;

    // ref_triad remains unchanged

    del_deriv() ;  // Deletes all derived quantities

}

          //---------------------------//
          //          Outputs          //
          //---------------------------//

// Save in a file
// --------------
void Bin_bhns_extr::sauve(FILE* fich) const {

    star.sauve(fich) ;

    fwrite_be(&omega, sizeof(double), 1, fich) ;
    fwrite_be(&separ, sizeof(double), 1, fich) ;
    fwrite_be(&mass_bh, sizeof(double), 1, fich) ;

}

// Printing
// --------
ostream& operator<<(ostream& ost, const Bin_bhns_extr& bibi) {

    bibi >> ost ;
    return ost ;

}

ostream& Bin_bhns_extr::operator>>(ostream& ost) const {

  using namespace Unites ;

    ost << endl ;
    ost << "Binary BH-NS system" << endl ;
    ost << "===================" << endl ;
    ost << endl ;
    if (star.in_kerrschild()) {
        ost << "Kerr-Schild background metric" << endl ;
	ost << "-----------------------------" << endl ;
    }
    else {
        ost << "Conformally flat background metric" << endl ;
	ost << "----------------------------------" << endl ;
    }
    if (star.with_multipole()) {
        ost << "Multipole falloff boundary condition" << endl ;
	ost << "------------------------------------" << endl ;
    }
    else {
        ost << "1/r falloff boundary condition" << endl ;
	ost << "------------------------------" << endl ;
    }
    ost << endl
	<< "Orbital angular velocity :                "
	<< omega * f_unit << " rad/s" << endl ;
    ost << "Coordinate separation between BH and NS : "
	<< separ / km << " km" << endl ;

    ost << endl << "Neutron star : " ;
    ost << endl << "============   " << endl ;
    ost << endl ;
    if (star.is_relativistic()) {
	ost << "Relativistic star" << endl ;
	ost << "-----------------" << endl ;
    }
    else {
        ost << "WARNING : BH-NS binary should be relativistic !!!" << endl ;
	ost << "-------------------------------------------------" << endl ;
    }
    ost << "Number of domains occupied by the star : " << star.get_nzet()
	<< endl ;
    ost << "Equation of state : " << endl ;
    ost << star.get_eos() << endl ;

    ost << endl
	<< "Enthalpy at the coordinate origin :              "
	<< star.get_ent()()(0,0,0,0) << " c^2" << endl ;
    ost << "Proper baryon density at the coordinate origin : "
	<< star.get_nbar()()(0,0,0,0) << " x 0.1 fm^-3" << endl ;
    ost << "Proper energy density at the coordinate origin : "
	<< star.get_ener()()(0,0,0,0) << " rho_nuc c^2" << endl ;
    ost << "Pressure at the coordinate origin :              "
	<< star.get_press()()(0,0,0,0) << " rho_nuc c^2" << endl ;
    ost << endl
	<< "Lapse N at the coordinate origin :              "
	<< star.get_nnn()()(0,0,0,0) << endl ;
    ost << "Conformal factor A^2 at the coordinate origin : "
	<< star.get_a_car()()(0,0,0,0) << endl ;

    ost << endl
	<< "Equatorial radius (to BH) a_to :       "
	<< star.ray_eq_pi()/km << " km" << endl ;
    ost << "Equatorial radius (opp. to BH) a_opp : "
	<< star.ray_eq()/km << " km" << endl ;

    ost << endl
	<< "Baryon mass in isolation :        " << star.mass_b() / msol
	<< " Mo" << endl
	<< "Gravitational mass in isolation : " << star.mass_g() / msol
	<< " Mo" << endl
	<< "Baryon mass in a binary system :  " << mass_b_extr() / msol
	<< " Mo" << endl ;

    ost << endl ;
    ost << "Star in a binary system" << endl ; 
    ost << "-----------------------" << endl ;

    if (star.is_irrotational()) {
        ost << "irrotational configuration" << endl ;
    }
    else {
        ost << "corotating configuration" << endl ;
    }

    ost << "Absolute abscidia of the stellar center: "
	<< star.get_mp().get_ori_x()/km << " km" << endl ;
    ost << "Orientation with respect to the absolute frame : "
	<< star.get_mp().get_rot_phi() << " rad" << endl ;

    double r0 = 0.5 * ( star.ray_eq() + star.ray_eq_pi() ) ;
                           // = 1 in Baumgarte et al.
    double d_ns = separ + star.ray_eq() - r0 ;
                           // Orbital separation of Baumgarte et al.
    double d_tilde = d_ns / r0 ; // Normalized separation of Baumgarte et al.

    ost << endl << "Comparison with those by Baumgarte et al. :" << endl ;
    ost << "     Radius r0 :              " << r0 / km << " km " << endl ;
    ost << "     Separation d :           " << d_ns / km << " km " << endl ;
    ost << "     Normalized sep. (d/r0) : " << d_tilde << endl ;

    ost << endl << "Black hole : " ;
    ost << endl << "==========   " << endl ;
    ost << "Gravitational mass of BH : "
	<< mass_bh / msol << " M_sol" << endl ;

    return ost ;

}

// Display in polytropic units
// ---------------------------
void Bin_bhns_extr::display_poly(ostream& ost) const {

  using namespace Unites ;

    const Eos* p_eos = &( star.get_eos() ) ;
    const Eos_poly* p_eos_poly = dynamic_cast<const Eos_poly*>( p_eos ) ;

    if (p_eos_poly != 0x0) {

        double kappa = p_eos_poly->get_kap() ;
	double gamma = p_eos_poly->get_gam() ;
	double kap_ns2 = pow( kappa,  0.5 /(gamma-1) ) ;

	// Polytropic unit of length in terms of r_unit :
	double r_poly = kap_ns2 / sqrt(ggrav) ;

	// Polytropic unit of time in terms of t_unit :
	double t_poly = r_poly ;

	// Polytropic unit of mass in terms of m_unit :
	double m_poly = r_poly / ggrav ;

	// Polytropic unit of angular momentum in terms of j_unit :
	//	  double j_poly = r_poly * r_poly / ggrav ;

	double r0 = 0.5 * ( star.ray_eq() + star.ray_eq_pi() ) ;
	            // = 1 in Baumgarte et al.
	double d_ns = separ + star.ray_eq() - r0 ;
	            // Orbital separation of Baumgarte et al.

	ost.precision(16) ;
	ost << endl << "Quantities in polytropic units : " ;
	ost << endl << "==============================" << endl ;
	ost << " ( r_poly = " << r_poly / km << " km )" << endl ;
	ost << "  d_e_max   :   " << separ / r_poly << endl ;
	ost << "  d_G_x :       " << xa_barycenter_extr() / r_poly << endl
	    << "  d_G_y :       " << ya_barycenter_extr() / r_poly << endl ;
	ost << "  d_bh/M_bh :   " << d_ns/r_poly / (mass_bh/m_poly) << endl ;
	ost << "  Omega :       " << omega * t_poly << endl ;
	ost << "  Omega M_bh :  " << omega * t_poly * mass_bh / m_poly
	    << endl ;
	ost << "  M_bar(NS) :   " << mass_b_extr() / m_poly << endl ;
	ost << "  M_bar(NS_0) : " << star.mass_b() / m_poly << endl ;
	ost << "  R_0(NS) :     "
	    << 0.5 * (star.ray_eq() + star.ray_eq_pi()) / r_poly << endl ;
	ost << "  M_grv(BH) :   " << mass_bh / m_poly << endl ;

    }

}
}
