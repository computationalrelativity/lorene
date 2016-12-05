/*
 *  Methods of class Bin_bhns
 *
 *    (see file bin_bhns.h for documentation).
 *
 */

/*
 *   Copyright (c) 2005-2007 Keisuke Taniguchi
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
 * Revision 1.5  2016/12/05 16:17:45  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.4  2014/10/13 08:52:41  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.3  2014/10/06 15:13:00  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2008/05/15 18:58:01  k_taniguchi
 * Change of some parameters and introduction of new
 * global quantities.
 *
 * Revision 1.1  2007/06/22 01:08:53  k_taniguchi
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// C++ headers
//#include <>

// C headers
#include <cmath>

// Lorene headers
#include "bin_bhns.h"
#include "eos.h"
#include "utilitaires.h"
#include "unites.h"

                    //---------------------//
                    //     Constructor     //
                    //---------------------//

// Standard constructor
// --------------------
namespace Lorene {
Bin_bhns::Bin_bhns(Map& mp_bh, Map& mp_ns, int nzet_i, const Eos& eos_i,
		   bool irrot_ns, bool kerrschild_i,
		   bool bc_lapconf_nd, bool bc_lapconf_fs, bool irrot_bh,
		   double massbh)
      : ref_triad(0., "Absolute frame Cartesian basis"),
	hole(mp_bh,kerrschild_i,bc_lapconf_nd,bc_lapconf_fs,irrot_bh, massbh),
	star(mp_ns, nzet_i, eos_i, irrot_ns)
{

    omega = 0. ;
    separ = 0. ;
    x_rot = 0. ;
    y_rot = 0. ;

    // Pointers of derived quantities are initialized to zero
    set_der_0x0() ;

}

// Copy constructor
// ----------------
Bin_bhns::Bin_bhns(const Bin_bhns& bibi)
      : ref_triad(0., "Absolute frame Cartesian basis"),
	hole(bibi.hole),
	star(bibi.star),
	omega(bibi.omega),
	separ(bibi.separ),
	x_rot(bibi.x_rot),
	y_rot(bibi.y_rot)
{

    // Pointers of derived quantities are initialized to zero
    set_der_0x0() ;

}

// Constructor from a file
// -----------------------
Bin_bhns::Bin_bhns(Map& mp_bh, Map& mp_ns, const Eos& eos_i, FILE* fich)
      : ref_triad(0., "Absolute frame Cartesian basis"),
	hole(mp_bh, fich),
	star(mp_ns, eos_i, fich)
{

    // omega, separ, and x_rot are read from the file
    fread_be(&omega, sizeof(double), 1, fich) ;
    fread_be(&separ, sizeof(double), 1, fich) ;
    fread_be(&x_rot, sizeof(double), 1, fich) ;
    fread_be(&y_rot, sizeof(double), 1, fich) ;

    // Pointers of derived quantities are initialized to zero
    set_der_0x0() ;

}


                    //--------------------//
                    //     Destructor     //
                    //--------------------//

Bin_bhns::~Bin_bhns()
{

    del_deriv() ;

}


                    //------------------------------------------//
                    //     Management of derived quantities     //
                    //------------------------------------------//

void Bin_bhns::del_deriv() const {

    if (p_mass_adm_bhns_surf != 0x0) delete p_mass_adm_bhns_surf ;
    if (p_mass_adm_bhns_vol != 0x0) delete p_mass_adm_bhns_vol ;
    if (p_mass_kom_bhns_surf != 0x0) delete p_mass_kom_bhns_surf ;
    if (p_mass_kom_bhns_vol != 0x0) delete p_mass_kom_bhns_vol ;
    if (p_line_mom_bhns != 0x0) delete p_line_mom_bhns ;
    if (p_angu_mom_bhns != 0x0) delete p_angu_mom_bhns ;
    if (p_virial_bhns_surf != 0x0) delete p_virial_bhns_surf ;
    if (p_virial_bhns_vol != 0x0) delete p_virial_bhns_vol ;
    if (p_xa_barycenter != 0x0) delete p_xa_barycenter ;
    if (p_ya_barycenter != 0x0) delete p_ya_barycenter ;
    if (p_omega_two_points != 0x0) delete p_omega_two_points ;
    //    if (p_ham_constr_bhns != 0x0) delete p_ham_constr_bhns ;
    //    if (p_mom_constr_bhns != 0x0) delete p_mom_constr_bhns ;

    set_der_0x0() ;

}

void Bin_bhns::set_der_0x0() const {

    p_mass_adm_bhns_surf = 0x0 ;
    p_mass_adm_bhns_vol = 0x0 ;
    p_mass_kom_bhns_surf = 0x0 ;
    p_mass_kom_bhns_vol = 0x0 ;
    p_line_mom_bhns = 0x0 ;
    p_angu_mom_bhns = 0x0 ;
    p_virial_bhns_surf = 0x0 ;
    p_virial_bhns_vol = 0x0 ;
    p_xa_barycenter = 0x0 ;
    p_ya_barycenter = 0x0 ;
    p_omega_two_points = 0x0 ;
    //    p_ham_constr_bhns = 0x0 ;
    //    p_mom_constr_bhns = 0x0 ;

}


                    //--------------------//
                    //     Assignment     //
                    //--------------------//

// Assignment to anothe Bin_bhns
// -----------------------------
void Bin_bhns::operator=(const Bin_bhns& bibi) {

    assert( bibi.ref_triad == ref_triad ) ;

    hole = bibi.hole ;
    star = bibi.star ;

    omega = bibi.omega ;
    separ = bibi.separ ;
    x_rot = bibi.x_rot ;
    y_rot = bibi.y_rot ;

    del_deriv() ;     // Deletes all derived quantities

}


                    //-----------------//
                    //     Outputs     //
                    //-----------------//

// Save in a file
// --------------
void Bin_bhns::sauve(FILE* fich) const {

    hole.sauve(fich) ;
    star.sauve(fich) ;

    fwrite_be(&omega, sizeof(double), 1, fich) ;
    fwrite_be(&separ, sizeof(double), 1, fich) ;
    fwrite_be(&x_rot, sizeof(double), 1, fich) ;
    fwrite_be(&y_rot, sizeof(double), 1, fich) ;

}

// Printing
// --------
ostream& operator<<(ostream& ost, const Bin_bhns& bibi) {

    bibi >> ost ;
    return ost ;

}

ostream& Bin_bhns::operator>>(ostream& ost) const {

    using namespace Unites ;

    ost << endl ;
    ost << "BH-NS binary system" << endl ;
    ost << "===================" << endl ;

    ost << endl
	<< "Orbital angular velocity :                "
	<< omega * f_unit << " [rad/s]" << endl ;
    ost << "Coordinate separation between BH and NS : "
	<< separ / km << " [km]" << endl ;
    ost << "Position of the rotation axis :           "
	<< x_rot / km << " [km]" << endl ;

    ost << endl << "Black hole : " ;
    ost << endl << "==========   " << endl ;

    int nt = (hole.get_mp()).get_mg()->get_nt(1) ;

    ost << "Irreducible mass of BH :         "
	<< hole.mass_irr_bhns() / msol << " [Mo]" << endl ;
    ost << "Mass in the background :         "
	<< hole.get_mass_bh() / msol << " [Mo]" << endl ;
    ost << "Radius of the apparent horison : "
	<< hole.rad_ah() / km << " [km]" << endl ;
    ost << "Lapse function on the AH :       "
	<< (hole.get_lapse_tot()).val_grid_point(1,0,nt-1,0) << endl ;
    ost << "Conformal factor on the AH :     "
	<< (hole.get_confo_tot()).val_grid_point(1,0,nt-1,0) << endl ;
    ost << "shift(1) on the AH :             "
	<< (hole.get_shift_tot()(1)).val_grid_point(1,0,nt-1,0) << endl ;
    ost << "shift(2) on the AH :             "
	<< (hole.get_shift_tot()(2)).val_grid_point(1,0,nt-1,0) << endl ;
    ost << "shift(3) on the AH :             "
	<< (hole.get_shift_tot()(3)).val_grid_point(1,0,nt-1,0) << endl ;

    ost << endl << "Neutron star : " ;
    ost << endl << "============   " << endl ;

    ost << "Baryon mass of NS in isolation :       "
    	<< star.mass_b()/msol << " [Mo]" << endl ;
    ost << "Gravitational mass of NS :             "
	<< star.mass_g_bhns()/msol << " [Mo]" << endl ;
    ost << "Baryon mass of NS in BH bg. :          "
	<< star.mass_b_bhns(hole.is_kerrschild(),hole.get_mass_bh(),separ)/msol
	<< " [Mo]" << endl ;
    ost << "Coordinate radius R_eq_tow :           "
	<< star.ray_eq_pi() / km << " [km]" << endl ;
    ost << "Coordinate radius R_eq_opp :           "
	<< star.ray_eq() / km << " [km]" << endl ;
    ost << "Coordinate radius R_eq_orb :           "
	<< star.ray_eq_pis2() / km << " [km]" << endl ;
    ost << "Coordinate radius R_eq_orb_opp :       "
	<< star.ray_eq_3pis2() / km << " [km]" << endl ;
    ost << "Coordinate radius R_pole :             "
	<< star.ray_pole() / km << " [km]" << endl ;
    ost << "Central enthalpy H_ent :               "
	<< (star.get_ent()).val_grid_point(0,0,0,0) << endl ;
    ost << "Lapse function at the center of NS :   "
	<< (star.get_lapse_tot()).val_grid_point(0,0,0,0) << endl ;
    ost << "Conformal factor at the center of NS : "
	<< (star.get_confo_tot()).val_grid_point(0,0,0,0) << endl ;
    ost << "shift(1) at the center of NS :         "
	<< (star.get_shift_tot()(1)).val_grid_point(0,0,0,0) << endl ;
    ost << "shift(2) at the center of NS :         "
	<< (star.get_shift_tot()(2)).val_grid_point(0,0,0,0) << endl ;
    ost << "shift(3) at the center of NS :         "
	<< (star.get_shift_tot()(3)).val_grid_point(0,0,0,0) << endl ;


    return ost ;

}

// Display in polytropic unites
// ----------------------------
void Bin_bhns::display_poly(ostream& ost) const {

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
	//	double j_poly = r_poly * r_poly / ggrav ;

	double r0 = 0.5 * ( star.ray_eq() + star.ray_eq_pi() ) ;
	            // = 1 in Baumgarte et al.
	//	double d_ns = separ + star.ray_eq() - r0 ;
	            // Orbital separation of Baumgarte et al.

	double y_ns = star.get_mp().get_ori_y() ;
	double d_coord = sqrt(separ*separ + y_ns*y_ns) ;

	ost.precision(16) ;
	ost << endl << "Quantities in polytropic units : " ;
	ost << endl << "==============================" << endl ;
	ost << " ( r_poly = " << r_poly / km << " km )" << endl ;
	ost << "  d_separ :      " << separ / r_poly << endl ;
	ost << "  d_coord :      " << d_coord / r_poly << endl ;
	ost << "  Omega :        " << omega * t_poly << endl ;
	ost << "  Omega M_irr :  "
	    << omega * t_poly * hole.mass_irr_bhns() / m_poly << endl ;
	ost << "  Omega_spin :   " << hole.get_omega_spin() * t_poly << endl ;
	ost << "  M_irr :        " << hole.mass_irr_bhns() / m_poly << endl ;
	ost << "  M_bh :         " << hole.get_mass_bh() / m_poly << endl ;
	ost << "  R_ah :         " << hole.rad_ah() / r_poly << endl ;
	ost << "  M_bar(NS)iso : " << star.mass_b() / m_poly
	    << endl ;
	ost << "  M_bar(NS) :    "
	    << star.mass_b_bhns(hole.is_kerrschild(),hole.get_mass_bh(),separ)
	  / m_poly << endl ;
	ost << "  M_g(NS)iso :   " << star.mass_g_bhns() / m_poly << endl ;
	ost << "  R_0(NS) :      " << r0 / r_poly << endl ;

    }

}
}
