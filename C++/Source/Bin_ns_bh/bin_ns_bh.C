/*
 *  Basic methods for class Bin_ns_bh
 *
 */

/*
 *   Copyright (c) 2002  Philippe Grandclement, Keisuke Taniguchi,
 *              Eric Gourgoulhon
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

char bin_ns_bh_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.7  2005/11/30 11:09:06  p_grandclement
 * Changes for the Bin_ns_bh project
 *
 * Revision 1.6  2005/08/29 15:10:15  p_grandclement
 * Addition of things needed :
 *   1) For BBH with different masses
 *   2) Provisory files for the mixted binaries (Bh and NS) : THIS IS NOT
 *   WORKING YET !!!
 *
 * Revision 1.5  2004/03/25 10:28:58  j_novak
 * All LORENE's units are now defined in the namespace Unites (in file unites.h).
 *
 * Revision 1.4  2003/02/13 16:40:25  p_grandclement
 * Addition of various things for the Bin_ns_bh project, non of them being
 * completely tested
 *
 * Revision 1.3  2002/12/19 14:51:19  e_gourgoulhon
 * Added the new functions set_omega and set_x_axe
 *
 * Revision 1.2  2002/12/18 10:31:15  e_gourgoulhon
 * irrot : int -> bool
 *
 * Revision 1.1  2002/12/17 13:10:11  e_gourgoulhon
 * Methods for class Bin_ns_bh
 *
 *
 *
 *
 * $Header$
 *
 */

// C++ headers
#include "headcpp.h"

// C headers
#include <math.h>

// Lorene headers
#include "map.h"
#include "bhole.h"
#include "bin_ns_bh.h"
#include "utilitaires.h"
#include "unites.h"
#include "graphique.h"

  			    //--------------//
			    // Constructors //
			    //--------------//

// Standard constructor
// --------------------
Bin_ns_bh::Bin_ns_bh(Map& mp_ns, int nzet, const Eos& eos, bool irrot_ns,
	        Map_af& mp_bh)
		: ref_triad(0., "Absolute frame Cartesian basis"),
		  star(mp_ns, nzet, true, eos, irrot_ns, ref_triad),
		  hole(mp_bh),
		  omega(0),
		  x_axe(0) {

     // Pointers of derived quantities initialized to zero :
    set_der_0x0() ;
}

// Copy constructor
// ----------------
Bin_ns_bh::Bin_ns_bh(const Bin_ns_bh& bibi)
                : ref_triad(0., "Absolute frame Cartesian basis"),
		 star(bibi.star),
		 hole(bibi.hole),
		 omega(bibi.omega),
		 x_axe(bibi.x_axe) {

     // Pointers of derived quantities initialized to zero :
    set_der_0x0() ;
}

// Constructor from a file
// -----------------------
Bin_ns_bh::Bin_ns_bh(Map& mp_ns, const Eos& eos, Map_af& mp_bh, FILE* fich)
		: ref_triad(0., "Absolute frame Cartesian basis"),
		  star(mp_ns, eos, ref_triad, fich),
		  hole(mp_bh, fich) {

    // omega and x_axe are read in the file:
    fread_be(&omega, sizeof(double), 1, fich) ;
    fread_be(&x_axe, sizeof(double), 1, fich) ;

    assert(hole.get_omega() == omega) ;

    // Pointers of derived quantities initialized to zero :
    set_der_0x0() ;
}

			    //------------//
			    // Destructor //
			    //------------//

Bin_ns_bh::~Bin_ns_bh(){

    del_deriv() ;

}

			//----------------------------------//
			// Management of derived quantities //
			//----------------------------------//

void Bin_ns_bh::del_deriv() const {

    if (p_mass_adm != 0x0) delete p_mass_adm ;
    if (p_mass_kom != 0x0) delete p_mass_kom ;
    if (p_angu_mom != 0x0) delete p_angu_mom ;
    if (p_total_ener != 0x0) delete p_total_ener ;
    if (p_virial != 0x0) delete p_virial ;
    if (p_virial_gb != 0x0) delete p_virial_gb ;
    if (p_virial_fus != 0x0) delete p_virial_fus ;
    if (p_ham_constr != 0x0) delete p_ham_constr ;
    if (p_mom_constr != 0x0) delete p_mom_constr ;

    set_der_0x0() ;
}




void Bin_ns_bh::set_der_0x0() const {

    p_mass_adm = 0x0 ;
    p_mass_kom = 0x0 ;
    p_angu_mom = 0x0 ;
    p_total_ener = 0x0 ;
    p_virial = 0x0 ;
    p_virial_gb = 0x0 ;
    p_virial_fus = 0x0 ;
    p_ham_constr = 0x0 ;
    p_mom_constr = 0x0 ;

}

			    //--------------//
			    //  Assignment  //
			    //--------------//

// Assignment to another Binaire
// -----------------------------

void Bin_ns_bh::operator=(const Bin_ns_bh& bibi) {

    assert( bibi.ref_triad == ref_triad ) ;

    star = bibi.star ;
    hole = bibi.hole ;

    omega = bibi.omega ;
    x_axe = bibi.x_axe ;

    // ref_triad remains unchanged

    del_deriv() ;  // Deletes all derived quantities

}

void Bin_ns_bh::set_omega(double omega_i) {

        omega = omega_i ;
	hole.set_omega(omega) ;

	del_deriv() ;

}

void Bin_ns_bh::set_x_axe(double x_axe_i) {

        x_axe = x_axe_i ;

	del_deriv() ;

}

double Bin_ns_bh::separation() const {
    return star.mp.get_ori_x() - hole.mp.get_ori_x() ;
}


			    //--------------//
			    //	  Outputs   //
			    //--------------//

// Save in a file
// --------------
void Bin_ns_bh::sauve(FILE* fich) const {

    star.sauve(fich) ;
    hole.sauve(fich) ;

    fwrite_be(&omega, sizeof(double), 1, fich) ;
    fwrite_be(&x_axe, sizeof(double), 1, fich) ;

}


void Bin_ns_bh::init_auto () {
	
	// On doit faire fonction pour assurer que tout va bien sur les trucs limites
	Cmp filtre_ns(star.get_mp()) ;
	Cmp radius_ns (star.get_mp()) ;
	radius_ns = star.get_mp().r ;
	double rlim_ns = star.get_mp().val_r (0, 1, 0, 0) ;
	filtre_ns = 0.5 * (1 + exp(-radius_ns*radius_ns/rlim_ns/rlim_ns)) ;
	filtre_ns.std_base_scal() ;
	
	Cmp filtre_bh(hole.get_mp()) ;
	Cmp radius_bh (hole.get_mp()) ;
	radius_bh = hole.get_mp().r ;
	double rlim_bh = hole.get_mp().val_r (0, 1, 0, 0) ;
	filtre_bh = 0.5 * (1 + exp(-radius_bh*radius_bh/rlim_bh/rlim_bh)) ;
	filtre_bh.std_base_scal() ;
	
	// Facteur conforme : pas de soucis
	star.set_confpsi_auto() = sqrt(exp(star.get_beta_auto()()-star.get_logn_auto()()))*filtre_ns ;
	star.set_confpsi_auto().set_std_base() ;
	hole.set_psi_auto() = hole.get_psi_auto()() * filtre_bh ;
	hole.set_psi_auto().std_base_scal() ;
	
	// Le lapse
	star.set_n_auto() = sqrt(exp(star.get_beta_auto()()-star.get_logn_auto()()))*filtre_ns ;
	star.set_n_auto().set_std_base() ;
		
	hole.set_n_auto() = hole.get_n_auto()() * filtre_bh ;
	hole.set_n_auto().std_base_scal() ;
	
	// On doit assurer que le lapse tot est bien zero sur l'horizon...
	Cmp soustrait ((filtre_bh-0.5)*2*exp(1.)) ;
	int nz = hole.get_mp().get_mg()->get_nzone() ;
	Mtbl xa_hole (hole.get_mp().get_mg()) ;
	xa_hole = hole.get_mp().xa ;
	Mtbl ya_hole (hole.get_mp().get_mg()) ;
	ya_hole = hole.get_mp().ya ;
	Mtbl za_hole (hole.get_mp().get_mg()) ;
	za_hole = hole.get_mp().za ;
	double xa_abs, ya_abs, za_abs ;
	double air, tet, phi ;

	int np = hole.get_mp().get_mg()->get_np(0) ;
	int nt = hole.get_mp().get_mg()->get_nt(0) ;
	for (int k=0 ; k<np ; k++)
	     for (int j=0 ; j<nt ; j++) {
	          double val_hole = hole.n_auto()(1, k,j,0) ;
		  xa_abs = xa_hole(1,k,j,0) ; 
		  ya_abs = ya_hole(1,k,j,0) ;
		  za_abs = za_hole(1,k,j,0) ;
		  star.get_mp().convert_absolute (xa_abs, ya_abs, za_abs, air, tet, phi) ;
		  double val_star  = star.get_n_auto()().val_point (air, tet, phi) ;
        	  for (int l=1 ; l<nz ; l++)
		      for (int i=0 ; i<hole.get_mp().get_mg()->get_nr(l) ; i++)
			   hole.set_n_auto().set(l,k,j,i) -= (val_star+val_hole)*soustrait(l,k,j,i) ;
	}
	hole.set_n_auto().std_base_scal() ;
	hole.set_n_auto().raccord(1) ;
}

	// *********************************
	// Affectation to another Bin_ns_bh
	//**********************************

void Bin_ns_bh::affecte(const Bin_ns_bh& so) {
        
	// Kinematic quantities :
	star.nzet = so.star.nzet ;
	set_omega(so.omega) ;
	x_axe = so.x_axe ;
	star.set_mp().set_ori (so.star.mp.get_ori_x(), 0., 0.) ;
        hole.set_mp().set_ori (so.hole.mp.get_ori_x(), 0., 0.) ;
   	// Faut gêrer le map_et :
	Map_et* map_et = dynamic_cast<Map_et*>(&star.mp) ;
	Map_et* map_et_so = dynamic_cast<Map_et*>(&so.star.mp) ;
	map_et->set_ff(map_et_so->get_ff()) ;
	map_et->set_gg(map_et_so->get_gg()) ;
	int l_min = (map_et->get_mg()->get_nzone() > map_et_so->get_mg()->get_nzone()) ?
	 	map_et->get_mg()->get_nzone() : map_et->get_mg()->get_nzone() ;
	for (int l=0 ; l<l_min+1 ; l++) {
	     map_et->set_alpha(map_et_so->get_alpha()[l], l) ;
	     map_et->set_beta(map_et_so->get_beta()[l],l) ;
        }
	
	cout << star.mp << endl ;
	cout << so.star.mp << endl ;
	abort() ;
	
	// The BH part :
	// Lapse :
	hole.n_auto.allocate_all() ;
	hole.n_auto.set().import(so.hole.n_auto()) ;
	hole.n_auto.set().std_base_scal() ;
	// Psi :
	hole.psi_auto.allocate_all() ;
	hole.psi_auto.set().import(so.hole.psi_auto()) ;
	hole.psi_auto.set().std_base_scal() ;
	// Shift :
	hole.shift_auto.allocate_all() ;
	for (int i=0 ; i<3 ; i++)
		hole.shift_auto.set(i).import (so.hole.shift_auto(i)) ;
	hole.shift_auto.set_std_base() ;
	
	// The NS part :
	star.n_auto.allocate_all() ;
	star.n_auto.set().import(so.star.n_auto()) ;
	star.n_auto.set().std_base_scal() ;
	// Psi :
	star.confpsi_auto.allocate_all() ;
	star.confpsi_auto.set().import(so.star.confpsi_auto()) ;
	star.confpsi_auto.set().std_base_scal() ;
	// Shift :
	star.shift.allocate_all() ;
	for (int i=0 ; i<3 ; i++)
		star.shift.set(i).import (so.star.shift(i)) ;
	star.shift_auto.set_std_base() ;
	
	// The matter : 
	star.ent.allocate_all() ;
	star.ent.set().import(star.nzet, so.star.ent()) ;
	star.ent.set_std_base() ;
	
	des_profile (star.ent(), 0, 3, 0, 0) ;
	des_profile (so.star.ent(), 0, 3, 0, 0) ;
	
	if (so.star.is_irrotational()) {
		star.d_psi.allocate_all() ;
		for (int i=0 ; i< 3 ; i++)
	        star.d_psi.set(i).import(star.nzet, so.star.d_psi(i)) ;
		star.d_psi.set_std_base() ;
	}
	
	star.kinematics(omega, x_axe) ;
        star.fait_d_psi() ;
        star.hydro_euler() ;

	// Reconstruction of the fields :
	hole.update_metric(star) ;    
        star.update_metric(hole) ;
	star.update_metric_der_comp(hole) ;
	fait_tkij() ;
}
	

// Printing
// --------
ostream& operator<<(ostream& ost, const Bin_ns_bh& bibi)  {
    bibi >> ost ;
    return ost ;
}


ostream& Bin_ns_bh::operator>>(ostream& ost) const {

  using namespace Unites ;

    ost << endl ;
    ost << "Neutron star - black hole binary system" << endl ;
    ost << "=======================================" << endl ;
    ost << endl <<
	"Orbital angular velocity : " << omega * f_unit << " rad/s" << endl ;
    ost <<
	"Absolute coordinate X of the rotation axis : " << x_axe / km
	    << " km" << endl ;
    ost << endl << "Neutron star : " << endl ;
    ost <<         "============   " << endl ;
    ost << star << endl ;

    ost << "Black hole : " << endl ;
    ost << "==========   " << endl ;
    ost << "Coordinate radius of the throat : " << hole.get_rayon() / km << " km" << endl ;
    ost << "Absolute abscidia of the throat center : " << (hole.get_mp()).get_ori_x() / km
        << " km" << endl ;
    return ost ;
}
