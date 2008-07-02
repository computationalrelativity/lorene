/*
 *  Methods of class Hole_bhns to compute global quantities
 *
 *    (see file hole_bhns.h for documentation).
 *
 */

/*
 *   Copyright (c) 2005,2007 Keisuke Taniguchi
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

char hole_bhns_global_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2008/07/02 21:10:15  k_taniguchi
 * A bug removed.
 *
 * Revision 1.2  2008/05/15 19:07:26  k_taniguchi
 * Introduction of the quasilocal spin angular momentum.
 *
 * Revision 1.1  2007/06/22 01:25:15  k_taniguchi
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// C++ headers
//#include <>

// C headers
#include <math.h>

// Lorene headers
#include "hole_bhns.h"
#include "unites.h"
#include "utilitaires.h"

                    //-----------------------------------------//
                    //          Irreducible mass of BH         //
                    //-----------------------------------------//

double Hole_bhns::mass_irr_bhns() const {

    // Fundamental constants and units
    // -------------------------------
    using namespace Unites ;

    if (p_mass_irr_bhns == 0x0) {   // a new computation is required

        Scalar psi4(mp) ;
	psi4 = pow(confo_tot, 4.) ;
	psi4.std_spectral_base() ;
	psi4.annule_domain(0) ;
	psi4.raccord(1) ;

	double radius_ah = mp.val_r(1,-1.,M_PI/2.,0.) ;

	Map_af& mp_aff= dynamic_cast<Map_af&>(mp) ;

	double a_ah = mp_aff.integrale_surface(psi4, radius_ah) ;
	double mirr = sqrt(a_ah/16./M_PI) / ggrav ;

	p_mass_irr_bhns = new double( mirr ) ;

    }

    return *p_mass_irr_bhns ;

}

          //----------------------------------------------------------//
          //          Quasilocal spin angular momentum of BH          //
          //----------------------------------------------------------//

double Hole_bhns::spin_am_bhns(const Tbl& xi_i, const double& phi_i,
			       const double& theta_i, const int& nrk_phi,
			       const int& nrk_theta) const {

    // Fundamental constants and units
    // -------------------------------
    using namespace Unites ;

    if (p_spin_am_bhns == 0x0) {   // a new computation is required

        double mass = ggrav * mass_bh ;

        Scalar rr(mp) ;
	rr = mp.r ;
	rr.std_spectral_base() ;

	Scalar st(mp) ;
	st = mp.sint ;
	st.std_spectral_base() ;
	Scalar ct(mp) ;
	ct = mp.cost ;
	ct.std_spectral_base() ;
	Scalar sp(mp) ;
	sp = mp.sinp ;
	sp.std_spectral_base() ;
	Scalar cp(mp) ;
	cp = mp.cosp ;
	cp.std_spectral_base() ;

	Vector ll(mp, CON, mp.get_bvect_cart()) ;
	ll.set_etat_qcq() ;
	ll.set(1) = st % cp ;
	ll.set(2) = st % sp ;
	ll.set(3) = ct ;
	ll.std_spectral_base() ;

	double radius_ah = mp.val_r(1,-1.,M_PI/2.,0.) ;

	if (kerrschild) {

	  cout << "Not yet prepared!!!" << endl ;
	  abort() ;

	}
	else {  // Isotropic coordinates

	  // Sets C/M^2 for each case of the lapse boundary condition
	  // --------------------------------------------------------
	  double cc ;

	  if (bc_lapconf_nd) {  // Neumann boundary condition
	    if (bc_lapconf_fs) {  // First condition
	      // d(\alpha \psi)/dr = 0
	      // ---------------------
	      cc = 2. * (sqrt(13.) - 1.) / 3. ;
	    }
	    else {  // Second condition
	      // d(\alpha \psi)/dr = (\alpha \psi)/(2 rah)
	      // -----------------------------------------
	      cc = 4. / 3. ;
	    }
	  }
	  else {  // Dirichlet boundary condition
	    if (bc_lapconf_fs) {  // First condition
	      // (\alpha \psi) = 1/2
	      // -------------------
	      cout << "!!!!! WARNING: Not yet prepared !!!!!" << endl ;
	      abort() ;
	    }
	    else {  // Second condition
	      // (\alpha \psi) = 1/sqrt(2.) \psi_KS
	      // ----------------------------------
	      cout << "!!!!! WARNING: Not yet prepared !!!!!" << endl ;
	      abort() ;
	      //	      cc = 2. * sqrt(2.) ;
	    }
	  }

	  Scalar r_are(mp) ;
	  r_are = r_coord(bc_lapconf_nd, bc_lapconf_fs) ;
	  r_are.std_spectral_base() ;

	  // Killing vector of the spherical components
	  Vector killing_spher(mp, COV, mp.get_bvect_spher()) ;
	  killing_spher.set_etat_qcq() ;
	  killing_spher = killing_vect(xi_i, phi_i, theta_i,
				       nrk_phi, nrk_theta) ;
	  killing_spher.std_spectral_base() ;

	  killing_spher.set(2) = confo_tot * confo_tot * radius_ah
	    * killing_spher(2) ;
	  killing_spher.set(3) = confo_tot * confo_tot * radius_ah
	    * killing_spher(3) ;
	  // killing_spher(3) is already divided by sin(theta)
	  killing_spher.std_spectral_base() ;

	  // Killing vector of the Cartesian components
	  Vector killing(mp, COV, mp.get_bvect_cart()) ;
	  killing.set_etat_qcq() ;
	  killing.set(1) = (killing_spher(2) * ct * cp - killing_spher(3) * sp)
	    / radius_ah ;
	  killing.set(2) = (killing_spher(2) * ct * sp + killing_spher(3) * cp)
	    / radius_ah ;
	  killing.set(3) = - killing_spher(2) * st / radius_ah ;
	  killing.std_spectral_base() ;

	  // Surface integral <- dzpuis should be 0
	  // --------------------------------------
	  // Source terms in the surface integral
	  Scalar source_1(mp) ;
	  source_1 = (ll(1) * (taij_tot_rs(1,1) * killing(1)
			       + taij_tot_rs(1,2) * killing(2)
			       + taij_tot_rs(1,3) * killing(3))
		      + ll(2) * (taij_tot_rs(2,1) * killing(1)
				 + taij_tot_rs(2,2) * killing(2)
				 + taij_tot_rs(2,3) * killing(3))
		      + ll(3) * (taij_tot_rs(3,1) * killing(1)
				 + taij_tot_rs(3,2) * killing(2)
				 + taij_tot_rs(3,3) * killing(3)))
	    / pow(confo_tot, 4.) ;
	  source_1.std_spectral_base() ;
	  source_1.dec_dzpuis(2) ;  // dzpuis : 2 -> 0

	  Scalar source_2(mp) ;
	  source_2 = -2. * pow(confo_tot, 3.) * mass * mass * cc
	    * sqrt(1. - 2.*mass/r_are/rr + cc*cc*pow(mass/r_are/rr,4.))
	    * (ll(1)*killing(1) + ll(2)*killing(2) + ll(3)*killing(3))
	    / lapconf_tot / pow(r_are*rr, 3.) ;
	  source_2.std_spectral_base() ;

	  Scalar source_surf(mp) ;
	  source_surf = source_1 + source_2 ;
	  source_surf.std_spectral_base() ;
	  source_surf.annule_domain(0) ;
	  source_surf.raccord(1) ;

	  Map_af& mp_aff= dynamic_cast<Map_af&>(mp) ;

	  double spin = mp_aff.integrale_surface(source_surf, radius_ah) ;
	  double spin_angmom = 0.5 * spin / qpig ;

	  p_spin_am_bhns = new double( spin_angmom ) ;

	}

    }

    return *p_spin_am_bhns ;

}
