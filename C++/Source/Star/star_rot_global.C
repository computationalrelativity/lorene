/*
 * Methods for computing global quantities within the class Star_rot
 *
 * (see file star_rot.h for documentation)
 */

/*
 *   Copyright (c) 2010 Eric Gourgoulhon
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


char star_rot_global_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.5  2015/05/19 09:30:56  j_novak
 * New methods for computing the area of the star and its mean radius.
 *
 * Revision 1.4  2014/10/13 08:53:39  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.3  2014/10/06 15:13:17  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2010/01/25 22:33:35  e_gourgoulhon
 * Debugging...
 *
 * Revision 1.1  2010/01/25 18:15:52  e_gourgoulhon
 * First version.
 * 
 *
 *
 * $Header$
 *
 */

// Headers C
#include <cstdlib>
#include <cmath>

// Headers Lorene
#include "star_rot.h"
#include "unites.h"

			//--------------------------//
			//	Stellar surface	    //
			//--------------------------//

namespace Lorene {
  const Itbl& Star_rot::l_surf() const {
    
    if (p_l_surf == 0x0) {    // a new computation is required
      
      assert(p_xi_surf == 0x0) ;  // consistency check
      
      int np = mp.get_mg()->get_np(0) ;   
      int nt = mp.get_mg()->get_nt(0) ;   
      
      p_l_surf = new Itbl(np, nt) ;
      p_xi_surf = new Tbl(np, nt) ;
      
      double ent0 = 0 ;	// definition of the surface
      double precis = 1.e-15 ; 
      int nitermax = 100 ; 
      int niter ; 
      
      (ent.get_spectral_va()).equipot(ent0, nzet, precis, nitermax, niter, *p_l_surf, 
				      *p_xi_surf) ; 
      
    }
    
    return *p_l_surf ; 
    
  }
  
                    //------------------------------//
                    //	       Baryon mass	    //
                    //------------------------------//
  
  double Star_rot::mass_b() const {
    
    if (p_mass_b == 0x0) {    // a new computation is required
      
      if (relativistic) {
	
	Scalar dens = a_car * bbb * gam_euler * nbar ;
	
	p_mass_b = new double( dens.integrale() ) ;	
      }
      else{  // Newtonian case 
	assert(nbar.get_etat() == ETATQCQ) ; 
	
	p_mass_b = new double( nbar.integrale() ) ;
	
      }
      
    }
    return *p_mass_b ; 
  } 


			//----------------------------//
			//	Gravitational mass    //
			//----------------------------//

  double Star_rot::mass_g() const {
    
    if (p_mass_g == 0x0) {    // a new computation is required
      
      if (relativistic) {
	
	Scalar source = nn * (ener_euler + s_euler) 
	  + 2 * bbb * (ener_euler + press)
	  * tnphi * uuu ; 
	source = a_car * bbb * source ;
	source.std_spectral_base() ; 
	
	p_mass_g = new double( source.integrale() ) ;	
      }
      else{  // Newtonian case 
	p_mass_g = new double( mass_b() ) ;   // in the Newtonian case
	                                      //  M_g = M_b
      }
    }
    return *p_mass_g ; 
} 
		
			//----------------------------//
			//	Angular momentum      //
			//----------------------------//

  double Star_rot::angu_mom() const {
    
    if (p_angu_mom == 0x0) {    // a new computation is required
      
      Scalar dens = uuu ; 
      
      dens.mult_r() ;			//  Multiplication by
      dens.set_spectral_va() = (dens.get_spectral_va()).mult_st() ;	//    r sin(theta)
      
      if (relativistic) {
	dens = a_car * b_car * (ener_euler + press) * dens ; 
      }
      else {    // Newtonian case 
	dens = nbar * dens ; 
      }
      
      p_angu_mom = new double( dens.integrale() ) ; 
    }
    return *p_angu_mom ; 
  }


			//----------------------------//
			//	     T/W	      //
			//----------------------------//

  double Star_rot::tsw() const {
    
    if (p_tsw == 0x0) {    // a new computation is required
      
      double tcin = 0.5 * omega * angu_mom() ;
      
      if (relativistic) {
	
	Scalar dens = a_car * bbb * gam_euler * ener ;
	double mass_p = dens.integrale() ; 
	
	p_tsw = new double( tcin / ( mass_p + tcin - mass_g() ) ) ;  	
	
      }
      else {	    // Newtonian case 
	Scalar dens = 0.5 * nbar * logn ;
	double wgrav = dens.integrale() ; 
	p_tsw = new double( tcin / fabs(wgrav) ) ;  
      } 
    }
    return *p_tsw ; 
}


			//----------------------------//
			//	     GRV2	      //
			//----------------------------//

double Star_rot::grv2() const {

      using namespace Unites ;	
      if (p_grv2 == 0x0) {    // a new computation is required
	
	Scalar sou_m(mp) ; 
	if (relativistic) {
	  sou_m =  2 * qpig * a_car * (press + (ener_euler+press)
				       * uuu*uuu ) ;
        }
	else {
	  sou_m = 2 * qpig * (press + nbar * uuu*uuu ) ;
	}
	
	Vector dlogn = logn.derive_cov( mp.flat_met_spher() ) ; 
        Scalar sou_q =  1.5 * ak_car
	  - dlogn(1)*dlogn(1) - dlogn(2)*dlogn(2) - dlogn(3)*dlogn(3) ;	
	
	p_grv2 = new double( double(1) - lambda_grv2(sou_m, sou_q) ) ; 
	
      }
    
      return *p_grv2 ; 
      
}


			//----------------------------//
			//	     GRV3	      //
			//----------------------------//

double Star_rot::grv3(ostream* ost) const {

  using namespace Unites ;	
  
  if (p_grv3 == 0x0) {    // a new computation is required

    Scalar source(mp) ; 
    
    // Gravitational term [cf. Eq. (43) of Gourgoulhon & Bonazzola
    // ------------------	    Class. Quantum Grav. 11, 443 (1994)]
    
    Vector dlogn = logn.derive_cov( mp.flat_met_spher() ) ;

    if (relativistic) {
      Scalar alpha = dzeta - logn ; 
      Scalar bet = log( bbb ) ; 
      bet.std_spectral_base() ; 
      
      Vector dalpha = alpha.derive_cov( mp.flat_met_spher() ) ;
      Vector dbet = bet.derive_cov( mp.flat_met_spher() ) ;

      source = 0.75 * ak_car 
	- dlogn(1)*dlogn(1) - dlogn(2)*dlogn(2) - dlogn(3)*dlogn(3) 
	+ 0.5 * ( dalpha(1)*dbet(1) + dalpha(2)*dbet(2) + dalpha(3)*dbet(3) ) ; 
      
      Scalar aa = alpha - 0.5 * bet ; 
      Scalar daadt = aa.srdsdt() ;	// 1/r d/dth
	    
      // What follows is valid only for a mapping of class Map_radial : 
      const Map_radial* mpr = dynamic_cast<const Map_radial*>(&mp) ; 
      if (mpr == 0x0) {
	cout << "Star_rot::grv3: the mapping does not belong"
	     << " to the class Map_radial !" << endl ; 
	abort() ; 
      }
      
      // Computation of 1/tan(theta) * 1/r daa/dtheta
      if (daadt.get_etat() == ETATQCQ) {
	Valeur& vdaadt = daadt.set_spectral_va() ; 
	vdaadt = vdaadt.ssint() ;	// division by sin(theta)
	vdaadt = vdaadt.mult_ct() ;	// multiplication by cos(theta)
      }
      
      Scalar temp = aa.dsdr() + daadt ; 
      temp = ( bbb - a_car/bbb ) * temp ; 
      temp.std_spectral_base() ; 
      
      // Division by r 
      Valeur& vtemp = temp.set_spectral_va() ; 
      vtemp = vtemp.sx() ;    // division by xi in the nucleus
      // Id in the shells
      // division by xi-1 in the ZEC
      vtemp = (mpr->xsr) * vtemp ; // multiplication by xi/r in the nucleus
      //		  by 1/r in the shells
      //		  by r(xi-1) in the ZEC
      
      // In the ZEC, a multiplication by r has been performed instead
      //   of the division: 			
      temp.set_dzpuis( temp.get_dzpuis() + 2 ) ;  
      
      source = bbb * source + 0.5 * temp ; 
      
    }
    else{
      source = - 0.5 * ( dlogn(1)*dlogn(1) + dlogn(2)*dlogn(2) + dlogn(3)*dlogn(3) ) ; 
    }
    
    source.std_spectral_base() ; 
    
    double int_grav = source.integrale() ; 
    
    // Matter term
    // -----------
    
    if (relativistic) {    
      source  = qpig * a_car * bbb * s_euler ;
    }
    else{
      source = qpig * ( 3 * press + nbar * uuu * uuu ) ; 
    }
    
    source.std_spectral_base() ; 

    double int_mat = source.integrale() ; 
    
    // Virial error
    // ------------
    if (ost != 0x0) {
      *ost << "Star_rot::grv3 : gravitational term : " << int_grav 
	   << endl ;
      *ost << "Star_rot::grv3 : matter term :        " << int_mat 
	   << endl ;
    }
    
    p_grv3 = new double( (int_grav + int_mat) / int_mat ) ; 	 
    
  }
  
  return *p_grv3 ; 
  
}


			//----------------------------//
			//	     R_circ	      //
			//----------------------------//

double Star_rot::r_circ() const {

    if (p_r_circ == 0x0) {    // a new computation is required
	
	// Index of the point at phi=0, theta=pi/2 at the surface of the star:
	const Mg3d* mg = mp.get_mg() ; 
	assert(mg->get_type_t() == SYM) ; 
	int l_b = nzet - 1 ; 
	int i_b = mg->get_nr(l_b) - 1 ; 
	int j_b = mg->get_nt(l_b) - 1 ; 
	int k_b = 0 ; 
    
	p_r_circ = new double( bbb.val_grid_point(l_b, k_b, j_b, i_b) * ray_eq() ) ; 

    }
    
    return *p_r_circ ; 

}

			//----------------------------//
			//	  Surface area	      //
			//----------------------------//

  double Star_rot::area() const {

    if (p_area == 0x0) {    // a new computation is required
      const Mg3d& mg = *(mp.get_mg()) ; 
      int np = mg.get_np(0) ;
      int nt = mg.get_nt(0) ;
      assert(np == 1) ; //Only valid for axisymmetric configurations
      
      const Map_radial* mp_rad = dynamic_cast<const Map_radial*>(&mp) ;
      assert(mp_rad != 0x0) ;

      Valeur va_r(mg.get_angu()) ;
      va_r.annule_hard() ;
      Itbl lsurf = l_surf() ;
      Tbl xisurf = xi_surf() ;
      for (int k=0; k<np; k++) {
	for (int j=0; j<nt; j++) {
	  int l_star = lsurf(k, j) ;
	  double xi_star = xisurf(k, j) ;
	  va_r.set(0, k, j, 0) = mp_rad->val_r_jk(l_star, xi_star, j, k) ;
	}
      }
      va_r.std_base_scal() ;
      
      Valeur integ(mg.get_angu()) ;
      Valeur dr = va_r.dsdt() ;
      integ = sqrt(va_r*va_r + dr*dr) ;
      Valeur a2 = get_a_car().get_spectral_va() ; a2.std_base_scal() ;
      Valeur b = get_bbb().get_spectral_va() ; b.std_base_scal() ;
      for (int k=0; k<np; k++) {
	for (int j=0; j<nt; j++) {
	  integ.set(0, k, j, 0) *= sqrt(a2.val_point_jk(lsurf(k, j), xisurf(k, j), j, k))
	    * b.val_point_jk(lsurf(k, j), xisurf(k, j), j, k) * va_r(0, k, j, 0) ;
	}
      }
      integ.std_base_scal() ;
      Valeur integ2 = integ.mult_st() ;
      double surftot = 0. ;
      for (int j=0; j<nt; j++) {
	surftot += (*integ2.c_cf)(0, 0, j, 0) / double(2*j+1) ;
      }
      
      p_area = new double( 4*M_PI*surftot) ;

    }
    
    return *p_area ; 

}

  double Star_rot::mean_radius() const {

    return sqrt(area()/(4*M_PI)) ;

  }

			//----------------------------//
			//	   Flattening	      //
			//----------------------------//

double Star_rot::aplat() const {

    if (p_aplat == 0x0) {    // a new computation is required
	
	p_aplat = new double( ray_pole() / ray_eq() ) ; 	 

    }
    
    return *p_aplat ; 

}


			//----------------------------//
			//	     Z_eq_f	      //
			//----------------------------//

double Star_rot::z_eqf() const {

    if (p_z_eqf == 0x0) {    // a new computation is required
	
	// Index of the point at phi=0, theta=pi/2 at the surface of the star:
	const Mg3d* mg = mp.get_mg() ; 
	assert(mg->get_type_t() == SYM) ; 
	int l_b = nzet - 1 ; 
	int i_b = mg->get_nr(l_b) - 1 ; 
	int j_b = mg->get_nt(l_b) - 1 ; 
	int k_b = 0 ; 
    
	double u_eq = uuu.val_grid_point(l_b, k_b, j_b, i_b) ; 
	double n_eq = nn.val_grid_point(l_b, k_b, j_b, i_b) ; 
	double nphi_eq = nphi.val_grid_point(l_b, k_b, j_b, i_b) ; 
	
	p_z_eqf = new double( sqrt((1.-u_eq)/(1.+u_eq)) 
				/ (n_eq + nphi_eq * r_circ() )
			      - 1. ) ;
    }
    
    return *p_z_eqf ; 

}
			//----------------------------//
			//	     Z_eq_b	      //
			//----------------------------//

double Star_rot::z_eqb() const {

    if (p_z_eqb == 0x0) {    // a new computation is required
	
	// Index of the point at phi=0, theta=pi/2 at the surface of the star:
	const Mg3d* mg = mp.get_mg() ; 
	assert(mg->get_type_t() == SYM) ; 
	int l_b = nzet - 1 ; 
	int i_b = mg->get_nr(l_b) - 1 ; 
	int j_b = mg->get_nt(l_b) - 1 ; 
	int k_b = 0 ; 
    
	double u_eq = uuu.val_grid_point(l_b, k_b, j_b, i_b) ; 
	double n_eq = nn.val_grid_point(l_b, k_b, j_b, i_b) ; 
	double nphi_eq = nphi.val_grid_point(l_b, k_b, j_b, i_b) ; 
	
	p_z_eqb = new double(  sqrt((1.+u_eq)/(1.-u_eq)) 
				/ (n_eq - nphi_eq * r_circ() )
			      - 1. )  ;

    }
    
    return *p_z_eqb ; 

}


			//----------------------------//
			//	     Z_pole	      //
			//----------------------------//

double Star_rot::z_pole() const {

    if (p_z_pole == 0x0) {    // a new computation is required
	
	double n_pole = nn.val_point(ray_pole(), 0., 0.) ; 
	
	p_z_pole = new double(  1. / n_pole - 1. ) ; 

    }
    
    return *p_z_pole ; 

}


			//----------------------------//
			//     Quadrupole moment      //
			//----------------------------//

double Star_rot::mom_quad() const {

  using namespace Unites ;
    if (p_mom_quad == 0x0) {    // a new computation is required
	
	// Source for of the Poisson equation for nu
	// -----------------------------------------

	Scalar source(mp) ; 
	
	if (relativistic) {
	    Scalar bet = log(bbb) ; 
	    bet.std_spectral_base() ; 
	    
	    Vector dlogn = logn.derive_cov( mp.flat_met_spher() ) ;
	    Vector dlogn_bet = dlogn + bet.derive_cov( mp.flat_met_spher() ) ;

	    source =  qpig * a_car *( ener_euler + s_euler ) + ak_car 
	      - dlogn(1)*dlogn_bet(1) - dlogn(2)*dlogn_bet(2) - dlogn(3)*dlogn_bet(3)   ; 
	}
	else {
	    source = qpig * nbar ; 
	}
	
	source.std_spectral_base() ;

	// Multiplication by -r^2 P_2(cos(theta))
	//  [cf Eq.(7) of Salgado et al. Astron. Astrophys. 291, 155 (1994) ]
	// ------------------------------------------------------------------
	
	// Multiplication by r^2 : 
	// ----------------------
	source.mult_r() ; 
	source.mult_r() ; 
	if (source.check_dzpuis(2)) {
	    source.inc_dzpuis(2) ; 
	}
		
	// Muliplication by cos^2(theta) :
	// -----------------------------
	Scalar temp = source ; 
	
	// What follows is valid only for a mapping of class Map_radial : 
	assert( dynamic_cast<const Map_radial*>(&mp) != 0x0 ) ; 
		
	if (temp.get_etat() == ETATQCQ) {
	    Valeur& vtemp = temp.set_spectral_va() ; 
	    vtemp = vtemp.mult_ct() ;	// multiplication by cos(theta)
	    vtemp = vtemp.mult_ct() ;	// multiplication by cos(theta)
	}
	
	// Muliplication by -P_2(cos(theta)) :
	// ----------------------------------
	source = 0.5 * source - 1.5 * temp ; 
	
	// Final result
	// ------------

	p_mom_quad = new double( source.integrale() / qpig ) ; 	 

    }
    
    return *p_mom_quad ; 

}




}
