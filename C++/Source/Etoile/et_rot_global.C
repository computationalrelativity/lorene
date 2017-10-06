/*
 * Methods for computing global quantities within the class Etoile_rot
 *
 * (see file etoile.h for documentation)
 */

/*
 *   Copyright (c) 2000-2001 Eric Gourgoulhon
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


 

/*
 * $Id$
 * $Log$
 * Revision 1.11  2017/10/06 12:36:34  a_sourie
 * Cleaning of tabulated 2-fluid EoS class + superfluid rotating star model.
 *
 * Revision 1.10  2016/12/05 16:17:54  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.9  2015/06/12 12:38:25  j_novak
 * Implementation of the corrected formula for the quadrupole momentum.
 *
 * Revision 1.8  2015/06/10 14:37:44  a_sourie
 * Corrected the formula for the quadrupole.
 *
 * Revision 1.7  2014/10/13 08:52:57  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.6  2014/10/06 15:13:09  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.5  2012/08/12 17:48:35  p_cerda
 * Magnetstar: New classes for magnetstar. Allowing for non-equatorial symmetry in Etoile et al. Adding B_phi in Et_rot_mag.
 *
 * Revision 1.4  2004/03/25 10:29:06  j_novak
 * All LORENE's units are now defined in the namespace Unites (in file unites.h).
 *
 * Revision 1.3  2003/11/03 16:47:13  e_gourgoulhon
 * Corrected error in grv2() in the Newtonian case.
 *
 * Revision 1.2  2002/04/05 09:09:37  j_novak
 * The inversion of the EOS for 2-fluids polytrope has been modified.
 * Some errors in the determination of the surface were corrected.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 1.5  2000/11/19  18:52:09  eric
 * grv2() operationnelle.
 *
 * Revision 1.4  2000/10/12  15:34:55  eric
 * Calcul de la masse grav, de GRV3 et du moment quadrupolaire.
 *
 * Revision 1.3  2000/08/31  11:25:58  eric
 * *** empty log message ***
 *
 * Revision 1.2  2000/08/25  12:28:16  eric
 * *** empty log message ***
 *
 * Revision 1.1  2000/07/20  15:32:56  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */

// Headers C
#include <cstdlib>
#include <cmath>

// Headers Lorene
#include "etoile.h"
#include "unites.h"

			//--------------------------//
			//	Stellar surface	    //
			//--------------------------//

namespace Lorene {
const Itbl& Etoile_rot::l_surf() const {

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
	
	(ent().va).equipot(ent0, nzet, precis, nitermax, niter, *p_l_surf, 
		    *p_xi_surf) ; 
    
    }
   
    return *p_l_surf ; 
    
}

			//--------------------------//
			//	Baryon mass	    //
			//--------------------------//

double Etoile_rot::mass_b() const {

    if (p_mass_b == 0x0) {    // a new computation is required
	
	if (relativistic) {

	    Cmp dens = a_car() * bbb() * gam_euler() * nbar() ;
	    
	    dens.std_base_scal() ; 

	    p_mass_b = new double( dens.integrale() ) ;


	}
	else{  // Newtonian case 
	    assert(nbar.get_etat() == ETATQCQ) ; 

	    p_mass_b = new double( nbar().integrale() ) ;

	}

    }
    
    return *p_mass_b ; 

} 


			//----------------------------//
			//	Gravitational mass    //
			//----------------------------//

double Etoile_rot::mass_g() const {

    if (p_mass_g == 0x0) {    // a new computation is required
	
	if (relativistic) {

	    Tenseur source = nnn * (ener_euler + s_euler) 
				+ 2 * bbb * (ener_euler + press)
				    * tnphi * uuu ; 
	    source = a_car * bbb * source ;

	    source.set_std_base() ;

	    p_mass_g = new double( source().integrale() ) ;


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

double Etoile_rot::angu_mom() const {

    if (p_angu_mom == 0x0) {    // a new computation is required
	
	Cmp dens = uuu() ; 

	dens.mult_r() ;			//  Multiplication by
	dens.va = (dens.va).mult_st() ;	//    r sin(theta)

	if (relativistic) {
	    dens = a_car() * b_car() * (ener_euler() + press()) 
			* dens ; 
	}
	else {    // Newtonian case 
	    dens = nbar() * dens ; 
	}

	dens.std_base_scal() ; 

	p_angu_mom = new double( dens.integrale() ) ;

    }
    
    return *p_angu_mom ; 

}


			//----------------------------//
			//	     T/W	      //
			//----------------------------//

double Etoile_rot::tsw() const {

    if (p_tsw == 0x0) {    // a new computation is required
	
	double tcin = 0.5 * omega * angu_mom() ;
	
	if (relativistic) {
	    
	    Cmp dens = a_car() * bbb() * gam_euler() * ener() ;
	    dens.std_base_scal() ; 
	    double mass_p = dens.integrale() ; 
	    
	    p_tsw = new double( tcin / ( mass_p + tcin - mass_g() ) ) ;  	
	   
	}
	else {	    // Newtonian case 
	    Cmp dens = 0.5 * nbar() * logn() ;
	    dens.std_base_scal() ; 
	    double wgrav = dens.integrale() ; 
	    p_tsw = new double( tcin / fabs(wgrav) ) ;  
	}


    }
    
    return *p_tsw ; 

}


			//----------------------------//
			//	     GRV2	      //
			//----------------------------//

double Etoile_rot::grv2() const {

      using namespace Unites ;	
      if (p_grv2 == 0x0) {    // a new computation is required
	
	Tenseur sou_m(mp) ; 
	if (relativistic) {
	  sou_m =  2 * qpig * a_car * (press + (ener_euler+press)
				       * uuu*uuu ) ;
        }
	else {
	  sou_m = 2 * qpig * (press + nbar * uuu*uuu ) ;
	}
	
        Tenseur sou_q =  1.5 * ak_car
	  - flat_scalar_prod(logn.gradient_spher(),
			     logn.gradient_spher() ) ;	
	
	p_grv2 = new double( double(1) - lambda_grv2(sou_m(), sou_q()) ) ; 
	
      }
    
      return *p_grv2 ; 
      
}


			//----------------------------//
			//	     GRV3	      //
			//----------------------------//

double Etoile_rot::grv3(ostream* ost) const {

  using namespace Unites ;	
  
  if (p_grv3 == 0x0) {    // a new computation is required

    Tenseur source(mp) ; 
    
    // Gravitational term [cf. Eq. (43) of Gourgoulhon & Bonazzola
    // ------------------	    Class. Quantum Grav. 11, 443 (1994)]
    
    if (relativistic) {
      Tenseur alpha = dzeta - logn ; 
      Tenseur beta = log( bbb ) ; 
      beta.set_std_base() ; 
      
      source = 0.75 * ak_car 
	- flat_scalar_prod(logn.gradient_spher(),
			   logn.gradient_spher() )
	+ 0.5 * flat_scalar_prod(alpha.gradient_spher(),
				 beta.gradient_spher() ) ; 
      
      Cmp aa = alpha() - 0.5 * beta() ; 
      Cmp daadt = aa.srdsdt() ;	// 1/r d/dth
	    
      // What follows is valid only for a mapping of class Map_radial : 
      const Map_radial* mpr = dynamic_cast<const Map_radial*>(&mp) ; 
      if (mpr == 0x0) {
	cout << "Etoile_rot::grv3: the mapping does not belong"
	     << " to the class Map_radial !" << endl ; 
	abort() ; 
      }
      
      // Computation of 1/tan(theta) * 1/r daa/dtheta
      if (daadt.get_etat() == ETATQCQ) {
	Valeur& vdaadt = daadt.va ; 
	vdaadt = vdaadt.ssint() ;	// division by sin(theta)
	vdaadt = vdaadt.mult_ct() ;	// multiplication by cos(theta)
      }
      
      Cmp temp = aa.dsdr() + daadt ; 
      temp = ( bbb() - a_car()/bbb() ) * temp ; 
      temp.std_base_scal() ; 
      
      // Division by r 
      Valeur& vtemp = temp.va ; 
      vtemp = vtemp.sx() ;    // division by xi in the nucleus
      // Id in the shells
      // division by xi-1 in the ZEC
      vtemp = (mpr->xsr) * vtemp ; // multiplication by xi/r in the nucleus
      //		  by 1/r in the shells
      //		  by r(xi-1) in the ZEC
      
      // In the ZEC, a multiplication by r has been performed instead
      //   of the division: 			
      temp.set_dzpuis( temp.get_dzpuis() + 2 ) ;  
      
      source = bbb() * source() + 0.5 * temp ; 
      
    }
    else{
      source = - 0.5 * flat_scalar_prod(logn.gradient_spher(),
					logn.gradient_spher() ) ; 
    }
    
    source.set_std_base() ; 
    
    double int_grav = source().integrale() ; 
    
    // Matter term
    // -----------
    
    if (relativistic) {    
      source  = qpig * a_car * bbb * s_euler ;
    }
    else{
      source = qpig * ( 3 * press + nbar * uuu * uuu ) ; 
    }
    
    source.set_std_base() ; 
    
    double int_mat = source().integrale() ; 
    
    // Virial error
    // ------------
    if (ost != 0x0) {
      *ost << "Etoile_rot::grv3 : gravitational term : " << int_grav 
	   << endl ;
      *ost << "Etoile_rot::grv3 : matter term :        " << int_mat 
	   << endl ;
    }
    
    p_grv3 = new double( (int_grav + int_mat) / int_mat ) ; 	 
    
  }
  
  return *p_grv3 ; 
  
}


			//----------------------------//
			//	     R_circ	      //
			//----------------------------//

double Etoile_rot::r_circ() const {

    if (p_r_circ == 0x0) {    // a new computation is required
	
	// Index of the point at phi=0, theta=pi/2 at the surface of the star:
	const Mg3d* mg = mp.get_mg() ; 


	int l_b = nzet - 1 ; 
	int i_b = mg->get_nr(l_b) - 1 ; 
	int j_b;
	if (mg->get_type_t() == SYM) {
	  j_b = mg->get_nt(l_b) - 1 ;
	}else{
	  j_b = (mg->get_nt(l_b) - 1)/2 ;
	}
	int k_b = 0 ; 
    
	p_r_circ = new double( bbb()(l_b, k_b, j_b, i_b) * ray_eq() ) ; 

    }
    
    return *p_r_circ ; 

}

			//----------------------------//
			//	  Surface area	      //
			//----------------------------//

  double Etoile_rot::area() const {

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
      Cmp aaaa = get_a_car()() ;
      Valeur a2 = aaaa.va ; a2.std_base_scal() ;
      Cmp bbbb = get_bbb()() ;
      Valeur b = bbbb.va ; b.std_base_scal() ;
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

  double Etoile_rot::mean_radius() const {

    return sqrt(area()/(4*M_PI)) ;

  }















			//----------------------------//
			//	   Flattening	      //
			//----------------------------//

double Etoile_rot::aplat() const {

    if (p_aplat == 0x0) {    // a new computation is required
	
	p_aplat = new double( ray_pole() / ray_eq() ) ; 	 

    }
    
    return *p_aplat ; 

}


			//----------------------------//
			//	     Z_eq_f	      //
			//----------------------------//

double Etoile_rot::z_eqf() const {

    if (p_z_eqf == 0x0) {    // a new computation is required
	
	// Index of the point at phi=0, theta=pi/2 at the surface of the star:
	const Mg3d* mg = mp.get_mg() ; 
	int l_b = nzet - 1 ; 
	int i_b = mg->get_nr(l_b) - 1 ; 
	int j_b;
	if (mg->get_type_t() == SYM) {
	  j_b = mg->get_nt(l_b) - 1 ; 
	}else{
	  j_b = (mg->get_nt(l_b) - 1)/2 ; 
	}
	int k_b = 0 ; 
    
	double u_eq = uuu()(l_b, k_b, j_b, i_b) ; 
	double n_eq = nnn()(l_b, k_b, j_b, i_b) ; 
	double nphi_eq = nphi()(l_b, k_b, j_b, i_b) ; 
	
	p_z_eqf = new double( sqrt((1.-u_eq)/(1.+u_eq)) 
				/ (n_eq + nphi_eq * r_circ() )
			      - 1. ) ;
    }
    
    return *p_z_eqf ; 

}
			//----------------------------//
			//	     Z_eq_b	      //
			//----------------------------//

double Etoile_rot::z_eqb() const {

    if (p_z_eqb == 0x0) {    // a new computation is required
	
	// Index of the point at phi=0, theta=pi/2 at the surface of the star:
	const Mg3d* mg = mp.get_mg() ; 
	int l_b = nzet - 1 ; 
	int i_b = mg->get_nr(l_b) - 1 ;
	int j_b;
	if (mg->get_type_t() == SYM) {
	  j_b = mg->get_nt(l_b) - 1 ;
	}else{
	  j_b = (mg->get_nt(l_b) - 1) / 2 ;
	}
	int k_b = 0 ; 
    
	double u_eq = uuu()(l_b, k_b, j_b, i_b) ; 
	double n_eq = nnn()(l_b, k_b, j_b, i_b) ; 
	double nphi_eq = nphi()(l_b, k_b, j_b, i_b) ; 
	
	p_z_eqb = new double(  sqrt((1.+u_eq)/(1.-u_eq)) 
				/ (n_eq - nphi_eq * r_circ() )
			      - 1. )  ;

    }
    
    return *p_z_eqb ; 

}


			//----------------------------//
			//	     Z_pole	      //
			//----------------------------//

double Etoile_rot::z_pole() const {

    if (p_z_pole == 0x0) {    // a new computation is required
	
	double n_pole = nnn().val_point(ray_pole(), 0., 0.) ; 
	
	p_z_pole = new double(  1. / n_pole - 1. ) ; 

    }
    
    return *p_z_pole ; 

}


			//----------------------------//
			//     Quadrupole moment      //
			//----------------------------//


double Etoile_rot::mom_quad() const {

  using namespace Unites ;

  if (p_mom_quad == 0x0) {    // a new computation is required
	
    p_mom_quad = new double( mom_quad_old() ) ;
    if (relativistic) {
      double b = mom_quad_Bo() / ( mass_g() * mass_g() ) ;
      *p_mom_quad -= 4./3. * ( 1./4. + b ) * pow(mass_g(), 3) * ggrav * ggrav  ;
    }
  }
    
  return *p_mom_quad ; 

}


double Etoile_rot::mom_quad_Bo() const {

  using namespace Unites ;

  if (p_mom_quad_Bo == 0x0) {    // a new computation is required
	
    Cmp dens(mp) ; 
   
   dens = press() ;
   dens = a_car() * bbb() * nnn() * dens ; 
   dens.mult_rsint() ;
   dens.std_base_scal() ; 
      
   p_mom_quad_Bo = new double( - 32. * dens.integrale() / qpig  ) ; 

  }
    
  return *p_mom_quad_Bo ; 

}



double Etoile_rot::mom_quad_old() const {

  using namespace Unites ;

  if (p_mom_quad_old == 0x0) {    // a new computation is required
	
    // Source for of the Poisson equation for nu
    // -----------------------------------------

    Tenseur source(mp) ; 
	
    if (relativistic) {
      Tenseur beta = log(bbb) ; 
      beta.set_std_base() ; 
      source =  qpig * a_car *( ener_euler + s_euler ) 
	+ ak_car - flat_scalar_prod(logn.gradient_spher(), 
				    logn.gradient_spher() + beta.gradient_spher()) ; 
    }
    else {
      source = qpig * nbar ; 
    }
    source.set_std_base() ; 	

    // Multiplication by -r^2 P_2(cos(theta))
    //  [cf Eq.(7) of Salgado et al. Astron. Astrophys. 291, 155 (1994) ]
    // ------------------------------------------------------------------
	
    // Multiplication by r^2 : 
    // ----------------------
    Cmp& csource = source.set() ; 
    csource.mult_r() ; 
    csource.mult_r() ; 
    if (csource.check_dzpuis(2)) {
      csource.inc2_dzpuis() ; 
    }
		
    // Muliplication by cos^2(theta) :
    // -----------------------------
    Cmp temp = csource ; 
	
    // What follows is valid only for a mapping of class Map_radial : 
    assert( dynamic_cast<const Map_radial*>(&mp) != 0x0 ) ; 
		
    if (temp.get_etat() == ETATQCQ) {
      Valeur& vtemp = temp.va ; 
      vtemp = vtemp.mult_ct() ;	// multiplication by cos(theta)
      vtemp = vtemp.mult_ct() ;	// multiplication by cos(theta)
    }
	
    // Muliplication by -P_2(cos(theta)) :
    // ----------------------------------
    source = 0.5 * source() - 1.5 * temp ; 
	
    // Final result
    // ------------

    p_mom_quad_old = new double(- source().integrale() / qpig ) ; 	 

  }
    
  return *p_mom_quad_old ; 

}




}
