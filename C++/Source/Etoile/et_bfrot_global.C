/*
 *   Copyright (c) 2001 Jerome Novak
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


char et_bfrot_global_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.13  2004/03/25 10:29:03  j_novak
 * All LORENE's units are now defined in the namespace Unites (in file unites.h).
 *
 * Revision 1.12  2003/11/20 14:01:26  r_prix
 * changed member names to better conform to Lorene coding standards:
 * J_euler -> j_euler, EpS_euler -> enerps_euler, Delta_car -> delta_car
 *
 * Revision 1.11  2003/11/18 18:38:11  r_prix
 * use of new member EpS_euler: matter sources in equilibrium() and global quantities
 * no longer distinguish Newtonian/relativistic, as all terms should have the right limit...
 *
 * Revision 1.10  2003/11/13 12:13:15  r_prix
 * *) removed all uses of Etoile-type specific quantities: u_euler, press
 *    and exclusively use ener_euler, s_euler, sphph_euler, J_euler instead
 * *) this also fixes a bug in previous expression for mass_g() in 2-fluid case
 *
 * NOTE: some current Newtonian expressions are still completely broken..
 *
 * Revision 1.9  2003/09/17 08:27:50  j_novak
 * New methods: mass_b1() and mass_b2().
 *
 * Revision 1.8  2003/02/07 10:47:43  j_novak
 * The possibility of having gamma5 xor gamma6 =0 has been introduced for
 * tests. The corresponding parameter files have been added.
 *
 * Revision 1.7  2002/10/16 14:36:36  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.6  2002/10/14 14:20:08  j_novak
 * Error corrected for angu_mom()
 *
 * Revision 1.5  2002/05/10 09:26:52  j_novak
 * Added new class Et_rot_mag for magnetized rotating neutron stars (under development)
 *
 * Revision 1.4  2002/04/05 09:09:36  j_novak
 * The inversion of the EOS for 2-fluids polytrope has been modified.
 * Some errors in the determination of the surface were corrected.
 *
 * Revision 1.3  2002/01/08 14:43:53  j_novak
 * better determination of surfaces for 2-fluid stars
 *
 * Revision 1.2  2002/01/03 15:30:28  j_novak
 * Some comments modified.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 1.4  2001/08/31  15:07:12  novak
 * Retour a la version 1.2, sans la routine prolonge_c1
 *
 * Revision 1.2  2001/08/28 15:32:17  novak
 * The surface is now defined by the baryonic density
 *
 * Revision 1.1  2001/06/22 15:39:46  novak
 * Initial revision
 *
 *
 * $Header$
 *
 */


// Headers C
#include <stdlib.h>
#include <math.h>

// Headers Lorene
#include "et_rot_bifluid.h"
#include "unites.h"

//--------------------------//
//	Baryon mass	    //
//--------------------------//

double Et_rot_bifluid::mass_b1() const {

  if (p_mass_b1 == 0x0) {    // a new computation is required
	
    assert(nbar.get_etat() == ETATQCQ);
    Cmp dens1 = a_car() * bbb() * gam_euler() * eos.get_m1()* nbar();
    dens1.std_base_scal() ; 

    p_mass_b1 = new double( dens1.integrale() ) ;

  }
    
  return *p_mass_b1 ; 

} 

double Et_rot_bifluid::mass_b2() const {

  if (p_mass_b2 == 0x0) {    // a new computation is required
      assert(nbar2.get_etat() == ETATQCQ);
	
      Cmp dens2 = a_car() * bbb() * gam_euler2() * eos.get_m2() * nbar2();
      dens2.std_base_scal() ; 

      p_mass_b2 = new double( dens2.integrale() ) ;

  }
    
  return *p_mass_b2 ; 

} 

double Et_rot_bifluid::mass_b() const {

  if (p_mass_b == 0x0) 
    p_mass_b = new double(mass_b1() + mass_b2() ) ;

  return *p_mass_b;

} 


//----------------------------//
//	Gravitational mass    //
//----------------------------//

double Et_rot_bifluid::mass_g() const {

  if (p_mass_g == 0x0) {    // a new computation is required
	
      Tenseur vtmp (j_euler);
      vtmp.change_triad( mp.get_bvect_spher() );  // change to spherical triad

      Tenseur tJphi (mp);
      tJphi = vtmp(2);  // get the phi tetrad-component: J^ph r sin(th)

      // relativistic AND Newtonian: enerps_euler has right limit and tnphi->0
      Tenseur source = nnn * enerps_euler + 2 * b_car * tnphi * tJphi;
      source = a_car * bbb * source ;
	  
      source.set_std_base() ;

      p_mass_g = new double( source().integrale() ) ;
	  
  }
    
  return *p_mass_g ; 

} 
		
//----------------------------//
//	Angular momentum      //
//----------------------------//

double Et_rot_bifluid::angu_mom() const {

  if (p_angu_mom == 0x0) {    // a new computation is required
	
    Cmp dens(mp) ;

    // this should work for both relativistic AND Newtonian, 
    // provided j_euler has the right limit...
    Tenseur tmp = j_euler;
    tmp.change_triad( mp.get_bvect_spher() ) ;  
    dens = tmp(2) ;
    dens.mult_rsint() ;

    dens = a_car() * b_car()* bbb() * dens ; 

    dens.std_base_scal() ; 
      
    p_angu_mom = new double( dens.integrale() ) ;
      
  }
    
  return *p_angu_mom ; 

}


//----------------------------//
//	     GRV2	      //
//----------------------------//

double Et_rot_bifluid::grv2() const {

  using namespace Unites ;

  if (p_grv2 == 0x0) {    // a new computation is required
	
    Tenseur sou_m(mp);
    // should work for both relativistic AND Newtonian, provided sphph_euler has right limit..
    sou_m =  2 * qpig * a_car * sphph_euler ;

    Tenseur sou_q =  1.5 * ak_car - flat_scalar_prod(logn.gradient_spher(),logn.gradient_spher() ) ;	

    p_grv2 = new double( double(1) - lambda_grv2(sou_m(), sou_q()) ) ; 	

  }
    
  return *p_grv2 ; 

}


//----------------------------//
//	     GRV3	      //
//----------------------------//

double Et_rot_bifluid::grv3(ostream* ost) const {

  using namespace Unites ;

  if (p_grv3 == 0x0) {    // a new computation is required

    Tenseur source(mp) ; 
	
    // Gravitational term [cf. Eq. (43) of Gourgoulhon & Bonazzola
    // ------------------	    Class. Quantum Grav. 11, 443 (1994)]

    if (relativistic) {
      Tenseur alpha = dzeta - logn ; 
      Tenseur beta = log( bbb ) ; 
      beta.set_std_base() ; 
	    
      source = 0.75 * ak_car - flat_scalar_prod(logn.gradient_spher(), logn.gradient_spher() )
	+ 0.5 * flat_scalar_prod(alpha.gradient_spher(), beta.gradient_spher() ) ; 
	    
      Cmp aa = alpha() - 0.5 * beta() ; 
      Cmp daadt = aa.srdsdt() ;	// 1/r d/dth
	    
      // What follows is valid only for a mapping of class Map_radial : 
      const Map_radial* mpr = dynamic_cast<const Map_radial*>(&mp) ; 
      if (mpr == 0x0) {
	cout << "Et_rot_bifluid::grv3: the mapping does not belong"
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
      source = - 0.5 * flat_scalar_prod(logn.gradient_spher(),logn.gradient_spher() ) ; 
    }
	
    source.set_std_base() ; 

    double int_grav = source().integrale() ; 

    // Matter term
    // -----------

    // should work for relativistic AND Newtonian, provided s_euler has the right limit...
    source  = qpig * a_car * bbb * s_euler ;

    source.set_std_base() ; 

    double int_mat = source().integrale() ; 

    // Virial error
    // ------------
    if (ost != 0x0) {
      *ost << "Et_rot_bifluid::grv3 : gravitational term : " << int_grav 
	   << endl ;
      *ost << "Et_rot_bifluid::grv3 : matter term :        " << int_mat 
	   << endl ;
    }

    p_grv3 = new double( (int_grav + int_mat) / int_mat ) ; 	 

  }
    
  return *p_grv3 ; 

}


//----------------------------//
//	     R_circ	      //
//----------------------------//

double Et_rot_bifluid::r_circ2() const {

  if (p_r_circ2 == 0x0) {    // a new computation is required
	
    // Index of the point at phi=0, theta=pi/2 at the surface of the star:
    const Mg3d* mg = mp.get_mg() ; 
    assert(mg->get_type_t() == SYM) ; 
    int l_b = nzet - 1 ; 
    int i_b = mg->get_nr(l_b) - 1 ; 
    int j_b = mg->get_nt(l_b) - 1 ; 
    int k_b = 0 ; 
    
    p_r_circ2 = new double( bbb()(l_b, k_b, j_b, i_b) * ray_eq2() ) ; 

  }
    
  return *p_r_circ2 ; 

}


//----------------------------//
//	   Flattening	      //
//----------------------------//

double Et_rot_bifluid::aplat2() const {

  if (p_aplat2 == 0x0) {    // a new computation is required
	
    p_aplat2 = new double( ray_pole2() / ray_eq2() ) ; 	 

  }
    
  return *p_aplat2 ; 

}



//----------------------------//
//     Quadrupole moment      //
//----------------------------//

double Et_rot_bifluid::mom_quad() const {

  using namespace Unites ;

  if (p_mom_quad == 0x0) {    // a new computation is required
	
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
      source = qpig * (nbar + nbar2); 
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

    p_mom_quad = new double( source().integrale() / qpig ) ; 	 

  }
    
  return *p_mom_quad ; 

}


//--------------------------//
//	Stellar surface	    //
//--------------------------//

const Itbl& Et_rot_bifluid::l_surf() const {

  if (p_l_surf == 0x0) {    // a new computation is required
    
    assert(p_xi_surf == 0x0) ;  // consistency check
	
    int np = mp.get_mg()->get_np(0) ;   
    int nt = mp.get_mg()->get_nt(0) ;   
	
    p_l_surf = new Itbl(np, nt) ;
    p_xi_surf = new Tbl(np, nt) ;
	
    double nb0 = 0 ;	// definition of the surface
    double precis = 1.e-15 ; 
    int nitermax = 100 ; 
    int niter ; 
	      
    // Cmp defining the surface of the star (via the density fields)
    // 
    Cmp surf(mp) ;
    surf = -0.2*nbar()(0,0,0,0) ;
    surf.annule(0, nzet-1) ;
    surf += nbar() ; ;
    surf = prolonge_c1(surf, nzet) ;
    
    (surf.va).equipot(nb0, nzet, precis, nitermax, niter, *p_l_surf, 
			*p_xi_surf) ; 
    
  }
   
  return *p_l_surf ; 
    
}
const Itbl& Et_rot_bifluid::l_surf2() const {

  if (p_l_surf2 == 0x0) {    // a new computation is required
    
    assert(p_xi_surf2 == 0x0) ;  // consistency check
	
    int np = mp.get_mg()->get_np(0) ;   
    int nt = mp.get_mg()->get_nt(0) ;   
	
    p_l_surf2 = new Itbl(np, nt) ;
    p_xi_surf2 = new Tbl(np, nt) ;
	
    double nb0 = 0 ;	// definition of the surface
    double precis = 1.e-15 ; 
    int nitermax = 100 ; 
    int niter ; 
	
    // Cmp defining the surface of the star (via the density fields)
    // 
    Cmp surf2(mp) ;
    surf2 = -0.2*nbar2()(0,0,0,0) ;
    surf2.annule(0, nzet-1) ;
    surf2 += nbar2() ; ;
    surf2 = prolonge_c1(surf2, nzet) ;
    
    (surf2.va).equipot(nb0, nzet, precis, nitermax, niter, *p_l_surf2, 
			*p_xi_surf2) ; 
    
  }
   
  return *p_l_surf2 ; 
    
}

const Tbl& Et_rot_bifluid::xi_surf2() const {

  if (p_xi_surf2 == 0x0) {    // a new computation is required
    
    assert(p_l_surf2 == 0x0) ;  // consistency check
	
    l_surf2() ;  // the computation is done by l_surf2()
    
  }
   
  return *p_xi_surf2 ; 
    
}

//--------------------------//
//	Coordinate radii    //
//--------------------------//

double Et_rot_bifluid::ray_eq2() const {

  if (p_ray_eq2 == 0x0) {    // a new computation is required
	
    const Mg3d& mg = *(mp.get_mg()) ;
	
    int type_t = mg.get_type_t() ; 
    int nt = mg.get_nt(0) ; 	
	
    if ( type_t == SYM ) {
      assert( ( mg.get_type_p() == SYM) || (mg.get_type_p() == NONSYM) ) ; 
      int k = 0 ; 
      int j = nt-1 ; 
      int l = l_surf2()(k, j) ; 
      double xi = xi_surf2()(k, j) ; 
      double theta = M_PI / 2 ; 
      double phi = 0 ; 
	    
      p_ray_eq2 = new double( mp.val_r(l, xi, theta, phi) ) ;

    }
    else {
      cout << "Et_rot_bifluid::ray_eq2 : the case type_t = " << type_t
	   << " is not contemplated yet !" << endl ;
      abort() ; 
    }

  }
    
  return *p_ray_eq2 ; 

} 


double Et_rot_bifluid::ray_eq2_pis2() const {

  if (p_ray_eq2_pis2 == 0x0) {    // a new computation is required
	
    const Mg3d& mg = *(mp.get_mg()) ;
	
    int type_t = mg.get_type_t() ; 
    int type_p = mg.get_type_p() ; 
    int nt = mg.get_nt(0) ; 	
    int np = mg.get_np(0) ; 	
    
    if ( type_t == SYM ) {
      
      int j = nt-1 ; 
      double theta = M_PI / 2 ; 
      double phi = M_PI / 2 ;
      
      switch (type_p) {
	
      case SYM : {
	int k = np / 2  ; 
	int l = l_surf2()(k, j) ; 
	double xi = xi_surf2()(k, j) ; 
	p_ray_eq2_pis2 = new double( mp.val_r(l, xi, theta, phi) ) ;
	break ; 
      }
      
      case NONSYM : {
	assert( np % 4 == 0 ) ; 
	int k = np / 4  ; 
	int l = l_surf2()(k, j) ; 
	double xi = xi_surf2()(k, j) ; 
	p_ray_eq2_pis2 = new double( mp.val_r(l, xi, theta, phi) ) ;
	break ; 
      }
	    
      default : {
	cout << "Et_rot_bifluid::ray_eq2_pis2 : the case type_p = " 
	     << type_p << " is not contemplated yet !" << endl ;
	abort() ; 
      }
      } 
      
    }
    else {
      cout << "Et_rot_bifluid::ray_eq2_pis2 : the case type_t = " << type_t
	   << " is not contemplated yet !" << endl ;
      abort() ; 
    }
    
  }
    
  return *p_ray_eq2_pis2 ; 

} 


double Et_rot_bifluid::ray_eq2_pi() const {

  if (p_ray_eq2_pi == 0x0) {    // a new computation is required
    
    const Mg3d& mg = *(mp.get_mg()) ;
    
    int type_t = mg.get_type_t() ; 
    int type_p = mg.get_type_p() ; 
    int nt = mg.get_nt(0) ; 	
    int np = mg.get_np(0) ; 	
    
    if ( type_t == SYM ) {
      
      switch (type_p) {
	
      case SYM : {
	p_ray_eq2_pi = new double( ray_eq2() ) ;
	break ; 
      }		
      
      case NONSYM : {
	int k = np / 2  ; 
	int j = nt-1 ; 
	int l = l_surf2()(k, j) ; 
	double xi = xi_surf2()(k, j) ; 
	double theta = M_PI / 2 ; 
	double phi = M_PI ; 
	
	p_ray_eq2_pi = new double( mp.val_r(l, xi, theta, phi) ) ;
	break ;
      }
      
      default : {
	
	cout << "Et_rot_bifluid::ray_eq2_pi : the case type_t = " << type_t
	     << " and type_p = " << type_p << endl ; 
	cout << " is not contemplated yet !" << endl ;
	abort() ; 
	break ; 
      }
      }
    }
    
  }
  
  return *p_ray_eq2_pi ; 
  
} 

double Et_rot_bifluid::ray_pole2() const {

  if (p_ray_pole2 == 0x0) {    // a new computation is required
	
    assert( ((mp.get_mg())->get_type_t() == SYM) 
	    || ((mp.get_mg())->get_type_t() == NONSYM) ) ;  
	
    int k = 0 ; 
    int j = 0 ; 
    int l = l_surf2()(k, j) ; 
    double xi = xi_surf2()(k, j) ; 
    double theta = 0 ; 
    double phi = 0 ; 
	    
    p_ray_pole2 = new double( mp.val_r(l, xi, theta, phi) ) ;

  }
    
  return *p_ray_pole2 ; 

} 



