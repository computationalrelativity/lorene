/*
 *  Methods for computing global quantities within the class Star_rot_CFC
 *
 *    (see file star_rot_cfc.h for documentation).
 *
 */

/*
 *   Copyright (c) 2024 Jerome Novak
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
#include "proto.h"

        //-------------------------------------------------------//
        //                  Baryon mass                          //
        //                                                       //
        //  Note: In Lorene units, mean particle mass is unity   //
        //-------------------------------------------------------//


namespace Lorene {
double Star_rot_CFC::mass_b() const {

  if (p_mass_b == 0x0) {    // a new computation is required

    Scalar dens = psi4*psi2 * gam_euler * nbar ;

    dens.std_spectral_base() ;

    p_mass_b = new double( dens.integrale() ) ;

  }

  return *p_mass_b ;

}

          //---------------------------------------------------------//
          //            Gravitational mass                           //   
          //                                                         //
          // Note: This is the Komar mass for stationary and         //
          //       asymptotically flat spacetime (see, eg, Wald)     //
          //---------------------------------------------------------//

double Star_rot_CFC::mass_g() const {

  if (p_mass_g == 0x0) {    // a new computation is required

    Scalar j_source = 2.* contract(contract(gamma.cov(), 0, j_euler, 0), 
				   0, beta, 0) ;
           
    Scalar dens = psi4*psi2 * ( nn * (ener_euler + s_euler) - j_source ) ;

    dens.std_spectral_base() ;

    p_mass_g = new double( dens.integrale() ) ;

  }

  return *p_mass_g ;

}

                //--------------------------------------//
                //    Angular momentum                  //
                //                                      // 
                // Komar-type integral (see, eg, Wald)  //
                //--------------------------------------// 

double Star_rot_CFC::angu_mom() const {

  if (p_angu_mom == 0x0) {    // a new computation is required

    // phi_kill = axial killing vector 

    Vector phi_kill(mp, CON, mp.get_bvect_spher()) ;

    phi_kill.set(1).set_etat_zero() ;
    phi_kill.set(2).set_etat_zero() ;
    phi_kill.set(3) = 1. ;
    phi_kill.set(3).std_spectral_base() ;
    phi_kill.set(3).mult_rsint() ;

    Scalar j_source = contract(contract(gamma.cov(), 0, j_euler, 0),
			       0, phi_kill, 0) ;

    Scalar dens = psi4*psi2 * j_source  ;

    dens.std_spectral_base() ;

    p_angu_mom = new double( dens.integrale() ) ;


  }

  return *p_angu_mom ;

}

                     //---------------------//
                     //        T/W          //
                     //---------------------//

double Star_rot_CFC::tsw() const {

  if (p_tsw == 0x0) {    // a new computation is required

    double tcin = 0.5 * omega * angu_mom() ;

    Scalar dens = psi4*psi2 * gam_euler * ener ;
    
    dens.std_spectral_base() ;
    
    double mass_p = dens.integrale() ;

    p_tsw = new double( tcin / ( mass_p + tcin - mass_g() ) ) ;

  }

  return *p_tsw ;

}


      //--------------------------------------------------------------//
      //                        GRV2                                  //   
      // cf. Eq. (28) of Bonazzola & Gourgoulhon CQG, 11, 1775 (1994) // 
      //                                                              //
      //--------------------------------------------------------------//

double Star_rot_CFC::grv2() const {

  using namespace Unites ;

  if (p_grv2 == 0x0) {    // a new computation is required

      // determinant of the 2-metric k_{ab}
      
      Scalar k_det = gamma.cov()(1,1)*gamma.cov()(2,2) - 
	  gamma.cov()(1,2)*gamma.cov()(1,2) ;
      
      
      //**
      // sou_m = 8\pi T_{\mu\nu} m^{\mu}m^{\nu}
      // => sou_m = 8\pi [ (E+P) U^2 + P ], where v^2 = v_i v^i 
      //
      //**
      
      Scalar sou_m = 2 * qpig * ( (ener_euler + press)*v2 + press ) ;
      
      sou_m = sqrt( k_det )*sou_m ;
      
      sou_m.std_spectral_base() ;
      
      
      // This is the term 3k_a k^a.
      
      Scalar sou_q = 3 *( hatA(1,3) * hatA(1,3) + hatA(2,3)*hatA(2,3) )
	/ (psi4*psi4*psi4)  ;
      
      
      // This is the term \nu_{|| a}\nu^{|| a}. 
      //
      
      Scalar sou_tmp = gamma.con()(1,1) * logn.dsdr() * logn.dsdr() ;
      
      Scalar term_2 = 2 * gamma.con()(1,2) * logn.dsdr() * logn.dsdt() ;
      
      term_2.div_r_dzpuis(4) ;
      
      Scalar term_3 = gamma.con()(2,2) * logn.dsdt() * logn.dsdt() ;
      
      term_3.div_r_dzpuis(2) ;
      term_3.div_r_dzpuis(4) ;
      
      sou_tmp += term_2 + term_3 ;
      
      
      // Source of the gravitational part
      
      sou_q -= sou_tmp ;
      
      sou_q = sqrt( k_det )*sou_q ;
      
      sou_q.std_spectral_base() ;
      
      p_grv2 = new double( double(1) + integrale2d(sou_m)/integrale2d(sou_q) ) ;

  }

  return *p_grv2 ;

}


     //-------------------------------------------------------------//
     //                 GRV3                                        //
     // cf. Eq. (29) of Gourgoulhon & Bonazzola CQG, 11, 443 (1994) //
     //-------------------------------------------------------------//

double Star_rot_CFC::grv3() const {

  using namespace Unites ;

  if (p_grv3 == 0x0) {    // a new computation is required

    // Gravitational term 
    // -------------------

    Scalar psi6 = psi4*psi2 ;
    Scalar sou_q = 0.75*hatA_quad/psi6
      - psi6*contract(logn.derive_cov(gamma), 0,logn.derive_con(gamma), 0) ;


    Tensor t_tmp = contract(gamma.connect().get_delta(), 2, 
     			    gamma.connect().get_delta(), 0) ;

    Scalar tmp_1 = 0.25* contract( gamma.con(), 0, 1, 
				   contract(t_tmp, 0, 3), 0, 1 ) ;

    Scalar tmp_2 = 0.25* contract( gamma.con(), 0, 1, 
		        contract( contract( gamma.connect().get_delta(), 0, 1), 
                        0, gamma.connect().get_delta(), 0), 0, 1)  ;

    sou_q = sou_q + psi6*(tmp_1 - tmp_2) ;

    // sou_q = psi4*psi2 * sou_q ; 

    sou_q.std_spectral_base() ;

    double int_grav = sou_q.integrale() ;


    // Matter term 
    // --------------

    Scalar sou_m = qpig*s_euler ;

    sou_m = psi6 * sou_m ;

    sou_m.std_spectral_base() ;

    double int_mat = sou_m.integrale() ;

    p_grv3 = new double( (int_grav + int_mat) / int_mat ) ;
    


  }

  return *p_grv3 ;

}


                //--------------------//
                //     R_circ         //
                //--------------------//

double Star_rot_CFC::r_circ() const {

  if (p_r_circ ==0x0) {  // a new computation is required

    // Index of the point at phi=0, theta=pi/2 at the surface of the star:
    const Mg3d* mg = mp.get_mg() ;
    int l_b = nzet - 1 ; 
    int i_b = mg->get_nr(l_b) - 1 ; 
    int j_b = (mg->get_type_t() == SYM ? mg->get_nt(l_b) - 1 : mg->get_nt(l_b) / 2) ; 
    int k_b = 0 ;

    double gamma_phi = gamma.cov()(3,3).val_grid_point(l_b, k_b, j_b, i_b) ;

    p_r_circ = new double( sqrt( gamma_phi ) * ray_eq() ) ;

  }

  return *p_r_circ ;

}


                //--------------------//
                //     R_circ         //
                //--------------------//

double Star_rot_CFC::rp_circ() const {

    if (p_rp_circ ==0x0) {  // a new computation is required
	const Mg3d& mg = *mp.get_mg() ;
	int nz = mg.get_nzone() ;
	assert(nz>2) ;
	int np = mg.get_np(0) ;
	if (np != 1) {
	    cout << "The polar circumferential radius is only well defined\n"
		 << "with np = 1!" << endl ;
	    abort() ;
	}
	int nt = mg.get_nt(0) ;
	Sym_tensor gam = gamma.cov() ;
	const Coord& tet = mp.tet ;
	Scalar rrr(mp) ;
	rrr.annule_hard() ;
	rrr.annule(nzet,nz-1) ;
	double phi = 0 ;
	for (int j=0; j<nt; j++) {
	    double theta = (+tet)(0, 0, j, 0) ;
	    double erre = 
		mp.val_r(l_surf()(0,j), xi_surf()(0, j), theta, phi) ;
	    for (int lz=0; lz<nzet; lz++) {
		int nrz = mg.get_nr(lz) ;
		for (int i=0; i<nrz; i++) {
		    rrr.set_grid_point(lz, 0, j, i) = erre ; 
		}
	    }
	}
	rrr.std_spectral_base() ;
	Scalar drrr = rrr.dsdt() ;

	Scalar fff(mp) ;
	fff.annule_hard() ;
	fff.annule(nzet,nz-1) ;
	for (int j=0; j<nt; j++) {
	    double theta = (+tet)(0, 0, j, 0) ;
	    int ls = l_surf()(0, j) ;
	    double xs = xi_surf()(0, j) ;
	    double grr = gam(1,1).get_spectral_va().val_point_jk(ls, xs, j, 0) ;
	    double grt = gam(1,2).get_spectral_va().val_point_jk(ls, xs, j, 0) ;
	    double gtt = gam(2,2).get_spectral_va().val_point_jk(ls, xs, j, 0) ;
	    double rr = mp.val_r(ls, xs, theta, phi) ;
	    double dr = drrr.get_spectral_va().val_point_jk(ls, xs, j, 0) ;
	    for (int i=0; i< mg.get_nr(nzet-1); i++) {
		fff.set_grid_point(nzet-1, 0, j, i) 
		    = sqrt(grr*dr*dr + 2*grt*rr*dr + gtt*rr*rr) ;
	    }
	}
	fff.std_spectral_base() ;
	fff.set_spectral_va().coef() ;
	p_rp_circ = new double((*fff.get_spectral_va().c_cf)(nzet-1, 0, 0, 0)) ;      
    }
    return *p_rp_circ ;
}

                //--------------------------//
                //       Flattening         //
                //--------------------------//

double Star_rot_CFC::aplat() const {

  return ray_pole() / ray_eq() ;

}

                //--------------------------//
                //       Ellipticity        //
                //--------------------------//

double Star_rot_CFC::ellipt() const {

  return sqrt(1. - (rp_circ()*rp_circ())/(r_circ()*r_circ())) ;

}
}
