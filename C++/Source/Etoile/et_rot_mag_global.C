/*
 * Methods for computing global quantities within the class Etoile_rot
 *
 * (see file etoile.h for documentation)
 */

/*
 *   Copyright (c) 2000-2001 Eric Gourgoulhon
 *   Copyright (c) 2002 Emmanuel Marcq
 *   Copyright (c) 2002 Jerome Novak
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


char et_rot_mag_global_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.15  2002/10/17 11:30:54  j_novak
 * Corrected mistake in angu_mom()
 *
 * Revision 1.14  2002/06/03 13:00:45  e_marcq
 *
 * conduc parameter read in parmag.d
 *
 * Revision 1.12  2002/05/22 12:20:17  j_novak
 * *** empty log message ***
 *
 * Revision 1.11  2002/05/20 15:44:55  e_marcq
 *
 * Dimension errors corrected, parmag.d input file created and read
 *
 * Revision 1.10  2002/05/20 10:31:59  j_novak
 * *** empty log message ***
 *
 * Revision 1.9  2002/05/20 08:27:59  j_novak
 * *** empty log message ***
 *
 * Revision 1.8  2002/05/17 15:08:01  e_marcq
 *
 * Rotation progressive plug-in, units corrected, Q and a_j new member data
 *
 * Revision 1.7  2002/05/16 13:27:11  j_novak
 * *** empty log message ***
 *
 * Revision 1.6  2002/05/16 10:02:09  j_novak
 * Errors in stress energy tensor corrected
 *
 * Revision 1.5  2002/05/15 09:53:59  j_novak
 * First operational version
 *
 * Revision 1.4  2002/05/14 13:45:30  e_marcq
 *
 * Correction de la formule du rapport gyromagnetique
 *
 * Revision 1.1  2002/05/10 09:26:52  j_novak
 * Added new class Et_rot_mag for magnetized rotating neutron stars (under development)
 *
 *
 * $Header$
 *
 */

// Headers C
#include <stdlib.h>
#include <math.h>

// Headers Lorene
#include "et_rot_mag.h"
#include "unites_mag.h"

// Definition des fonctions membres differentes ou nouvelles

void Et_rot_mag::MHD_comput() {
  // Calcule les grandeurs du tenseur impulsion-energie EM a partir des champs
  #include"unites_mag.h"
  if (mu0<0.) { //to avoid compilation warnings
    cout << qpig << f_unit << msol << km << mevpfm3 << endl ; 
    cout << mag_unit << elec_unit << endl ;
  }    
  
  Tenseur ATTENS(A_t) ;
  Tenseur APTENS(A_phi) ;
  
  Tenseur ApAp ( flat_scalar_prod_desal(APTENS.gradient_spher(),
					APTENS.gradient_spher())() );
  Tenseur ApAt ( flat_scalar_prod_desal(APTENS.gradient_spher(),
					ATTENS.gradient_spher())() );
  Tenseur AtAt ( flat_scalar_prod_desal(ATTENS.gradient_spher(),
					ATTENS.gradient_spher())() );

  if (ApAp.get_etat() != ETATZERO) {
    ApAp.set().div_rsint() ;
    ApAp.set().div_rsint() ;
  }
  if (ApAt.get_etat() != ETATZERO) 
    ApAt.set().div_rsint() ;
  
  E_em = 0.5*mu0 * ( 1/(a_car*nnn*nnn) * (AtAt + 2*tnphi*ApAt)
	      + ( (tnphi*tnphi/(a_car*nnn*nnn)) + 1/(a_car*b_car) )*ApAp );
  Jp_em = -mu0 * (ApAt + tnphi*ApAp) /(a_car*nnn) ;
  if (Jp_em.get_etat() != ETATZERO) Jp_em.set().mult_rsint() ;
  Srr_em = 0 ;
  // Stt_em = -Srr_em
  Spp_em = E_em ;
}

Tenseur Et_rot_mag::Elec() const {

  Cmp E_r(mp); Cmp E_t(mp);
  E_r = 1/(sqrt(a_car())*nnn())*(A_t.dsdr()+nphi()*A_phi.dsdr()) ;
  E_t = 1/(sqrt(a_car())*nnn())*(A_t.srdsdt()+nphi()*A_phi.srdsdt()) ;
  E_r.va.set_base((A_t.dsdr()).va.base) ;
  E_t.va.set_base((A_t.srdsdt()).va.base) ;
  Tenseur Elect(mp, 1, CON, mp.get_bvect_spher()) ;
  Elect.set_etat_qcq() ;
  Elect.set(0) = E_r ;
  Elect.set(1) = E_t ;
  Elect.set(2) = 0. ;
  
    return Elect*elec_unit ;

}

Tenseur Et_rot_mag::Magn() const {

  Cmp B_r(mp); Cmp B_t(mp);
  B_r = 1/(sqrt(a_car())*bbb())*A_phi.srdsdt();
  B_r.va.set_base((A_phi.srdsdt()).va.base) ;
  B_r.div_rsint();
  B_t = 1/(sqrt(a_car())*bbb())*A_phi.dsdr();
  B_t.va.set_base((A_phi.dsdr()).va.base) ;
  B_t.div_rsint();

  Tenseur Bmag(mp, 1, CON, mp.get_bvect_spher()) ;
  Bmag.set_etat_qcq() ;
  Bmag.set(0) = B_r ;
  Bmag.set(1) = B_t ;
  Bmag.set(2) = 0. ;

    return Bmag*mag_unit ;

}

double Et_rot_mag::MagMom() const {

  int Z = mp.get_mg()->get_nzone();   
  double mm ;

  if(A_phi.get_etat()==ETATZERO) {

    mm = 0 ;
  }else{

  Valeur** asymp = A_phi.asymptot(1) ;
  mm = 4*M_PI*(*asymp[1])(Z-1,0,mp.get_mg()->get_nt(Z-1)-1,0) ;

  delete asymp[0] ;
  delete asymp[1] ;

  delete [] asymp ;
  }

  return mm*mag_unit*r_unit*r_unit*r_unit/mu_si*1.e9 ;

}

double Et_rot_mag::Q_comput() const {

  int Z = mp.get_mg()->get_nzone();

  if(A_t.get_etat()==ETATZERO) {
    return 0 ;
  }else{
  Valeur** asymp = A_t.asymptot(1) ;

  double Q_c = -4*M_PI*(*asymp[1])(Z-1,0,0,0) ;
  delete asymp[0] ;
  delete asymp[1] ;

  delete [] asymp ;

  return Q_c *(j_unit/v_unit*pow(r_unit,3)) ;}
  }

double Et_rot_mag::Q_int() const {

  double Qi = 0. ;

	if (relativistic) {

	    Cmp dens = a_car() * bbb() * nnn() * j_t ;
	    
	    dens.std_base_scal() ; 

	    Qi = dens.integrale() ;


	}
	else{  // Newtonian case 
	    assert(nbar.get_etat() == ETATQCQ) ; 

	    Qi = ( j_t.integrale() ) ;

	}

    
    
    return Qi * (j_unit/v_unit*pow(r_unit,3)) ; 


}

double Et_rot_mag::GyroMag() const {
  return 2*MagMom()*mass_g()/(Q_comput()*angu_mom()*v_unit*r_unit); 
}
			//----------------------------//
			//	Gravitational mass    //
			//----------------------------//

double Et_rot_mag::mass_g() const {

    if (p_mass_g == 0x0) {    // a new computation is required
	
	if (relativistic) {

	    Tenseur source = nnn * (ener_euler + E_em + s_euler + Spp_em) + 
	      nphi * Jp_em + 2 * bbb * (ener_euler + press) * tnphi * uuu ;
 
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

double Et_rot_mag::angu_mom() const {

    if (p_angu_mom == 0x0) {    // a new computation is required
	
	Cmp dens = uuu() ; 

	dens.mult_r() ;			//  Multiplication by
	dens.va = (dens.va).mult_st() ;	//    r sin(theta)

	if (relativistic) {
	    dens = a_car() * (b_car() * (ener_euler() + press()) 
			* dens + bbb() * Jp_em()) ; 
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

// Redefini en virtual dans le .h : A CHANGER

double Et_rot_mag::tsw() const {

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

double Et_rot_mag::grv2() const {

    if (p_grv2 == 0x0) {    // a new computation is required
	
		// To get qpig:	
		#include "unites.h"	
		// To avoid some compilation warnings
		if (p_grv2 != 0x0) {
	    	cout << f_unit << msol << km << mevpfm3 << endl ;
		}

        Tenseur sou_m =  2 * qpig * a_car * (press + (ener_euler+press)
        						* uuu*uuu ) ;
        						
        Tenseur sou_q =   2 * qpig * a_car * Spp_em + 1.5 * ak_car
        				 - flat_scalar_prod(logn.gradient_spher(),
						     				logn.gradient_spher() ) ;	

		p_grv2 = new double( double(1) - lambda_grv2(sou_m(), sou_q()) ) ; 	

    }
    
    return *p_grv2 ; 

}


			//----------------------------//
			//	     GRV3	      //
			//----------------------------//

double Et_rot_mag::grv3(ostream* ost) const {

    if (p_grv3 == 0x0) {    // a new computation is required

	// To get qpig:	
	#include "unites.h"	    
	// To avoid some compilation warnings
	if (p_grv3 != 0x0) {
	    cout << f_unit << msol << km << mevpfm3 << endl ; 
	}    


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
	    source  = qpig * a_car * bbb * ( s_euler + Spp_em ) ;
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
			//     Quadrupole moment      //
			//----------------------------//

double Et_rot_mag::mom_quad() const {

    if (p_mom_quad == 0x0) {    // a new computation is required
	
	// To get qpig:	
	#include "unites.h"	    
	// To avoid some compilation warnings
	if (p_mom_quad != 0x0) {
	    cout << f_unit << msol << km << mevpfm3 << endl ; 
	}    

	// Source for of the Poisson equation for nu
	// -----------------------------------------

	Tenseur source(mp) ; 
	
	if (relativistic) {
	    Tenseur beta = log(bbb) ; 
	    beta.set_std_base() ; 
	    source =  qpig * a_car *( ener_euler + s_euler + Spp_em ) 
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

	p_mom_quad = new double( source().integrale() / qpig ) ; 	 

    }
    
    return *p_mom_quad ; 

}
