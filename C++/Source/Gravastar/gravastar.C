/*
 *  Methods of class Gravastar
 *
 *    (see file gravastar.h for documentation).
 *
 */

/*
 *   Copyright (c) 2010 Frederic Vincent
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

 


// C++ headers
//#include <>

// C headers
#include <cmath>

// Lorene headers
#include "gravastar.h"
#include "eos.h"
#include "utilitaires.h"
#include "param.h"
#include "unites.h"
#include "nbr_spx.h"


                    //---------------------//
                    //     Constructor     //
                    //---------------------//

namespace Lorene {
Gravastar::Gravastar(Map& mpi, int nzet_i, const Eos& eos_i, const double rho_core_i) 
  : Star_rot(mpi,nzet_i,true,eos_i), rho_core(rho_core_i)
{} 

                            //------------//
			    // Destructor //
			    //------------//

Gravastar::~Gravastar(){
    del_deriv() ; 
}

                    //----------------------//
                    //     Eq. of state     //
                    //----------------------//

void Gravastar::equation_of_state() {

	Scalar ent_eos = ent ;


    // Slight rescale of the enthalpy field in case of 2 domains inside the
    //  star


        double epsilon = 1.e-12 ;

	const Mg3d* mg = mp.get_mg() ;
        int nz = mg->get_nzone() ;

        Mtbl xi(mg) ;
        xi.set_etat_qcq() ;
        for (int l=0; l<nz; l++) {
        	xi.t[l]->set_etat_qcq() ;
        	for (int k=0; k<mg->get_np(l); k++) {
        		for (int j=0; j<mg->get_nt(l); j++) {
        			for (int i=0; i<mg->get_nr(l); i++) {
        				xi.set(l,k,j,i) =
        					mg->get_grille3d(l)->x[i] ;
        			}
        		}
        	}

        }

     	Scalar fact_ent(mp) ;
     	fact_ent.allocate_all() ;
     	
     	fact_ent.set_domain(0) = 1 + epsilon * xi(0) * xi(0) ;
     	fact_ent.set_domain(1) = 1 - 0.25 * epsilon * (xi(1) - 1) * (xi(1) - 1) ;
     	
     	for (int l=nzet; l<nz; l++) {
     		fact_ent.set_domain(l) = 1 ;
     	}

    if (nzet > 1) {

    	if (nzet > 2) {
    	
    		cout << "Gravastar::equation_of_state: not ready yet for nzet > 2 !"
    		     << endl ;    	
    	}

    	ent_eos = fact_ent * ent_eos ;
    	ent_eos.std_spectral_base() ;
    }



    /*
      Call to the Eos
      There is no Eos inside the core where p=-rho=cst. Thus the most internal domain in treated separately.
      Inside the crust, the Eos is called domain by domain.
     */

    Scalar tempo(mp) ;

    //nbar.set_etat_qcq() ;
    nbar = 0 ;
    for (int l=1; l<nzet; l++) {

        Param par ;       // Paramater for multi-domain equation of state
        par.add_int(l) ;

	

        tempo =  eos.nbar_ent(ent_eos, 1, l, &par) ;
	//cout << "tempo eos= " << tempo << endl;
        nbar = nbar + tempo ;

    }

    //ener.set_etat_qcq() ;
    ener = rho_core ;
    ener.annule(1,nz-1);
    for (int l=1; l<nzet; l++) {

        Param par ;    // Paramater for multi-domain equation of state
        par.add_int(l) ;

        tempo =  eos.ener_ent(ent_eos, 1, l, &par) ;

        ener = ener + tempo ;

    }

    //press.set_etat_qcq() ;
    press = -rho_core ;
    press.annule(1,nz-1);
    for (int l=1; l<nzet; l++) {

        Param par ;     // Paramater for multi-domain equation of state
        par.add_int(l) ;

        tempo =  eos.press_ent(ent_eos, 1, l, &par) ;

        press = press + tempo ;

    }


    // Set the bases for spectral expansion
    nbar.std_spectral_base() ; 
    ener.std_spectral_base() ; 
    press.std_spectral_base() ; 

    // The derived quantities are obsolete
    del_deriv() ; 
    
}

// Printing
// --------

ostream& Gravastar::operator>>(ostream& ost) const {

  using namespace Unites ;
  
  ost << endl ; 
  
  ost << "Number of domains occupied by the star : " << nzet << endl ; 
  
  ost << "Equation of state : " << endl ; 
  ost << eos << endl ; 
  
  const Mg3d* mg = mp.get_mg() ; 
  int l_cr = 0, i_cr = mg->get_nr(l_cr) - 1, j_cr = mg->get_nt(l_cr) - 1, k_cr = 0 ; 

  ost << endl << "Inner crust enthalpy : " << ent.val_grid_point(l_cr, k_cr, j_cr, i_cr) << " c^2" << endl ; 
  
  ost << "Central proper energy density : " << ener.val_grid_point(0,0,0,0) 
      << " rho_nuc c^2" << endl ; 
  ost << "Central pressure : " << press.val_grid_point(0,0,0,0) 
      << " rho_nuc c^2" << endl ; 
  
  ost << endl ;
  ost << "Central lapse N :      " << nn.val_grid_point(0,0,0,0) <<  endl ; 
  //    ost << "Central value of lnq : " << lnq.val_grid_point(0,0,0,0) <<  endl ; 
  
  ost << endl 
      << "Coordinate equatorial radius (phi=0) a1 =    " 
      << ray_eq()/km << " km" << endl ;  
  ost << "Coordinate equatorial radius (phi=pi/2) a2 = " 
      << ray_eq_pis2()/km << " km" << endl ;  
  ost << "Coordinate equatorial radius (phi=pi):       " 
      << ray_eq_pi()/km << " km" << endl ;  
  ost << "Coordinate polar radius a3 =                 " 
      << ray_pole()/km << " km" << endl ;  
  ost << "Axis ratio a2/a1 = " << ray_eq_pis2() / ray_eq() 
      << "  a3/a1 = " << ray_pole() / ray_eq() << endl ; 	
  
  ///////////////////
  
  double omega_c = get_omega_c() ; 
  
  ost << endl ; 
  
  if (omega != __infinity) {
    ost << "Uniformly rotating star" << endl ; 
    ost << "-----------------------" << endl ; 
     
    double freq = omega / (2.*M_PI) ;  
    ost << "Omega : " << omega * f_unit 
	<< " rad/s     f : " << freq * f_unit << " Hz" << endl ; 
    ost << "Rotation period : " << 1000. / (freq * f_unit) << " ms"
	<< endl ;
    
  }
  else {
    ost << "Differentially rotating star" << endl ; 
    ost << "----------------------------" << endl ; 
    
    double freq = omega_c / (2.*M_PI) ;  
    ost << "Central value of Omega : " << omega_c * f_unit 
	<< " rad/s     f : " << freq * f_unit << " Hz" << endl ; 
    ost << "Central rotation period : " << 1000. / (freq * f_unit) << " ms"
	<< endl ;
  }

  // double compact = qpig/(4.*M_PI) * mass_g() / r_circ() ; 
  //ost << "Compactness G M_g /(c^2 R_circ) : " << compact << endl ;     
  
  double nphi_c = nphi.val_grid_point(0, 0, 0, 0) ;
  if ( (omega_c==0) && (nphi_c==0) ) {
    ost << "Central N^phi :               " << nphi_c << endl ;
  }
  else{
    ost << "Central N^phi/Omega :         " << nphi_c / omega_c << endl ;
  }
  
  ost << "Error on the virial identity GRV2 : " <<  grv2() << endl ; 
  double xgrv3 = grv3(&ost) ; 
  ost << "Error on the virial identity GRV3 : " << xgrv3 << endl ; 
  
  // cout << "ICI" << endl;
  double mom_quad_38si = mom_quad() * rho_unit * (pow(r_unit, double(5.)) 
  					  / double(1.e38) ) ;
//cout << "LA" << endl;
  ost << "Quadrupole moment Q : " << mom_quad_38si << " 10^38 kg m^2"
      << endl ; 

  /*ost << "Q / (M R_circ^2) :       " 
    << mom_quad() / ( mass_g() * pow( r_circ(), 2. ) ) << endl ; 
  ost << "c^4 Q / (G^2 M^3) :      " 
      << mom_quad() / ( pow(qpig/(4*M_PI), 2.) * pow(mass_g(), 3.) ) 
      << endl ; */
  
  ost << "Angular momentum J :      " 
      << angu_mom()/( qpig / (4* M_PI) * msol*msol) << " G M_sol^2 / c"
      << endl ; 
  /*  ost << "c J / (G M^2) :           " 
      << angu_mom()/( qpig / (4* M_PI) * pow(mass_g(), 2.) ) << endl ;    	 	*/
  
if (omega != __infinity) {
   double mom_iner = angu_mom() / omega ; 
   double mom_iner_38si = mom_iner * rho_unit * (pow(r_unit, double(5.)) 
						  / double(1.e38) ) ; 
   ost << "Moment of inertia:       " << mom_iner_38si << " 10^38 kg m^2"
	<< endl ; 
 }
  
  ost << "Ratio T/W :            " << tsw() << endl ; 
  ost << "Circumferential equatorial radius R_circ :     " 
      << r_circ()/km << " km" << endl ;  
  ost << "Coordinate equatorial radius r_eq : " << ray_eq()/km << " km" 
      << endl ;  
  ost << "Flattening r_pole/r_eq :        " << aplat() << endl ; 

  //cout << "ICI" << endl;
  
  int lsurf = nzet - 1; 
  int nt = mp.get_mg()->get_nt(lsurf) ; 
  int nr = mp.get_mg()->get_nr(lsurf) ; 
  ost << "Equatorial value of the velocity U:         " 
      << uuu.val_grid_point(nzet-1, 0, nt-1, nr-1) << " c" << endl ; 
  ost << "Redshift at the equator (forward) : " << z_eqf() << endl ; 
  ost << "Redshift at the equator (backward): " << z_eqb() << endl ; 
  ost << "Redshift at the pole              : " << z_pole() << endl ; 
  
  
  ost << "Central value of log(N)        : " 
      << logn.val_grid_point(0, 0, 0, 0) << endl ; 
  
  ost << "Central value of dzeta=log(AN) : " 
      << dzeta.val_grid_point(0, 0, 0, 0) << endl ; 
  
  if ( (omega_c==0) && (nphi_c==0) ) {
    ost << "Central N^phi :               " << nphi_c << endl ;
  }
  else{
    ost << "Central N^phi/Omega :         " << nphi_c / omega_c << endl ;
  }
  
  ost << "  ... w_shift (NB: components in the star Cartesian frame) [c] :  "
      << endl  
      << w_shift(1).val_grid_point(0, 0, 0, 0) << "  " 
      << w_shift(2).val_grid_point(0, 0, 0, 0) << "  " 
      << w_shift(3).val_grid_point(0, 0, 0, 0) << endl ; 
  
  ost << "Central value of khi_shift [km c] : " 
      << khi_shift.val_grid_point(0, 0, 0, 0) / km << endl ; 
  
  ost << "Central value of B^2 : " << b_car.val_grid_point(0,0,0,0) <<  endl ; 
  
  Tbl diff_a_b = diffrel( a_car, b_car ) ;
  ost << 
    "Relative discrepancy in each domain between the metric coef. A^2 and B^2 : "
      << endl ;
  for (int l=0; l<diff_a_b.get_taille(); l++) {
    ost << diff_a_b(l) << "  " ;
  }
  ost << endl ;
  
  // Approximate formula for R_isco = 6 R_g (1-(2/3)^1.5 j )
  // up to the first order in j
  
  /*double jdimless = angu_mom() / ( ggrav * pow(mass_g(), 2.) ) ;
  double r_grav = ggrav * mass_g() ;
  double r_isco_appr = 6. * r_grav * ( 1. - pow(2./3.,1.5) * jdimless ) ;*/
  
  // Approximate formula for the ISCO frequency
  //   freq_ms = 6^{-1.5}/2pi/R_g (1+11*6^(-1.5) j )
  // up to the first order in j

  /*  double f_isco_appr = ( 1. + 11. /6. /sqrt(6.) * jdimless ) / r_grav /
      (12. * M_PI ) / sqrt(6.) ;*/
  
  /*  ost << endl << "Innermost stable circular orbit (ISCO) : " << endl ;
  double xr_isco = r_isco(&ost) ;
  ost <<"    circumferential radius r_isco = "
      << xr_isco / km << " km" << endl ;
  ost <<"     (approx. 6M + 1st order in j : "
      << r_isco_appr / km << " km)" << endl ;
  ost <<"                      (approx. 6M : "
      << 6. * r_grav / km << " km)" << endl ;
  ost <<"    orbital frequency f_isco = "
      << f_isco() * f_unit << " Hz" << endl ;
  ost <<"     (approx. 1st order in j : "
  << f_isco_appr * f_unit << " Hz)" << endl ;*/
  
  
  return ost ;
  
}
}
