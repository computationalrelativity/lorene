/*
 *  Methods of class Star_rot
 *
 *    (see file star_rot.h for documentation).
 *
 */

/*
 *   Copyright (c) 2010 Eric Gourgoulhon
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

char star_rot_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.5  2010/02/08 10:56:30  j_novak
 * Added a few things missing for the reading from resulting file.
 *
 * Revision 1.4  2010/02/02 12:45:16  e_gourgoulhon
 * Improved the display (operator>>)
 *
 * Revision 1.3  2010/01/25 22:33:35  e_gourgoulhon
 * Debugging...
 *
 * Revision 1.2  2010/01/25 18:15:32  e_gourgoulhon
 * Added member unsurc2
 *
 * Revision 1.1  2010/01/24 16:09:39  e_gourgoulhon
 * New class Star_rot.
 *
 *
 * $Header$
 *
 */


// C headers
#include <math.h>
#include <assert.h>

// Lorene headers
#include "star_rot.h"
#include "eos.h"
#include "unites.h" 
#include "utilitaires.h"
#include "nbr_spx.h"



                   //--------------//
                   // Constructors //
                   //--------------//

// Standard constructor
// --------------------
Star_rot::Star_rot(Map& mpi, int nzet_i, bool relat, const Eos& eos_i)
		       : Star(mpi, nzet_i, eos_i),
			 relativistic(relat),
			 a_car(mpi), 
			 bbb(mpi), 
			 b_car(mpi), 
			 nphi(mpi), 
			 tnphi(mpi), 
			 uuu(mpi), 
			 nuf(mpi), 
			 nuq(mpi), 
			 dzeta(mpi), 
			 tggg(mpi), 
			 w_shift(mpi, CON, mp.get_bvect_cart()), 
			 khi_shift(mpi), 
			 tkij(mpi, COV, mp.get_bvect_cart()), 
			 ak_car(mpi), 
			 ssjm1_nuf(mpi), 
			 ssjm1_nuq(mpi), 
			 ssjm1_dzeta(mpi), 
			 ssjm1_tggg(mpi), 
			 ssjm1_khi(mpi), 
			 ssjm1_wshift(mpi, CON, mp.get_bvect_cart())			 
{

    // Parameter 1/c^2 is deduced from relativistic:
    unsurc2 = relativistic ? double(1) : double(0) ; 

    // Initialization to a static state : 
    omega = 0 ; 
    uuu = 0 ; 

    // Initialization to a flat metric : 
    a_car = 1 ;
    bbb = 1 ;
    bbb.std_spectral_base() ; 
    b_car = 1 ;
    nphi = 0 ;   
    tnphi = 0 ;   
    nuf = 0 ; 
    nuq = 0 ;
    dzeta = 0 ; 
    tggg = 0 ;   

    w_shift.set_etat_zero() ; 
    khi_shift =  0 ; 

    beta.set_etat_zero() ; 
    beta.set_triad( mp.get_bvect_cart() ) ;

    tkij.set_etat_zero() ; 

    ak_car = 0 ; 

    ssjm1_nuf = 0 ; 
    ssjm1_nuq = 0 ; 
    ssjm1_dzeta = 0 ; 
    ssjm1_tggg = 0 ; 
    ssjm1_khi = 0 ; 

    ssjm1_wshift.set_etat_zero() ; 
    
    // Pointers of derived quantities initialized to zero : 
    set_der_0x0() ;
    
}

// Copy constructor
// ----------------

Star_rot::Star_rot(const Star_rot& et)
		       : Star(et), 
			 relativistic(et.relativistic),
			 unsurc2(et.unsurc2),
			 omega(et.omega),
			 a_car(et.a_car), 
			 bbb(et.bbb), 
			 b_car(et.b_car), 
			 nphi(et.nphi), 
			 tnphi(et.tnphi), 
			 uuu(et.uuu), 
			 nuf(et.nuf), 
			 nuq(et.nuq),
			 dzeta(et.dzeta), 
			 tggg(et.tggg), 
			 w_shift(et.w_shift), 
			 khi_shift(et.khi_shift), 
			 tkij(et.tkij),
			 ak_car(et.ak_car), 
			 ssjm1_nuf(et.ssjm1_nuf), 
			 ssjm1_nuq(et.ssjm1_nuq), 
			 ssjm1_dzeta(et.ssjm1_dzeta), 
			 ssjm1_tggg(et.ssjm1_tggg), 
			 ssjm1_khi(et.ssjm1_khi), 
			 ssjm1_wshift(et.ssjm1_wshift)			 
{
    // Pointers of derived quantities initialized to zero : 
    set_der_0x0() ;
}    


// Constructor from a file
// -----------------------
Star_rot::Star_rot(Map& mpi, const Eos& eos_i, FILE* fich)
		       : Star(mpi, eos_i, fich), 
			 a_car(mpi), 
			 bbb(mpi), 
			 b_car(mpi), 
			 nphi(mpi), 
			 tnphi(mpi), 
			 uuu(mpi), 
			 nuf(mpi), 
			 nuq(mpi), 
			 dzeta(mpi), 
			 tggg(mpi), 
			 w_shift(mpi, CON, mp.get_bvect_cart()), 
			 khi_shift(mpi), 
			 tkij(mpi, COV, mp.get_bvect_cart()), 
			 ak_car(mpi), 
			 ssjm1_nuf(mpi), 
			 ssjm1_nuq(mpi), 
			 ssjm1_dzeta(mpi), 
			 ssjm1_tggg(mpi), 
			 ssjm1_khi(mpi), 
			 ssjm1_wshift(mpi, CON, mp.get_bvect_cart())			 
{

    // Star parameters
    // -----------------

    // relativistic is read in the file: 
    fread(&relativistic, sizeof(bool), 1, fich) ;	//## to be checked !	

    // Parameter 1/c^2 is deduced from relativistic:
    unsurc2 = relativistic ? double(1) : double(0) ; 

    // omega is read in the file:     
    fread_be(&omega, sizeof(double), 1, fich) ;		
    	  
   
    // Read of the saved fields:
    // ------------------------

    Scalar nuf_file(mp, *(mp.get_mg()), fich) ; 
    nuf = nuf_file ; 
        
    Scalar nuq_file(mp, *(mp.get_mg()), fich) ; 
    nuq = nuq_file ; 
        
    Scalar dzeta_file(mp, *(mp.get_mg()), fich) ; 
    dzeta = dzeta_file ; 
        
    Scalar tggg_file(mp, *(mp.get_mg()), fich) ; 
    tggg = tggg_file ; 
        
    Vector w_shift_file(mp, mp.get_bvect_cart(), fich) ; 
    w_shift = w_shift_file ;
    
    Scalar khi_shift_file(mp, *(mp.get_mg()), fich) ; 
    khi_shift = khi_shift_file ;
    
    fait_shift() ;	    // constructs shift from w_shift and khi_shift
    fait_nphi() ;       // constructs N^phi from (N^x,N^y,N^z)
    
    Scalar ssjm1_nuf_file(mp, *(mp.get_mg()), fich) ; 
    ssjm1_nuf = ssjm1_nuf_file ; 

    Scalar ssjm1_nuq_file(mp, *(mp.get_mg()), fich) ; 
    ssjm1_nuq = ssjm1_nuq_file ; 

    Scalar ssjm1_dzeta_file(mp, *(mp.get_mg()), fich) ; 
    ssjm1_dzeta = ssjm1_dzeta_file ; 

    Scalar ssjm1_tggg_file(mp, *(mp.get_mg()), fich) ; 
    ssjm1_tggg = ssjm1_tggg_file ; 

    Scalar ssjm1_khi_file(mp, *(mp.get_mg()), fich) ; 
    ssjm1_khi = ssjm1_khi_file ; 

    Vector ssjm1_wshift_file(mp, mp.get_bvect_cart(), fich) ; 
    ssjm1_wshift = ssjm1_wshift_file ; 

    // All other fields are initialized to zero : 
    // ----------------------------------------
    a_car = 0 ; 
    bbb = 0 ; 
    b_car = 0 ; 
    uuu = 0 ; 
    tkij.set_etat_zero() ; 
    ak_car = 0 ;  

    // Pointers of derived quantities initialized to zero 
    // --------------------------------------------------
    set_der_0x0() ;
    
}

			    //------------//
			    // Destructor //
			    //------------//

Star_rot::~Star_rot(){

    del_deriv() ; 

}

		//----------------------------------//
		// Management of derived quantities //
		//----------------------------------//

void Star_rot::del_deriv() const {

    Star::del_deriv() ; 

    if (p_angu_mom != 0x0) delete p_angu_mom ; 
    if (p_tsw != 0x0) delete p_tsw ; 
    if (p_grv2 != 0x0) delete p_grv2 ; 
    if (p_grv3 != 0x0) delete p_grv3 ; 
    if (p_r_circ != 0x0) delete p_r_circ ; 
    if (p_aplat != 0x0) delete p_aplat ; 
    if (p_z_eqf != 0x0) delete p_z_eqf ; 
    if (p_z_eqb != 0x0) delete p_z_eqb ; 
    if (p_z_pole != 0x0) delete p_z_pole ; 
    if (p_mom_quad != 0x0) delete p_mom_quad ;
    if (p_r_isco != 0x0) delete p_r_isco ;
    if (p_f_isco != 0x0) delete p_f_isco ;
    if (p_lspec_isco != 0x0) delete p_lspec_isco ;
    if (p_espec_isco != 0x0) delete p_espec_isco ;
    if (p_f_eq != 0x0) delete p_f_eq ;
    
    Star_rot::set_der_0x0() ; 
}			    


void Star_rot::set_der_0x0() const {

    Star::set_der_0x0() ;

    p_angu_mom = 0x0 ; 
    p_tsw = 0x0 ;
    p_grv2 = 0x0 ;
    p_grv3 = 0x0 ;
    p_r_circ = 0x0 ;
    p_aplat = 0x0 ;
    p_z_eqf = 0x0 ;
    p_z_eqb = 0x0 ;
    p_z_pole = 0x0 ;
    p_mom_quad = 0x0 ;
    p_r_isco = 0x0 ;
    p_f_isco = 0x0 ;
    p_lspec_isco = 0x0 ;
    p_espec_isco = 0x0 ;
    p_f_eq = 0x0 ;
    
}			    

void Star_rot::del_hydro_euler() {

    Star::del_hydro_euler() ; 

    del_deriv() ; 

}			    

			    //--------------//
			    //  Assignment  //
			    //--------------//

// Assignment to another Star_rot
// --------------------------------
void Star_rot::operator=(const Star_rot& et) {

    // Assignment of quantities common to all the derived classes of Star
    Star::operator=(et) ;	    

    // Assignement of proper quantities of class Star_rot
    relativistic = et.relativistic ; 
    unsurc2 = et.unsurc2 ; 
    omega = et.omega ; 

    a_car = et.a_car ; 
    bbb = et.bbb ; 
    b_car = et.b_car ; 
    nphi = et.nphi ; 
    tnphi = et.tnphi ; 
    uuu = et.uuu ; 
    nuf = et.nuf ; 
    nuq = et.nuq ; 
    dzeta = et.dzeta ; 
    tggg = et.tggg ; 
    w_shift = et.w_shift ;
    khi_shift = et.khi_shift ;
    tkij = et.tkij ; 
    ak_car = et.ak_car ;
    ssjm1_nuf = et.ssjm1_nuf ;
    ssjm1_nuq = et.ssjm1_nuq ;
    ssjm1_dzeta = et.ssjm1_dzeta ; 
    ssjm1_tggg = et.ssjm1_tggg ;
    ssjm1_khi = et.ssjm1_khi ;
    ssjm1_wshift = et.ssjm1_wshift ; 
    
    del_deriv() ;  // Deletes all derived quantities

}	


			    //--------------//
			    //	  Outputs   //
			    //--------------//

// Save in a file
// --------------
void Star_rot::sauve(FILE* fich) const {
    
    Star::sauve(fich) ; 
    
    fwrite(&relativistic, sizeof(bool), 1, fich) ;		
    fwrite_be(&omega, sizeof(double), 1, fich) ;		
    
    nuf.sauve(fich) ; 
    nuq.sauve(fich) ; 
    dzeta.sauve(fich) ; 
    tggg.sauve(fich) ; 
    w_shift.sauve(fich) ; 
    khi_shift.sauve(fich) ; 
    
    ssjm1_nuf.sauve(fich) ; 
    ssjm1_nuq.sauve(fich) ; 
    ssjm1_dzeta.sauve(fich) ; 
    ssjm1_tggg.sauve(fich) ; 
    ssjm1_khi.sauve(fich) ; 
    ssjm1_wshift.sauve(fich) ; 
    
    
}

// Printing
// --------

ostream& Star_rot::operator>>(ostream& ost) const {
    
  using namespace Unites ;

    Star::operator>>(ost) ; 

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
    if (relativistic) {
	ost << "Relativistic star" << endl ; 
    }
    else {
	ost << "Newtonian star" << endl ; 
    }
    double compact = qpig/(4.*M_PI) * mass_g() / r_circ() ; 
    ost << "Compactness G M_g /(c^2 R_circ) : " << compact << endl ;     
       
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

    double mom_quad_38si = mom_quad() * rho_unit * (pow(r_unit, double(5.)) 
 							/ double(1.e38) ) ;
    ost << "Quadrupole moment Q : " << mom_quad_38si << " 10^38 kg m^2"
         << endl ; 
    ost << "Q / (M R_circ^2) :       " 
	 << mom_quad() / ( mass_g() * pow( r_circ(), 2. ) ) << endl ; 
    ost << "c^4 Q / (G^2 M^3) :      " 
	 << mom_quad() / ( pow(qpig/(4*M_PI), 2.) * pow(mass_g(), 3.) ) 
	 << endl ; 

    ost << "Angular momentum J :      " 
	 << angu_mom()/( qpig / (4* M_PI) * msol*msol) << " G M_sol^2 / c"
	 << endl ; 
    ost << "c J / (G M^2) :           " 
	 << angu_mom()/( qpig / (4* M_PI) * pow(mass_g(), 2.) ) << endl ;    	 	

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
    double jdimless = angu_mom() / ( ggrav * pow(mass_g(), 2.) ) ;
    double r_grav = ggrav * mass_g() ;
    double r_isco_appr = 6. * r_grav * ( 1. - pow(2./3.,1.5) * jdimless ) ;

	// Approximate formula for the ISCO frequency
	//   freq_ms = 6^{-1.5}/2pi/R_g (1+11*6^(-1.5) j )
	// up to the first order in j
	double f_isco_appr = ( 1. + 11. /6. /sqrt(6.) * jdimless ) / r_grav /
      						(12. * M_PI ) / sqrt(6.) ;

    ost << endl << "Innermost stable circular orbit (ISCO) : " << endl ;
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
	   	<< f_isco_appr * f_unit << " Hz)" << endl ;
    	

    return ost ;
    
}


void Star_rot::partial_display(ostream& ost) const {
    
  using namespace Unites ;

    double omega_c = get_omega_c() ; 
    double freq = omega_c / (2.*M_PI) ;  
    ost << "Central Omega : " << omega_c * f_unit 
        << " rad/s     f : " << freq * f_unit << " Hz" << endl ; 
    ost << "Rotation period : " << 1000. / (freq * f_unit) << " ms"
	 << endl ;
    ost << endl << "Central enthalpy : " << ent.val_grid_point(0,0,0,0) << " c^2" << endl ; 
    ost << "Central proper baryon density : " << nbar.val_grid_point(0,0,0,0) 
	<< " x 0.1 fm^-3" << endl ; 
    ost << "Central proper energy density : " << ener.val_grid_point(0,0,0,0) 
	<< " rho_nuc c^2" << endl ; 
    ost << "Central pressure : " << press.val_grid_point(0,0,0,0) 
	<< " rho_nuc c^2" << endl ; 

    ost << "Central value of log(N)        : " 
	<< logn.val_grid_point(0, 0, 0, 0) << endl ; 
    ost << "Central lapse N :      " << nn.val_grid_point(0,0,0,0) <<  endl ; 
    ost << "Central value of dzeta=log(AN) : " 
	<< dzeta.val_grid_point(0, 0, 0, 0) << endl ; 
    ost << "Central value of A^2 : " << a_car.val_grid_point(0,0,0,0) <<  endl ; 
    ost << "Central value of B^2 : " << b_car.val_grid_point(0,0,0,0) <<  endl ; 

    double nphi_c = nphi.val_grid_point(0, 0, 0, 0) ;
    if ( (omega_c==0) && (nphi_c==0) ) {
		ost << "Central N^phi :               " << nphi_c << endl ;
    }
    else{
		ost << "Central N^phi/Omega :         " << nphi_c / omega_c 
							<< endl ;
    }
	    

    int lsurf = nzet - 1; 
    int nt = mp.get_mg()->get_nt(lsurf) ; 
    int nr = mp.get_mg()->get_nr(lsurf) ; 
    ost << "Equatorial value of the velocity U:         " 
	 << uuu.val_grid_point(nzet-1, 0, nt-1, nr-1) << " c" << endl ; 

    ost << endl 
	<< "Coordinate equatorial radius r_eq =    " 
	<< ray_eq()/km << " km" << endl ;  
    ost << "Flattening r_pole/r_eq :        " << aplat() << endl ; 

}


double Star_rot::get_omega_c() const {
    
    return omega ; 
    
}

// display_poly
// ------------

void Star_rot::display_poly(ostream& ost) const {

  using namespace Unites ;

  const Eos_poly* p_eos_poly = dynamic_cast<const Eos_poly*>( &eos ) ; 	  

    if (p_eos_poly != 0x0) {

	double kappa = p_eos_poly->get_kap() ; 

	// kappa^{n/2}
	double kap_ns2 = pow( kappa,  0.5 /(p_eos_poly->get_gam()-1) ) ; 
    
	// Polytropic unit of length in terms of r_unit : 
	double r_poly = kap_ns2 / sqrt(ggrav) ; 
    
	// Polytropic unit of time in terms of t_unit :
	double t_poly = r_poly ; 

	// Polytropic unit of mass in terms of m_unit :
	double m_poly = r_poly / ggrav ; 
    
	// Polytropic unit of angular momentum in terms of j_unit :
	double j_poly = r_poly * r_poly / ggrav ; 
    
	// Polytropic unit of density in terms of rho_unit :
	double rho_poly = 1. / (ggrav * r_poly * r_poly) ;  

	ost.precision(10) ; 
	ost << endl << "Quantities in polytropic units : " << endl ; 
	ost	 << "==============================" << endl ; 
	ost << " ( r_poly = " << r_poly / km << " km )" << endl ; 
	ost << "  n_c     : " << nbar.val_grid_point(0, 0, 0, 0) / rho_poly << endl ; 
	ost << "  e_c     : " << ener.val_grid_point(0, 0, 0, 0) / rho_poly << endl ; 
	ost << "  Omega_c : " << get_omega_c() * t_poly << endl ; 
	ost << "  P_c     : " << 2.*M_PI / get_omega_c() / t_poly << endl ; 
	ost << "  M_bar   : " << mass_b() / m_poly << endl ; 
	ost << "  M       : " << mass_g() / m_poly << endl ; 
	ost << "  J	  : " << angu_mom() / j_poly << endl ; 
	ost << "  r_eq	  : " << ray_eq() / r_poly << endl ; 
	ost << "  R_circ  : " << r_circ() / r_poly << endl ; 
	
    
    }
    

} 



			    //-------------------------//
			    //	Computational methods  //
			    //-------------------------//
			    
void Star_rot::fait_shift() {

    Vector d_khi = khi_shift.derive_con( mp.flat_met_cart() ) ; 
    
    d_khi.dec_dzpuis(2) ;   // divide by r^2 in the external compactified domain  
 
    // x_k dW^k/dx_i
      Scalar xx(mp) ; 
      Scalar yy(mp) ; 
      Scalar zz(mp) ; 
      xx = mp.x ; 
      yy = mp.y ; 
      zz = mp.z ; 

      Tensor d_w = w_shift.derive_con( mp.flat_met_cart() ) ;
      
      Vector x_d_w(mp, CON, mp.get_bvect_cart()) ;
      x_d_w.set(1) = xx*d_w(1,1) + yy*d_w(2,1) + zz*d_w(3,1) ; 
      x_d_w.set(2) = xx*d_w(1,2) + yy*d_w(2,2) + zz*d_w(3,2) ; 
      x_d_w.set(3) = xx*d_w(1,3) + yy*d_w(2,3) + zz*d_w(3,3) ; 
      x_d_w.std_spectral_base() ;
      x_d_w.dec_dzpuis(2) ;

      double lambda = double(1) / double(3) ; 

      beta = - (lambda+2)/2./(lambda+1) * w_shift
		+ (lambda/2./(lambda+1))  * (d_khi + x_d_w) ;      
    
} 

void Star_rot::fait_nphi() {

    if ( (beta(1).get_etat() == ETATZERO) && (beta(2).get_etat() == ETATZERO) ) {
	tnphi = 0 ; 
	nphi = 0 ;
	return ;  
    }

    // Computation of tnphi
    // --------------------
    
    mp.comp_p_from_cartesian( -beta(1), -beta(2), tnphi ) ; 

    // Computation of nphi
    // -------------------
    
    nphi = tnphi ; 
    nphi.div_rsint() ; 
    
}

