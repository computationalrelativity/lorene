/*
 * Methods for the class Etoile_rot
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
 * Revision 1.10  2023/07/04 08:59:53  j_novak
 * Added method "r_circ_merid()" to compute the circumferential meridional radius.
 *
 * Revision 1.9  2021/04/13 11:24:53  j_novak
 * Corrected a bug in Etoile_rot::fait_shift() which was missing the np=1 case.
 *
 * Revision 1.8  2016/12/05 16:17:55  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.7  2015/12/03 14:17:24  j_novak
 * Check added for the computation of area (thanks S. Koeppel).
 *
 * Revision 1.6  2015/06/10 14:37:44  a_sourie
 * Corrected the formula for the quadrupole.
 *
 * Revision 1.5  2014/10/13 08:52:59  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.4  2004/03/25 10:29:07  j_novak
 * All LORENE's units are now defined in the namespace Unites (in file unites.h).
 *
 * Revision 1.3  2001/12/06 15:11:43  jl_zdunik
 * Introduction of the new function f_eq() in the class Etoile_rot
 *
 * Revision 1.2  2001/12/04 21:27:53  e_gourgoulhon
 *
 * All writing/reading to a binary file are now performed according to
 * the big endian convention, whatever the system is big endian or
 * small endian, thanks to the functions fwrite_be and fread_be
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 2.17  2001/10/24  15:36:20  eric
 * Ajout de la fonction display_poly.
 *
 * Revision 2.16  2001/10/16  14:49:02  eric
 * Appel de get_omega_c() pour avoir la valeur centrale de Omega.
 * Affichage different si rotation differentielle.
 *
 * Revision 2.15  2001/09/13  08:32:01  eric
 * Ajout du facteur de compacite M/R dans l'affichage.
 *
 * Revision 2.14  2001/06/20  14:20:56  novak
 * Appel a Etoile_rot::set_der0x0 dans del_deriv (au lieu de set_der0x0
 * tout court).
 *
 * Revision 2.13  2001/03/26 09:30:58  jlz
 * New members p_espec_isco and p_lspec_isco.
 *
 * Revision 2.12  2000/11/20  21:42:02  eric
 * Appel de fait_nphi() dans le constructeur par lecture de fichier.
 *
 * Revision 2.11  2000/11/18  23:18:30  eric
 * Modifs affichage.
 *
 * Revision 2.10  2000/11/18  21:09:57  eric
 * Ajout des membres  p_r_isco et p_f_isco.
 *
 * Revision 2.9  2000/11/07  16:33:08  eric
 * Modif affichage.
 *
 * Revision 2.8  2000/10/12  15:37:01  eric
 * Ajout de la fonction fait_nphi().
 *
 * Revision 2.7  2000/09/18  16:15:12  eric
 * Ajout du membre tkij.
 *
 * Revision 2.6  2000/08/31  15:38:00  eric
 * Bases spectrales standards pour bbb et b_car dans le constructeur
 * standard (initialisation a la metrique plate).
 *
 * Revision 2.5  2000/08/31  11:25:45  eric
 * Ajout des membres tnphi et ak_car.
 *
 * Revision 2.4  2000/08/25  12:28:29  eric
 * Modif affichage.
 *
 * Revision 2.3  2000/08/18  14:01:59  eric
 * Ajout de partial_display
 *
 * Revision 2.2  2000/08/17  12:40:04  eric
 * *** empty log message ***
 *
 * Revision 2.1  2000/07/21  16:31:26  eric
 * *** empty log message ***
 *
 * Revision 1.1  2000/07/20  15:32:37  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */

// Headers C
#include "math.h"

// Headers Lorene
#include "etoile.h"
#include "eos.h"
#include "nbr_spx.h"
#include "utilitaires.h"
#include "unites.h"	    

			    //--------------//
			    // Constructors //
			    //--------------//

// Standard constructor
// --------------------
namespace Lorene {
Etoile_rot::Etoile_rot(Map& mpi, int nzet_i, bool relat, const Eos& eos_i)
		       : Etoile(mpi, nzet_i, relat, eos_i), 
			 bbb(mpi), 
			 b_car(mpi), 
			 nphi(mpi), 
			 tnphi(mpi), 
			 uuu(mpi), 
			 logn(logn_auto), 
			 nuf(mpi), 
			 nuq(mpi), 
			 dzeta(beta_auto), 
			 tggg(mpi), 
			 w_shift(mpi, 1, CON, mp.get_bvect_cart()), 
			 khi_shift(mpi), 
			 tkij(mpi, 2, COV, mp.get_bvect_cart()), 
			 ak_car(mpi), 
			 ssjm1_nuf(mpi), 
			 ssjm1_nuq(mpi), 
			 ssjm1_dzeta(mpi), 
			 ssjm1_tggg(mpi), 
			 ssjm1_khi(mpi), 
			 ssjm1_wshift(mpi, 1, CON, mp.get_bvect_cart())			 
{

    // Initialization to a static state : 
    omega = 0 ; 
    uuu = 0 ; 

    // Initialization to a flat metric : 
    bbb = 1 ;
    bbb.set_std_base() ; 
    b_car = 1 ;
    b_car.set_std_base() ;
    nphi = 0 ;   
    tnphi = 0 ;   
    nuf = 0 ; 
    nuq = 0 ; 
    tggg = 0 ;   

    w_shift.set_etat_qcq() ; 
    for (int i=0; i<3; i++) {
	w_shift.set(i) = 0 ; 
    }

    khi_shift.set_etat_qcq() ; 
    khi_shift.set() = 0 ; 

    tkij.set_etat_zero() ; 

    ak_car = 0 ; 

    ssjm1_nuf = 0 ; 
    ssjm1_nuq = 0 ; 
    ssjm1_dzeta = 0 ; 
    ssjm1_tggg = 0 ; 
    ssjm1_khi = 0 ; 

    ssjm1_wshift.set_etat_qcq() ; 
    for (int i=0; i<3; i++) {
	ssjm1_wshift.set(i) = 0 ; 
    }
    
    // Pointers of derived quantities initialized to zero : 
    set_der_0x0() ;
    
}

// Copy constructor
// ----------------

Etoile_rot::Etoile_rot(const Etoile_rot& et)
		       : Etoile(et), 
			 bbb(et.bbb), 
			 b_car(et.b_car), 
			 nphi(et.nphi), 
			 tnphi(et.tnphi), 
			 uuu(et.uuu), 
			 logn(logn_auto), 
			 nuf(et.nuf), 
			 nuq(et.nuq),
			 dzeta(beta_auto), 
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
    omega = et.omega ; 

    // Pointers of derived quantities initialized to zero : 
    set_der_0x0() ;
}    


// Constructor from a file
// -----------------------
Etoile_rot::Etoile_rot(Map& mpi, const Eos& eos_i, FILE* fich)
		       : Etoile(mpi, eos_i, fich), 
			 bbb(mpi), 
			 b_car(mpi), 
			 nphi(mpi), 
			 tnphi(mpi), 
			 uuu(mpi), 
			 logn(logn_auto), 
			 nuf(mpi), 
			 nuq(mpi), 
			 dzeta(beta_auto), 
			 tggg(mpi), 
			 w_shift(mpi, 1, CON, mp.get_bvect_cart()), 
			 khi_shift(mpi), 
			 tkij(mpi, 2, COV, mp.get_bvect_cart()), 
			 ak_car(mpi), 
			 ssjm1_nuf(mpi), 
			 ssjm1_nuq(mpi), 
			 ssjm1_dzeta(mpi), 
			 ssjm1_tggg(mpi), 
			 ssjm1_khi(mpi), 
			 ssjm1_wshift(mpi, 1, CON, mp.get_bvect_cart())			 
{

    // Etoile parameters
    // -----------------

    // omega is read in the file:     
    fread_be(&omega, sizeof(double), 1, fich) ;		
    	  
   
    // Read of the saved fields:
    // ------------------------

    Tenseur nuf_file(mp, fich) ; 
    nuf = nuf_file ; 
        
    Tenseur nuq_file(mp, fich) ; 
    nuq = nuq_file ; 
        
    Tenseur tggg_file(mp, fich) ; 
    tggg = tggg_file ; 
        
    Tenseur w_shift_file(mp, mp.get_bvect_cart(), fich) ; 
    w_shift = w_shift_file ;
    
    Tenseur khi_shift_file(mp, fich) ; 
    khi_shift = khi_shift_file ;
    
    fait_shift() ;	    // constructs shift from w_shift and khi_shift
    fait_nphi() ;       // constructs N^phi from (N^x,N^y,N^z)
    
    Cmp ssjm1_nuf_file(mp, *(mp.get_mg()), fich) ; 
    ssjm1_nuf = ssjm1_nuf_file ; 

    Cmp ssjm1_nuq_file(mp, *(mp.get_mg()), fich) ; 
    ssjm1_nuq = ssjm1_nuq_file ; 

    Cmp ssjm1_dzeta_file(mp, *(mp.get_mg()), fich) ; 
    ssjm1_dzeta = ssjm1_dzeta_file ; 

    Cmp ssjm1_tggg_file(mp, *(mp.get_mg()), fich) ; 
    ssjm1_tggg = ssjm1_tggg_file ; 

    Cmp ssjm1_khi_file(mp, *(mp.get_mg()), fich) ; 
    ssjm1_khi = ssjm1_khi_file ; 

    Tenseur ssjm1_wshift_file(mp, mp.get_bvect_cart(), fich) ; 
    ssjm1_wshift = ssjm1_wshift_file ; 

    // All other fields are initialized to zero : 
    // ----------------------------------------
    bbb = 0 ; 
    b_car = 0 ; 
    uuu = 0 ; 

    // Pointers of derived quantities initialized to zero 
    // --------------------------------------------------
    set_der_0x0() ;
    
}

			    //------------//
			    // Destructor //
			    //------------//

Etoile_rot::~Etoile_rot(){

    del_deriv() ; 

}

		//----------------------------------//
		// Management of derived quantities //
		//----------------------------------//

void Etoile_rot::del_deriv() const {

    Etoile::del_deriv() ; 

    if (p_angu_mom != 0x0) delete p_angu_mom ; 
    if (p_tsw != 0x0) delete p_tsw ; 
    if (p_grv2 != 0x0) delete p_grv2 ; 
    if (p_grv3 != 0x0) delete p_grv3 ; 
    if (p_r_circ != 0x0) delete p_r_circ ; 
    if (p_r_circ_merid != 0x0) delete p_r_circ_merid ; 
    if (p_area != 0x0) delete p_area ;
    if (p_aplat != 0x0) delete p_aplat ; 
    if (p_z_eqf != 0x0) delete p_z_eqf ; 
    if (p_z_eqb != 0x0) delete p_z_eqb ; 
    if (p_z_pole != 0x0) delete p_z_pole ; 
    if (p_mom_quad != 0x0) delete p_mom_quad ;
    if (p_mom_quad_old != 0x0) delete p_mom_quad_old ;  
    if (p_mom_quad_Bo != 0x0) delete p_mom_quad_Bo ;
    if (p_r_isco != 0x0) delete p_r_isco ;
    if (p_f_isco != 0x0) delete p_f_isco ;
    if (p_lspec_isco != 0x0) delete p_lspec_isco ;
    if (p_espec_isco != 0x0) delete p_espec_isco ;
    if (p_f_eq != 0x0) delete p_f_eq ;
    
    Etoile_rot::set_der_0x0() ; 
}			    




void Etoile_rot::set_der_0x0() const {

    Etoile::set_der_0x0() ;

    p_angu_mom = 0x0 ; 
    p_tsw = 0x0 ;
    p_grv2 = 0x0 ;
    p_grv3 = 0x0 ;
    p_r_circ = 0x0 ;
    p_r_circ_merid = 0x0 ;
    p_area = 0x0 ;
    p_aplat = 0x0 ;
    p_z_eqf = 0x0 ;
    p_z_eqb = 0x0 ;
    p_z_pole = 0x0 ;
    p_mom_quad = 0x0 ;
    p_mom_quad_old = 0x0 ;
    p_mom_quad_Bo = 0x0 ;
    p_r_isco = 0x0 ;
    p_f_isco = 0x0 ;
    p_lspec_isco = 0x0 ;
    p_espec_isco = 0x0 ;
    p_f_eq = 0x0 ;
    
}			    

void Etoile_rot::del_hydro_euler() {

    Etoile::del_hydro_euler() ; 

    del_deriv() ; 

}			    


			    //--------------//
			    //  Assignment  //
			    //--------------//

// Assignment to another Etoile_rot
// --------------------------------
void Etoile_rot::operator=(const Etoile_rot& et) {

    // Assignment of quantities common to all the derived classes of Etoile
    Etoile::operator=(et) ;	    

    // Assignement of proper quantities of class Etoile_rot
    omega = et.omega ; 

    bbb = et.bbb ; 
    b_car = et.b_car ; 
    nphi = et.nphi ; 
    tnphi = et.tnphi ; 
    uuu = et.uuu ; 
    nuf = et.nuf ; 
    nuq = et.nuq ; 
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
void Etoile_rot::sauve(FILE* fich) const {
    
    Etoile::sauve(fich) ; 
    
    fwrite_be(&omega, sizeof(double), 1, fich) ;		
    
    nuf.sauve(fich) ; 
    nuq.sauve(fich) ; 
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

ostream& Etoile_rot::operator>>(ostream& ost) const {
    
  using namespace Unites ;

    Etoile::operator>>(ost) ; 

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
    
       
    double nphi_c = nphi()(0, 0, 0, 0) ;
    if ( (omega_c==0) && (nphi_c==0) ) {
	 	ost << "Central N^phi :               " << nphi_c << endl ;
    }
    else{
		ost << "Central N^phi/Omega :         " << nphi_c / omega_c << endl ;
    }
	    
    ost << "Error on the virial identity GRV2 : " << endl ; 
    ost << "GRV2 = " << grv2() << endl ; 
    ost << "Error on the virial identity GRV3 : " << endl ; 
    double xgrv3 = grv3(&ost) ; 
    ost << "GRV3 = " << xgrv3 << endl ; 

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
    ost << "Circumferential equatorial radius R_circ_eq :     " 
	<< r_circ()/km << " km" << endl ; 
    if (mp.get_mg()->get_np(0) == 1) {
      ost << "Circumferential meridional radius R_circ_merid :     " 
	  << r_circ_merid()/km << " km" << endl ;
      ost << "Flattening r_circ_merid/r_circ_eq :        "
	  << r_circ_merid() / r_circ() << endl ; 
      ost << "Surface area :   " << area()/(km*km) << " km^2" << endl ;
      ost << "Mean radius R_mean :    " 
	  << mean_radius()/km << " km" << endl ;
    } else {
      ost << 
    "Skipping polar radius / surface statements due to number of points in phi direction np > 1" 
	  << endl;
    }
    ost << "Coordinate equatorial radius r_eq : " << ray_eq()/km << " km" 
	 << endl ;  
    ost << "Flattening r_pole/r_eq :        " << aplat() << endl ; 

    double compact = qpig/(4.*M_PI) * mass_g() / r_circ() ; 
    ost << "Compaction parameter M_g / R_circ : " << compact << endl ; 

    int lsurf = nzet - 1; 
    int nt = mp.get_mg()->get_nt(lsurf) ; 
    int nr = mp.get_mg()->get_nr(lsurf) ; 
    ost << "Equatorial value of the velocity U:         " 
	 << uuu()(nzet-1, 0, nt-1, nr-1) << " c" << endl ; 
    ost << "Redshift at the equator (forward) : " << z_eqf() << endl ; 
    ost << "Redshift at the equator (backward): " << z_eqb() << endl ; 
    ost << "Redshift at the pole              : " << z_pole() << endl ; 


    ost << "Central value of log(N)        : " 
	<< logn()(0, 0, 0, 0) << endl ; 

    ost << "Central value of dzeta=log(AN) : " 
	<< dzeta()(0, 0, 0, 0) << endl ; 

    if ( (omega_c==0) && (nphi_c==0) ) {
		ost << "Central N^phi :               " << nphi_c << endl ;
    }
    else{
		ost << "Central N^phi/Omega :         " << nphi_c / omega_c << endl ;
    }

    ost << "  ... w_shift (NB: components in the star Cartesian frame) [c] :  "
	<< endl  
	<< w_shift(0)(0, 0, 0, 0) << "  " 
	<< w_shift(1)(0, 0, 0, 0) << "  " 
	<< w_shift(2)(0, 0, 0, 0) << endl ; 

    ost << "Central value of khi_shift [km c] : " 
        << khi_shift()(0, 0, 0, 0) / km << endl ; 

    ost << "Central value of B^2 : " << b_car()(0,0,0,0) <<  endl ; 

    Tbl diff_a_b = diffrel( a_car(), b_car() ) ;
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


void Etoile_rot::partial_display(ostream& ost) const {
    
  using namespace Unites ;

    double omega_c = get_omega_c() ; 
    double freq = omega_c / (2.*M_PI) ;  
    ost << "Central Omega : " << omega_c * f_unit 
        << " rad/s     f : " << freq * f_unit << " Hz" << endl ; 
    ost << "Rotation period : " << 1000. / (freq * f_unit) << " ms"
	 << endl ;
    ost << endl << "Central enthalpy : " << ent()(0,0,0,0) << " c^2" << endl ; 
    ost << "Central proper baryon density : " << nbar()(0,0,0,0) 
	<< " x 0.1 fm^-3" << endl ; 
    ost << "Central proper energy density : " << ener()(0,0,0,0) 
	<< " rho_nuc c^2" << endl ; 
    ost << "Central pressure : " << press()(0,0,0,0) 
	<< " rho_nuc c^2" << endl ; 

    ost << "Central value of log(N)        : " 
	<< logn()(0, 0, 0, 0) << endl ; 
    ost << "Central lapse N :      " << nnn()(0,0,0,0) <<  endl ; 
    ost << "Central value of dzeta=log(AN) : " 
	<< dzeta()(0, 0, 0, 0) << endl ; 
    ost << "Central value of A^2 : " << a_car()(0,0,0,0) <<  endl ; 
    ost << "Central value of B^2 : " << b_car()(0,0,0,0) <<  endl ; 

    double nphi_c = nphi()(0, 0, 0, 0) ;
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
	 << uuu()(nzet-1, 0, nt-1, nr-1) << " c" << endl ; 

    ost << endl 
	<< "Coordinate equatorial radius r_eq =    " 
	<< ray_eq()/km << " km" << endl ;  
    ost << "Flattening r_pole/r_eq :        " << aplat() << endl ; 
    if (mp.get_mg()->get_np(0) == 1)
      ost << "Flattening r_circ_pole/r_circ_eq :        "
	  << r_circ_merid() / r_circ() << endl ; 

}


double Etoile_rot::get_omega_c() const {
    
    return omega ; 
    
}


// display_poly
// ------------

void Etoile_rot::display_poly(ostream& ost) const {

  using namespace Unites ;

  const Eos_poly* p_eos_poly = dynamic_cast<const Eos_poly*>( &eos ) ; 	  

    if (p_eos_poly != 0x0) {

	double kappa = p_eos_poly->get_kap() ; 
	double gamma = p_eos_poly->get_gam() ;  ; 

	// kappa^{n/2}
	double kap_ns2 = pow( kappa,  0.5 /(gamma-1) ) ; 
    
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
	ost << "  n_c     : " << nbar()(0, 0, 0, 0) / rho_poly << endl ; 
	ost << "  e_c     : " << ener()(0, 0, 0, 0) / rho_poly << endl ; 
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
			    //	Computational routines //
			    //-------------------------//
			    
void Etoile_rot::fait_shift() {

    Tenseur d_khi = khi_shift.gradient() ; 
    
    if (d_khi.get_etat() == ETATQCQ) { 
	d_khi.dec2_dzpuis() ;   // divide by r^2 in the external compactified
				// domain
    }
 
    // x_k dW^k/dx_i
    
    Tenseur x_d_w = skxk( w_shift.gradient() ) ;
    x_d_w.dec_dzpuis() ;
    
    double lambda = double(1) / double(3) ;

    if ( mp.get_mg()->get_np(0) == 1 ) {
      lambda = 0 ; 
    }

    // The final computation is done component by component because
    // d_khi and x_d_w are covariant comp. whereas w_shift is
    // contravariant
    
    shift.set_etat_qcq() ; 
    
    for (int i=0; i<3; i++) {
	shift.set(i) = (lambda+2)/2./(lambda+1) * w_shift(i)
		- (lambda/2./(lambda+1))  * (d_khi(i) + x_d_w(i)) ;      
    }
    
    shift.set_triad( *(w_shift.get_triad()) ) ; 
    
} 



void Etoile_rot::fait_nphi() {

    if ( shift.get_etat() == ETATZERO ) {
	tnphi = 0 ; 
	nphi = 0 ;
	return ;  
    }

    assert( shift.get_etat() == ETATQCQ ) ; 
    
    // Computation of tnphi
    // --------------------
    tnphi.set_etat_qcq() ; 
    
    mp.comp_p_from_cartesian( shift(0), shift(1), tnphi.set() ) ; 

    // Computation of nphi
    // -------------------
    
    nphi = tnphi ; 
    (nphi.set()).div_rsint() ; 
    
}
}
