/*
 * Methods for two fluids rotating relativistic stars.
 *
 * See the file et_rot_bifluid.h for documentation
 *
 */

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


 

/*
 * $Id$
 * $Log$
 * Revision 1.17  2017/10/06 12:36:34  a_sourie
 * Cleaning of tabulated 2-fluid EoS class + superfluid rotating star model.
 *
 * Revision 1.16  2016/12/05 16:17:53  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.15  2015/06/11 13:50:19  j_novak
 * Minor corrections
 *
 * Revision 1.14  2015/06/10 14:39:17  a_sourie
 * New class Eos_bf_tabul for tabulated 2-fluid EoSs and associated functions for the computation of rotating stars with such EoSs.
 *
 * Revision 1.13  2014/10/13 08:52:57  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.12  2004/03/25 10:29:04  j_novak
 * All LORENE's units are now defined in the namespace Unites (in file unites.h).
 *
 * Revision 1.11  2003/12/04 14:28:26  r_prix
 * allow for the case of "slow-rot-style" EOS inversion, in which we need to adapt
 * the inner domain to n_outer=0 instead of mu_outer=0 ...
 * (this should only be used for comparison to analytic slow-rot solution!)
 *
 * Revision 1.10  2003/11/20 14:01:26  r_prix
 * changed member names to better conform to Lorene coding standards:
 * J_euler -> j_euler, EpS_euler -> enerps_euler, Delta_car -> delta_car
 *
 * Revision 1.9  2003/11/18 18:38:11  r_prix
 * use of new member EpS_euler: matter sources in equilibrium() and global quantities
 * no longer distinguish Newtonian/relativistic, as all terms should have the right limit...
 *
 * Revision 1.8  2003/11/17 13:49:43  r_prix
 * - moved superluminal check into hydro_euler()
 * - removed some warnings
 *
 * Revision 1.7  2003/11/13 12:07:57  r_prix
 * *) changed xxx2 -> Delta_car
 * *) added (non 2-fluid specific!) members sphph_euler J_euler
 * *) more or less rewritten hydro_euler() to see if I understand it ;)
 *   - somewhat simplified and more adapted to the notation used in our notes/paper.
 *   - Main difference: u_euler is no longer used!!, the "output" instead
 *     consists of ener_euler, s_euler, sphph_euler and J_euler, which are
 *     the general 3+1 components for Tmunu.
 *
 * Revision 1.6  2003/09/17 08:27:50  j_novak
 * New methods: mass_b1() and mass_b2().
 *
 * Revision 1.5  2002/10/18 08:42:58  j_novak
 * Take into account the sign for uuu and uuu2
 *
 * Revision 1.4  2002/01/16 15:03:28  j_novak
 * *** empty log message ***
 *
 * Revision 1.3  2002/01/11 14:09:34  j_novak
 * Added newtonian version for 2-fluid stars
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
 * Revision 1.3  2001/08/28  16:04:22  novak
 * Use of new definition of relative velocity and new declarations for EOS
 *
 * Revision 1.2  2001/08/27 09:58:43  novak
 * *** empty log message ***
 *
 * Revision 1.1  2001/06/22 15:39:17  novak
 * Initial revision
 *
 *
 * $Header$
 *
 */
// Headers C
#include "math.h"

// Headers Lorene
#include "et_rot_bifluid.h"
#include "nbr_spx.h"
#include "utilitaires.h"
#include "unites.h"	    
 
			    //--------------//
			    // Constructors //
			    //--------------//
// Standard constructor
// --------------------
namespace Lorene {
Et_rot_bifluid::Et_rot_bifluid(Map& mpi, int nzet_i, bool relat, const Eos_bifluid& eos_i):
  Etoile_rot(mpi, nzet_i, relat, *eos_i.trans2Eos()), 
  eos(eos_i),
  ent2(mpi),
  nbar2(mpi),
  K_nn(mpi),
  K_np(mpi),
  K_pp(mpi),
  alpha_eos(mpi),
  sphph_euler(mpi),
  j_euler(mpi, 1, CON, mp.get_bvect_cart()), 
  j_euler1 (mpi, 1, CON, mp.get_bvect_cart()), 
  j_euler2(mpi, 1, CON, mp.get_bvect_cart()),  
  enerps_euler(mpi),
  uuu2(mpi),
  gam_euler2(mpi),
  delta_car(mpi)
{
  // All the matter quantities are initialized to zero :
  nbar2 = 0 ;
  ent2 = 0 ; 
  K_nn = 0 ;
  K_np = 0 ;
  K_pp = 0 ; 
  alpha_eos = 0.;
  sphph_euler = 0;
  j_euler = 0;
  j_euler1 = 0 ; 
  j_euler2 = 0; 
  enerps_euler = 0;
  gam_euler.set_std_base() ; 
  
  // Initialization to a static state : 
  omega2 = 0 ; 
  uuu2 = 0 ; 
  gam_euler2 = 1 ; 
  delta_car = 0 ;
  
  // Pointers of derived quantities initialized to zero : 
  set_der_0x0() ;
  
}

// Copy constructor
// ----------------

Et_rot_bifluid::Et_rot_bifluid(const Et_rot_bifluid& et):
  Etoile_rot(et), 
  eos(et.eos),
  ent2(et.ent2),
  nbar2(et.nbar2),
  K_nn(et.K_nn),
  K_np(et.K_np),
  K_pp(et.K_pp),
  alpha_eos(et.alpha_eos),
  sphph_euler(et.sphph_euler),
  j_euler(et.j_euler),
  j_euler1(et.j_euler1),
  j_euler2(et.j_euler2), 
  enerps_euler(et.enerps_euler),
  uuu2(et.uuu2),
  gam_euler2(et.gam_euler2),
  delta_car(et.delta_car)
{
  omega2 = et.omega2 ; 

  // Pointers of derived quantities initialized to zero : 
  set_der_0x0() ;
}    


// Constructor from a file 
// ------------------------
Et_rot_bifluid::Et_rot_bifluid(Map& mpi, const Eos_bifluid& eos_i, FILE* fich):
  Etoile_rot(mpi, *eos_i.trans2Eos(), fich),
  eos(eos_i),
  ent2(mpi),
  nbar2(mpi),
  K_nn(mpi),
  K_np(mpi),
  K_pp(mpi),
  alpha_eos(mpi),
  sphph_euler(mpi),
  j_euler(mpi, 1, CON, mp.get_bvect_cart()), 
  j_euler1(mpi, 1, CON, mp.get_bvect_cart()),  
  j_euler2(mpi, 1, CON, mp.get_bvect_cart()),  
  enerps_euler(mpi),
  uuu2(mpi),
  gam_euler2(mpi),
  delta_car(mpi)
{

  // Etoile parameters
  // -----------------
  // omega2 is read in the file:     
  fread_be(&omega2, sizeof(double), 1, fich) ;		
  
  
  // Read of the saved fields:
  // ------------------------
  
  Tenseur ent2_file(mp, fich) ; 
  ent2 = ent2_file ; 
        
  // All other fields are initialized to zero : 
  // ----------------------------------------
  uuu2 = 0 ;
  delta_car = 0 ;

  // Pointers of derived quantities initialized to zero 
  // --------------------------------------------------
  set_der_0x0() ;
  
}

			    //------------//
			    // Destructor //
			    //------------//

Et_rot_bifluid::~Et_rot_bifluid(){

  del_deriv() ; 

}

		//----------------------------------//
		// Management of derived quantities //
		//----------------------------------//

void Et_rot_bifluid::del_deriv() const {

  Etoile_rot::del_deriv() ; 
  
  if (p_ray_eq2 != 0x0) delete p_ray_eq2 ; 
  if (p_ray_eq2_pis2 != 0x0) delete p_ray_eq2_pis2 ; 
  if (p_ray_eq2_pi != 0x0) delete p_ray_eq2_pi ; 
  if (p_ray_pole2 != 0x0) delete p_ray_pole2 ; 
  if (p_l_surf2 != 0x0) delete p_l_surf2 ; 
  if (p_xi_surf2 != 0x0) delete p_xi_surf2 ;
  if (p_r_circ2 != 0x0) delete p_r_circ2 ;
  if (p_area2 != 0x0) delete p_area2 ;
  if (p_aplat2 != 0x0) delete p_aplat2 ; 
  if (p_mass_b1 != 0x0) delete p_mass_b1 ;
  if (p_mass_b2 != 0x0) delete p_mass_b2 ;
  if (p_angu_mom_1 != 0x0) delete p_angu_mom_1 ; 
  if (p_angu_mom_2 != 0x0) delete p_angu_mom_2 ; 
  if (p_coupling_mominert_1 != 0x0) delete p_coupling_mominert_1 ;
  if (p_coupling_mominert_2 != 0x0) delete p_coupling_mominert_2 ; 
  if (p_coupling_entr != 0x0) delete p_coupling_entr ;
  if (p_coupling_LT_1 != 0x0) delete p_coupling_LT_1 ;
  if (p_coupling_LT_2 != 0x0) delete p_coupling_LT_2 ;
   
  set_der_0x0() ; 
}			    




void Et_rot_bifluid::set_der_0x0() const {

  Etoile_rot::set_der_0x0() ;
  
  p_ray_eq2 = 0x0 ;
  p_ray_eq2_pis2 = 0x0 ; 
  p_ray_eq2_pi = 0x0 ; 
  p_ray_pole2 = 0x0 ; 
  p_l_surf2 = 0x0 ; 
  p_xi_surf2 = 0x0 ; 
  p_r_circ2 = 0x0 ;
  p_area2 = 0x0 ;
  p_aplat2 = 0x0 ;
  p_mass_b1 = 0x0;
  p_mass_b2 = 0x0;
  p_angu_mom_1 = 0x0; 
  p_angu_mom_2 = 0x0; 
  p_coupling_mominert_1 = 0x0;
  p_coupling_mominert_2 = 0x0;
  p_coupling_entr = 0x0;
  p_coupling_LT_1 = 0x0;
  p_coupling_LT_2 = 0x0;

}			    

void Et_rot_bifluid::del_hydro_euler() {

  Etoile_rot::del_hydro_euler() ; 
  sphph_euler.set_etat_nondef();
  j_euler.set_etat_nondef();
  j_euler1.set_etat_nondef(); 
  j_euler2.set_etat_nondef();  
  enerps_euler.set_etat_nondef();
  uuu2.set_etat_nondef();
  gam_euler2.set_etat_nondef() ; 
  delta_car.set_etat_nondef();
  K_nn.set_etat_nondef(); 
  K_np.set_etat_nondef();
  K_pp.set_etat_nondef();
  alpha_eos.set_etat_nondef() ;
  del_deriv() ; 
}			    

// Assignment of the enthalpy field
// --------------------------------

void Et_rot_bifluid::set_enthalpies(const Cmp& ent_i, const Cmp& ent2_i) {
    
  ent = ent_i ; 
  ent2 = ent2_i ;
    
  // Update of (nbar, ener, press) :
  equation_of_state() ; 
    
  // The derived quantities are obsolete:
  del_deriv() ; 
    
}


			    //--------------//
			    //  Assignment  //
			    //--------------//

// Assignment to another Et_rot_bifluid
// --------------------------------
void Et_rot_bifluid::operator=(const Et_rot_bifluid& et) {

    // Assignment of quantities common to all the derived classes of Etoile
    Etoile_rot::operator=(et) ;	    

    assert( &(et.eos) == &eos ) ;	    // Same EOS
    // Assignement of proper quantities of class Et_rot_bifluid
    omega2 = et.omega2 ; 

    ent2 = et.ent2 ;
    nbar2 = et.nbar2 ;
    K_nn = et.K_nn ;
    K_np = et.K_np ;
    K_pp = et.K_pp ;
    alpha_eos = et.alpha_eos ;
    sphph_euler = et.sphph_euler;
    j_euler = et.j_euler;
    j_euler1 = et.j_euler1; 
    j_euler2 = et.j_euler2; 
    enerps_euler = et.enerps_euler;
    uuu2 = et.uuu2 ;
    gam_euler2 = et.gam_euler2 ;
    delta_car = et.delta_car ;

    del_deriv() ;  // Deletes all derived quantities

}	

			    //--------------//
			    //	  Outputs   //
			    //--------------//

// Save in a file
// --------------
void Et_rot_bifluid::sauve(FILE* fich) const {
    
    Etoile_rot::sauve(fich) ; 
    
    fwrite_be(&omega2, sizeof(double), 1, fich) ;		
    
    ent2.sauve(fich) ; 
    
    
}

// Printing
// --------

ostream& Et_rot_bifluid::operator>>(ostream& ost) const {
    
  using namespace Unites ;

    Etoile::operator>>(ost) ; 
    
    ost << endl ; 
    ost << "Bifluid rotating star" << endl ; 
    ost << "-------------" << endl ; 
    ost << setprecision(16);
    double freq = omega / (2.*M_PI) ;  
    ost << "Omega1 : " << omega * f_unit 
        << " rad/s     f : " << freq * f_unit << " Hz" << endl ; 
    ost << "Rotation period 1: " << 1000. / (freq * f_unit) << " ms"
	    << endl ;
       
    double freq2 = omega2 / (2.*M_PI) ;  
    ost << "Omega2 : " << omega2 * f_unit 
        << " rad/s     f : " << freq2 * f_unit << " Hz" << endl ; 
    ost << "Rotation period 2: " << 1000. / (freq2 * f_unit) << " ms"
	    << endl ;
      

	ost << "Total angular momentum J :      " 
	 << angu_mom()/( qpig / (4* M_PI) * msol*msol) << " G M_sol^2 / c"
	 << endl ; 
   	ost << "c J / (G M^2) :           " 
	 << angu_mom()/( qpig / (4* M_PI) * pow(mass_g(), 2.) ) << endl ;   
	
	double mom_iner = fabs(angu_mom() / omega2) ;
	ost <<  "Total moment of inertia I = J/Omega2 :      "  
	    <<  mom_iner << " Lorene units" 
	    << endl; 
	double mom_iner_38si = mom_iner * rho_unit * (pow(r_unit, double(5.)) / double(1.e38) ) ; 	 
	ost << "Total moment of inertia I = J/Omega2 :      "  
	<< mom_iner_38si << " 10^38 kg m^2"
	    << endl ; 
	
   // partial angular momenta
	ost << "Angular momentum J of fluid 1 :      " 
	 << angu_mom_1()/( qpig / (4* M_PI) * msol*msol) << " G M_sol^2 / c"
	 << endl ; 	
	ost << "Angular momentum J of fluid 2 :      " 
	 << angu_mom_2()/( qpig / (4* M_PI) * msol*msol) << " G M_sol^2 / c"
	 << endl ; 


 	//	partial moments of inertia defined as Jn/Omega and Jp/Omega
	// Remark : such definitions only makes sense in corotation
	ost.precision(16) ;     
	if (omega == omega2) {
	
		double mom_iner_1 = 0.; //< In
		double mom_iner_2 = 0.; //< Ip
  
      if (omega != __infinity) { 
			
			ost << "Partial moments of inertia (defined in corotation only) :" << endl;	

			mom_iner_1 = fabs(angu_mom_1() / omega) ; 
			ost 	<< "Moment of inertia of fluid 1 :       " << mom_iner_1 << " Lorene units " << endl ;
			double mom_iner_1_38si = mom_iner_1 * rho_unit * (pow(r_unit, double(5.)) / double(1.e38) ) ; 
			ost 	<< "Moment of inertia of fluid 1 :       " << mom_iner_1_38si  << " 10^38 kg m^2"
	    			<< endl ; 

			mom_iner_2 = fabs(angu_mom_2() / omega) ; 
			ost 	<< "Moment of inertia of fluid 2 :       " << mom_iner_2 << " Lorene units " << endl ;
			double mom_iner_2_38si = mom_iner_2 * rho_unit * (pow(r_unit, double(5.)) / double(1.e38) ) ; 
 			ost 	<< "Moment of inertia of fluid 2 :       " << mom_iner_2_38si << " 10^38 kg m^2"
	    			<< endl;	    
	    }
	} 

	ost << "***** Fluid coupling quantities *****" << endl;
	
	double mominert_1_tilde_si = coupling_mominert_1() * rho_unit * (pow(r_unit, double(5.)) / double(1.e38) ) ; 
	double mominert_2_tilde_si = coupling_mominert_2() * rho_unit * (pow(r_unit, double(5.)) / double(1.e38) ) ; 

	ost << "tilde{I}_n : " << coupling_mominert_1() << "    SI: " << mominert_1_tilde_si << " 10^38 kg m^2" 
	    << endl;
	ost << "tilde{I}_p : " << coupling_mominert_2() << "    SI: " << mominert_2_tilde_si << " 10^38 kg m^2" 
	    << endl;	
	ost << "tilde{varepsilon}_n : " << coupling_entr() / coupling_mominert_1() << endl;
	ost << "tilde{varepsilon}_p : " << coupling_entr() / coupling_mominert_2() << endl;
	ost << "tilde{omega}_n : " << coupling_LT_1() / coupling_mominert_1() << endl;
	ost << "tilde{omega}_p : " << coupling_LT_2() / coupling_mominert_2() << endl;
	ost << " Verif :  Jn = " << coupling_mominert_1() * omega - coupling_LT_1() +  coupling_entr() * (omega2 - omega) 
	    << "         Jp = "  << coupling_mominert_2() * omega2 - coupling_LT_2() +  coupling_entr() * (omega - omega2) 
	    << endl ;
	ost << " Num :    Jn = " <<  angu_mom_1() << "         Jp = " << angu_mom_2()
	    << endl;
	
    double nphi_c = nphi()(0, 0, 0, 0) ;
    if ( (omega==0) && (nphi_c==0) ) {
	 	ost << "Central N^phi :               " << nphi_c << endl ;
    }
    else{
		ost << "Central N^phi/Omega :    " << nphi_c / omega << endl ;
    }
    if (omega2!=0) 
      ost << "Central N^phi/Omega2 :     " << nphi_c / omega2 << endl ;
    
    ost << "Error on the virial identity GRV2 : " << endl ; 
    ost << "GRV2 = " << grv2() << endl ; 
    ost << "Error on the virial identity GRV3 : " << endl ; 
    double xgrv3 = grv3(&ost) ; 
    ost << "GRV3 = " << xgrv3 << endl ; 

    ost << "Circumferential equatorial radius R_circ :     " 
	<< r_circ()/km << " km" << endl ;  
    ost << "Mean radius R_mean :    " 
	<< mean_radius()/km << " km" << endl ;
    ost << "Coordinate equatorial radius r_eq : " << ray_eq()/km << " km" 
	 << endl ;  
    ost << "Flattening r_pole/r_eq :        " << aplat() << endl ; 
    ost << "Circumferential equatorial radius R_circ2 :     " 
	<< r_circ2()/km << " km" << endl ;  
    ost << "Mean radius R_mean2 :    " 
	<< mean_radius2()/km << " km" << endl ;	
    ost << "Coordinate equatorial radius r_eq2 : " << ray_eq2()/km << " km" 
	 << endl ;  
    ost << "Flattening r_pole2/r_eq2 :        " << aplat2() << endl ; 

    int lsurf = nzet - 1; 
    int nt = mp.get_mg()->get_nt(lsurf) ; 
    int nr = mp.get_mg()->get_nr(lsurf) ; 
    ost << "Equatorial value of the velocity U:         " 
	 << uuu()(nzet-1, 0, nt-1, nr-1) << " c" << endl ; 
    ost << "Equatorial value of the velocity U2:         " 
	 << uuu2()(nzet-1, 0, nt-1, nr-1) << " c" << endl ; 
    ost << "Redshift at the equator (forward) : " << z_eqf() << endl ; 
    ost << "Redshift at the equator (backward): " << z_eqb() << endl ; 
    ost << "Redshift at the pole              : " << z_pole() << endl ; 


    ost << "Central value of log(N)        : " 
	<< logn()(0, 0, 0, 0) << endl ; 

    ost << "Central value of dzeta=log(AN) : " 
	<< dzeta()(0, 0, 0, 0) << endl ; 

    ost << "Central value of B^2 : " << b_car()(0,0,0,0) <<  endl ; 

    Tbl diff_a_b = diffrel( a_car(), b_car() ) ;
    ost << 
    "Relative discrepancy in each domain between the metric coef. A^2 and B^2 : "
       << endl ;
    for (int l=0; l<diff_a_b.get_taille(); l++) {
    	ost << diff_a_b(l) << "  " ;
    }
    ost << endl;
     ost << "Quadrupole moment  : " <<  mom_quad() << endl ;
    double mom_quad_38si = mom_quad() * rho_unit * (pow(r_unit, double(5.))  / double(1.e38) ) ;
     ost << "Quadrupole moment Q : " << mom_quad_38si << " 10^38 kg m^2"
 	<< endl ; 
    ost << "Old quadrupole moment  : " << mom_quad_old() << endl ;
    ost << "Coefficient b  : " << mom_quad_Bo() /  pow(mass_g(), 2.) << endl ;
    ost << "q = c^4 Q / (G^2 M^3) :  " 
        <<  mom_quad() / ( ggrav * ggrav *  pow(mass_g(), 3.) )  << endl ;
    ost << "j = c J / (G M^2) :  "         
	<< angu_mom()/( ggrav * pow(mass_g(), 2.) ) << endl ;   
    

    ost << "Baryon mass 1  : " << mass_b1() / msol << "  Msol" <<  endl ;
    ost << "Baryon mass 2  : " << mass_b2() / msol << "  Msol" <<  endl ;

    ost << endl ;    	

    return ost ;
    
}


void Et_rot_bifluid::partial_display(ostream& ost) const {
    
  using namespace Unites ;

    double freq = omega / (2.*M_PI) ;  
    ost << "Omega : " << omega * f_unit 
        << " rad/s     f : " << freq * f_unit << " Hz" << endl ; 
    ost << "Rotation period : " << 1000. / (freq * f_unit) << " ms"
	 << endl ;
    ost << endl << "Central enthalpy : " << ent()(0,0,0,0) << " c^2" << endl ; 
    ost << "Central proper baryon density : " << nbar()(0,0,0,0) 
	<< " x 0.1 fm^-3" << endl ; 
    double freq2 = omega2 / (2.*M_PI) ;  
    ost << "Omega2 : " << omega2 * f_unit 
        << " rad/s     f : " << freq2 * f_unit << " Hz" << endl ; 
    ost << "Rotation period 2: " << 1000. / (freq2 * f_unit) << " ms"
	 << endl ;
    ost << endl << "Central enthalpy 2: " << ent2()(0,0,0,0) << " c^2" << endl ; 
    ost << "Central proper baryon density 2: " << nbar2()(0,0,0,0) 
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
    if ( (omega==0) && (nphi_c==0) ) {
		ost << "Central N^phi :               " << nphi_c << endl ;
    }
    else{
		ost << "Central N^phi/Omega :         " << nphi_c / omega << endl ;
    }
	    

    int lsurf = nzet - 1; 
    int nt = mp.get_mg()->get_nt(lsurf) ; 
    int nr = mp.get_mg()->get_nr(lsurf) ; 
    ost << "Equatorial value of the velocity U:         " 
	 << uuu()(nzet-1, 0, nt-1, nr-1) << " c" << endl ; 

    ost << "Equatorial value of the velocity U2:         " 
	 << uuu2()(nzet-1, 0, nt-1, nr-1) << " c" << endl ; 

    ost << endl 
	<< "Coordinate equatorial radius r_eq =    " 
	<< ray_eq()/km << " km" << endl ;  
    ost << "Flattening r_pole/r_eq :        " << aplat() << endl ; 
    ost << endl 
	<< "Coordinate equatorial radius r_eq2 =    " 
	<< ray_eq2()/km << " km" << endl ;  
    ost << "Flattening r_pole2/r_eq2 :        " << aplat2() << endl ; 

}

//
//   Equation of state
//

void Et_rot_bifluid::equation_of_state() {

  Cmp ent_eos = ent() ;
  Cmp ent2_eos = ent2() ;
  Tenseur rel_vel(delta_car) ;

  if (nzet > 1) {
    // Slight rescale of the enthalpy field in case of 2 domains inside the
    //  star
  
    if (nzet > 2) {
      cout << "Et_rot_bifluid::equation_of_state: not ready yet for nzet > 2 !" << endl ;    	
    }
    
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
  
    Cmp fact_ent(mp) ;
    fact_ent.allocate_all() ;
    
    fact_ent.set(0) = 1 + epsilon * xi(0) * xi(0) ;
    fact_ent.set(1) = 1 - 0.25 * epsilon * (xi(1) - 1) * (xi(1) - 1) ;
  
    for (int l=nzet; l<nz; l++) {
      fact_ent.set(l) = 1 ;
    }

    ent_eos = fact_ent * ent_eos ;
    ent2_eos = fact_ent * ent2_eos ;
    ent_eos.std_base_scal() ;
    ent2_eos.std_base_scal() ;

  } // if nzet > 1
  
  
  // Call to the EOS
  nbar.set_etat_qcq() ;
  nbar2.set_etat_qcq() ;
  ener.set_etat_qcq() ;
  press.set_etat_qcq() ;
  K_nn.set_etat_qcq() ;
  K_np.set_etat_qcq() ;
  K_pp.set_etat_qcq() ;
  alpha_eos.set_etat_qcq() ;
 
  const Eos_bf_tabul* eos_t = dynamic_cast<const Eos_bf_tabul*>(&eos) ;
  
  if (eos_t != 0x0) {  	// The EoS is tabulated
  		eos_t->calcule_interpol(ent_eos, ent2_eos, rel_vel(), nbar.set(), nbar2.set(), 
  		   ener.set(), press.set(), K_nn.set(), K_np.set(), K_pp.set(), alpha_eos.set(), 
		   nzet)  ; 
  }
  else {	      			// The EoS is analytic
  		eos.calcule_tout(ent_eos, ent2_eos, rel_vel(), nbar.set(), nbar2.set(), 
  		   ener.set(), press.set(), nzet)  ; 	  
		   
 		K_nn.set() = eos.get_Knn(nbar(), nbar2(), delta_car(), nzet);
		K_pp.set() = eos.get_Kpp(nbar(), nbar2(), delta_car(), nzet);
		K_np.set() = eos.get_Knp(nbar(), nbar2(), delta_car(), nzet);
		alpha_eos.set() = eos.get_Knp(nbar(), nbar2(), delta_car(), nzet) * nbar() * nbar2() * pow(1. - unsurc2 * rel_vel(), -1.5) / 2. ; 
  }
  // Set the bases for spectral expansion 
  nbar.set_std_base() ; 
  nbar2.set_std_base() ; 
  ener.set_std_base() ; 
  press.set_std_base() ; 
  K_pp.set_std_base() ; 
  K_nn.set_std_base() ; 
  K_np.set_std_base() ; 
  alpha_eos.set_std_base() ;
  // The derived quantities are obsolete
  del_deriv() ; 
  
}

//
// Computation of hydro quantities
//

void Et_rot_bifluid::hydro_euler(){

  const Mg3d* mg = mp.get_mg(); 
  int nz = mg->get_nzone() ; 
  int nzm1 = nz - 1 ; 

    // RP: I prefer to use the 3-vector j_euler instead of u_euler
    // for better physical "encapsulation" 
    // (i.e. --> use same form of Poisson-equations for all etoile sub-classes!)
    u_euler.set_etat_nondef(); // make sure it's not used

    // (Norm of) Euler-velocity of the first fluid
    //------------------------------
    uuu.set_etat_qcq();

    uuu.set() = bbb() * (omega - nphi() ) / nnn();
    uuu.annule(nzm1) ; 

    // gosh, we have to exclude the thing being zero here... :(
    if( uuu.get_etat() != ETATZERO )
      {
	(uuu.set()).std_base_scal() ;
	(uuu.set()).mult_rsint();
      }
    uuu.set_std_base();
    

    // (Norm of) Euler-velocity of the second fluid
    //----------------------------------------
    uuu2.set_etat_qcq();
    
    uuu2.set() = bbb() * (omega2 - nphi() ) / nnn();
    uuu2.annule(nzm1) ; 

    if( uuu2.get_etat() != ETATZERO )
      {
	(uuu2.set()).std_base_scal();
	(uuu2.set()).mult_rsint();
      }
    uuu2.set_std_base();

    // Sanity check:
    // Is one of the new velocities larger than c in the equatorial plane ?
    //----------------------------------------

    // Index of the point at phi=0, theta=pi/2 at the surface of the star:
    int l_b = nzet - 1 ; 
    int j_b = mg->get_nt(l_b) - 1 ; 

    bool superlum = false ; 
    if (relativistic) {
      for (int l=0; l<nzet; l++) {
	for (int i=0; i<mg->get_nr(l); i++) {
	  
	  double u1 = uuu()(l, 0, j_b, i) ; 
	  double u2 = uuu2()(l, 0, j_b, i) ;
	  if ((u1 >= 1.) || (u2>=1.)) {	    // superluminal velocity
	    superlum = true ; 
	    cout << "U > c  for l, i : " << l << "  " << i 
		 << "   U1 = " << u1 << endl ;
	    cout << "   U2 = " << u2 << endl ;
	  }
	}
      }
      if ( superlum ) {
	cout << "**** VELOCITY OF LIGHT REACHED ****" << endl ; 
	abort() ;
      }
    }


    Tenseur uuu_car = uuu * uuu;
    Tenseur uuu2_car = uuu2 * uuu2;


    // Lorentz factors
    // --------------
    Tenseur gam_car = 1.0 / (1.0 - unsurc2 * uuu_car) ;
    gam_euler = sqrt(gam_car) ;
    gam_euler.set_std_base() ;  // sets the standard spectral bases for a scalar field


    Tenseur gam2_car = 1.0 / (1.0 - unsurc2 * uuu2_car) ;
    gam_euler2 = sqrt(gam2_car) ;
    gam_euler2.set_std_base() ;

    // Update of "relative velocity" squared: $\Delta^2$
    // ---------------------------

    delta_car =  (uuu - uuu2)*(uuu - uuu2) / ( (1 - unsurc2* uuu*uuu2) *(1 - unsurc2* uuu*uuu2 ) ) ;
    delta_car.set_std_base() ;
 
 	 Tenseur Ann(ent) ;
	 Tenseur Anp(ent) ;
  	 Tenseur App(ent) ;

	 Ann = gam_car() * nbar() * nbar() * K_nn() ;   			
	 Anp = gam_euler() * gam_euler2() * nbar() * nbar2() * K_np();	
	 App = gam2_car() * nbar2() * nbar2() * K_pp();		
	
   
    //  Energy density E with respect to the Eulerian observer
    //------------------------------------
    // use of ener_euler is deprecated, because it's useless in Newtonian limit!
    ener_euler.set_etat_nondef(); // make sure, it's not used

    Tenseur E_euler(mp);
    E_euler =  Ann + 2. * Anp + App - press ;
    E_euler.set_std_base() ; 
    

    // S^phi_phi component of stress-tensor S^i_j
    //------------------------------------
    sphph_euler = press() + Ann() * uuu_car() + 2. * Anp() * uuu() * uuu2() + App() * uuu2_car();
    sphph_euler.set_std_base(); 


    // Trace of the stress tensor with respect to the Eulerian observer
    //------------------------------------
    s_euler = 2. * press() + sphph_euler();
    s_euler.set_std_base() ; 

    // The combination enerps_euler := (E + S_i^i) which has Newtonian limit -> rho
    if (relativistic)
      enerps_euler = E_euler + s_euler;
    else
      enerps_euler = eos.get_m1() * nbar() + eos.get_m2() * nbar2();


    // the (flat-space) angular-momentum 3-vector j_euler^i
    //-----------------------------------
    Tenseur Jph(mp);   // the normalized phi-component of J^i: Sqrt[g_phiphi]*J^phi
    Jph = Ann*uuu + Anp*(uuu + uuu2) + App*uuu2 ;

    j_euler.set_etat_qcq();

    j_euler.set(0) = 0;						// r tetrad component
    j_euler.set(1) = 0;						// theta tetrad component
    j_euler.set(2) = Jph()/ bbb(); 		// phi tetrad component ... = J^phi r sin(th)
    j_euler.set_triad (mp.get_bvect_spher());
    j_euler.set_std_base();

    // RP: it seems that j_euler _HAS_ to have cartesian triad set on exit from here...!!
    j_euler.change_triad( mp.get_bvect_cart() ) ;	// Triad = Cartesian triad

    if( (j_euler(0).get_etat() == ETATZERO)&&(j_euler(1).get_etat() == ETATZERO)&&(j_euler(2).get_etat()==ETATZERO))
      j_euler = 0;

    // the (flat-space) angular-momentum 3-vector j_euler^i for fluid 1 
    //-----------------------------------
    Tenseur Jph1(mp);   // the normalized phi-component of J_n^i: Sqrt[g_phiphi]*J_n^phi
    Jph1 = Ann*uuu + Anp*uuu2 ;

    j_euler1.set_etat_qcq();

    j_euler1.set(0) = 0;					 
    j_euler1.set(1) = 0;					  
    j_euler1.set(2) = Jph1()/ bbb(); 	 
    j_euler1.set_triad (mp.get_bvect_spher());
    j_euler1.set_std_base();

    j_euler1.change_triad( mp.get_bvect_cart() ) ;	 

    if( (j_euler1(0).get_etat() == ETATZERO)&&(j_euler1(1).get_etat() == ETATZERO)&&(j_euler1(2).get_etat()==ETATZERO))
      j_euler1 = 0;

    // the (flat-space) angular-momentum 3-vector j_euler^i for fluid 2 
    //-----------------------------------
    Tenseur Jph2(mp);   // the normalized phi-component of J_p^i: Sqrt[g_phiphi]*J_p^phi
    Jph2 = Anp*uuu + App*uuu2 ;

    j_euler2.set_etat_qcq();

    j_euler2.set(0) = 0;			 
    j_euler2.set(1) = 0;		 
    j_euler2.set(2) = Jph2()/ bbb(); 		 
    j_euler2.set_triad (mp.get_bvect_spher());
    j_euler2.set_std_base();

     j_euler2.change_triad( mp.get_bvect_cart() ) ;	// Triad = Cartesian triad

    if( (j_euler2(0).get_etat() == ETATZERO)&&(j_euler2(1).get_etat() == ETATZERO)&&(j_euler2(2).get_etat()==ETATZERO))
      j_euler2 = 0;
    
    // The derived quantities are obsolete
    // -----------------------------------
    del_deriv() ;                

} // hydro_euler()

}

