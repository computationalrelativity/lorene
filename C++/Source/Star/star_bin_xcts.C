 /*
 * Methods for the class Star_bin_xcts
 * (see file star.h for documentation)
 */

/*
 *   Copyright (c) 2010 Michal Bejger
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

char star_bin_xcts_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.8  2010/12/09 10:45:26  m_bejger
 * Added psi4 in the constructor, removed fait_decouple
 *
 * Revision 1.7  2010/10/28 13:50:27  m_bejger
 * Added mass-shedding estimation to Star_bin_xcts::operator>>
 *
 * Revision 1.6  2010/10/26 20:18:34  m_bejger
 * Correction to stdin output
 *
 * Revision 1.5  2010/10/26 19:57:02  m_bejger
 * Various cleanups
 *
 * Revision 1.4  2010/06/17 15:06:25  m_bejger
 * N, Psi output corrected
 *
 * Revision 1.3  2010/06/15 08:18:19  m_bejger
 * Method set_chi_comp() added
 *
 * Revision 1.2  2010/06/04 19:59:56  m_bejger
 * Corrected definitions of lapse and Psi4
 *
 * Revision 1.1  2010/05/04 07:51:05  m_bejger
 * Initial version
 *
 * $Header$
 *
 */

// Headers C
#include "math.h"

// Headers Lorene
#include "etoile.h"
#include "star.h"
#include "eos.h"
#include "unites.h"

// Local prototype
Cmp raccord_c1(const Cmp& uu, int l1) ;

			    //--------------//
			    // Constructors //
			    //--------------//

// Standard constructor
// --------------------
Star_bin_xcts::Star_bin_xcts(Map& mpi,
							int nzet_i,
							const Eos& eos_i,
		     				bool irrot)
		: Star(mpi, nzet_i, eos_i),
		irrotational(irrot),
      	psi0(mpi),
      	d_psi(mpi, COV, mpi.get_bvect_cart()),
      	wit_w(mpi, CON, mpi.get_bvect_cart()),
      	loggam(mpi),
      	bsn(mpi, CON, mpi.get_bvect_cart()),
      	pot_centri(mpi),
      	chi_auto(mpi),
      	chi_comp(mpi),
      	Psi_auto(mpi),
      	Psi_comp(mpi),
      	Psi(mpi),
      	chi(mpi),
      	psi4(mpi),
      	w_beta(mpi, CON, mpi.get_bvect_cart()),
      	khi(mpi),
      	dcov_Psi(mpi, COV, mpi.get_bvect_cart()),
      	dcov_chi(mpi, COV, mpi.get_bvect_cart()),
      	flat(mpi, mpi.get_bvect_cart()),
      	beta_auto(mpi, CON, mpi.get_bvect_cart()),
      	beta_comp(mpi, CON, mpi.get_bvect_cart()),
      	haij_auto(mpi, CON, mpi.get_bvect_cart()),
      	haij_comp(mpi, CON, mpi.get_bvect_cart()),
      	hacar_auto(mpi),
      	hacar_comp(mpi),
      	ssjm1_chi(mpi),
      	ssjm1_psi(mpi),
      	ssjm1_khi(mpi),
      	ssjm1_wbeta(mpi, CON, mpi.get_bvect_cart()) {

    // Pointers of derived quantities initialized to zero :
    set_der_0x0() ;

    // Quantities defined on a spherical triad in star.C are put on
    // a cartesian one
    u_euler.change_triad(mpi.get_bvect_cart()) ;
    stress_euler.change_triad(mpi.get_bvect_cart()) ;
    beta.change_triad(mpi.get_bvect_cart()) ;
    Metric temp_met (mp.flat_met_cart()) ;
    gamma = temp_met ;

    // Initialization of all quantities :
    psi0 = 0 ;
    d_psi.set_etat_zero() ;
    wit_w.set_etat_zero() ;
    loggam = 0 ;
    bsn.set_etat_zero() ;
    pot_centri = 0 ;

    Psi_auto = 0 ;
    Psi_comp = 0 ;

    beta_auto.set_etat_zero() ;
    beta_comp.set_etat_zero() ;
    chi_auto = 0 ;
    chi_comp = 0 ;

    Psi = 1. ;
    chi = 1. ;

    w_beta.set_etat_zero() ;
    khi = 0 ;

    dcov_Psi.set_etat_zero() ;
    dcov_chi.set_etat_zero() ;

    haij_auto.set_etat_zero() ;
    haij_comp.set_etat_zero() ;
    hacar_auto = 0 ;
    hacar_comp = 0 ;
    ssjm1_psi = 0 ;
    ssjm1_chi = 0 ;
    ssjm1_khi = 0 ;
    ssjm1_wbeta.set_etat_zero() ;

}

// Copy constructor
// ----------------
Star_bin_xcts::Star_bin_xcts(const Star_bin_xcts& star)
		: Star(star),
		irrotational(star.irrotational),
		psi0(star.psi0),
		d_psi(star.d_psi),
		wit_w(star.wit_w),
		loggam(star.loggam),
		bsn(star.bsn),
		pot_centri(star.pot_centri),
		chi_auto(star.chi_auto),
		chi_comp(star.chi_comp),
		Psi_auto(star.Psi_auto),
		Psi_comp(star.Psi_comp),
		Psi(star.Psi),
		chi(star.chi),
		psi4(star.psi4),
		w_beta(star.w_beta),
		khi(star.khi),
		dcov_Psi(star.dcov_Psi),
		dcov_chi(star.dcov_chi),
		flat(star.flat),
		beta_auto(star.beta_auto),
		beta_comp(star.beta_comp),
		haij_auto(star.haij_auto),
		haij_comp(star.haij_comp),
		hacar_auto(star.hacar_auto),
		hacar_comp(star.hacar_comp),
		ssjm1_chi(star.ssjm1_chi),
		ssjm1_psi(star.ssjm1_psi),
		ssjm1_khi(star.ssjm1_khi),
		ssjm1_wbeta(star.ssjm1_wbeta) {

    set_der_0x0() ;

}

// Constructor from a file
// -----------------------
Star_bin_xcts::Star_bin_xcts(Map& mpi,
							const Eos& eos_i,
							FILE* fich)
		: Star(mpi, eos_i, fich),
		psi0(mpi),
		d_psi(mpi, COV, mpi.get_bvect_cart()),
		wit_w(mpi, CON, mpi.get_bvect_cart()),
		loggam(mpi),
		bsn(mpi, CON, mpi.get_bvect_cart()),
		pot_centri(mpi),
		chi_auto(mpi, *(mpi.get_mg()), fich),
		chi_comp(mpi),
		Psi_auto(mpi, *(mpi.get_mg()), fich),
		Psi_comp(mpi),
		Psi(mpi),
		chi(mpi),
		psi4(mpi),
		w_beta(mpi, mpi.get_bvect_cart(), fich),
		khi(mpi, *(mpi.get_mg()), fich),
		dcov_Psi(mpi, COV, mpi.get_bvect_cart()),
		dcov_chi(mpi, COV, mpi.get_bvect_cart()),
		flat(mpi, mpi.get_bvect_cart()),
		beta_auto(mpi, mpi.get_bvect_cart(), fich),
		beta_comp(mpi, CON, mpi.get_bvect_cart()),
     	haij_auto(mpi, CON, mpi.get_bvect_cart()),
	    haij_comp(mpi, CON, mpi.get_bvect_cart()),
		hacar_auto(mpi),
		hacar_comp(mpi),
		ssjm1_chi(mpi, *(mpi.get_mg()), fich),
		ssjm1_psi(mpi, *(mpi.get_mg()), fich),
		ssjm1_khi(mpi, *(mpi.get_mg()), fich),
		ssjm1_wbeta(mpi, mpi.get_bvect_cart(), fich) {

    // Star parameters
    // -----------------

    // irrotational is read in the file:
    bool status = fread(&irrotational, sizeof(bool), 1, fich) ;
    if(!status)
    	cout << "Star_bin_xcts::Constructor from a file: Problem with reading ! " << endl ;

    // Read of the saved fields:
    // ------------------------

    if (irrotational) {
		Scalar psi0_file(mp, *(mp.get_mg()), fich) ;
		psi0 = psi0_file ;

		Scalar gam_euler_file(mp, *(mp.get_mg()), fich) ;
		gam_euler = gam_euler_file ;
    }

    // Quantities defined on a spherical triad in star.C are put on
    // a cartesian one

    u_euler.change_triad(mpi.get_bvect_cart()) ;
    stress_euler.change_triad(mpi.get_bvect_cart()) ;
    beta.change_triad(mpi.get_bvect_cart()) ;
    Metric temp_met (mp.flat_met_cart()) ;
    gamma = temp_met ;

    // All other fields are initialized to initial values :
    // ----------------------------------------------------

    d_psi.set_etat_zero() ;
    wit_w.set_etat_zero() ;
    loggam = 0 ;
    bsn.set_etat_zero() ;
    pot_centri = 0 ;
    Psi_comp = 0 ;
    dcov_Psi.set_etat_zero() ;
    beta_comp.set_etat_zero() ;
    chi_comp = 0 ;
    w_beta.set_etat_zero() ;
    khi = 0 ;
    dcov_chi.set_etat_zero() ;
    haij_auto.set_etat_zero() ;
    haij_comp.set_etat_zero() ;
    hacar_auto = 0 ;
    hacar_comp = 0 ;

    // Pointers of derived quantities initialized to zero
    // --------------------------------------------------
    set_der_0x0() ;

}

			    //------------//
			    // Destructor //
			    //------------//

Star_bin_xcts::~Star_bin_xcts() {

    del_deriv() ;

}

			//----------------------------------//
			// Management of derived quantities //
			//----------------------------------//

void Star_bin_xcts::del_deriv() const {

    Star::del_deriv() ;

    if (p_xa_barycenter != 0x0) delete p_xa_barycenter ;

    set_der_0x0() ;
}

void Star_bin_xcts::set_der_0x0() const {

    Star::set_der_0x0() ;

    p_xa_barycenter = 0x0 ;

}

void Star_bin_xcts::del_hydro_euler() {

    Star::del_hydro_euler() ;

    del_deriv() ;

}


			    //--------------//
			    //  Assignment  //
			    //--------------//

// Assignment to another Star_bin_xcts
// --------------------------------
void Star_bin_xcts::operator=(const Star_bin_xcts& star) {

    // Assignment of quantities common to the derived classes of Star
    Star::operator=(star) ;

    // Assignement of proper quantities of class Star_bin_xcts
    irrotational = star.irrotational ;
    psi0 = star.psi0 ;
    d_psi = star.d_psi ;
    wit_w = star.wit_w ;
    loggam = star.loggam ;
    bsn = star.bsn ;
    pot_centri = star.pot_centri ;
    Psi_auto = star.Psi_auto ;
    Psi_comp = star.Psi_comp ;
    chi_auto = star.chi_auto ;
    chi_comp = star.chi_comp ;
    Psi = star.Psi ;
    chi = star.chi ;
    psi4 = star.psi4 ;
    w_beta = star.w_beta ;
    khi = star.khi ;
    dcov_Psi = star.dcov_Psi ;
    dcov_chi = star.dcov_chi ;
    flat = star.flat ;
    beta_auto = star.beta_auto ;
    beta_comp = star.beta_comp ;
    haij_auto = star.haij_auto ;
    haij_comp = star.haij_comp ;
    hacar_auto = star.hacar_auto ;
    hacar_comp = star.hacar_comp ;
    ssjm1_psi = star.ssjm1_psi ;
    ssjm1_chi = star.ssjm1_chi ;
    ssjm1_khi = star.ssjm1_khi ;
    ssjm1_wbeta = star.ssjm1_wbeta ;

    del_deriv() ;  // Deletes all derived quantities

}

Scalar& Star_bin_xcts::set_pot_centri() {

    del_deriv() ;	// sets to 0x0 all the derived quantities
    return pot_centri ;

}

Scalar& Star_bin_xcts::set_Psi_comp() {

    del_deriv() ;	// sets to 0x0 all the derived quantities
    return Psi_comp ;

}

Scalar& Star_bin_xcts::set_Psi_auto() {

    del_deriv() ;	// sets to 0x0 all the derived quantities
    return Psi_auto ;

}

Scalar& Star_bin_xcts::set_chi_comp() {

    del_deriv() ;	// sets to 0x0 all the derived quantities
    return chi_comp ;

}

Scalar& Star_bin_xcts::set_chi_auto() {

    del_deriv() ;	// sets to 0x0 all the derived quantities
    return chi_auto ;

}

Vector& Star_bin_xcts::set_beta_auto() {

    del_deriv() ;	// sets to 0x0 all the derived quantities
    return beta_auto ;

}

Vector& Star_bin_xcts::set_beta() {

    del_deriv() ;	// sets to 0x0 all the derived quantities
   return beta ;

}


			    //--------------//
			    //	  Outputs   //
			    //--------------//

// Save in a file
// --------------
void Star_bin_xcts::sauve(FILE* fich) const {

    Star::sauve(fich) ;

    chi_auto.sauve(fich) ;
    Psi_auto.sauve(fich) ;

    w_beta.sauve(fich) ;
    khi.sauve(fich) ;

    beta_auto.sauve(fich) ;

    ssjm1_chi.sauve(fich) ;
    ssjm1_psi.sauve(fich) ;
    ssjm1_khi.sauve(fich) ;
    ssjm1_wbeta.sauve(fich) ;

    fwrite(&irrotational, sizeof(bool), 1, fich) ;

    if (irrotational) {

		psi0.sauve(fich) ;
		gam_euler.sauve(fich) ; // required to construct d_psi from psi0

   }

}

// Printing
// --------

ostream& Star_bin_xcts::operator>>(ostream& ost) const {

  using namespace Unites ;

//    Star::operator>>(ost) ;

    ost << endl ;

    ost << "Number of domains occupied by the star : " << nzet << endl ;

    ost << "Equation of state : " << endl ;
    ost << eos << endl ;

    ost << endl << "Central enthalpy : " << ent.val_grid_point(0,0,0,0) << " c^2" << endl ;
    ost << "Central proper baryon density : " << nbar.val_grid_point(0,0,0,0)
	<< " x 0.1 fm^-3" << endl ;
    ost << "Central proper energy density : " << ener.val_grid_point(0,0,0,0)
	<< " rho_nuc c^2" << endl ;
    ost << "Central pressure : " << press.val_grid_point(0,0,0,0)
	<< " rho_nuc c^2" << endl ;

    ost << endl ;

	Scalar psi4_local = pow(Psi_auto + Psi_comp + 1., 4.) ;
	psi4_local.std_spectral_base() ;

    ost << "Central lapse N            :  " << nn.val_grid_point(0,0,0,0) <<  endl ;
    ost << "Central value of Psi^4     :  " << psi4_local.val_grid_point(0,0,0,0) <<  endl ;

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

	double dent_eq   = ent.dsdr().val_point(ray_eq(),M_PI/2.,0.) ;
  	double dent_pole = ent.dsdr().val_point(ray_pole(),0.,0.) ;
  	double mass_shedd_chi = fabs( dent_eq / dent_pole ) ;

    ost << "Mass-shedding estimator = " << mass_shedd_chi << endl ;

    ost << endl << "Baryon mass         : " << mass_b() / msol << " M_sol" << endl ;
    ost << "Gravitational mass  : " << mass_g() / msol << " M_sol" << endl ;



    ost << endl ;
    ost << "Star in a binary system" << endl ;
    ost << "-----------------------" << endl ;

    if (irrotational) {
	ost << "irrotational configuration" << endl ;
    }
    else {
	ost << "corotating configuration" << endl ;
    }

    ost << "Absolute abscidia of the stellar center: " <<
	mp.get_ori_x() / km << " km" << endl ;

    ost << "Absolute abscidia of the barycenter of the baryon density : " <<
	xa_barycenter() / km << " km" << endl ;

    double r_0 = 0.5 * ( ray_eq() + ray_eq_pi() ) ;
    double d_ns = fabs( mp.get_ori_x() ) + ray_eq_pi() - r_0 ;
    double d_tilde = 2 * d_ns / r_0 ;

    ost << "d_tilde : " << d_tilde << endl ;

    ost << "Central value of gam_euler : "
        << gam_euler.val_grid_point(0, 0, 0, 0)  << endl ;

    ost << "Central u_euler (U^r, U^t, U^p) [c] : "
	<< u_euler(1).val_grid_point(0, 0, 0, 0) << "  "
	<< u_euler(2).val_grid_point(0, 0, 0, 0) << "  "
	<< u_euler(3).val_grid_point(0, 0, 0, 0) << endl ;

    if (irrotational) {
    ost << "Central d_psi (r, t, p) [c] :         "
	    << d_psi(1).val_grid_point(0, 0, 0, 0) << "  "
	    << d_psi(2).val_grid_point(0, 0, 0, 0) << "  "
	    << d_psi(3).val_grid_point(0, 0, 0, 0) << endl ;

	ost << "Central vel. / co-orb. (W^r, W^t, W^p) [c] : "
	    << wit_w(1).val_grid_point(0, 0, 0, 0) << "  "
	    << wit_w(2).val_grid_point(0, 0, 0, 0) << "  "
	    << wit_w(3).val_grid_point(0, 0, 0, 0) << endl ;

	ost << "Max vel. / co-orb. (W^r, W^t, W^p) [c] : "
	    << max(max(wit_w(1))) << "  "
	    << max(max(wit_w(2))) << "  "
	    << max(max(wit_w(3))) << endl ;

	ost << "Min vel. / co-orb. (W^r, W^t, W^p) [c] : "
	    << min(min(wit_w(1))) << "  "
	    << min(min(wit_w(2))) << "  "
	    << min(min(wit_w(3))) << endl ;

	double r_surf = mp.val_r(0,1.,M_PI/4,M_PI/4) ;

	ost << "Velocity at (r_surf,pi/4,pi/4) / co-orb. [c] : "
	    << wit_w(1).val_point(r_surf,M_PI/4,M_PI/4) << "  "
	    << wit_w(2).val_point(r_surf,M_PI/4,M_PI/4) << "  "
	    << wit_w(3).val_point(r_surf,M_PI/4,M_PI/4) << endl ;

	ost << "Central value of loggam : "
	    << loggam.val_grid_point(0, 0, 0, 0)  << endl ;
    }

    ost << "Central value of Psi auto, comp :         "
	<< Psi_auto.val_grid_point(0, 0, 0, 0) << "  "
	<< Psi_comp.val_grid_point(0, 0, 0, 0) << endl ;

    ost << "Central value of beta (N^r, N^t, N^p) [c] : "
	<< beta(1).val_grid_point(0, 0, 0, 0) << "  "
	<< beta(2).val_grid_point(0, 0, 0, 0) << "  "
	<< beta(3).val_grid_point(0, 0, 0, 0) << endl ;

    ost << "  ... beta_auto part of it [c] :            "
	<< beta_auto(1).val_grid_point(0, 0, 0, 0) << "  "
	<< beta_auto(2).val_grid_point(0, 0, 0, 0) << "  "
	<< beta_auto(3).val_grid_point(0, 0, 0, 0) << endl ;

    ost << endl << "Central value of (B^r, B^t, B^p)/N [c] : "
	<< bsn(1).val_grid_point(0, 0, 0, 0) << "  "
	<< bsn(2).val_grid_point(0, 0, 0, 0) << "  "
	<< bsn(3).val_grid_point(0, 0, 0, 0) << endl ;


    ost << endl << "Central \\hat{A}^{ij} [c/km] : " << endl ;
    ost << "  \\hat{A}^{xx} auto, comp : "
	<< haij_auto(1, 1).val_grid_point(0, 0, 0, 0) * km  << "  "
	<< haij_comp(1, 1).val_grid_point(0, 0, 0, 0) * km << endl ;
    ost << "  A^{xy} auto, comp : "
	<< haij_auto(1, 2).val_grid_point(0, 0, 0, 0) * km  << "  "
	<< haij_comp(1, 2).val_grid_point(0, 0, 0, 0) * km << endl ;
    ost << "  A^{xz} auto, comp : "
	<< haij_auto(1, 3).val_grid_point(0, 0, 0, 0) * km  << "  "
	<< haij_comp(1, 3).val_grid_point(0, 0, 0, 0) * km << endl ;
    ost << "  A^{yy} auto, comp : "
	<< haij_auto(2, 2).val_grid_point(0, 0, 0, 0) * km  << "  "
	<< haij_comp(2, 2).val_grid_point(0, 0, 0, 0) * km << endl ;
    ost << "  A^{yz} auto, comp : "
	<< haij_auto(2, 3).val_grid_point(0, 0, 0, 0) * km  << "  "
	<< haij_comp(2, 3).val_grid_point(0, 0, 0, 0) * km << endl ;
    ost << "  A^{zz} auto, comp : "
	<< haij_auto(3, 3).val_grid_point(0, 0, 0, 0) * km  << "  "
	<< haij_comp(3, 3).val_grid_point(0, 0, 0, 0) * km << endl ;

    ost << endl << "Central \\hat{A}_{ij}\\hat{A}^{ij} [c^2/km^2] : "
    	<< endl ;
    ost << "   \\hat{A}_{ij}\\hat{A}^{ij}  auto, comp : "
	<< hacar_auto.val_grid_point(0, 0, 0, 0) * km*km  << "  "
	<< hacar_comp.val_grid_point(0, 0, 0, 0) * km*km << endl ;

    return ost ;
}

			    //-------------------------//
			    //	Computational routines //
			    //-------------------------//

void Star_bin_xcts::fait_d_psi() {

    if (!irrotational) {
		d_psi.set_etat_nondef() ;
	return ;
    }

    // Specific relativistic enthalpy		    ---> hhh
    //----------------------------------

    Scalar hhh = exp(ent) ;  // = 1 at the Newtonian limit

    //  Computation of W^i = - h Gamma_n B^i/N
    //----------------------------------------------

    Vector www = hhh * gam_euler * bsn * psi4 ;

    // Constant value of W^i at the center of the star
    //-------------------------------------------------

    Vector v_orb(mp, COV, mp.get_bvect_cart()) ;

    for (int i=1; i<=3; i++)
		v_orb.set(i) = www(i).val_grid_point(0, 0, 0, 0) ;

    // Gradient of psi
    //----------------

    Vector d_psi0 = psi0.derive_cov(flat) ;

    d_psi0.change_triad( mp.get_bvect_cart() ) ;
    d_psi0.std_spectral_base() ;

    d_psi = d_psi0 + v_orb ;

    for (int i=1; i<=3; i++) {
	if (d_psi(i).get_etat() == ETATZERO)
	    d_psi.set(i).annule_hard() ;
    }

    // C^1 continuation of d_psi outside the star
    // (to ensure a smooth enthalpy field accross the stellar surface)
    // ----------------------------------------------------------------

    int nzm1 = mp.get_mg()->get_nzone() - 1 ;
    d_psi.annule(nzet, nzm1) ;

    for (int i=1; i<=3; i++) {
		Cmp d_psi_i (d_psi(i)) ;
		d_psi_i.va.set_base( d_psi0(i).get_spectral_va().base ) ;
		d_psi_i = raccord_c1(d_psi_i, nzet) ;
		d_psi.set(i) = d_psi_i ;
    }

    d_psi.std_spectral_base() ;

}

void Star_bin_xcts::relaxation(const Star_bin_xcts& star_jm1,
							double relax_ent,
			    			double relax_met,
			    		  	int mer,
			    			int fmer_met) {

    double relax_ent_jm1 = 1. - relax_ent ;
    double relax_met_jm1 = 1. - relax_met ;

    ent = relax_ent * ent + relax_ent_jm1 * star_jm1.ent ;

    if ( (mer != 0) && (mer % fmer_met == 0)) {

	Psi_auto = relax_met * Psi_auto
				+ relax_met_jm1 * star_jm1.Psi_auto ;

	chi_auto = relax_met * chi_auto
				+ relax_met_jm1 * star_jm1.chi_auto ;

	beta_auto = relax_met * beta_auto
				+ relax_met_jm1 * star_jm1.beta_auto ;

    }

    del_deriv() ;

    equation_of_state() ;

}
