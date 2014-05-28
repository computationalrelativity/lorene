/*
 * Methods for axisymmetric rotating neutron stars with magnetisation in the stress-energy tensor. 
 *
 * See the file et_rot_mag.h for documentation
 *
 */

/*
 *   Copyright (c) 2013-2014 Jerome Novak, Debarti Chatterjee, Micaela Oertel
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

char et_magnetisation_C[] = "$Header $" ;

/*
 * $Id$
 * $Log$
 * Revision 1.8  2014/05/28 07:46:06  j_novak
 * Minor modifications.
 *
 * Revision 1.7  2014/05/27 12:32:28  j_novak
 * Added possibility to converge to a given magnetic moment.
 *
 * Revision 1.6  2014/05/14 15:19:05  j_novak
 * The magnetisation field is now filtered.
 *
 * Revision 1.5  2014/04/29 13:46:07  j_novak
 * Addition of switches 'use_B_in_eos' and 'include_magnetisation' to control the model.
 *
 * Revision 1.4  2014/04/28 12:48:13  j_novak
 * Minor modifications.
 *
 * Revision 1.3  2013/12/19 17:05:40  j_novak
 * Corrected a dzpuis problem.
 *
 * Revision 1.2  2013/12/13 16:36:51  j_novak
 * Addition and computation of magnetisation terms in the Einstein equations.
 *
 * Revision 1.1  2013/11/25 13:52:11  j_novak
 * New class Et_magnetisation to include magnetization terms in the stress energy tensor.
 *
 *
 * $Header$
 *
 */

// Headers C
#include <cmath>

// Headers Lorene
#include "et_rot_mag.h"
#include "utilitaires.h"
#include "unites.h"
#include "param.h"
#include "eos.h"

			    //--------------//
			    // Constructors //
			    //--------------//
// Standard constructor
// --------------------


Et_magnetisation::Et_magnetisation(Map& mp_i, int nzet_i, bool relat, 
				   const Eos& eos_i, bool magnetisation, 
				   bool B_eos)
  : Et_rot_mag(mp_i, nzet_i, relat, eos_i, 1),
    use_B_in_eos(B_eos), include_magnetisation(magnetisation),
    xmag(mp_i),
    E_I(mp_i),
    J_I(mp_i, COV, mp_i.get_bvect_spher()),
    Sij_I(mp_i, COV, mp_i.get_bvect_spher())
{
  const Eos_mag* tmp_p_eos = dynamic_cast<const Eos_mag*>(&eos_i) ;
  if (tmp_p_eos == 0x0) {
    cerr << "Et_magnetisation::Et_magnetisation : " << endl ;
    cerr << "Only magnetised EoS is admitted for this class!" << endl ;
    cerr << "Aborting ... " << endl ;
    abort() ;
  }
  if ( include_magnetisation && (!use_B_in_eos) ) {
    cerr << "Et_magnetisation::Et_magnetisation : " << endl ;
    cerr << "Magnetisation terms can be included only if " 
	 << "the magnetic field is used in the EoS!" << endl;
    cerr << "Aborting ... " << endl ;
    abort() ;
  }
 xmag = 0 ;
  E_I = 0 ;
  J_I.set_etat_zero() ;
  Sij_I.set_etat_zero() ;

}



Et_magnetisation::Et_magnetisation(Map& mp_i, const Eos& eos_i, FILE* fich)
  : Et_rot_mag(mp_i, eos_i, fich, 0), 
    use_B_in_eos(true), include_magnetisation(true),
    xmag(mp_i),
    E_I(mp_i),
    J_I(mp_i, COV, mp_i.get_bvect_spher()),
    Sij_I(mp_i, COV, mp_i.get_bvect_spher())
{

  // Read of the saved fields:
  // ------------------------

  fread(&use_B_in_eos, sizeof(bool), 1, fich) ;

  fread(&include_magnetisation, sizeof(bool), 1, fich) ;
    
  Scalar xmag_file(mp, *mp.get_mg(), fich) ;
  xmag = xmag_file ;
  
  Scalar E_I_file(mp, *mp.get_mg(), fich) ;
  E_I = E_I_file ;
  
  Vector J_I_file(mp, mp.get_bvect_spher(), fich) ;
  J_I = J_I_file ;
  
  Sym_tensor Sij_I_file(mp, mp.get_bvect_spher(), fich) ;
  Sij_I = Sij_I_file ;
  
}


// Copy constructor
// ----------------

Et_magnetisation::Et_magnetisation(const Et_magnetisation& et)
  : Et_rot_mag(et),
    use_B_in_eos(et.use_B_in_eos), include_magnetisation(et.include_magnetisation),
    xmag(et.xmag),
    E_I(et.E_I),
    J_I(et.J_I),
    Sij_I(et.Sij_I)
{  }


			    //------------//
			    // Destructor //
			    //------------//

Et_magnetisation::~Et_magnetisation(){ }


// Assignment to another Et_magnetisation
// --------------------------------

void Et_magnetisation::operator=(const Et_magnetisation& et) {

  // Assignement of quantities common to all the derived classes of Et_rot_mag
  Et_rot_mag::operator=(et) ;
  use_B_in_eos = et.use_B_in_eos;
  include_magnetisation = et.include_magnetisation ;
  xmag   = et.xmag  ;
  E_I    = et.E_I   ;
  J_I    = et.J_I   ;
  Sij_I  = et.Sij_I ;

}


void Et_magnetisation::equation_of_state() {

  Cmp ent_eos = ent() ;
  const Eos_mag* mageos = dynamic_cast<const Eos_mag*>(&eos) ;
  assert (mageos != 0x0) ;

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
	  xi.set(l,k,j,i) = mg->get_grille3d(l)->x[i] ;
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
  
  if (nzet > 1) {
    
    if (nzet == 3) {
      fact_ent.set(1) = 1 - 0.5 * epsilon * (xi(1) - 0.5) * (xi(1) - 0.5) ;
      fact_ent.set(2) = 1 - 0.25 * epsilon * (xi(2) - 1) * (xi(2) - 1) ;	
    }
    
    if (nzet > 3) {
      cerr << "Et_magnetisation::equation_of_state: "
	   << "not ready yet for nzet > 3 !" << endl ;    	
    }
    
    ent_eos = fact_ent * ent_eos ;
    ent_eos.std_base_scal() ;
  }
  
  
  
  // Call to a magnetized EOS
  // with the norm of the magnetic field passed as an argument to the EoS.
  
  double magb0 = 0 ;
  Cmp norm_b(mp) ;
  norm_b.set_etat_zero() ;
  if (use_B_in_eos) {
    Tenseur Bmag = Magn() ;
    norm_b = sqrt(a_car() * ( Bmag(0)*Bmag(0) + Bmag(1)*Bmag(1) )) 
      / gam_euler() ; // Use of b^\mu in the fluid frame
    norm_b.std_base_scal() ;
  }
  Param par ;
  par.add_double_mod(magb0) ;
  
  nbar.set_etat_qcq() ; nbar.set().set_etat_qcq() ;
  nbar.set().va.set_etat_c_qcq() ; nbar.set().va.c->set_etat_qcq() ;
  ener.set_etat_qcq() ; ener.set().set_etat_qcq() ;
  ener.set().va.set_etat_c_qcq() ; ener.set().va.c->set_etat_qcq() ;
  press.set_etat_qcq() ; press.set().set_etat_qcq() ;
  press.set().va.set_etat_c_qcq() ; press.set().va.c->set_etat_qcq() ;
  xmag.allocate_all() ;
  
  for (int l=0; l< nz; l++) {
    Tbl* tent = ent_eos.va.c->t[l] ;
    if ( (tent->get_etat() == ETATZERO) || (l >= nzet) ) {
      nbar.set().set(l).set_etat_zero() ;
      ener.set().set(l).set_etat_zero() ;
      press.set().set(l).set_etat_zero() ;
      xmag.annule_domain(l) ;
    }
    else {
      nbar.set().set(l).annule_hard() ;
      ener.set().set(l).annule_hard() ;
      press.set().set(l).annule_hard() ;
      xmag.set_domain(l).annule_hard() ;
      for (int k=0; k < mg->get_np(l); k++) {
	for (int j=0; j < mg->get_nt(l); j++) {
	  for (int i=0; i < mg->get_nr(l); i++) {
	    magb0 = norm_b(l, k, j, i) ;
	    double ent0 = ent_eos(l, k, j, i) ;
	    nbar.set().set(l, k, j, i) = mageos->nbar_ent_p(ent0, &par) ;
	    ener.set().set(l, k, j, i) = mageos->ener_ent_p(ent0, &par) ;
	    press.set().set(l, k, j, i) = mageos->press_ent_p(ent0, &par) ;
	    xmag.set_grid_point(l, k, j, i) = mageos->mag_ent_p(ent0, &par) ;
	  }
	}
      }
    }
  }
  
  if (!include_magnetisation)
    xmag.set_etat_zero() ;

  // Set the bases for spectral expansion
  nbar.set_std_base() ; 
  ener.set_std_base() ; 
  press.set_std_base() ; 
  xmag.std_spectral_base() ;

  Scalar tmp_scal(xmag) ;
  tmp_scal.exponential_filter_r(0, 0, 1) ;
  xmag = tmp_scal ;

  // The derived quantities are obsolete
  del_deriv() ; 
    
}



			    //--------------//
			    //	  Outputs   //
			    //--------------//

// Save in a file
// --------------
void Et_magnetisation::sauve(FILE* fich) const {
    
    Et_rot_mag::sauve(fich) ; 

    fwrite(&use_B_in_eos, sizeof(bool), 1, fich) ;		
    
    fwrite(&include_magnetisation, sizeof(bool), 1, fich) ;		
    
    xmag.sauve(fich) ;
    E_I.sauve(fich) ;
    J_I.sauve(fich) ;
    Sij_I.sauve(fich) ;
    
}


// Printing
// --------


ostream& Et_magnetisation::operator>>(ostream& ost) const {

  using namespace Unites_mag ;

  Et_rot_mag::operator>>(ost) ;
  //  int theta_eq = mp.get_mg()->get_nt(nzet-1)-1 ;
  ost << endl ;
  ost << "Rotating magnetized neutron star"
      << endl ;
  if (use_B_in_eos)
    ost << "Using magnetic field in the EoS" << endl ;
  if (include_magnetisation) {
    ost << "Including magnetisation terms in the equations" << endl ;
    ost << "Maximal value of the magnetization scalar x : " 
	<< max(maxabs(xmag)) << endl ;
  }
  return ost ;
}

//=================================================================================
//
//                Computation of equilibrium configuration
//
//=================================================================================


void Et_magnetisation::equilibrium_mag(double ent_c, double omega0, 
     double fact_omega, int nzadapt, const Tbl& ent_limit, 
     const Itbl& icontrol, const Tbl& control, double mbar_wanted,
				       double magmom_wanted,
     double aexp_mass, Tbl& diff, double Q0, double a_j0, 
     Cmp (*f_j)(const Cmp&, const double), 
     Cmp (*M_j)(const Cmp& x, const double)) {
			     
  // Fundamental constants and units
  // -------------------------------
  using namespace Unites_mag ;
    
  // For the display 
  // ---------------
  char display_bold[]="x[1m" ; display_bold[0] = 27 ;
  char display_normal[] = "x[0m" ; display_normal[0] = 27 ;

  // Grid parameters
  // ---------------
  
  const Mg3d* mg = mp.get_mg() ; 
  int nz = mg->get_nzone() ;	    // total number of domains
  int nzm1 = nz - 1 ; 
  
  // The following is required to initialize mp_prev as a Map_et:
  Map_et& mp_et = dynamic_cast<Map_et&>(mp) ; 
  
  // Index of the point at phi=0, theta=pi/2 at the surface of the star:
  assert(mg->get_type_t() == SYM) ; 
  int l_b = nzet - 1 ; 
  int i_b = mg->get_nr(l_b) - 1 ; 
  int j_b = mg->get_nt(l_b) - 1 ; 
  int k_b = 0 ; 
    
  // Value of the enthalpy defining the surface of the star
  double ent_b = ent_limit(nzet-1) ;
    
  // Parameters to control the iteration
  // -----------------------------------
    
  int mer_max = icontrol(0) ; 
  int mer_rot = icontrol(1) ;
  int mer_change_omega = icontrol(2) ; 
  int mer_fix_omega = icontrol(3) ; 
  int mer_mass = icontrol(4) ; 
  int mermax_poisson = icontrol(5) ; 
  int delta_mer_kep = icontrol(6) ; 
  int mer_mag = icontrol(7) ;
  int mer_change_mag = icontrol(8) ;
  int mer_fix_mag = icontrol(9) ;
  int mer_magmom = icontrol(10) ;

  // Protections:
  if (mer_change_omega < mer_rot) {
    cerr << "Et_magnetisation::equilibrium_mag: mer_change_omega < mer_rot !" 
	 << endl ;
    cerr << " mer_change_omega = " << mer_change_omega << endl ; 
    cerr << " mer_rot = " << mer_rot << endl ; 
    abort() ; 
  }
  if (mer_fix_omega < mer_change_omega) {
    cerr << "Et_magnetisation::equilibrium_mag: mer_fix_omega < mer_change_omega !" 
	 << endl ;
    cerr << " mer_fix_omega = " << mer_fix_omega << endl ; 
    cerr << " mer_change_omega = " << mer_change_omega << endl ; 
    abort() ; 
  }

  // In order to converge to a given baryon mass, shall the central
  // enthalpy be varied or Omega ?
  bool change_ent = true ; 
  if (mer_mass < 0) {
    change_ent = false ; 
    mer_mass = abs(mer_mass) ;
  }

  double precis = control(0) ; 
  double omega_ini = control(1) ; 
  double relax = control(2) ;
  double relax_prev = double(1) - relax ;  
  double relax_poisson = control(3) ; 
  double thres_adapt = control(4) ; 
  double precis_adapt = control(5) ; 
  double Q_ini = control(6) ;
  double a_j_ini = control (7) ;

  // Error indicators
  // ----------------
    
  diff.set_etat_qcq() ; 
  double& diff_ent = diff.set(0) ; 

  // Parameters for the function Map_et::adapt
  // -----------------------------------------
    
  Param par_adapt ; 
  int nitermax = 100 ;  
  int niter ; 
  int adapt_flag = 1 ;    //  1 = performs the full computation, 
                          //  0 = performs only the rescaling by 
                          //      the factor alpha_r
  int nz_search = nzet + 1 ;  // Number of domains for searching the enthalpy
			      //  isosurfaces
  double alpha_r ; 
  double reg_map = 1. ; // 1 = regular mapping, 0 = contracting mapping

  par_adapt.add_int(nitermax, 0) ; // maximum number of iterations to
				   // locate zeros by the secant method
  par_adapt.add_int(nzadapt, 1) ; // number of domains where the adjustment 
				  // to the isosurfaces of ent is to be 
				  // performed
  par_adapt.add_int(nz_search, 2) ; // number of domains to search for
				    // the enthalpy isosurface
  par_adapt.add_int(adapt_flag, 3) ; //  1 = performs the full computation, 
				     //  0 = performs only the rescaling by 
				     //      the factor alpha_r
  par_adapt.add_int(j_b, 4) ; //  theta index of the collocation point 
			      //  (theta_*, phi_*)
  par_adapt.add_int(k_b, 5) ; //  theta index of the collocation point 
			      //  (theta_*, phi_*)

  par_adapt.add_int_mod(niter, 0) ;  //  number of iterations actually used in 
				     //  the secant method
    
  par_adapt.add_double(precis_adapt, 0) ; // required absolute precision in 
					  // the determination of zeros by 
					  // the secant method
  par_adapt.add_double(reg_map, 1)	;  // 1. = regular mapping, 
                                           // 0 = contracting mapping
    
  par_adapt.add_double(alpha_r, 2) ;	// factor by which all the radial 
					// distances will be multiplied 
    	   
  par_adapt.add_tbl(ent_limit, 0) ;	// array of values of the field ent 
				        // to define the isosurfaces.


  // Parameters for the function Map_et::poisson for nuf
  // ----------------------------------------------------
  
  double precis_poisson = 1.e-16 ;     
  
  Param par_poisson_nuf ; 
  par_poisson_nuf.add_int(mermax_poisson,  0) ;  // maximum number of iterations
  par_poisson_nuf.add_double(relax_poisson,  0) ; // relaxation parameter
  par_poisson_nuf.add_double(precis_poisson, 1) ; // required precision
  par_poisson_nuf.add_int_mod(niter, 0) ;  //  number of iterations actually used 
  par_poisson_nuf.add_cmp_mod( ssjm1_nuf ) ; 
  
  Param par_poisson_nuq ; 
  par_poisson_nuq.add_int(mermax_poisson,  0) ;  // maximum number of iterations
  par_poisson_nuq.add_double(relax_poisson,  0) ; // relaxation parameter
  par_poisson_nuq.add_double(precis_poisson, 1) ; // required precision
  par_poisson_nuq.add_int_mod(niter, 0) ;  //  number of iterations actually used 
  par_poisson_nuq.add_cmp_mod( ssjm1_nuq ) ; 
  
  Param par_poisson_tggg ; 
  par_poisson_tggg.add_int(mermax_poisson,  0) ;  // maximum number of iterations
  par_poisson_tggg.add_double(relax_poisson,  0) ; // relaxation parameter
  par_poisson_tggg.add_double(precis_poisson, 1) ; // required precision
  par_poisson_tggg.add_int_mod(niter, 0) ;  // number of iterations actually used 
  par_poisson_tggg.add_cmp_mod( ssjm1_tggg ) ; 
  double lambda_tggg ;
  par_poisson_tggg.add_double_mod( lambda_tggg ) ; 
  
  Param par_poisson_dzeta ; 
  double lbda_grv2 ;
  par_poisson_dzeta.add_double_mod( lbda_grv2 ) ; 

  // Parameters for the function Tenseur::poisson_vect
  // -------------------------------------------------
  
  Param par_poisson_vect ; 
  
  par_poisson_vect.add_int(mermax_poisson,  0) ;  // maximum number of iterations
  par_poisson_vect.add_double(relax_poisson,  0) ; // relaxation parameter
  par_poisson_vect.add_double(precis_poisson, 1) ; // required precision
  par_poisson_vect.add_cmp_mod( ssjm1_khi ) ; 
  par_poisson_vect.add_tenseur_mod( ssjm1_wshift ) ; 
  par_poisson_vect.add_int_mod(niter, 0) ;   
					   
  // Parameters for the Maxwell equations
  // -------------------------------------

  Param par_poisson_At ; // For scalar At Poisson equation
  Cmp ssjm1_At(mp) ;
  ssjm1_At.set_etat_zero() ;
  par_poisson_At.add_int(mermax_poisson,  0) ;  // maximum number of iterations
  par_poisson_At.add_double(relax_poisson,  0) ; // relaxation parameter
  par_poisson_At.add_double(precis_poisson, 1) ; // required precision
  par_poisson_At.add_int_mod(niter, 0) ;  //  number of iterations actually used 
  par_poisson_At.add_cmp_mod( ssjm1_At ) ; 
  
  Param par_poisson_Avect ;  // For vector Aphi Poisson equation
  
  Cmp ssjm1_khi_mag(ssjm1_khi) ;
  Tenseur ssjm1_w_mag(ssjm1_wshift) ;
  
  par_poisson_Avect.add_int(mermax_poisson,  0) ;  // maximum number of iterations
  par_poisson_Avect.add_double(relax_poisson,  0) ; // relaxation parameter
  par_poisson_Avect.add_double(precis_poisson, 1) ; // required precision
  par_poisson_Avect.add_cmp_mod( ssjm1_khi_mag ) ; 
  par_poisson_Avect.add_tenseur_mod( ssjm1_w_mag ) ; 
  par_poisson_Avect.add_int_mod(niter, 0) ;   
  				   
  // Initializations
  // ---------------
  
  // Initial angular velocity / magnetic quantities
  omega = 0 ; 
  Q = 0 ;
  a_j = 0 ;

  double accrois_omega = (omega0 - omega_ini) 
    / double(mer_fix_omega - mer_change_omega) ; 
  double accrois_Q = (Q0 - Q_ini) 
    / double(mer_fix_mag - mer_change_mag);
  double accrois_a_j = (a_j0 - a_j_ini) 
    / double(mer_fix_mag - mer_change_mag); 

  update_metric() ;	// update of the metric coefficients

  equation_of_state() ;	// update of the density, pressure, etc...
    
  hydro_euler() ;	// update of the hydro quantities relative to the 
			//  Eulerian observer

  MHD_comput() ; // update of EM contributions to stress-energy tensor


  // Quantities at the previous step : 	
  Map_et mp_prev = mp_et ; 
  Tenseur ent_prev = ent ;	    
  Tenseur logn_prev = logn ;	    
  Tenseur dzeta_prev = dzeta ;	    
  
  // Creation of uninitialized tensors:
  Tenseur source_nuf(mp) ;    // source term in the equation for nuf
  Tenseur source_nuq(mp) ;    // source term in the equation for nuq
  Tenseur source_dzf(mp) ;	// matter source term in the eq. for dzeta
  Tenseur source_dzq(mp) ;	// quadratic source term in the eq. for dzeta
  Tenseur source_tggg(mp) ;	// source term in the eq. for tggg
  Tenseur source_shift(mp, 1, CON, mp.get_bvect_cart()) ; // source term for shift
  Tenseur mlngamma(mp) ;	// centrifugal potential
    
  // Preparations for the Poisson equations:
  // --------------------------------------
  if (nuf.get_etat() == ETATZERO) {
    nuf.set_etat_qcq() ; 
    nuf.set() = 0 ; 
  }
  
  if (relativistic) {
    if (nuq.get_etat() == ETATZERO) {
      nuq.set_etat_qcq() ; 
      nuq.set() = 0 ; 
    }
    
    if (tggg.get_etat() == ETATZERO) {
      tggg.set_etat_qcq() ; 
      tggg.set() = 0 ; 
    }
    
    if (dzeta.get_etat() == ETATZERO) {
      dzeta.set_etat_qcq() ; 
      dzeta.set() = 0 ; 
    }
  }
		    
  ofstream fichconv("convergence.d") ;    // Output file for diff_ent
  fichconv << "#     diff_ent     GRV2    " << endl ; 
  
  ofstream fichfreq("frequency.d") ;    // Output file for  omega
  fichfreq << "#       f [Hz]" << endl ; 
  
  ofstream fichevol("evolution.d") ;    // Output file for various quantities
  fichevol << 
    "#       |dH/dr_eq/dH/dr_pole|      r_pole/r_eq	ent_c" 
	   << endl ; 
  
  diff_ent = 1 ; 
  double err_grv2 = 1 ; 
    
  //=========================================================================
  // 			Start of iteration
  //=========================================================================

  for(int mer=0 ; (diff_ent > precis) && (mer<mer_max) ; mer++ ) {

    cout << "-----------------------------------------------" << endl ;
    cout << "step: " << mer << endl ;
    cout << "diff_ent = " << display_bold << diff_ent << display_normal
	 << endl ;    
    cout << "err_grv2 = " << err_grv2 << endl ;    
    fichconv << mer ;
    fichfreq << mer ;
    fichevol << mer ;
    
    if (mer >= mer_rot) {
      
      if (mer < mer_change_omega) {
	omega = omega_ini ; 
      }
      else {
	if (mer <= mer_fix_omega) {
	  omega = omega_ini + accrois_omega * 
	    (mer - mer_change_omega) ;
	}
      }
      
    }

    if (mer >= mer_mag) {
      if (mer < mer_change_mag) {
	Q   = Q_ini ;
	a_j = a_j_ini ;
      }
      else {
	if (mer <= mer_fix_mag) {
	  Q = Q_ini + accrois_Q * (mer - mer_change_mag) ;
	  a_j = a_j_ini + accrois_a_j * (mer - mer_change_mag) ;
	}
      }
    }

    //-----------------------------------------------
    // Computation of electromagnetic potentials :
    // -------------------------------------------
    magnet_comput(adapt_flag, f_j, par_poisson_At, par_poisson_Avect) ;

    MHD_comput() ; // computes EM contributions to T_{mu,nu}

    // S_{rr} + S_{\theta\theta}
    Tenseur SrrplusStt( Cmp(Sij_I(1, 1) + Sij_I(2, 2)) ) ; 
    SrrplusStt = SrrplusStt / a_car ; // S^r_r + S^\theta_\theta

    Tenseur Spp (Cmp(Sij_I(3, 3))) ; //S_{\phi\phi}
    Spp = Spp / b_car ; // S^\phi_\phi

    //-----------------------------------------------
    //  Sources of the Poisson equations
    //-----------------------------------------------
    
    // Source for nu
    // -------------
    Tenseur beta = log(bbb) ; 
    beta.set_std_base() ; 
    
    if (relativistic) {
      source_nuf =  
	qpig * a_car *( ener_euler + s_euler + SrrplusStt + Spp ) ; 
      
      source_nuq = ak_car 
	- flat_scalar_prod(logn.gradient_spher(), logn.gradient_spher() 
			   + beta.gradient_spher()) 
	+ qpig * a_car * 2*E_em ;
    }
    else {
      source_nuf = qpig * nbar ; 
      
      source_nuq = 0 ; 
    }
    source_nuf.set_std_base() ; 	
    source_nuq.set_std_base() ; 	

    // Source for dzeta
    // ----------------
    source_dzf = 2 * qpig * a_car 
      * (press + (ener_euler+press) * uuu*uuu + Spp ) ;
    source_dzf.set_std_base() ; 
 
    source_dzq = 2 * qpig * a_car * E_em + 1.5 * ak_car - 
      flat_scalar_prod(logn.gradient_spher(), logn.gradient_spher() ) ;  
    source_dzq.set_std_base() ; 	
	

    // Source for tggg
    // ---------------
    
    source_tggg = 2 * qpig * nnn * a_car * bbb 
      * ( 2*press + SrrplusStt ) ;
    source_tggg.set_std_base() ; 
	
    (source_tggg.set()).mult_rsint() ; 
    
    
    // Source for shift
    // ----------------
	    
    // Matter term: 
    
    Cmp tjpem(J_I(3)) ;
    tjpem += Jp_em() ;
    tjpem.div_rsint() ;
    
    source_shift = (-4*qpig) * nnn * a_car * (ener_euler + press) * u_euler ;

    // Quadratic terms:
    Tenseur vtmp =  6 * beta.gradient_spher() - 2 * logn.gradient_spher() ;
    Tenseur mtmp(mp, 1, COV, mp.get_bvect_spher()) ;
    if (tjpem.get_etat() == ETATZERO) mtmp.set_etat_zero() ;
    else {
      mtmp.set_etat_qcq() ;
      mtmp.set(0) = 0 ;
      mtmp.set(1) = 0 ;
      mtmp.set(2) = (-4*qpig)*tjpem*nnn()*a_car()/b_car() ;
    }
    mtmp.change_triad(mp.get_bvect_cart()) ; 
    
    vtmp.change_triad(mp.get_bvect_cart()) ; 
    
    Tenseur squad  = nnn * flat_scalar_prod(tkij, vtmp) ;     
    
    // The addition of matter terms and quadratic terms is performed
    //  component by component because u_euler is contravariant,
    //  while squad is covariant. 
    
    if (squad.get_etat() == ETATQCQ) {
      for (int i=0; i<3; i++) {
	source_shift.set(i) += squad(i) ; 
      }
    }
    if (mtmp.get_etat() == ETATQCQ) {
      if (source_shift.get_etat() == ETATZERO) {
	source_shift.set_etat_qcq() ;
	for (int i=0; i<3; i++) {
	  source_shift.set(i) = mtmp(i) ;
	  source_shift.set(i).va.coef_i() ;
	}
      }
      else
	for (int i=0; i<3; i++) 
	  source_shift.set(i) += mtmp(i) ; 
    }
    
    source_shift.set_std_base() ; 	
    
    //----------------------------------------------
    // Resolution of the Poisson equation for nuf 
    //----------------------------------------------
    
    source_nuf().poisson(par_poisson_nuf, nuf.set()) ; 
		
    if (relativistic) {
      
      //----------------------------------------------
      // Resolution of the Poisson equation for nuq 
      //----------------------------------------------
      
      source_nuq().poisson(par_poisson_nuq, nuq.set()) ; 
      
      //---------------------------------------------------------
      // Resolution of the vector Poisson equation for the shift
      //---------------------------------------------------------
      
      
      if (source_shift.get_etat() != ETATZERO) {
	
	for (int i=0; i<3; i++) {
	  if(source_shift(i).dz_nonzero()) {
	    assert( source_shift(i).get_dzpuis() == 4 ) ; 
	  }
	  else{
	    (source_shift.set(i)).set_dzpuis(4) ; 
	  }
	}
	
      }
      //##
      // source_shift.dec2_dzpuis() ;    // dzpuis 4 -> 2
      
      double lambda_shift = double(1) / double(3) ; 
      
      if ( mg->get_np(0) == 1 ) {
	lambda_shift = 0 ; 
      }
      
      source_shift.poisson_vect(lambda_shift, par_poisson_vect, 
				shift, w_shift, khi_shift) ;      
	    
      // Computation of tnphi and nphi from the Cartesian components
      //  of the shift
      // -----------------------------------------------------------
      
      fait_nphi() ; 
	
    }
    
    //-----------------------------------------
    // Determination of the fluid velociy U
    //-----------------------------------------
    
    if (mer > mer_fix_omega + delta_mer_kep) {
      
      omega *= fact_omega ;  // Increase of the angular velocity if 
    }			   //  fact_omega != 1
    
    bool omega_trop_grand = false ; 
    bool kepler = true ; 
    
    while ( kepler ) {
      
      // Possible decrease of Omega to ensure a velocity < c 
      
      bool superlum = true ; 
      
      while ( superlum ) {
	
	// New fluid velocity U :
	
	Cmp tmp = omega - nphi() ; 
	tmp.annule(nzm1) ; 
	tmp.std_base_scal() ;
	
	tmp.mult_rsint() ;	    //  Multiplication by r sin(theta)
	
	uuu = bbb() / nnn() * tmp ; 
	
	if (uuu.get_etat() == ETATQCQ) {
	  // Same basis as (Omega -N^phi) r sin(theta) :
	  ((uuu.set()).va).set_base( (tmp.va).base ) ;   
	}
	
	// Is the new velocity larger than c in the equatorial plane ?
	
	superlum = false ; 
	
	for (int l=0; l<nzet; l++) {
	  for (int i=0; i<mg->get_nr(l); i++) {
	    
	    double u1 = uuu()(l, 0, j_b, i) ; 
	    if (u1 >= 1.) {	    // superluminal velocity
	      superlum = true ; 
	      cout << "U > c  for l, i : " << l << "  " << i 
		   << "   U = " << u1 << endl ;  
	    }
	  }
	}
	if ( superlum ) {
	  cout << "**** VELOCITY OF LIGHT REACHED ****" << endl ; 
	  omega /= fact_omega ;    // Decrease of Omega
	  cout << "New rotation frequency : " 
	       << omega/(2.*M_PI) * f_unit <<  " Hz" << endl ; 
	  omega_trop_grand = true ;  
	}
      }	// end of while ( superlum )
      
      
      // New computation of U (which this time is not superluminal)
      //  as well as of gam_euler, ener_euler, etc...
      // -----------------------------------
      
      hydro_euler() ; 
            
      //------------------------------------------------------
      //	First integral of motion 
      //------------------------------------------------------
      
      // Centrifugal potential : 
      if (relativistic) {
	mlngamma = - log( gam_euler ) ;
      }
      else {
	mlngamma = - 0.5 * uuu*uuu ; 
      }
      
      Tenseur mag(mp) ;
      if (is_conduct()) {
	mag = mu0*M_j(A_phi, a_j) ;}
      else{
	mag = mu0*M_j(omega*A_phi-A_t, a_j) ;}
      
      // Equatorial values of various potentials :
      double nuf_b  = nuf()(l_b, k_b, j_b, i_b) ; 
      double nuq_b  = nuq()(l_b, k_b, j_b, i_b) ; 
      double mlngamma_b  = mlngamma()(l_b, k_b, j_b, i_b) ; 
      double mag_b = mag()(l_b, k_b, j_b, i_b) ; 
      
      // Central values of various potentials :
      double nuf_c = nuf()(0,0,0,0) ; 
      double nuq_c = nuq()(0,0,0,0) ; 
      double mlngamma_c = 0 ;
      double mag_c = mag()(0,0,0,0) ;
      
      // Scale factor to ensure that the enthalpy is equal to ent_b at 
      //  the equator
      double alpha_r2 = ( ent_c - ent_b + mlngamma_c - mlngamma_b
			  + nuq_c - nuq_b + mag_c - mag_b) 
	/ ( nuf_b - nuf_c  ) ;
      alpha_r = sqrt(alpha_r2) ;
      cout << "alpha_r = " << alpha_r << endl ; 
      
      // Readjustment of nu :
      // -------------------
      
      logn = alpha_r2 * nuf + nuq ;
      double nu_c =  logn()(0,0,0,0) ;
      
      // First integral	--> enthalpy in all space
      //-----------------
      ent = (ent_c + nu_c + mlngamma_c + mag_c) - logn - mlngamma - mag ;
      
      // Test: is the enthalpy negative somewhere in the equatorial plane
      //  inside the star ? If yes, this means that the Keplerian velocity
      //  has been overstep.
      
      kepler = false ; 
      for (int l=0; l<nzet; l++) {
	int imax = mg->get_nr(l) - 1 ;
	if (l == l_b) imax-- ;	// The surface point is skipped
	for (int i=0; i<imax; i++) { 
	  if ( ent()(l, 0, j_b, i) < 0. ) {
	    kepler = true ;
	    cout << "ent < 0 for l, i : " << l << "  " << i 
		 << "   ent = " << ent()(l, 0, j_b, i) << endl ;  
	  } 
	}
      }
      
      if ( kepler ) {
	cout << "**** KEPLERIAN VELOCITY REACHED ****" << endl ; 
	omega /= fact_omega ;    // Omega is decreased
	cout << "New rotation frequency : " 
	     << omega/(2.*M_PI) * f_unit << " Hz" << endl ; 
	omega_trop_grand = true ;  
      }
      
    }   // End of while ( kepler )
    
    if ( omega_trop_grand ) {	// fact_omega is decreased for the
      //  next step 
      fact_omega = sqrt( fact_omega ) ; 
      cout << "**** New fact_omega : " << fact_omega << endl ; 
    }
    
    //----------------------------------------------------
    // Adaptation of the mapping to the new enthalpy field
    //----------------------------------------------------
    
    // Shall the adaptation be performed (cusp) ?
    // ------------------------------------------
    
    double dent_eq = ent().dsdr()(l_b, k_b, j_b, i_b) ; 
    double dent_pole = ent().dsdr()(l_b, k_b, 0, i_b) ;
    double rap_dent = fabs( dent_eq / dent_pole ) ; 
    cout << "| dH/dr_eq / dH/dr_pole | = " << rap_dent << endl ; 
    
    if ( rap_dent < thres_adapt ) {
      adapt_flag = 0 ;	// No adaptation of the mapping 
      cout << "******* FROZEN MAPPING  *********" << endl ; 
    }
    else{
      adapt_flag = 1 ;	// The adaptation of the mapping is to be
      //  performed
    }
    
    mp_prev = mp_et ; 
    
    mp.adapt(ent(), par_adapt) ; 
    
    //----------------------------------------------------
    // Computation of the enthalpy at the new grid points
    //----------------------------------------------------
    
    mp_prev.homothetie(alpha_r) ; 
    
    mp.reevaluate(&mp_prev, nzet+1, ent.set()) ; 
    
    //----------------------------------------------------
    // Equation of state  
    //----------------------------------------------------
    
    equation_of_state() ; 	// computes new values for nbar (n), ener (e) 
				// and press (p) from the new ent (H)
    
    //---------------------------------------------------------
    // Matter source terms in the gravitational field equations	
    //---------------------------------------------------------
    
    //## Computation of tnphi and nphi from the Cartesian components
    //  of the shift for the test in hydro_euler():
    
    fait_nphi() ; 
    
    hydro_euler() ;		// computes new values for ener_euler (E), 
				// s_euler (S) and u_euler (U^i)
    
    if (relativistic) {
      
      //-------------------------------------------------------
      //	2-D Poisson equation for tggg
      //-------------------------------------------------------
      
      mp.poisson2d(source_tggg(), mp.cmp_zero(), par_poisson_tggg, tggg.set()) ; 
	    
      //-------------------------------------------------------
      //	2-D Poisson equation for dzeta
      //-------------------------------------------------------
      
      mp.poisson2d(source_dzf(), source_dzq(), par_poisson_dzeta, dzeta.set()) ; 
	    
      err_grv2 = lbda_grv2 - 1; 
      cout << "GRV2: " << err_grv2 << endl ; 
      
    }
    else {
      err_grv2 = grv2() ; 
    }
        
    //---------------------------------------
    // Computation of the metric coefficients (except for N^phi)
    //---------------------------------------
    
    // Relaxations on nu and dzeta :  
    
    if (mer >= 10) {
      logn = relax * logn + relax_prev * logn_prev ;
      
      dzeta = relax * dzeta + relax_prev * dzeta_prev ; 
    }
    
    // Update of the metric coefficients N, A, B and computation of K_ij :
    
    update_metric() ; 
    
    //-----------------------
    //  Informations display
    //-----------------------
    
    //	partial_display(cout) ; 
    fichfreq << "  " << omega / (2*M_PI) * f_unit ; 
    fichevol << "  " << rap_dent ; 
    fichevol << "  " << ray_pole() / ray_eq() ; 
    fichevol << "  " << ent_c ; 
    
    //-----------------------------------------
    // Convergence towards a given baryon mass 
    //-----------------------------------------
    
    if (mer > mer_mass) {
      
      double xx ; 
      if (mbar_wanted > 0.) {
	xx = mass_b() / mbar_wanted - 1. ;
	cout << "Discrep. baryon mass <-> wanted bar. mass : " << xx 
	     << endl ; 
      }
      else{
	xx = mass_g() / fabs(mbar_wanted) - 1. ;
	cout << "Discrep. grav. mass <-> wanted grav. mass : " << xx 
	     << endl ; 
      }
      double xprog = ( mer > 2*mer_mass) ? 1. : 
	double(mer-mer_mass)/double(mer_mass) ; 
      xx *= xprog ; 
      double ax = .5 * ( 2. + xx ) / (1. + xx ) ; 
      double fact = pow(ax, aexp_mass) ; 
      cout << "  xprog, xx, ax, fact : " << xprog << "  " <<
	xx << "  " << ax << "  " << fact << endl ; 
      
      if ( change_ent ) {
	ent_c *= fact ; 
      }
      else {
	if (mer%4 == 0) omega *= fact ; 
      }
    }
        
    //---------------------------------------------
    // Convergence towards a given magnetic moment 
    //---------------------------------------------
    
    if (mer > mer_magmom) {
      
      // if(mer > mer_mass) {
      // 	cerr << "et_magnetisation::equilibrium()" << endl ;
      // 	cerr << "Impossible to converge to a given baryon mass" << endl ;
      // 	cerr << "and a given magnetic moment!" << endl ;
      // 	cerr << "Aborting..." << endl ;
      // 	abort() ;
      // }
      double xx = MagMom() / magmom_wanted - 1. ;
	cout << "Discrep. mag. moment <-> wanted mag. moment : " << xx 
	     << endl ; 

      double xprog = ( mer > 2*mer_magmom) ? 1. : 
	double(mer-mer_magmom)/double(mer_magmom) ; 
      xx *= xprog ; 
      double ax = .5 * ( 2. + xx ) / (1. + xx ) ; 
      double fact = pow(ax, aexp_mass) ; 
      cout << "  xprog, xx, ax, fact : " << xprog << "  " <<
	xx << "  " << ax << "  " << fact << endl ; 
      
      a_j *= fact ; 
    }
        
    //------------------------------------------------------------
    //  Relative change in enthalpy with respect to previous step 
    //------------------------------------------------------------
    
    Tbl diff_ent_tbl = diffrel( ent(), ent_prev() ) ; 
    diff_ent = diff_ent_tbl(0) ; 
    for (int l=1; l<nzet; l++) {
      diff_ent += diff_ent_tbl(l) ; 
    }
    diff_ent /= nzet ; 
    
    fichconv << "  " << log10( fabs(diff_ent) + 1.e-16 ) ;
    fichconv << "  " << log10( fabs(err_grv2) + 1.e-16 ) ;
    
    //------------------------------
    //  Recycling for the next step
    //------------------------------
    
    ent_prev = ent ; 
    logn_prev = logn ; 
    dzeta_prev = dzeta ; 
    
    fichconv << endl ;
    fichfreq << endl ;
    fichevol << endl ;
    fichconv.flush() ; 
    fichfreq.flush() ; 
    fichevol.flush() ; 
    
  } // End of main loop
  
    //=========================================================================
    // 			End of iteration
    //=========================================================================

  fichconv.close() ; 
  fichfreq.close() ; 
  fichevol.close() ; 
}














