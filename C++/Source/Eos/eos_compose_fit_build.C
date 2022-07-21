/*
 *  Methods of class Eos_compose_fit building.
 *
 *  (see file eos_compose_fit.h for documentation).
 *
 */

/*
 *   Copyright (c) 2022 Jerome Novak
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
 * Revision 1.2  2022/07/21 12:34:18  j_novak
 * Corrected units
 *
 * Revision 1.1  2022/04/15 13:39:24  j_novak
 * New class Eos_compose_fit to generate fitted EoSs from CompOSE tables.
 *
 *
 * $Header$
 *
 */

// Headers Lorene
#include "eos_compose_fit.h"
#include "scalar.h"
#include "utilitaires.h"
#include "unites.h"

namespace Lorene {
void Eos_compose_fit::read_compose_data(int& nbp, Tbl*& logh_read, Tbl*& logp_read,
					Tbl*& loge_read, Tbl*& lognb_read,
					Tbl*& gam1_read)
{
  // Files containing data and a test
  //---------------------------------
  string file_nb = tablename + "eos.nb" ;
  string file_thermo = tablename + "eos.thermo" ;

  cout << "opening " << file_nb << " and " << file_thermo << " ... " << flush ;
  ifstream in_nb(file_nb.data()) ;
  if (!in_nb) {
    cerr << "Eos_compose_fit::read_compose_data:" << endl ;
    cerr << "Problem in opening the EOS file!" << endl ;
    cerr << "While trying to open " << file_nb << endl ;
    cerr << "Aborting..." << endl ;
    abort() ;
  }
  
  // obtaining the size of the tables for memory allocation
  //-------------------------------------------------------
  int index1, index2 ;
  in_nb >> index1 >> index2 ;
  nbp = index2 - index1 + 1 ;
  assert(nbp > 0) ;

  Tbl press(nbp) ; press.set_etat_qcq() ;
  Tbl nb(nbp) ; nb.set_etat_qcq() ;
  Tbl ener(nbp) ; ener.set_etat_qcq() ; 
  
  if (logh_read != nullptr) delete logh_read ;
  logh_read = new Tbl(nbp) ;
  logh_read->set_etat_qcq() ;
  if (logp_read != nullptr) delete logp_read ;
  logp_read = new Tbl(nbp) ;
  logp_read->set_etat_qcq() ;
  if (loge_read != nullptr) delete loge_read ;
  loge_read = new Tbl(nbp) ;
  loge_read->set_etat_qcq() ;
  if (lognb_read != nullptr) delete lognb_read ;
  lognb_read = new Tbl(nbp) ;
  lognb_read->set_etat_qcq() ;
  if (gam1_read != nullptr) delete gam1_read ;
  gam1_read = new Tbl(nbp) ;
  gam1_read->set_etat_qcq() ;

  // Variables and conversion
  //-------------------------
  double nb_fm3, rho_cgs, p_cgs, p_over_nb_comp, eps_comp ;
  double dummy_x ;
  int dummy_n ;

  //# Dealing with units could be optimized...
  double rhonuc_cgs = Unites::rhonuc_si * 1e-3 ;
  double c2_cgs = Unites::c_si * Unites::c_si * 1e4 ;
  double m_neutron_MeV, m_proton_MeV ;
  
  ifstream in_p_rho (file_thermo.data()) ;
  if (!in_p_rho) {
    cerr << "Reading CompOSE data: " << endl ;
    cerr << "Problem in opening the EOS file!" << endl ;
    cerr << "While trying to open " << file_thermo << endl ;
    cerr << "Aborting..." << endl ;
    abort() ;
  }
  in_p_rho >> m_neutron_MeV >> m_proton_MeV ; //Neutron and proton masses
  in_p_rho.ignore(1000, '\n') ;

  // Conversion from MeV/fm^3 to cgs
  double p_convert = Unites::mev_si * 1.e45 * 10. ;
  // From meV/fm^3 to g/cm^3
  double eps_convert = Unites::mev_si * 1.e42 / (Unites::c_si*Unites::c_si) ;

  cout << " done. " << endl ;
  cout << "Reading the table ... " << flush ;
  
  //--------------------------------------
  // Main loop reading the CompOSE tables
  //--------------------------------------
  for (int i=0; i<nbp; i++) {
    in_nb >> nb_fm3 ;
    in_p_rho >> dummy_n >> dummy_n >> dummy_n >> p_over_nb_comp ;
    in_p_rho >> dummy_x >> dummy_x >> dummy_x >> dummy_x >> dummy_x
	     >> eps_comp ;
    in_p_rho.ignore(1000, '\n') ;
    p_cgs = p_over_nb_comp * nb_fm3 * p_convert ;
    rho_cgs = ( eps_comp + 1. ) * m_neutron_MeV * nb_fm3 * eps_convert ;
    
    if ( (nb_fm3 < 0) || (rho_cgs < 0) || (p_cgs < 0) ){
      cerr << "Eos_compose_fit::read_compose_data: " << endl ;
      cerr << "Negative value in table!" << endl ;
      cerr << "nb = " << nb_fm3 << ", rho = " << rho_cgs <<
	", p = " << p_cgs << endl ;
      cerr << "Aborting..." << endl ;
      abort() ;
    }

    press.set(i) = p_cgs / c2_cgs ;
    nb.set(i)    = nb_fm3 ;
    ener.set(i)    = rho_cgs ;

    // log-quantities in cgs units
    logp_read->set(i) = log( press(i) / rhonuc_cgs ) ;
    loge_read->set(i) = log( ener(i) / rhonuc_cgs ) ; 
    lognb_read->set(i) = log( 10.* nb(i) ) ;
  }
  cout << " done." << endl ;  
  // -----------------------------------------
  // End of the reading of the CompOSE tables
  // -----------------------------------------

  cout << "Computation of derived quantities and file outputs.. " << flush ;
  
  // Computation of various derived quantities
  //-------------------------------------------
  
  // log-enthallpy, with a shift ensuring h=10^{-14} for the lowest density
  double ww = 0. ;
  for (int i=0; i<nbp; i++) {
    double h = log( (ener(i) + press(i)) / (10 * nb(i) * rhonuc_cgs) ) ;
      
    if (i==0) { ww = h ; }
    h = h - ww + 1.e-14 ;
    
    logh_read->set(i) = log( h ) ;
  }

  // d log(p) / d log(n) -> adiabatic index \Gamma_1
  compute_derivative((*lognb_read), (*logp_read), (*gam1_read) ) ;

  cout << "done" << endl ;

}


void Eos_compose_fit::adiabatic_index_fit(int i_min, int i_max,
			 const Tbl& logh_read, const Tbl& logp_read,
			 const Tbl& loge_read, const Tbl& lognb_read,
			 const Tbl& gam1_read) {

  int nbp = logh_read.get_dim(0) ;
  double p_bound = exp(logp_read(nbp-1)) ;
  double e_bound = exp(loge_read(nbp-1)) ;
  double nb_bound = exp(lognb_read(nbp-1)) ;

  assert((mp != nullptr) && (mg != nullptr)) ;
  const Map_af& mpd = *mp ;
  int nz = mg->get_nzone() ;
  int nr = mg->get_nr(0) ;
#ifndef NDEBUG 
    for (int l=1; l<nz; l++)
      assert(mg->get_nr(l) == nr) ;
#endif
  const Coord& xx = mpd.r ; // xx = log(h) = log(log(e+p/m_B n_B))
  double x_end = (+xx)(nz-1, 0, 0, nr-1) ;
  double x_lim_max = (+xx)(nz-2, 0, 0, nr-1) ;
  double x_lim_min = (+xx)(1, 0, 0, 0) ;
  double x_0 = (+xx)(0, 0, 0, 0) ;

  hmin = exp(x_0) ;
  hmax = exp(x_end) ;
  
  if (log_p != nullptr) delete log_p ;
  if (log_e != nullptr) delete log_e ;
  if (log_nb != nullptr) delete log_nb ;
  if (log_cs2 != nullptr) delete log_cs2 ;
  
  int np_gam = 0 ;
  double mean_gam = 0. ;
  for (int i=0; i<i_min; i++) {
    mean_gam += gam1_read(i) ;
    np_gam++ ;
  }

  mean_gam /= double(np_gam) ;
  cout << "Mean Gamma_1 = " << mean_gam << endl ;
  
  // Fitting of the adiabatic index \Gamma_1
  Scalar gam1_spec(mpd) ;

  // Low density part : polytrope, \Gamma_1 = const.
  //-------------------------------------------------
  gam1_spec = mean_gam ;
  gam1_spec.std_spectral_base() ;

    // High density part : polynomial fit (polynomial regression)
  //--------------------------------------------------------------
  int ndata = nbp - i_max ;
  Tbl xdata(ndata) ; xdata.set_etat_qcq() ;
  Tbl ydata(ndata) ; ydata.set_etat_qcq() ;
  for (int i=i_max; i<nbp; i++) {
    xdata.set(i-i_max) = 2.*(logh_read(i) - x_lim_max)
      / (x_end - x_lim_max) - 1. ;
    ydata.set(i-i_max) = gam1_read(i) ;
  }

  Tbl tcoefs = poly_regression(xdata, ydata, n_coefs) ;

  // Putting coefficients to the 'Scalar' object 
  gam1_spec.set_spectral_va().coef() ;
  for (int i=0; i<=n_coefs; i++) {
    (*gam1_spec.set_spectral_va().c_cf).set(nz-1, 0, 0, i) = tcoefs(i) ;
  }
  
  // Not nice, but only simple way to update values on grid points
  if (gam1_spec.set_spectral_va().c != nullptr) {
    delete gam1_spec.set_spectral_va().c ;
    gam1_spec.set_spectral_va().c = nullptr ;
  }
  gam1_spec.set_spectral_va().coef_i() ;

  // Intermediate part : linear interpolation
  //------------------------------------------
  cout << "gam1_spec at lower boundary of fit = "
       << gam1_spec.val_grid_point(nz-1, 0, 0, 0) << endl ;
  Scalar interm(mpd) ;
  // Linear formula
  interm = mean_gam
    + (gam1_spec.val_grid_point(nz-1, 0, 0, 0) - mean_gam)
    * (xx - x_lim_min) / (x_lim_max - x_lim_min) ;
  
  gam1_spec.set_domain(nz-2) = interm.domain(nz-2) ;

  // Solution of the differential equations to compute pressure, energy, etc
  //--------------------------------------------------------------------------
  Scalar expexpx(mpd) ; // e^(e^x) , or simply (e+p)/(m_B n_B)
  expexpx = exp(exp(xx)) ;
  expexpx.std_spectral_base() ;
  Scalar expx(mpd) ; // e^x or simply h = log((e+p)/(m_B n_B))
  expx = exp(xx) ;
  expx.std_spectral_base() ;

   // rhs: f = (\Gamma_1 - 1.)/\Gamma_1 * e^x * e^(e^x)
  Scalar fff = ((gam1_spec-1.)/gam1_spec)*expx*expexpx ;
  fff.std_spectral_base() ;

  // First integration
  //-------------------
  // F = primitive(f)
  Scalar FFF = fff.primr() ;
  FFF.std_spectral_base() ;
  // Integration constant
  double Y0 = p_bound / (e_bound + p_bound) ;
  double A0 = Y0*expexpx.val_grid_point(nz-1, 0, 0, nr-1)
    - FFF.val_grid_point(nz-1, 0, 0, nr-1) ;
  // # the parameter '1.e-4' below is to avoid negative values for Y
  // # baryon density depends on it at low densities
  double diff = (1.+1.e-4)*(A0 + FFF.val_grid_point(0, 0, 0, 0)) ;
  A0 -= diff ;

  // Y is solution of Y' + Y e^x = f 
  Scalar YYY = (A0 + FFF)/expexpx ;

  // log(cs^2)
  log_cs2 = new Scalar(log(YYY*gam1_spec)) ;
  log_cs2->std_spectral_base();

  // Second integration
  //--------------------
  // \Pi' = e^x / Y
  Scalar Pidot = expx / YYY ;
  log_p = new Scalar(Pidot.primr()) ; // log_p_spec = log(p)
  // Integration constant
  double lnp0 = log(p_bound) -  log_p->val_grid_point(nz-1, 0, 0, nr-1) ;
  (*log_p) = (*log_p) + lnp0 ; 
  
  // Pressure
  Scalar spec_press = exp((*log_p)) ;
  spec_press.std_spectral_base( );

  // Energy density 
  Scalar spec_ener = spec_press*(1./YYY - 1.) ; spec_ener.std_spectral_base() ;
  log_e = new Scalar(log(spec_ener)) ;
  log_e->std_spectral_base() ; 


  // Baryon density
  Scalar spec_nbar = spec_press / (YYY * expexpx) ; 
  cout << "Relative difference in baryon density at the last point" << endl ;
  cout << "(fitted computed / original tabulated)  : "
           << fabs(1. - spec_nbar.val_grid_point(nz-1, 0, 0, nr-1) / nb_bound)
       << endl ;
  log_nb = new Scalar(log(spec_nbar)) ;
  log_nb->std_spectral_base() ; // log_nb_spec = log(n_B)

  // Matching to the low-density part polytrope
  //--------------------------------------------

  // Determining kappa and mu_0 constants 
  double kappa = spec_press.val_grid_point(nz-2, 0, 0, 0)
    / pow(spec_nbar.val_grid_point(nz-2, 0, 0, 0), mean_gam) ;

  double mu_0 = (spec_ener.val_grid_point(nz-2, 0, 0, 0)
		 - kappa/(mean_gam - 1.)
		 * pow(spec_nbar.val_grid_point(nz-2, 0, 0, 0), mean_gam))
    / spec_nbar.val_grid_point(nz-2, 0, 0, 0);

  log_p->set_spectral_va().coef_i() ;
  log_e->set_spectral_va().coef_i() ;
  log_nb->set_spectral_va().coef_i() ;
  log_cs2->set_spectral_va().coef_i() ;

  for (int i=0; i<nr; i++) {
    double h0 = exp((+xx)(0, 0, 0, i)) ;
    double n0gamm1 = (exp(h0) - mu_0)*(mean_gam - 1.) / (kappa*mean_gam) ;
    double n0 = pow(n0gamm1, 1./(mean_gam-1.)) ;
    double p0 = kappa*n0gamm1*n0 ;
    double e0 = mu_0*n0 + p0/(mean_gam - 1.) ;
    double cs2_0 = mean_gam*p0 / (p0*mean_gam/(mean_gam - 1.) + mu_0*n0 ) ;
    log_p->set_grid_point(0, 0, 0, i) = log(p0) ;
    log_e->set_grid_point(0, 0, 0, i) = log(e0) ;
    log_nb->set_grid_point(0, 0, 0, i) = log(n0) ;
    log_cs2->set_grid_point(0, 0, 0, i) = log(cs2_0) ;
  }

  cout << "... done (analytic representation)" << endl ;

}

}
