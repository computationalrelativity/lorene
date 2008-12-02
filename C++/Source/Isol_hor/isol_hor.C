/*
 *  Methods of class Isol_hor
 *
 *    (see file isol_hor.h for documentation).
 *
 */

/*
 *   Copyright (c) 2004 Jose Luis Jaramillo
 *                      Francois Limousin
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

char isol_hor_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.32  2008/12/02 15:02:21  j_novak
 * Implementation of the new constrained formalism, following Cordero et al. 2009
 * paper. The evolution eqs. are solved as a first-order system. Not tested yet!
 *
 * Revision 1.31  2006/05/24 16:55:31  f_limousin
 * Improvement of dn_comp() and dpsi_comp()
 *
 * Revision 1.30  2005/10/21 16:20:55  jl_jaramillo
 * Version for the paper JaramL05
 *
 * Revision 1.29  2005/09/13 18:33:17  f_limousin
 * New function vv_bound_cart_bin(double) for computing binaries with
 * berlin condition for the shift vector.
 * Suppress all the symy and asymy in the importations.
 *
 * Revision 1.28  2005/09/12 12:33:54  f_limousin
 * Compilation Warning - Change of convention for the angular velocity
 * Add Berlin boundary condition in the case of binary horizons.
 *
 * Revision 1.27  2005/07/08 13:15:23  f_limousin
 * Improvements of boundary_vv_cart(), boundary_nn_lapl().
 * Add a fonction to compute the departure of axisymmetry.
 *
 * Revision 1.26  2005/05/12 14:48:07  f_limousin
 * New boundary condition for the lapse : boundary_nn_lapl().
 *
 * Revision 1.25  2005/04/15 09:34:16  jl_jaramillo
 * Function "adapt_hor" for adapting the the excised surface to
 * a given surface. The zero expansion surface is not properly implemented
 *
 * Revision 1.24  2005/04/08 12:16:52  f_limousin
 * Function set_psi(). And dependance in phi.
 *
 * Revision 1.23  2005/04/03 19:48:22  f_limousin
 * Implementation of set_psi(psi_in). And minor changes to avoid warnings.
 *
 * Revision 1.22  2005/04/02 15:49:21  f_limousin
 * New choice (Lichnerowicz) for aaquad. New member data nz.
 *
 * Revision 1.21  2005/03/31 09:45:31  f_limousin
 * New functions compute_ww(...) and aa_kerr_ww().
 *
 * Revision 1.20  2005/03/30 12:08:20  f_limousin
 * Implementation of K^{ij} (Eq.(13) Of Sergio (2002)).
 *
 * Revision 1.19  2005/03/28 19:42:39  f_limousin
 * Implement the metric and A^{ij}A_{ij} of Sergio for pertubations
 * of Kerr black holes.
 *
 * Revision 1.18  2005/03/24 17:05:34  f_limousin
 * Small change
 *
 * Revision 1.17  2005/03/24 16:50:28  f_limousin
 * Add parameters solve_shift and solve_psi in par_isol.d and in function
 * init_dat(...). Implement Isolhor::kerr_perturb().
 *
 * Revision 1.16  2005/03/10 10:19:42  f_limousin
 * Add the regularisation of the shift in the case of a single black hole
 * and lapse zero on the horizon.
 *
 * Revision 1.15  2005/03/09 10:29:53  f_limousin
 * New function update_aa().
 *
 * Revision 1.14  2005/03/06 16:59:14  f_limousin
 * New function Isol_hor::aa() (the one belonging to the class
 * Time_slice_conf need to compute the time derivative of hh and thus
 * cannot work in the class Isol_hor).
 *
 * Revision 1.13  2005/03/03 15:12:17  f_limousin
 * Implement function operator>>
 *
 * Revision 1.12  2005/03/03 10:05:36  f_limousin
 * Introduction of members boost_x and boost_z.
 *
 * Revision 1.11  2005/02/07 10:35:05  f_limousin
 * Add the regularisation of the shift for the case N=0 on the horizon.
 *
 * Revision 1.10  2004/12/31 15:36:43  f_limousin
 * Add the constructor from a file and change the standard constructor.
 *
 * Revision 1.9  2004/12/29 16:14:22  f_limousin
 * Add new function beta_comp(const Isol_hor& comp).
 *
 * Revision 1.7  2004/11/05 10:57:03  f_limousin
 * Delete argument partial_save in the function sauve.
 *
 * Revision 1.6  2004/11/05 10:10:21  f_limousin
 * Construction of an isolhor with the Metric met_gamt instead
 * of a Sym_tensor.
 *
 * Revision 1.5  2004/11/03 17:16:06  f_limousin
 * Change the standart constructor. Add 4 memebers : trK, trK_point,
 * gamt and gamt_point.
 * Add also a constructor from a file.
 *
 * Revision 1.3  2004/10/29 15:44:45  jl_jaramillo
 * Remove two members
 *
 * Revision 1.2  2004/09/28 16:07:16  f_limousin
 * Remove all unused functions.
 *
 * Revision 1.1  2004/09/09 14:07:26  jl_jaramillo
 * First version
 *
 * Revision 1.1  2004/03/30 14:00:31  jl_jaramillo
 * New class Isol_hor (first version).
 *
 *
 * $Header$
 *
 */

// C headers
#include <stdlib.h>
#include <assert.h>

// Lorene headers
#include "param.h"
#include "utilitaires.h"
#include "time_slice.h"
#include "isol_hor.h"
#include "tensor.h"
#include "metric.h"
#include "evolution.h"
//#include "graphique.h"

//--------------//
// Constructors //
//--------------//
// Standard constructor
// --------------------

Isol_hor::Isol_hor(Map_af& mpi, int depth_in) : 
    Time_slice_conf(mpi, mpi.get_bvect_spher(), mpi.flat_met_spher()),
    mp(mpi), nz(mpi.get_mg()->get_nzone()), radius ((mpi.get_alpha())[0]), 
	      omega(0), boost_x(0), boost_z(0), regul(0),
  n_auto_evol(depth_in), n_comp_evol(depth_in), 
  psi_auto_evol(depth_in), psi_comp_evol(depth_in),
  dn_evol(depth_in), dpsi_evol(depth_in),
  beta_auto_evol(depth_in), beta_comp_evol(depth_in),
  aa_auto_evol(depth_in), aa_comp_evol(depth_in),
  aa_nn(depth_in), aa_quad_evol(depth_in),
  met_gamt(mpi.flat_met_spher()), gamt_point(mpi, CON, mpi.get_bvect_spher()),
  trK(mpi), trK_point(mpi), decouple(mpi){
}		  

// Constructor from conformal decomposition
// ----------------------------------------

Isol_hor::Isol_hor(Map_af& mpi, const Scalar& lapse_in, 
		   const Scalar& psi_in, const Vector& shift_in,
		   const Sym_tensor& aa_in, 
		   const Metric& metgamt, const Sym_tensor& gamt_point_in, 
		   const Scalar& trK_in, const Scalar& trK_point_in,
		   const Metric_flat& ff_in, int depth_in) 	  
    : Time_slice_conf(lapse_in, shift_in, ff_in, psi_in, metgamt.con() -
		      ff_in.con(), pow(psi_in, 6)*aa_in, trK_in, depth_in),
      mp(mpi), nz(mpi.get_mg()->get_nzone()), radius ((mpi.get_alpha())[0]), 
      omega(0), boost_x(0), boost_z(0), regul(0),
      n_auto_evol(depth_in), n_comp_evol(depth_in), 
      psi_auto_evol(depth_in), psi_comp_evol(depth_in),
      dn_evol(depth_in), dpsi_evol(depth_in),
      beta_auto_evol(depth_in), beta_comp_evol(depth_in),
      aa_auto_evol(depth_in), aa_comp_evol(depth_in), 
      aa_nn(depth_in), aa_quad_evol(depth_in),
      met_gamt(metgamt), gamt_point(gamt_point_in),
      trK(trK_in), trK_point(trK_point_in), decouple(lapse_in.get_mp()){

    // hh_evol, trk_evol
    hh_evol.update(met_gamt.con() - ff.con(), jtime, the_time[jtime]) ;
    trk_evol.update(trK, jtime, the_time[jtime]) ;
 
}

// Copy constructor
// ----------------

Isol_hor::Isol_hor(const Isol_hor& isolhor_in) 
    : Time_slice_conf(isolhor_in),
      mp(isolhor_in.mp),
      nz(isolhor_in.nz),
      radius(isolhor_in.radius),
      omega(isolhor_in.omega),
      boost_x(isolhor_in.boost_x),
      boost_z(isolhor_in.boost_z),
      regul(isolhor_in.regul),
      n_auto_evol(isolhor_in.n_auto_evol),
      n_comp_evol(isolhor_in.n_comp_evol),
      psi_auto_evol(isolhor_in.psi_auto_evol),
      psi_comp_evol(isolhor_in.psi_comp_evol),
      dn_evol(isolhor_in.dn_evol),
      dpsi_evol(isolhor_in.dpsi_evol),
      beta_auto_evol(isolhor_in.beta_auto_evol),
      beta_comp_evol(isolhor_in.beta_comp_evol),
      aa_auto_evol(isolhor_in.aa_auto_evol),
      aa_comp_evol(isolhor_in.aa_comp_evol),
      aa_nn(isolhor_in.aa_nn),
      aa_quad_evol(isolhor_in.aa_quad_evol),
      met_gamt(isolhor_in.met_gamt),
      gamt_point(isolhor_in.gamt_point),
      trK(isolhor_in.trK),
      trK_point(isolhor_in.trK_point),
      decouple(isolhor_in.decouple){
}

// Constructor from a file
// -----------------------

Isol_hor::Isol_hor(Map_af& mpi, FILE* fich, 
		   bool partial_read, int depth_in)
    : Time_slice_conf(mpi, mpi.get_bvect_spher(), mpi.flat_met_spher(), 
		      fich, partial_read, depth_in),
      mp(mpi), nz(mpi.get_mg()->get_nzone()), radius ((mpi.get_alpha())[0]), 
      omega(0), boost_x(0), boost_z(0), regul(0),
      n_auto_evol(depth_in), n_comp_evol(depth_in), 
      psi_auto_evol(depth_in), psi_comp_evol(depth_in),
      dn_evol(depth_in), dpsi_evol(depth_in),
      beta_auto_evol(depth_in), beta_comp_evol(depth_in),
      aa_auto_evol(depth_in), aa_comp_evol(depth_in), 
      aa_nn(depth_in), aa_quad_evol(depth_in),
      met_gamt(mpi.flat_met_spher()), 
      gamt_point(mpi, CON, mpi.get_bvect_spher()),
      trK(mpi), trK_point(mpi), decouple(mpi){

    fread_be(&omega, sizeof(double), 1, fich) ;
    fread_be(&boost_x, sizeof(double), 1, fich) ;
    fread_be(&boost_z, sizeof(double), 1, fich) ;
  
    int jmin = jtime - depth + 1 ; 
    int indicator ; 

    // psi_evol
    for (int j=jmin; j<=jtime; j++) {
	fread_be(&indicator, sizeof(int), 1, fich) ;	
	if (indicator == 1) {
	    Scalar psi_file(mp, *(mp.get_mg()), fich) ; 
	    psi_evol.update(psi_file, j, the_time[j]) ; 
	}
    }

    // n_auto_evol
    for (int j=jmin; j<=jtime; j++) {
	fread_be(&indicator, sizeof(int), 1, fich) ;	
	if (indicator == 1) {
	    Scalar nn_auto_file(mp, *(mp.get_mg()), fich) ; 
	    n_auto_evol.update(nn_auto_file, j, the_time[j]) ; 
	}
    }

  // psi_auto_evol
  for (int j=jmin; j<=jtime; j++) {
      fread_be(&indicator, sizeof(int), 1, fich) ;	
      if (indicator == 1) {
	  Scalar psi_auto_file(mp, *(mp.get_mg()), fich) ; 
	  psi_auto_evol.update(psi_auto_file, j, the_time[j]) ; 
      }
  }

  // beta_auto_evol
  for (int j=jmin; j<=jtime; j++) {
      fread_be(&indicator, sizeof(int), 1, fich) ;	
      if (indicator == 1) {
	  Vector beta_auto_file(mp, mpi.get_bvect_spher(), fich) ; 
	  beta_auto_evol.update(beta_auto_file, j, the_time[j]) ; 
      }
  }
  
  // met_gamt, gamt_point, trK, trK_point

  Sym_tensor met_file (mp, mp.get_bvect_spher(), fich) ;
  met_gamt = met_file ;

  Sym_tensor gamt_point_file (mp, mp.get_bvect_spher(), fich) ;
  gamt_point = gamt_point_file ;
  
  Scalar trK_file (mp, *(mp.get_mg()), fich) ;
  trK = trK_file ;
  
  Scalar trK_point_file (mp, *(mp.get_mg()), fich) ;
  trK_point = trK_point_file ;
  
  // hh_evol, trk_evol
  hh_evol.update(met_gamt.con() - ff.con(), jtime, the_time[jtime]) ;
  trk_evol.update(trK, jtime, the_time[jtime]) ;

}

			    //--------------//
			    //  Destructor  //
			    //--------------//

Isol_hor::~Isol_hor(){}


                    //-----------------------//
                    // Mutators / assignment //
                    //-----------------------//

void Isol_hor::operator=(const Isol_hor& isolhor_in) {

Time_slice_conf::operator=(isolhor_in) ;
 mp = isolhor_in.mp ;
 nz = isolhor_in.nz ;
 radius = isolhor_in.radius ;
 omega = isolhor_in.omega ;
 boost_x = isolhor_in.boost_x ;
 boost_z = isolhor_in.boost_z ;
 regul = isolhor_in.regul ;
 n_auto_evol = isolhor_in.n_auto_evol ;
 n_comp_evol = isolhor_in.n_comp_evol ;
 psi_auto_evol = isolhor_in.psi_auto_evol ;
 psi_comp_evol = isolhor_in.psi_comp_evol ;
 dn_evol = isolhor_in.dn_evol ;
 dpsi_evol = isolhor_in.dpsi_evol ;
 beta_auto_evol = isolhor_in.beta_auto_evol ;
 beta_comp_evol = isolhor_in.beta_comp_evol ;
 aa_auto_evol = isolhor_in.aa_auto_evol ;
 aa_comp_evol = isolhor_in.aa_comp_evol ;
 aa_nn = isolhor_in.aa_nn ;
 aa_quad_evol = isolhor_in.aa_quad_evol ;
 met_gamt = isolhor_in.met_gamt ;
 gamt_point = isolhor_in.gamt_point ;
 trK = isolhor_in.trK ;
 trK_point = isolhor_in.trK_point ;
 decouple = isolhor_in.decouple ;
}


                //------------------//
                //      output      //
                //------------------//


ostream& Isol_hor::operator>>(ostream& flux) const {
    
    Time_slice_conf::operator>>(flux) ; 
    
    flux << '\n' << "radius of the horizon  : " << radius << '\n' ;
    flux << "boost in x-direction   : " << boost_x << '\n' ;
    flux << "boost in z-direction   : " << boost_z << '\n' ;
    flux << "angular velocity omega : " << omega_hor() << '\n' ;
    flux << "area of the horizon    : " << area_hor() << '\n' ;
    flux << "ang. mom. of horizon   : " << ang_mom_hor() << '\n' ;
    flux << "ADM ang. mom.          : " << ang_mom_adm() << '\n' ;
    flux << "Mass of the horizon    : " << mass_hor() << '\n' ;
    flux << "ADM Mass               : " << adm_mass() << '\n' ;
   
    return flux ;     
}


                //--------------------------//
                //      Save in a file      //
                //--------------------------//


void Isol_hor::sauve(FILE* fich, bool partial_save) const {


    // Writing of quantities common to all derived classes of Time_slice
    // -----------------------------------------------------------------
    
    Time_slice_conf::sauve(fich, partial_save) ; 

    fwrite_be (&omega, sizeof(double), 1, fich) ;
    fwrite_be (&boost_x, sizeof(double), 1, fich) ;
    fwrite_be (&boost_z, sizeof(double), 1, fich) ;
    
    // Writing of quantities common to Isol_hor
    // -----------------------------------------

    int jmin = jtime - depth + 1 ; 

    // psi_evol
    for (int j=jmin; j<=jtime; j++) {
	int indicator = (psi_evol.is_known(j)) ? 1 : 0 ; 
        fwrite_be(&indicator, sizeof(int), 1, fich) ;
        if (indicator == 1) psi_evol[j].sauve(fich) ; 
    }
    
    // n_auto_evol
    for (int j=jmin; j<=jtime; j++) {
	int indicator = (n_auto_evol.is_known(j)) ? 1 : 0 ; 
        fwrite_be(&indicator, sizeof(int), 1, fich) ;
        if (indicator == 1) n_auto_evol[j].sauve(fich) ; 
    }
	 
    // psi_auto_evol
    for (int j=jmin; j<=jtime; j++) {
	int indicator = (psi_auto_evol.is_known(j)) ? 1 : 0 ; 
        fwrite_be(&indicator, sizeof(int), 1, fich) ;
        if (indicator == 1) psi_auto_evol[j].sauve(fich) ; 
    }
	 
    // beta_auto_evol
    for (int j=jmin; j<=jtime; j++) {
	int indicator = (beta_auto_evol.is_known(j)) ? 1 : 0 ; 
        fwrite_be(&indicator, sizeof(int), 1, fich) ;
        if (indicator == 1) beta_auto_evol[j].sauve(fich) ; 
    }
	
    met_gamt.con().sauve(fich) ;
    gamt_point.sauve(fich) ;    
    trK.sauve(fich) ;
    trK_point.sauve(fich) ;
}

// Accessors
// ---------

const Scalar& Isol_hor::n_auto() const {

    assert( n_auto_evol.is_known(jtime) ) ; 
    return n_auto_evol[jtime] ;   
} 

const Scalar& Isol_hor::n_comp() const {

    assert( n_comp_evol.is_known(jtime) ) ; 
    return n_comp_evol[jtime] ;   
} 

const Scalar& Isol_hor::psi_auto() const {

    assert( psi_auto_evol.is_known(jtime) ) ; 
    return psi_auto_evol[jtime] ;   
} 

const Scalar& Isol_hor::psi_comp() const {

    assert( psi_comp_evol.is_known(jtime) ) ; 
    return psi_comp_evol[jtime] ;   
} 

const Vector& Isol_hor::dnn() const {

    assert( dn_evol.is_known(jtime) ) ; 
    return dn_evol[jtime] ;   
} 

const Vector& Isol_hor::dpsi() const {

    assert( dpsi_evol.is_known(jtime) ) ; 
    return dpsi_evol[jtime] ;   
} 

const Vector& Isol_hor::beta_auto() const {

    assert( beta_auto_evol.is_known(jtime) ) ; 
    return beta_auto_evol[jtime] ;   
} 

const Vector& Isol_hor::beta_comp() const {

    assert( beta_comp_evol.is_known(jtime) ) ; 
    return beta_comp_evol[jtime] ;   
} 

const Sym_tensor& Isol_hor::aa_auto() const {

    assert( aa_auto_evol.is_known(jtime) ) ; 
    return aa_auto_evol[jtime] ;   
} 

const Sym_tensor& Isol_hor::aa_comp() const {

    assert( aa_comp_evol.is_known(jtime) ) ; 
    return aa_comp_evol[jtime] ;   
} 

void Isol_hor::set_psi(const Scalar& psi_in) {

    psi_evol.update(psi_in, jtime, the_time[jtime]) ;

    // Reset of quantities depending on Psi:
    if (p_psi4 != 0x0) {
        delete p_psi4 ;
        p_psi4 = 0x0 ;
    }
    if (p_ln_psi != 0x0) {
        delete p_ln_psi ;
        p_ln_psi = 0x0 ;
    }
    if (p_gamma != 0x0) {
        delete p_gamma ;
        p_gamma = 0x0 ;
    }
    gam_dd_evol.downdate(jtime) ;
    gam_uu_evol.downdate(jtime) ;
    k_dd_evol.downdate(jtime) ;
    k_uu_evol.downdate(jtime) ;
    adm_mass_evol.downdate(jtime) ;
    
}

void Isol_hor::set_nn(const Scalar& nn_in) {

    n_evol.update(nn_in, jtime, the_time[jtime]) ;

    hata_evol.downdate(jtime) ;
    aa_quad_evol.downdate(jtime) ;
    k_dd_evol.downdate(jtime) ;
    k_uu_evol.downdate(jtime) ;
}

void Isol_hor::set_gamt(const Metric& gam_tilde) {

    if (p_tgamma != 0x0) {
        delete p_tgamma ;
        p_tgamma = 0x0 ;
    }
    if (p_gamma != 0x0) {
        delete p_gamma ;
        p_gamma = 0x0 ;
    }

    met_gamt = gam_tilde ;
 
    gam_dd_evol.downdate(jtime) ;
    gam_uu_evol.downdate(jtime) ;
    k_dd_evol.downdate(jtime) ;
    k_uu_evol.downdate(jtime) ;
    hh_evol.downdate(jtime) ;

    hh_evol.update(gam_tilde.con() - ff.con(), jtime, the_time[jtime]) ;

}



// Import the lapse from the companion (Bhole case)

void Isol_hor::n_comp(const Isol_hor& comp) {

    double ttime = the_time[jtime] ;    

    Scalar temp (mp) ;
    temp.import(comp.n_auto()) ;
    temp.std_spectral_base() ;
    n_comp_evol.update(temp, jtime, ttime) ;
    n_evol.update(temp + n_auto(), jtime, ttime) ;
     
    Vector dn_comp (mp, COV, mp.get_bvect_cart()) ;
    dn_comp.set_etat_qcq() ;
    Vector auxi (comp.n_auto().derive_cov(comp.ff)) ;
    auxi.dec_dzpuis(2) ;
    auxi.change_triad(auxi.get_mp().get_bvect_cart()) ;
    for (int i=1 ; i<=3 ; i++){
      if (auxi(i).get_etat() != ETATZERO)
    	auxi.set(i).raccord(3) ;
    }

    auxi.change_triad(mp.get_bvect_cart()) ;
    assert ( *(auxi.get_triad()) == *(dn_comp.get_triad())) ;

    for (int i=1 ; i<=3 ; i++){
    dn_comp.set(i).import(auxi(i)) ;
    dn_comp.set(i).set_spectral_va().set_base(auxi(i).get_spectral_va().
					      get_base()) ;
    }
    dn_comp.inc_dzpuis(2) ;
    dn_comp.change_triad(mp.get_bvect_spher()) ;
    /*    
    Vector dn_comp_zec (n_comp().derive_cov(ff)) ;
    for (int i=1 ; i<=3 ; i++)
      for (int l=nz-1 ; l<=nz-1 ; l++) {
	  if (dn_comp.set(i).get_etat() == ETATQCQ)
	    dn_comp.set(i).set_domain(l) = dn_comp_zec(i).domain(l) ;
      }
    */
    dn_evol.update(n_auto().derive_cov(ff) + dn_comp, jtime, ttime) ;


    /*
    Scalar tr_K (mp) ;
 
    Mtbl mxabs (mp.xa) ;
    Mtbl myabs (mp.ya) ;
    Mtbl mzabs (mp.za) ;
    Scalar comp_r (mp) ;
    comp_r.annule_hard() ;
    for (int l=1 ; l<mp.get_mg()->get_nzone() ; l++) {
	int nr = mp.get_mg()->get_nr (l) ;
	if (l==mp.get_mg()->get_nzone()-1)
	    nr -- ;	
	int np = mp.get_mg()->get_np (l) ;
	int nt = mp.get_mg()->get_nt (l) ;
	double xabs, yabs, zabs, air, theta, phi ;
    
    for (int k=0 ; k<np ; k++)
	for (int j=0 ; j<nt ; j++)
	    for (int i=0 ; i<nr ; i++) {
		
		xabs = mxabs (l, k, j, i) ;
		yabs = myabs (l, k, j, i) ;
		zabs = mzabs (l, k, j, i) ;

		// coordinates of the point in 2 :
		comp.mp.convert_absolute 
		    (xabs, yabs, zabs, air, theta, phi) ;
		comp_r.set_grid_point(l,k,j,i) = air ;
	    }
    }

    Scalar fact(mp) ;
    fact = 0.0000000000000001 ;

//    Scalar fact(mp) ;
//    fact = 0.7*gam().radial_vect().divergence(gam()) ;
//    fact.dec_dzpuis(2) ;

    tr_K = 1/mp.r/mp.r ;
    tr_K = tr_K * fact ;
    tr_K += fact/comp_r/comp_r ;
    tr_K.std_spectral_base() ;
    tr_K.annule(0, 0) ;
    tr_K.raccord(1) ;
    tr_K.inc_dzpuis(2) ;
    trk_evol.update(tr_K, jtime, the_time[jtime]) ;
*/
}

// Import the conformal factor from the companion (Bhole case)

void Isol_hor::psi_comp (const Isol_hor& comp) {
  
    double ttime = the_time[jtime] ;    
    
    Scalar temp (mp) ;
    temp.import(comp.psi_auto()) ;
    temp.std_spectral_base() ;
    psi_comp_evol.update(temp, jtime, ttime) ;
    psi_evol.update(temp + psi_auto(), jtime, ttime) ;
    
    Vector dpsi_comp (mp, COV, mp.get_bvect_cart()) ;
    dpsi_comp.set_etat_qcq() ;
    Vector auxi (comp.psi_auto().derive_cov(comp.ff)) ;
    auxi.dec_dzpuis(2) ;
    auxi.change_triad(auxi.get_mp().get_bvect_cart()) ;
    for (int i=1 ; i<=3 ; i++){
      if (auxi(i).get_etat() != ETATZERO)
        auxi.set(i).raccord(3) ;
    }

    auxi.change_triad(mp.get_bvect_cart()) ;
    assert ( *(auxi.get_triad()) == *(dpsi_comp.get_triad())) ;

    for (int i=1 ; i<=3 ; i++){
      dpsi_comp.set(i).import(auxi(i)) ;
      dpsi_comp.set(i).set_spectral_va().set_base(auxi(i).get_spectral_va().
						  get_base()) ;
    }
    dpsi_comp.inc_dzpuis(2) ;
    dpsi_comp.change_triad(mp.get_bvect_spher()) ;
    /*    
    Vector dpsi_comp_zec (psi_comp().derive_cov(ff)) ;
    for (int i=1 ; i<=3 ; i++)
      for (int l=nz-1 ; l<=nz-1 ; l++) {
	if (dpsi_comp.set(i).get_etat() == ETATQCQ)
	  dpsi_comp.set(i).set_domain(l) = dpsi_comp_zec(i).domain(l) ;
      }
    */
    
    dpsi_evol.update(psi_auto().derive_cov(ff) + dpsi_comp, jtime, ttime) ;

}

void Isol_hor::beta_comp (const Isol_hor& comp) {
  
    double ttime = the_time[jtime] ;    
    
    Vector tmp_vect (mp, CON, mp.get_bvect_cart()) ;
    Vector shift_comp (comp.beta_auto()) ;
    shift_comp.change_triad(comp.mp.get_bvect_cart()) ;
    shift_comp.change_triad(mp.get_bvect_cart()) ;
    assert (*(shift_comp.get_triad()) == *(tmp_vect.get_triad())) ;

    tmp_vect.set(1).import(shift_comp(1)) ;
    tmp_vect.set(2).import(shift_comp(2)) ;
    tmp_vect.set(3).import(shift_comp(3)) ;
    tmp_vect.std_spectral_base() ;
    tmp_vect.change_triad(mp.get_bvect_spher()) ;
    
    beta_comp_evol.update(tmp_vect, jtime,ttime) ;
    beta_evol.update(beta_auto() + beta_comp(), jtime, ttime) ;
}

//Initialisation to Schwartzchild
void Isol_hor::init_bhole () {
    
    double ttime = the_time[jtime] ;    
    Scalar auxi(mp) ;
    
    // Initialisation of the lapse different of zero on the horizon
    // at the first step
    auxi = 0.5 - 0.5/mp.r ;
    auxi.annule(0, 0);
    auxi.set_dzpuis(0) ;
    
    Scalar temp(mp) ;
    temp = auxi;
    temp.std_spectral_base() ;
    temp.raccord(1) ;
    n_auto_evol.update(temp, jtime, ttime) ;

    temp = 0.5 ;
    temp.std_spectral_base() ;
    n_comp_evol.update(temp, jtime, ttime) ;
    n_evol.update(n_auto() + n_comp(), jtime, ttime) ;
  
    auxi = 0.5 + radius/mp.r ;
    auxi.annule(0, 0);
    auxi.set_dzpuis(0) ;
    temp = auxi;
    temp.std_spectral_base() ;
    temp.raccord(1) ;
    psi_auto_evol.update(temp, jtime, ttime) ;

    temp = 0.5 ;
    temp.std_spectral_base() ;
    psi_comp_evol.update(temp, jtime, ttime) ;
    psi_evol.update(psi_auto() + psi_comp(), jtime, ttime) ;
    
    dn_evol.update(nn().derive_cov(ff), jtime, ttime) ;
    dpsi_evol.update(psi().derive_cov(ff), jtime, ttime) ;
    
    Vector temp_vect1(mp, CON, mp.get_bvect_spher()) ;
    temp_vect1.set(1) = 0.0/mp.r/mp.r ;
    temp_vect1.set(2) = 0. ;
    temp_vect1.set(3) = 0. ;
    temp_vect1.std_spectral_base() ;

    Vector temp_vect2(mp, CON, mp.get_bvect_spher()) ;
    temp_vect2.set_etat_zero() ;    

    beta_auto_evol.update(temp_vect1, jtime, ttime) ;
    beta_comp_evol.update(temp_vect2, jtime, ttime) ;
    beta_evol.update(temp_vect1, jtime, ttime) ;    
}

void Isol_hor::init_met_trK() {
 
  Metric flat (mp.flat_met_spher()) ;
  met_gamt = flat ;

  gamt_point.set_etat_zero() ;
  trK.set_etat_zero() ;
  trK_point.set_etat_zero() ;
 
}


void Isol_hor::init_bhole_seul () {
    
    double ttime = the_time[jtime] ;    
    Scalar auxi(mp) ;
    
    auxi = (1-radius/mp.r)/(1+radius/mp.r) ;
    auxi.annule(0, 0);
    auxi.set_outer_boundary((*mp.get_mg()).get_nzone(), 1.) ;
    auxi.set_dzpuis(0) ;

    Scalar temp(mp) ;
    temp = auxi;
    temp.std_spectral_base() ;
    temp.raccord(1) ;
    n_auto_evol.update(temp, jtime, ttime) ;

    temp.set_etat_zero() ;
    n_comp_evol.update(temp, jtime, ttime) ;
    n_evol.update(temp, jtime, ttime) ;
 
    
    auxi = 1 + radius/mp.r ;
    auxi.annule(0, 0);
    auxi.set_outer_boundary((*mp.get_mg()).get_nzone(), 1.) ;
    auxi.set_dzpuis(0) ;
  
    temp = auxi;
    temp.std_spectral_base() ;
    temp.raccord(1) ;
    psi_auto_evol.update(temp, jtime, ttime) ;
    temp.set_etat_zero() ;
    psi_comp_evol.update(temp, jtime, ttime) ;
    psi_evol.update(temp, jtime, ttime) ;
    
    dn_evol.update(nn().derive_cov(ff), jtime, ttime) ;
    dpsi_evol.update(psi().derive_cov(ff), jtime, ttime) ;

    Vector temp_vect(mp, CON, mp.get_bvect_spher()) ;
    temp_vect.set_etat_zero() ;
    beta_auto_evol.update(temp_vect, jtime, ttime) ;
    beta_comp_evol.update(temp_vect, jtime, ttime) ;
    beta_evol.update(temp_vect, jtime, ttime) ;    	

}		   

void Isol_hor::update_aa() {
	
  Sym_tensor aa_new (mp, CON, mp.get_bvect_spher()) ;
  int nnt = mp.get_mg()->get_nt(1) ;
  int nnp = mp.get_mg()->get_np(1) ;

  int check ;
  check = 0 ;
  for (int k=0; k<nnp; k++)
    for (int j=0; j<nnt; j++){
      if (nn().val_grid_point(1, k, j , 0) < 1e-12){
	check = 1 ;
	break ;
      }
    }
  
  if (check == 0)
    aa_new = ( beta().ope_killing_conf(met_gamt) + gamt_point ) 
      / (2.* nn()) ;            
  else {
    regul = regularise_one() ;
    cout << "regul = " << regul << endl ;
    Scalar nn_sxpun (division_xpun (Cmp(nn()), 0)) ;
    aa_new = beta().ope_killing_conf(met_gamt) + gamt_point ;
    
    Scalar auxi (mp) ;
    for (int i=1 ; i<=3 ; i++)
	for (int j=i ; j<=3 ; j++) {
   	    auxi = aa_new(i, j) ;
	    auxi = division_xpun (auxi, 0) ;
	    aa_new.set(i,j) = auxi / nn_sxpun / 2. ;
	}
  }
  
  Sym_tensor hata_new = aa_new*psi4()*psi()*psi() ;
  set_hata(hata_new) ;  // set aa to aa_new and delete values of 
                    // k_uu and k_dd.
  Sym_tensor aa_dd (aa_new.up_down(met_gamt)) ;
  Scalar aquad (contract(aa_dd, 0, 1, aa_new, 0, 1)*psi4()*psi4()*psi4()) ;
  aa_quad_evol.update(aquad, jtime, the_time[jtime]) ;
}

const Scalar& Isol_hor::aa_quad() const {

  if (!aa_quad_evol.is_known(jtime) ) {
    Sym_tensor aa_dd (aa().up_down(met_gamt)) ;
    Scalar aquad (contract(aa_dd, 0, 1, aa(), 0, 1)*psi4()*psi4()*psi4()) ;
    aa_quad_evol.update(aquad, jtime, the_time[jtime]) ;
  } 

    return aa_quad_evol[jtime] ;   
} 

void Isol_hor::met_kerr_perturb() {

    Sym_tensor gamm (gam().cov()) ;
    Sym_tensor gamt (gamm / gamm(3,3)) ;

    Metric metgamt (gamt) ;
    met_gamt = metgamt ;

    Scalar psi_perturb (pow(gamm(3,3), 0.25)) ;
    psi_perturb.std_spectral_base() ;
    set_psi(psi_perturb) ;

    cout << "met_gamt" << endl << norme(met_gamt.cov()(1,1)) << endl 
	 << norme(met_gamt.cov()(2,1)) << endl << norme(met_gamt.cov()(3,1)) 
	 << endl << norme(met_gamt.cov()(2,2)) << endl 
	 << norme(met_gamt.cov()(3,2)) << endl << norme(met_gamt.cov()(3,3)) 
	 << endl ;
    cout << "determinant" << norme(met_gamt.determinant()) << endl ;

    hh_evol.update(met_gamt.con() - ff.con(), jtime, the_time[jtime]) ;
}


void Isol_hor::aa_kerr_ww(double mm, double aaa) {
  
  Scalar rr(mp) ;
  rr = mp.r ;
  Scalar cost (mp) ;
  cost = mp.cost ;
  Scalar sint (mp) ;
  sint = mp.sint ;
  
  // rbl
  Scalar rbl = rr + mm + (mm*mm - aaa*aaa) / (4*rr) ;
  
  // sigma inverse
  Scalar sigma_inv = 1. / (rbl*rbl + aaa*aaa*cost*cost) ;
  sigma_inv.set_domain(0) = 1. ;

  // ww perturbation
  Scalar ww_pert (mp) ;
  ww_pert = - 1*(mm*aaa*aaa*aaa*pow(sint, 4.)*cost) * sigma_inv ;
  ww_pert.set_spectral_va().set_base_r(0,R_CHEBPIM_P) ;
  for (int l=1; l<nz-1; l++)
    ww_pert.set_spectral_va().set_base_r(l,R_CHEB) ;
  ww_pert.set_spectral_va().set_base_r(nz-1,R_CHEBU) ;
  ww_pert.set_spectral_va().set_base_t(T_COSSIN_CI) ;
  ww_pert.set_spectral_va().set_base_p(P_COSSIN) ;

  // Quadratic part A^{ij]A_{ij}
  // Lichnerowicz choice
  //----------------------------

  // BY - BY
  Vector dw_by (mp, COV, mp.get_bvect_spher()) ;
  dw_by.set(1) = 0. ;
  dw_by.set(2) = 3 * aaa*mm*sint*sint*sint / rr ;
  dw_by.set(3) = 0. ;
  dw_by.set(2).set_spectral_va().set_base_r(0,R_CHEBPIM_P) ;
  for (int l=1; l<nz-1; l++)
    dw_by.set(2).set_spectral_va().set_base_r(l,R_CHEB) ;
  dw_by.set(2).set_spectral_va().set_base_r(nz-1,R_CHEBU) ;
  dw_by.set(2).set_spectral_va().set_base_t(T_COSSIN_SI) ;
  dw_by.set(2).set_spectral_va().set_base_p(P_COSSIN) ;

  Scalar aquad_1 = 2*contract(dw_by, 0, dw_by.up_down(ff), 0) * 
      gam_dd()(3,3) / gam_dd()(1,1) ;
  aquad_1.div_rsint_dzpuis(1) ;
  aquad_1.div_rsint_dzpuis(2) ;
  aquad_1.div_rsint_dzpuis(3) ;
  aquad_1.div_rsint_dzpuis(4) ;

  // BY - dw_pert
  Vector dw_pert(ww_pert.derive_con(ff)) ;
  Scalar aquad_2 = 4*contract(dw_by, 0, dw_pert, 0) *
      gam_dd()(3,3) / gam_dd()(1,1) ;

  aquad_2.div_rsint_dzpuis(3) ;
  aquad_2.div_rsint_dzpuis(4) ;
  aquad_2.div_rsint() ;
  aquad_2.div_rsint() ;

  // dw_pert - dw_pert
  Scalar aquad_3 = 2*contract(dw_pert, 0, dw_pert.up_down(ff), 0) *
      gam_dd()(3,3) / gam_dd()(1,1) ;

  aquad_3.div_rsint() ;
  aquad_3.div_rsint() ;
  aquad_3.div_rsint() ;
  aquad_3.div_rsint() ;

  // Total
  Scalar aquad = aquad_1 + aquad_2 + aquad_3 ;

  aquad.set_domain(0) = 0. ;
  Base_val sauve_base (aquad.get_spectral_va().get_base()) ;
  
  aquad = aquad * pow(gam_dd()(1,1), 2.) * pow(gam_dd()(3,3), -2.) ;
  aquad.set_spectral_va().set_base(sauve_base) ;
 
/*
  cout << "norme de aquad" << endl << norme(aquad) << endl ;
  cout << "norme de aa_quad" << endl << norme(aa_quad()) << endl ;

  des_meridian (aquad, 0, 4, "aquad", 1) ;
  des_meridian (aa_quad(), 0, 4, "aa_quad()", 2) ;
  des_meridian (aa_quad()-aquad, 0, 4, "diff aa_quad", 3) ;
  arrete() ;
*/

  aa_quad_evol.update(aquad, jtime, the_time[jtime]) ;
  

  // Extrinsic curvature A^{ij} and A_{ij}
  // Dynamical choice
  // -------------------------------------

    Scalar s_r (ww_pert.derive_cov(ff)(2)) ;
    s_r = - s_r * gam().cov()(3,3) / gam().cov()(1,1) ;
    s_r.div_rsint() ;

    Scalar temp = dw_by(2) ;
    temp = - temp *  gam().cov()(3,3) / gam().cov()(1,1) ;
    temp.div_rsint_dzpuis(2) ;

    s_r = s_r + temp ;
    s_r.annule_domain(0) ;

    Scalar s_t (ww_pert.derive_cov(ff)(1)) ;
    s_t = s_t * gam().cov()(3,3) / gam().cov()(1,1)  ;
    s_t.div_rsint() ;

    temp = dw_by(1) ;
    temp = temp *  gam().cov()(3,3) / gam().cov()(1,1) ;
    temp.div_rsint_dzpuis(2) ;

    s_t = s_t + temp ;
    s_t.annule_domain(0) ;


    Vector ss (mp, CON, mp.get_bvect_spher()) ;
    ss.set(1) = s_r ;
    ss.set(2) = s_t ;
    ss.set(3) = 0. ;

    Sym_tensor aij (mp, CON, mp.get_bvect_spher()) ;
    aij.set(1,1) = 0. ;
    aij.set(2,1) = 0. ;
    aij.set(2,2) = 0. ;
    aij.set(3,3) = 0. ;
    aij.set(3,1) = s_r ;
    aij.set(3,1).div_rsint() ;
    aij.set(3,2) = s_t ;
    aij.set(3,2).div_rsint() ;

    Base_val base_31 (aij(3,1).get_spectral_va().get_base()) ;
    Base_val base_32 (aij(3,2).get_spectral_va().get_base()) ;

    aij.set(3,1) = aij(3,1) * pow(gam_dd()(1,1), 5./3.) 
	                    * pow(gam_dd()(3,3), -5./3.) ;
    aij.set(3,1) = aij(3,1) * pow(psi(), -6.) ;
    aij.set(3,1).set_spectral_va().set_base(base_31) ;
    aij.set(3,2) = aij(3,2) * pow(gam_dd()(1,1), 5./3.) 
	                    * pow(gam_dd()(3,3), -5./3.) ;
    aij.set(3,2) = aij(3,2) * pow(psi(), -6.) ;
    aij.set(3,2).set_spectral_va().set_base(base_32) ;

    /*
    cout << "norme de A(3,1)" << endl << norme(aij(3,1)) << endl ;
    cout << "norme de A(3,2)" << endl << norme(aij(3,2)) << endl ;
	
    cout << "norme de A_init(3,1)" << endl << norme(aa()(3,1)) << endl ;
    cout << "norme de A_init(3,2)" << endl << norme(aa()(3,2)) << endl ;

    des_meridian(aij(3,1), 0., 4., "aij(3,1)", 0) ;
    des_meridian(aa()(3,1), 0., 4., "aa_init(3,1)", 1) ;
    des_meridian(aa()(3,1)-aij(3,1), 0., 4., "diff_aa(3,1)", 2) ;
    des_meridian(aij(3,2), 0., 4., "aij(3,2)", 3) ;
    des_meridian(aa()(3,2), 0., 4., "aa_init(3,2)", 4) ;
    des_meridian(aa()(3,2)-aij(3,2), 0., 4., "diff_aa(3,2)", 5) ;
    arrete() ;
    */
    Sym_tensor hataij = aij*psi4()*psi()*psi() ;
    hata_evol.update(hataij, jtime, the_time[jtime]) ;
    Sym_tensor kij (aij) ;
    kij = kij * pow(gam().determinant(), -1./3.) ;
    kij.std_spectral_base() ;
    k_uu_evol.update(kij, jtime, the_time[jtime]) ;
    k_dd_evol.update(kij.up_down(gam()), jtime, the_time[jtime]) ;

}

double Isol_hor::axi_break() const {

    Vector phi (ff.get_mp(), CON, *(ff.get_triad()) ) ;

    Scalar tmp (ff.get_mp() ) ;
    tmp = 1 ;
    tmp.std_spectral_base() ;
    tmp.mult_rsint() ;

    phi.set(1) = 0. ;
    phi.set(2) = 0. ;
    phi.set(3) = tmp ; 
    
    Sym_tensor q_uu ( gam_uu() - gam().radial_vect() * gam().radial_vect() ) ;
    Sym_tensor q_dd ( q_uu.up_down(gam()) ) ;
		          
    Sym_tensor L_phi_q_dd ( q_dd.derive_lie( phi) ) ;
    Sym_tensor L_phi_q_uu ( contract(contract(L_phi_q_dd, 0, q_uu, 0), 0,  q_uu,0) ) ;
    
 
    Scalar integrand ( contract(  L_phi_q_dd, 0, 1, L_phi_q_uu, 0, 1 ) * darea_hor()  ) ;
    
    double axibreak = mp.integrale_surface(integrand, 1.0000000001)/ 
	(4 * M_PI * radius_hor()* radius_hor() ) ;
    
    return axibreak ;

}

double fonc_expansion(double rr, const Param& par_expansion) {

  Scalar expa = par_expansion.get_scalar(0) ;
  double theta = par_expansion.get_double(0) ;
  double phi = par_expansion.get_double(1) ;
 
  return expa.val_point(rr, theta, phi) ;
    
}
void Isol_hor::adapt_hor(double c_min, double c_max) {

  Scalar expa (expansion()) ;
  Scalar app_hor(mp) ;
  app_hor.annule_hard() ;
  int nitmax = 200 ; 
  int nit ; 

  double precis = 1.e-13 ;

  // Calculation of the radius of the apparent horizon for each (theta, phi)
  // -----------------------------------------------------------------------

  for (int nt=0; nt<mp.get_mg()->get_nt(1); nt++)
    for (int np=0; np<mp.get_mg()->get_np(1); np++) {

      double theta = mp.get_mg()->get_grille3d(1)->tet[nt] ;
      double phi = mp.get_mg()->get_grille3d(1)->phi[np] ;

      Param par_expansion ;
      par_expansion.add_scalar(expa, 0) ;
      par_expansion.add_double(theta, 0) ;
      par_expansion.add_double(phi, 1) ;
      double r_app_hor = zerosec_b(fonc_expansion, par_expansion, c_min, c_max, 
				 precis, nitmax, nit) ;
   
      app_hor.set_grid_point(1, np, nt, 0) = r_app_hor ;
    }
  
  // Transformation of the 3-metric and extrinsic curvature 
  // ------------------------------------------------------
  
  Scalar rr (mp) ;
  rr = mp.r ;
  
  Scalar trans_11 (mp) ;
  Scalar r_new (mp) ;
  r_new.annule_hard() ;
  //  trans_11.annule_hard() ;
  for (int l=1; l<nz; l++)
    for (int nr=0; nr<mp.get_mg()->get_nr(1); nr++)
      for (int nt=0; nt<mp.get_mg()->get_nt(1); nt++)
	for (int np=0; np<mp.get_mg()->get_np(1); np++) {
	  r_new.set_grid_point(l, np, nt, nr) = rr.val_grid_point(l, np, nt, nr) -
	    app_hor.val_grid_point(1, np, nt, 0) + 1 ;
	  //	  trans_11.set_grid_point(l, np, nt, nr) = 1. / 
	  //  app_hor.val_grid_point(1, np, nt, 0) ;                 // !
	}
  r_new.std_spectral_base() ;
      
  Itbl comp(2) ;
  comp.set(0) = CON ;
  comp.set(1) = COV ;

  Scalar trans_12 (r_new.dsdt()) ;
  trans_12.div_r() ;
  Scalar trans_13 (r_new.stdsdp()) ;
  trans_13.div_r() ;                                        
  for (int nr=0; nr<mp.get_mg()->get_nr(1); nr++)
    for (int nt=0; nt<mp.get_mg()->get_nt(1); nt++)
      for (int np=0; np<mp.get_mg()->get_np(1); np++) {
	trans_12.set_grid_point(nz-1, np, nt, nr) = trans_12.val_grid_point(1, np, nt, nr) ;  
	trans_13.set_grid_point(nz-1, np, nt, nr) = trans_13.val_grid_point(1, np, nt, nr) ;  
      }
      
  // Transformation matrix
  Tensor trans (mp, 2, comp, mp.get_bvect_spher()) ;
  trans.set(1,1) = 1 ;
  trans.set(1,2) = 0;//trans_12 ;
  trans.set(1,3) = 0;//trans_13 ;
  trans.set(2,2) = 1. ;  
  trans.set(3,3) = 1. ;
  trans.set(2,1) = 0. ;
  trans.set(3,1) = 0. ;
  trans.set(3,2) = 0. ;
  trans.set(2,3) = 0. ;
  trans.std_spectral_base() ;

  cout << "trans(1,3)" << endl << norme(trans(1,3)) << endl ;

  Sym_tensor gamma_uu (gam().con()) ;
  Sym_tensor kk_uu (k_uu()) ;

  gamma_uu = contract(gamma_uu, 0, 1, trans * trans, 1, 3) ;
  kk_uu = contract(kk_uu, 0, 1, trans * trans, 1, 3) ;
  
  Sym_tensor copie_gamma (gamma_uu) ;
  Sym_tensor copie_kk (kk_uu) ;

  // dz_puis set to zero
  kk_uu.dec_dzpuis(2) ;
  for(int i=1; i<=3; i++)
    for(int j=i; j<=3; j++){
      kk_uu.set(i,j).annule_hard() ;
      gamma_uu.set(i,j).annule_hard() ;
    }

  copie_kk.dec_dzpuis(2) ;

  Scalar expa_trans(mp) ;
  expa_trans.annule_hard() ;
  expa.dec_dzpuis(2) ;
  
  /*
  copie_gamma.set(2,2).div_r() ;
  copie_gamma.set(2,2).div_r() ;
  copie_gamma.set(3,3).div_rsint() ;
  copie_gamma.set(3,3).div_rsint() ;
  copie_gamma.set(1,2).div_r() ;
  copie_gamma.set(1,3).div_rsint() ;
  //  gamma_uu.set(2,3).div_r() ;           
  //  gamma_uu.set(2,3).div_rsint() ;       
  */  

  //Importation
  for(int i=1; i<=3; i++)
    for(int j=i; j<=3; j++)
      for (int l=1; l<nz; l++)
	for (int nr=0; nr<mp.get_mg()->get_nr(1); nr++)
	  for (int nt=0; nt<mp.get_mg()->get_nt(1); nt++)
	    for (int np=0; np<mp.get_mg()->get_np(1); np++) {
	      
	      double theta = mp.get_mg()->get_grille3d(1)->tet[nt] ;
	      double phi = mp.get_mg()->get_grille3d(1)->phi[np] ;
	      double r_inv = rr.val_grid_point(l, np, nt, nr) +
		app_hor.val_grid_point(1, np, nt, 0) - 1. ;
	      	      
	      gamma_uu.set(i,j).set_grid_point(l, np, nt, nr) = 
		copie_gamma(i,j).val_point(r_inv, theta, phi) ;
	      kk_uu.set(i,j).set_grid_point(l, np, nt, nr) = 
		copie_kk(i,j).val_point(r_inv, theta, phi) ;

	      expa_trans.set_grid_point(l, np, nt, nr) = expa.val_point(r_inv, theta, phi) ;
	    }
  kk_uu.std_spectral_base() ;                 // Save the base?
  gamma_uu.std_spectral_base() ;
  expa_trans.std_spectral_base() ;
  
  for (int l=1; l<nz; l++)
    for (int nr=0; nr<mp.get_mg()->get_nr(1); nr++)
      for (int nt=0; nt<mp.get_mg()->get_nt(1); nt++)
	for (int np=0; np<mp.get_mg()->get_np(1); np++) {
	  gamma_uu.set(1,2).set_grid_point(l,np,nt,nr) =  gamma_uu.set(1,2).val_grid_point(l,np,nt,nr) 
	     / (1 +app_hor.val_grid_point(1, np, nt, 0)/ rr.val_grid_point(l, np, nt, nr) 
		- 1 / rr.val_grid_point(l, np, nt, nr) ) ;
	  gamma_uu.set(1,3).set_grid_point(l,np,nt,nr) =  gamma_uu.set(1,3).val_grid_point(l,np,nt,nr) 
	     / (1 +app_hor.val_grid_point(1, np, nt, 0)/ rr.val_grid_point(l, np, nt, nr) 
		- 1 / rr.val_grid_point(l, np, nt, nr) ) ;
	  gamma_uu.set(2,2).set_grid_point(l,np,nt,nr) =  gamma_uu.set(2,2).val_grid_point(l,np,nt,nr) 
	     / (1 +app_hor.val_grid_point(1, np, nt, 0)/ rr.val_grid_point(l, np, nt, nr) 
		- 1 / rr.val_grid_point(l, np, nt, nr) ) 
	      / (1 +app_hor.val_grid_point(1, np, nt, 0)/ rr.val_grid_point(l, np, nt, nr) 
		 - 1 / rr.val_grid_point(l, np, nt, nr) ) ;
	  gamma_uu.set(2,3).set_grid_point(l,np,nt,nr) =  gamma_uu.set(2,3).val_grid_point(l,np,nt,nr) 
	     / (1 +app_hor.val_grid_point(1, np, nt, 0)/ rr.val_grid_point(l, np, nt, nr) 
		- 1 / rr.val_grid_point(l, np, nt, nr) ) 
	      / (1 +app_hor.val_grid_point(1, np, nt, 0)/ rr.val_grid_point(l, np, nt, nr) 
		 - 1 / rr.val_grid_point(l, np, nt, nr) ) ;
	  gamma_uu.set(3,3).set_grid_point(l,np,nt,nr) =  gamma_uu.set(3,3).val_grid_point(l,np,nt,nr) 
	     / (1 +app_hor.val_grid_point(1, np, nt, 0)/ rr.val_grid_point(l, np, nt, nr) 
		- 1 / rr.val_grid_point(l, np, nt, nr) ) 
	      / (1 +app_hor.val_grid_point(1, np, nt, 0)/ rr.val_grid_point(l, np, nt, nr) 
		 - 1 / rr.val_grid_point(l, np, nt, nr) ) ;
	}
	  
  /*
  gamma_uu.set(2,2).mult_r() ;
  gamma_uu.set(2,2).mult_r() ;
  gamma_uu.set(3,3).mult_rsint() ;
  gamma_uu.set(3,3).mult_rsint() ;
  gamma_uu.set(1,2).mult_r() ;
  gamma_uu.set(1,3).mult_rsint() ;
  //  gamma_uu.set(2,3).mult_r() ;                
  //  gamma_uu.set(2,3).mult_rsint() ;            
  */
  


  cout << "gamma_uu(1,1)" << endl << norme(gamma_uu(1,1)) << endl ;
  cout << "gamma_uu(2,1)" << endl << norme(gamma_uu(2,1)) << endl ;
  cout << "gamma_uu(3,1)" << endl << norme(gamma_uu(3,1)) << endl ;
  cout << "gamma_uu(2,2)" << endl << norme(gamma_uu(2,2)) << endl ;
  cout << "gamma_uu(3,2)" << endl << norme(gamma_uu(3,2)) << endl ;
  cout << "gamma_uu(3,3)" << endl << norme(gamma_uu(3,3)) << endl ;

  kk_uu.inc_dzpuis(2) ;
  expa_trans.inc_dzpuis(2) ;

  Metric gamm (gamma_uu) ;

  // Updates
  gam_uu_evol.update(gamma_uu, jtime, the_time[jtime]) ;
  gam_dd_evol.update(gamm.cov(), jtime, the_time[jtime]) ;
  k_uu_evol.update(kk_uu, jtime, the_time[jtime]) ;
  
    if (p_psi4 != 0x0) {
        delete p_psi4 ;
        p_psi4 = 0x0 ;
    }
    if (p_ln_psi != 0x0) {
        delete p_ln_psi ;
        p_ln_psi = 0x0 ;
    }
    if (p_gamma != 0x0) {
        delete p_gamma ;
        p_gamma = 0x0 ;
    }
    if (p_tgamma != 0x0) {
        delete p_tgamma ;
        p_tgamma = 0x0 ;
    }
    if (p_hdirac != 0x0) {
      delete p_hdirac ;
      p_hdirac = 0x0 ;
    }

  k_dd_evol.downdate(jtime) ;
  psi_evol.downdate(jtime) ;
  hata_evol.downdate(jtime) ;
  aa_quad_evol.downdate(jtime) ;
  beta_evol.downdate(jtime) ;
  n_evol.downdate(jtime) ;
  hh_evol.downdate(jtime) ;
  

  Scalar new_expa (expansion()) ;
  //des_meridian(expa_trans, 1., 6., "Expansion trans", 1) ;
  //des_meridian(new_expa, 1.000000001, 6., "Expansion new", 2) ;
  //des_meridian(expa_trans- new_expa, 1.000000001, 4., "diff Expansion trans", 3) ;
  


}

