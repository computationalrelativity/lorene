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

#include "utilitaires.h"
#include "time_slice.h"
#include "isol_hor.h"
#include "tensor.h"
#include "metric.h"
#include "evolution.h"



			    //--------------//
			    // Constructors //
			    //--------------//
// Standard constructor
// --------------------

Isol_hor::Isol_hor(Map_af& mpi, int depth_in) : 
  Time_slice_conf(mpi, mpi.get_bvect_spher(), mpi.flat_met_spher()),
  mp(mpi), radius ((mpi.get_alpha())[0]), omega(0), boost_x(0),
  boost_z(0), regul(0),
  n_auto_evol(depth_in), n_comp_evol(depth_in), 
  psi_auto_evol(depth_in), psi_comp_evol(depth_in),
  dn_evol(depth_in), dpsi_evol(depth_in),
  beta_auto_evol(depth_in), beta_comp_evol(depth_in),
  aa_auto_evol(depth_in), aa_comp_evol(depth_in), aa_nn(depth_in),
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
		      ff_in.con(), aa_in, trK_in, depth_in),
      mp(mpi), radius ((mpi.get_alpha())[0]), 
      omega(0), boost_x(0), boost_z(0), regul(0),
      n_auto_evol(depth_in), n_comp_evol(depth_in), 
      psi_auto_evol(depth_in), psi_comp_evol(depth_in),
      dn_evol(depth_in), dpsi_evol(depth_in),
      beta_auto_evol(depth_in), beta_comp_evol(depth_in),
      aa_auto_evol(depth_in), aa_comp_evol(depth_in), aa_nn(depth_in),
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
      mp(mpi), radius ((mpi.get_alpha())[0]), omega(0), boost_x(0),
      boost_z(0), regul(0),
      n_auto_evol(depth_in), n_comp_evol(depth_in), 
      psi_auto_evol(depth_in), psi_comp_evol(depth_in),
      dn_evol(depth_in), dpsi_evol(depth_in),
      beta_auto_evol(depth_in), beta_comp_evol(depth_in),
      aa_auto_evol(depth_in), aa_comp_evol(depth_in), aa_nn(depth_in),
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

// Import the lapse from the companion (Bhole case)

void Isol_hor::n_comp(const Isol_hor& comp) {

    double ttime = the_time[jtime] ;    

    Scalar temp (mp) ;
    temp.import_symy(comp.n_auto()) ;
    temp.std_spectral_base() ;
    n_comp_evol.update(temp, jtime, ttime) ;
    n_evol.update(temp + n_auto(), jtime, ttime) ;
     
    Vector dn_comp (mp, COV, mp.get_bvect_cart()) ;
    dn_comp.set_etat_qcq() ;
    Vector auxi (comp.n_auto().derive_cov(comp.ff)) ;
    auxi.dec_dzpuis(2) ;
    auxi.change_triad(auxi.get_mp().get_bvect_cart()) ;
    auxi.change_triad(mp.get_bvect_cart()) ;
    assert ( *(auxi.get_triad()) == *(dn_comp.get_triad())) ;

    dn_comp.set(1).import_symy(auxi(1)) ;
    dn_comp.set(2).import_asymy(auxi(2)) ;
    dn_comp.set(3).import_symy(auxi(3)) ;
    dn_comp.std_spectral_base() ;
    dn_comp.inc_dzpuis(2) ;
    dn_comp.change_triad(mp.get_bvect_spher()) ;

    dn_evol.update(n_auto().derive_cov(ff) + dn_comp, jtime, ttime) ;
}

// Import the conformal factor from the companion (Bhole case)

void Isol_hor::psi_comp (const Isol_hor& comp) {
  
    double ttime = the_time[jtime] ;    
    
    Scalar temp (mp) ;
    temp.import_symy(comp.psi_auto()) ;
    temp.std_spectral_base() ;
    psi_comp_evol.update(temp, jtime, ttime) ;
    psi_evol.update(temp + psi_auto(), jtime, ttime) ;
    
    Vector dpsi_comp (mp, COV, mp.get_bvect_cart()) ;
    dpsi_comp.set_etat_qcq() ;
    Vector auxi (comp.psi_auto().derive_cov(comp.ff)) ;
    auxi.dec_dzpuis(2) ;
    auxi.change_triad(auxi.get_mp().get_bvect_cart()) ;
    auxi.change_triad(mp.get_bvect_cart()) ;
    assert ( *(auxi.get_triad()) == *(dpsi_comp.get_triad())) ;

    dpsi_comp.set(1).import_symy(auxi(1)) ;
    dpsi_comp.set(2).import_asymy(auxi(2)) ;
    dpsi_comp.set(3).import_symy(auxi(3)) ;
    dpsi_comp.std_spectral_base() ;
    dpsi_comp.inc_dzpuis(2) ;
    dpsi_comp.change_triad(mp.get_bvect_spher()) ;
    
    dpsi_evol.update(psi_auto().derive_cov(ff) + dpsi_comp, jtime, ttime) ;
}

void Isol_hor::beta_comp (const Isol_hor& comp) {
  
    double ttime = the_time[jtime] ;    
    
    Vector tmp_vect (mp, CON, mp.get_bvect_cart()) ;
    Vector shift_comp (comp.beta_auto()) ;
    shift_comp.change_triad(comp.mp.get_bvect_cart()) ;
    shift_comp.change_triad(mp.get_bvect_cart()) ;
    assert (*(shift_comp.get_triad()) == *(tmp_vect.get_triad())) ;

    tmp_vect.set(1).import_asymy(shift_comp(1)) ;
    tmp_vect.set(2).import_symy(shift_comp(2)) ;
    tmp_vect.set(3).import_asymy(shift_comp(3)) ;
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
    auxi = 0.5 - radius/mp.r ;
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
    
    Vector temp_vect(mp, CON, mp.get_bvect_spher()) ;
    temp_vect.set_etat_zero() ;
    beta_auto_evol.update(temp_vect, jtime, ttime) ;
    beta_comp_evol.update(temp_vect, jtime, ttime) ;
    beta_evol.update(temp_vect, jtime, ttime) ;    
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
  
  set_aa(aa_new) ;

  return ;
  
}

void Isol_hor::kerr_perturb() {

    Sym_tensor gam (gam().cov()) ;
    Sym_tensor gamt (gam / gam(3,3)) ;

    Metric metgamt (gamt) ;
    met_gamt = metgamt ;

    cout << "met_gamt" << endl << norme(met_gamt.cov()(1,1)) << endl << norme(met_gamt.cov()(2,1)) << endl << norme(met_gamt.cov()(3,1)) << endl << norme(met_gamt.cov()(2,2)) << endl << norme(met_gamt.cov()(3,2)) << endl << norme(met_gamt.cov()(3,3)) << endl ;
    cout << "determinant" << norme(met_gamt.determinant()) << endl ;


    hh_evol.update(met_gamt.con() - ff.con(), jtime, the_time[jtime]) ;

}
