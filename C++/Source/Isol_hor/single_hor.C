/*
 *  Methods of class single_hor
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

 

/*
 * $Id$
 * $Log$
 * Revision 1.4  2016/12/05 16:17:56  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.3  2014/10/13 08:53:01  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2014/10/06 15:13:11  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.1  2007/04/13 15:28:35  f_limousin
 * Lots of improvements, generalisation to an arbitrary state of
 * rotation, implementation of the spatial metric given by Samaya.
 *
 *
 * $Header$
 *
 */

// C headers
#include <cstdlib>
#include <cassert>

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

namespace Lorene {
Single_hor::Single_hor(Map_af& mpi) : 
  mp(mpi), nz(mpi.get_mg()->get_nzone()), radius ((mpi.get_alpha())[0]), 
  omega(0),regul(0),
  n_auto(mpi), n_comp(mpi), nn(mpi), 
  psi_auto(mpi), psi_comp(mpi), psi(mpi),
  dn(mpi, COV, mpi.get_bvect_spher()), 
  dpsi(mpi, COV, mpi.get_bvect_spher()),
  beta_auto(mpi, CON, mpi.get_bvect_spher()), 
  beta_comp(mpi, CON, mpi.get_bvect_spher()),
  beta(mpi, CON, mpi.get_bvect_spher()),
  aa_auto(mpi, CON, mpi.get_bvect_spher()), 
  aa_comp(mpi, CON, mpi.get_bvect_spher()),
  aa(mpi, CON, mpi.get_bvect_spher()),
  tgam(mpi.flat_met_spher()), 
  ff(mpi.flat_met_spher()), 
  hh(mpi, CON, mpi.get_bvect_spher()),
  gamt_point(mpi, CON, mpi.get_bvect_spher()),
  trK(mpi), trK_point(mpi), decouple(mpi){
  
  hh.set_etat_zero() ;
  set_der_0x0() ;
}		  

// Copy constructor
// ----------------

Single_hor::Single_hor(const Single_hor& singlehor_in) 
    : mp(singlehor_in.mp),
      nz(singlehor_in.nz),
      radius(singlehor_in.radius),
      omega(singlehor_in.omega),
      regul(singlehor_in.regul),      
      n_auto(singlehor_in.n_auto),
      n_comp(singlehor_in.n_comp),
      nn(singlehor_in.nn),
      psi_auto(singlehor_in.psi_auto),
      psi_comp(singlehor_in.psi_comp),
      psi(singlehor_in.psi),
      dn(singlehor_in.dn),
      dpsi(singlehor_in.dpsi),
      beta_auto(singlehor_in.beta_auto),
      beta_comp(singlehor_in.beta_comp),
      beta(singlehor_in.beta),
      aa_auto(singlehor_in.aa_auto),
      aa_comp(singlehor_in.aa_comp),
      aa(singlehor_in.aa),
      tgam(singlehor_in.tgam),
      ff(singlehor_in.ff),
      hh(singlehor_in.hh),
      gamt_point(singlehor_in.gamt_point),
      trK(singlehor_in.trK),
      trK_point(singlehor_in.trK_point),
      decouple(singlehor_in.decouple){

  set_der_0x0() ;
}

// Constructor from a file
// -----------------------

Single_hor::Single_hor(Map_af& mpi, FILE* fich)
  : mp(mpi), nz(mpi.get_mg()->get_nzone()), radius ((mpi.get_alpha())[0]), 
    omega(0),regul(0), 
    n_auto(mpi, *(mpi.get_mg()), fich), n_comp(mpi), 
    nn(mpi), 
    psi_auto(mpi, *(mpi.get_mg()), fich), psi_comp(mpi), 
    psi(mpi),
    dn(mpi, COV, mpi.get_bvect_spher()), 
    dpsi(mpi, COV, mpi.get_bvect_spher()),
    beta_auto(mpi, mpi.get_bvect_spher(), fich), 
    beta_comp(mpi, CON, mpi.get_bvect_spher()),
    beta(mpi, CON, mpi.get_bvect_spher()),
    aa_auto(mpi, CON, mpi.get_bvect_spher()), 
    aa_comp(mpi, CON, mpi.get_bvect_spher()),
    aa(mpi, CON, mpi.get_bvect_spher()),
    tgam(mpi.flat_met_spher()), 
    ff(mpi.flat_met_spher()), 
    hh(mpi, CON, mpi.get_bvect_spher()),
    gamt_point(mpi, CON, mpi.get_bvect_spher()),
    trK(mpi), trK_point(mpi), decouple(mpi){
  
  fread_be(&omega, sizeof(double), 1, fich) ;

 // tgam, gamt_point, trK, trK_point

  Sym_tensor met_file (mp, mp.get_bvect_spher(), fich) ;
  tgam = met_file ;

  Sym_tensor gamt_point_file (mp, mp.get_bvect_spher(), fich) ;
  gamt_point = gamt_point_file ;
  
  Scalar trK_file (mp, *(mp.get_mg()), fich) ;
  trK = trK_file ;
  
  Scalar trK_point_file (mp, *(mp.get_mg()), fich) ;
  trK_point = trK_point_file ;
  
  set_der_0x0() ;
  hh = tgam.con() - ff.con() ;

}

			    //--------------//
			    //  Destructor  //
			    //--------------//

Single_hor::~Single_hor(){

  Single_hor::del_deriv() ;
}


                    //-----------------------//
                    // Mutators / assignment //
                    //-----------------------//

void Single_hor::operator=(const Single_hor& singlehor_in) {

 mp = singlehor_in.mp ;
 nz = singlehor_in.nz ;
 radius = singlehor_in.radius ;
 omega = singlehor_in.omega ;
 regul = singlehor_in.regul ;
 n_auto = singlehor_in.n_auto ;
 n_comp = singlehor_in.n_comp ;
 nn = singlehor_in.nn ;
 psi_auto = singlehor_in.psi_auto ;
 psi_comp = singlehor_in.psi_comp ;
 psi = singlehor_in.psi ;
 dn = singlehor_in.dn ;
 dpsi = singlehor_in.dpsi ;
 beta_auto = singlehor_in.beta_auto ;
 beta_comp = singlehor_in.beta_comp ;
 beta = singlehor_in.beta ;
 aa_auto = singlehor_in.aa_auto ;
 aa_comp = singlehor_in.aa_comp ;
 aa = singlehor_in.aa ;
 tgam = singlehor_in.tgam ;
 ff = singlehor_in.ff ;
 hh = singlehor_in.hh ;
 gamt_point = singlehor_in.gamt_point ;
 trK = singlehor_in.trK ;
 trK_point = singlehor_in.trK_point ;
 decouple = singlehor_in.decouple ;
}


                    //---------------------//
                    //  Memory management  //
                    //---------------------//

void Single_hor::del_deriv() const {

    if (p_psi4 != 0x0) delete p_psi4 ; 
    if (p_gam != 0x0) delete p_gam ; 
    if (p_k_dd != 0x0) delete p_k_dd ; 
    
    set_der_0x0() ;
}


void Single_hor::set_der_0x0() const {

    p_psi4 = 0x0 ; 
    p_gam = 0x0 ; 
    p_k_dd = 0x0 ; 
    
}

                //--------------------------//
                //      Save in a file      //
                //--------------------------//


void Single_hor::sauve(FILE* fich) const {

  n_auto.sauve(fich) ;
  psi_auto.sauve(fich) ;
  beta_auto.sauve(fich) ;
  
  fwrite_be (&omega, sizeof(double), 1, fich) ;
  
  tgam.con().sauve(fich) ;
  gamt_point.sauve(fich) ;    
  trK.sauve(fich) ;
  trK_point.sauve(fich) ;
}

// Accessors
// ---------

const Scalar& Single_hor::get_n_auto() const {

    return n_auto ;   
} 

const Scalar& Single_hor::get_n_comp() const {

    return n_comp ;   
} 
const Scalar& Single_hor::get_nn() const {

    return nn ;   
} 

const Scalar& Single_hor::get_psi_auto() const {

    return psi_auto ;   
} 

const Scalar& Single_hor::get_psi_comp() const {

    return psi_comp ;   
} 
const Scalar& Single_hor::get_psi() const {

    return psi ;   
} 
const Scalar& Single_hor::get_psi4() const {

  if (p_psi4 == 0x0)  {
    
    p_psi4 = new Scalar( pow( psi, 4.) ) ; 
    p_psi4->std_spectral_base() ;
  }
  
  return *p_psi4 ;
} 

const Vector& Single_hor::get_dn() const {

  return dn ;   
} 

const Vector& Single_hor::get_dpsi() const {

  return dpsi ;   
} 

const Vector& Single_hor::get_beta_auto() const {

  return beta_auto ;   
} 

const Vector& Single_hor::get_beta_comp() const {

  return beta_comp ;   

} 
const Vector& Single_hor::get_beta() const {

  return beta ;   
} 

const Sym_tensor& Single_hor::get_aa_auto() const {

  return aa_auto ;   
} 

const Sym_tensor& Single_hor::get_aa_comp() const {

  return aa_comp ;   

} 
const Sym_tensor& Single_hor::get_aa() const {

  return aa ;   

} 
const Metric& Single_hor::get_gam() const {

  if (p_gam == 0x0) {
    p_gam = new Metric( get_psi4()*tgam.cov() ) ; 
  }
  
  return *p_gam ; 
} 

const Sym_tensor& Single_hor::get_k_dd() const {

  if (p_k_dd == 0x0)  {
    
    Sym_tensor temp (aa/get_psi4()+1./3.*trK*get_gam().con()) ;
    
    p_k_dd = new Sym_tensor( temp.up_down(get_gam()) ) ; 
        p_k_dd->std_spectral_base() ;
  }
  
  return *p_k_dd ;
} 


void Single_hor::set_psi_auto(const Scalar& psi_in) {

  psi_auto = psi_in ;    
}
void Single_hor::set_n_auto(const Scalar& n_in) {

  n_auto = n_in ;    
}
void Single_hor::set_beta_auto(const Scalar& beta_in) {

  beta_auto = beta_in ;    
}
void Single_hor::set_aa_auto(const Scalar& aa_in) {

  aa_auto = aa_in ;    
}
void Single_hor::set_aa_comp(const Scalar& aa_in) {

  aa_comp = aa_in ;    
}
void Single_hor::set_aa(const Scalar& aa_in) {

  aa = aa_in ;    
}


// Import the lapse from the companion (Bhole case)

void Single_hor::n_comp_import(const Single_hor& comp) {

    Scalar temp (mp) ;
    temp.import(comp.n_auto) ;
    temp.std_spectral_base() ;
    n_comp = temp ;
    nn = temp + n_auto ;
     
    Vector dn_comp (mp, COV, mp.get_bvect_cart()) ;
    dn_comp.set_etat_qcq() ;
    Vector auxi (comp.n_auto.derive_cov(comp.ff)) ;
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

    dn = n_auto.derive_cov(ff) + dn_comp ;
}

// Import the conformal factor from the companion (Bhole case)

void Single_hor::psi_comp_import(const Single_hor& comp) {

    Scalar temp (mp) ;
    temp.import(comp.psi_auto) ;
    temp.std_spectral_base() ;
    psi_comp = temp ;
    psi = temp + psi_auto ;
    
    Vector dpsi_comp (mp, COV, mp.get_bvect_cart()) ;
    dpsi_comp.set_etat_qcq() ;
    Vector auxi (comp.psi_auto.derive_cov(comp.ff)) ;
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
    
    dpsi = psi_auto.derive_cov(ff) + dpsi_comp ;

}

void Single_hor::beta_comp_import(const Single_hor& comp) {

    Vector tmp_vect (mp, CON, mp.get_bvect_cart()) ;
    Vector shift_comp (comp.beta_auto) ;
    shift_comp.change_triad(comp.mp.get_bvect_cart()) ;
    shift_comp.change_triad(mp.get_bvect_cart()) ;
    assert (*(shift_comp.get_triad()) == *(tmp_vect.get_triad())) ;

    tmp_vect.set(1).import(shift_comp(1)) ;
    tmp_vect.set(2).import(shift_comp(2)) ;
    tmp_vect.set(3).import(shift_comp(3)) ;
    tmp_vect.std_spectral_base() ;
    tmp_vect.change_triad(mp.get_bvect_spher()) ;
    
    beta_comp = tmp_vect ;
    beta = beta_auto + beta_comp ;
}

//Initialisation to Schwartzchild
void Single_hor::init_bhole () {
    
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
    n_auto = temp ;

    temp = 0.5 ;
    temp.std_spectral_base() ;
    n_comp = temp ;
    nn = n_auto  ;
  
    auxi = 0.5 + radius/mp.r ;
    auxi.annule(0, 0);
    auxi.set_dzpuis(0) ;
    temp = auxi;
    temp.std_spectral_base() ;
    temp.raccord(1) ;
    psi_auto = temp ;

    temp = 0.5 ;
    temp.std_spectral_base() ;
    psi_comp = temp ;
    psi = psi_auto + psi_comp ;

    dn = nn.derive_cov(ff) ;
    dpsi = psi.derive_cov(ff) ;
    
    Vector temp_vect1(mp, CON, mp.get_bvect_spher()) ;
    temp_vect1.set(1) = 0.0/mp.r/mp.r ;
    temp_vect1.set(2) = 0. ;
    temp_vect1.set(3) = 0. ;
    temp_vect1.std_spectral_base() ;

    Vector temp_vect2(mp, CON, mp.get_bvect_spher()) ;
    temp_vect2.set_etat_zero() ;    

    beta_auto = temp_vect1 ;
    beta_comp = temp_vect2 ;
    beta = temp_vect1 ;    

}

void Single_hor::init_met_trK() {
 
  Metric flat (mp.flat_met_spher()) ;
  tgam = flat ;

  gamt_point.set_etat_zero() ;
  trK.set_etat_zero() ;
  trK_point.set_etat_zero() ;
 
}


void Single_hor::init_bhole_seul () {
    
    Scalar auxi(mp) ;
    
    auxi = (1-radius/mp.r)/(1+radius/mp.r) ;
    auxi.annule(0, 0);
    auxi.set_outer_boundary((*mp.get_mg()).get_nzone(), 1.) ;
    auxi.set_dzpuis(0) ;

    Scalar temp(mp) ;
    temp = auxi;
    temp.std_spectral_base() ;
    temp.raccord(1) ;
    n_auto = temp;

    temp.set_etat_zero() ;
    n_comp = temp ;
    nn = temp ; 
    
    auxi = 1 + radius/mp.r ;
    auxi.annule(0, 0);
    auxi.set_outer_boundary((*mp.get_mg()).get_nzone(), 1.) ;
    auxi.set_dzpuis(0) ;
  
    temp = auxi;
    temp.std_spectral_base() ;
    temp.raccord(1) ;
    psi_auto = temp ;
    temp.set_etat_zero() ;
    psi_comp = temp ;
    psi = temp ;
    
    dn = nn.derive_cov(ff) ;
    dpsi = psi.derive_cov(ff) ;

    Vector temp_vect(mp, CON, mp.get_bvect_spher()) ;
    temp_vect.set_etat_zero() ;
    beta_auto = temp_vect ;
    beta_comp = temp_vect ;
    beta = temp_vect ;    	

}		   



}
