/*
 *  Methods of class Metconf
 *
 *   (see file metconf.h for documentation)
 *
 */

/*
 *   Copyright (c) 2002 Jerome Novak
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


char metconf_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.4  2002/08/14 15:31:17  j_novak
 * The ZEC is now correctly treated in Metconf. A test code is added in
 * Codes/Test/Metrique
 *
 * Revision 1.3  2002/08/14 13:46:15  j_novak
 * Derived quantities of a Tenseur can now depend on several Metrique's
 *
 * Revision 1.2  2002/08/13 08:02:45  j_novak
 * Handling of spherical vector/tensor components added in the classes
 * Mg3d and Tenseur. Minor corrections for the class Metconf.
 *
 * Revision 1.1  2002/08/09 15:41:09  j_novak
 * New class Metconf added for conformal metric handling.
 *
 *
 * $Header$ 
 *
 */

// Headers C++
 #include <iostream.h>

// Headers C
 #include <stdlib.h>
 #include <assert.h>
 #include <math.h>

// Headers Lorene
 #include "metconf.h"
 #include "utilitaires.h"

//Constructeur standard (ne fait pas grand chose) :
Metconf::Metconf (const Map& mapping, const Metrique& metric, const Metrique& 
		  metplat, bool jauge, bool plate) : 
  Metrique(mapping, plate), gamij(&metric), fij(&metplat), dirac(jauge) {
  assert (metplat.is_flat()) ;
  assert(metplat.get_etat() != ETATNONDEF) ;
  assert(metplat.cov().get_poids() == 0.) ;
  assert(metric.cov().get_poids() == 0.) ;
  assert(metric.get_etat() != ETATNONDEF) ;
  assert( *(metplat.cov().get_triad()) == *(metric.cov().get_triad()) ) ;
  set_der_0x0() ;
}

// COPY :
Metconf::Metconf (const Metconf& source) : Metrique(source),
					   gamij(source.gamij),
					   fij(source.fij),
					   dirac(source.dirac) {
  if (source.p_delta != 0x0)
    p_delta = new Tenseur_sym(*source.p_delta) ;
  else 
    p_delta = 0x0 ;
    
  if (source.p_Hi != 0x0)
    p_Hi = new Tenseur_sym(*source.p_Hi) ;
  else
    p_Hi = 0x0 ;
}


// Constructeur from Tensor d'ordre 2 symetrique  :
Metconf::Metconf (const Tenseur_sym& source, const Metrique& metplat, 
	         bool jauge, bool plate) : 
  Metrique(source, plate), gamij(source.get_metric()), fij(&metplat),
  dirac(jauge) {
  
  assert(metplat.is_flat()) ;
  assert(metplat.get_etat() != ETATNONDEF) ;
  assert(metplat.cov().get_poids() == 0.) ;
  assert(gamij->cov().get_poids() == 0.) ;
  assert(gamij->get_etat() != ETATNONDEF) ;
  assert( metplat.cov().get_triad()->identify() == 
	  gamij->cov().get_triad()->identify() ) ;
  int tipe = source.get_type_indice(0) ;
  
  if (tipe == CON) assert(source.get_poids() == 2./3.) ;
  else assert(source.get_poids() == -2./3.) ;
  set_der_0x0() ;
}
	
// From a file and a mapping :
Metconf::Metconf (const Map& mapping, const Base_vect& triad, const Metrique&
		metric, const Metrique& metplat, FILE* fd) :
  Metrique(mapping, triad, fd), gamij(&metric), fij(&metplat) {
    
  assert (metplat.is_flat()) ;
  assert(metplat.get_etat() != ETATNONDEF) ;
  assert(metplat.cov().get_poids() == 0.) ;
  assert(metric.cov().get_poids() == 0.) ;
  assert(metric.get_etat() != ETATNONDEF) ;
  assert( *(metplat.cov().get_triad()) == *(metric.cov().get_triad()) ) ;
  int jauge ;
  fread_be (&jauge, sizeof(int), 1, fd) ;
  dirac = jauge ;
  set_der_0x0() ;
    
}

//Destructor :
Metconf::~Metconf() {
  del_deriv() ;
}

void Metconf::del_deriv() {
  if (p_delta != 0x0)
    delete p_delta ;
  
  if (p_Hi != 0x0)
    delete p_Hi ;
    
  Metrique::del_deriv() ;
}

void Metconf::set_der_0x0() {
  Metrique::set_der_0x0() ;
  p_delta = 0x0 ;
  p_Hi = 0x0 ;
}

//AFFECTATIONS :

void Metconf::operator= (const Metconf& source) {
    
  Metrique::operator=(source) ;
  del_deriv() ;
   
  dirac = source.dirac ;
  gamij = source.gamij ;
  fij = source.fij ;
  if (source.etat == ETATZERO)
    set_etat_zero() ;
  else {
    if (source.p_delta != 0x0)
      p_delta = new Tenseur_sym (*source.p_delta) ;
    if (source.p_Hi != 0x0)
      p_Hi = new Tenseur_sym(*source.p_Hi) ;
  }
}


void Metconf::operator= (const Tenseur_sym& source) {
    
  int tipe = source.get_type_indice(0) ;
  if (tipe == COV) assert(source.get_poids() == -2./3.) ;
  else assert(source.get_poids() == 2./3.) ;
  
  Metrique::operator=(source) ;
  gamij = source.get_metric() ;
  assert(gamij->cov().get_poids() == 0.) ;
  assert(gamij->get_etat() != ETATNONDEF) ;
}

void Metconf::set_dirac_gauge(bool gauge) {
  dirac = gauge ;
  del_deriv() ;
}

void Metconf::set_flat_metric(const Metrique& metflat) {
  assert(metflat.is_flat()) ;
  assert(metflat.get_etat() != ETATNONDEF) ;
  assert( *(metflat.cov().get_triad()) == *(gamij->cov().get_triad()) ) ;
  fij = &metflat ;
  del_deriv() ;
}

void Metconf::sauve(FILE* fd) const {

  Metrique::sauve(fd) ;
  int jauge = dirac ;
  fwrite_be (&jauge, sizeof(int), 1, fd) ;

}

// Le calcul des Delta's
void Metconf::fait_delta() const {
    
  assert (etat != ETATNONDEF) ;
  if (p_delta != 0x0)
    return ;
  
  else {   // Calcul a faire :
    Itbl tipe (3) ;
    tipe.set_etat_qcq() ;
    tipe.set(0) = CON ; tipe.set(1) = COV ; tipe.set(2) = COV ; 
    p_delta = new Tenseur_sym (*mp, 3, tipe, mp->get_bvect_cart() ) ;
    bool cart = cov().get_triad()->identify() == 
      (mp->get_bvect_cart()).identify() ;
    
    if ( (etat == ETATZERO) || (plat && cart) )
      p_delta->set_etat_zero() ;
    else {
      p_delta->set_etat_qcq() ;
      assert(fij != 0x0) ;
      Tenseur t1 (contract(con(), 1, cov().derive_cov(*fij), 2)) ;
      Tenseur t2 (contract(con(), 1, cov().derive_cov(*fij), 0)) ;
      
      Cmp auxi(mp) ;
      
      // Boucle sur les composantes :
      for (int i=0 ; i<3 ; i++)
	for (int j=0 ; j<3 ; j++)
	  for (int k=j ; k<3 ; k++) {
	    auxi = 0.5*( t1(i, j, k) + t1(i, k, j) - t2(i, j, k) ) ;
	    p_delta->set(i, j, k) = auxi ;
	  }
    }
  }
}

// Le calcul des Christoffel, cas general :
void Metconf::fait_gamma() const {
    
  assert (etat != ETATNONDEF) ;
  if (p_gamma != 0x0)
    return ;
  
  else {   // Calcul a faire :
    
    Itbl tipe (3) ;
    tipe.set_etat_qcq() ;
    tipe.set(0) = CON ; tipe.set(1) = COV ; tipe.set(2) = COV ; 
    p_gamma = new Tenseur_sym (*mp, 3, tipe, mp->get_bvect_cart() ) ;
    if (etat == ETATZERO) p_gamma->set_etat_zero() ;
    else {
      assert(fij != 0x0) ;
      p_gamma->set_etat_qcq() ;
      *p_gamma = delta() + fij->gamma() ;
    }
  }
}

// Le calcul de H^i
void Metconf::fait_Hi() const {
    
  assert (etat != ETATNONDEF) ;
  if (p_Hi != 0x0)
    return ;
  
  else    // Calcul a faire :
    {
      p_Hi = new Tenseur(*mp, 1, CON, mp->get_bvect_cart(), gamij,
			 2./3.) ;
      
      if ((etat == ETATZERO) || (dirac)) p_Hi->set_etat_zero() ;
      else {
	p_Hi->set_etat_qcq() ;
	*p_Hi = - contract(contract(con(), 0, delta(), 1), 0, 2) ;
      }
    }
}

// Calcul de ricci :
void Metconf::fait_ricci() const {
    
  assert(etat != ETATNONDEF) ;
  if (p_ricci != 0x0)
    return ;
  
  else {
    p_ricci = new Tenseur_sym (*mp, 2, COV, mp->get_bvect_cart() ) ;
    if ( (etat == ETATZERO) || (plat) )	     
      p_ricci->set_etat_zero() ;
    else {
	    
      p_ricci->set_etat_qcq() ;
	    
      Tenseur_sym dcov(cov().derive_cov(*fij)) ;
      Tenseur_sym dcov2(dcov) ;
      dcov2.dec2_dzpuis() ;
      Tenseur_sym dd(dcov2.derive_cov(*fij)) ;
      dd.inc_dzpuis() ;
      Tenseur_sym dcon(con().derive_cov(*fij)) ;
      Tenseur Hi2(Hi()) ;
      Hi2.dec2_dzpuis() ;
      Tenseur dH(Hi2.derive_cov(*fij)) ;
      dH.inc_dzpuis() ;

      Tenseur_sym T1( contract(contract(delta(), 0, delta(), 2) , 1,2) ) ;
      T1.dec_dzpuis() ;
	    
      Tenseur_sym T2( contract(contract(con(), 0, dd, 0) , 0, 1) ) ;
	    
      Tenseur auxi(contract(contract(dcon, 1, dcov, 0), 1, 2)  ) ;
      auxi.dec_dzpuis() ;
      auxi = auxi + contract(cov(), 1, dH, 1) ;
      Tenseur_sym T3(*mp, 2, COV, mp->get_bvect_cart() ) ;
      T3.set_etat_qcq() ;
      // Boucle sur les composantes :
      for (int i=0 ; i<3 ; i++)
	for (int j=i ; j<3 ; j++) 
	  T3.set(i, j) = auxi(i, j) + auxi(j, i) ;
      Tenseur_sym T4 (contract(Hi(), 0, dcov, 0) );
      T4.dec_dzpuis() ;
	    
      *p_ricci = -0.5*(T2 + T3 + T4) - T1 ;
    }
  }
}

// Calcul du scalaire de ricci :
void Metconf::fait_ricci_scal() const {
    
  assert(etat != ETATNONDEF) ;
  if (p_ricci_scal != 0x0)
    return ;
  
  else {
    
    // Il s'agit d'une densite de  scalaire ...
    p_ricci_scal = new Tenseur(*mp, gamij, 2./3.) ;	    
    if ( (etat == ETATZERO) || (plat) )
      p_ricci_scal->set_etat_zero() ;
    else {
      
      p_ricci_scal->set_etat_qcq() ;
      Tenseur Hi2(Hi()) ;
      Hi2.dec2_dzpuis() ;
      Tenseur T0(contract(Hi2.derive_cov(*fij), 0, 1) ) ;
      T0.inc_dzpuis() ;
      Tenseur_sym dcov(cov().derive_cov(*fij)) ;
      Tenseur_sym dcon(con().derive_cov(*fij)) ;
      Tenseur auxi(contract(dcon, 1, dcov, 1)) ;
      Tenseur auxi1(contract(auxi, 1, 3)) ;
      Tenseur auxi2(contract(auxi, 1, 2)) ;
      Tenseur T1(contract(contract(con(), 0, auxi1, 0), 0, 1)) ;
      T1.dec_dzpuis() ;
      Tenseur T2(contract(contract(con(), 0, auxi2, 0), 0, 1)) ;
      T2.dec_dzpuis() ;

      *p_ricci_scal = -T0 + 0.25*T1 - 0.5*T2 ;
      if (dirac) p_ricci_scal->inc_dzpuis() ;
    }
  }
}


// Calcul du determinant :
void Metconf::fait_determinant() const {
    
  assert(etat != ETATNONDEF) ;
  if (p_determinant != 0x0)
    return ;
  
  else {
    
    p_determinant = new Tenseur(*mp, this, 2) ;	   
    if (etat == ETATZERO)
      p_determinant->set_etat_zero() ;
    else {
      
      p_determinant->set_etat_qcq() ;
      
      p_determinant->set() = 1. ;
    }
  }
}

const Tenseur_sym& Metconf::delta() const{
    if (p_delta == 0x0)
	fait_delta() ;
    return *p_delta ;
}

const Tenseur& Metconf::Hi() const{
    if (p_Hi == 0x0)
	fait_Hi() ;
    return *p_Hi ;
}

