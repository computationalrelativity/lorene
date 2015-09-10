/*
 *  Methods of the class Hot_eos.
 *
 *    (see file hoteos.h for documentation).
 *
 */

/*
 *   Copyright (c) 2015  Jerome Novak
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

char hoteos_C[] = "$Header$" ;

/*
 * $Id$
 * $Log $
 *
 * $Header $
 *
 */

// C headers
#include <cassert>

// Lorene headers
#include "hoteos.h"
#include "eos.h"
#include "tensor.h"
#include "utilitaires.h"

namespace Lorene {

			//--------------//
			// Constructors //
			//--------------//


// Standard constructor without name
// ---------------------------------
  Hot_eos::Hot_eos(){}
  
  // Standard constructor with name
  // ---------------------------------
  Hot_eos::Hot_eos(const string& name_i):name(name_i){ 
    set_der_0x0() ;
  }

  Hot_eos::Hot_eos(const char* name_i):name(name_i){ 
    set_der_0x0() ;
  }
  
  // Copy constructor
  // ----------------
  Hot_eos::Hot_eos(const Hot_eos& eos_i):name(eos_i.name){
    set_der_0x0() ;
  }
  
  // Constructor from a binary file
  // ------------------------------
  Hot_eos::Hot_eos(FILE* fich){
    
    int taille ;
    fread(&taille, sizeof(int), 1, fich) ;		
    assert(taille > 0) ;
    char* t_name = new char[taille] ;
    fread(t_name, sizeof(char), taille, fich) ;		
    set_name(t_name) ;
    delete [] t_name ;
    set_der_0x0() ;
  }
  
  // Constructor from a formatted file
  // ---------------------------------
  Hot_eos::Hot_eos(ifstream& fich){
    
    char t_name[100] ;
    fich.getline(t_name, 100) ;
    set_name(t_name) ;
    set_der_0x0() ;
  }

			//--------------//
			//  Destructor  //
			//--------------//
  Hot_eos::~Hot_eos(){
    Hot_eos::del_deriv() ;
  }

               //----------------------------------//
               // Management of derived quantities //
               //----------------------------------//

  void Hot_eos::del_deriv() const {
    if (p_cold_eos != 0x0)  delete p_cold_eos ;
    set_der_0x0() ;
  }

  void Hot_eos::set_der_0x0() const {
    p_cold_eos = 0x0 ;
  }

void Hot_eos::set_name(const char* name_i) {

  name.assign(name_i) ; 
    
}


			//------------//
			//  Outputs   //
			//------------//

void Hot_eos::sauve(FILE* fich) const {

    int ident = identify() ; 
    fwrite_be(&ident, sizeof(int), 1, fich) ;

    int taille = int(name.size()) ;
    fwrite_be(&taille, sizeof(int), 1, fich) ;
    fwrite(name.c_str(), sizeof(char), name.size(), fich) ;		
   
}
    



ostream& operator<<(ostream& ost, const Hot_eos& eqetat)  {
    ost << eqetat.get_name() << endl ; 
    eqetat >> ost ;
    return ost ;
}

			//-------------------------------//
			// Generic computational routine //
			//-------------------------------//

  void Hot_eos::calcule(const Scalar& ent, const Scalar& sb, int nzet, int l_min,  
			double (Hot_eos::*fait)(double, double) const,
			Scalar& resu) const {
    
    assert(ent.get_etat() != ETATNONDEF) ; 
    assert(sb.get_etat() != ETATNONDEF) ; 
    
    const Map* mp = &(ent.get_mp()) ;	// Mapping
    
    const Mg3d* mg = mp->get_mg() ;	// Multi-grid
    
    int nz = mg->get_nzone() ;		// total number of domains
        
    if (ent.get_etat() == ETATZERO) {
      resu.set_etat_zero() ; 
      return ; 
    }

    assert(ent.get_etat() == ETATQCQ) ; 
    const Valeur& vent = ent.get_spectral_va() ;
    vent.coef_i() ;	// the values in the configuration space are required
   
    const Valeur* vsb = &sb.get_spectral_va() ;
    Valeur vzero(mg) ;
    if (sb.get_etat() == ETATZERO) {
      vzero.annule_hard() ;
      vsb  = &vzero ;
    }

    assert(vsb->get_mg() == vent.get_mg()) ;

    // Preparations for a point by point computation:
    resu.set_etat_qcq() ;
    Valeur& vresu = resu.set_spectral_va() ; 
    vresu.set_etat_c_qcq() ;
    vresu.c->set_etat_qcq() ;

    // Loop on domains where the computation has to be done :
    for (int l = l_min; l< l_min + nzet; l++) {
      
      assert(l>=0) ; 
      assert(l<nz) ;
      
      bool tsb0 = false ;
      Tbl* tent = vent.c->t[l] ; 
      Tbl* tsb = vsb->c->t[l] ; 
      Tbl* tresu = vresu.c->t[l] ; 
	
      if (tent->get_etat() == ETATZERO) {
	tresu->set_etat_zero() ; 
      }
      else {
	assert( tent->get_etat() == ETATQCQ ) ; 
	tresu->set_etat_qcq() ;
	if (tsb->get_etat() == ETATZERO) {
	  tsb0 = true ;
	  tsb = new Tbl(tent->dim) ;
	  tsb->annule_hard() ;
	}
	  
	for (int i=0; i<tent->get_taille(); i++) {
	  
	  tresu->t[i] = (this->*fait)( tent->t[i], tsb->t[i] ) ;
	}  
	    
      }  // End of the case where ent != 0 in the considered domain  
      if (tsb0) delete tsb ;
    }  // End of the loop on domains where the computation had to be done

    // resu is set to zero in the other domains :
    
    if (l_min > 0) {
      resu.annule(0, l_min-1) ; 
    }
    
    if (l_min + nzet < nz) {
      resu.annule(l_min + nzet, nz - 1) ; 
    }
  }

  // Baryon density from enthalpy 
  //------------------------------
  
  Scalar Hot_eos::nbar_Hs(const Scalar& ent, const Scalar& sb, int nzet, int l_min)
    const {
    
    Scalar resu(ent.get_mp()) ; 
    
    calcule(ent, sb, nzet, l_min, &Hot_eos::nbar_Hs_p, resu) ;
    
    return resu ; 
    
  }
  
  
  
  // Energy density from enthalpy 
  //------------------------------
  
  Scalar Hot_eos::ener_Hs(const Scalar& ent, const Scalar& sb, int nzet, int l_min) 
    const {
    
    Scalar resu(ent.get_mp()) ; 
    
    calcule(ent, sb, nzet, l_min, &Hot_eos::ener_Hs_p, resu) ;
    
    return resu ; 
    
  }

  // Pressure from enthalpy 
  //-----------------------
  
  Scalar Hot_eos::press_Hs(const Scalar& ent, const Scalar& sb, int nzet, int l_min) 
    const {
    
    Scalar resu(ent.get_mp()) ; 
    
    calcule(ent, sb, nzet, l_min, &Hot_eos::press_Hs_p, resu) ;
    
    return resu ;
    
  }
  
  // Temperature from enthalpy 
  //--------------------------
  
  Scalar Hot_eos::temp_Hs(const Scalar& ent, const Scalar& sb, int nzet, int l_min) 
    const {
    
    Scalar resu(ent.get_mp()) ; 
    
    calcule(ent, sb, nzet, l_min, &Hot_eos::temp_Hs_p, resu) ;
    
    return resu ;
    
  }
  
}
