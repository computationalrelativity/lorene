/*
 * Methods of class Dyn_eos.
 *
 * (see file dyneos.h for documentation).
 */

/*
 *   Copyright (c) 2019 Jerome Novak
 *             (c) 2000 Eric Gourgoulhon for Eos classes
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
 * Revision 1.3  2021/05/06 14:33:17  j_novak
 * New conversion function from Eos to Dyn_eos.
 *
 * Revision 1.2  2020/12/17 17:00:27  j_novak
 * Output of sound speed squared, instead of sound speed.
 *
 * Revision 1.1  2019/12/06 14:30:50  j_novak
 * New classes Dyn_eos... for cold Eos's with baryon density as input.
 *
 *
 *
 * $Header$
 *
 */

// Headers Lorene
#include "dyneos.h"
#include "eos.h"
#include "scalar.h"
#include "utilitaires.h"
#include "param.h"


namespace Lorene {

			//--------------//
			// Constructors //
			//--------------//


  // Standard constructor without name
  // ---------------------------------
  Dyn_eos::Dyn_eos(){}

  // Standard constructor with name
  // ---------------------------------
  Dyn_eos::Dyn_eos(const string& name_i):name(name_i){}

  // Copy constructor
  // ----------------
  Dyn_eos::Dyn_eos(const Dyn_eos& eos_i):name(eos_i.name){}

  // Constructor from a binary file
  // ------------------------------
  Dyn_eos::Dyn_eos(FILE* fich)
  {
    const int nchar = 100 ;
    char tmp_c[nchar] ;
    size_t ret = fread(tmp_c, sizeof(char), nchar, fich) ;
    if (int(ret) == nchar)
      name = tmp_c ;
  }

  // Constructor from a formatted file
  // ---------------------------------
  Dyn_eos::Dyn_eos(ifstream& fich)
  {      
    getline(fich, name, '\n') ;
  }

			//--------------//
			//  Destructor  //
			//--------------//

  Dyn_eos::~Dyn_eos(){}   // does nothing

			//-------------------------//
			//  Manipulation of name   //
			//-------------------------//
			
  void Dyn_eos::set_name(const string& name_i)
  {
    name = name_i ;  
  }

  const string& Dyn_eos::get_name() const
  {  
    return name ; 
  }

  
			//-------------------------//
			//  Conversion from Eos    //
			//-------------------------//
			
  Dyn_eos* Dyn_eos::convert_from_Eos(const Eos& eos_in) {
    int eos_type = eos_in.identify() ;
    Dyn_eos* resu ;
    switch (eos_type) {
    case 1: {
      const Eos_poly* p_poly = dynamic_cast<const Eos_poly*>(&eos_in) ;
      assert(p_poly != nullptr) ;
      resu = new Dyn_eos_poly(p_poly->get_gam(), p_poly->get_kap(),
			      p_poly->get_m_0(), p_poly->get_mu_0() ) ;
      break ;
    }
    case 10:
    case 11:
    case 12:
    case 13:
    case 14:
    case 15:
    case 16: {
      const Eos_tabul* p_tabul = dynamic_cast<const Eos_tabul*>(&eos_in) ;
      assert(p_tabul != nullptr) ;
      string eos_name = p_tabul->get_name() ;
      resu = new Dyn_eos_tab(eos_name, p_tabul->get_tablename(), false) ;
      break ;
    }
    case 17: {
      const Eos_CompOSE* p_compose
	= dynamic_cast<const Eos_CompOSE*>(&eos_in) ;
      assert(p_compose != nullptr) ;
      string eos_name = p_compose->get_name() ;
      resu = new Dyn_eos_tab(eos_name, p_compose->get_tablename(),
			     p_compose->get_format()) ;
      break ;
    }
    case 20: {
      const Eos_consistent* p_consistent
	= dynamic_cast<const Eos_consistent*>(&eos_in) ;
      assert(p_consistent != nullptr) ;
      string eos_name = p_consistent->get_name() ;
      resu = new Dyn_eos_cons(eos_name, p_consistent->get_tablename(),
			     p_consistent->get_format()) ;
      break ;
    }
    default: {
      cerr << "Dyn_eos::convert_from_Eos : " ;
      cerr << "Unknown input Eos!" << endl ;
      abort() ;
      break ;
    }
    }

    return resu ;
    
  }
  
			//------------//
			//  Outputs   //
			//------------//

  void Dyn_eos::sauve(FILE* fich) const
  {
    int ident = identify() ; 
    fwrite_be(&ident, sizeof(int), 1, fich) ;	
    fwrite(name.c_str(), sizeof(char), 100, fich) ;		
  }
    
  ostream& operator<<(ostream& ost, const Dyn_eos& eqetat)
  {
    ost << eqetat.get_name() << endl ; 
    eqetat >> ost ;
    return ost ;
  }

			//-------------------------------//
			// Generic computational routine //
			//-------------------------------//

  void Dyn_eos::calcule(const Scalar& nbar, int nzet, int l_min,  
			double (Dyn_eos::*fait)(double, const Param*) const,
			Param* par, Scalar& resu) const
  {    
    assert(nbar.get_etat() != ETATNONDEF) ;     
    const Map& mp = nbar.get_mp() ;	// Mapping
    
    if (nbar.get_etat() == ETATZERO) {
      resu.set_etat_zero() ; 
      return ; 
    }
    
    assert(nbar.get_etat() == ETATQCQ) ;
    const Valeur& vnbar = nbar.get_spectral_va() ;
    vnbar.coef_i() ;	// the values in the configuration space are required
   
    const Mg3d* mg = mp.get_mg() ;	// Multi-grid
    int nz = mg->get_nzone() ;		// total number of domains
    
    // Preparations for a point by point computation:
    resu.set_etat_qcq() ;
    Valeur& vresu = resu.set_spectral_va() ; 
    vresu.set_etat_c_qcq() ;
    vresu.c->set_etat_qcq() ;

    // Loop on domains where the computation has to be done :
    for (int l = l_min; l< l_min + nzet; l++)
      {	
	assert(l>=0) ; 
	assert(l<nz) ;
	
	Tbl* tnbar = vnbar.c->t[l] ; 
	Tbl* tresu = vresu.c->t[l] ; 
	
	if (tnbar->get_etat() == ETATZERO)
	  {
	    tresu->set_etat_zero() ; 
	  }
	else
	  {
	    assert( tnbar->get_etat() == ETATQCQ ) ; 
	    tresu->set_etat_qcq() ;
	    for (int i=0; i<tnbar->get_taille(); i++)
	      {	    
		tresu->t[i] = (this->*fait)( tnbar->t[i], par ) ;
	      }  	    
	  }  // End of the case where nbar != 0 in the considered domain  
      }  // End of the loop on domains where the computation had to be done
    
    // resu is set to zero in the other domains :
    if (l_min > 0)
      {
	resu.annule(0, l_min-1) ; 
      }
    if (l_min + nzet < nz)
      {
	resu.annule(l_min + nzet, nz - 1) ; 
      }
  }
			//---------------------------------//
			//	Public  functions	   //
			//---------------------------------//
			

  // Enthalpy from baryon density
  //------------------------------

  Scalar Dyn_eos::ent_nbar(const Scalar& nbar, int nzet, int l_min, Param* par) const
  { 
    Scalar resu(nbar.get_mp()) ; 
    
    calcule(nbar, nzet, l_min, &Dyn_eos::ent_nbar_p, par, resu) ;
    
    return resu ;    
  }


  // Energy density from baryon density 
  //------------------------------------
  
  Scalar Dyn_eos::ener_nbar(const Scalar& nbar, int nzet, int l_min, Param* par) const
  {
    Scalar resu(nbar.get_mp()) ; 
    
    calcule(nbar, nzet, l_min, &Dyn_eos::ener_nbar_p, par, resu) ;
    
    return resu ;    
  }

  // Pressure from baryon density 
  //------------------------------

  Scalar Dyn_eos::press_nbar(const Scalar& nbar, int nzet, int l_min, Param* par) const
  {    
    Scalar resu(nbar.get_mp()) ; 
    
    calcule(nbar, nzet, l_min, &Dyn_eos::press_nbar_p, par, resu) ;

    return resu ;
  }
  
  // Sound speed from baryon density
  //---------------------------------
  
  Scalar Dyn_eos::csound_square_nbar(const Scalar& nbar, int nzet,
				      int l_min, Param* par) const {
    Scalar resu(nbar.get_mp()) ;
    
    calcule(nbar, nzet, l_min, &Dyn_eos::csound_square_nbar_p, par, resu) ;
    
    return resu ;
  }

  		//--------------------------------------//
		//  Identification virtual functions	//
		//--------------------------------------//

  int Dyn_eos_poly::identify() const		{ return 1; }

  int Dyn_eos_tab::identify() const	{ return 17; }

  int Dyn_eos_cons::identify() const	{ return 20; }

		//---------------------------------------------//
		//    EOS construction from a binary file      //
		//---------------------------------------------//

  Dyn_eos* Dyn_eos::eos_from_file(FILE* fich) {
    
    Dyn_eos* p_eos ; 
    
    // Type (class) of EOS :
    int identificator ;     
    fread_be(&identificator, sizeof(int), 1, fich) ;		
    
    switch(identificator) {
      
    case 1 : {
      p_eos = new Dyn_eos_poly(fich) ; 
      break ; 
    }
      
    case 17 : {
      p_eos = new Dyn_eos_tab(fich) ;
      break ;
    }
      
    case 20 : {
      p_eos = new Dyn_eos_cons(fich) ;
      break ;
    }
      
    default : {
      cout << "Dyn_eos::eos_from_file : unknown type of EOS !" << endl ; 
      cout << " identificator = " << identificator << endl ; 
      abort() ; 
      break ; 
    } 
    }
    return p_eos ; 
  }

		//----------------------------------------------//
		//    EOS construction from a formatted file    //
		//----------------------------------------------//

Dyn_eos* Dyn_eos::eos_from_file(ifstream& fich) {
    
    int identificator ; 

    // EOS identificator : 
    fich >> identificator ; fich.ignore(1000, '\n') ;

    Dyn_eos* p_eos ; 
    
    switch(identificator) {
	
	case 1 : {
	    p_eos = new Dyn_eos_poly(fich) ; 
	    break ; 
	}
	
 	case 17 : {
	  p_eos = new Dyn_eos_tab(fich) ; 
	  break ;
	}

 	case 20 : {
	  p_eos = new Dyn_eos_cons(fich) ; 
	  break ;
	}
    default : {
      cout << "Dyn_eos::eos_from_file : unknown type of EOS !" << endl ; 
      cout << " identificator = " << identificator << endl ; 
      abort() ; 
      break ; 
    }
    }
    return p_eos ; 
}


}
