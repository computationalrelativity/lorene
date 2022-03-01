/*
 * Methods for Hot_eos and file manipulation
 *
 * (see file hoteos.h for documentation)
 */

/*
 *   Copyright (c) 2015 Jerome Novak
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
 * Revision 1.4  2022/03/01 10:03:06  g_servignat
 * Corrected all pure virtual interference between Hot_eos derived classes ; amended use of interpol_linear_2D
 *
 * Revision 1.3  2016/12/05 16:17:52  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.2  2015/12/08 10:52:18  j_novak
 * New class Hoteos_tabul for tabulated temperature-dependent EoSs.
 *
 * Revision 1.1  2015/03/17 14:20:00  j_novak
 * New class Hot_eos to deal with temperature-dependent EOSs.
 *
 * $Header$
 *
 */
 
// Headers C
#include <cstdlib>

// Header Lorene
#include "headcpp.h"
#include "hoteos.h"
#include "utilitaires.h"

namespace Lorene {

		//--------------------------------------//
		//  Identification virtual functions	//
		//--------------------------------------//

  int Ideal_gas::identify() const		{ return 1; }
  
  int Hoteos_tabul::identify() const	        { return 2; }
  
  int Ye_eos_tabul::identify() const {return 3;}


		//-------------------------------------------------//
		//    Hot EOS construction from a binary file      //
		//-------------------------------------------------//

  Hot_eos* Hot_eos::hoteos_from_file(FILE* fich) {
    
    Hot_eos* p_eos ; 
    
    // Type (class) of EOS :
    int identificator ;     
    fread_be(&identificator, sizeof(int), 1, fich) ;		
    
    switch(identificator) {

    case 1 : {
      p_eos = new Ideal_gas(fich) ; 
      break ; 
    }
	
    case 2 : {
      p_eos = new Hoteos_tabul(fich) ; 
      break ; 
    }
	
    default : {
      cout << "Hot_eos::hoteos_from_file : unknown type of EOS !" << endl ; 
      cout << " identificator = " << identificator << endl ; 
      abort() ; 
      break ; 
    }
      
    }
    
    return p_eos ; 
    
  }

		//--------------------------------------------------//
		//    Hot EOS construction from a formatted file    //
		//--------------------------------------------------//

  Hot_eos* Hot_eos::hoteos_from_file(ifstream& fich) {
    
    int identificator ; 
    if (!fich) {
      cerr << "Hot_eos::hoteos_from_file: file cannot be opened!" << endl ;
      abort() ;
    }
    
    // EOS identificator : 
    fich >> identificator ; fich.ignore(1000, '\n') ;
    
    Hot_eos* p_eos ; 
    
    switch(identificator) {

    case 1 : {
      p_eos = new Ideal_gas(fich) ; 
      break ; 
    }
	
    case 2 : {
      p_eos = new Hoteos_tabul(fich) ; 
      break ; 
    }
	
    default : {
      cout << "Hot_eos::hoteos_from_file : unknown type of EOS !" << endl ; 
      cout << " identificator = " << identificator << endl ; 
      abort() ; 
      break ; 
    }
      
    }
    
    return p_eos ; 
    
  }


}
