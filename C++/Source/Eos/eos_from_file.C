/*
 * Methods for Eos and file manipulation
 *
 * (see file eos.h for documentation)
 */

/*
 *   Copyright (c) 2000-2001 Eric Gourgoulhon
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
 * Revision 1.19  2023/06/15 08:04:05  g_servignat
 * Added identifier for Pseudo_polytrope_1D and Piecewise_polytrope_1D
 *
 * Revision 1.18  2022/04/15 13:39:24  j_novak
 * New class Eos_compose_fit to generate fitted EoSs from CompOSE tables.
 *
 * Revision 1.17  2016/12/05 16:17:51  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.16  2015/08/04 14:41:29  j_novak
 * Back to previous version for Eos_CompOSE. Enthalpy-consistent EoS can be accessed using Eos_consistent class (derived from Eos_CompOSE).
 *
 * Revision 1.15  2015/01/27 14:22:38  j_novak
 * New methods in Eos_tabul to correct for EoS themro consistency (optional).
 *
 * Revision 1.14  2014/10/13 08:52:53  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.13  2014/06/30 16:13:18  j_novak
 * New methods for reading directly from CompOSE files.
 *
 * Revision 1.12  2014/03/06 15:53:35  j_novak
 * Eos_compstar is now Eos_compOSE. Eos_tabul uses strings and contains informations about authors.
 *
 * Revision 1.11  2012/10/26 14:09:33  e_gourgoulhon
 * Added new class Eos_Fermi
 *
 * Revision 1.10  2011/06/16 10:49:18  j_novak
 * New class Eos_mag for EOSs depending on density and magnetic field.
 *
 * Revision 1.9  2010/02/02 13:22:16  j_novak
 * New class Eos_Compstar.
 *
 * Revision 1.8  2005/05/22 20:51:41  k_taniguchi
 * Add a new Eos Eos_fit_AkmalPR.
 *
 * Revision 1.7  2004/09/26 18:53:08  k_taniguchi
 * Introduction of new EOSs: Eos_fit_SLy4 and Eos_fit_FPS
 *
 * Revision 1.6  2004/05/07 08:06:45  k_taniguchi
 * Add the case of Eos_multi_poly.C
 *
 * Revision 1.5  2003/12/08 15:47:03  m_bejger
 * GlendNH3 EOS (Glendenning 1985, case 3) added
 *
 * Revision 1.4  2002/10/16 14:36:35  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.3  2002/04/09 14:32:15  e_gourgoulhon
 * 1/ Added extra parameters in EOS computational functions (argument par)
 * 2/ New class MEos for multi-domain EOS
 *
 * Revision 1.2  2001/12/04 21:27:53  e_gourgoulhon
 *
 * All writing/reading to a binary file are now performed according to
 * the big endian convention, whatever the system is big endian or
 * small endian, thanks to the functions fwrite_be and fread_be
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.6  2001/09/11  16:23:08  eric
 * Ajout de Eos_AkmalPR, Eos_BBB2 et Eos_BalbN1H1.
 *
 * Revision 2.5  2000/11/23  22:34:10  eric
 * Ajout de Eos_BPAL12.
 *
 * Revision 2.4  2000/11/23  14:46:16  eric
 * Ajout de Eos_strange_cr.
 *
 * Revision 2.3  2000/11/22  19:30:55  eric
 * Ajout des Eos_SLy4 et Eos_FPS
 *
 * Revision 2.2  2000/10/24  15:29:22  eric
 * Ajout de l'EOS matiere etrange (Eos_strange).
 *
 * Revision 2.1  2000/02/14  14:33:41  eric
 * Ajout du constructeur par lecture de fichier formate.
 *
 * Revision 2.0  2000/01/21  15:18:08  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */
 
// Headers C
#include <cstdlib>

// Header Lorene
#include "headcpp.h"
#include "eos.h"
#include "eos_multi_poly.h"
#include "eos_compose_fit.h"
#include "eos_fitting.h"
#include "utilitaires.h"

		//--------------------------------------//
		//  Identification virtual functions	//
		//--------------------------------------//


namespace Lorene {
int Eos_poly::identify() const		{ return 1; }

int Eos_poly_newt::identify() const	{ return 2; }

int Eos_incomp::identify() const	{ return 3; }

int Eos_incomp_newt::identify() const	{ return 4; }

int Eos_strange::identify() const	{ return 5; }

int Eos_strange_cr::identify() const	{ return 6; }

int Eos_compose_fit::identify() const	{ return 7; }

int Eos_SLy4::identify() const		{ return 10; }

int Eos_FPS::identify() const		{ return 11; }

int Eos_BPAL12::identify() const	{ return 12; }

int Eos_AkmalPR::identify() const	{ return 13; }

int Eos_BBB2::identify() const		{ return 14; }

int Eos_BalbN1H1::identify() const	{ return 15; }

int Eos_GlendNH3::identify() const	{ return 16; }

int Eos_CompOSE::identify() const	{ return 17; }

int Eos_mag::identify() const	{ return 18; }

int Eos_Fermi::identify() const	{ return 19; }

int Eos_consistent::identify() const	{ return 20; }

int Pseudo_polytrope_1D::identify() const	{ return 21; }

int Piecewise_polytrope_1D::identify() const	{ return 22; }

int MEos::identify() const	{ return 100; }

int Eos_multi_poly::identify() const	{ return 110; }

int Eos_fit_SLy4::identify() const      { return 120; }

int Eos_fit_FPS::identify() const       { return 121; }

int Eos_fit_AkmalPR::identify() const   { return 122; }

		//---------------------------------------------//
		//    EOS construction from a binary file      //
		//---------------------------------------------//

Eos* Eos::eos_from_file(FILE* fich) {
    
    Eos* p_eos ; 
    
    // Type (class) of EOS :
    int identificator ;     
    fread_be(&identificator, sizeof(int), 1, fich) ;		

    switch(identificator) {
	
	case 1 : {
	    p_eos = new Eos_poly(fich) ; 
	    break ; 
	}
	
	case 2 : {
	    p_eos = new Eos_poly_newt(fich) ; 
	    break ; 
	}
	
	case 3 : {
	    p_eos = new Eos_incomp(fich) ; 
	    break ; 
	}
	
	case 4 : {
	    p_eos = new Eos_incomp_newt(fich) ; 
	    break ; 
	}
	
	case 5 : {
	    p_eos = new Eos_strange(fich) ;
	    break ;
	}
	
	case 6 : {
	    p_eos = new Eos_strange_cr(fich) ;
	    break ;
	}
	
	case 7 : {
	    p_eos = new Eos_compose_fit(fich) ;
	    break ;
	}
	
	case 10 : {
	    p_eos = new Eos_SLy4(fich) ;
	    break ;
	}
	
	case 11 : {
	    p_eos = new Eos_FPS(fich) ;
	    break ;
	}
	
	case 12 : {
	    p_eos = new Eos_BPAL12(fich) ;
	    break ;
	}
	
	case 13 : {
	    p_eos = new Eos_AkmalPR(fich) ;
	    break ;
	}
	
	case 14 : {
	    p_eos = new Eos_BBB2(fich) ;
	    break ;
	}
	
	case 15 : {
	    p_eos = new Eos_BalbN1H1(fich) ;
	    break ;
	}

	case 16 : {
	    p_eos = new Eos_GlendNH3(fich) ;
	    break ;
	}

	case 17 : {
	    p_eos = new Eos_CompOSE(fich) ;
	    break ;
	}

	case 18 : {
	    p_eos = new Eos_mag(fich) ;
	    break ;
	}

	case 19 : {
	    p_eos = new Eos_Fermi(fich) ;
	    break ;
	}

	case 20 : {
	    p_eos = new Eos_consistent(fich) ;
	    break ;
	}

	case 21 : {
	    p_eos = new Pseudo_polytrope_1D(fich) ;
	    break ;
	}

	case 22 : {
	    p_eos = new Piecewise_polytrope_1D(fich) ;
	    break ;
	}

	case 100 : {
	    p_eos = new MEos(fich) ;
	    break ;
	}

	case 110 : {
	    p_eos = new Eos_multi_poly(fich) ;
	    break ;
	}

	case 120 : {
	    p_eos = new Eos_fit_SLy4(fich) ;
	    break ;
	}

	case 121 : {
	    p_eos = new Eos_fit_FPS(fich) ;
	    break ;
	}

	case 122 : {
	    p_eos = new Eos_fit_AkmalPR(fich) ;
	    break ;
	}

	default : {
	    cout << "Eos::eos_from_file : unknown type of EOS !" << endl ; 
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

Eos* Eos::eos_from_file(ifstream& fich) {
    
    int identificator ; 

    // EOS identificator : 
    fich >> identificator ; fich.ignore(1000, '\n') ;

    Eos* p_eos ; 
    
    switch(identificator) {
	
	case 1 : {
	    p_eos = new Eos_poly(fich) ; 
	    break ; 
	}
	
	case 2 : {
	    p_eos = new Eos_poly_newt(fich) ; 
	    break ; 
	}
	
	case 3 : {
	    p_eos = new Eos_incomp(fich) ; 
	    break ; 
	}
	
	case 4 : {
	    p_eos = new Eos_incomp_newt(fich) ; 
	    break ; 
	}
	
	case 5 : {
	    p_eos = new Eos_strange(fich) ;
	    break ;
	}
	
	case 6 : {
	    p_eos = new Eos_strange_cr(fich) ;
	    break ;
	}
	
	case 7 : {
	    p_eos = new Eos_compose_fit(fich) ;
	    break ;
	}
	
	case 10 : {
	    p_eos = new Eos_SLy4(fich) ;
	    break ;
	}
	
	case 11 : {
	    p_eos = new Eos_FPS(fich) ;
	    break ;
	}
	
	case 12 : {
	    p_eos = new Eos_BPAL12(fich) ;
	    break ;
	}
	
	case 13 : {
	    p_eos = new Eos_AkmalPR(fich) ;
	    break ;
	}
	
	case 14 : {
	    p_eos = new Eos_BBB2(fich) ;
	    break ;
	}
	
	case 15 : {
	    p_eos = new Eos_BalbN1H1(fich) ;
	    break ;
	}

	case 16 : {
	    p_eos = new Eos_GlendNH3(fich) ;
	    break ;
	}

	case 17 : {
	  int format ;
	  fich >> format ;
	  fich.ignore(1000, '\n') ;
#ifndef NDEBUG
	  cout << "Reading tabulated EoS, with "
	       << ( (format == 0) ? "standard LORENE " : "original CompOSE ")
	       << "format." << endl ;
#endif
	  if (format == 1) {
	    fich.ignore(1000, '\n') ;
	    string files_path ;
	    fich >> files_path ;	    
	    p_eos = new Eos_CompOSE(files_path ) ;
	  }
	  else 
	    p_eos = new Eos_CompOSE(fich ) ;
	  break ;
	}

	case 18 : {
	    p_eos = new Eos_mag(fich) ;
	    break ;
	}

	case 19 : {
	    p_eos = new Eos_Fermi(fich) ;
	    break ;
	}

	case 20 : {
	  int format ;
	  fich >> format ;
	  fich.ignore(1000, '\n') ;
#ifndef NDEBUG
	  cout << "Reading tabulated EoS, with "
	       << ( (format == 0) ? "standard LORENE " : "original CompOSE ")
	       << "format." << endl ;
#endif
	  if (format == 1) {
	    fich.ignore(1000, '\n') ;
	    string files_path ;
	    fich >> files_path ;	    
	    p_eos = new Eos_consistent(files_path ) ;
	  }
	  else 
	    p_eos = new Eos_consistent(fich ) ;
	  break ;
	}

	case 21 : {
	    p_eos = new Pseudo_polytrope_1D(fich) ;
	    break ;
	}

	case 22 : {
	    p_eos = new Piecewise_polytrope_1D(fich) ;
	    break ;
	}

	case 100 : {
	    p_eos = new MEos(fich) ;
	    break ;
	}

	case 110 : {
	    p_eos = new Eos_multi_poly(fich) ;
	    break ;
	}

	case 120 : {
	    p_eos = new Eos_fit_SLy4(fich) ;
	    break ;
	}

	case 121 : {
	    p_eos = new Eos_fit_FPS(fich) ;
	    break ;
	}

	case 122 : {
	    p_eos = new Eos_fit_AkmalPR(fich) ;
	    break ;
	}

	default : {
	    cout << "Eos::eos_from_file : unknown type of EOS !" << endl ; 
	    cout << " identificator = " << identificator << endl ; 
	    abort() ; 
	    break ; 
	}
	
    }
    
    return p_eos ; 
    
}






}
