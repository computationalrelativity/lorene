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


char eos_from_file_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:27  e_gourgoulhon
 * Initial revision
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
 
// Headers C++
#include <iostream.h>
#include <fstream.h>

// Headers C
#include <stdlib.h>

// Header Lorene
#include "eos.h"

		//--------------------------------------//
		//  Identification virtual functions	//
		//--------------------------------------//


int Eos_poly::identify() const		{ return 1; }

int Eos_poly_newt::identify() const	{ return 2; }

int Eos_incomp::identify() const	{ return 3; }

int Eos_incomp_newt::identify() const	{ return 4; }

int Eos_strange::identify() const	{ return 5; }

int Eos_strange_cr::identify() const	{ return 6; }

int Eos_SLy4::identify() const		{ return 10; }

int Eos_FPS::identify() const		{ return 11; }

int Eos_BPAL12::identify() const	{ return 12; }

int Eos_AkmalPR::identify() const	{ return 13; }

int Eos_BBB2::identify() const		{ return 14; }

int Eos_BalbN1H1::identify() const	{ return 15; }


		//---------------------------------------------//
		//    EOS construction from a binary file      //
		//---------------------------------------------//

Eos* Eos::eos_from_file(FILE* fich) {
    
    Eos* p_eos ; 
    
    // Type (class) of EOS :
    int identificator ;     
    fread(&identificator, sizeof(int), 1, fich) ;		

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
    char blabla[80] ;

    // EOS identificator : 
    fich >> identificator ; fich.getline(blabla, 80) ;

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
	
	default : {
	    cout << "Eos::eos_from_file : unknown type of EOS !" << endl ; 
	    cout << " identificator = " << identificator << endl ; 
	    abort() ; 
	    break ; 
	}
	
    }
    
    return p_eos ; 
    
}






