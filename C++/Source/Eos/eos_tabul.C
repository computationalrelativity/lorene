/*
 *  Methods of class Eos_tabul
 *
 *  (see file eos.h for documentation).
 *
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


char eos_tabul_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2002/04/09 14:32:15  e_gourgoulhon
 * 1/ Added extra parameters in EOS computational functions (argument par)
 * 2/ New class MEos for multi-domain EOS
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.3  2001/09/13  13:35:49  eric
 * Suppression des affichages dans read_table().
 *
 * Revision 2.2  2001/02/07  09:48:05  eric
 * Suppression de la fonction derent_ent_p.
 * Ajout des fonctions donnant les derivees de l'EOS:
 *      der_nbar_ent_p
 *      der_ener_ent_p
 *      der_press_ent_p
 *
 * Revision 2.1  2000/11/23  00:10:20  eric
 * Enthalpie minimale fixee a 1e-14.
 *
 * Revision 2.0  2000/11/22  19:31:30  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers C++
#include <iostream.h>
#include <fstream.h>

// headers C
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Headers Lorene
#include "eos.h"
#include "tbl.h"
#include "utilitaires.h"


void interpol_herm(const Tbl& , const Tbl&, const Tbl&, double, int&,
		   double&, double& ) ;


			//----------------------------//
			//   	Constructors	      //
			//----------------------------//

// Standard constructor
// --------------------			
Eos_tabul::Eos_tabul(const char* name_i, const char* table,
		     const char* path) : Eos(name_i) {	

	strcpy(tablename, path) ;
	strcat(tablename, "/") ;
	strcat(tablename, table) ;
	
	read_table() ; 	

}


// Constructor from binary file
// ----------------------------
Eos_tabul::Eos_tabul(FILE* fich) : Eos(fich) {

       fread(tablename, sizeof(char), 160, fich) ;		

       read_table() ;

}



// Constructor from a formatted file
// ---------------------------------
Eos_tabul::Eos_tabul(ifstream& fich, const char* table) : Eos(fich) {

	char path[160] ; 	
	
	fich.getline(path, 160) ;
	
	strcpy(tablename, path) ;
	strcat(tablename, "/") ;
	strcat(tablename, table) ;
	
	read_table() ; 	

}


			//--------------//
			//  Destructor  //
			//--------------//

Eos_tabul::~Eos_tabul(){
	delete logh ;
	delete logp ;
	delete dlpsdlh ;
}

			//------------//
			//  Outputs   //
			//------------//

void Eos_tabul::sauve(FILE* fich) const {

	Eos::sauve(fich) ;

    	fwrite(tablename, sizeof(char), 160, fich) ;		

}

			//------------------------//
			//  Reading of the table  //
			//------------------------//
			
void Eos_tabul::read_table() {

#include "unites.h"
	if (this == 0x0) {
		cout << mevpfm3 << km << msol << qpig << f_unit << endl ;
	}
    	
    	char blabla[120] ;
    	
	ifstream fich(tablename) ;
	
	for (int i=0; i<5; i++) {		//  jump over the file
    		fich.getline(blabla, 120) ;     //  header
    	}                                       //

    	int nbp ;
    	fich >> nbp ; fich.getline(blabla, 120) ;   // number of data
    	    	
	for (int i=0; i<3; i++) {		//  jump over the table
    		fich.getline(blabla, 120) ;     //  header
    	}                                       //

    	logh = new Tbl(nbp) ;
    	logp = new Tbl(nbp) ;
    	dlpsdlh = new Tbl(nbp) ;
    	
    	logh->set_etat_qcq() ;
    	logp->set_etat_qcq() ;
    	dlpsdlh->set_etat_qcq() ;
    	
    	
	double rhonuc_cgs = rhonuc_si * 1e-3 ;
	double c2_cgs = c_si * c_si * 1e4 ;    	
    	
    	int no ;
    	double nb_fm3, rho_cgs, p_cgs ;
    	
    	double ww = 0 ;
    	
    	for (int i=0; i<nbp; i++) {
    	
    		fich >> no ;
    		fich >> nb_fm3 ;
    		fich >> rho_cgs ;
    		fich >> p_cgs ; fich.getline(blabla, 120) ;    		
    		double psc2_cgs = p_cgs / c2_cgs ;
    		double h = log( (rho_cgs + psc2_cgs) /
    		                    (10 * nb_fm3 * rhonuc_cgs) ) ;
    		
    		if (i==0) {
    			ww = h ;
    		}
    		
    		h = h - ww + 1.e-14 ;    		
    		
//##
//	cout.precision(8) ;
// 	cout << i << "  " << rho_cgs << "  " << p_cgs << "  " << h
// 		  << "  "  << nb_fm3 << endl ;
//## 		
    		
    		logh->set(i) = log10( h ) ;
    		logp->set(i) = log10( psc2_cgs / rhonuc_cgs ) ;
    		dlpsdlh->set(i) = h * (rho_cgs + psc2_cgs) / psc2_cgs ;
    	
    	}
     	
     	hmin = pow( double(10), (*logh)(0) ) ;
     	hmax = pow( double(10), (*logh)(nbp-1) ) ;
     	
//##     	   	
//     	cout << "logh : " << endl ;
//     	cout << *logh << endl ;
//     	
//     	cout << "logp : " << endl ;
//     	cout << *logp << endl ;
//     	
//     	cout << "dlpsdlh : " << endl ;
//     	cout << *dlpsdlh << endl ;
//     	
//     	cout << "hmin, hmax : " << hmin << "  " << hmax << endl ;
//##     	
     	
     	fich.close();
   	
}


			//------------------------------//
			//    Computational routines    //
			//------------------------------//

// Baryon density from enthalpy
//------------------------------

double Eos_tabul::nbar_ent_p(double ent, const Param* ) const {

    static int i_near = logh->get_taille() / 2 ;

    if ( ent > hmin ) {
           if (ent > hmax) {
           	cout << "Eos_tabul::nbar_ent_p : ent > hmax !" << endl ;
           	abort() ;
           }
           double logh0 = log10( ent ) ;
           double logp0 ;
           double dlpsdlh0 ;
           interpol_herm(*logh, *logp, *dlpsdlh, logh0, i_near, logp0,
           		 dlpsdlh0) ;

           double pp = pow(double(10), logp0) ;

           return pp / ent * dlpsdlh0 * exp(-ent) ;
    }
    else{
	return 0 ;
    }
}

// Energy density from enthalpy
//------------------------------

double Eos_tabul::ener_ent_p(double ent, const Param* ) const {

    static int i_near = logh->get_taille() / 2 ;

    if ( ent > hmin ) {
           if (ent > hmax) {
           	cout << "Eos_tabul::ener_ent_p : ent > hmax !" << endl ;
           	abort() ;
           }
           double logh0 = log10( ent ) ;
           double logp0 ;
           double dlpsdlh0 ;
           interpol_herm(*logh, *logp, *dlpsdlh, logh0, i_near, logp0,
           		 dlpsdlh0) ;

           double pp = pow(double(10), logp0) ;

           return pp / ent * dlpsdlh0 - pp ;
    }
    else{
	return 0 ;
    }
}

// Pressure from enthalpy
//------------------------

double Eos_tabul::press_ent_p(double ent, const Param* ) const {

    static int i_near = logh->get_taille() / 2 ;

    if ( ent > hmin ) {
           if (ent > hmax) {
           	cout << "Eos_tabul::press_ent_p : ent > hmax !" << endl ;
           	abort() ;
           }

           double logh0 = log10( ent ) ;
           double logp0 ;
           double dlpsdlh0 ;
           interpol_herm(*logh, *logp, *dlpsdlh, logh0, i_near, logp0,
           		 dlpsdlh0) ;

           return pow(double(10), logp0) ;
    }
    else{
	return 0 ;
    }
}

// dln(n)/ln(H) from enthalpy 
//---------------------------

double Eos_tabul::der_nbar_ent_p(double ent, const Param* ) const {

    if ( ent > hmin ) {
           if (ent > hmax) {
           	cout << "Eos_tabul::der_nbar_ent_p : ent > hmax !" << endl ;
           	abort() ;
           }
	   
           cout << "Eos_tabul::der_nbar_ent_p : ent > hmax !" << endl ;
           abort() ;	       
	   return 0 ;
    }
    else{
	return 0 ;
    }
}


// dln(e)/ln(H) from enthalpy 
//---------------------------

double Eos_tabul::der_ener_ent_p(double ent, const Param* ) const {

    if ( ent > hmin ) {
           if (ent > hmax) {
           	cout << "Eos_tabul::der_ener_ent_p : ent > hmax !" << endl ;
           	abort() ;
           }

           cout << "Eos_tabul::der_ener_ent_p : not ready yet !" << endl ;
           abort() ;
	   return 0 ;
    }
    else{
	return 0 ;
    }
}


// dln(p)/ln(H) from enthalpy 
//---------------------------

double Eos_tabul::der_press_ent_p(double ent, const Param* ) const {

    if ( ent > hmin ) {
           if (ent > hmax) {
           	cout << "Eos_tabul::der_press_ent_p : ent > hmax !" << endl ;
           	abort() ;
           }

           cout << "Eos_tabul::der_press_ent_p : not ready yet !" << endl ;
           abort() ;
	   return 0 ;

    }
    else{
	return 0 ;
    }
}
