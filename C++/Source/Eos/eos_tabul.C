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
 * Revision 1.8  2004/03/25 10:29:02  j_novak
 * All LORENE's units are now defined in the namespace Unites (in file unites.h).
 *
 * Revision 1.7  2003/11/25 13:42:50  m_bejger
 * read_table written in more ordered way
 *
 * Revision 1.6  2003/11/21 16:11:16  m_bejger
 * added the computation dlnp/dlnn_b, dlnn/dlnH
 *
 * Revision 1.5  2003/05/30 07:50:06  e_gourgoulhon
 *
 * Reformulate the computation of der_nbar_ent
 * Added computation of der_press_ent_p.
 *
 * Revision 1.4  2003/05/15 09:42:47  e_gourgoulhon
 *
 * Now computes d ln / dln H.
 *
 * Revision 1.3  2002/10/16 14:36:35  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
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

// headers C
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Headers Lorene
#include "eos.h"
#include "tbl.h"
#include "utilitaires.h"
#include "unites.h"


void interpol_herm(const Tbl& , const Tbl&, const Tbl&, double, int&,
		   double&, double& ) ;

void interpol_linear(const Tbl&, const Tbl&, double, int&, double&) ;

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
        delete lognb ;
        delete dlpsdlnb ;
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

  using namespace Unites ;
    	
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

        press = new double[nbp] ;
        nb    = new double[nbp] ;
        ro    = new double[nbp] ; 
 
    	logh = new Tbl(nbp) ;
    	logp = new Tbl(nbp) ;
    	dlpsdlh = new Tbl(nbp) ;
	lognb = new Tbl(nbp) ;
        dlpsdlnb = new Tbl(nbp) ;
    	
    	logh->set_etat_qcq() ;
    	logp->set_etat_qcq() ;
    	dlpsdlh->set_etat_qcq() ;
        lognb->set_etat_qcq() ;
        dlpsdlnb->set_etat_qcq() ;   	
    	
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
    		
                press[i] = psc2_cgs ; 
                nb[i]    = nb_fm3 ;
                ro[i]    = rho_cgs ; 

    		if (i==0) { ww = h ; }
    		
    		h = h - ww + 1.e-14 ;    		
    		
//##
//	cout.precision(8) ;
// 	cout << i << "  " << rho_cgs << "  " << p_cgs << "  " << h
// 		  << "  "  << nb_fm3 << endl ;
//## 		
    		
    		logh->set(i) = log10( h ) ;
    		logp->set(i) = log10( psc2_cgs / rhonuc_cgs ) ;
    		dlpsdlh->set(i) = h * (rho_cgs + psc2_cgs) / psc2_cgs ;
		lognb->set(i) = log10(nb_fm3) ;

    	}

        double p0, p1, p2, n0, n1, n2, dpdnb; 

	// special case: i=0

          p0 = log(press[0]);
          p1 = log(press[1]);
          p2 = log(press[2]);

          n0 = log(nb[0]);
          n1 = log(nb[1]);
          n2 = log(nb[2]);

          dpdnb = p0*(2*n0-n1-n2)/(n0-n1)/(n0-n2) +
	                 p1*(n0-n2)/(n1-n0)/(n1-n2) +
                         p2*(n0-n1)/(n2-n0)/(n2-n1) ;

          dlpsdlnb->set(0) = dpdnb ; 

/*
	  cout << exp(n0) << "  " << ro[0] << "  " 
               << c2_cgs*exp(p0) << "  " << dpdnb << endl ;
*/
	for(int i=1;i<nbp-1;i++) { 

          p0 = log(press[i-1]);
          p1 = log(press[i]);
          p2 = log(press[i+1]);

          n0 = log(nb[i-1]);
          n1 = log(nb[i]);
          n2 = log(nb[i+1]);

          dpdnb = p0*(n1-n2)/(n0-n1)/(n0-n2) +
	                 p1*(2*n1-n0-n2)/(n1-n0)/(n1-n2) +
                         p2*(n1-n0)/(n2-n0)/(n2-n1) ;

	  dlpsdlnb->set(i) = dpdnb ;

/*
	  cout << exp(n1) << "  " << ro[i] << "  " 
               << c2_cgs*exp(p1) << "  " << dpdnb << endl ;
*/

        } 
     	
	// special case: i=nbp-1

          p0 = log(press[nbp-3]);
          p1 = log(press[nbp-2]);
          p2 = log(press[nbp-1]);

          n0 = log(nb[nbp-3]);
          n1 = log(nb[nbp-2]);
          n2 = log(nb[nbp-1]);

          dpdnb = p0*(n2-n1)/(n0-n1)/(n0-n2) +
	                 p1*(n2-n0)/(n1-n0)/(n1-n2) +
                         p2*(2*n2-n0-n1)/(n2-n0)/(n2-n1) ;

	  dlpsdlnb->set(nbp-1) = dpdnb ;

/*
	  cout << exp(n2) << "  " << ro[nbp-1] << "  " 
               << c2_cgs*exp(p2) << "  " << dpdnb << endl ;
*/
                
     	hmin = pow( double(10), (*logh)(0) ) ;
     	hmax = pow( double(10), (*logh)(nbp-1) ) ;
     	
//##     	   	
//     	cout << "logh : " << endl ;
//     	cout << *logh << endl ;
//     	
//
//     	cout << "logp : " << endl ;
//     	cout << *logp << endl ;
//     	
//     	cout << "dlpsdlh : " << endl ;
//     	cout << *dlpsdlh << endl ;
//
//     	cout << "dlpsdlnb : " << endl ;
//     	cout << *dlpsdlnb << endl ;
//
//     	
//     	cout << "hmin, hmax : " << hmin << "  " << hmax << endl ;
//##     	
     	
     	fich.close();
 
        delete [] press ; 
        delete [] nb ; 
        delete [] ro ; 
   	
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

	   double zeta = der_press_ent_p(ent) / der_press_nbar_p(ent) ; 

	   return zeta ; 
	  
    }
    else 

    return 1./(der_press_nbar_p(hmin)-1.) ;  // to ensure continuity

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

    static int i_near = logh->get_taille() / 2 ;

    if ( ent > hmin ) {
           if (ent > hmax) {
           	cout << "Eos_tabul::der_press_ent_p : ent > hmax !" << endl ;
           	abort() ;
           }

           double logh0 = log10( ent ) ;
           double logp0 ;
           double dlpsdlh0 ;
           interpol_herm(*logh, *logp, *dlpsdlh, logh0, i_near, logp0,
           		 dlpsdlh0) ;

           return dlpsdlh0 ;

    }
    else{
        
        return 0 ; 
	// return der_press_ent_p(hmin) ; // to ensure continuity
    }
}


// dln(p)/dln(n) from enthalpy 
//---------------------------

double Eos_tabul::der_press_nbar_p(double ent, const Param*) const {

    static int i_near = logh->get_taille() / 2 ;

    if ( ent > hmin ) {
           if (ent > hmax) {
           	cout << "Eos_tabul::der_press_nbar_p : ent > hmax !" << endl ;
           	abort() ;
           }

           double logh0 = log10( ent ) ;
           double dlpsdlnb0 ;

           interpol_linear(*logh, *dlpsdlnb, logh0, i_near, dlpsdlnb0) ;
 
//	   cout << "gamma: " << dlpsdlnb0 << " ent: " << ent << endl; 

          return dlpsdlnb0 ;

    }
    else{
        
         return (*dlpsdlnb)(0) ; 
	 // return der_press_nbar_p(hmin) ; // to ensure continuity

    }

          
}
