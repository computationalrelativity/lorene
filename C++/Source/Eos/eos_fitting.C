/*
 *  Method of class Eos_fitting
 *
 *    (see file eos_fitting.h for documentation).
 *
 */

/*
 *   Copyright (c) 2004 Keisuke Taniguchi
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

char eos_fitting_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2004/09/26 18:53:53  k_taniguchi
 * Initial revision
 *
 *
 * $Header$
 *
 */

// C headers
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Lorene headers
#include "headcpp.h"
#include "eos_fitting.h"
#include "eos.h"
#include "utilitaires.h"
#include "unites.h"

double fs(double) ;
double fd(double) ;
double dfs(double) ;
double dfd(double) ;

//************************************************************************

                    //--------------------------------//
                    //          Constructors          //
                    //--------------------------------//

// Standard constructor
// --------------------
Eos_fitting::Eos_fitting(const char* name_i, const char* data,
			 const char* path) : Eos(name_i) {

    strcpy(dataname, path) ;
    strcat(dataname, "/") ;
    strcat(dataname, data) ;

    read_coef() ;

}

// Constructor from a binary file
// ------------------------------
Eos_fitting::Eos_fitting(FILE* fich) : Eos(fich) {

    fread(dataname, sizeof(char), 160, fich) ;

    read_coef() ;

}

// Constructor from a formatted file
// ---------------------------------
Eos_fitting::Eos_fitting(ifstream& fich, const char* data) : Eos(fich) {

    char path[160] ;

    fich.getline(path, 160) ;

    strcpy(dataname, path) ;
    strcat(dataname, "/") ;
    strcat(dataname, data) ;

    read_coef() ;

}

// Destructor
Eos_fitting::~Eos_fitting() {

    delete [] pp ;

}


                     //---------------------------------------//
                     //              Outputs                  //
                     //---------------------------------------//

void Eos_fitting::sauve(FILE* fich) const {

    Eos::sauve(fich) ;

    fwrite(dataname, sizeof(char), 160, fich) ;

}
		  //-----------------------//
		  //	Miscellaneous	   //
		  //-----------------------//

void Eos_fitting::read_coef() {

    char blabla[120] ;

    ifstream fich(dataname) ;

    for (int i=0; i<3; i++) {      // Jump over the file header
        fich.getline(blabla, 120) ;
    }

    int nb_coef ;
    fich >> nb_coef ; fich.getline(blabla, 120) ;  // Number of coefficients

    for (int i=0; i<3; i++) {      // Jump over the table header
        fich.getline(blabla, 120) ;
    }

    pp = new double[nb_coef] ;

    for (int i=0; i<nb_coef; i++) {
        fich >> pp[i] ; fich.getline(blabla, 120) ;
    }

    fich.close() ;

}


			//------------------------------//
			//    Computational routines    //
			//------------------------------//

// Baryon density from enthalpy
//------------------------------

double Eos_fitting::nbar_ent_p(double ent, const Param* ) const {

  using namespace Unites ;

    if ( ent > double(0) ) {

        double aa = 0. ;
	double xx = 1. ;
	int m ;
	double yy ;
	double ent_value ;
	double rhob ;      /// Baryon density in the unit of c=G=Msol=1
	double nb ;        /// Number density in the unit of n_nuc=0.1fm^{-3}
	double trans_dens = msol_si / pow(g_si*msol_si/c_si/c_si,3.)
	  / rhonuc_si ;

	while (xx > 1.e-15) {

	    ent_value = 1. ;   // Initialization
	    xx = 0.1 * xx ;
	    m = 0 ;

	    while (ent_value > 1.e-15) {

	        m++ ;
		yy = aa + m * xx ;

		double aaa = 1.+pp[0]*pow(yy,pp[1])+pp[2]*pow(yy,pp[3]) ;
		double bbb = pow(1.+pp[4]*pow(yy,pp[5]),pp[6]) ;
		double ccc = pp[0]*pp[1]*pow(yy,pp[1])
		  +pp[2]*pp[3]*pow(yy,pp[3]) ;
		double ddd = pow(1.+pp[4]*pow(yy,pp[5]),pp[18]) ;
		double eee = -pp[7]*yy+pp[19] ;
		double fff = -pp[8]*yy+pp[20] ;
		double ggg = pp[4]*pp[5]*pp[6]*pow(yy,pp[5]) ;

		ent_value = exp(ent) - 1.0
		  -( aaa*bbb - 1. ) * fs(eee)
		  -pp[10]*pow(yy,pp[11])*fs(-eee)*fs(fff)
		  -pp[14]*pow(yy,pp[15])*fs(-fff)
		  -( ccc*bbb + aaa*ggg*ddd )*fs(eee)
		  +( aaa*bbb - 1.)*pp[7]*fd(eee)*yy
		  -pp[10]*pow(yy,pp[11])*(pp[11]*fs(-eee)*fs(fff)
					  +yy*(pp[7]*fd(-eee)*fs(fff)
					       -pp[8]*fs(-eee)*fd(fff)))
		  -pp[14]*pow(yy,pp[15])*(pp[15]*fs(-fff)
					  +pp[8]*yy*fd(-fff)) ;

	    }
	    aa += (m - 1) * xx ;
	}
	rhob = aa ;

	// The transformation from rhob to nb
	nb = rhob * trans_dens ;

	return nb ;
    }
    else {
        return 0 ;
    }

}

// Energy density from enthalpy
//------------------------------

double Eos_fitting::ener_ent_p(double ent, const Param* ) const {

  using namespace Unites ;

    if ( ent > double(0) ) {

        // Number density in the unit of [n_nuc]
        double nb = nbar_ent_p(ent) ;

	// The transformation from nb to yy
	// --------------------------------

	double trans_dens = msol_si / pow(g_si*msol_si/c_si/c_si,3.)
	  / rhonuc_si ; /// rho_b -> n_b

	// Baryon density in the unit of c=G=Msol=1
	double yy = nb / trans_dens ;

	double aaa = 1.+pp[0]*pow(yy,pp[1])+pp[2]*pow(yy,pp[3]) ;
	double bbb = pow(1.+pp[4]*pow(yy,pp[5]),pp[6]) ;
	double eee = -pp[7]*yy+pp[19] ;
	double fff = -pp[8]*yy+pp[20] ;

	double epsil = ( aaa*bbb - 1. ) * fs(eee)
	  +pp[10]*pow(yy,pp[11])*fs(-eee)*fs(fff)
	  +pp[14]*pow(yy,pp[15])*fs(-fff) ;

	// The transformation from epsil to ee
	// -----------------------------------

	// Energy density in the unit of [rho_nuc * c^2]
	double ee = nb * (1. + epsil) ;

	return ee ;
    }
    else {
        return 0 ;
    }

}

// Pressure from enthalpy
//------------------------

double Eos_fitting::press_ent_p(double ent, const Param* ) const {

  using namespace Unites ;

    if ( ent > double(0) ) {

        // Number density in the unit of [n_nuc]
        double nb = nbar_ent_p(ent) ;

	// The transformation from nb to yy
	// --------------------------------

	double trans_dens = msol_si / pow(g_si*msol_si/c_si/c_si,3.)
	  / rhonuc_si ; /// rho_b -> n_b

	// Baryon density in the unit of c=G=Msol=1
	double yy = nb / trans_dens ;

	double aaa = 1.+pp[0]*pow(yy,pp[1])+pp[2]*pow(yy,pp[3]) ;
	double bbb = pow(1.+pp[4]*pow(yy,pp[5]),pp[6]) ;
	double ccc = pp[0]*pp[1]*pow(yy,pp[1])+pp[2]*pp[3]*pow(yy,pp[3]) ;
	double ddd = pow(1.+pp[4]*pow(yy,pp[5]),pp[18]) ;
	double eee = -pp[7]*yy+pp[19] ;
	double fff = -pp[8]*yy+pp[20] ;
	double ggg = pp[4]*pp[5]*pp[6]*pow(yy,pp[5]) ;

	double ppp = yy*( ccc*bbb + aaa*ggg*ddd )*fs(eee)
	  -( aaa*bbb - 1. )*pp[7]*fd(eee)*yy*yy
	  +pp[10]*(pp[11]*pow(yy,pp[12])*fs(-eee)*fs(fff)
		   +pow(yy,pp[13])*(pp[7]*fd(-eee)*fs(fff)
				    -pp[8]*fs(-eee)*fd(fff)))
	  +pp[14]*(pp[15]*pow(yy,pp[16])*fs(-fff)
		   +pp[8]*pow(yy,pp[17])*fd(-fff)) ;

	// The transformation from ppp to pres
	// -----------------------------------

	// Pressure in the unit of [rho_nuc * c^2]
	double pres = ppp * trans_dens ;

	return pres ;
    }
    else {
        return 0 ;
    }

}

// dln(n)/dln(H) from enthalpy
//----------------------------

double Eos_fitting::der_nbar_ent_p(double ent, const Param* ) const {

  using namespace Unites ;

    if ( ent > double(0) ) {

        // Number density in the unit of [n_nuc]
        double nb = nbar_ent_p(ent) ;

	// The transformation from nb to yy
	// --------------------------------

	double trans_dens = msol_si / pow(g_si*msol_si/c_si/c_si,3.)
	  / rhonuc_si ; /// rho_b -> n_b

	// Baryon density in the unit of c=G=Msol=1
	double yy = nb / trans_dens ;

	double aaa = 1.+pp[0]*pow(yy,pp[1])+pp[2]*pow(yy,pp[3]) ;
	double bbb = pow(1.+pp[4]*pow(yy,pp[5]),pp[6]) ;
	double ccc = pp[0]*pp[1]*pow(yy,pp[1]) + pp[2]*pp[3]*pow(yy,pp[3]) ;
	double ddd = pow(1.+pp[4]*pow(yy,pp[5]),pp[18]) ;
	double eee = -pp[7]*yy+pp[19] ;
	double fff = -pp[8]*yy+pp[20] ;
	double ggg = pp[4]*pp[5]*pp[6]*pow(yy,pp[5]) ;
	double jjj = pp[0]*pp[1]*pp[1]*pow(yy,pp[1])
	  +pp[2]*pp[3]*pp[3]*pow(yy,pp[3]) ;

        double dlnsdlh = exp(ent) * ent /
	  ( ( ccc*bbb + aaa*ggg*ddd )*fs(eee)
	    -( aaa*bbb - 1. )*yy*pp[7]*dfs(eee)
	    +pp[10]*pow(yy,pp[11])*(pp[11]*fs(-eee)*fs(fff)
				    +yy*(pp[7]*dfs(-eee)*fs(fff)
					 -pp[8]*fs(-eee)*dfs(fff)))
	    +pp[14]*pow(yy,pp[15])*(pp[15]*fs(-fff)
				    +yy*pp[8]*dfs(-fff))  // xdf1/dx
	    +( jjj*bbb + 2.*ccc*ggg*ddd + aaa*ggg*pp[5]*ddd
	       + aaa*ggg*ggg*pp[18]/pp[6]
	       *pow(1.+pp[4]*pow(yy,pp[5]),pp[18]-1.) )*fs(eee)
	    -( ccc*bbb + aaa*ggg*ddd )*yy*pp[7]*dfs(eee)
	    -( ccc*bbb + aaa*ggg*ddd )*pp[7]*yy*fd(eee)
	    -pp[7]*( aaa*bbb - 1. )*yy*(-pp[7]*dfd(eee)*yy+fd(eee))
	    +pp[10]*(pp[11]*pow(yy,pp[11])*(pp[11]*fs(-eee)*fs(fff)
					    +yy*(pp[7]*dfs(-eee)*fs(fff)
						 -pp[8]*fs(-eee)*dfs(fff)))
		     +pp[12]*pow(yy,pp[12])*(pp[7]*fd(-eee)*fs(fff)
					     -pp[8]*fs(-eee)*fd(fff))
		      +pow(yy,pp[13])*(pp[7]*(pp[7]*dfd(-eee)*fs(fff)
					      -pp[8]*fd(-eee)*dfs(fff))
				       -pp[8]*(pp[7]*dfs(-eee)*fd(fff)
					       -pp[8]*fs(-eee)*dfd(fff))))
	    +pp[14]*(pp[15]*pow(yy,pp[15])*(pp[15]*fs(-fff)
					    +yy*pp[8]*dfs(-fff))
		     +pp[8]*pow(yy,pp[16])*(pp[16]*fd(-fff)
					    +yy*pp[8]*dfd(-fff)))
	    ) ;  // xd(f2/x)/dx

	return dlnsdlh ;

    }
    else {

        return double(1) / pp[11] ;  // To ensure continuity at ent=0

    }

}

// dln(e)/dln(H) from enthalpy
//----------------------------

double Eos_fitting::der_ener_ent_p(double ent, const Param* ) const {

  using namespace Unites ;

    if ( ent > double(0) ) {

         // Number density in the unit of [n_nuc]
        double nb = nbar_ent_p(ent) ;

	// The transformation from nb to yy
	// --------------------------------

	double trans_dens = msol_si / pow(g_si*msol_si/c_si/c_si,3.)
	  / rhonuc_si ; /// rho_b -> n_b

	// Baryon density in the unit of c=G=Msol=1
	double yy = nb / trans_dens ;

	double aaa = 1.+pp[0]*pow(yy,pp[1])+pp[2]*pow(yy,pp[3]) ;
	double bbb = pow(1.+pp[4]*pow(yy,pp[5]),pp[6]) ;
	double ccc = pp[0]*pp[1]*pow(yy,pp[1]) + pp[2]*pp[3]*pow(yy,pp[3]) ;
	double ddd = pow(1.+pp[4]*pow(yy,pp[5]),pp[18]) ;
	double eee = -pp[7]*yy+pp[19] ;
	double fff = -pp[8]*yy+pp[20] ;
	double ggg = pp[4]*pp[5]*pp[6]*pow(yy,pp[5]) ;

	double dlnsdlh = der_nbar_ent_p(ent) ;

	double dlesdlh = dlnsdlh
	  * (1. + ( (ccc*bbb + aaa*ggg*ddd)*fs(eee)
		    - (aaa*bbb-1.)*yy*pp[7]*dfs(eee)
		    +pp[10]*pow(yy,pp[11])*(pp[11]*fs(-eee)*fs(fff)
					    +yy*(pp[7]*dfs(-eee)*fs(fff)
						 -pp[8]*fs(-eee)*dfs(fff)))
		    +pp[14]*pow(yy,pp[15])*(pp[15]*fs(-fff)
					    +yy*pp[8]*dfs(-fff)) )
	     / ( 1. + (aaa*bbb-1.)*fs(eee)
		 + pp[10]*pow(yy,pp[11])*fs(-eee)*fs(fff)
		 + pp[14]*pow(yy,pp[15])*fs(-fff) ) ) ;

	return dlesdlh ;

    }
    else {

        return double(1) / pp[11] ;  // To ensure continuity at ent=0

    }

}

// dln(p)/dln(H) from enthalpy
//----------------------------

double Eos_fitting::der_press_ent_p(double ent, const Param* ) const {

  using namespace Unites ;

    if ( ent > double(0) ) {

         // Number density in the unit of [n_nuc]
        double nb = nbar_ent_p(ent) ;

	// The transformation from nb to yy
	// --------------------------------

	double trans_dens = msol_si / pow(g_si*msol_si/c_si/c_si,3.)
	  / rhonuc_si ; /// rho_b -> n_b

	// Baryon density in the unit of c=G=Msol=1
	double yy = nb / trans_dens ;

	double aaa = 1.+pp[0]*pow(yy,pp[1])+pp[2]*pow(yy,pp[3]) ;
	double bbb = pow(1.+pp[4]*pow(yy,pp[5]),pp[6]) ;
	double ccc = pp[0]*pp[1]*pow(yy,pp[1]) + pp[2]*pp[3]*pow(yy,pp[3]) ;
	double ddd = pow(1.+pp[4]*pow(yy,pp[5]),pp[18]) ;
	double eee = -pp[7]*yy+pp[19] ;
	double fff = -pp[8]*yy+pp[20] ;
	double ggg = pp[4]*pp[5]*pp[6]*pow(yy,pp[5]) ;
	double jjj = pp[0]*pp[1]*pp[1]*pow(yy,pp[1])
	  +pp[2]*pp[3]*pp[3]*pow(yy,pp[3]) ;

	double dlnsdlh = der_nbar_ent_p(ent) ;

	double dlpsdlh = dlnsdlh
	  * ( (ccc*bbb+aaa*ggg*ddd)*fs(eee)
	      +( jjj*bbb + 2.*ccc*ggg*ddd + aaa*ggg*pp[5]*ddd
		 + aaa*ggg*ggg*pp[18]/pp[6]
		 *pow(1.+pp[4]*pow(yy,pp[5]),pp[18]-1.) )*fs(eee)
	      -(ccc*bbb + aaa*ggg*ddd)*yy*pp[7]*dfs(eee)
	      -(ccc*bbb + aaa*ggg*ddd)*pp[7]*yy*fd(eee)
	      -(aaa*bbb-1.)*pp[7]*yy*(-pp[7]*dfd(eee)*yy+2.*fd(eee))
	      +pp[10]*(pp[11]*pow(yy,pp[11])*(pp[12]*fs(-eee)*fs(fff)
					      +yy*pp[7]*dfs(-eee)*fs(fff)
					      +-pp[8]*yy*fs(-eee)*dfs(fff))
		       +pp[13]*pow(yy,pp[12])*(pp[7]*fd(-eee)*fs(fff)
						  -pp[8]*fs(-eee)*fd(fff))
		       +pow(yy,pp[13])*(pp[7]*(pp[7]*dfd(-eee)*fs(fff)
					       -pp[8]*fd(-eee)*dfs(fff))
					-pp[8]*(pp[7]*dfs(-eee)*fd(fff)
						-pp[8]*fs(-eee)*dfd(fff))) )
	      +pp[14]*(pp[15]*pow(yy,pp[11])*(pp[16]*fs(-fff)
					      +yy*pp[8]*dfs(-fff))
		       +pp[8]*pow(yy,pp[16])*(pp[17]*pp[8]*dfd(-fff)
					      +yy*pp[8]*dfd(-fff)) ) )
	  / ( (ccc*bbb + aaa*ggg*ddd)*fs(eee) - (aaa*bbb-1.)*pp[7]*yy*fd(eee)
	      +pp[10]*pow(yy,pp[11])*(pp[11]*fs(-eee)*fs(fff)
				      +yy*(pp[7]*fd(-eee)*fs(fff)
					   -pp[8]*fs(-eee)*fd(fff)))
	      +pp[14]*pow(yy,pp[15])*(pp[15]*fs(-fff)
				      +pp[8]*yy*fd(-fff)) ) ;

	return dlpsdlh ;

    }
    else {

        return pp[12] / pp[11] ;  // To ensure continuity at ent=0

    }

}

//********************************************************
//     Functions which appear in the fitting formula
//********************************************************

 double fs(double xx) {

     double resu = double(1) / (double(1) + exp(xx)) ;

     return resu ;

 }

 double fd(double xx) {

     double resu = double(1) / pow(double(1)+exp(xx),double(2))
       - double(1) / (double(1) + exp(xx)) ;

     return resu ;

 }

 double dfs(double xx) {

     double resu ;

     if (xx > 0.) {
         resu = - exp(-xx) / pow(exp(-xx)+double(1),double(2)) ;
     }
     else {
         resu = - exp(xx) / pow(double(1)+exp(xx),double(2)) ;
     }

     return resu ;

 }

 double dfd(double xx) {

     double resu ;

     if (xx > 0.) {
         resu = - double(2) * exp(-2.*xx) /
	   pow(exp(-xx)+double(1),double(3))
	   + exp(-xx) / pow(exp(-xx)+double(1),double(2)) ;
     }
     else {
         resu = - double(2) * exp(xx) / 
	   pow(double(1)+exp(xx),double(3))
	   + exp(xx) / pow(double(1)+exp(xx),double(2)) ;
     }

     return resu ;

 }
