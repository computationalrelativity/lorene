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
 * Revision 1.5  2014/10/13 08:52:53  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.4  2014/10/06 15:13:06  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.3  2005/05/23 14:14:00  k_taniguchi
 * Minor modification.
 *
 * Revision 1.2  2005/05/22 20:53:06  k_taniguchi
 * Modify the method to calculate baryon number density from enthalpy.
 *
 * Revision 1.1  2004/09/26 18:53:53  k_taniguchi
 * Initial revision
 *
 *
 * $Header$
 *
 */

// C headers
#include <cstdlib>
#include <cstring>
#include <cmath>

// Lorene headers
#include "headcpp.h"
#include "eos_fitting.h"
#include "eos.h"
#include "utilitaires.h"
#include "unites.h"

namespace Lorene {
double fc(double) ;
double gc(double) ;
double hc(double) ;

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
	double xx = 0.01 ;
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
		double ddd = pow(1.+pp[4]*pow(yy,pp[5]),pp[6]-1.) ;
		double eee = -pp[7]*yy+pp[9] ;
		double fff = -pp[8]*yy+pp[10] ;
		double ggg = pp[4]*pp[5]*pp[6]*pow(yy,pp[5]) ;

		ent_value = exp(ent) - 1.0
		  -( aaa*bbb - 1. ) * fc(eee)
		  -pp[11]*pow(yy,pp[12])*fc(-eee)*fc(fff)
		  -pp[13]*pow(yy,pp[14])*fc(-fff)
		  -( ccc*bbb + aaa*ddd*ggg )*fc(eee)
		  -yy*( aaa*bbb - 1. )*pp[7]*gc(eee)
		  -pp[11]*pp[12]*pow(yy,pp[12])*fc(-eee)*fc(fff)
		  +pp[11]*pow(yy,pp[12]+1.)*(pp[7]*gc(-eee)*fc(fff)
					     -pp[8]*fc(-eee)*gc(fff))
		  -pp[13]*pow(yy,pp[14])*(pp[14]*fc(-fff)
					  -pp[8]*yy*gc(-fff)) ;

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
	double eee = -pp[7]*yy+pp[9] ;
	double fff = -pp[8]*yy+pp[10] ;

	double epsil = ( aaa*bbb - 1. ) * fc(eee)
	  +pp[11]*pow(yy,pp[12])*fc(-eee)*fc(fff)
	  +pp[13]*pow(yy,pp[14])*fc(-fff) ;

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
	double ddd = pow(1.+pp[4]*pow(yy,pp[5]),pp[6]-1.) ;
	double eee = -pp[7]*yy+pp[9] ;
	double fff = -pp[8]*yy+pp[10] ;
	double ggg = pp[4]*pp[5]*pp[6]*pow(yy,pp[5]) ;

	double ppp = yy*( ccc*bbb + aaa*ddd*ggg )*fc(eee)
	  +yy*yy*( aaa*bbb - 1. )*pp[7]*gc(eee)
	  +pp[11]*pp[12]*pow(yy,pp[12]+1.)*fc(-eee)*fc(fff)
	  -pp[11]*pow(yy,pp[12]+2.)*(pp[7]*gc(-eee)*fc(fff)
				     -pp[8]*fc(-eee)*gc(fff))
	  +pp[13]*pow(yy,pp[14]+1.)*(pp[14]*fc(-fff)-pp[8]*yy*gc(-fff)) ;

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
	double ddd = pow(1.+pp[4]*pow(yy,pp[5]),pp[6]-1.) ;
	double eee = -pp[7]*yy+pp[9] ;
	double fff = -pp[8]*yy+pp[10] ;
	double ggg = pp[4]*pp[5]*pp[6]*pow(yy,pp[5]) ;
	double jjj = pp[0]*pp[1]*pp[1]*pow(yy,pp[1])
	  +pp[2]*pp[3]*pp[3]*pow(yy,pp[3]) ;

        double dlnsdlh = exp(ent) * ent /
	  ( ( ccc*bbb + aaa*ddd*ggg )*( fc(eee) + 2.*yy*pp[7]*gc(eee) )
	    +yy*pp[7]*( aaa*bbb - 1. )*( 2.*gc(eee) - yy*pp[7]*hc(eee) )
	    +pp[11]*pp[12]*(pp[12]+1.)*pow(yy,pp[12])*fc(-eee)*fc(fff)
	    +2.*pp[11]*(pp[12]+1.)*pow(yy,pp[12]+1.)
	      *( -pp[7]*gc(-eee)*fc(fff) + pp[8]*fc(-eee)*gc(fff) )
	    -pp[11]*pow(yy,pp[12]+2.)*( pp[7]*pp[7]*hc(-eee)*fc(fff)
					+2.*pp[7]*pp[8]*gc(-eee)*gc(fff)
					+pp[8]*pp[8]*fc(-eee)*hc(fff) )
	    +pp[13]*pp[14]*(pp[14]+1.)*pow(yy,pp[14])*fc(-fff)
	    -2.*pp[8]*pp[13]*(pp[14]+1.)*pow(yy,pp[14]+1.)*gc(-fff)
	    -pp[8]*pp[8]*pp[13]*pow(yy,pp[14]+2.)*hc(-fff)
	    +( jjj*bbb + 2.*ccc*ddd*ggg
	       + aaa*pow(1.+pp[4]*pow(yy,pp[5]),pp[6]-2.)
	       *ggg*(ggg+pp[5]) )*fc(eee)
	    ) ;

	return dlnsdlh ;

    }
    else {

        return double(1) / pp[12] ;  // To ensure continuity at ent=0

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
	double ddd = pow(1.+pp[4]*pow(yy,pp[5]),pp[6]-1.) ;
	double eee = -pp[7]*yy+pp[9] ;
	double fff = -pp[8]*yy+pp[10] ;
	double ggg = pp[4]*pp[5]*pp[6]*pow(yy,pp[5]) ;

	double dlnsdlh = der_nbar_ent_p(ent) ;

	double dlesdlh = dlnsdlh
	  * (1. + ( ( ccc*bbb + aaa*ddd*ggg )*fc(eee)
		    +yy*pp[7]*( aaa*bbb - 1. )*gc(eee)
		    +pp[11]*pp[12]*pow(yy,pp[12])*fc(-eee)*fc(fff)
		    +pp[11]*pow(yy,pp[12]+1.)*( -pp[7]*gc(-eee)*fc(fff)
						+pp[8]*fc(-eee)*gc(fff) )
		    +pp[13]*pp[14]*pow(yy,pp[14])*fc(-fff)
		    -pp[8]*pp[13]*pow(yy,pp[14]+1.)*gc(-fff) )
	     / ( 1. + ( aaa*bbb - 1. )*fc(eee)
		 + pp[11]*pow(yy,pp[12])*fc(-eee)*fc(fff)
		 + pp[13]*pow(yy,pp[14])*fc(-fff) ) ) ;

	return dlesdlh ;

    }
    else {

        return double(1) / pp[12] ;  // To ensure continuity at ent=0

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
	double ddd = pow(1.+pp[4]*pow(yy,pp[5]),pp[6]-1.) ;
	double eee = -pp[7]*yy+pp[9] ;
	double fff = -pp[8]*yy+pp[10] ;
	double ggg = pp[4]*pp[5]*pp[6]*pow(yy,pp[5]) ;
	double jjj = pp[0]*pp[1]*pp[1]*pow(yy,pp[1])
	  +pp[2]*pp[3]*pp[3]*pow(yy,pp[3]) ;

	double dlnsdlh = der_nbar_ent_p(ent) ;

	double dlpsdlh = dlnsdlh
	  * ( ( ccc*bbb + aaa*ddd*ggg )*fc(eee)
	      +( jjj*bbb + 2.*ccc*ddd*ggg
		 + aaa*pow(1.+pp[4]*pow(yy,pp[5]),pp[6]-2.)*ggg*(ggg+pp[5])
		 )*fc(eee)
	      +2.*yy*pp[7]*( ccc*bbb + aaa*ddd*ggg )*gc(eee)
	      +yy*pp[7]*( aaa*bbb - 1. )*(2.*gc(eee) - yy*pp[7]*hc(eee))
	      +pp[11]*pp[12]*(pp[12]+1.)*pow(yy,pp[12])*fc(-eee)*fc(fff)
	      +2.*pp[11]*(pp[12]+1.)*pow(yy,pp[12]+1.)
	        *( -pp[7]*gc(-eee)*fc(fff) + pp[8]*fc(-eee)*gc(fff) )
	      -pp[11]*pow(yy,pp[12]+2.)*( pp[7]*pp[7]*hc(-eee)*fc(fff)
					  +2.*pp[7]*pp[8]*gc(-eee)*gc(fff)
					  +pp[8]*pp[8]*fc(-eee)*hc(fff) )
	      +pp[13]*(pp[14]+1.)*pow(yy,pp[14])*( pp[14]*fc(-fff)
						   -2.*pp[8]*yy*gc(-fff) )
	      -pp[8]*pp[8]*pp[13]*pow(yy,pp[14]+2.)*hc(-fff) )
	  / ( ( ccc*bbb + aaa*ddd*ggg )*fc(eee)
	      +yy*pp[7]*( aaa*bbb - 1. )*gc(eee)
	      +pp[11]*pow(yy,pp[12])*( pp[12]*fc(-eee)*fc(fff)
				       -yy*pp[7]*gc(-eee)*fc(fff)
				       +yy*pp[8]*fc(-eee)*gc(fff) )
	      +pp[13]*pow(yy,pp[14])*( pp[14]*fc(-fff)
				       -yy*pp[8]*gc(-fff) ) ) ;

	return dlpsdlh ;

    }
    else {

        return (pp[12] + 1.) / pp[12] ;  // To ensure continuity at ent=0

    }

}

//********************************************************
//     Functions which appear in the fitting formula
//********************************************************

double fc(double xx) {

    double resu = double(1) / (double(1) + exp(xx)) ;

    return resu ;

}

double gc(double xx) {

    double resu ;

    if (xx > 0.) {
        resu = exp(-xx) / pow(exp(-xx)+double(1),double(2)) ;
    }
    else {
        resu = exp(xx) / pow(double(1)+exp(xx),double(2)) ;
    }

    return resu ;

}

 double hc(double xx) {

     double resu ;

     if (xx > 0.) {
         resu = exp(-xx) * (exp(-xx)-double(1)) /
	     pow(exp(-xx)+double(1),double(3)) ;
     }
     else {
         resu = exp(xx) * (double(1)-exp(xx)) / 
	     pow(double(1)+exp(xx),double(3)) ;
     }

     return resu ;

 }
}
