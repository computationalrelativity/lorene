/*
 *  Methods of class Eos_multi_poly
 *
 *    (see file eos_multi_poly.h for documentation).
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

char eos_multi_poly_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2004/05/07 11:55:59  k_taniguchi
 * Change the searching procedure of the baryon density.
 *
 * Revision 1.1  2004/05/07 08:10:58  k_taniguchi
 * Initial revision
 *
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
#include "eos_multi_poly.h"
#include "eos.h"
#include "utilitaires.h"
#include "param.h"


                     //--------------------------------------//
                     //              Constructors            //
                     //--------------------------------------//

Eos_multi_poly::Eos_multi_poly(int ndom_p, double* gamma_p, double* kappa_p,
			       int ndom_e, double* gamma_e, double* kappa_e)
    : Eos("Multi-polytropic EOS"),
      ndom_press(ndom_p), ndom_epsil(ndom_e),
      m_0(double(1)), mu_0(double(1)) {

    gam_p = new double [ndom_press] ;

    for (int l=0; l<ndom_press; l++) {
        gam_p[l] = gamma_p[l] ;
    }

    gam_e = new double [ndom_epsil] ;

    for (int l=0; l<ndom_epsil; l++) {
        gam_e[l] = gamma_e[l] ;
    }

    kap_p = new double [ndom_press] ;

    for (int l=0; l<ndom_press; l++) {
        kap_p[l] = kappa_p[l] ;
    }

    kap_e = new double [ndom_epsil] ;

    for (int l=0; l<ndom_epsil; l++) {
        kap_e[l] = kappa_e[l] ;
    }

    set_auxiliary() ;

}

// Copy constructor
Eos_multi_poly::Eos_multi_poly(const Eos_multi_poly& eosmp)
    : Eos(eosmp),
      ndom_press(eosmp.ndom_press), ndom_epsil(eosmp.ndom_epsil),
      m_0(eosmp.m_0), mu_0(eosmp.mu_0) {

    gam_p = new double [ndom_press] ;

    for (int l=0; l<ndom_press; l++) {
        gam_p[l] = eosmp.gam_p[l] ;
    }

    gam_e = new double [ndom_epsil] ;

    for (int l=0; l<ndom_epsil; l++) {
        gam_e[l] = eosmp.gam_e[l] ;
    }

    kap_p = new double [ndom_press] ;

    for (int l=0; l<ndom_press; l++) {
        kap_p[l] = eosmp.kap_p[l] ;
    }

    kap_e = new double [ndom_epsil] ;

    for (int l=0; l<ndom_epsil; l++) {
        kap_e[l] = eosmp.kap_e[l] ;
    }

    set_auxiliary() ;

}

//  Constructor from a binary file
Eos_multi_poly::Eos_multi_poly(FILE* fich) : Eos(fich) {

    fread_be(&ndom_press, sizeof(int), 1, fich) ;
    fread_be(&ndom_epsil, sizeof(int), 1, fich) ;

    gam_p = new double [ndom_press] ;

    for (int l=0; l<ndom_press; l++) {
        fread_be(&gam_p[l], sizeof(double), 1, fich) ;
    }

    gam_e = new double [ndom_epsil] ;

    for (int l=0; l<ndom_epsil; l++) {
        fread_be(&gam_e[l], sizeof(double), 1, fich) ;
    }

    kap_p = new double [ndom_press] ;

    for (int l=0; l<ndom_press; l++) {
        fread_be(&kap_p[l], sizeof(double), 1, fich) ;
    }

    kap_e = new double [ndom_epsil] ;

    for (int l=0; l<ndom_epsil; l++) {
        fread_be(&kap_e[l], sizeof(double), 1, fich) ;
    }

    set_auxiliary() ;

}

//  Constructor from a formatted file
Eos_multi_poly::Eos_multi_poly(ifstream& fich) : Eos(fich) {

    char blabla[80] ;

    fich >> ndom_press ; fich.getline(blabla, 80) ;
    fich >> ndom_epsil ; fich.getline(blabla, 80) ;

    gam_p = new double [ndom_press] ;

    for (int l=0; l<ndom_press; l++) {
        fich >> gam_p[l] ; fich.getline(blabla, 80) ;
    }

    gam_e = new double [ndom_epsil] ;

    for (int l=0; l<ndom_epsil; l++) {
        fich >> gam_e[l] ; fich.getline(blabla, 80) ;
    }

    kap_p = new double [ndom_press] ;

    for (int l=0; l<ndom_press; l++) {
        fich >> kap_p[l] ; fich.getline(blabla, 80) ;
    }

    kap_e = new double [ndom_epsil] ;

    for (int l=0; l<ndom_epsil; l++) {
        fich >> kap_e[l] ; fich.getline(blabla, 80) ;
    }

    set_auxiliary() ;

}

// Destructor
Eos_multi_poly::~Eos_multi_poly() {

    delete [] nb_press ;
    delete [] nb_epsil ;
    delete [] ent_crit_p ;
    delete [] ent_crit_e ;
    delete [] gam_p ;
    delete [] gam_e ;
    delete [] kap_p ;
    delete [] kap_e ;

}

			//--------------//
			//  Assignment  //
			//--------------//

void Eos_multi_poly::operator=(const Eos_multi_poly& ) {

        cout << "Eos_multi_poly::operator=  : not implemented yet !" << endl ;
                abort() ;

}

		  //-----------------------//
		  //	Miscellaneous	   //
		  //-----------------------//

void Eos_multi_poly::set_auxiliary() {

    nb_press = new double [ndom_press-1] ;

    for (int l=0; l<ndom_press-1; l++)
	nb_press[l] = pow( kap_p[l+1]/kap_p[l],
			   double(1)/(gam_p[l]-gam_p[l+1]) ) ;

    nb_epsil = new double [ndom_epsil-1] ;

    for (int l=0; l<ndom_epsil-1; l++)
	nb_epsil[l] = pow( kap_e[l+1]/kap_e[l],
			   double(1)/(gam_e[l]-gam_e[l+1]) ) ;

    ent_crit_p = new double [ndom_press-1] ;

    int i = 0 ;
    for (int l=0; l<ndom_press-1; l++) {
        if (nb_press[l] > nb_epsil[ndom_epsil-2]) {
	    i = ndom_epsil - 1 ;
	}
	else {
	    while (nb_press[l] > nb_epsil[i]) {
	        i++ ;
	    }
	}
	assert(i<=ndom_epsil-1) ;

        double kp_e = kap_e[i] ;
	double gm_e = gam_e[i] ;
        ent_crit_p[l] = log(double(1)+kp_e*pow(nb_press[l],gm_e)
			    +kap_p[l]*pow(nb_press[l],gam_p[l]-double(1))) ;
    }

    ent_crit_e = new double [ndom_epsil-1] ;

    int j = 0 ;
    for (int l=0; l<ndom_epsil-1; l++) {
        if (nb_epsil[l] > nb_press[ndom_press-2]) {
	    j = ndom_press - 1 ;
	}
	else {
	    while (nb_epsil[l] > nb_press[j]) {
	        j++ ;
	    }
	}
	assert(j<=ndom_press-1) ;

        double kp_p = kap_p[j] ;
	double gm_p = gam_p[j] ;
        ent_crit_e[l] = log(double(1)+kap_e[l]*pow(nb_epsil[l],gam_e[l])
			    +kp_p*pow(nb_epsil[l],gm_p-double(1))) ;
    }
}


			//------------------------//
			//  Comparison operators  //
			//------------------------//

bool Eos_multi_poly::operator==(const Eos& eos_i) const {
    
    bool resu = true ; 
    
    if ( eos_i.identify() != identify() ) {
	cout << "The second EOS is not of type Eos_multi_poly !" << endl ; 
	resu = false ; 
    }
    else{
	
	const Eos_multi_poly& eos = dynamic_cast<const Eos_multi_poly&>(eos_i) ; 

	if (eos.get_ndom_press() != ndom_press) {
	    cout << "The two Eos_multi_poly have "
		 << "different number of domains for the pressure : "
		 << ndom_press << " <-> " << eos.get_ndom_press() << endl ; 
	    resu = false ; 
	}

	if (eos.get_ndom_epsil() != ndom_epsil) {
	    cout << "The two Eos_multi_poly have "
		 << "different number of domains for the energy density : "
		 << ndom_epsil << " <-> " << eos.get_ndom_epsil() << endl ; 
	    resu = false ; 
	}

	for (int l=0; l<ndom_press; l++) {
	    if (eos.get_gam_press(l) != gam_p[l]) {
	        cout << "The two Eos_poly have different gamma "
		     << "for the pressure : " << gam_p[l] << " <-> " 
		     << eos.get_gam_press(l) << endl ;
	    resu = false ;
	    }
	}

	for (int l=0; l<ndom_epsil; l++) {
	    if (eos.get_gam_epsil(l) != gam_e[l]) {
	        cout << "The two Eos_poly have different gamma "
		     << "for the energy density : " << gam_e[l] << " <-> " 
		     << eos.get_gam_epsil(l) << endl ;
	    resu = false ;
	    }
	}

	for (int l=0; l<ndom_press; l++) {
	    if (eos.get_kap_press(l) != kap_p[l]) {
	        cout << "The two Eos_poly have different kappa "
		     << "for the pressure : " << kap_p[l] << " <-> " 
		     << eos.get_kap_press(l) << endl ;
	    resu = false ;
	    }
	}

	for (int l=0; l<ndom_epsil; l++) {
	    if (eos.get_kap_epsil(l) != kap_e[l]) {
	        cout << "The two Eos_poly have different kappa "
		     << "for the energy density : " << kap_e[l] << " <-> " 
		     << eos.get_kap_epsil(l) << endl ;
	    resu = false ;
	    }
	}

    }
    
    return resu ; 
    
}

bool Eos_multi_poly::operator!=(const Eos& eos_i) const {
 
    return !(operator==(eos_i)) ; 
       
}

                     //---------------------------------------//
                     //              Outputs                  //
                     //---------------------------------------//

void Eos_multi_poly::sauve(FILE* fich) const {

    Eos::sauve(fich) ;

    fwrite_be(&ndom_press, sizeof(int), 1, fich) ;
    fwrite_be(&ndom_epsil, sizeof(int), 1, fich) ;

    for (int l=0; l<ndom_press; l++) {
        fwrite_be(&gam_p[l], sizeof(double), 1, fich) ;
    }

    for (int l=0; l<ndom_epsil; l++) {
        fwrite_be(&gam_e[l], sizeof(double), 1, fich) ;
    }

    for (int l=0; l<ndom_press; l++) {
        fwrite_be(&kap_p[l], sizeof(double), 1, fich) ;
    }

    for (int l=0; l<ndom_epsil; l++) {
        fwrite_be(&kap_e[l], sizeof(double), 1, fich) ;
    }

}

ostream& Eos_multi_poly::operator>>(ostream & ost) const {

    ost << "EOS of class Eos_multi_poly "
	<< "(multi-domain polytropic equation of state) : " << endl ;

    ost << "   Number of domains for the pressure       :      "
	<< ndom_press << endl ;

    ost << "   Number of domains for the energy density :      "
	<< ndom_epsil << endl ;

    for (int l=0; l<ndom_press; l++) {
        ost << "EOS for the pressure in domain " << l << " : " << endl ;
	ost << "----------------------------------" << endl ;
	ost << "   gamma_press :" << gam_p[l] << endl ;
	ost << "   kappa_press :" << kap_p[l] << endl ;
    }

    for (int l=0; l<ndom_epsil; l++) {
        ost << "EOS for the energy density in domain " << l << " : " << endl ;
	ost << "----------------------------------------" << endl ;
	ost << "   gamma_epsil :" << gam_e[l] << endl ;
	ost << "   kappa_epsil :" << kap_e[l] << endl ;
    }

    return ost ;

}


			//------------------------------//
			//    Computational routines    //
			//------------------------------//

// Baryon density from enthalpy
//------------------------------

double Eos_multi_poly::nbar_ent_p(double ent, const Param* ) const {

    if ( ent > double(0) ) {

        int i = 0 ;

	if ( ent >= ent_crit_p[ndom_press-2] ) {
	    i = ndom_press - 1 ;
	}
	else {
	    while ( ent > ent_crit_p[i] ) {
	        i++ ;
	    }
	}
	assert(i <= ndom_press-1) ;

	double kap_pr = kap_p[i] ;
	double gam_pr = gam_p[i] ;

	int j = 0 ;

	if ( ent >= ent_crit_e[ndom_epsil-2] ) {
	    j = ndom_epsil - 1 ;
	}
	else {
	    while ( ent > ent_crit_e[j] ) {
	        j++ ;
	    }
	}
	assert(j <= ndom_epsil-1) ;

        double kap_ep = kap_e[j] ;
	double gam_ep = gam_e[j] ;

	double nb = 1. ;
	double nb_m1 = 0. ;

	// Switch the method of computation of the baryon density
	double diff_gamma = gam_ep - gam_pr + double(1) ;

	if (diff_gamma <= 0.) {
	    while (abs(1.-nb_m1/nb) > 1.e-15) {
	        nb_m1 = nb ;
		nb = pow( (exp(ent) - double(1))
			  / (kap_pr + kap_ep*pow(nb_m1, diff_gamma)),
			  double(1)/(gam_pr-double(1)) ) ;
	    }
	}
	else {
	    double aa = 0. ;
	    double xx = 1. ;
	    int m ;
	    double yy ;
	    double ent_value ;

	    while (xx > 1.e-15) {

	      ent_value = 1. ;   // Initialization
	      xx = 0.1 * xx ;
	      m = 0 ;

	      while (ent_value > 1.e-15) {

		m++ ;
		yy = aa + m * xx ;
		ent_value = exp(ent) - double(1) - kap_ep*pow(yy,gam_ep)
		  - kap_pr*pow(yy,gam_pr-double(1)) ;

	      }
	      aa += (m - 1) * xx ;
	    }
	    nb = aa ;
	}

	return nb ;
    }
    else {
        return 0 ;
    }

}

// Energy density from enthalpy
//------------------------------

double Eos_multi_poly::ener_ent_p(double ent, const Param* ) const {

    if ( ent > double(0) ) {

        int i = 0 ;

	if ( ent >= ent_crit_p[ndom_press-2] ) {
	    i = ndom_press - 1 ;
	}
	else {
	    while ( ent > ent_crit_p[i] ) {
	        i++ ;
	    }
	}
	assert(i <= ndom_press-1) ;

	double kap_pr = kap_p[i] ;
	double gam_pr = gam_p[i] ;

	int j = 0 ;

	if ( ent >= ent_crit_e[ndom_epsil-2] ) {
	    j = ndom_epsil - 1 ;
	}
	else {
	    while ( ent > ent_crit_e[j] ) {
	        j++ ;
	    }
	}
	assert(j <= ndom_epsil-1) ;

        double kap_ep = kap_e[j] ;
	double gam_ep = gam_e[j] ;

	double nb = 1. ;
	double nb_m1 = 0. ;

	// Switch the method of computation of the baryon density
	double diff_gamma = gam_ep - gam_pr + double(1) ;

	if (diff_gamma <= 0.) {
	    while (abs(1.-nb_m1/nb) > 1.e-15) {
	        nb_m1 = nb ;
		nb = pow( (exp(ent) - double(1))
			  / (kap_pr + kap_ep*pow(nb_m1, diff_gamma)),
			  double(1)/(gam_pr-double(1)) ) ;
	    }
	}
	else {
	    double aa = 0. ;
	    double xx = 1. ;
	    int m ;
	    double yy ;
	    double ent_value ;

	    while (xx > 1.e-15) {

	      ent_value = 1. ;   // Initialization
	      xx = 0.1 * xx ;
	      m = 0 ;

	      while (ent_value > 1.e-15) {

		m++ ;
		yy = aa + m * xx ;
		ent_value = exp(ent) - double(1) - kap_ep*pow(yy,gam_ep)
		  - kap_pr*pow(yy,gam_pr-double(1)) ;

	      }
	      aa += (m - 1) * xx ;
	    }
	    nb = aa ;
	}
	return nb * exp(ent) - kap_pr * pow(nb, gam_pr) ;
    }
    else {
        return 0 ;
    }

}

// Pressure from enthalpy
//------------------------

double Eos_multi_poly::press_ent_p(double ent, const Param* ) const {

    if ( ent > double(0) ) {

        int i = 0 ;

	if ( ent >= ent_crit_p[ndom_press-2] ) {
	    i = ndom_press - 1 ;
	}
	else {
	    while ( ent > ent_crit_p[i] ) {
	        i++ ;
	    }
	}
	assert(i <= ndom_press-1) ;

	double kap_pr = kap_p[i] ;
	double gam_pr = gam_p[i] ;

	int j = 0 ;

	if ( ent >= ent_crit_e[ndom_epsil-2] ) {
	    j = ndom_epsil - 1 ;
	}
	else {
	    while ( ent > ent_crit_e[j] ) {
	        j++ ;
	    }
	}
	assert(j <= ndom_epsil-1) ;

        double kap_ep = kap_e[j] ;
	double gam_ep = gam_e[j] ;

	double nb = 1. ;
	double nb_m1 = 0. ;

	// Switch the method of computation of the baryon density
	double diff_gamma = gam_ep - gam_pr + double(1) ;

	if (diff_gamma <= 0.) {
	    while (abs(1.-nb_m1/nb) > 1.e-15) {
	        nb_m1 = nb ;
		nb = pow( (exp(ent) - double(1))
			  / (kap_pr + kap_ep*pow(nb_m1, diff_gamma)),
			  double(1)/(gam_pr-double(1)) ) ;
	    }
	}
	else {
	    double aa = 0. ;
	    double xx = 1. ;
	    int m ;
	    double yy ;
	    double ent_value ;

	    while (xx > 1.e-15) {

	      ent_value = 1. ;   // Initialization
	      xx = 0.1 * xx ;
	      m = 0 ;

	      while (ent_value > 1.e-15) {

		m++ ;
		yy = aa + m * xx ;
		ent_value = exp(ent) - double(1) - kap_ep*pow(yy,gam_ep)
		  - kap_pr*pow(yy,gam_pr-double(1)) ;

	      }
	      aa += (m - 1) * xx ;
	    }
	    nb = aa ;
	}
	return kap_pr * pow(nb, gam_pr) ;
    }
    else {
        return 0 ;
    }

}

// dln(n)/dln(H) from enthalpy
//----------------------------

double Eos_multi_poly::der_nbar_ent_p(double ent, const Param* ) const {

    if ( ent > double(0) ) {

        if ( ent < 1.e-13 ) {
	    double gam_ep = gam_e[0] ;
	    return (double(1) + ent/double(2) + ent*ent/double(12)) / gam_ep ;
	}
	else {

	    int i = 0 ;

	    if ( ent >= ent_crit_p[ndom_press-2] ) {
	        i = ndom_press - 1 ;
	    }
	    else {
	        while ( ent > ent_crit_p[i] ) {
		    i++ ;
		}
	    }
	    assert(i <= ndom_press-1) ;

	    double kap_pr = kap_p[i] ;
	    double gam_pr = gam_p[i] ;

	    int j = 0 ;

	    if ( ent >= ent_crit_e[ndom_epsil-2] ) {
	        j = ndom_epsil - 1 ;
	    }
	    else {
	        while ( ent > ent_crit_e[j] ) {
		    j++ ;
		}
	    }
	    assert(j <= ndom_epsil-1) ;

	    double kap_ep = kap_e[j] ;
	    double gam_ep = gam_e[j] ;

	    double nb = 1. ;
	    double nb_m1 = 0. ;

	    // Switch the method of computation of the baryon density
	    double diff_gamma = gam_ep - gam_pr + double(1) ;

	    if (diff_gamma <= 0.) {
	        while (abs(1.-nb_m1/nb) > 1.e-15) {
		    nb_m1 = nb ;
		    nb = pow( (exp(ent) - double(1))
			      / (kap_pr + kap_ep*pow(nb_m1, diff_gamma)),
			      double(1)/(gam_pr-double(1)) ) ;
		}
	    }
	    else {
	        double aa = 0. ;
		double xx = 1. ;
		int m ;
		double yy ;
		double ent_value ;

		while (xx > 1.e-15) {

		    ent_value = 1. ;   // Initialization
		    xx = 0.1 * xx ;
		    m = 0 ;

		    while (ent_value > 1.e-15) {

		        m++ ;
			yy = aa + m * xx ;
			ent_value = exp(ent) - double(1)
			  - kap_ep*pow(yy,gam_ep)
			  - kap_pr*pow(yy,gam_pr-double(1)) ;

		    }
		    aa += (m - 1) * xx ;
		}
		nb = aa ;
	    }

	    return ent * exp(ent) / ( gam_ep * (exp(ent)-double(1))
				      + kap_pr * (gam_pr-gam_ep-double(1))
				      * pow(nb, gam_pr-double(1)) ) ;
	}
    }
    else {
        double gam_ep = gam_e[0] ;
        return double(1) / gam_ep ;	//  to ensure continuity at ent=0
    }

}

// dln(e)/dln(H) from enthalpy
//----------------------------

double Eos_multi_poly::der_ener_ent_p(double ent, const Param* ) const {

    if ( ent > double(0) ) {

        if ( ent < 1.e-13 ) {
	    double gam_ep = gam_e[0] ;
	    double kap_ep = kap_e[0] ;
	    double kap_pr = kap_p[0] ;
	    double nn = pow( (exp(ent)-double(1))/(kap_pr+kap_ep),
			     double(1)/gam_ep ) ;

	    return (double(1) + ent/double(2) + ent*ent/double(12)) / gam_ep
	      * ( double(1) + kap_ep*gam_ep*pow(nn, gam_ep)
		  / ( double(1) + kap_ep*pow(nn, gam_ep) ) ) ;
	}
	else {

	    int i = 0 ;

	    if ( ent >= ent_crit_p[ndom_press-2] ) {
	        i = ndom_press - 1 ;
	    }
	    else {
	        while ( ent > ent_crit_p[i] ) {
		    i++ ;
		}
	    }
	    assert(i <= ndom_press-1) ;

	    double kap_pr = kap_p[i] ;
	    double gam_pr = gam_p[i] ;

	    int j = 0 ;

	    if ( ent >= ent_crit_e[ndom_epsil-2] ) {
	        j = ndom_epsil - 1 ;
	    }
	    else {
	        while ( ent > ent_crit_e[j] ) {
		    j++ ;
		}
	    }
	    assert(j <= ndom_epsil-1) ;

	    double kap_ep = kap_e[j] ;
	    double gam_ep = gam_e[j] ;

	    double nb = 1. ;
	    double nb_m1 = 0. ;

	    // Switch the method of computation of the baryon density
	    double diff_gamma = gam_ep - gam_pr + double(1) ;

	    if (diff_gamma <= 0.) {
	        while (abs(1.-nb_m1/nb) > 1.e-15) {
		    nb_m1 = nb ;
		    nb = pow( (exp(ent) - double(1))
			      / (kap_pr + kap_ep*pow(nb_m1, diff_gamma)),
			      double(1)/(gam_pr-double(1)) ) ;
		}
	    }
	    else {
	        double aa = 0. ;
		double xx = 1. ;
		int m ;
		double yy ;
		double ent_value ;

		while (xx > 1.e-15) {

		    ent_value = 1. ;   // Initialization
		    xx = 0.1 * xx ;
		    m = 0 ;

		    while (ent_value > 1.e-15) {

		        m++ ;
			yy = aa + m * xx ;
			ent_value = exp(ent) - double(1)
			  - kap_ep*pow(yy,gam_ep)
			  - kap_pr*pow(yy,gam_pr-double(1)) ;

		    }
		    aa += (m - 1) * xx ;
		}
		nb = aa ;
	    }

	    return ent * exp(ent) / ( gam_ep * (exp(ent)-double(1))
				     + kap_pr * (gam_pr-gam_ep-double(1))
				     * pow(nb, gam_pr-double(1)) )
	      * ( double(1) + kap_ep*gam_ep*pow(nb, gam_ep)
		  / ( double(1) + kap_ep*pow(nb, gam_ep) ) ) ;
	}
    }
    else {
        double gam_ep = gam_e[0] ;
        return double(1) / gam_ep ;	//  to ensure continuity at ent=0
    }

}

// dln(p)/dln(H) from enthalpy
//----------------------------

double Eos_multi_poly::der_press_ent_p(double ent, const Param* ) const {

    if ( ent > double(0) ) {

        if ( ent < 1.e-13 ) {
	    double gam_ep = gam_e[0] ;
	    double gam_pr = gam_p[0] ;
	    return gam_pr * (double(1) + ent/double(2) + ent*ent/double(12))
	      / gam_ep ;
	}
	else {

	    int i = 0 ;

	    if ( ent >= ent_crit_p[ndom_press-2] ) {
	        i = ndom_press - 1 ;
	    }
	    else {
	        while ( ent > ent_crit_p[i] ) {
		    i++ ;
		}
	    }
	    assert(i <= ndom_press-1) ;

	    double kap_pr = kap_p[i] ;
	    double gam_pr = gam_p[i] ;

	    int j = 0 ;

	    if ( ent >= ent_crit_e[ndom_epsil-2] ) {
	        j = ndom_epsil - 1 ;
	    }
	    else {
	        while ( ent > ent_crit_e[j] ) {
		    j++ ;
		}
	    }
	    assert(j <= ndom_epsil-1) ;

	    double kap_ep = kap_e[j] ;
	    double gam_ep = gam_e[j] ;

	    double nb = 1. ;
	    double nb_m1 = 0. ;

	    // Switch the method of computation of the baryon density
	    double diff_gamma = gam_ep - gam_pr + double(1) ;

	    if (diff_gamma <= 0.) {
	        while (abs(1.-nb_m1/nb) > 1.e-15) {
		    nb_m1 = nb ;
		    nb = pow( (exp(ent) - double(1))
			      / (kap_pr + kap_ep*pow(nb_m1, diff_gamma)),
			      double(1)/(gam_pr-double(1)) ) ;
		}
	    }
	    else {
	        double aa = 0. ;
		double xx = 1. ;
		int m ;
		double yy ;
		double ent_value ;

		while (xx > 1.e-15) {

		    ent_value = 1. ;   // Initialization
		    xx = 0.1 * xx ;
		    m = 0 ;

		    while (ent_value > 1.e-15) {

		        m++ ;
			yy = aa + m * xx ;
			ent_value = exp(ent) - double(1)
			  - kap_ep*pow(yy,gam_ep)
			  - kap_pr*pow(yy,gam_pr-double(1)) ;

		    }
		    aa += (m - 1) * xx ;
		}
		nb = aa ;
	    }

	    return gam_pr * ent * exp(ent)
	      / ( gam_ep * (exp(ent)-double(1))
		  + kap_pr * (gam_pr-gam_ep-double(1))
		  * pow(nb, gam_pr-double(1)) ) ;
	}
    }
    else {
        double gam_ep = gam_e[0] ;
	double gam_pr = gam_p[0] ;
        return gam_pr / gam_ep ;	//  to ensure continuity at ent=0
    }

}
