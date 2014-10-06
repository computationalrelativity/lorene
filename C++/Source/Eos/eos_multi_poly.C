/*
 *  Methods of class Eos_multi_poly
 *
 *    (see file eos_multi_poly.h for documentation).
 *
 */

/*
 *   Copyright (c) 2009 Keisuke Taniguchi
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
 * Revision 1.8  2014/10/06 15:13:06  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.7  2012/10/09 16:15:26  j_novak
 * Corrected a bug in the constructor and save into a file
 *
 * Revision 1.6  2009/06/23 14:34:04  k_taniguchi
 * Completely revised.
 *
 * Revision 1.5  2004/06/23 15:42:08  e_gourgoulhon
 * Replaced all "abs" by "fabs".
 *
 * Revision 1.4  2004/05/13 07:38:57  k_taniguchi
 * Change the procedure for searching the baryon density from enthalpy.
 *
 * Revision 1.3  2004/05/09 10:43:52  k_taniguchi
 * Change the searching method of the baryon density again.
 *
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
#include <cstdlib>
#include <cstring>
#include <cmath>

// Lorene headers
#include "headcpp.h"
#include "eos_multi_poly.h"
#include "eos.h"
#include "utilitaires.h"
#include "unites.h"

double logp(double, double, double, double, double, double) ;
double dlpsdlh(double, double, double, double, double, double) ;
double dlpsdlnb(double, double, double, double, double) ;

//*********************************************************************

               //--------------------------------------//
               //              Constructors            //
               //--------------------------------------//

// Standard constructor
Eos_multi_poly::Eos_multi_poly(int npoly, double* gamma_i, double kappa0_i,
			       double logP1_i, double* logRho_i,
			       double* decInc_i)
    : Eos("Multi-polytropic EOS"),
      npeos(npoly), kappa0(kappa0_i), logP1(logP1_i), m0(double(1)) {

    assert(npeos > 1) ;

    gamma = new double [npeos] ;

    for (int l=0; l<npeos; l++) {
        gamma[l] = gamma_i[l] ;
    }

    logRho = new double [npeos-1] ;

    for (int l=0; l<npeos-1; l++) {
        logRho[l] = logRho_i[l] ;
    }

    decInc = new double [npeos-1] ;

    for (int l=0; l<npeos-1; l++) {
        decInc[l] = decInc_i[l] ;
    }

    set_auxiliary() ;

}


// Copy constructor
Eos_multi_poly::Eos_multi_poly(const Eos_multi_poly& eosmp)
    : Eos(eosmp),
      npeos(eosmp.npeos), kappa0(eosmp.kappa0), logP1(eosmp.logP1),
      m0(eosmp.m0) {

    gamma = new double [npeos] ;

    for (int l=0; l<npeos; l++) {
        gamma[l] = eosmp.gamma[l] ;
    }

    logRho = new double [npeos-1] ;

    for (int l=0; l<npeos-1; l++) {
        logRho[l] = eosmp.logRho[l] ;
    }

    kappa = new double [npeos] ;

    for (int l=0; l<npeos; l++) {
        kappa[l] = eosmp.kappa[l] ;
    }

    nbCrit = new double [npeos-1] ;

    for (int l=0; l<npeos-1; l++) {
        nbCrit[l] = eosmp.nbCrit[l] ;
    }

    entCrit = new double [npeos-1] ;

    for (int l=0; l<npeos-1; l++) {
        entCrit[l] = eosmp.entCrit[l] ;
    }

    decInc = new double [npeos-1] ;

    for (int l=0; l<npeos-1; l++) {
        decInc[l] = eosmp.decInc[l] ;
    }

    mu0 = new double [npeos] ;

    for (int l=0; l<npeos; l++) {
        mu0[l] = eosmp.mu0[l] ;
    }

}

//  Constructor from a binary file
Eos_multi_poly::Eos_multi_poly(FILE* fich) : Eos(fich) {

    fread_be(&npeos, sizeof(int), 1, fich) ;

    gamma = new double [npeos] ;

    for (int l=0; l<npeos; l++) {
        fread_be(&gamma[l], sizeof(double), 1, fich) ;
    }

    fread_be(&kappa0, sizeof(double), 1, fich) ;
    fread_be(&logP1, sizeof(double), 1, fich) ;

    logRho = new double [npeos-1] ;

    for (int l=0; l<npeos-1; l++) {
        fread_be(&logRho[l], sizeof(double), 1, fich) ;
    }

    decInc = new double [npeos-1] ;

    for (int l=0; l<npeos-1; l++) {
        fread_be(&decInc[l], sizeof(double), 1, fich) ;
    }

    m0 = double(1) ;

    set_auxiliary() ;

}

//  Constructor from a formatted file
Eos_multi_poly::Eos_multi_poly(ifstream& fich) : Eos(fich) {

    char blabla[80] ;

    fich >> npeos ; fich.getline(blabla, 80) ;

    gamma = new double [npeos] ;

    for (int l=0; l<npeos; l++) {
        fich >> gamma[l] ; fich.getline(blabla, 80) ;
    }

    fich >> kappa0 ; fich.getline(blabla, 80) ;
    fich >> logP1 ; fich.getline(blabla, 80) ;

    logRho = new double [npeos-1] ;

    for (int l=0; l<npeos-1; l++) {
        fich >> logRho[l] ; fich.getline(blabla, 80) ;
    }

    decInc = new double [npeos-1] ;

    for (int l=0; l<npeos-1; l++) {
        fich >> decInc[l] ; fich.getline(blabla, 80) ;
    }

    m0 = double(1) ;

    set_auxiliary() ;

}

// Destructor
Eos_multi_poly::~Eos_multi_poly() {

    delete [] gamma ;
    delete [] logRho ;
    delete [] kappa ;
    delete [] nbCrit ;
    delete [] entCrit ;
    delete [] decInc ;
    delete [] mu0 ;

}

			//--------------//
			//  Assignment  //
			//--------------//

void Eos_multi_poly::operator=(const Eos_multi_poly& ) {

    cout << "Eos_multi_poly::operator=  : not implemented yet !" << endl ;
    abort() ;

}


		     //-----------------------//
		     //	    Miscellaneous     //
		     //-----------------------//

void Eos_multi_poly::set_auxiliary() {

    using namespace Unites ;

    double* kappa_cgs = new double [npeos] ;

    kappa_cgs[0] = kappa0 ;

    kappa_cgs[1] = pow(10., logP1-logRho[0]*gamma[1]) ;

    if (npeos > 2) {

        kappa_cgs[2] = kappa_cgs[1]
	  * pow(10., logRho[0]*(gamma[1]-gamma[2])) ;

	if (npeos > 3) {

	    for (int l=3; l<npeos; l++) {

	        kappa_cgs[l] = kappa_cgs[l-1]
		  * pow(10., logRho[l-2]*(gamma[l-1]-gamma[l])) ;

	    }

	}

    }

    kappa = new double [npeos] ;

    double rhonuc_cgs = rhonuc_si * 1.e-3 ;

    for (int l=0; l<npeos; l++) {
        kappa[l] = kappa_cgs[l] * pow( rhonuc_cgs, gamma[l] - double(1) ) ;
	// Conversion from cgs units to Lorene units
    }

    delete [] kappa_cgs ;

    mu0 = new double [npeos] ;
    mu0[0] = double(1) ;  // We define

    entCrit = new double [npeos-1] ;

    nbCrit = new double [npeos-1] ;

    for (int l=0; l<npeos-1; l++) {

        nbCrit[l] =
	  pow(kappa[l]/kappa[l+1], double(1)/(gamma[l+1]-gamma[l])) ;

	mu0[l+1] = mu0[l]
	  + ( kappa[l] * pow(nbCrit[l], gamma[l]-double(1))
	      / (gamma[l]-double(1))
	      - kappa[l+1] * pow(nbCrit[l], gamma[l+1]-double(1))
	      / (gamma[l+1]-double(1)) ) ;

	entCrit[l] = log ( mu0[l] / m0
			   + kappa[l] * gamma[l]
			   * pow(nbCrit[l], gamma[l]-double(1))
			   / (gamma[l]-double(1)) / m0 ) ;

    }

}


               //------------------------------//
               //     Comparison operators     //
               //------------------------------//

bool Eos_multi_poly::operator==(const Eos& eos_i) const {
    
    bool resu = true ; 
    
    if ( eos_i.identify() != identify() ) {
	cout << "The second EOS is not of type Eos_multi_poly !" << endl ; 
	resu = false ; 
    }
    else{
	
	const Eos_multi_poly& eos
	  = dynamic_cast<const Eos_multi_poly&>(eos_i) ; 

	if (eos.get_npeos() != npeos) {
	    cout << "The two Eos_multi_poly have "
		 << "different number of polytropes : "
		 << npeos << " <-> " << eos.get_npeos() << endl ; 
	    resu = false ; 
	}

	for (int l=0; l<npeos; l++) {
	    if (eos.get_gamma(l) != gamma[l]) {
	        cout << "The two Eos_multi_poly have different gamma "
		     << gamma[l] << " <-> " 
		     << eos.get_gamma(l) << endl ;
	    resu = false ;
	    }
	}

	for (int l=0; l<npeos; l++) {
	    if (eos.get_kappa(l) != kappa[l]) {
	        cout << "The two Eos_multi_poly have different kappa "
		     << kappa[l] << " <-> " 
		     << eos.get_kappa(l) << endl ;
	    resu = false ;
	    }
	}

    }
    
    return resu ; 
    
}

bool Eos_multi_poly::operator!=(const Eos& eos_i) const {
 
    return !(operator==(eos_i)) ; 
       
}

                     //--------------------------//
                     //         Outputs          //
                     //--------------------------//

void Eos_multi_poly::sauve(FILE* fich) const {

    Eos::sauve(fich) ;

    fwrite_be(&npeos, sizeof(int), 1, fich) ;

    for (int l=0; l<npeos; l++) {
        fwrite_be(&gamma[l], sizeof(double), 1, fich) ;
    }

    fwrite_be(&kappa0, sizeof(double), 1, fich) ;
    fwrite_be(&logP1, sizeof(double), 1, fich) ;

    for (int l=0; l<npeos-1; l++) {
        fwrite_be(&logRho[l], sizeof(double), 1, fich) ;
    }

    for (int l=0; l<npeos-1; l++) {
        fwrite_be(&decInc[l], sizeof(double), 1, fich) ;
    }

}


ostream& Eos_multi_poly::operator>>(ostream & ost) const {

    using namespace Unites ;

    ost << "EOS of class Eos_multi_poly "
	<< "(multiple polytropic equation of state) : " << endl ;

    ost << "  Number of polytropes : "
	<< npeos << endl << endl ;

    double rhonuc_cgs = rhonuc_si * 1.e-3 ;

    ost.precision(16) ;
    for (int l=0; l<npeos; l++) {
        ost << "  EOS in region " << l << " : " << endl ;
	ost << "  ---------------" << endl ;
	ost << "    gamma  : " << gamma[l] << endl ;
	ost << "    kappa  : " << kappa[l]
	    << " [Lorene units: rho_nuc c^2 / n_nuc^gamma]" << endl ;

	double kappa_cgs = kappa[l]
	  * pow( rhonuc_cgs, double(1) - gamma[l] ) ;

	ost << "           : " << kappa_cgs
	    << " [(g/cm^3)^{1-gamma}]" << endl ;
    }

    ost << endl ;
    ost << "  Exponent of the pressure at the fiducial density rho_1"
	<< endl ;
    ost << "  ------------------------------------------------------"
	<< endl ;
    ost << "    log P1 : " << logP1 << endl ;

    ost << endl ;
    ost << "  Exponent of fiducial densities" << endl ;
    ost << "  ------------------------------" << endl ;
    for (int l=0; l<npeos-1; l++) {
      ost << "    log rho[" << l << "] : " << logRho[l] << endl ;
    }

    ost << endl ;
    for (int l=0; l<npeos-1; l++) {
        ost << "  Critical density and enthalpy between domains "
	    << l << " and " << l+1 << " : " << endl ;
	ost << "  -----------------------------------------------------"
	    << endl ;
	ost << "    num. dens. : " << nbCrit[l] << " [Lorene units: n_nuc]"
	    << endl ;
	ost << "    density :    " << nbCrit[l] * rhonuc_cgs << " [g/cm^3]"
	    << endl ;

	ost << "    ln(ent) :    " << entCrit[l] << endl ;
    }

    ost << endl ;
    for (int l=0; l<npeos; l++) {
        ost << "  Relat. chem. pot. in region " << l << " : " << endl ;
	ost << "  -----------------------------" << endl ;
	ost << "    mu : " << mu0[l] << " [m_B c^2]" << endl ;
    }

    return ost ;

}


               //------------------------------------------//
               //          Computational routines          //
	       //------------------------------------------//

// Baryon rest-mass density from enthalpy
//----------------------------------------

double Eos_multi_poly::nbar_ent_p(double ent, const Param* ) const {

    int i = 0 ; // "i" corresponds to the number of domain
                // i=0: gamma[0], kappa[0], i=1: gamma[1], kappa[1], .....
                // The buffer zone is included in the next zone.
    for (int l=0; l<npeos-1; l++) {
        if ( ent > entCrit[l]*(double(1)-decInc[l]) ) {
	    i++ ;
	}
    }

    double mgam1skapgam = 1. ;  // Initialization
    if (i == 0) {

        if ( ent > double(0) ) {

	    mgam1skapgam = m0*(gamma[0]-double(1))/kappa[0]/gamma[0] ;

	    return pow( mgam1skapgam*(exp(ent)-double(1)),
			double(1)/(gamma[0]-double(1)) ) ; // mu0[0]/m0=1

	}
	else {
	    return double(0) ;
	}

    }
    else {

        double entSmall = entCrit[i-1]*(double(1)-decInc[i-1]) ;
        double entLarge = entCrit[i-1]*(double(1)+decInc[i-1]) ;

        if ( ent < entLarge ) {

	    double log10H = log10( ent ) ;
	    double log10HSmall = log10( entSmall ) ;
	    double log10HLarge = log10( entLarge ) ;
	    double dH = log10HLarge - log10HSmall ;
	    double uu = (log10H - log10HSmall) / dH ;

	    double mgam1skapgamSmall = m0*(gamma[i-1]-double(1))
	      /kappa[i-1]/gamma[i-1] ;
	    double mgam1skapgamLarge = m0*(gamma[i]-double(1))
	      /kappa[i]/gamma[i] ;

	    double nnSmall = pow( mgam1skapgamSmall
				  *(exp(entSmall)-mu0[i-1]/m0),
				  double(1)/(gamma[i-1]-double(1)) ) ;
	    double nnLarge = pow( mgam1skapgamLarge
				  *(exp(entLarge)-mu0[i]/m0),
				  double(1)/(gamma[i]-double(1)) ) ;

	    double ppSmall = kappa[i-1] * pow( nnSmall, gamma[i-1] ) ;
	    double ppLarge = kappa[i] * pow( nnLarge, gamma[i] ) ;

	    double eeSmall = mu0[i-1] * nnSmall
	      + ppSmall / (gamma[i-1] - double(1)) ;
	    double eeLarge = mu0[i] * nnLarge
	      + ppLarge / (gamma[i] - double(1)) ;

	    double log10PSmall = log10( ppSmall ) ;
	    double log10PLarge = log10( ppLarge ) ;

	    double dlpsdlhSmall = entSmall
	      * (double(1) + eeSmall / ppSmall) ;
	    double dlpsdlhLarge = entLarge
	      * (double(1) + eeLarge / ppLarge) ;

	    double log10PInterpolate = logp(log10PSmall, log10PLarge,
					    dlpsdlhSmall, dlpsdlhLarge,
					    dH, uu) ;

	    double dlpsdlhInterpolate = dlpsdlh(log10PSmall, log10PLarge,
						dlpsdlhSmall, dlpsdlhLarge,
						dH, uu) ;

	    double pp = pow(double(10), log10PInterpolate) ;

	    return pp / ent * dlpsdlhInterpolate * exp(-ent) / m0 ;
	    // Is m0 necessary?

	}
	else {

	    mgam1skapgam = m0*(gamma[i]-double(1))/kappa[i]/gamma[i] ;

	    return pow( mgam1skapgam*(exp(ent)-mu0[i]/m0),
			double(1)/(gamma[i]-double(1)) ) ;

	}

    }

}

// Energy density from enthalpy
//------------------------------

double Eos_multi_poly::ener_ent_p(double ent, const Param* ) const {

    int i = 0 ; // "i" corresponds to the number of domain
                // i=0: gamma[0], kappa[0], i=1: gamma[1], kappa[1], .....
                // The buffer zone is included in the next zone.
    for (int l=0; l<npeos-1; l++) {
        if ( ent > entCrit[l]*(double(1)-decInc[l]) ) {
	    i++ ;
	}
    }

    double mgam1skapgam = 1. ;  // Initialization
    double nn = 0. ;  // Initialization
    double pp = 0. ;  // Initialization

    if (i == 0) {

        if ( ent > double(0) ) {

	    mgam1skapgam = m0*(gamma[0]-double(1))/kappa[0]/gamma[0] ;

	    nn = pow( mgam1skapgam*(exp(ent)-double(1)),
		      double(1)/(gamma[0]-double(1)) ) ; // mu0[0]/m0=1

	    pp = kappa[0] * pow( nn, gamma[0] ) ;

	    return mu0[0] * nn + pp / (gamma[0] - double(1)) ;

	}
	else {
	    return double(0) ;
	}

    }
    else {

        double entSmall = entCrit[i-1]*(double(1)-decInc[i-1]) ;
        double entLarge = entCrit[i-1]*(double(1)+decInc[i-1]) ;

        if ( ent < entLarge ) {

	    double log10H = log10( ent ) ;
	    double log10HSmall = log10( entSmall ) ;
	    double log10HLarge = log10( entLarge ) ;
	    double dH = log10HLarge - log10HSmall ;
	    double uu = (log10H - log10HSmall) / dH ;

	    double mgam1skapgamSmall = m0*(gamma[i-1]-double(1))
	      /kappa[i-1]/gamma[i-1] ;
	    double mgam1skapgamLarge = m0*(gamma[i]-double(1))
	      /kappa[i]/gamma[i] ;

	    double nnSmall = pow( mgam1skapgamSmall
				  *(exp(entSmall)-mu0[i-1]/m0),
				  double(1)/(gamma[i-1]-double(1)) ) ;
	    double nnLarge = pow( mgam1skapgamLarge
				  *(exp(entLarge)-mu0[i]/m0),
				  double(1)/(gamma[i]-double(1)) ) ;

	    double ppSmall = kappa[i-1] * pow( nnSmall, gamma[i-1] ) ;
	    double ppLarge = kappa[i] * pow( nnLarge, gamma[i] ) ;

	    double eeSmall = mu0[i-1] * nnSmall
	      + ppSmall / (gamma[i-1] - double(1)) ;
	    double eeLarge = mu0[i] * nnLarge
	      + ppLarge / (gamma[i] - double(1)) ;

	    double log10PSmall = log10( ppSmall ) ;
	    double log10PLarge = log10( ppLarge ) ;

	    double dlpsdlhSmall = entSmall
	      * (double(1) + eeSmall / ppSmall) ;
	    double dlpsdlhLarge = entLarge
	      * (double(1) + eeLarge / ppLarge) ;

	    double log10PInterpolate = logp(log10PSmall, log10PLarge,
					    dlpsdlhSmall, dlpsdlhLarge,
					    dH, uu) ;

	    double dlpsdlhInterpolate = dlpsdlh(log10PSmall, log10PLarge,
						dlpsdlhSmall, dlpsdlhLarge,
						dH, uu) ;

	    pp = pow(double(10), log10PInterpolate) ;

	    return pp / ent * dlpsdlhInterpolate - pp ;

	}
	else {

	    mgam1skapgam = m0*(gamma[i]-double(1))/kappa[i]/gamma[i] ;

	    nn = pow( mgam1skapgam*(exp(ent)-mu0[i]/m0),
		      double(1)/(gamma[i]-double(1)) ) ;

	    pp = kappa[i] * pow( nn, gamma[i] ) ;

	    return mu0[i] * nn + pp / (gamma[i] - double(1)) ;

	}

    }

}


// Pressure from enthalpy
//------------------------

double Eos_multi_poly::press_ent_p(double ent, const Param* ) const {

    int i = 0 ; // "i" corresponds to the number of domain
                // i=0: gamma[0], kappa[0], i=1: gamma[1], kappa[1], .....
                // The buffer zone is included in the next zone.
    for (int l=0; l<npeos-1; l++) {
        if ( ent > entCrit[l]*(double(1)-decInc[l]) ) {
	    i++ ;
	}
    }

    double mgam1skapgam = 1. ;  // Initialization
    double nn = 0. ;  // Initialization

    if (i == 0) {

        if ( ent > double(0) ) {

	    mgam1skapgam = m0*(gamma[0]-double(1))/kappa[0]/gamma[0] ;

	    nn = pow( mgam1skapgam*(exp(ent)-double(1)),
		      double(1)/(gamma[0]-double(1)) ) ; // mu0[0]/m0=1

	    return kappa[0] * pow( nn, gamma[0] ) ;

	}
	else {
	    return double(0) ;
	}

    }
    else {

        double entSmall = entCrit[i-1]*(double(1)-decInc[i-1]) ;
        double entLarge = entCrit[i-1]*(double(1)+decInc[i-1]) ;

        if ( ent < entLarge ) {

	    double log10H = log10( ent ) ;
	    double log10HSmall = log10( entSmall ) ;
	    double log10HLarge = log10( entLarge ) ;
	    double dH = log10HLarge - log10HSmall ;
	    double uu = (log10H - log10HSmall) / dH ;

	    double mgam1skapgamSmall = m0*(gamma[i-1]-double(1))
	      /kappa[i-1]/gamma[i-1] ;
	    double mgam1skapgamLarge = m0*(gamma[i]-double(1))
	      /kappa[i]/gamma[i] ;

	    double nnSmall = pow( mgam1skapgamSmall
				  *(exp(entSmall)-mu0[i-1]/m0),
				  double(1)/(gamma[i-1]-double(1)) ) ;
	    double nnLarge = pow( mgam1skapgamLarge
				  *(exp(entLarge)-mu0[i]/m0),
				  double(1)/(gamma[i]-double(1)) ) ;

	    double ppSmall = kappa[i-1] * pow( nnSmall, gamma[i-1] ) ;
	    double ppLarge = kappa[i] * pow( nnLarge, gamma[i] ) ;

	    double eeSmall = mu0[i-1] * nnSmall
	      + ppSmall / (gamma[i-1] - double(1)) ;
	    double eeLarge = mu0[i] * nnLarge
	      + ppLarge / (gamma[i] - double(1)) ;

	    double log10PSmall = log10( ppSmall ) ;
	    double log10PLarge = log10( ppLarge ) ;

	    double dlpsdlhSmall = entSmall
	      * (double(1) + eeSmall / ppSmall) ;
	    double dlpsdlhLarge = entLarge
	      * (double(1) + eeLarge / ppLarge) ;

	    double log10PInterpolate = logp(log10PSmall, log10PLarge,
					    dlpsdlhSmall, dlpsdlhLarge,
					    dH, uu) ;

	    return pow(double(10), log10PInterpolate) ;

	}
	else {

	    mgam1skapgam = m0*(gamma[i]-double(1))/kappa[i]/gamma[i] ;

	    nn = pow( mgam1skapgam*(exp(ent)-mu0[i]/m0),
		      double(1)/(gamma[i]-double(1)) ) ;

	    return kappa[i] * pow( nn, gamma[i] ) ;

	}

    }

}


// dln(n)/dln(H) from enthalpy
//----------------------------

double Eos_multi_poly::der_nbar_ent_p(double ent, const Param* ) const {

    int i = 0 ; // "i" corresponds to the number of domain
                // i=0: gamma[0], kappa[0], i=1: gamma[1], kappa[1], .....
                // The buffer zone is included in the next zone.
    for (int l=0; l<npeos-1; l++) {
        if ( ent > entCrit[l]*(double(1)-decInc[l]) ) {
	    i++ ;
	}
    }

    if (i == 0) {

        if ( ent > double(0) ) {

	    if ( ent < 1.e-13 ) {

	        return ( double(1) + ent/double(2) + ent*ent/double(12) )
		  / (gamma[0] - double(1)) ;

	    }
	    else {

	        return ent / (double(1) - exp(-ent))
		  / (gamma[0] - double(1)) ; // mu0[0]/m0=1

	    }

	}
	else {

	  return double(1) / (gamma[0] - double(1)) ;

	}

    }
    else {

        if ( ent < entCrit[i-1]*(double(1)+decInc[i-1]) ) {

	    double zeta = der_press_ent_p(ent) / der_press_nbar_p(ent) ;

	    return zeta ;

	}
	else {

	    return ent / (double(1) - (mu0[i]/m0) * exp(-ent))
	      / (gamma[i] - double(1)) ;

	}

    }
}


// dln(e)/dln(H) from enthalpy
//----------------------------

double Eos_multi_poly::der_ener_ent_p(double ent, const Param* ) const {

    int i = 0 ; // "i" corresponds to the number of domain
                // i=0: gamma[0], kappa[0], i=1: gamma[1], kappa[1], .....
                // The buffer zone is included in the next zone.
    for (int l=0; l<npeos-1; l++) {
        if ( ent > entCrit[l]*(double(1)-decInc[l]) ) {
	    i++ ;
	}
    }

    double mgam1skapgam = 1. ;  // Initialization
    double nn = 0. ;  // Initialization
    double pp = 0. ;  // Initialization
    double ee = 0. ;  // Initialization

    if (i == 0) {

        if ( ent > double(0) ) {

	    mgam1skapgam = m0*(gamma[0]-double(1))/kappa[0]/gamma[0] ;

	    nn = pow( mgam1skapgam*(exp(ent)-double(1)),
		      double(1)/(gamma[0]-double(1)) ) ; // mu0[0]/m0=1

	    pp = kappa[0] * pow( nn, gamma[0] ) ;

	    ee = mu0[0] * nn + pp / (gamma[0] - double(1)) ;

	    if ( ent < 1.e-13 ) {

	        return ( double(1) + ent/double(2) + ent*ent/double(12) )
		  / (gamma[0] - double(1)) * (double(1) + pp / ee) ;

	    }
	    else {

	        return ent / (double(1) - exp(-ent))
		  / (gamma[0] - double(1)) * (double(1) + pp / ee) ;
		// mu0[0]/m0=1

	    }

	}
	else {

	  return double(1) / (gamma[0] - double(1)) ;

	}

    }
    else {

        double entSmall = entCrit[i-1]*(double(1)-decInc[i-1]) ;
        double entLarge = entCrit[i-1]*(double(1)+decInc[i-1]) ;

        if ( ent < entLarge ) {

	    double log10H = log10( ent ) ;
	    double log10HSmall = log10( entSmall ) ;
	    double log10HLarge = log10( entLarge ) ;
	    double dH = log10HLarge - log10HSmall ;
	    double uu = (log10H - log10HSmall) / dH ;

	    double mgam1skapgamSmall = m0*(gamma[i-1]-double(1))
	      /kappa[i-1]/gamma[i-1] ;
	    double mgam1skapgamLarge = m0*(gamma[i]-double(1))
	      /kappa[i]/gamma[i] ;

	    double nnSmall = pow( mgam1skapgamSmall
				  *(exp(entSmall)-mu0[i-1]/m0),
				  double(1)/(gamma[i-1]-double(1)) ) ;
	    double nnLarge = pow( mgam1skapgamLarge
				  *(exp(entLarge)-mu0[i]/m0),
				  double(1)/(gamma[i]-double(1)) ) ;

	    double ppSmall = kappa[i-1] * pow( nnSmall, gamma[i-1] ) ;
	    double ppLarge = kappa[i] * pow( nnLarge, gamma[i] ) ;

	    double eeSmall = mu0[i-1] * nnSmall
	      + ppSmall / (gamma[i-1] - double(1)) ;
	    double eeLarge = mu0[i] * nnLarge
	      + ppLarge / (gamma[i] - double(1)) ;

	    double log10PSmall = log10( ppSmall ) ;
	    double log10PLarge = log10( ppLarge ) ;

	    double dlpsdlhSmall = entSmall
	      * (double(1) + eeSmall / ppSmall) ;
	    double dlpsdlhLarge = entLarge
	      * (double(1) + eeLarge / ppLarge) ;

	    double log10PInterpolate = logp(log10PSmall, log10PLarge,
					    dlpsdlhSmall, dlpsdlhLarge,
					    dH, uu) ;

	    double dlpsdlhInterpolate = dlpsdlh(log10PSmall, log10PLarge,
						dlpsdlhSmall, dlpsdlhLarge,
						dH, uu) ;

	    pp = pow(double(10), log10PInterpolate) ;

	    ee = pp / ent * dlpsdlhInterpolate - pp ;

	    double zeta = (double(1) + pp / ee) * der_press_ent_p(ent)
	      / der_press_nbar_p(ent) ;

	    return zeta ;

	}
	else {

	    mgam1skapgam = m0*(gamma[i]-double(1))/kappa[i]/gamma[i] ;

	    nn = pow( mgam1skapgam*(exp(ent)-mu0[i]/m0),
		      double(1)/(gamma[i]-double(1)) ) ;

	    pp = kappa[i] * pow( nn, gamma[i] ) ;

	    ee = mu0[i] * nn + pp / (gamma[i] - double(1)) ;

	    return ent / (double(1) - (mu0[i]/m0) * exp(-ent))
	      / (gamma[i] - double(1)) * (double(1) + pp / ee) ;

	}

    }
}


// dln(p)/dln(H) from enthalpy
//----------------------------

double Eos_multi_poly::der_press_ent_p(double ent, const Param* ) const {

    int i = 0 ; // "i" corresponds to the number of domain
                // i=0: gamma[0], kappa[0], i=1: gamma[1], kappa[1], .....
                // The buffer zone is included in the next zone.
    for (int l=0; l<npeos-1; l++) {
        if ( ent > entCrit[l]*(double(1)-decInc[l]) ) {
	    i++ ;
	}
    }

    if (i == 0) {

        if ( ent > double(0) ) {

	    if ( ent < 1.e-13 ) {

	        return gamma[0]
		  * ( double(1) + ent/double(2) + ent*ent/double(12) )
		  / (gamma[0] - double(1)) ;

	    }
	    else {

	        return gamma[0] * ent / (double(1) - exp(-ent))
		  / (gamma[0] - double(1)) ; // mu0[0]/m0=1

	    }

	}
	else {

	  return gamma[0] / (gamma[0] - double(1)) ;

	}

    }
    else {

        double entSmall = entCrit[i-1]*(double(1)-decInc[i-1]) ;
        double entLarge = entCrit[i-1]*(double(1)+decInc[i-1]) ;

        if ( ent < entLarge ) {

	    double log10H = log10( ent ) ;
	    double log10HSmall = log10( entSmall ) ;
	    double log10HLarge = log10( entLarge ) ;
	    double dH = log10HLarge - log10HSmall ;
	    double uu = (log10H - log10HSmall) / dH ;

	    double mgam1skapgamSmall = m0*(gamma[i-1]-double(1))
	      /kappa[i-1]/gamma[i-1] ;
	    double mgam1skapgamLarge = m0*(gamma[i]-double(1))
	      /kappa[i]/gamma[i] ;

	    double nnSmall = pow( mgam1skapgamSmall
				  *(exp(entSmall)-mu0[i-1]/m0),
				  double(1)/(gamma[i-1]-double(1)) ) ;
	    double nnLarge = pow( mgam1skapgamLarge
				  *(exp(entLarge)-mu0[i]/m0),
				  double(1)/(gamma[i]-double(1)) ) ;

	    double ppSmall = kappa[i-1] * pow( nnSmall, gamma[i-1] ) ;
	    double ppLarge = kappa[i] * pow( nnLarge, gamma[i] ) ;

	    double eeSmall = mu0[i-1] * nnSmall
	      + ppSmall / (gamma[i-1] - double(1)) ;
	    double eeLarge = mu0[i] * nnLarge
	      + ppLarge / (gamma[i] - double(1)) ;

	    double log10PSmall = log10( ppSmall ) ;
	    double log10PLarge = log10( ppLarge ) ;

	    double dlpsdlhSmall = entSmall
	      * (double(1) + eeSmall / ppSmall) ;
	    double dlpsdlhLarge = entLarge
	      * (double(1) + eeLarge / ppLarge) ;

	    double dlpsdlhInterpolate = dlpsdlh(log10PSmall, log10PLarge,
						dlpsdlhSmall, dlpsdlhLarge,
						dH, uu) ;
	    return dlpsdlhInterpolate ;

	}
	else {

	    return gamma[i] * ent / (double(1) - (mu0[i]/m0) * exp(-ent))
	      / (gamma[i] - double(1)) ;

	}

    }
}


// dln(p)/dln(n) from enthalpy
//----------------------------

double Eos_multi_poly::der_press_nbar_p(double ent, const Param* ) const {

    int i = 0 ; // "i" corresponds to the number of domain
                // i=0: gamma[0], kappa[0], i=1: gamma[1], kappa[1], .....
                // The buffer zone is included in the next zone.
    for (int l=0; l<npeos-1; l++) {
        if ( ent > entCrit[l]*(double(1)-decInc[l]) ) {
	    i++ ;
	}
    }

    if (i == 0) {

        return gamma[0] ;

    }
    else {

        double entSmall = entCrit[i-1]*(double(1)-decInc[i-1]) ;
        double entLarge = entCrit[i-1]*(double(1)+decInc[i-1]) ;

        if ( ent < entLarge ) {

	    double log10H = log10( ent ) ;
	    double log10HSmall = log10( entSmall ) ;
	    double log10HLarge = log10( entLarge ) ;

	    double dlpsdlnbInterpolate = dlpsdlnb(log10HSmall, log10HLarge,
						  gamma[i-1], gamma[i],
						  log10H) ;

	    return dlpsdlnbInterpolate ;

	}
	else {

	    return gamma[i] ;

	}

    }
}


//***************************************************//
//     Functions which appear in the computation     //
//***************************************************//

double logp(double log10PressSmall, double log10PressLarge,
	    double dlpsdlhSmall, double dlpsdlhLarge,
	    double dx, double u) {

    double resu = log10PressSmall * (double(2) * u * u * u
				     - double(3) * u * u + double(1))
      + log10PressLarge * (double(3) * u * u - double(2) * u * u * u)
      + dlpsdlhSmall * dx * (u * u * u - double(2) * u * u + u)
      - dlpsdlhLarge * dx * (u * u - u * u * u) ;

    return resu ;

}

double dlpsdlh(double log10PressSmall, double log10PressLarge,
	       double dlpsdlhSmall, double dlpsdlhLarge,
	       double dx, double u) {

    double resu = double(6) * (log10PressSmall - log10PressLarge)
      * (u * u - u) / dx
      + dlpsdlhSmall * (double(3) * u * u - double(4) * u + double(1))
      + dlpsdlhLarge * (double(3) * u * u - double(2) * u) ;

    return resu ;

}

double dlpsdlnb(double log10HSmall, double log10HLarge,
		double dlpsdlnbSmall, double dlpsdlnbLarge,
		double log10H) {

    double resu = log10H * (dlpsdlnbSmall - dlpsdlnbLarge)
      / (log10HSmall - log10HLarge)
      + (log10HSmall * dlpsdlnbLarge - log10HLarge * dlpsdlnbSmall)
      / (log10HSmall - log10HLarge) ;

    return resu ;

}
