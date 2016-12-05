/*
 * Methods for the class Eos_strange_cr_cr
 *
 * (see file eos.h for documentation)
 *
 */

/*
 *   Copyright (c) 2000 J. Leszek Zdunik
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
 * Revision 1.8  2016/12/05 16:17:51  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.7  2014/10/13 08:52:54  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.6  2014/10/06 15:13:07  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.5  2004/03/25 10:29:02  j_novak
 * All LORENE's units are now defined in the namespace Unites (in file unites.h).
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
 * Revision 2.2  2001/02/07  09:49:26  eric
 * Suppression de la fonction derent_ent_p.
 * Ajout des fonctions donnant les derivees de l'EOS:
 *      der_nbar_ent_p
 *      der_ener_ent_p
 *      der_press_ent_p
 * /
 *
 * Revision 2.1  2000/11/24  13:26:50  eric
 * Correction erreur dans set_auxillary(): rho_nd_nucl est un membre !
 *
 * Revision 2.0  2000/11/23  14:46:35  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */


// Headers C
#include <cstdlib>
#include <cstring>
#include <cmath>

// Headers Lorene
#include "eos.h"
#include "cmp.h"
#include "utilitaires.h"
#include "unites.h"

		    //------------------------------------//
		    //		Constructors		  //
		    //------------------------------------//

// Standard constructor
// --------------------
namespace Lorene {
Eos_strange_cr::Eos_strange_cr(double n0_b60_i, double b60_i, double ent0_i,
		    double eps_fit_i, double rho0_b60_i,
		    double ent_nd_i, double rho_nd_i, double gam_i) :
	     Eos("Strange matter EOS from Zdunik (2000) with crust"),
	     n0_b60(n0_b60_i), 
	     b60(b60_i), 
	     ent0(ent0_i), 
	     eps_fit(eps_fit_i), 
	     rho0_b60(rho0_b60_i),
	     ent_nd(ent_nd_i),
	     rho_nd(rho_nd_i),
	     gam(gam_i) {
	
    set_auxiliary() ; 
		    
}		    

// Copy constructor
// -----------------

Eos_strange_cr::Eos_strange_cr(const Eos_strange_cr& eos_i) :
	    Eos(eos_i), 
	    n0_b60(eos_i.n0_b60), 
	    b60(eos_i.b60), 
	    ent0(eos_i.ent0), 
	    eps_fit(eos_i.eps_fit), 
	    rho0_b60(eos_i.rho0_b60),
	    ent_nd(eos_i.ent_nd),
	    rho_nd(eos_i.rho_nd),
	    gam(eos_i.gam)		{

    set_auxiliary() ; 

}


// Constructor from binary file
// ----------------------------
Eos_strange_cr::Eos_strange_cr(FILE* fich) :
	Eos(fich) {
        
    fread_be(&n0_b60, sizeof(double), 1, fich) ;		
    fread_be(&b60, sizeof(double), 1, fich) ;		
    fread_be(&ent0, sizeof(double), 1, fich) ;		
    fread_be(&eps_fit, sizeof(double), 1, fich) ;		
    fread_be(&rho0_b60, sizeof(double), 1, fich) ;		
    fread_be(&ent_nd, sizeof(double), 1, fich) ;		
    fread_be(&rho_nd, sizeof(double), 1, fich) ;		
    fread_be(&gam, sizeof(double), 1, fich) ;		

    set_auxiliary() ; 

}

// Constructor from a formatted file
// ---------------------------------
Eos_strange_cr::Eos_strange_cr(ifstream& fich) :
	Eos(fich) {

    char blabla[80] ;

    fich >> n0_b60 ; fich.getline(blabla, 80) ;
    fich >> b60 ; fich.getline(blabla, 80) ;
    fich >> ent0 ; fich.getline(blabla, 80) ;
    fich >> eps_fit ; fich.getline(blabla, 80) ;
    fich >> rho0_b60 ; fich.getline(blabla, 80) ;
    fich >> ent_nd ; fich.getline(blabla, 80) ;
    fich >> rho_nd ; fich.getline(blabla, 80) ;
    fich >> gam ; fich.getline(blabla, 80) ;

    set_auxiliary() ;

}
			//--------------//
			//  Destructor  //
			//--------------//

Eos_strange_cr::~Eos_strange_cr(){

    // does nothing

}

			//--------------//
			//  Assignment  //
			//--------------//

void Eos_strange_cr::operator=(const Eos_strange_cr& eosi) {

    set_name(eosi.name) ;

    n0_b60 = eosi.n0_b60 ;
    b60 = eosi.b60 ;
    ent0 = eosi.ent0 ;
    eps_fit = eosi.eps_fit ;
    rho0_b60 = eosi.rho0_b60 ;
    ent_nd = eosi.ent_nd ;
    rho_nd = eosi.rho_nd ;
    gam = eosi.gam ;

    set_auxiliary() ;

}


		  //-----------------------//
		  //	Miscellaneous	   //
		  //-----------------------//

void Eos_strange_cr::set_auxiliary() {

  using namespace Unites ;

    rho_nd_nucl = rho_nd * mevpfm3 ;

    rho0 = b60 * rho0_b60 * mevpfm3 ;

    b34 = pow(b60, double(0.75)) ;

    n0 = b34 * n0_b60 * double(10) ; 	// 10 : fm^{-3} --> 0.1 fm^{-3}

    fach = (double(4) + eps_fit) / (double(1) + eps_fit) ;

    x_nd = (exp (ent_nd) - 1. )/( 1. +exp(ent_nd)/( gam - 1. ) );

    delent = log(1+fach*x_nd*rho_nd_nucl/rho0) / fach - ent_nd;

    ncr_nd = (1. + x_nd) * n0 * rho_nd_nucl/rho0
    			/ exp(ent_nd) / exp(delent);

    unsgam1 = double(1) / ( gam - double(1) ) ;

    gam1sx = ( gam - double(1) - x_nd) / gam / x_nd ;


    cout << "ncr_nd =" << ncr_nd << endl ;

    cout << "eos ent_nd, x_nd, delent,  ent_nd +delent " << ent_nd << "    "<< x_nd << "    "
    << delent <<"      "<< (ent_nd+delent) << endl ;

    double pcr = x_nd * rho_nd_nucl *
    		 pow(( gam - 1. - x_nd)/gam/x_nd *
    		       fabs((exp(ent_nd)-1.)), gam/(gam-1.)) ;
    		
    double psq = ( exp(fach * (ent_nd + delent)) - 1) / fach
   		* rho0 ;
   		
    cout << " eos sqm fach =" << fach << "      PNS_s  " <<
    		psq << "      PNS_cr = " << pcr << endl ;
    //double relp=psq/pcr ;
    cout << " 1+fach ... " << (1+fach*x_nd*rho_nd_nucl/rho0) << endl;
    cout << " log(miu) " <<
    	log(pow((1+fach*x_nd*rho_nd_nucl/rho0),1./fach)) << endl;
    double miu_rate = rho0/n0*pow((1+fach*x_nd*rho_nd_nucl/rho0),1./fach)/
                      rho_nd_nucl/(1-x_nd/(gam-1.))*ncr_nd/exp(ent_nd) ;
    cout << " miu_rate -1 " << miu_rate -1. ;

    double aa0 = log(1+fach*x_nd*rho_nd_nucl/rho0) / fach ;
    double aa1 = log(pow(1+fach*x_nd*rho_nd_nucl/rho0,  1./fach)) ;
    cout << "aa0, aa1,  aa0-aa1 : " << aa0 << "  " <<  aa1 << "  " <<
	    aa0-aa1 << endl ;
    cout << "fach : " << fach << endl ;
    cout << "1+fach*x_nd*rho_nd_nucl/rho0 : " << 1+fach*x_nd*rho_nd_nucl/rho0 << endl ;	
	
}


			//------------------------//
			//  Comparison operators  //
			//------------------------//


bool Eos_strange_cr::operator==(const Eos& eos_i) const {

    bool resu = true ;

    if ( eos_i.identify() != identify() ) {
	cout << "The second EOS is not of type Eos_strange_cr !" << endl ;
	resu = false ;
    }
    else{
	
	const Eos_strange_cr& eos = dynamic_cast<const Eos_strange_cr&>( eos_i ) ;

	if (eos.n0_b60 != n0_b60) {
	    cout
	    << "The two Eos_strange_cr have different n0_b60 : " << n0_b60 << " <-> "
		<< eos.n0_b60 << endl ;
	    resu = false ;
	}

	if (eos.b60 != b60) {
	    cout
	    << "The two Eos_strange_cr have different b60 : " << b60 << " <-> "
		<< eos.b60 << endl ;
	    resu = false ;
	}

	if (eos.ent0 != ent0) {
	    cout
	    << "The two Eos_strange_cr have different ent0 : " << ent0 << " <-> "
		<< eos.ent0 << endl ;
	    resu = false ;
	}

	if (eos.eps_fit != eps_fit) {
	    cout
	    << "The two Eos_strange_cr have different eps_fit : " << eps_fit
		<< " <-> " << eos.eps_fit << endl ;
	    resu = false ;
	}

	if (eos.rho0_b60 != rho0_b60) {
	    cout
	    << "The two Eos_strange_cr have different rho0_b60 : " << rho0_b60
		<< " <-> " << eos.rho0_b60 << endl ;
	    resu = false ;
	}

	
	if (eos.ent_nd != ent_nd) {
	    cout
	    << "The two Eos_strange_cr have different ent_nd : "
	    << ent_nd
	    << " <-> " << eos.ent_nd << endl ;
	    resu = false ;
	}

	
	if (eos.rho_nd != rho_nd) {
	    cout
	    << "The two Eos_strange_cr have different rho_nd : "
	    << rho_nd
	    << " <-> " << eos.rho_nd << endl ;
	    resu = false ;
	}

	
	if (eos.gam != gam) {
	    cout
	    << "The two Eos_strange_cr have different gam : " << gam
		<< " <-> " << eos.gam << endl ;
	    resu = false ;
	}

	
	
    }

    return resu ;

}

bool Eos_strange_cr::operator!=(const Eos& eos_i) const {

    return !(operator==(eos_i)) ;

}

			//------------//
			//  Outputs   //
			//------------//

void Eos_strange_cr::sauve(FILE* fich) const {

    Eos::sauve(fich) ;

    fwrite_be(&n0_b60, sizeof(double), 1, fich) ;	
    fwrite_be(&b60, sizeof(double), 1, fich) ;	
    fwrite_be(&ent0, sizeof(double), 1, fich) ;	
    fwrite_be(&eps_fit, sizeof(double), 1, fich) ;	
    fwrite_be(&rho0_b60, sizeof(double), 1, fich) ;	
    fwrite_be(&ent_nd, sizeof(double), 1, fich) ;	
    fwrite_be(&rho_nd, sizeof(double), 1, fich) ;	
    fwrite_be(&gam, sizeof(double), 1, fich) ;	

}		

ostream& Eos_strange_cr::operator>>(ostream & ost) const {

    ost <<
    "EOS of class Eos_strange_cr " << endl
      << " (Strange matter EOS from Zdunik (2000) with crust) "
    << endl ;
    ost << "   Baryon density at zero pressure : " << n0_b60
	<< " * B_{60}^{3/4}" << endl ;
    ost << "   Bag constant B :  " << b60 << " * 60 MeV/fm^3"<< endl ;
    ost <<
    "   Log-enthalpy threshold for setting the energy density to non-zero: "
	<< endl << "              " << ent0 << endl ;
    ost << "   Fitting parameter eps_fit : " << eps_fit << endl ;
    ost << "   Energy density at zero pressure : " << rho0_b60
	<< " * B_{60}  MeV/fm^3" << endl ;
    ost << "   Log-enthalpy at neutron drip point : " << ent_nd << endl ;
    ost << "   Energy density at neutron drip point : " << rho_nd
        << " MeV/fm^3" << endl ;
    ost << "   Adiabatic index for the crust : " << gam << endl ;

    return ost ;

}


			//------------------------------//
			//    Computational routines    //
			//------------------------------//

// Baryon density from enthalpy
//------------------------------

double Eos_strange_cr::nbar_ent_p(double ent, const Param* ) const {

    if ( ent > ent0 ) {

    	if (ent > ent_nd) {	// Strange quark matter

    		double entsqm = ent + delent ;

		return n0 * exp( double(3) * entsqm / (double(1) + eps_fit))  ;
	}
	else {		// crust

		double fn_au = pow( gam1sx * ( exp(ent)- double(1) ), unsgam1) ;

		return ncr_nd * fn_au ;

	}

    }
    else{
	return 0 ;
    }
}

// Energy density from enthalpy
//------------------------------

double Eos_strange_cr::ener_ent_p(double ent, const Param* ) const {


    if ( ent > ent0 ) {

    	if (ent > ent_nd) {	// Strange quark matter

    		double entsqm = ent + delent ;

		double pp = ( exp(fach * entsqm) - double(1)) / fach  * rho0 ;

		return rho0 + double(3) * pp / (double(1) + eps_fit) ;

	}
	else {		// crust

		double fn_au = pow( gam1sx *(exp(ent) - double(1) ), unsgam1) ;

		return  rho_nd_nucl * ( (double(1) - x_nd * unsgam1 ) * fn_au
			 + x_nd * unsgam1 * pow(fn_au,  gam)) ;
	}

    }
    else{
	return 0 ;
    }
}

// Pressure from enthalpy
//------------------------

double Eos_strange_cr::press_ent_p(double ent, const Param* ) const {

    if ( ent > ent0 ) {

    	if (ent > ent_nd) {	// Strange quark matter

    		double entsqm = ent + delent ;

		return ( exp(fach * entsqm) - double(1) ) / fach  * rho0 ;
    	}
    	else{		// crust

 		double fn_au = pow( gam1sx *(exp(ent) - double(1) ), unsgam1) ;

 		return x_nd * rho_nd_nucl * pow(fn_au, gam) ;

    	}


    }
    else{
	return 0 ;
    }
}


// dln(n)/ln(H) from enthalpy
//---------------------------

double Eos_strange_cr::der_nbar_ent_p(double ent, const Param* ) const {

    if ( ent > ent0 ) {

   	if (ent > ent_nd) {	// Strange quark matter

    		double entsqm = ent + delent ;

		return double(3) * entsqm / ( double(1) +  eps_fit ) ;
	}
       	else{		// crust

	    cout << "Eos_strange_cr::der_nbar_ent_p not ready yet !" << endl ;
	    abort() ;
	    return 0 ;
      	}

    }
    else{
	return 0 ;
    }
}

// dln(e)/ln(H) from enthalpy
//---------------------------

double Eos_strange_cr::der_ener_ent_p(double ent, const Param* ) const {

    if ( ent > ent0 ) {

   	if (ent > ent_nd) {	// Strange quark matter

    		double entsqm = ent + delent ;

		double xx = fach * entsqm ;

		return xx / ( double(1) +
		    ( double(1) + eps_fit ) / double(3) * exp(-xx) ) ;
	}
      	else{		// crust

	    cout << "Eos_strange_cr::der_ener_ent_p not ready yet !" << endl ;
	    abort() ;
	    return 0 ;
       	}
    }
    else{
	return 0 ;
    }
}

// dln(p)/ln(H) from enthalpy
//---------------------------

double Eos_strange_cr::der_press_ent_p(double ent, const Param* ) const {
    
    if ( ent > ent0 ) {


   	if (ent > ent_nd) {	// Strange quark matter
    		
    		double entsqm = ent + delent ;

		double xx = fach * entsqm ; 

		return xx / ( double(1) - exp(-xx) ) ;
	}
      	else{		// crust

	    cout << "Eos_strange_cr::der_press_ent_p not ready yet !" << endl ;
	    abort() ;
 	    return 0 ;      
       	}

    }
    else{
	return 0 ;
    }
}



	    
}
