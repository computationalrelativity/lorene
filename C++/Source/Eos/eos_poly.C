/*
 * Methods of the class Eos_poly.
 *
 * (see file eos.h for documentation).
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
 * Revision 1.14  2023/03/13 15:48:46  j_novak
 * Use of function expm1.
 *
 * Revision 1.13  2023/01/11 15:30:00  j_novak
 * Change the computation at low enthalpies for better accuracy.
 *
 * Revision 1.12  2023/01/10 15:46:11  j_novak
 * Improvement computation of derivatives, including sound speed.
 *
 * Revision 1.11  2021/05/14 15:39:22  g_servignat
 * Added sound speed computation from enthalpy to Eos class and tabulated+polytropic derived classes
 *
 * Revision 1.10  2016/12/05 16:17:51  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.9  2014/10/13 08:52:53  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.8  2014/10/06 15:13:06  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.7  2009/05/25 06:52:27  k_taniguchi
 * Allowed the case of mu_0 != 1 for der_ener_ent_p and der_press_ent_p.
 *
 * Revision 1.6  2003/12/10 08:58:20  r_prix
 * - added new Eos_bifluid paramter for eos-file: bool slow_rot_style
 *  to indicate if we want this particular kind of EOS-inversion (only works for
 *  the  Newtonian 'analytic' polytropes) --> replaces former dirty hack with gamma1<0
 *
 * Revision 1.5  2002/10/16 14:36:35  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.4  2002/04/11 13:28:40  e_gourgoulhon
 * Added the parameter mu_0
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
 * Revision 2.9  2001/02/23  15:17:30  eric
 * Methodes der_nbar_ent_p, der_ener_ent_p, der_press_ent_p :
 *   traitement du cas ent<1.e-13 par un DL
 *   continuite des quantites pour ent<=0.
 *
 * Revision 2.8  2001/02/07  09:50:30  eric
 * Suppression de la fonction derent_ent_p.
 * Ajout des fonctions donnant les derivees de l'EOS:
 *      der_nbar_ent_p
 *      der_ener_ent_p
 *      der_press_ent_p
 *
 * Revision 2.7  2000/06/20  08:34:46  eric
 * Ajout des fonctions get_gam(), etc...
 *
 * Revision 2.6  2000/02/14  14:49:34  eric
 * Modif affichage.
 *
 * Revision 2.5  2000/02/14  14:33:15  eric
 * Ajout du constructeur par lecture de fichier formate.
 *
 * Revision 2.4  2000/01/21  16:05:47  eric
 * Corrige erreur dans set_auxiliary: calcul de gam1.
 *
 * Revision 2.3  2000/01/21  15:18:45  eric
 * Ajout des operateurs de comparaison == et !=
 *
 * Revision 2.2  2000/01/18  14:26:37  eric
 * *** empty log message ***
 *
 * Revision 2.1  2000/01/18  13:47:17  eric
 * Premiere version operationnelle
 *
 * Revision 2.0  2000/01/18  10:46:28  eric
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

			//--------------//
			// Constructors //
			//--------------//

// Standard constructor with m_0 = 1 and mu_0 = 1
// -----------------------------------------------
namespace Lorene {
Eos_poly::Eos_poly(double gam0, double kappa) : 
	Eos("Relativistic polytropic EOS"), 
	gam(gam0), kap(kappa), m_0(double(1)), mu_0(double(1)) {

    set_auxiliary() ; 

}  

// Standard constructor with mu_0 = 1
// ----------------------------------
Eos_poly::Eos_poly(double gam0, double kappa, double mass) :
	Eos("Relativistic polytropic EOS"),
	gam(gam0), kap(kappa), m_0(mass), mu_0(double(1)) {

    set_auxiliary() ;

}

// Standard constructor with mu_0 = 1
// ----------------------------------
Eos_poly::Eos_poly(double gam0, double kappa, double mass, double mu_zero) :
	Eos("Relativistic polytropic EOS"),
	gam(gam0), kap(kappa), m_0(mass), mu_0(mu_zero) {

    set_auxiliary() ;

}

// Copy constructor
// ----------------
Eos_poly::Eos_poly(const Eos_poly& eosi) : 
	Eos(eosi), 
	gam(eosi.gam), kap(eosi.kap), m_0(eosi.m_0), mu_0(eosi.mu_0) {

    set_auxiliary() ; 

}  
  

// Constructor from binary file
// ----------------------------
Eos_poly::Eos_poly(FILE* fich) : 
	Eos(fich) {
        
    fread_be(&gam, sizeof(double), 1, fich) ;		
    fread_be(&kap, sizeof(double), 1, fich) ;		
    fread_be(&m_0, sizeof(double), 1, fich) ;

    if (m_0 < 0) {       // to ensure compatibility with previous version (revision <= 1.2)
                        //  of Eos_poly
        m_0 = fabs( m_0 ) ;
        fread_be(&mu_0, sizeof(double), 1, fich) ;
    }
    else {
        mu_0 = double(1) ;
    }

    set_auxiliary() ; 

}


// Constructor from a formatted file
// ---------------------------------
Eos_poly::Eos_poly(ifstream& fich) : 
	Eos(fich) {

    char blabla[80] ;
        
    fich >> gam ; fich.getline(blabla, 80) ;
    fich >> kap ; fich.getline(blabla, 80) ;
    fich >> m_0 ; fich.getline(blabla, 80) ;

    if (m_0 < 0) {       // to ensure compatibility with previous version (revision <= 1.2)
                        //  of Eos_poly
        m_0 = fabs( m_0 ) ;
        fich >> mu_0 ; fich.getline(blabla, 80) ;
    }
    else {
        mu_0 = double(1) ;
    }

    set_auxiliary() ; 

}
			//--------------//
			//  Destructor  //
			//--------------//

Eos_poly::~Eos_poly(){
    
    // does nothing
        
}
			//--------------//
			//  Assignment  //
			//--------------//

void Eos_poly::operator=(const Eos_poly& eosi) {
    
    set_name(eosi.name) ; 
    
    gam = eosi.gam ; 
    kap = eosi.kap ; 
    m_0 = eosi.m_0 ;
    mu_0 = eosi.mu_0 ;

    set_auxiliary() ;

}


		  //-----------------------//
		  //	Miscellaneous	   //
		  //-----------------------//

void Eos_poly::set_auxiliary() {

    gam1 = gam - double(1) ;

    unsgam1 = double(1) / gam1 ;

    gam1sgamkap = m_0 * gam1 / (gam * kap) ;

    rel_mu_0 = mu_0 / m_0 ;

    ent_0 = log( rel_mu_0 ) ;

}

double Eos_poly::get_gam() const {
    return gam ;
}

double Eos_poly::get_kap() const {
    return kap ; 
}

double Eos_poly::get_m_0() const {
    return m_0 ;
}

double Eos_poly::get_mu_0() const {
    return mu_0 ;
}


			//------------------------//
			//  Comparison operators  //
			//------------------------//


bool Eos_poly::operator==(const Eos& eos_i) const {
    
    bool resu = true ; 
    
    if ( eos_i.identify() != identify() ) {
	cout << "The second EOS is not of type Eos_poly !" << endl ; 
	resu = false ; 
    }
    else{
	
	const Eos_poly& eos = dynamic_cast<const Eos_poly&>( eos_i ) ; 

	if (eos.gam != gam) {
	    cout 
	    << "The two Eos_poly have different gamma : " << gam << " <-> " 
		<< eos.gam << endl ; 
	    resu = false ; 
	}

	if (eos.kap != kap) {
	    cout 
	    << "The two Eos_poly have different kappa : " << kap << " <-> " 
		<< eos.kap << endl ; 
	    resu = false ; 
	}

	if (eos.m_0 != m_0) {
	    cout
	    << "The two Eos_poly have different m_0 : " << m_0 << " <-> "
		<< eos.m_0 << endl ;
	    resu = false ;
	}

	if (eos.mu_0 != mu_0) {
	    cout
	    << "The two Eos_poly have different mu_0 : " << mu_0 << " <-> "
		<< eos.mu_0 << endl ;
	    resu = false ;
	}

    }
    
    return resu ; 
    
}

bool Eos_poly::operator!=(const Eos& eos_i) const {
 
    return !(operator==(eos_i)) ; 
       
}


			//------------//
			//  Outputs   //
			//------------//

void Eos_poly::sauve(FILE* fich) const {

    Eos::sauve(fich) ; 
    
    fwrite_be(&gam, sizeof(double), 1, fich) ;	
    fwrite_be(&kap, sizeof(double), 1, fich) ;
    double tempo = - m_0 ;
    fwrite_be(&tempo, sizeof(double), 1, fich) ;
    fwrite_be(&mu_0, sizeof(double), 1, fich) ;

}

ostream& Eos_poly::operator>>(ostream & ost) const {
    
    ost << "EOS of class Eos_poly (relativistic polytrope) : " << endl ; 
    ost << "   Adiabatic index gamma :      " << gam << endl ; 
    ost << "   Pressure coefficient kappa : " << kap << 
	   " rho_nuc c^2 / n_nuc^gamma" << endl ; 
    ost << "   Mean particle mass : " << m_0 << " m_B" << endl ;
    ost << "   Relativistic chemical potential at zero pressure : " << mu_0 << " m_B c^2" << endl ;

    return ost ;

}


			//------------------------------//
			//    Computational routines    //
			//------------------------------//

// Baryon density from enthalpy
//------------------------------

double Eos_poly::nbar_ent_p(double ent, const Param* ) const {

    if ( ent > ent_0 ) {

      double xx = ent - ent_0 ;
      return pow( gam1sgamkap * rel_mu_0 * expm1(xx), unsgam1 ) ;
    }
    else{
	return 0 ;
    }
}

// Energy density from enthalpy
//------------------------------

double Eos_poly::ener_ent_p(double ent, const Param* ) const {

  double nn = nbar_ent_p(ent) ;
  double pp = kap * pow( nn, gam ) ;

  return  unsgam1 * pp + mu_0 * nn ;

}

// Pressure from enthalpy
//------------------------

double Eos_poly::press_ent_p(double ent, const Param* ) const {

  double nn = nbar_ent_p(ent) ;

  return kap * pow( nn, gam ) ;

}

// dln(n)/ln(H) from enthalpy
//---------------------------

double Eos_poly::der_nbar_ent_p(double ent, const Param* ) const {

    if ( ent > ent_0 )
      {
	double xx = ent - ent_0 ;
	return ent / (- expm1(-xx)) / gam1 ;
      }
    else{
      return double(1) / gam1 ;	//  to ensure continuity at ent=0
    }
}

// dln(e)/ln(H) from enthalpy
//---------------------------

double Eos_poly::der_ener_ent_p(double ent, const Param* ) const {

  if ( ent > ent_0 ) {
    
    double nn = nbar_ent_p(ent) ;
    double pp = kap * pow( nn, gam ) ;
    double ee =  unsgam1 * pp + mu_0 * nn ;

    double xx = ent - ent_0 ;
    return ent / (- expm1(-xx)) / gam1 * ( double(1) + pp / ee) ;
  }
  else{
    return double(1) / gam1 ;   //  to ensure continuity at ent=0
  }
}

// dln(p)/ln(H) from enthalpy
//---------------------------

double Eos_poly::der_press_ent_p(double ent, const Param* ) const {
    
    if ( ent > ent_0 )
      {
	double xx = ent - ent_0 ;
	return gam * ent / (- expm1(-xx)) / gam1 ;
      }
    else{
      return gam / gam1 ;	//  to ensure continuity at ent=0
    }
}
  
// Sound speed from enthalpy
//---------------------------------
double Eos_poly::csound_square_ent_p(double ent, const Param* ) const {

  if ( ent > ent_0 ) {
    double xx = ent - ent_0 ;
    return gam1*( - expm1(-xx) ) ;
  }
  else
    return 0. ;
}

}
