/*
 * Methods of class Binaire
 *
 * (see file binaire.h for documentation)
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


char binaire_bin_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2001/12/20 13:03:24  k_taniguchi
 * Addition of the Komar mass, the virial error by Gourgoulhon and Bonazzola, and the virial error by Friedman, Uryu, and Shibata.
 *
 * Revision 1.2  2001/12/04 21:27:52  e_gourgoulhon
 *
 * All writing/reading to a binary file are now performed according to
 * the big endian convention, whatever the system is big endian or
 * small endian, thanks to the functions fwrite_be and fread_be
 *
 * Revision 1.1.1.1  2001/11/20 15:19:30  e_gourgoulhon
 * LORENE
 *
 * Revision 2.7  2001/06/25  12:55:46  eric
 * Traitement des etoiles compagnons (appel de Etoile_bin::set_companion dans
 * les constructeurs).
 *
 * Revision 2.6  2000/07/07  14:10:17  eric
 * AJout de display_poly.
 *
 * Revision 2.5  2000/03/13  14:25:52  eric
 *  Ajout des membres p_ham_constr et p_mom_constr.
 *
 * Revision 2.4  2000/02/18  14:52:44  eric
 * Ajout des membres p_virial et p_total_ener.
 * p_masse_adm --> p_mass_adm.
 *
 * Revision 2.3  2000/02/12  17:36:25  eric
 * Ajout de la fonction separation().
 *
 * Revision 2.2  2000/02/04  17:14:47  eric
 * Ajout du membre ref_triad.
 *
 * Revision 2.1  2000/02/02  10:40:36  eric
 * Modif affichage.
 *
 * Revision 2.0  2000/01/31  15:57:49  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers C
#include "math.h"

// Headers Lorene
#include "binaire.h"
#include "eos.h"
#include "utilitaires.h"

			    //--------------//
			    // Constructors //
			    //--------------//

// Standard constructor
// --------------------

Binaire::Binaire(Map& mp1, int nzet1, const Eos& eos1, int irrot1, 
		 Map& mp2, int nzet2, const Eos& eos2, int irrot2,
		 int relat) 
		 : ref_triad(0., "Absolute frame Cartesian basis"),  
		   star1(mp1, nzet1, relat, eos1, irrot1, ref_triad), 
		   star2(mp2, nzet2, relat, eos2, irrot2, ref_triad)
{

    et[0] = &star1 ; 
    et[1] = &star2 ; 
    
    omega = 0 ; 
    x_axe = 0 ; 

    // Companion stars:
    star1.set_companion( star2 ) ; 
    star2.set_companion( star1 ) ; 

    // Pointers of derived quantities initialized to zero : 
    set_der_0x0() ;
}

// Copy constructor
// ----------------
Binaire::Binaire(const Binaire& bibi) 
		: ref_triad(0., "Absolute frame Cartesian basis"), 
		  star1(bibi.star1), 
		  star2(bibi.star2),
		  omega(bibi.omega), 
		  x_axe(bibi.x_axe) 
{
    et[0] = &star1 ; 
    et[1] = &star2 ; 

    // Companion stars:
    star1.set_companion( star2 ) ; 
    star2.set_companion( star1 ) ; 

    // Pointers of derived quantities initialized to zero : 
    set_der_0x0() ;    
}

// Constructor from a file
// -----------------------
Binaire::Binaire(Map& mp1, const Eos& eos1, Map& mp2, const Eos& eos2, 
		 FILE* fich)
		: ref_triad(0., "Absolute frame Cartesian basis"), 
		  star1(mp1, eos1, ref_triad, fich), 
		  star2(mp2, eos2, ref_triad, fich) 
{
    et[0] = &star1 ; 
    et[1] = &star2 ; 

    // Companion stars:
    star1.set_companion( star2 ) ; 
    star2.set_companion( star1 ) ; 

    // omega and x_axe are read in the file:     
    fread_be(&omega, sizeof(double), 1, fich) ;		
    fread_be(&x_axe, sizeof(double), 1, fich) ;		

    // Pointers of derived quantities initialized to zero : 
    set_der_0x0() ;    
    
}			

			    //------------//
			    // Destructor //
			    //------------//

Binaire::~Binaire(){

    del_deriv() ; 

}

			//----------------------------------//
			// Management of derived quantities //
			//----------------------------------//

void Binaire::del_deriv() const {

    if (p_mass_adm != 0x0) delete p_mass_adm ; 
    if (p_mass_kom != 0x0) delete p_mass_kom ; 
    if (p_angu_mom != 0x0) delete p_angu_mom ; 
    if (p_total_ener != 0x0) delete p_total_ener ; 
    if (p_virial != 0x0) delete p_virial ; 
    if (p_virial_gb != 0x0) delete p_virial_gb ; 
    if (p_virial_fus != 0x0) delete p_virial_fus ; 
    if (p_ham_constr != 0x0) delete p_ham_constr ; 
    if (p_mom_constr != 0x0) delete p_mom_constr ; 

    set_der_0x0() ; 
}			    




void Binaire::set_der_0x0() const {

    p_mass_adm = 0x0 ; 
    p_mass_kom = 0x0 ; 
    p_angu_mom = 0x0 ; 
    p_total_ener = 0x0 ; 
    p_virial = 0x0 ; 
    p_virial_gb = 0x0 ; 
    p_virial_fus = 0x0 ; 
    p_ham_constr = 0x0 ; 
    p_mom_constr = 0x0 ; 

}			    


			    //--------------//
			    //  Assignment  //
			    //--------------//

// Assignment to another Binaire
// -----------------------------

void Binaire::operator=(const Binaire& bibi) {

    assert( bibi.ref_triad == ref_triad ) ; 
    
    star1 = bibi.star1 ; 
    star2 = bibi.star2 ; 
    
    omega = bibi.omega ; 
    x_axe = bibi.x_axe ; 
    
    // ref_triad remains unchanged 
    
    del_deriv() ;  // Deletes all derived quantities
    
}

			    //--------------//
			    //	  Outputs   //
			    //--------------//

// Save in a file
// --------------
void Binaire::sauve(FILE* fich) const {
    
    star1.sauve(fich) ; 
    star2.sauve(fich) ; 
    
    fwrite_be(&omega, sizeof(double), 1, fich) ;		
    fwrite_be(&x_axe, sizeof(double), 1, fich) ;		
    
}

// Printing
// --------
ostream& operator<<(ostream& ost, const Binaire& bibi)  {
    bibi >> ost ;
    return ost ;
}
    

ostream& Binaire::operator>>(ostream& ost) const {

    #include "unites.h"	    
    // To avoid some compilation warnings
    if (&ost == 0x0) {
	cout << qpig << msol << mevpfm3 << endl ; 
    }    

    ost << endl ; 
    ost << "Binary system" << endl ; 
    ost << "=============" << endl ; 
    ost << endl << 
	"Orbital angular velocity : " << omega * f_unit << " rad/s" << endl ; 
    ost << endl << 
	"Coordinate separation between the two stellar centers : " 
	<< separation() / km  << " km" << endl ; 
    ost << 
	"Absolute coordinate X of the rotation axis : " << x_axe / km 
	    << " km" << endl ; 
    ost << endl << "Star 1 : " << endl ; 
    ost << "======   " << endl ; 
    ost << star1 << endl ; 
    ost << "Star 2 : " << endl ; 
    ost << "======   " << endl ; 
    ost << star2 << endl ; 
    return ost ;
}

// Display in polytropic units
// ---------------------------

void Binaire::display_poly(ostream& ost) const {

    #include "unites.h"	    
    // To avoid some compilation warnings
    if (&ost == 0x0) {
	cout << f_unit << qpig << msol << mevpfm3 << endl ; 
    }    

    const Eos* p_eos1 = &( star1.get_eos() ) ; 
    const Eos_poly* p_eos_poly = dynamic_cast<const Eos_poly*>( p_eos1 ) ; 	  

    if (p_eos_poly != 0x0) {

	assert( star1.get_eos() == star2.get_eos() ) ; 

	double kappa = p_eos_poly->get_kap() ; 
	double gamma = p_eos_poly->get_gam() ;  ; 
	double kap_ns2 = pow( kappa,  0.5 /(gamma-1) ) ; 
    
	// Polytropic unit of length in terms of r_unit : 
	double r_poly = kap_ns2 / sqrt(ggrav) ; 
    
	// Polytropic unit of time in terms of t_unit :
	double t_poly = r_poly ; 

	// Polytropic unit of mass in terms of m_unit :
	double m_poly = r_poly / ggrav ; 
    
	// Polytropic unit of angular momentum in terms of j_unit :
	double j_poly = r_poly * r_poly / ggrav ; 
    
	ost.precision(10) ; 
	ost << endl << "Quantities in polytropic units : " << endl ; 
	ost	 << "==============================" << endl ; 
	ost << " ( r_poly = " << r_poly / km << " km )" << endl ; 
	ost << "  d_e_max	: " << separation() / r_poly << endl ; 
	ost << "  d_G		: " 
	     << ( star2.xa_barycenter() - star1.xa_barycenter() ) / r_poly 
	     << endl ; 
	ost << "  Omega	  : " << omega * t_poly << endl ; 
	ost << "  J	  : " << angu_mom()(2) / j_poly << endl ; 
	ost << "  M_ADM   : " << mass_adm() / m_poly << endl ; 
	ost << "  M_Komar : " << mass_kom() / m_poly << endl ; 
	ost << "  E	  : " << total_ener() / m_poly << endl ; 
	ost << "  M_bar(star 1) : " << star1.mass_b() / m_poly << endl ; 
	ost << "  M_bar(star 2) : " << star2.mass_b() / m_poly << endl ; 
	ost << "  R_0(star 1)	: " << 
	0.5 * ( star1.ray_eq() + star1.ray_eq_pi() ) / r_poly << endl ;  
	ost << "  R_0(star 2)	: " << 
	0.5 * ( star2.ray_eq() + star2.ray_eq_pi() ) / r_poly << endl ;  
    
    }
    

} 


   
		    //-------------------------------//
		    //		Miscellaneous	     //
		    //-------------------------------//

double Binaire::separation() const {
    
    double dx = star1.get_mp().get_ori_x() - star2.get_mp().get_ori_x() ; 
    double dy = star1.get_mp().get_ori_y() - star2.get_mp().get_ori_y() ; 
    double dz = star1.get_mp().get_ori_z() - star2.get_mp().get_ori_z() ; 
    
    return sqrt( dx*dx + dy*dy + dz*dz ) ; 
    
}
