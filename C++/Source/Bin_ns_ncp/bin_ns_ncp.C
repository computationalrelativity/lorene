/*
 * Methods of class Bin_ns_ncp.C
 *
 */

/*
 *   Copyright (c) 2002 Francois Limousin
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
char Bin_ns_ncp_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.6  2003/06/20 13:49:53  f_limousin
 * Add a new argument conf_flat in the constructors and a new function fait_decouple().
 *
 * Revision 1.5  2003/03/03 19:32:06  f_limousin
 * Suppression of the member ref_triad.
 *
 * Revision 1.4  2003/02/12 18:46:59  f_limousin
 * Change the arguments of the standard constructor.
 *
 * Revision 1.3  2003/01/20 17:13:25  j_novak
 * Modif des include <math.h> pour eviter les warning sous SGI.
 *
 * Revision 1.2  2003/01/20 09:38:17  f_limousin
 * Modification of the standard constructor
 *
 * Revision 1.1  2003/01/14 13:43:38  f_limousin
 * Methods of class Bin_ns_ncp.
 *
 * Revision 1.2  2001/12/11 06:44:41  e_gourgoulhon
 * template files
 *
 *
 *
 * $Header$
 *
 */

// Headers C
#include <math.h>

// Headers Lorene
#include "bin_ns_ncp.h"
#include "eos.h"
#include "utilitaires.h"
#include "graphique.h"
#include "param.h"

			    //--------------//
			    // Constructors //
			    //--------------//

// Standard constructor
// --------------------

Bin_ns_ncp::Bin_ns_ncp(Map& mp1, int nzet1, const Eos& eos1, int irrot1, 
		 Map& mp2, int nzet2, const Eos& eos2, int irrot2, int relat,
		 int conf_flat, const Metrique& flat1, const Metrique& flat2,
		       const Tenseur_sym &source1, const Tenseur_sym &source2) 
                 : star1(mp1, nzet1, relat, eos1, irrot1, conf_flat,
			 mp1.get_bvect_cart(), flat1, source1), 
		   star2(mp2, nzet2, relat, eos2, irrot2, conf_flat,
			 mp2.get_bvect_cart(), flat2, source2)
{

    et[0] = &star1 ; 
    et[1] = &star2 ; 
    
    omega = 0 ; 
    x_axe = 0 ; 

    // Pointers of derived quantities initialized to zero :
    set_der_0x0() ;
}

// Copy constructor
// ----------------
Bin_ns_ncp::Bin_ns_ncp(const Bin_ns_ncp& bibi) 
		: star1(bibi.star1), 
		  star2(bibi.star2),
		  omega(bibi.omega), 
		  x_axe(bibi.x_axe) 
{
    et[0] = &star1 ; 
    et[1] = &star2 ; 

    // Pointers of derived quantities initialized to zero :
    set_der_0x0() ;    
}

// Constructor from a file
// -----------------------
Bin_ns_ncp::Bin_ns_ncp(Map& mp1, const Eos& eos1, Map& mp2, const Eos& eos2, 
		 const Metrique& flat1, const Metrique& flat2, FILE* fich)
		: star1(mp1, eos1, mp1.get_bvect_cart(), flat1, fich), 
		  star2(mp2, eos2, mp2.get_bvect_cart(), flat2, fich) 
{
    et[0] = &star1 ; 
    et[1] = &star2 ; 

    // omega and x_axe are read in the file:
    fread_be(&omega, sizeof(double), 1, fich) ;		
    fread_be(&x_axe, sizeof(double), 1, fich) ;		

    // Pointers of derived quantities initialized to zero : 
    set_der_0x0() ;    
    
}			

			    //------------//
			    // Destructor //
			    //------------//

Bin_ns_ncp::~Bin_ns_ncp(){

    del_deriv() ; 

}

			//----------------------------------//
			// Management of derived quantities //
			//----------------------------------//

void Bin_ns_ncp::del_deriv() const {

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




void Bin_ns_ncp::set_der_0x0() const {

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

// Assignment to another Bin_ns_ncp
// --------------------------------

void Bin_ns_ncp::operator=(const Bin_ns_ncp& bibi) {

    star1 = bibi.star1 ; 
    star2 = bibi.star2 ; 
    
    omega = bibi.omega ; 
    x_axe = bibi.x_axe ; 
    
    del_deriv() ;  // Deletes all derived quantities
    
}

			    //--------------//
			    //	  Outputs   //
			    //--------------//

// Save in a file
// --------------
void Bin_ns_ncp::sauve(FILE* fich) const {
    
    star1.sauve(fich) ; 
    star2.sauve(fich) ; 
    
    fwrite_be(&omega, sizeof(double), 1, fich) ;		
    fwrite_be(&x_axe, sizeof(double), 1, fich) ;		
    
}

// Printing
// --------
ostream& operator<<(ostream& ost, const Bin_ns_ncp& bibi)  {
    bibi >> ost ;
    return ost ;
}
    

ostream& Bin_ns_ncp::operator>>(ostream& ost) const {

    #include "unites.h"	    
    // To avoid some compilation warnings
    if (&ost == 0x0) {
	cout << qpig << msol << mevpfm3 << endl ; 
    }    

    ost << endl ; 
    ost << "Binary neutron stars with non comformally flat metric" << endl ; 
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

void Bin_ns_ncp::display_poly(ostream& ost) const {

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
	//	double j_poly = r_poly * r_poly / ggrav ; 
    
	ost.precision(10) ; 
	ost << endl << "Quantities in polytropic units : " << endl ; 
	ost	 << "==============================" << endl ; 
	ost << " ( r_poly = " << r_poly / km << " km )" << endl ; 
	ost << "  d_e_max	: " << separation() / r_poly << endl ; 
	ost << "  d_G		: " 
	     << ( star2.xa_barycenter() - star1.xa_barycenter() ) / r_poly 
	     << endl ; 
	ost << "  Omega	  : " << omega * t_poly << endl ; 
	//	ost << "  J	  : " << angu_mom()(2) / j_poly << endl ; 
	//      ost << "  M_ADM   : " << mass_adm() / m_poly << endl ;      
	//      ost << "  M_Komar : " << mass_kom() / m_poly << endl ; 
	//	ost << "  E	  : " << total_ener() / m_poly << endl ; 
	ost << "  M_bar(star 1) : " << star1.mass_b() / m_poly << endl ; 
	ost << "  M_bar(star 2) : " << star2.mass_b() / m_poly << endl ; 
	ost << "  R_0(star 1)	: " << 
	0.5 * ( star1.ray_eq() + star1.ray_eq_pi() ) / r_poly << endl ;  
	ost << "  R_0(star 2)	: " << 
	0.5 * ( star2.ray_eq() + star2.ray_eq_pi() ) / r_poly << endl ;  
    
    }
    

} 


void Bin_ns_ncp::fait_decouple () {
    
    int nz_un = star1.mp.get_mg()->get_nzone() ;
    int nz_deux = star2.mp.get_mg()->get_nzone() ;
    
    // On determine R_limite (pour le moment en tout cas...) :
    double distance = fabs(star1.mp.get_ori_x() - star2.mp.get_ori_x()) ;
    double lim_un = -1*distance/2. ;
    double lim_deux = -1*distance/2. ;
    double int_un = 0*distance/6. ;
    double int_deux = 0*distance/6. ;
    

    /*
    // Les fonctions de base
    Cmp fonction_f_un (star1.mp) ;
    //   fonction_f_un = (exp(-pow(star1.mp.r/lim_un, 2)) - exp(-1.)) / (1-exp(-1.))/2. + 0.5 ;
    
    fonction_f_un = 0.5*pow(
      cos((star1.mp.r-int_un)*M_PI/2./(lim_un-int_un)), 2.)+0.5 ;
    fonction_f_un.std_base_scal();

    des_coupe_z(fonction_f_un, 0, 2) ;
    des_profile(fonction_f_un, 0, 10, 0, 0) ;
    des_coef_xi(fonction_f_un.va, 0, 0, 0) ;
    des_coef_xi(fonction_f_un.va, 1, 0, 0) ;
    des_coef_xi(fonction_f_un.va, 2, 0, 0) ;
    
    Cmp fonction_g_un (star1.mp) ;
    //   fonction_g_un = (1 - exp(-pow(star1.mp.r/lim_un, 2))) /
    //(1-exp(-1.))/2. ;

    fonction_g_un = 0.5*pow
      (sin((star1.mp.r-int_un)*M_PI/2./(lim_un-int_un)), 2.) ;
    fonction_g_un.std_base_scal();
    
    Cmp fonction_f_deux (star2.mp) ;
    fonction_f_deux = 0.5*pow(
 cos((star2.mp.r-int_deux)*M_PI/2./(lim_deux-int_deux)), 2.)+0.5 ;
    fonction_f_deux.std_base_scal();
    
    Cmp fonction_g_deux (star2.mp) ;
    fonction_g_deux = 0.5*pow
 (sin((star2.mp.r-int_deux)*M_PI/2./(lim_deux-int_deux)), 2.) ;
    fonction_g_deux.std_base_scal();
    */   


     // Les fonctions totales :
    Cmp decouple_un (star1.mp) ;
    decouple_un.allocate_all() ;
    Cmp decouple_deux (star2.mp) ;
    decouple_deux.allocate_all() ;
    
    Mtbl xabs_un (star1.mp.xa) ;
    Mtbl yabs_un (star1.mp.ya) ;
    Mtbl zabs_un (star1.mp.za) ;
	    
    Mtbl xabs_deux (star2.mp.xa) ;
    Mtbl yabs_deux (star2.mp.ya) ;
    Mtbl zabs_deux (star2.mp.za) ;
	    
    double xabs, yabs, zabs, air_un, air_deux, theta, phi ;
	    
    // On boucle sur les autres zones :
    for (int l=0 ; l<nz_un ; l++) {
	int nr = star1.mp.get_mg()->get_nr (l) ;
		
	if (l==nz_un-1)
	    nr -- ;
		
	int np = star1.mp.get_mg()->get_np (l) ;
	int nt = star1.mp.get_mg()->get_nt (l) ;
		
	for (int k=0 ; k<np ; k++)
	    for (int j=0 ; j<nt ; j++)
		for (int i=0 ; i<nr ; i++) {
			    
		    xabs = xabs_un (l, k, j, i) ;
		    yabs = yabs_un (l, k, j, i) ;
		    zabs = zabs_un (l, k, j, i) ;
			    
		    // les coordonnees du point :
		    star1.mp.convert_absolute 
			(xabs, yabs, zabs, air_un, theta, phi) ;
		    star2.mp.convert_absolute 
			(xabs, yabs, zabs, air_deux, theta, phi) ;

		    if (air_un <= lim_un)
			if (air_un < int_un)
			    decouple_un.set(l, k, j, i) = 1 ;
			else
			// pres de l'etoile une :
			decouple_un.set(l, k, j, i) =  0.5*pow(
      cos((air_un-int_un)*M_PI/2./(lim_un-int_un)), 2.)+0.5 ;

		    else 
			if (air_deux <= lim_deux)
			    if (air_deux < int_deux)
				decouple_un.set(l, k, j, i) = 0 ;
			    else
			// On est pres de l'etoile deux :
			     decouple_un.set(l, k, j, i) = 0.5*pow
      (sin((air_deux-int_deux)*M_PI/2./(lim_deux-int_deux)), 2.) ;
						    
			else
			    // On est loin des deux etoiles :
			    decouple_un.set(l, k, j, i) = 0.5 ;
		}
	
    
	        // Cas infini :
		if (l==nz_un-1)
		    for (int k=0 ; k<np ; k++)
			for (int j=0 ; j<nt ; j++)
			    decouple_un.set(nz_un-1, k, j, nr) = 0.5 ;
    }


    for (int l=0 ; l<nz_deux ; l++) {
	int nr = star2.mp.get_mg()->get_nr (l) ;
		
	if (l==nz_deux-1)
	    nr -- ;
		
	int np = star2.mp.get_mg()->get_np (l) ;
	int nt = star2.mp.get_mg()->get_nt (l) ;
		
	for (int k=0 ; k<np ; k++)
	    for (int j=0 ; j<nt ; j++)
		for (int i=0 ; i<nr ; i++) {
			    
		    xabs = xabs_deux (l, k, j, i) ;
		    yabs = yabs_deux (l, k, j, i) ;
		    zabs = zabs_deux (l, k, j, i) ;
			    
		    // les coordonnees du point  :
		    star1.mp.convert_absolute 
			(xabs, yabs, zabs, air_un, theta, phi) ;
		    star2.mp.convert_absolute 
			(xabs, yabs, zabs, air_deux, theta, phi) ;
		    
		    if (air_deux <= lim_deux)
			if (air_deux < int_deux)
			    decouple_deux.set(l, k, j, i) = 1 ;
			else
			  // pres de l'etoile deux :
			decouple_deux.set(l, k, j, i) =  0.5*pow(
	       cos((air_deux-int_deux)*M_PI/2./(lim_deux-int_deux)), 2.)+0.5 ;
			
		    else 
			if (air_un <= lim_un)
			    if (air_un < int_un)
				decouple_deux.set(l, k, j, i) = 0 ;
			    else
			// On est pres de l'etoile une :
			     decouple_deux.set(l, k, j, i)=0.5*pow
               (sin((air_un-int_un)*M_PI/2./(lim_un-int_un)), 2.) ;
		   
		
			else
			    // On est loin des deux etoiles :
			    decouple_deux.set(l, k, j, i) = 0.5 ;
		}
			    
		// Cas infini :
		if (l==nz_deux-1)
		    for (int k=0 ; k<np ; k++)
			for (int j=0 ; j<nt ; j++)
			    decouple_deux.set(nz_un-1, k, j, nr) = 0.5 ;
   }
   
    double somme = 0 ;

    int nr = star2.mp.get_mg()->get_nr (2) ;
    int np = star2.mp.get_mg()->get_np (2) ;
    int nt = star2.mp.get_mg()->get_nt (2) ;
    for (int k=0 ; k<np ; k++)
      for (int j=0 ; j<nt ; j++)
	for (int i=0 ; i<nr ; i++) {
	  somme += decouple_un(2,k,j,5) ;
	}
 
    cout << "decouple_un"  << endl << norme(decouple_un) << endl ;
    cout << "decouple_deux"  << endl << norme(decouple_deux) << endl ;
    /*
    decouple_un.std_base_scal() ;

    des_coef_xi(decouple_un.va, 0, 0, 0) ;
    des_coef_xi(decouple_un.va, 1, 0, 0) ;
    des_coef_xi(decouple_un.va, 2, 0, 0) ;

    decouple_deux.std_base_scal() ;

    des_coef_xi(decouple_deux.va, 0, 0, 0) ;
    des_coef_xi(decouple_deux.va, 1, 0, 0) ;
    des_coef_xi(decouple_deux.va, 2, 0, 0) ;
    */
    star1.decouple = decouple_un ;
    star2.decouple = decouple_deux ;

     
}



   
		    //-------------------------------//
		    //		Miscellaneous	     //
		    //-------------------------------//

double Bin_ns_ncp::separation() const {
    
    double dx = star1.mp.get_ori_x() - star2.mp.get_ori_x() ; 
    double dy = star1.mp.get_ori_y() - star2.mp.get_ori_y() ; 
    double dz = star1.mp.get_ori_z() - star2.mp.get_ori_z() ; 
    
    return sqrt( dx*dx + dy*dy + dz*dz ) ; 
    
}
