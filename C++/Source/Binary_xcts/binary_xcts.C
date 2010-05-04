/*
 * Methods of class Binary_xcts
 * (see file binary_xcts.h for documentation)
 */

/*
 *   Copyright (c) 2010 Michal Bejger
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

char binary_xcts_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2010/05/04 07:35:54  m_bejger
 * Initial version
 *
 * $Header$
 *
 */
 
// Headers C
#include <math.h>

// Headers Lorene
#include "binary_xcts.h"
#include "eos.h"
#include "utilitaires.h"
#include "graphique.h"
#include "param.h"
#include "unites.h"	    

			    //--------------//
			    // Constructors //
			    //--------------//

// Standard constructor
// --------------------

Binary_xcts::Binary_xcts(Map& mp1, 
						int nzet1, 
						const Eos& eos1, 
						int irrot1, 
	       				Map& mp2, 
	       				int nzet2, 
	       				const Eos& eos2, 
	       				int irrot2) 
                 : star1(mp1, nzet1, eos1, irrot1), 
		   		   star2(mp2, nzet2, eos2, irrot2)
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
Binary_xcts::Binary_xcts(const Binary_xcts& bibi) 
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
Binary_xcts::Binary_xcts(Map& mp1, 
						 const Eos& eos1, 
						 Map& mp2, 
						 const Eos& eos2, 
	       				 FILE* fich)
						: star1(mp1, eos1, fich), 
		  				  star2(mp2, eos2, fich) 
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
			    //------- -----//

Binary_xcts::~Binary_xcts(){

    del_deriv() ; 

}

			//----------------------------------//
			// Management of derived quantities //
			//----------------------------------//

void Binary_xcts::del_deriv() const {

    if (p_mass_adm != 0x0) delete p_mass_adm ; 
    if (p_mass_kom != 0x0) delete p_mass_kom ; 
    if (p_angu_mom != 0x0) delete p_angu_mom ; 
    if (p_total_ener != 0x0) delete p_total_ener ; 
    if (p_virial != 0x0) delete p_virial ; 
    if (p_ham_constr != 0x0) delete p_ham_constr ; 
    if (p_mom_constr != 0x0) delete p_mom_constr ; 

    set_der_0x0() ; 
}			    




void Binary_xcts::set_der_0x0() const {

    p_mass_adm = 0x0 ; 
    p_mass_kom = 0x0 ; 
    p_angu_mom = 0x0 ; 
    p_total_ener = 0x0 ; 
    p_virial = 0x0 ; 
    p_ham_constr = 0x0 ; 
    p_mom_constr = 0x0 ; 

}			    


			    //--------------//
			    //  Assignment  //
			    //--------------//

// Assignment to another Binary_xcts
// --------------------------------

void Binary_xcts::operator=(const Binary_xcts& bibi) {

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
void Binary_xcts::sauve(FILE* fich) const {
    
    star1.sauve(fich) ; 
    star2.sauve(fich) ; 
    
    fwrite_be(&omega, sizeof(double), 1, fich) ;		
    fwrite_be(&x_axe, sizeof(double), 1, fich) ;		
    
}

// Printing
// --------
ostream& operator<<(ostream& ost, const Binary_xcts& bibi)  {
	
    bibi >> ost ;
    return ost ;
}
    

ostream& Binary_xcts::operator>>(ostream& ost) const {

  using namespace Unites ;

    ost << endl ; 
    ost << "Binary neutron stars" << endl ; 
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

void Binary_xcts::display_poly(ostream& ost) const {

  using namespace Unites ;

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
	ost << "  M_bar(star 1) : " << star1.mass_b() / m_poly << endl ; 
	ost << "  M_bar(star 2) : " << star2.mass_b() / m_poly << endl ; 
	ost << "  R_0(star 1)	: " << 
	0.5 * ( star1.ray_eq() + star1.ray_eq_pi() ) / r_poly << endl ;  
	ost << "  R_0(star 2)	: " << 
	0.5 * ( star2.ray_eq() + star2.ray_eq_pi() ) / r_poly << endl ;  
    
    }
    
} 

void Binary_xcts::fait_decouple () {

    int nz_un = star1.get_mp().get_mg()->get_nzone() ;
    int nz_deux = star2.get_mp().get_mg()->get_nzone() ;
    
    // We determine R_limite :
    double distance = star2.get_mp().get_ori_x() - star1.get_mp().get_ori_x() ;
    cout << "distance = " << distance << endl ;
    double lim_un = distance/2. ;
    double lim_deux = distance/2. ;
    double int_un = distance/6. ;
    double int_deux = distance/6. ;
    
    // The functions used.
    Scalar fonction_f_un (star1.get_mp()) ;
    fonction_f_un = 0.5*pow(
	cos((star1.get_mp().r-int_un)*M_PI/2./(lim_un-int_un)), 2.)+0.5 ;
    fonction_f_un.std_spectral_base();
    
    Scalar fonction_g_un (star1.get_mp()) ;
    fonction_g_un = 0.5*pow
	(sin((star1.get_mp().r-int_un)*M_PI/2./(lim_un-int_un)), 2.) ;
    fonction_g_un.std_spectral_base();
    
    Scalar fonction_f_deux (star2.get_mp()) ;
    fonction_f_deux = 0.5*pow(
	cos((star2.get_mp().r-int_deux)*M_PI/2./(lim_deux-int_deux)), 2.)+0.5 ;
    fonction_f_deux.std_spectral_base();
    
    Scalar fonction_g_deux (star2.get_mp()) ;
    fonction_g_deux = 0.5*pow(
	sin((star2.get_mp().r-int_deux)*M_PI/2./(lim_un-int_deux)), 2.) ;
    fonction_g_deux.std_spectral_base();
    
    // The functions total :
    Scalar decouple_un (star1.get_mp()) ;
    decouple_un.allocate_all() ;
    Scalar decouple_deux (star2.get_mp()) ;
    decouple_deux.allocate_all() ;
    
    Mtbl xabs_un (star1.get_mp().xa) ;
    Mtbl yabs_un (star1.get_mp().ya) ;
    Mtbl zabs_un (star1.get_mp().za) ;
	    
    Mtbl xabs_deux (star2.get_mp().xa) ;
    Mtbl yabs_deux (star2.get_mp().ya) ;
    Mtbl zabs_deux (star2.get_mp().za) ;
	    
    double xabs, yabs, zabs, air_un, air_deux, theta, phi ;
	    
    for (int l=0 ; l<nz_un ; l++) {
	int nr = star1.get_mp().get_mg()->get_nr (l) ;
		
	if (l==nz_un-1)
	    nr -- ;
		
	int np = star1.get_mp().get_mg()->get_np (l) ;
	int nt = star1.get_mp().get_mg()->get_nt (l) ;
		
	for (int k=0 ; k<np ; k++)
	    for (int j=0 ; j<nt ; j++)
		for (int i=0 ; i<nr ; i++) {
			    
		    xabs = xabs_un (l, k, j, i) ;
		    yabs = yabs_un (l, k, j, i) ;
		    zabs = zabs_un (l, k, j, i) ;
			    
		    // Coordinates of the point
		    star1.get_mp().convert_absolute 
			(xabs, yabs, zabs, air_un, theta, phi) ;
		    star2.get_mp().convert_absolute 
			(xabs, yabs, zabs, air_deux, theta, phi) ;
		
		    if (air_un <= lim_un)
			if (air_un < int_un)
			    decouple_un.set_grid_point(l, k, j, i) = 1 ;
			else
			// Close to star 1 :
			decouple_un.set_grid_point(l, k, j, i) = 
			    fonction_f_un.val_grid_point(l, k, j, i) ;
		    else 
			if (air_deux <= lim_deux)
			    if (air_deux < int_deux)
				decouple_un.set_grid_point(l, k, j, i) = 0 ;
			    else
			// Close to star 2 :
				decouple_un.set_grid_point(l, k, j, i) = 
		fonction_g_deux.val_point (air_deux, theta, phi) ;
		
			else
			    // Far from each star :
			    decouple_un.set_grid_point(l, k, j, i) = 0.5 ;
		}
			    
		// Case infinity :
		if (l==nz_un-1)
		    for (int k=0 ; k<np ; k++)
			for (int j=0 ; j<nt ; j++)
			    decouple_un.set_grid_point(nz_un-1, k, j, nr)=0.5 ;
	    }
    
    for (int l=0 ; l<nz_deux ; l++) {
	int nr = star2.get_mp().get_mg()->get_nr (l) ;
		
	if (l==nz_deux-1)
	    nr -- ;
		
	int np = star2.get_mp().get_mg()->get_np (l) ;
	int nt = star2.get_mp().get_mg()->get_nt (l) ;
		
	for (int k=0 ; k<np ; k++)
	    for (int j=0 ; j<nt ; j++)
		for (int i=0 ; i<nr ; i++) {
			    
		    xabs = xabs_deux (l, k, j, i) ;
		    yabs = yabs_deux (l, k, j, i) ;
		    zabs = zabs_deux (l, k, j, i) ;
			    
		    // coordinates of the point  :
		    star1.get_mp().convert_absolute 
			(xabs, yabs, zabs, air_un, theta, phi) ;
		    star2.get_mp().convert_absolute 
			(xabs, yabs, zabs, air_deux, theta, phi) ;
		    
		    if (air_deux <= lim_deux)
			if (air_deux < int_deux)
			    decouple_deux.set_grid_point(l, k, j, i) = 1 ;
			else
			// close to star two :
			decouple_deux.set_grid_point(l, k, j, i) = 
			    fonction_f_deux.val_grid_point(l, k, j, i) ;
		    else 
			if (air_un <= lim_un)
			    if (air_un < int_un)
				decouple_deux.set_grid_point(l, k, j, i) = 0 ;
			    else
			// close to star one :
				decouple_deux.set_grid_point(l, k, j, i) = 
			 fonction_g_un.val_point (air_un, theta, phi) ;
		
			else
			    // Far from each star :
			    decouple_deux.set_grid_point(l, k, j, i) = 0.5 ;
		}
			    
		// Case infinity :
		if (l==nz_deux-1)
		    for (int k=0 ; k<np ; k++)
			for (int j=0 ; j<nt ; j++)
			 decouple_deux.set_grid_point(nz_un-1, k, j, nr)=0.5 ;
   }
   
    decouple_un.std_spectral_base() ;
    decouple_deux.std_spectral_base() ;

    int nr = star1.get_mp().get_mg()->get_nr (1) ;
    int nt = star1.get_mp().get_mg()->get_nt (1) ;
    int np = star1.get_mp().get_mg()->get_np (1) ;
    cout << "decouple_un"  << endl << norme(decouple_un/(nr*nt*np)) << endl ;
    cout << "decouple_deux"  << endl << norme(decouple_deux/(nr*nt*np)) 
	 << endl ;

    star1.decouple = decouple_un ;
    star2.decouple = decouple_deux ;
}

void Binary_xcts::write_global(ostream& ost) const {

  using namespace Unites ;

  const Map&  mp1 = star1.get_mp() ;
  const Mg3d* mg1 = mp1.get_mg() ;
  int nz1 = mg1->get_nzone() ; 

  ost.precision(5) ;
  ost << "# Grid 1 : " << nz1 << "x"
      << mg1->get_nr(0) << "x" << mg1->get_nt(0) << "x" << mg1->get_np(0) 
      << "  R_out(l) [km] : " ;
  for (int l=0; l<nz1; l++) {
    ost << " " << mp1.val_r(l, 1., M_PI/2, 0) / km ; 
  }
  ost << endl ; 

  ost << "#     VE(M)	   " << endl ;
  
  
  ost.setf(ios::scientific) ; 
  ost.width(14) ; 
  ost << virial() << endl ;
  
  ost << "#      d [km]         "  
      << "       d_G [km]       "
      << "     d/(a1 +a1')      "
      << "       f [Hz]         "
      << "    M_ADM [M_sol]     "     
      << "    M_ADM_vol [M_sol]     "     
      << "    M_Komar [M_sol]     "     
      << "    M_Komar_vol [M_sol]     "     
      << "   J [G M_sol^2/c]    "  << endl ;   
  
  ost.precision(14) ;
  ost.width(20) ; 
  ost << separation() / km ; ost.width(22) ;
  ost	<< ( star2.xa_barycenter() - star1.xa_barycenter() ) / km ; ost.width(22) ;
  ost	<< separation() / (star1.ray_eq() + star2.ray_eq()) ; ost.width(22) ;
  ost	<< omega / (2*M_PI)* f_unit ; ost.width(22) ;
  ost	<< mass_adm() / msol ; ost.width(22) ; 
  ost	<< mass_adm_vol() / msol ; ost.width(22) ; 
  ost	<< mass_kom() / msol ; ost.width(22) ; 
  ost	<< mass_kom_vol() / msol ; ost.width(22) ; 
  ost	<< angu_mom()(2)/ ( qpig / (4* M_PI) * msol*msol) << endl ; 
  
  ost << "#     H_c(1)[c^2]     "
      << "    e_c(1)[rho_nuc]   " 
      << "    M_B(1) [M_sol]    "
      << "     r_eq(1) [km]     "
      << "        a2/a1(1)	  " 
	    << "        a3/a1(1)	  " << endl ; 
  
  ost.width(20) ; 
  ost << star1.get_ent().val_grid_point(0,0,0,0) ; ost.width(22) ;
  ost	<< star1.get_ener().val_grid_point(0,0,0,0) ; ost.width(22) ;
  ost	<< star1.mass_b() / msol ; ost.width(22) ;	
  ost << star1.ray_eq() / km ; ost.width(22) ; 
  ost	<< star1.ray_eq_pis2() / star1.ray_eq() ; ost.width(22) ;
  ost	<< star1.ray_pole() / star1.ray_eq() << endl ;
  
  ost << "#     H_c(2)[c^2]     "
      << "    e_c(2)[rho_nuc]   " 
      << "    M_B(2) [M_sol]    "
      << "     r_eq(2) [km]     "
      << "        a2/a1(2)	  " 
      << "        a3/a1(2)	  " << endl ; 
  
  ost.width(20) ; 
  ost << star2.get_ent().val_grid_point(0,0,0,0) ; ost.width(22) ;
  ost	<< star2.get_ener().val_grid_point(0,0,0,0) ; ost.width(22) ;
  ost	<< star2.mass_b() / msol ; ost.width(22) ;	
  ost << star2.ray_eq() / km ; ost.width(22) ; 
  ost	<< star2.ray_eq_pis2() / star1.ray_eq() ; ost.width(22) ;
  ost	<< star2.ray_pole() / star1.ray_eq() << endl ;
  
  // Quantities in polytropic units if the EOS is a polytropic one
  // -------------------------------------------------------------
  const Eos* p_eos1 = &( star1.get_eos() ) ; 
  const Eos_poly* p_eos_poly = dynamic_cast<const Eos_poly*>( p_eos1 ) ; 	  
  
  if ((p_eos_poly != 0x0) && ( star1.get_eos() == star2.get_eos() )) {
      
      double kappa = p_eos_poly->get_kap() ; 
      double gamma = p_eos_poly->get_gam() ;  ; 
      double kap_ns2 = pow( kappa,  0.5 /(gamma-1.) ) ; 
      
      // Polytropic unit of length in terms of r_unit : 
      double r_poly = kap_ns2 / sqrt(ggrav) ; 
      
      // Polytropic unit of time in terms of t_unit :
      double t_poly = r_poly ; 
      
      // Polytropic unit of mass in terms of m_unit :
      double m_poly = r_poly / ggrav ; 
      
      // Polytropic unit of angular momentum in terms of j_unit :
      double j_poly = r_poly * r_poly / ggrav ; 
      
      ost << "#      d [poly]       "  
	  << "       d_G [poly]     "
	  << "     Omega [poly]     "
	  << "     M_ADM [poly]     "     
	  << "       J [poly]       "  
	  << "    M_B(1) [poly]     "
	  << "    M_B(2) [poly]     " << endl ; 
      
      ost.width(20) ; 
      ost << separation() / r_poly ; ost.width(22) ;
      ost << ( star2.xa_barycenter() - star1.xa_barycenter() ) / r_poly ; ost.width(22) ; 
      ost << omega * t_poly ; ost.width(22) ;
      ost << mass_adm() / m_poly ; ost.width(22) ;
      ost << angu_mom()(2) / j_poly ; ost.width(22) ;
      ost << star1.mass_b() / m_poly ; ost.width(22) ;
      ost << star2.mass_b() / m_poly << endl ; 
      
  }
  
}
   
		    //-------------------------------//
		    //		Miscellaneous	     //
		    //-------------------------------//

double Binary_xcts::separation() const {
    
    double dx = star1.mp.get_ori_x() - star2.mp.get_ori_x() ; 
    double dy = star1.mp.get_ori_y() - star2.mp.get_ori_y() ; 
    double dz = star1.mp.get_ori_z() - star2.mp.get_ori_z() ; 
    
    return sqrt( dx*dx + dy*dy + dz*dz ) ; 
    
}
