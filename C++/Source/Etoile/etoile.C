/*
 * Methods of class Etoile
 *
 * (see file etoile.h for documentation)
 */

/*
 *   Copyright (c) 2000-2001 Eric Gourgoulhon
 *   Copyright (c) 2000-2001 Keisuke Taniguchi
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
 * Revision 1.11  2016/12/05 16:17:54  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.10  2014/10/13 08:52:58  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.9  2012/08/12 17:48:35  p_cerda
 * Magnetstar: New classes for magnetstar. Allowing for non-equatorial symmetry in Etoile et al. Adding B_phi in Et_rot_mag.
 *
 * Revision 1.8  2005/01/18 22:36:50  k_taniguchi
 * Delete a pointer for ray_eq(int kk).
 *
 * Revision 1.7  2005/01/18 20:35:05  k_taniguchi
 * Addition of ray_eq(int kk).
 *
 * Revision 1.6  2005/01/17 20:40:25  k_taniguchi
 * Addition of ray_eq_3pis2().
 *
 * Revision 1.5  2004/03/25 10:29:06  j_novak
 * All LORENE's units are now defined in the namespace Unites (in file unites.h).
 *
 * Revision 1.4  2003/10/13 15:23:56  f_limousin
 * *** empty log message ***
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
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 2.14  2000/11/24  13:27:44  eric
 * Dans eqution_of_state(): changement leger de ent dans le cas ou l'on a
 * deux domaine avant d'appeler l'EOS.
 *
 * Revision 2.13  2000/09/25  12:22:02  keisuke
 * *** empty log message ***
 *
 * Revision 2.12  2000/09/22  15:50:58  keisuke
 * Ajout du membre d_logn_auto_div.
 *
 * Revision 2.11  2000/09/07  14:34:09  keisuke
 * Ajout du membre logn_auto_regu.
 *
 * Revision 2.10  2000/08/31  15:36:54  eric
 * Bases spectrales standards pour nnn, a_car et gam_euler dans le
 * constructeur (initialisation a la metrique plate).
 *
 * Revision 2.9  2000/08/29  11:37:49  eric
 * Ajout des membres k_div et logn_auto_div.
 *
 * Revision 2.8  2000/07/21  12:01:11  eric
 * Modif dans Etoile::del_deriv() :
 *   appel de Etoile::set_der_0x0() et non de la fonction virtuelle set_der_0x0().
 *
 * Revision 2.7  2000/03/21  12:39:34  eric
 * Le constructeur standard teste la compatibilite de l'EOS avec le
 * caractere relativiste de l'etoile.
 *
 * Revision 2.6  2000/02/21  14:32:40  eric
 * gam_euler est initialise a 1 dans le constructeur standard.
 * Suppression de l'appel a del_hydro_euler dans equation_of_state().
 *
 * Revision 2.5  2000/02/09  19:30:47  eric
 * La triade de decomposition doit desormais figurer en argument des
 *  constructeurs de Tenseur.
 *
 * Revision 2.4  2000/02/02  09:23:34  eric
 * Affichage de la masse.
 *
 * Revision 2.3  2000/01/28  17:18:10  eric
 * Modifs noms des quantites globales.
 * Affichage.
 *
 * Revision 2.2  2000/01/24  17:13:36  eric
 * Le mapping mp n'est plus constant.
 *
 * Revision 2.1  2000/01/24  13:37:22  eric
 * *** empty log message ***
 *
 * Revision 2.0  2000/01/20  17:04:45  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers C
#include "math.h"

// Headers Lorene
#include "etoile.h"
#include "eos.h"
#include "utilitaires.h"
#include "param.h"
#include "unites.h"

			    //--------------//
			    // Constructors //
			    //--------------//

// Standard constructor
// --------------------
namespace Lorene {
Etoile::Etoile(Map& mpi, int nzet_i, bool relat, const Eos& eos_i)
		 : mp(mpi), 
		   nzet(nzet_i), 
		   relativistic(relat), 
		   k_div(0), 
		   eos(eos_i), 
		   ent(mpi), 
		   nbar(mpi), 
		   ener(mpi), 
		   press(mpi),  
		   ener_euler(mpi), 
		   s_euler(mpi), 
		   gam_euler(mpi), 
		   u_euler(mpi, 1, CON, mp.get_bvect_cart()), 
		   logn_auto(mpi), 
		   logn_auto_regu(mpi), 
		   logn_auto_div(mpi),
		   d_logn_auto_div(mpi, 1, COV, mp.get_bvect_spher()), 
		   beta_auto(mpi), 
		   nnn(mpi), 
		   shift(mpi, 1, CON, mp.get_bvect_cart()),
		   a_car(mpi) {

    // Check of the EOS
    const Eos_poly* p_eos_poly = dynamic_cast<const Eos_poly*>( &eos ) ; 	  
    const Eos_poly_newt* p_eos_poly_newt = 
			    dynamic_cast<const Eos_poly_newt*>( &eos ) ; 	  
    const Eos_incomp* p_eos_incomp = dynamic_cast<const Eos_incomp*>( &eos ) ; 	  
    const Eos_incomp_newt* p_eos_incomp_newt = 
			    dynamic_cast<const Eos_incomp_newt*>( &eos ) ; 	  

    if (relativistic) {

	if (p_eos_poly_newt != 0x0) {
	    cout << 
	    "Etoile::Etoile : the EOS Eos_poly_newt must not be employed"
		<< " for a relativistic star ! " << endl ; 
	    cout << "(Use Eos_poly instead)" << endl ; 
	    abort() ; 
	}
	if (p_eos_incomp_newt != 0x0) {
	    cout << 
	    "Etoile::Etoile : the EOS Eos_incomp_newt must not be employed"
		<< " for a relativistic star ! " << endl ; 
	    cout << "(Use Eos_incomp instead)" << endl ; 
	    abort() ; 
	}

    }
    else{

	if ( (p_eos_poly != 0x0) && (p_eos_poly_newt == 0x0) ) {
	    cout << 
	    "Etoile::Etoile : the EOS Eos_poly must not be employed"
		<< " for a Newtonian star ! " << endl ; 
	    cout << "(Use Eos_poly_newt instead)" << endl ; 
	    abort() ; 
	}
	if ( (p_eos_incomp != 0x0) && (p_eos_incomp_newt == 0x0) ) {
	    cout << 
	    "Etoile::Etoile : the EOS Eos_incomp must not be employed"
		<< " for a relativistic star ! " << endl ; 
	    cout << "(Use Eos_incomp_newt instead)" << endl ; 
	    abort() ; 
	}
	
    }


    // Parameter 1/c^2      
    unsurc2 = relativistic ? double(1) : double(0) ; 

    // Pointers of derived quantities initialized to zero : 
    set_der_0x0() ;

    // All the matter quantities are initialized to zero :
    nbar = 0 ; 
    ener = 0 ; 
    press = 0 ; 
    ent = 0 ; 
    ener_euler = 0 ; 
    s_euler = 0 ; 
    gam_euler = 1 ; 
    gam_euler.set_std_base() ; 
    u_euler = 0 ; 

    // The metric is initialized to the flat one : 
    logn_auto = 0 ; 
    logn_auto_regu = 0 ; 
    logn_auto_div = 0 ; 
    d_logn_auto_div = 0 ; 
    beta_auto = 0 ; 
    nnn = 1 ; 
    nnn.set_std_base() ; 
    shift = 0 ; 
    a_car = 1 ; 
    a_car.set_std_base() ; 
    
}

// Copy constructor
// ----------------
Etoile::Etoile(const Etoile& et) 
		 : mp(et.mp), 
		   nzet(et.nzet), 
		   relativistic(et.relativistic), 
		   unsurc2(et.unsurc2), 
		   k_div(et.k_div), 
		   eos(et.eos), 
		   ent(et.ent), 
		   nbar(et.nbar), 
		   ener(et.ener), 
		   press(et.press),  
		   ener_euler(et.ener_euler), 
		   s_euler(et.s_euler), 
		   gam_euler(et.gam_euler), 
		   u_euler(et.u_euler), 
		   logn_auto(et.logn_auto), 
		   logn_auto_regu(et.logn_auto_regu),
		   logn_auto_div(et.logn_auto_div), 
		   d_logn_auto_div(et.d_logn_auto_div), 
		   beta_auto(et.beta_auto), 
		   nnn(et.nnn), 
		   shift(et.shift), 
		   a_car(et.a_car) {
	       
    set_der_0x0() ;

}

// Constructor from a file
// -----------------------
Etoile::Etoile(Map& mpi, const Eos& eos_i, FILE* fich)
		 : mp(mpi), 
		   eos(eos_i), 
		   ent(mpi), 
		   nbar(mpi), 
		   ener(mpi), 
		   press(mpi),  
		   ener_euler(mpi), 
		   s_euler(mpi), 
		   gam_euler(mpi), 
		   u_euler(mpi, 1, CON, mp.get_bvect_cart()), 
		   logn_auto(mpi), 
		   logn_auto_regu(mpi), 
		   logn_auto_div(mpi),
		   d_logn_auto_div(mpi, 1, COV, mp.get_bvect_spher()), 
		   beta_auto(mpi), 
		   nnn(mpi), 
		   shift(mpi, 1, CON, mp.get_bvect_cart()),
		   a_car(mpi) {

    // Etoile parameters
    // -----------------

    // nzet and relativistic are read in the file:     
    int xx ; 
    fread_be(&xx, sizeof(int), 1, fich) ;	
    k_div = xx / 1000 ;	    // integer part
    nzet = xx - k_div * 1000  ;
    		
    fread(&relativistic, sizeof(bool), 1, fich) ;		
    	  
    // Parameter 1/c^2 is deduced from relativistic:
    unsurc2 = relativistic ? double(1) : double(0) ; 


    // Equation of state
    // -----------------
    
    // Read of the saved EOS
    Eos* p_eos_file = Eos::eos_from_file(fich) ; 
    
    // Comparison with the assigned EOS:
    if (eos != *p_eos_file) {
	cout << 
	"Etoile::Etoile(const Map&, const Eos&, FILE*) : the EOS given in "
	<< endl << 
	" argument and that read in the file are different !" << endl ; 
	abort() ;  
    }
    
    // p_eos_file is no longer required (it was used only for checking the
    //  EOS compatibility)
    delete p_eos_file ;
    
    // Read of the saved fields:
    // ------------------------
    Tenseur ent_file(mp, fich) ; 
    ent = ent_file ; 
    
    Tenseur logn_auto_file(mp, fich) ; 
    logn_auto = logn_auto_file ; 
    
    Tenseur beta_auto_file(mp, fich) ; 
    beta_auto = beta_auto_file ; 

    if (k_div == 0) {
	logn_auto_div = 0 ; 
	d_logn_auto_div = 0 ; 
    }
    else {

	Tenseur logn_auto_div_file(mp, fich) ; 
	logn_auto_div = logn_auto_div_file ; 

	Tenseur d_logn_auto_div_file(mp, mp.get_bvect_spher(), fich) ; 
	d_logn_auto_div = d_logn_auto_div_file ; 
    }

    logn_auto_regu = logn_auto - logn_auto_div ;

     shift = 0 ;    

    // Pointers of derived quantities initialized to zero 
    // --------------------------------------------------
    set_der_0x0() ;
    
}

			    //------------//
			    // Destructor //
			    //------------//

Etoile::~Etoile(){

    del_deriv() ; 

}


			//----------------------------------//
			// Management of derived quantities //
			//----------------------------------//

void Etoile::del_deriv() const {

    if (p_mass_b != 0x0) delete p_mass_b ; 
    if (p_mass_g != 0x0) delete p_mass_g ; 
    if (p_ray_eq != 0x0) delete p_ray_eq ; 
    if (p_ray_eq_pis2 != 0x0) delete p_ray_eq_pis2 ; 
    if (p_ray_eq_pi != 0x0) delete p_ray_eq_pi ; 
    if (p_ray_eq_3pis2 != 0x0) delete p_ray_eq_3pis2 ; 
    if (p_ray_pole != 0x0) delete p_ray_pole ; 
    if (p_l_surf != 0x0) delete p_l_surf ; 
    if (p_xi_surf != 0x0) delete p_xi_surf ; 

    Etoile::set_der_0x0() ; 
}			    




void Etoile::set_der_0x0() const {

    p_mass_b = 0x0 ; 
    p_mass_g = 0x0 ; 
    p_ray_eq = 0x0 ; 
    p_ray_eq_pis2 = 0x0 ; 
    p_ray_eq_pi = 0x0 ; 
    p_ray_eq_3pis2 = 0x0 ; 
    p_ray_pole = 0x0 ; 
    p_l_surf = 0x0 ; 
    p_xi_surf = 0x0 ; 

}			    

void Etoile::del_hydro_euler() {

    ener_euler.set_etat_nondef() ; 
    s_euler.set_etat_nondef() ; 
    gam_euler.set_etat_nondef() ; 
    u_euler.set_etat_nondef() ; 

    del_deriv() ; 

}			    




			    //--------------//
			    //  Assignment  //
			    //--------------//

// Assignment to another Etoile
// ----------------------------
void Etoile::operator=(const Etoile& et) {

    assert( &(et.mp) == &mp ) ;		    // Same mapping
    assert( &(et.eos) == &eos ) ;	    // Same EOS
    
    nzet = et.nzet ; 
    relativistic = et.relativistic ; 
    k_div = et.k_div ; 
    unsurc2 = et.unsurc2 ; 
    
    ent = et.ent ;
    nbar = et.nbar ; 
    ener = et.ener ;
    press = et.press ;
    ener_euler = et.ener_euler ;
    s_euler = et.s_euler ;
    gam_euler = et.gam_euler ;
    u_euler = et.u_euler ;
    logn_auto = et.logn_auto ;
    logn_auto_regu = et.logn_auto_regu ;
    logn_auto_div = et.logn_auto_div ;
    d_logn_auto_div = et.d_logn_auto_div ;
    beta_auto = et.beta_auto ;
    nnn = et.nnn ;
    shift = et.shift ;
    a_car = et.a_car ;
    
    
    del_deriv() ;  // Deletes all derived quantities

}	

// Assignment of the enthalpy field
// --------------------------------

void Etoile::set_enthalpy(const Cmp& ent_i) {
    
    ent = ent_i ; 
    
    // Update of (nbar, ener, press) :
    equation_of_state() ; 
    
    // The derived quantities are obsolete:
    del_deriv() ; 
    
}

			    //--------------//
			    //	  Outputs   //
			    //--------------//

// Save in a file
// --------------
void Etoile::sauve(FILE* fich) const {
    
    int xx = nzet + k_div * 1000 ;     
    fwrite_be(&xx, sizeof(int), 1, fich) ;			

    fwrite(&relativistic, sizeof(bool), 1, fich) ;		

    eos.sauve(fich) ; 
    
    ent.sauve(fich) ;     
    logn_auto.sauve(fich) ;     
    beta_auto.sauve(fich) ;  
    
    if (k_div != 0) {
	logn_auto_div.sauve(fich) ; 
	d_logn_auto_div.sauve(fich) ; 
    }   
    
}

// Printing
// --------

ostream& operator<<(ostream& ost, const Etoile& et)  {
    et >> ost ;
    return ost ;
}
    
ostream& Etoile::operator>>(ostream& ost) const {
    
  using namespace Unites ;

    ost << endl ; 
    if (relativistic) {
	ost << "Relativistic star" << endl ; 
	ost << "-----------------" << endl ; 
    }
    else {
	ost << "Newtonian star" << endl ; 
	ost << "--------------" << endl ; 
    }
    
    ost << "Number of domains occupied by the star : " << nzet << endl ; 
    
    ost << "Equation of state : " << endl ; 
    ost << eos << endl ; 
    
    ost << endl << "Central enthalpy : " << ent()(0,0,0,0) << " c^2" << endl ; 
    ost << "Central proper baryon density : " << nbar()(0,0,0,0) 
	<< " x 0.1 fm^-3" << endl ; 
    ost << "Central proper energy density : " << ener()(0,0,0,0) 
	<< " rho_nuc c^2" << endl ; 
    ost << "Central pressure : " << press()(0,0,0,0) 
	<< " rho_nuc c^2" << endl ; 

    ost << endl 
	<< "Regularization index of the gravitational potential : k_div = "
	<< k_div << endl ; 
    ost << "Central lapse N :      " << nnn()(0,0,0,0) <<  endl ; 
    ost <<	   "Central value of A^2 : " << a_car()(0,0,0,0) <<  endl ; 

    ost << endl 
	<< "Coordinate equatorial radius (phi=0) a1 =    " 
	<< ray_eq()/km << " km" << endl ;  
    ost << "Coordinate equatorial radius (phi=pi/2) a2 = " 
	<< ray_eq_pis2()/km << " km" << endl ;  
    ost << "Coordinate equatorial radius (phi=pi):       " 
	<< ray_eq_pi()/km << " km" << endl ;  
    ost << "Coordinate polar radius a3 =                 " 
	<< ray_pole()/km << " km" << endl ;  
    ost << "Axis ratio a2/a1 = " << ray_eq_pis2() / ray_eq() 
	<< "  a3/a1 = " << ray_pole() / ray_eq() << endl ; 	

    ost << endl << "Baryon mass :        " << mass_b() / msol << " M_sol" << endl ; 
    ost << "Gravitational mass : " << mass_g() / msol << " M_sol" << endl ; 
    
    return ost ; 
}

		//-----------------------------------------//
		//	Computation of hydro quantities	   //
		//-----------------------------------------//

void Etoile::equation_of_state() {

	Cmp ent_eos = ent() ;


    // Slight rescale of the enthalpy field in case of 2 domains inside the
    //  star


        double epsilon = 1.e-12 ;

	const Mg3d* mg = mp.get_mg() ;
        int nz = mg->get_nzone() ;

        Mtbl xi(mg) ;
        xi.set_etat_qcq() ;
        for (int l=0; l<nz; l++) {
        	xi.t[l]->set_etat_qcq() ;
        	for (int k=0; k<mg->get_np(l); k++) {
        		for (int j=0; j<mg->get_nt(l); j++) {
        			for (int i=0; i<mg->get_nr(l); i++) {
        				xi.set(l,k,j,i) =
        					mg->get_grille3d(l)->x[i] ;
        			}
        		}
        	}

        }

     	Cmp fact_ent(mp) ;
     	fact_ent.allocate_all() ;
     	
     	fact_ent.set(0) = 1 + epsilon * xi(0) * xi(0) ;
     	fact_ent.set(1) = 1 - 0.25 * epsilon * (xi(1) - 1) * (xi(1) - 1) ;
     	
     	for (int l=nzet; l<nz; l++) {
     		fact_ent.set(l) = 1 ;
     	}

    if (nzet > 1) {

      if(nzet == 3) {
	fact_ent.set(1) = 1 - 0.5 * epsilon * (xi(1) - 0.5) * (xi(1) - 0.5) ;
	fact_ent.set(2) = 1 - 0.25 * epsilon * (xi(2) - 1) * (xi(2) - 1) ;
      }
      
      if (nzet > 3) {
    	
    		cout << "Etoile::equation_of_state: not ready yet for nzet > 3 !"
    		     << endl ;    	
    	}

    	ent_eos = fact_ent * ent_eos ;
    	ent_eos.std_base_scal() ;
    }





    // Call to the EOS (the EOS is called domain by domain in order to
    //          allow for the use of MEos)

    Cmp tempo(mp) ;

    nbar.set_etat_qcq() ;
    nbar.set() = 0 ;
    for (int l=0; l<nzet; l++) {

        Param par ;       // Paramater for multi-domain equation of state
        par.add_int(l) ;

        tempo =  eos.nbar_ent(ent_eos, 1, l, &par) ;

        nbar = nbar() + tempo ;

    }

    ener.set_etat_qcq() ;
    ener.set() = 0 ;
    for (int l=0; l<nzet; l++) {

        Param par ;    // Paramater for multi-domain equation of state
        par.add_int(l) ;

        tempo =  eos.ener_ent(ent_eos, 1, l, &par) ;

        ener = ener() + tempo ;

    }

    press.set_etat_qcq() ;
    press.set() = 0 ;
    for (int l=0; l<nzet; l++) {

        Param par ;     // Paramater for multi-domain equation of state
        par.add_int(l) ;

        tempo =  eos.press_ent(ent_eos, 1, l, &par) ;

        press = press() + tempo ;

    }


    // Set the bases for spectral expansion
    nbar.set_std_base() ; 
    ener.set_std_base() ; 
    press.set_std_base() ; 

    // The Eulerian quantities are obsolete
    //## del_hydro_euler() ; 
    
    // The derived quantities are obsolete
    del_deriv() ; 
    
}

void Etoile::hydro_euler() {
    
    cout << 
    "Etoile::hydro_euler : hydro_euler must be called via a derived class"
    << endl << " of Etoile !" << endl ; 
    
    abort() ;        
    
}
}
