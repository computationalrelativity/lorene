/*
 * Methods for two fluids rotating relativistic stars.
 *
 * See the file et_rot_bifluid.h for documentation
 *
 */

/*
 *   Copyright (c) 2001 Jerome Novak
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


char et_rot_bifluid_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.7  2003/11/13 12:07:57  r_prix
 * *) changed xxx2 -> Delta_car
 * *) added (non 2-fluid specific!) members sphph_euler J_euler
 * *) more or less rewritten hydro_euler() to see if I understand it ;)
 *   - somewhat simplified and more adapted to the notation used in our notes/paper.
 *   - Main difference: u_euler is no longer used!!, the "output" instead
 *     consists of ener_euler, s_euler, sphph_euler and J_euler, which are
 *     the general 3+1 components for Tmunu.
 *
 * Revision 1.6  2003/09/17 08:27:50  j_novak
 * New methods: mass_b1() and mass_b2().
 *
 * Revision 1.5  2002/10/18 08:42:58  j_novak
 * Take into account the sign for uuu and uuu2
 *
 * Revision 1.4  2002/01/16 15:03:28  j_novak
 * *** empty log message ***
 *
 * Revision 1.3  2002/01/11 14:09:34  j_novak
 * Added newtonian version for 2-fluid stars
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
 * Revision 1.3  2001/08/28  16:04:22  novak
 * Use of new definition of relative velocity and new declarations for EOS
 *
 * Revision 1.2  2001/08/27 09:58:43  novak
 * *** empty log message ***
 *
 * Revision 1.1  2001/06/22 15:39:17  novak
 * Initial revision
 *
 *
 * $Header$
 *
 */
// Headers C
#include "math.h"

// Headers Lorene
#include "et_rot_bifluid.h"
#include "utilitaires.h"

			    //--------------//
			    // Constructors //
			    //--------------//
// Standard constructor
// --------------------
Et_rot_bifluid::Et_rot_bifluid(Map& mpi, int nzet_i, bool relat, 
			      const Eos_bifluid& eos_i)
  : Etoile_rot(mpi, nzet_i, relat, *eos_i.trans2Eos()), 
  eos(eos_i),
  ent2(mpi),
  nbar2(mpi),
  Delta_car(mpi),
  gam_euler2(mpi),
  sphph_euler(mpi),
  J_euler(mpi, 1, CON, mp.get_bvect_cart()), 
  uuu2(mpi)
{
  // All the matter quantities are initialized to zero :
  nbar2 = 0 ;
  Delta_car = 0 ;
  ent2 = 0 ; 
  gam_euler2 = 1 ; 
  sphph_euler = 0;
  J_euler = 0;
  gam_euler.set_std_base() ; 
  
  // Initialization to a static state : 
  omega2 = 0 ; 
  uuu2 = 0 ; 
  
  // Pointers of derived quantities initialized to zero : 
  set_der_0x0() ;
  
}

// Copy constructor
// ----------------

Et_rot_bifluid::Et_rot_bifluid(const Et_rot_bifluid& et)
  : Etoile_rot(et), 
  eos(et.eos),
  ent2(et.ent2),
  nbar2(et.nbar2),
  Delta_car(et.Delta_car),
  gam_euler2(et.gam_euler2),
  sphph_euler(et.sphph_euler),
  J_euler(et.J_euler),
  uuu2(et.uuu2)
{
  omega2 = et.omega2 ; 
  
  // Pointers of derived quantities initialized to zero : 
  set_der_0x0() ;
}    


// Constructor from a file 
// ------------------------
Et_rot_bifluid::Et_rot_bifluid(Map& mpi, const Eos_bifluid& eos_i, FILE* fich):
  Etoile_rot(mpi, *eos_i.trans2Eos(), fich),
  eos(eos_i),
  ent2(mpi),
  nbar2(mpi),
  Delta_car(mpi),
  gam_euler2(mpi),
  sphph_euler(mpi),
  J_euler(mpi, 1, CON, mp.get_bvect_cart()), 
  uuu2(mpi) {

  // Etoile parameters
  // -----------------
  // omega2 is read in the file:     
  fread_be(&omega2, sizeof(double), 1, fich) ;		
  
  
  // Read of the saved fields:
  // ------------------------
  
  Tenseur ent2_file(mp, fich) ; 
  ent2 = ent2_file ; 
        
  // All other fields are initialized to zero : 
  // ----------------------------------------
  uuu2 = 0 ;
  Delta_car = 0 ;
  
  // Pointers of derived quantities initialized to zero 
  // --------------------------------------------------
  set_der_0x0() ;
  
}

			    //------------//
			    // Destructor //
			    //------------//

Et_rot_bifluid::~Et_rot_bifluid(){

  del_deriv() ; 

}

		//----------------------------------//
		// Management of derived quantities //
		//----------------------------------//

void Et_rot_bifluid::del_deriv() const {

  Etoile_rot::del_deriv() ; 
  
  if (p_ray_eq2 != 0x0) delete p_ray_eq2 ; 
  if (p_ray_eq2_pis2 != 0x0) delete p_ray_eq2_pis2 ; 
  if (p_ray_eq2_pi != 0x0) delete p_ray_eq2_pi ; 
  if (p_ray_pole2 != 0x0) delete p_ray_pole2 ; 
  if (p_l_surf2 != 0x0) delete p_l_surf2 ; 
  if (p_xi_surf2 != 0x0) delete p_xi_surf2 ;
  if (p_r_circ2 != 0x0) delete p_r_circ2 ;
  if (p_aplat2 != 0x0) delete p_aplat2 ; 
  if (p_mass_b1 != 0x0) delete p_mass_b1 ;
  if (p_mass_b2 != 0x0) delete p_mass_b2 ;
  
  set_der_0x0() ; 
}			    




void Et_rot_bifluid::set_der_0x0() const {

  Etoile_rot::set_der_0x0() ;
  
  p_ray_eq2 = 0x0 ;
  p_ray_eq2_pis2 = 0x0 ; 
  p_ray_eq2_pi = 0x0 ; 
  p_ray_pole2 = 0x0 ; 
  p_l_surf2 = 0x0 ; 
  p_xi_surf2 = 0x0 ; 
  p_r_circ2 = 0x0 ;
  p_aplat2 = 0x0 ;
  p_mass_b1 = 0x0;
  p_mass_b2 = 0x0;
}			    

void Et_rot_bifluid::del_hydro_euler() {

  Etoile_rot::del_hydro_euler() ; 
  gam_euler2.set_etat_nondef() ; 
  sphph_euler.set_etat_nondef();
  J_euler.set_etat_nondef();
  
  del_deriv() ; 

}			    

// Assignment of the enthalpy field
// --------------------------------

void Et_rot_bifluid::set_enthalpies(const Cmp& ent_i, const Cmp& ent2_i) {
    
  ent = ent_i ; 
  ent2 = ent2_i ;
    
  // Update of (nbar, ener, press) :
  equation_of_state() ; 
    
  // The derived quantities are obsolete:
  del_deriv() ; 
    
}


			    //--------------//
			    //  Assignment  //
			    //--------------//

// Assignment to another Et_rot_bifluid
// --------------------------------
void Et_rot_bifluid::operator=(const Et_rot_bifluid& et) {

    // Assignment of quantities common to all the derived classes of Etoile
    Etoile_rot::operator=(et) ;	    

    assert( &(et.eos) == &eos ) ;	    // Same EOS
    // Assignement of proper quantities of class Et_rot_bifluid
    omega2 = et.omega2 ; 

    ent2 = et.ent2 ;
    nbar2 = et.nbar2 ;
    Delta_car = et.Delta_car ;
    gam_euler2 = et.gam_euler2 ;
    sphph_euler = et.sphph_euler;
    J_euler = et.J_euler;
    uuu2 = et.uuu2 ;
    
    del_deriv() ;  // Deletes all derived quantities

}	

			    //--------------//
			    //	  Outputs   //
			    //--------------//

// Save in a file
// --------------
void Et_rot_bifluid::sauve(FILE* fich) const {
    
    Etoile_rot::sauve(fich) ; 
    
    fwrite_be(&omega2, sizeof(double), 1, fich) ;		
    
    ent2.sauve(fich) ; 
    
    
}

// Printing
// --------

ostream& Et_rot_bifluid::operator>>(ostream& ost) const {
    
    #include "unites.h"	    
    // To avoid some compilation warnings
    if (&ost == 0x0) {
	cout << qpig << msol << f_unit << mevpfm3 << endl ; 
    }    

    Etoile::operator>>(ost) ; 
    
    ost << endl ; 
    ost << "Bifluid rotating star" << endl ; 
    ost << "-------------" << endl ; 
    
    double freq = omega / (2.*M_PI) ;  
    ost << "Omega1 : " << omega * f_unit 
        << " rad/s     f : " << freq * f_unit << " Hz" << endl ; 
    ost << "Rotation period 1: " << 1000. / (freq * f_unit) << " ms"
	    << endl ;
       
    double freq2 = omega2 / (2.*M_PI) ;  
    ost << "Omega2 : " << omega2 * f_unit 
        << " rad/s     f : " << freq2 * f_unit << " Hz" << endl ; 
    ost << "Rotation period 2: " << 1000. / (freq2 * f_unit) << " ms"
	    << endl ;
       
    double nphi_c = nphi()(0, 0, 0, 0) ;
    if ( (omega==0) && (nphi_c==0) ) {
	 	ost << "Central N^phi :               " << nphi_c << endl ;
    }
    else{
		ost << "Central N^phi/Omega :    " << nphi_c / omega << endl ;
    }
    if (omega2!=0) 
      ost << "Central N^phi/Omega2 :     " << nphi_c / omega2 << endl ;
    
    ost << "Error on the virial identity GRV2 : " << endl ; 
    ost << "GRV2 = " << grv2() << endl ; 
    ost << "Error on the virial identity GRV3 : " << endl ; 
    double xgrv3 = grv3(&ost) ; 
    ost << "GRV3 = " << xgrv3 << endl ; 

    ost << "Circumferential equatorial radius R_circ :     " 
	<< r_circ()/km << " km" << endl ;  
    ost << "Coordinate equatorial radius r_eq : " << ray_eq()/km << " km" 
	 << endl ;  
    ost << "Flattening r_pole/r_eq :        " << aplat() << endl ; 
    ost << "Circumferential equatorial radius R_circ2 :     " 
	<< r_circ2()/km << " km" << endl ;  
    ost << "Coordinate equatorial radius r_eq2 : " << ray_eq2()/km << " km" 
	 << endl ;  
    ost << "Flattening r_pole2/r_eq2 :        " << aplat2() << endl ; 

    int lsurf = nzet - 1; 
    int nt = mp.get_mg()->get_nt(lsurf) ; 
    int nr = mp.get_mg()->get_nr(lsurf) ; 
    ost << "Equatorial value of the velocity U:         " 
	 << uuu()(nzet-1, 0, nt-1, nr-1) << " c" << endl ; 
    ost << "Equatorial value of the velocity U2:         " 
	 << uuu2()(nzet-1, 0, nt-1, nr-1) << " c" << endl ; 
    ost << "Redshift at the equator (forward) : " << z_eqf() << endl ; 
    ost << "Redshift at the equator (backward): " << z_eqb() << endl ; 
    ost << "Redshift at the pole              : " << z_pole() << endl ; 


    ost << "Central value of log(N)        : " 
	<< logn()(0, 0, 0, 0) << endl ; 

    ost << "Central value of dzeta=log(AN) : " 
	<< dzeta()(0, 0, 0, 0) << endl ; 

    ost << "Central value of B^2 : " << b_car()(0,0,0,0) <<  endl ; 

    Tbl diff_a_b = diffrel( a_car(), b_car() ) ;
    ost << 
    "Relative discrepancy in each domain between the metric coef. A^2 and B^2 : "
       << endl ;
    for (int l=0; l<diff_a_b.get_taille(); l++) {
    	ost << diff_a_b(l) << "  " ;
    }
    ost << endl ;    	

    return ost ;
    
}


void Et_rot_bifluid::partial_display(ostream& ost) const {
    
    #include "unites.h"	    
    // To avoid some compilation warnings
    if (&ost == 0x0) {
	cout << qpig << msol << f_unit << mevpfm3 << endl ; 
    }    

    double freq = omega / (2.*M_PI) ;  
    ost << "Omega : " << omega * f_unit 
        << " rad/s     f : " << freq * f_unit << " Hz" << endl ; 
    ost << "Rotation period : " << 1000. / (freq * f_unit) << " ms"
	 << endl ;
    ost << endl << "Central enthalpy : " << ent()(0,0,0,0) << " c^2" << endl ; 
    ost << "Central proper baryon density : " << nbar()(0,0,0,0) 
	<< " x 0.1 fm^-3" << endl ; 
    double freq2 = omega2 / (2.*M_PI) ;  
    ost << "Omega2 : " << omega2 * f_unit 
        << " rad/s     f : " << freq2 * f_unit << " Hz" << endl ; 
    ost << "Rotation period 2: " << 1000. / (freq2 * f_unit) << " ms"
	 << endl ;
    ost << endl << "Central enthalpy 2: " << ent2()(0,0,0,0) << " c^2" << endl ; 
    ost << "Central proper baryon density 2: " << nbar2()(0,0,0,0) 
	<< " x 0.1 fm^-3" << endl ; 
   ost << "Central proper energy density : " << ener()(0,0,0,0) 
	<< " rho_nuc c^2" << endl ; 
    ost << "Central pressure : " << press()(0,0,0,0) 
	<< " rho_nuc c^2" << endl ; 

    ost << "Central value of log(N)        : " 
	<< logn()(0, 0, 0, 0) << endl ; 
    ost << "Central lapse N :      " << nnn()(0,0,0,0) <<  endl ; 
    ost << "Central value of dzeta=log(AN) : " 
	<< dzeta()(0, 0, 0, 0) << endl ; 
    ost << "Central value of A^2 : " << a_car()(0,0,0,0) <<  endl ; 
    ost << "Central value of B^2 : " << b_car()(0,0,0,0) <<  endl ; 

    double nphi_c = nphi()(0, 0, 0, 0) ;
    if ( (omega==0) && (nphi_c==0) ) {
		ost << "Central N^phi :               " << nphi_c << endl ;
    }
    else{
		ost << "Central N^phi/Omega :         " << nphi_c / omega << endl ;
    }
	    

    int lsurf = nzet - 1; 
    int nt = mp.get_mg()->get_nt(lsurf) ; 
    int nr = mp.get_mg()->get_nr(lsurf) ; 
    ost << "Equatorial value of the velocity U:         " 
	 << uuu()(nzet-1, 0, nt-1, nr-1) << " c" << endl ; 

    ost << "Equatorial value of the velocity U2:         " 
	 << uuu2()(nzet-1, 0, nt-1, nr-1) << " c" << endl ; 

    ost << endl 
	<< "Coordinate equatorial radius r_eq =    " 
	<< ray_eq()/km << " km" << endl ;  
    ost << "Flattening r_pole/r_eq :        " << aplat() << endl ; 
    ost << endl 
	<< "Coordinate equatorial radius r_eq2 =    " 
	<< ray_eq2()/km << " km" << endl ;  
    ost << "Flattening r_pole2/r_eq2 :        " << aplat2() << endl ; 

}

//
//   Equation of state
//

void Et_rot_bifluid::equation_of_state() {

  Cmp ent_eos = ent() ;
  Cmp ent2_eos = ent2() ;
  Tenseur rel_vel(Delta_car) ;

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

    if (nzet > 2) {
      
      cout << "Et_rot_bifluid::equation_of_state: not ready yet for nzet > 2 !"
	   << endl ;    	
    }
    
    ent_eos = fact_ent * ent_eos ;
    ent2_eos = fact_ent * ent2_eos ;
    ent_eos.std_base_scal() ;
    ent2_eos.std_base_scal() ;
  }
  
  
  // Call to the EOS
  nbar.set_etat_qcq() ;
  nbar2.set_etat_qcq() ;
  ener.set_etat_qcq() ;
  press.set_etat_qcq() ;

  eos.calcule_tout(ent_eos, ent2_eos, rel_vel(), nbar.set(), nbar2.set(), 
  		   ener.set(), press.set(), nzet)  ; 
  
  // Set the bases for spectral expansion 
  nbar.set_std_base() ; 
  nbar2.set_std_base() ; 
  ener.set_std_base() ; 
  press.set_std_base() ; 
  
  // The derived quantities are obsolete
  del_deriv() ; 
  
}

//
// Computation of hydro quantities
//

void Et_rot_bifluid::hydro_euler(){

    int nz = mp.get_mg()->get_nzone() ; 
    int nzm1 = nz - 1 ; 

    // RP: I prefer to use the 3-vector J_euler instead of u_euler
    // for better physical "encapsulation" 
    // (i.e. --> use same form of Poisson-equations for all etoile sub-classes!)
    u_euler.set_etat_nondef(); // make sure it's not used

    // (Norm of) Euler-velocity of the first fluid
    //------------------------------
    uuu.set_etat_qcq();

    uuu.set() = bbb() * (omega - nphi() ) / nnn();
    uuu.annule(nzm1) ; 

    // gosh, we have to exclude the thing being zero here... :(
    if( uuu.get_etat() != ETATZERO )
      {
	(uuu.set()).std_base_scal() ;
	(uuu.set()).mult_rsint();
      }
    uuu.set_std_base();
    

    // (Norm of) Euler-velocity of the second fluid
    //----------------------------------------
    uuu2.set_etat_qcq();
    
    uuu2.set() = bbb() * (omega2 - nphi() ) / nnn();
    uuu2.annule(nzm1) ; 

    if( uuu2.get_etat() != ETATZERO )
      {
	(uuu2.set()).std_base_scal();
	(uuu2.set()).mult_rsint();
      }
    uuu2.set_std_base();

    Tenseur uuu_car = uuu * uuu;
    Tenseur uuu2_car = uuu2 * uuu2;

    // Lorentz factors
    // --------------
    Tenseur gam_car = 1.0 / (1.0 - unsurc2 * uuu_car) ; 
    gam_euler = sqrt(gam_car) ; 
    gam_euler.set_std_base() ;  // sets the standard spectral bases for a scalar field

    Tenseur gam2_car = 1.0 /(1.0 - unsurc2 * uuu2_car) ;
    gam_euler2 = sqrt(gam2_car) ;
    gam_euler2.set_std_base() ;

    // Update of "relative velocity" squared: $\Delta^2$
    // ---------------------------
    Delta_car = (uuu - uuu2)*(uuu - uuu2) / ( (1 - unsurc2* uuu*uuu2) *(1 - unsurc2* uuu*uuu2 ) ) ;
    Delta_car.set_std_base() ;

    // some auxiliary EOS variables
    Tenseur Ann(ent) ;
    Ann = gam_car()*nbar()*nbar() * eos.get_Knn(nbar(), nbar2(), Delta_car(), nzet);
    Tenseur Anp(ent) ;
    Anp = gam_euler()*gam_euler2()* nbar()*nbar2()* eos.get_Knp(nbar(),nbar2(),Delta_car(),nzet);
    Tenseur App(ent) ;
    App = gam2_car()*nbar2()*nbar2() *eos.get_Kpp(nbar(), nbar2(), Delta_car(), nzet);
  

    //  Energy density E with respect to the Eulerian observer
    //------------------------------------
    ener_euler =  Ann + 2*Anp + App - press ;
    ener_euler.set_std_base() ; 
    

    // S^phi_phi component of stress-tensor S^i_j
    //------------------------------------
    sphph_euler = press() + Ann()*uuu_car() + 2*Anp()*uuu()*uuu2() + App()* uuu2_car();
    sphph_euler.set_std_base();


    // Trace of the stress tensor with respect to the Eulerian observer
    //------------------------------------
    s_euler = 2*press() + sphph_euler();
    s_euler.set_std_base() ; 


    // the (flat-space) angular-momentum 3-vector J_euler^i
    //-----------------------------------
    Tenseur Jph(mp);   // the normalized phi-component of J^i: Sqrt[g_phiphi]*J^phi
    Jph = Ann*uuu + Anp*(uuu + uuu2) + App*uuu2 ;

    J_euler.set_etat_qcq();

    J_euler.set(0) = 0;		// r tetrad component
    J_euler.set(1) = 0;		// theta tetrad component
    J_euler.set(2) = Jph()/ bbb(); // phi tetrad component ... = J^phi r sin(th)
    J_euler.set_triad (mp.get_bvect_spher());
    J_euler.set_std_base();

    // RP: it seems that J_euler _HAS_ to have cartesian triad set on exit from here...!!
    J_euler.change_triad( mp.get_bvect_cart() ) ;	// Triad = Cartesian triad

    if( (J_euler(0).get_etat() == ETATZERO)&&(J_euler(1).get_etat() == ETATZERO)&&(J_euler(2).get_etat()==ETATZERO))
      J_euler = 0;

    // The derived quantities are obsolete
    // -----------------------------------
    del_deriv() ;                

} // hydro_euler()


