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
  : Etoile_rot(mpi, nzet_i, relat, *eos_i.trans2Eos(relat)), 
  eos(eos_i),
  ent2(mpi),
  nbar2(mpi),
  xxx2(mpi),
  gam_euler2(mpi),
  uuu2(mpi)
  
{   //Here, there should be a check of the EOS...

    // All the matter quantities are initialized to zero :
    nbar2 = 0 ;
    xxx2 = 0 ;
    ent2 = 0 ; 
    gam_euler2 = 1 ; 
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
  xxx2(et.xxx2),
  gam_euler2(et.gam_euler2),
  uuu2(et.uuu2)
{
    omega2 = et.omega2 ; 

    // Pointers of derived quantities initialized to zero : 
    set_der_0x0() ;
}    


// Constructor from a file (works only for relativistic stars!)
// ------------------------------------------------------------
Et_rot_bifluid::Et_rot_bifluid(Map& mpi, const Eos_bifluid& eos_i, FILE* fich):
  Etoile_rot(mpi, *eos_i.trans2Eos(true), fich),
  eos(eos_i),
  ent2(mpi),
  nbar2(mpi),
  xxx2(mpi),
  gam_euler2(mpi),
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
}			    

void Et_rot_bifluid::del_hydro_euler() {

    Etoile_rot::del_hydro_euler() ; 
    gam_euler2.set_etat_nondef() ; 

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
    xxx2 = et.xxx2 ;
    gam_euler2 = et.gam_euler2 ;
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
    ost << "Omega : " << omega * f_unit 
        << " rad/s     f : " << freq * f_unit << " Hz" << endl ; 
    ost << "Rotation period : " << 1000. / (freq * f_unit) << " ms"
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
		ost << "Central N^phi/Omega2 :         " << nphi_c / omega2 << endl ;
    }
	    
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
  Tenseur rel_vel(*ent.get_mp()) ;
  if ((uuu.get_etat() == ETATNONDEF)||(uuu2.get_etat() == ETATNONDEF)) {
    rel_vel = 0 ; }
  else
    rel_vel = (uuu - uuu2)*(uuu - uuu2) / ((1 - uuu*uuu2)*(1-uuu*uuu2)) ;
  
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
  
  xxx2 = rel_vel ;
  
  // Set the bases for spectral expansion 
  nbar.set_std_base() ; 
  nbar2.set_std_base() ; 
  xxx2.set_std_base() ; 
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

    // Computation of u_euler
    // WARNING!!: at the end of this routine u_euler is NOT the same 
    // type of quantity as in Etoile_rot, it represents the (normalized) 
    // components of the momentum 3-vector (3+1 decomposition of $T^\mu_\nu$). 
    // ---------------------------------------------------------------------
    
    const Coord& x = mp.x ; 
    const Coord& y = mp.y ; 
    
    u_euler.set_etat_qcq() ; 
    
    u_euler.set(0) = - omega * y ;	    // Cartesian components of solid rotation
    u_euler.set(1) =   omega * x ;
    u_euler.set(2) = 0 ;
    u_euler.annule(nzm1) ; 
    
    u_euler.set_triad( mp.get_bvect_cart() ) ;	// Triad = Cartesian triad
    
    u_euler.set_std_base() ;	// sets the standard bases for spectral expansions

    u_euler = ( u_euler - shift ) / nnn ; 

    u_euler.set_std_base() ;	// sets the standard bases for spectral expansions
    // At this point u_euler is the equivalent of Etoile_rot.u_euler
    
//## Test
    Tenseur utest(mp, 1, CON, mp.get_bvect_spher()) ; 
    utest.set_etat_qcq() ; 
    
    utest.set(0) = 0 ;	    // Spherical components of solid rotation
    utest.set(1) = 0 ;
    utest.set(2) = ( omega - nphi() ) / nnn();

    utest.set(2).annule(nzm1) ; 
    utest.set(2).std_base_scal() ;
    utest.set(2).mult_rsint() ;	    //  Multiplication by r sin(theta)
    
    utest.set_triad( mp.get_bvect_spher() ) ; 

    utest.change_triad( mp.get_bvect_cart() ) ; 
    
    for (int i=0; i<3; i++) {
	Valeur& uu = u_euler.set(i).va ;
	Valeur& ut = utest.set(i).va ;
	
	if (uu.get_etat() != ETATZERO) {
	    uu.coef() ; 
	    
	    if (ut.get_etat() == ETATZERO) {
		ut.set_etat_cf_qcq() ; 
		*(ut.c_cf) = 0 ; 
		ut.c_cf->base = uu.c_cf->base ; 
	    }
	    else {
		ut.coef() ; 
	    }
	    
	    Mtbl_cf diff = *(uu.c_cf) - *(ut.c_cf) ;
	    cout << "Etoile_rot::hydro_euler: test u_euler(" << i << ") : " 
		 << max( abs(diff) )(0) << endl ; 
	
	}
    }
//##

    if ( (u_euler(0).get_etat() == ETATZERO) &&
	 (u_euler(1).get_etat() == ETATZERO) &&
	 (u_euler(2).get_etat() == ETATZERO) )    {
	
	u_euler = 0 ;    
    }



    // Computation of uuu (norme of u_euler)
    // ------------------

    // The scalar product is performed on the spherical components: 

    Tenseur us = u_euler ; 
    us.change_triad( mp.get_bvect_spher() ) ; 

    Cmp uu2 =	a_car() * ( us(0) * us(0) + us(1) * us(1) ) 
	     +	b_car() * us(2) * us(2) ; 

    uuu = sqrt( uu2 ) ; 
    
    if (uuu.get_etat() == ETATQCQ) {
	((uuu.set()).va).set_base( us(2).va.base ) ;   // Same basis as 
    }						   // (Omega -N^phi) r sin(theta)

    us.set_etat_qcq() ; 
    
    // Computation of uuu2 (without test!)
    // ------------------------------------------------
    
    us.set(0) = - omega2 * y ;	    // Cartesian components of solid rotation
    us.set(1) =   omega2 * x ;
    us.set(2) = 0 ;
    us.annule(nzm1) ; 
    
    us.set_triad( mp.get_bvect_cart() ) ;	// Triad = Cartesian triad
    
    us.set_std_base() ;	// sets the standard bases for spectral expansions

    us = ( us - shift ) / nnn ; 

    us.set_std_base() ;	// sets the standard bases for spectral expansions 
    if ( (us(0).get_etat() == ETATZERO) &&
	 (us(1).get_etat() == ETATZERO) &&
	 (us(2).get_etat() == ETATZERO) )    {
	
	us = 0 ;    
    }

    Tenseur u_euler2 = us ;
    us.change_triad( mp.get_bvect_spher() ) ; 

    Cmp uu22 =	a_car() * ( us(0) * us(0) + us(1) * us(1) ) 
	     +	b_car() * us(2) * us(2) ; 

    uuu2 = sqrt( uu22 ) ; 
    
    if (uuu2.get_etat() == ETATQCQ) {
	((uuu2.set()).va).set_base( us(2).va.base ) ;   // Same basis as 
    }						   // (Omega -N^phi) r sin(theta)


    // Lorentz factor
    // --------------
    
    Tenseur u2(mp) ; 
    u2 = unsurc2 * uu2 ; 
    
    Tenseur gam2 = 1 / (1 - u2) ; 
    
    gam_euler = sqrt(gam2) ; 

    gam_euler.set_std_base() ;  // sets the standard spectral bases for
				    // a scalar field

    // Second fluid Lorentz factor
    // ---------------------------
    u2 = unsurc2 * uu22 ;

    Tenseur gam22 = 1 /(1 - u2) ;

    gam_euler2 = sqrt(gam22) ;

    gam_euler2.set_std_base() ;


    //  Energy density E with respect to the Eulerian observer
    //------------------------------------
    
    Tenseur Ann(ent) ;
    Ann = gam2()*eos.get_Knn(nbar(), nbar2(), xxx2(), nzet)*nbar()*nbar() ;
    Tenseur Anp(ent) ;
    Anp = 2*gam_euler()*gam_euler2()*eos.get_Knp(nbar(), nbar2(), 
				      xxx2(),nzet)*nbar()*nbar2() ;
    Tenseur App(ent) ;
    App = gam22()*eos.get_Kpp(nbar(), nbar2(), xxx2(), nzet)
    *nbar2()*nbar2() ;
    
    ener_euler =  Ann + Anp + App - press ;
    ener_euler.set_std_base() ; 
    
    // Trace of the stress tensor with respect to the Eulerian observer
    //------------------------------------
    
    s_euler = 3*press() + Ann()*uu2 + Anp()*uuu()*uuu2() + App()*uu22  ;

    s_euler.set_std_base() ; 

    // momentum vector (E+p)U in cartesian components
    //------------------------------------------------

    u_euler = Ann*u_euler + 0.5*Anp*(u_euler + u_euler2) + App*u_euler2 ;
    u_euler.set_std_base() ;
    
    // The derived quantities are obsolete
    // -----------------------------------
    
    del_deriv() ;                
    

}


