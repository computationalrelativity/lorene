/*
 * Methods of class Isol_hole
 *
 * (see file isol_hole.h for documentation)
 */

/*
 *   Copyright (c) 2009 Nicolas Vasset

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


// Headers C
#include "math.h"

// Headers Lorene
#include "isol_hole.h"
#include "spheroid.h"
#include "excision_surf.h"
#include "excision_hor.h"
#include "utilitaires.h"
#include "param.h"
#include "unites.h"
#include "proto.h"
#include "map.h"
#include "scalar.h"

    // Fundamental constants and units
    // -------------------------------
    using namespace Unites ;

			    //--------------//
			    // Constructors //
			    //--------------//


// Standard constructor
// --------------------
Isol_hole::Isol_hole (Map& mpi, double Omega_i, bool NorKappa_i, Scalar NoK_i, bool isCF_i)
  : mp(mpi), 
                   Omega(Omega_i),
		   NorKappa(NorKappa_i),
		   boundNoK(NoK_i),
		   isCF (isCF_i),
		   lapse(mpi),
		   conf_fact(mpi),
		   shift(mpi, CON, mpi.get_bvect_spher()),
		   hij(mpi, CON, mpi.get_bvect_spher()),
		   hatA(mpi, CON, mpi.get_bvect_spher()),
		   hor(*(dynamic_cast<const Map_af*>(&mpi)),1.){
		   




  assert (boundNoK.get_mp() == mpi); // Check if this is not too strong a condition
    
      // Pointers of derived quantities initialized to zero : 
    set_der_0x0() ;

  // Initializing primary quantities.
 
    lapse = 1. ;
    conf_fact = 1. ;
    shift.set_etat_zero() ;
    hij.set_etat_zero();
    hatA.set_etat_zero();
}
     
// Copy constructor
// ----------------
Isol_hole::Isol_hole(const Isol_hole& ih) 
		 : mp(ih.mp), 
		   Omega(ih.Omega),
		   NorKappa(ih.NorKappa),
		   boundNoK(ih.boundNoK),
		   isCF (ih.isCF),
		   lapse(ih.lapse),
		   conf_fact(ih.conf_fact),
		   shift(ih.shift),
		   hij(ih.hij),
		   hatA(ih.hatA),
		   hor(ih.hor){
	       
    set_der_0x0() ;

}

// Constructor from a file
// -----------------------
Isol_hole::Isol_hole(Map& mpi, double Omega_i, bool NorKappa_i, Scalar NoK_i, bool isCF_i, FILE* fich)
		 : mp(mpi),
		   Omega(Omega_i),
		   NorKappa(NorKappa_i),
		   boundNoK(NoK_i),
		   isCF(isCF_i),
                   lapse(mpi,*(mpi.get_mg()), fich), 
		   conf_fact(mpi,*(mpi.get_mg()), fich), 
		   shift(mpi, mpi.get_bvect_spher(), fich),
		   hij(mpi, mpi.get_bvect_spher(), fich),
		   hatA(mpi, mpi.get_bvect_spher(), fich),
		   hor(*(dynamic_cast<const Map_af*>(&mpi)),1.){

  cout << "WARNING: constructor from a FILE for Isol_hole is not fully operational yet, and does not assign the correct value to the associated spheroid" << endl;


    // Pointers of derived quantities initialized to zero 
    // --------------------------------------------------
    set_der_0x0() ;
    
}

			    //------------//
			    // Destructor //
			    //------------//

Isol_hole::~Isol_hole(){

    del_deriv() ; 

}


			//----------------------------------//
			// Management of derived quantities //
			//----------------------------------//

void Isol_hole::del_deriv() const {

  if (p_adm_mass != 0x0) delete p_adm_mass ;
  if (p_komar_angmom != 0x0) delete p_komar_angmom ;
  if (p_virial_residue != 0x0) delete p_virial_residue ;
  
    Isol_hole::set_der_0x0() ; 
}			    




void Isol_hole::set_der_0x0() const {

    p_adm_mass = 0x0 ; 
    p_komar_angmom = 0x0 ;
    p_virial_residue = 0x0 ;
  
}			    



			    //--------------//
			    //  Assignment  //
			    //--------------//

// Assignment to another Isol_hole
// ----------------------------
void Isol_hole::operator=(const Isol_hole& ih) {

    assert( &(ih.mp) == &mp ) ;		    // Same mapping  
    Omega = ih.Omega;
    NorKappa = ih.NorKappa ;
    boundNoK = ih.boundNoK ;
    isCF = ih.isCF ;
    lapse = ih.lapse ;
    conf_fact = ih.conf_fact ;
    shift = ih.shift ;
    hij = ih.hij ;
    hatA = ih.hatA ;
    hor = ih.hor ;

    del_deriv() ;  // Deletes all derived quantities

}	

			    //--------------//
			    //	  Outputs   //
			    //--------------//

// Save in a file
// --------------
void Isol_hole::sauve(FILE* fich) const {

    lapse.sauve(fich) ;
    conf_fact.sauve(fich) ;
    shift.sauve(fich);
    hij.sauve(fich);
    hatA.sauve(fich);
}


                  //----------------------------//
                  //  Accessors/ Derived data   //
                  //----------------------------//

// Computation of the ADM mass of the BH spacetime
double Isol_hole::adm_mass() {
 const Map_af* map = dynamic_cast<const Map_af*>(&mp) ;
 const Metric_flat& mets = (*map).flat_met_spher() ; 
  
  Sym_tensor gamtcon = mets.con() + hij;
  Metric gamt(gamtcon);
  Sym_tensor gamcon = gamtcon/(conf_fact*conf_fact*conf_fact*conf_fact);
  Metric gam(gamcon);

  Scalar detgam = sqrt((gam.cov())(2,2)*(gam.cov())(3,3) - (gam.cov())(2,3)*(gam.cov())(3,2));
    detgam.std_spectral_base();  
    Vector corr_adm =  - (0.125*contract(gamt.cov().derive_con(mets),1,2));
    Scalar admintegr = conf_fact.dsdr() + corr_adm(1);
 
    double M_ADM = - (1/(2.*3.1415927))*(*map).integrale_surface_infini(admintegr*detgam);
  return M_ADM;
}


// Computation of the Komar angular momentum w.r.t. assumed rotational symmetry
double Isol_hole:: komar_angmom() {
const Map_af* map = dynamic_cast<const Map_af*>(&mp) ;
 const Metric_flat& mets = (*map).flat_met_spher() ; 
  
  Sym_tensor gamtcon = mets.con() + hij;
  Metric gamt(gamtcon);
  Sym_tensor gamcon = gamtcon/(conf_fact*conf_fact*conf_fact*conf_fact);
  Metric gam(gamcon);
  Sym_tensor k_uu = hatA/(conf_fact*conf_fact*conf_fact*conf_fact) ;
  k_uu.dec_dzpuis(k_uu(1,1).get_dzpuis()); 
  Sym_tensor k_dd = k_uu.up_down(gam);
 
  Scalar detgam = sqrt((gam.cov())(2,2)*(gam.cov())(3,3) - (gam.cov())(2,3)*(gam.cov())(3,2));
    Scalar contraction = k_dd(1,3); contraction.mult_r_dzpuis(2); contraction.mult_sint();
    double angu_komar = - (1/(8.*3.1415927*ggrav))*(*map).integrale_surface_infini(detgam*contraction);

  return angu_komar;
}


// Computation of the Virial residual, as rescaled difference at infinity
// between the ADM mass and the Komar integral associated to the mass

double Isol_hole::virial_residue() {
  const Mg3d* mgrid = mp.get_mg();
  const int nz = (*mgrid).get_nzone(); 	// Number of domains
  Valeur** devel_psi (conf_fact.asymptot(1)) ;
    Valeur** devel_n (conf_fact.asymptot(1)) ;
    
    double erreur = (2*(*devel_psi[1])(nz-1, 0, 0, 0)
		     + (*devel_n[1])(nz-1, 0, 0, 0))/fabs ((*devel_n[1])(nz-1, 0, 0, 0)) ;

    return erreur;
}



