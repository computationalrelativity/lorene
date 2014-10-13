/*
 * Methods of class Excised_slice
 *
 * (see file excised_slice.h for documentation)
 * Disclaimer: the class Isol_hole() is redundant with this class under a set of parameters;
 * therefore it has to go at some point.
 */
 

/*
 *   Copyright (c) 2010 Nicolas Vasset

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
#include "excised_slice.h"
#include "spheroid.h"
#include "excision_surf.h"
#include "excision_hor.h"
#include "utilitaires.h"
#include "param.h"
#include "unites.h"
#include "proto.h"

namespace Lorene {

    // Fundamental constants and units
    // -------------------------------
    using namespace Unites ;

			    //--------------//
			    // Constructors //
			    //--------------//


// Standard constructor
// --------------------
Excised_slice::Excised_slice (const Map& mpi, int hor_type, int scheme_type)
  : mp(mpi),   
    type_hor(hor_type),
    field_set(scheme_type),
    lapse(mpi),
    conf_fact(mpi),
    shift(mpi, CON, mpi.get_bvect_spher()),
    hij(mpi, CON, mpi.get_bvect_spher()),
    hatA(mpi, CON, mpi.get_bvect_spher()),
    Xx(mpi, CON, mpi.get_bvect_spher()){

      // Pointers of derived quantities initialized to zero : 
    set_der_0x0() ;

  // Initializing primary quantities.
 
    lapse = 1. ;
    conf_fact = 1. ;
    shift.set_etat_zero() ;
    hij.set_etat_zero();
    hatA.set_etat_zero();
    Xx.set_etat_zero();
}
     
// Copy constructor
// ----------------
Excised_slice::Excised_slice(const Excised_slice& ih) 
		 : mp(ih.mp),
		   type_hor(ih.type_hor),
		   field_set(ih.field_set),
		   lapse(ih.lapse),
		   conf_fact(ih.conf_fact),
		   shift(ih.shift),
		   hij(ih.hij),
		   hatA(ih.hatA),
                   Xx(ih.Xx){
	       
    set_der_0x0() ;

}

// Constructor from a file
// -----------------------
Excised_slice::Excised_slice(const Map& mpi, int hor_type, int scheme_type, FILE* fich)
		 : mp(mpi),
		   type_hor(hor_type),
		   field_set(scheme_type),
                   lapse(mpi,*(mpi.get_mg()), fich), 
		   conf_fact(mpi,*(mpi.get_mg()), fich), 
		   shift(mpi, mpi.get_bvect_spher(), fich),
		   hij(mpi, mpi.get_bvect_spher(), fich),
		   hatA(mpi, mpi.get_bvect_spher(), fich), 
                   Xx(mpi, mpi.get_bvect_spher(), fich){



    // Pointers of derived quantities initialized to zero 
    // --------------------------------------------------
    set_der_0x0() ;
    
}

			    //------------//
			    // Destructor //
			    //------------//

Excised_slice::~Excised_slice(){

    del_deriv() ; 

}


			//----------------------------------//
			// Management of derived quantities //
			//----------------------------------//

void Excised_slice::del_deriv() const {
  if (p_hor != 0x0) delete p_hor ;
  if (p_adm_mass != 0x0) delete p_adm_mass ;
  if (p_komar_angmom != 0x0) delete p_komar_angmom ;
  if (p_virial_residue != 0x0) delete p_virial_residue ;
    Excised_slice::set_der_0x0() ; 
}			    


void Excised_slice::set_der_0x0() const {
  p_hor = 0x0;
    p_adm_mass = 0x0 ; 
    p_komar_angmom = 0x0 ;
    p_virial_residue = 0x0 ;
}			    



			    //--------------//
			    //  Assignment  //
			    //--------------//

// Assignment to another Excised_slice
// ----------------------------
void Excised_slice::operator=(const Excised_slice& ih) {

    assert( &(ih.mp) == &mp ) ;		    // Same mapping  
    type_hor = ih.type_hor ;
    field_set = ih.field_set ;
    lapse = ih.lapse ;
    conf_fact = ih.conf_fact ;
    shift = ih.shift ;
    hij = ih.hij ;
    hatA = ih.hatA ;
    Xx = ih.Xx;

    del_deriv() ;  // Deletes all derived quantities

}	

			    //--------------//
			    //	  Outputs   //
			    //--------------//

// Save in a file
// --------------
void Excised_slice::sauve(FILE* fich) const {

    lapse.sauve(fich) ;
    conf_fact.sauve(fich) ;
    shift.sauve(fich);
    hij.sauve(fich);
    hatA.sauve(fich);
    Xx.sauve(fich);}


// Prints out maximal errors in Einstein equations for the obtained metric fields 

void Excised_slice::Einstein_errors() {

const Map_af* map = dynamic_cast<const Map_af*>(&mp) ;
 const Metric_flat& mets = (*map).flat_met_spher() ; 
  
  Sym_tensor gamtcon = mets.con() + hij;
  Metric gamt(gamtcon);
  Sym_tensor gamcon = gamtcon/(conf_fact*conf_fact*conf_fact*conf_fact);
  gamcon.std_spectral_base();
  Metric gam(gamcon);
  Sym_tensor k_uu = hatA/(pow(conf_fact,10)) ;
  k_uu.std_spectral_base();
  k_uu.dec_dzpuis(k_uu(1,1).get_dzpuis()); //WTF?
  Sym_tensor k_dd = k_uu.up_down(gam);

  Scalar TrK3 = k_uu.trace(gam);
 
  // Hamiltonian constraint
  //-----------------------
  Scalar ham_constr = gam.ricci_scal() ;
  ham_constr.dec_dzpuis(3) ; // To check
  ham_constr +=  TrK3*TrK3 - contract(k_uu, 0, 1, k_dd, 0, 1) ;
  maxabs(ham_constr, "Hamiltonian constraint: ") ;
 
  // Momentum constraint
  //-------------------
  Vector mom_constr = k_uu.divergence(gam)  - TrK3.derive_con(gam) ;
  mom_constr.dec_dzpuis(2) ; // To check
  maxabs(mom_constr, "Momentum constraint: ") ;

  // Evolution equations
  //--------------------
  Sym_tensor evol_eq = lapse*gam.ricci() 
    - lapse.derive_cov(gam).derive_cov(gam);
  evol_eq.dec_dzpuis() ;
  evol_eq += k_dd.derive_lie(shift) ;
  evol_eq.dec_dzpuis(2) ; // To check
  evol_eq += lapse*(TrK3*k_dd - 2*contract(k_dd, 1, k_dd.up(0, gam), 0) ) ;
  maxabs(evol_eq, "Evolution equations: ") ;
    
  return; 
}





                  //----------------------------//
                  //  Accessors/ Derived data   //
                  //----------------------------//

// Computation of the Spheroid corresponding to the black hole excision surface.

Spheroid Excised_slice::hor() {
 const Map_af* map = dynamic_cast<const Map_af*>(&mp) ;
  const Mg3d* mgrid = (*map).get_mg();
	
  // Construct angular grid for h(theta,phi) 
  const Mg3d* g_angu = (*mgrid).get_angu_1dom() ;

  const Coord& rr = (*map).r;
   Scalar rrr (*map) ; 
  rrr = rr ; 
  rrr.std_spectral_base();  
  assert((rrr.val_grid_point(1,0,0,0) - 1.) <= 1.e-9); // For now the code handles only horizons at r=1, corresponding to the first shell inner boundary. This test assures this is the case with our mapping.
 
  // Angular mapping defined for one domain (argument of spheroid Class)
  //--------------------------------------------------------------------

  double r_limits2[] = {rrr.val_grid_point(1,0,0,0), rrr.val_grid_point(2,0,0,0)} ; 
  const Map_af map_2(*g_angu, r_limits2);

  //Full 3-metric and extrrinsic curvature 
 const Metric_flat& mets = (*map).flat_met_spher() ; 
  
  Sym_tensor gamtcon = mets.con() + hij;
  Metric gamt(gamtcon);
  Sym_tensor gamcon = gamtcon/(conf_fact*conf_fact*conf_fact*conf_fact);
  Metric gam(gamcon);


  Sym_tensor kuu = hatA/pow(conf_fact,10) ;
  kuu.std_spectral_base();
  Sym_tensor kdd = kuu.up_down(gam);
 
  //---------------------------------------------------------
  // Construction of the spheroid associated with those data 
  //--------------------------------------------------------
  double hor_posd = rrr.val_grid_point(1,0,0,0);
  Scalar hor_pos(map_2); hor_pos = hor_posd; hor_pos.std_spectral_base();
  Spheroid hor_loc(hor_pos, gam, kdd);
  return hor_loc; 
}

// Computation of the ADM mass of the BH spacetime
double Excised_slice::adm_mass() {
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
 
    double M_ADM = - (1/(2.*3.1415927*ggrav))*(*map).integrale_surface_infini(admintegr*detgam);
  return M_ADM;
}


// Computation of the Komar angular momentum w.r.t. assumed rotational symmetry
double Excised_slice:: komar_angmom() {
const Map_af* map = dynamic_cast<const Map_af*>(&mp) ;
 const Metric_flat& mets = (*map).flat_met_spher() ; 
  
  Sym_tensor gamtcon = mets.con() + hij;
  Sym_tensor gamcon = gamtcon/(conf_fact*conf_fact*conf_fact*conf_fact);
  gamcon.std_spectral_base();
  Metric gam(gamcon);
  Sym_tensor k_uu = hatA/(pow(conf_fact,10));
  k_uu.std_spectral_base();
  k_uu.dec_dzpuis(k_uu(1,1).get_dzpuis()); //WTF?
  Sym_tensor k_dd = k_uu.up_down(gam);
 
  Scalar detgam = sqrt((gam.cov())(2,2)*(gam.cov())(3,3) - (gam.cov())(2,3)*(gam.cov())(3,2));
  detgam.std_spectral_base();
    Scalar contraction = k_dd(1,3); contraction.mult_r_dzpuis(2); contraction.mult_sint();
    double angu_komar = - (1/(8.*3.1415927*ggrav))*(*map).integrale_surface_infini(detgam*contraction);

  return angu_komar;
}


// Computation of the Virial residual, as rescaled difference at infinity
// between the ADM mass and the Komar integral associated to the mass

double Excised_slice::virial_residue() {
  const Mg3d* mgrid = mp.get_mg();
  const int nz = (*mgrid).get_nzone(); 	// Number of domains
  Valeur** devel_psi (conf_fact.asymptot(1)) ;
    Valeur** devel_n (lapse.asymptot(1)) ;
    

    double erreur = (2*(*devel_psi[1])(nz-1, 0, 0, 0)
		     + (*devel_n[1])(nz-1, 0, 0, 0))/fabs ((*devel_n[1])(nz-1, 0, 0, 0)) ;

    return erreur;
}





}
