/*
 *   Copyright (c) 2005 Francois Limousin
 *                      Jose Luis Jaramillo
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


char binhor_glob_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.7  2006/05/24 16:56:37  f_limousin
 * Many small modifs.
 *
 * Revision 1.6  2005/09/24 08:31:38  f_limousin
 * Improve the computation of moment_adm() and moment_hor().
 *
 * Revision 1.5  2005/09/13 18:33:15  f_limousin
 * New function vv_bound_cart_bin(double) for computing binaries with
 * berlin condition for the shift vector.
 * Suppress all the symy and asymy in the importations.
 *
 * Revision 1.4  2005/06/09 16:12:04  f_limousin
 * Implementation of amg_mom_adm().
 *
 * Revision 1.3  2005/04/29 14:02:44  f_limousin
 * Important changes : manage the dependances between quantities (for
 * instance psi and psi4). New function write_global(ost).
 *
 * Revision 1.2  2005/03/04 17:09:57  jl_jaramillo
 * Change to avoid warnings
 *
 * Revision 1.1  2005/03/03 13:48:56  f_limousin
 * First version
 *
 *
 * $Header$
 *
 */



//standard
#include <stdlib.h>
#include <math.h>

// Lorene
#include "nbr_spx.h"
#include "tensor.h"
#include "isol_hor.h"
#include "proto.h"
#include "utilitaires.h"
//#include "graphique.h"

double Bin_hor::adm_mass() const {
 
   Vector dpsi_un (hole1.psi_auto().derive_con(hole1.ff)) ;
    Vector dpsi_deux (hole2.psi_auto().derive_con(hole2.ff)) ;

    Vector ww (0.125*(hole1.hdirac() - (hole1.hh().trace(hole1.ff)).
		      derive_con(hole1.ff))) ;

    double inf = hole1.mp.val_r(hole1.mp.get_mg()->get_nzone()-1, 1., 0., 0.) ;
    
    double masse = dpsi_un.flux(inf, hole1.ff) + 
	           dpsi_deux.flux(inf, hole2.ff) +
	           ww.flux(inf, hole1.ff) ;
    masse /= -2*M_PI ;
    return masse ;
}

double Bin_hor::komar_mass() const {

    Vector dnn_un (hole1.n_auto().derive_con(hole1.ff)) ;
    Vector dnn_deux (hole2.n_auto().derive_con(hole2.ff)) ;
    
    Vector ww (contract(hole1.hh(), 1, hole1.nn().derive_cov(hole1.ff), 0)) ;
	       
    double inf = hole1.mp.val_r(hole1.mp.get_mg()->get_nzone()-1, 1., 0., 0.) ;

    double mass = dnn_un.flux(inf, hole1.ff) + 
	dnn_deux.flux(inf, hole2.ff) + 
	ww.flux(inf, hole1.ff) ;
    
    mass /= 4*M_PI ;
    return mass ;
}
    
double Bin_hor::ang_mom_hor() const {
    
    if (omega == 0)
	return 0 ;
    else {
	// Alignement
	double orientation_un = hole1.mp.get_rot_phi() ;
	assert ((orientation_un==0) || (orientation_un==M_PI)) ;
	double orientation_deux = hole2.mp.get_rot_phi() ;
	assert ((orientation_deux==0) || (orientation_deux==M_PI)) ;
	int same_orient = (orientation_un == orientation_deux) ? 1 : -1 ;
	
	// Integral on the first horizon :
	Scalar xa_un (hole1.mp) ;
	xa_un = hole1.mp.xa ;
	xa_un.std_spectral_base() ;
	
	Scalar ya_un (hole1.mp) ;
	ya_un = hole1.mp.ya ;
	ya_un.std_spectral_base() ;
	
	Sym_tensor tkij_un (hole1.aa()) ;
	tkij_un.change_triad(hole1.mp.get_bvect_cart()) ;

	Vector vecteur_un (hole1.mp, CON, hole1.mp.get_bvect_cart()) ;
	for (int i=1 ; i<=3 ; i++)
	    vecteur_un.set(i) = -ya_un*tkij_un(1, i)+
		xa_un * tkij_un(2, i) ;
	vecteur_un.annule_domain(hole1.mp.get_mg()->get_nzone()-1) ;
	vecteur_un.change_triad (hole1.mp.get_bvect_spher()) ;
	
	Scalar integrant_un (hole1.psi4()*hole1.psi()*hole1.psi()
			     *vecteur_un(1)) ;
	double moment_un = hole1.mp.integrale_surface
	    (integrant_un, hole1.radius+1e-12)/8/M_PI ;
	
	//Integral on the second horizon :
	Scalar xa_deux (hole2.mp) ;
	xa_deux = hole2.mp.xa ;
	xa_deux.std_spectral_base() ;
	
	Scalar ya_deux (hole2.mp) ;
	ya_deux = hole2.mp.ya ;
	ya_deux.std_spectral_base() ;
	
	Sym_tensor tkij_deux (hole2.aa()) ;
	tkij_deux.change_triad(hole2.mp.get_bvect_cart()) ;

	Vector vecteur_deux (hole2.mp, CON, hole2.mp.get_bvect_cart()) ;
	for (int i=1 ; i<=3 ; i++)
	    vecteur_deux.set(i) = -ya_deux*tkij_deux(1, i)+
		xa_deux * tkij_deux(2, i) ;
	vecteur_deux.annule_domain(hole2.mp.get_mg()->get_nzone()-1) ;
	vecteur_deux.change_triad (hole2.mp.get_bvect_spher()) ;
	
	Scalar integrant_deux (hole2.psi4()*hole2.psi()*hole2.psi()
			       *vecteur_deux(1)) ;
	double moment_deux = hole2.mp.integrale_surface
	    (integrant_deux, hole2.radius+1e-12)/8/M_PI ;
	
	return moment_un+same_orient*moment_deux ;
	}
}



double Bin_hor::ang_mom_adm() const {
/*    
    Scalar integrand_un (hole1.k_dd()(1,3) - hole1.gam_dd()(1,3) * hole1.trK) ;
    
    integrand_un.mult_rsint() ;  // in order to pass from the triad components
    
    double mom = hole1.mp.integrale_surface_infini(integrand_un) ;    
    mom /= 8*M_PI ;
    return mom ;
*/

    if (omega == 0)
	return 0 ;
    else {
	// Alignement 
	double orientation_un = hole1.mp.get_rot_phi() ;
	assert ((orientation_un==0) || (orientation_un==M_PI)) ;
	double orientation_deux = hole2.mp.get_rot_phi() ;
	assert ((orientation_deux==0) || (orientation_deux==M_PI)) ;
	int same_orient = (orientation_un == orientation_deux) ? 1 : -1 ;
	
	// Construction of an auxiliar grid and mapping
	int nzones = hole1.mp.get_mg()->get_nzone() ;
	double* bornes = new double [nzones+1] ;
	double courant = (hole1.mp.get_ori_x()-hole2.mp.get_ori_x())+1 ;
	for (int i=nzones-1 ; i>0 ; i--) {
	bornes[i] = courant ;
	courant /= 2. ;
	}
	bornes[0] = 0 ;
	bornes[nzones] = __infinity ;
	
	Map_af mapping (*hole1.mp.get_mg(), bornes) ;
	
	delete [] bornes ; 
	
	// Construction of k_total
	Sym_tensor k_total (mapping, CON, mapping.get_bvect_cart()) ;
	
	Vector shift_un (mapping, CON, mapping.get_bvect_cart()) ;
	Vector shift_deux (mapping, CON, mapping.get_bvect_cart()) ;
	
	Vector beta_un (hole1.beta_auto()) ;
	Vector beta_deux (hole2.beta_auto()) ;
	beta_un.change_triad(hole1.mp.get_bvect_cart()) ;
	beta_deux.change_triad(hole2.mp.get_bvect_cart()) ;
	
	shift_un.set(1).import(beta_un(1)) ;
	shift_un.set(2).import(beta_un(2)) ;
	shift_un.set(3).import(beta_un(3)) ;
	
	shift_deux.set(1).import(same_orient*beta_deux(1)) ;
	shift_deux.set(2).import(same_orient*beta_deux(2)) ;
	shift_deux.set(3).import(beta_deux(3)) ;
	
	Vector shift_tot (shift_un+shift_deux) ;
	shift_tot.std_spectral_base() ;
	shift_tot.annule(0, nzones-2) ;

	// Substract the residuals
	shift_tot.inc_dzpuis(2) ;
	shift_tot.dec_dzpuis(2) ;
	
	const Metric_flat& flat0 (mapping.flat_met_cart()) ;

	k_total = shift_tot.ope_killing_conf(flat0) / 2. ;
	//- flat0.con() * hole1.trK;

	for (int lig=1 ; lig<=3 ; lig++)
	    for (int col=lig ; col<=3 ; col++)
		k_total.set(lig, col).mult_r_ced() ;
	
	Vector vecteur_un (mapping, CON, mapping.get_bvect_cart()) ;
	for (int i=1 ; i<=3 ; i++)
	    vecteur_un.set(i) = k_total(1, i) ;
	vecteur_un.change_triad (mapping.get_bvect_spher()) ;
	Scalar integrant_un (vecteur_un(1)) ;
	
	Vector vecteur_deux (mapping, CON, mapping.get_bvect_cart()) ;
	for (int i=1 ; i<=3 ; i++)
	    vecteur_deux.set(i) = k_total(2, i) ;
	vecteur_deux.change_triad (mapping.get_bvect_spher()) ;
	Scalar integrant_deux (vecteur_deux(1)) ;
	
	// Multiplication by y and x :
	integrant_un.set_spectral_va() = integrant_un.get_spectral_va()
	    .mult_st() ;
	integrant_un.set_spectral_va() = integrant_un.get_spectral_va()
	    .mult_sp() ;
	
	integrant_deux.set_spectral_va() = integrant_deux.get_spectral_va()
	    .mult_st() ;
	integrant_deux.set_spectral_va() = integrant_deux.get_spectral_va()
	    .mult_cp() ;
	
	double moment = mapping.integrale_surface_infini (-integrant_un
							  +integrant_deux) ;
	
	moment /= 8*M_PI ;
	
	return moment ;
    }
}


/*
double Bin_hor::proper_distance(const int nr) const {
    

    // On determine les rayons coordonnes des points limites de l'integrale :
    double x_un = hole1.mp.get_ori_x() - hole1.rayon ;
    double x_deux = hole2.mp.get_ori_x() + hole2.rayon ;
    
    // Les coefficients du changement de variable :
    double pente = 2./(x_un-x_deux) ;
    double constante = - (x_un+x_deux)/(x_un-x_deux) ;
    
    
    double ksi ; // variable d'integration.
    double xabs ; // x reel.
    double air_un, theta_un, phi_un ; // coordonnee spheriques 1
    double air_deux, theta_deux, phi_deux ; // coordonnee spheriques 2
    
    double* coloc = new double[nr] ;
    double* coef = new double[nr] ;
    int* deg = new int[3] ;
    deg[0] = 1 ; deg[1] = 1 ; deg[2] = nr ;
    
    for (int i=0 ; i<nr ; i++) {
	ksi = -cos (M_PI*i/(nr-1)) ;
	xabs = (ksi-constante)/pente ;
	
	hole1.mp.convert_absolute (xabs, 0, 0, air_un, theta_un, phi_un) ;
	hole2.mp.convert_absolute (xabs, 0, 0, air_deux, theta_deux, phi_deux) ;
	
	coloc[i] = pow (hole1.psi_auto().val_point (air_un, theta_un, phi_un) +
		   hole2.psi_auto().val_point (air_deux, theta_deux, phi_deux), 2.) ;
    }
    
    // On prend les coefficients de la fonction
    cfrcheb(deg, deg, coloc, deg, coef) ;
    
    // On integre
    double* som = new double[nr] ;
    som[0] = 2 ;
    for (int i=2 ; i<nr ; i+=2)
	som[i] = 1./(i+1)-1./(i-1) ;
    for (int i=1 ; i<nr ; i+=2)
	som[i] = 0 ;
    
    double res = 0 ;
    for (int i=0 ; i<nr ; i++)
	res += som[i]*coef[i] ;
	
    res /= pente ;
    
    delete [] deg ;
    delete [] coef ;
    delete [] coloc ;
    delete [] som ;
   
    return res ;

}
*/
