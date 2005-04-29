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
#include "graphique.h"

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
    
double Bin_hor::ang_mom_adm() const {
    
    Scalar integrand_un (hole1.aa_auto().up_down(hole1.met_gamt)(1,3) 
			 - hole1.gam_dd()(1,3) * hole1.trK) ;
    Scalar integrand_deux (hole2.aa_auto().up_down(hole2.met_gamt)(1,3)) ;
    
    integrand_un.mult_rsint() ;  // in order to pass from the triad components
    integrand_deux.mult_rsint() ; // to the coordinate basis
    
    double mom = hole1.mp.integrale_surface_infini(integrand_un)
	+ hole2.mp.integrale_surface_infini(integrand_deux) ;
    
    mom /= 8*M_PI ;
    return mom ;
    
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
