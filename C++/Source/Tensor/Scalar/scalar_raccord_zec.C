/*
 *   Copyright (c) 2003 Eric Gourgoulhon & Jerome Novak
 *
 *   Copyright (c) 2000-2001 Philippe Grandclement (for preceding Cmp version)
 *
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


char scalar_raccord_zec_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2003/09/25 08:58:10  e_gourgoulhon
 * First version.
 *
 *
 * $Header$
 *
 */

//standard
#include <stdlib.h>

// LORENE
#include "matrice.h"
#include "tensor.h"
#include "proto.h"

// Fait le raccord C1 dans la zec ...
// Suppose (pour le moment, le meme nbre de points sur les angles ...)
// et que la zone precedente est une coquille

void Scalar::raccord_c1_zec(int puis, int nbre, int lmax) {
    
    assert (nbre>0) ;
    assert (etat != ETATNONDEF) ;
    if (etat == ETATZERO)
	return ;

    // Le mapping doit etre affine :
    const Map_af* map  = dynamic_cast<const Map_af*>(mp) ;
    if (map == 0x0) {
	cout << "Le mapping doit etre affine" << endl ;
	abort() ;
    }
    
    int nz = map->get_mg()->get_nzone() ;
    int nr = map->get_mg()->get_nr (nz-1) ;
    int nt = map->get_mg()->get_nt (nz-1) ;
    int np = map->get_mg()->get_np (nz-1) ;
    
    double alpha = map->get_alpha()[nz-1] ;
    double r_cont = -1./2./alpha ;	//Rayon de debut de la zec.
  
    // On calcul les coefficients des puissances de 1./r
    Tbl coef (nbre+2*lmax, nr) ;
    coef.set_etat_qcq() ;
    
    int* deg = new int[3] ;
    deg[0] = 1 ; deg[1] = 1 ; deg[2] = nr ;
    double* auxi = new double[nr] ;
    
    for (int conte=0 ; conte<nbre+2*lmax ; conte++) {
	for (int i=0 ; i<nr ; i++)
	    auxi[i] = pow(-1-cos(M_PI*i/(nr-1)), (conte+puis)) ;
	    
	cfrcheb(deg, deg, auxi, deg, auxi) ;
	for (int i=0 ; i<nr ; i++)
	    coef.set(conte, i) = auxi[i]*pow (alpha, conte+puis) ;
	}
     
    delete[] deg ;
    // Maintenant on va calculer les valeurs de la ieme derivee :
    Tbl valeurs (nbre, nt, np+1) ;
    valeurs.set_etat_qcq() ;
    
    Scalar courant (*this) ;
    double* res_val = new double[1] ;
    
    for (int conte=0 ; conte<nbre ; conte++) {
	
	courant.va.coef() ;
	courant.va.ylm() ;
	courant.va.c_cf->t[nz-1]->annule_hard() ;
	
	int base_r = courant.va.base.get_base_r(nz-2) ;
	for (int k=0 ; k<np+1 ; k++)
	    for (int j=0 ; j<nt ; j++) 
		if (nullite_plm(j, nt, k, np, courant.va.base) == 1) {
	    
		    for (int i=0 ; i<nr ; i++)
			auxi[i] = (*courant.va.c_cf)(nz-2, k, j, i) ;

		    switch (base_r) {
			case R_CHEB :
			    som_r_cheb (auxi, nr, 1, 1, 1, res_val) ;
			    break ;
			default :
			    cout << "Cas non prevu dans raccord_zec" << endl ;
			    abort() ;
			    break ;
		    }
		    valeurs.set(conte, k, j) = res_val[0] ;
		}
	Scalar copie (courant) ;
	copie.dec2_dzpuis() ;
	courant = copie.dsdr() ;
    }
 
    delete [] auxi ;
    delete [] res_val ; 
   
    // On boucle sur les harmoniques : construction de la matrice 
    // et du second membre
    va.coef() ;
    va.ylm() ;
    va.c_cf->t[nz-1]->annule_hard() ;
    va.set_etat_cf_qcq() ;
    
    const Base_val& base = va.base ;
    int base_r, l_quant, m_quant ;
    for (int k=0 ; k<np+1 ; k++)
	for (int j=0 ; j<nt ; j++) 
	    if (nullite_plm(j, nt, k, np, va.base) == 1) {
	    
	    donne_lm (nz, nz-1, j, k, base, m_quant, l_quant, base_r) ;
	    
	    if (l_quant<=lmax) {
	    
		Matrice systeme (nbre, nbre) ;
		systeme.set_etat_qcq() ;
	    
		for (int col=0 ; col<nbre ; col++)
		    for (int lig=0 ; lig<nbre ; lig++) {
			
			int facteur = (lig%2==0) ? 1 : -1 ;
			for (int conte=0 ; conte<lig ; conte++)
			    facteur *= puis+col+conte+2*l_quant ;
			systeme.set(lig, col) = facteur/pow(r_cont, puis+col+lig+2*l_quant) ;
		    }
		
		systeme.set_band(nbre, nbre) ;
		systeme.set_lu() ;
		
	        Tbl sec_membre (nbre) ;
		sec_membre.set_etat_qcq() ;
		for (int conte=0 ; conte<nbre ; conte++)
		    sec_membre.set(conte) = valeurs(conte, k, j) ;
		
		Tbl inv (systeme.inverse(sec_membre)) ;
		
		for (int conte=0 ; conte<nbre ; conte++)
		    for (int i=0 ; i<nr ; i++)
			va.c_cf->set(nz-1, k, j, i)+= 
			    inv(conte)*coef(conte+2*l_quant, i) ;    
	    }
	else for (int i=0 ; i<nr ; i++)
		va.c_cf->set(nz-1, k, j, i)
		    = 0 ;
	}
	
    va.ylm_i() ;
    set_dzpuis (0) ;
}
