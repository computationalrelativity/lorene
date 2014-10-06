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
 * Revision 1.7  2014/10/06 15:16:16  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.6  2004/06/04 16:14:18  j_novak
 * In smooth_decay, the configuration space was not up-to-date.
 *
 * Revision 1.5  2004/03/05 15:09:59  e_gourgoulhon
 * Added method smooth_decay.
 *
 * Revision 1.4  2003/10/29 11:04:34  e_gourgoulhon
 * dec2_dpzuis() replaced by dec_dzpuis(2).
 * inc2_dpzuis() replaced by inc_dzpuis(2).
 *
 * Revision 1.3  2003/10/10 15:57:29  j_novak
 * Added the state one (ETATUN) to the class Scalar
 *
 * Revision 1.2  2003/09/25 09:22:33  j_novak
 * Added a #include<math.h>
 *
 * Revision 1.1  2003/09/25 08:58:10  e_gourgoulhon
 * First version.
 *
 *
 * $Header$
 *
 */

//standard
#include <cstdlib>
#include <cmath>

// LORENE
#include "matrice.h"
#include "tensor.h"
#include "proto.h"
#include "nbr_spx.h"
#include "utilitaires.h"

// Fait le raccord C1 dans la zec ...
// Suppose (pour le moment, le meme nbre de points sur les angles ...)
// et que la zone precedente est une coquille

void Scalar::raccord_c1_zec(int puis, int nbre, int lmax) {
    
    assert (nbre>0) ;
    assert (etat != ETATNONDEF) ;
    if ((etat == ETATZERO) || (etat == ETATUN))
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
	copie.dec_dzpuis(2) ;
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


//***************************************************************************

void Scalar::smooth_decay(int kk, int nn) {

    assert(kk >= 0) ; 

    if (etat == ETATZERO) return ; 
    if (va.get_etat() == ETATZERO) return ; 
    
    const Mg3d& mgrid = *(mp->get_mg()) ; 
    
    int nz = mgrid.get_nzone() ; 
    int nzm1 = nz - 1 ; 
    int np = mgrid.get_np(nzm1) ; 
    int nt = mgrid.get_nt(nzm1) ; 
    int nr_ced = mgrid.get_nr(nzm1) ; 
    int nr_shell = mgrid.get_nr(nzm1-1) ; 
    
    assert(mgrid.get_np(nzm1-1) == np) ; 
    assert(mgrid.get_nt(nzm1-1) == nt) ; 
    
    assert(mgrid.get_type_r(nzm1) == UNSURR) ; 

    // In the present version, the mapping must be affine :
    const Map_af* mapaf  = dynamic_cast<const Map_af*>(mp) ;
    if (mapaf == 0x0) {
	    cout << "Scalar::smooth_decay: present version supports only \n" 
         << "  affine mappings !" << endl ;
	    abort() ;
    }


    assert(va.get_etat() == ETATQCQ) ; 
    
    
    // Spherical harmonics are required
    va.ylm() ; 
    
    // Array for the spectral coefficients in the CED:
    assert( va.c_cf->t[nzm1] != 0x0) ; 
    Tbl& coefresu = *( va.c_cf->t[nzm1] ) ; 
    coefresu.set_etat_qcq() ; 
    
    // 1-D grid 
    //----------
    int nbr1[] = {nr_shell, nr_ced} ; 
    int nbt1[] = {1, 1} ; 
    int nbp1[] = {1, 1} ; 
    int typr1[] = {mgrid.get_type_r(nzm1-1), UNSURR} ; 
    Mg3d grid1d(2, nbr1, typr1, nbt1, SYM, nbp1, SYM) ; 
    
    // 1-D mapping
    //------------
    double rbound = mapaf->get_alpha()[nzm1-1] 
                    + mapaf->get_beta()[nzm1-1] ; 
    double rmin = - mapaf->get_alpha()[nzm1-1] 
                    + mapaf->get_beta()[nzm1-1] ; 
    double bound1[] = {rmin, rbound, __infinity} ; 

    Map_af map1d(grid1d, bound1) ;     
    
    // 1-D scalars
    // -----------
    Scalar uu(map1d) ; 
    uu.std_spectral_base() ; 
    Scalar duu(map1d) ; 
    Scalar vv(map1d) ; 
    Scalar tmp(map1d) ; 

    const Base_val& base = va.get_base() ; 
    
    // Loop on the spherical harmonics
    // -------------------------------
    for (int k=0; k<=np; k++) {
        for (int j=0; j<nt; j++) {
            
	        if (nullite_plm(j, nt, k, np, base) != 1) {
                for (int i=0 ; i<nr_ced ; i++) {
                    coefresu.set(k, j, i) = 0 ;
                }
            }
            else {
                int base_r_ced, base_r_shell , l_quant, m_quant ;
	            donne_lm(nz, nzm1, j, k, base, m_quant, l_quant, base_r_ced) ;
	            donne_lm(nz, nzm1-1, j, k, base, m_quant, l_quant, base_r_shell) ;
               
                int nn0 = l_quant + nn ;    // lowerst power of decay
                
                uu.set_etat_qcq() ; 
                uu.va.set_etat_cf_qcq() ; 
                uu.va.c_cf->set_etat_qcq() ; 
                uu.va.c_cf->t[0]->set_etat_qcq() ;
                uu.va.c_cf->t[1]->set_etat_qcq() ;
                for (int i=0; i<nr_shell; i++) {
                    uu.va.c_cf->t[0]->set(0, 0, i) = 
                        va.c_cf->operator()(nzm1-1, k, j, i) ; 
                    uu.va.c_cf->t[0]->set(1, 0, i) = 0. ; 
                    uu.va.c_cf->t[0]->set(2, 0, i) = 0. ;
                }
 
                uu.va.c_cf->t[1]->set_etat_zero() ; 

                uu.va.set_base_r(0, base_r_shell) ; 
                uu.va.set_base_r(1, base_r_ced) ; 
 
                // Computation of the derivatives up to order k at the outer
                // of the last shell:
                // ----------------------------------------------------------
                Tbl derivb(kk+1) ;
                derivb.set_etat_qcq() ; 
                duu = uu ; 
                derivb.set(0) = uu.val_grid_point(0, 0, 0, nr_shell-1) ; 
                for (int p=1; p<=kk; p++) {
                    tmp = duu.dsdr() ; 
                    duu = tmp ; 
                    derivb.set(p) = duu.val_grid_point(0, 0, 0, nr_shell-1) ; 
                }
               
                // Matrix of the linear system to be solved
                // ----------------------------------------
                
                Matrice mat(kk+1,kk+1) ; 
                mat.set_etat_qcq() ; 
                
                for (int im=0; im<=kk; im++) {
                    
                    double fact = (im%2 == 0) ? 1. : -1. ; 
                    fact /= pow(rbound, nn0 + im) ; 
                
                    for (int jm=0; jm<=kk; jm++) {
                        
                        double prod = 1 ; 
                        for (int ir=0; ir<im; ir++) {
                            prod *= nn0 + jm + ir ; 
                        } 
                        
                        mat.set(im, jm) = fact * prod / pow(rbound, jm) ; 
                        
                    }
                }
           
                // Resolution of the linear system to get the coefficients
                // of the 1/r^i functions
                //---------------------------------------------------------
                mat.set_band(kk+1, kk+1) ;
		        mat.set_lu() ;
                Tbl aa = mat.inverse( derivb ) ; 

                // Decaying function 
                // -----------------
                vv = 0 ; 
                const Coord& r = map1d.r ; 
                for (int p=0; p<=kk; p++) {
                    tmp = aa(p) / pow(r, nn0 + p) ; 
                    vv += tmp ; 
                }
                vv.va.set_base( uu.va.get_base() ) ; 

                // The coefficients of the decaying function 
                // are set to the result
                //-------------------------------------------
                vv.va.coef() ; 
                
                if (vv.get_etat() == ETATZERO) {
                     for (int i=0; i<nr_ced; i++) {
                        coefresu.set(k,j,i) = 0 ; 
                    }
                }
                else {
                    assert( vv.va.c_cf != 0x0 ) ; 
                    for (int i=0; i<nr_ced; i++) {
                        coefresu.set(k,j,i) = vv.va.c_cf->operator()(1,0,0,i) ; 
                    }
                }
                
               
            }
         
            
        }
    } 
    
    if (va.c != 0x0) {
      delete va.c ;
      va.c = 0x0 ;
    }
    va.ylm_i() ;     
    
    // Test of the computation
    // -----------------------
        
    Scalar tmp1(*this) ; 
    Scalar tmp2(*mp) ; 
    for (int p=0; p<=kk; p++) {
        double rd = pow(rbound, tmp1.get_dzpuis()) ; 
        double err = 0 ; 
        for (int k=0; k<np; k++) {
            for (int j=0; j<nt; j++) {
                double diff = fabs( tmp1.val_grid_point(nzm1, k, j, 0) / rd - 
                  tmp1.val_grid_point(nzm1-1, k, j, nr_shell-1) ) ;
                if (diff > err) err = diff ;  
            }
        }
        cout << 
        "Scalar::smooth_decay: Max error matching of " << p 
            << "-th derivative: " << err << endl ; 
        tmp2 = tmp1.dsdr() ; 
        tmp1 = tmp2 ; 
    }
    

}
