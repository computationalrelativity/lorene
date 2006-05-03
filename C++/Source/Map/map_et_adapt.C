/*
 *  Method of the class Map_et to compute the mapping parameters to let the
 *  domain boundaries coincide with isosurfaces of a given scalar field.
 *
 *  (see file map.h for documentation).
 */

/*
 *   Copyright (c) 2000-2001 Eric Gourgoulhon
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


char map_et_adapt_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.4  2006/05/03 07:07:52  p_grandclement
 * Petite correction
 *
 * Revision 1.3  2006/04/25 07:21:59  p_grandclement
 * Various changes for the NS_BH project
 *
 * Revision 1.2  2003/10/03 15:58:48  j_novak
 * Cleaning of some headers
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.6  2000/11/10  15:19:30  eric
 * Appel de Valeur::equipot_outward plutot que Valeur::equipot
 *   dans la determination des isopotentielles.
 *
 * Revision 2.5  2000/10/23  14:01:17  eric
 * Changement des arguments (en autre, ajout de nz_search, qui est
 * desormais a priori independant de nzadapt).
 *
 * Revision 2.4  2000/10/20  13:13:01  eric
 * nzet --> nzadapt.
 * nz_search = nzadapt --> nz_search = nzadapt + 1
 *
 * Revision 2.3  2000/02/17  19:00:53  eric
 * nz_search = nzet  et non plus  nz_search = nz-1.
 *
 * Revision 2.2  2000/02/15  15:09:12  eric
 * Changement du Param : fact_echelle est desormais passe en double_mod.
 *
 * Revision 2.1  2000/01/03  15:46:59  eric
 * *** empty log message ***
 *
 * Revision 2.0  2000/01/03  13:30:59  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers C
#include <assert.h>

// Headers Lorene
#include "cmp.h"
#include "itbl.h"
#include "param.h"
#include "proto.h"

void Map_et::adapt(const Cmp& ent, const Param& par, int nbr_filtre) {
    
    // Parameters of the computation
    // -----------------------------
    
    int nitermax = par.get_int(0) ;
    int& niter = par.get_int_mod(0) ; 
    int nzadapt = par.get_int(1) ;
    int nz_search = par.get_int(2)  ;
    int adapt_flag = par.get_int(3) ;
    int j_bord = par.get_int(4) ;	// index of theta_*
    int k_bord = par.get_int(5) ;	// index of phi_*
    double precis = par.get_double(0) ;
    double fact_lamu = par.get_double(1) ;
    double fact_echelle = par.get_double(2) ;

    const Tbl& ent_limit = par.get_tbl(0) ; 
    
    int nz = mg->get_nzone() ;	
    int nt = mg->get_nt(0) ; 
    int np = mg->get_np(0) ; 
    
    
    // Protections 
    // -----------
    assert(ent.get_mp() == this) ; 
    assert(nzadapt < nz) ;
    assert(nzadapt > 0) ;
    assert(nz_search >= nzadapt) ; 
    for (int l=1; l<nz; l++) {
	assert( mg->get_nt(l) == nt ) ; 
	assert( mg->get_np(l) == np ) ; 
    }
    assert(ent_limit.get_ndim() == 1) ; 
    assert(ent_limit.get_taille() >= nzadapt) ; 
    assert(ent_limit.get_etat() == ETATQCQ) ; 

    // Initializations 
    // ---------------
    
    niter = 0 ;
    const double xi_max = 1 ; 
    
    //=======================================================================
    //   Special case of a fixed mapping : only a rescaling is performed
    //=======================================================================

    if ( (adapt_flag == 0) || (nzadapt == 0) ) { 

	homothetie(fact_echelle) ; 
	
	return ; 
    }        
    
    //=======================================================================
    //   General case 
    //=======================================================================
    
    // New mapping parameters
    // ----------------------
    
    double* nalpha = new double[nz] ;
    double* nbeta = new double[nz] ;

    Itbl l_ext(np, nt) ;
    Tbl  x_ext(np, nt) ;
    
    Tbl xtgg(np, nt, 1) ;	// working array constructed from the isosurface
				//  equation
    Tbl xtff(np, nt, 1) ;	// idem 

    //----------------------------------------------------------------
    //   Adaptation of the nucleus
    //----------------------------------------------------------------
    
    double ent0 = ent_limit(0) ; 

    // Search for the equipotential surface ent = ent0 
    //    --->  l = l_ext(theta, phi) 
    //		xi = x_ext(theta, phi)
    // -----------------------------------------------

    int niter0 ; 
    (ent.va).equipot_outward(ent0, nz_search, precis, nitermax, niter0, 
			     l_ext, x_ext) ; 

    niter = ( niter0 > niter ) ? niter0 : niter ;

    // The new mapping must be such that the found isosurface corresponds
    //  to xi = 1
    
    // Computation of r_ext(theta, phi) - r_eq ---> xtgg
    // -------------------------------------------------
    
    xtgg.set_etat_qcq() ; 
    assert(l_ext.get_etat() == ETATQCQ) ; 

    double r_eq = val_r_jk(0, xi_max, j_bord, k_bord) ;

    double* pxtgg = xtgg.t ; 
    int* pl_ext = l_ext.t ;
    double* px_ext = x_ext.t ;
    
    for (int k=0; k<np; k++) {
	for (int j=0; j<nt; j++) {

	    *pxtgg = val_r_jk(*pl_ext, *px_ext, j, k) - r_eq  ;
	    			    
	    // ... next one
	    pl_ext++ ;
	    px_ext++ ;
	    pxtgg++ ;
	}	   
    }	

    // Decomposition into even and odd harmonics in phi of xtgg = r_ext - r_eq :
    //	    even phi harmonics --> xtgg
    //	    odd  phi harmonics --> xtff
    // ------------------------------------------------------------------------
    
    int base_p = ( ((ent.va).base).b[0] ) & MSQ_P ;

    switch( base_p ) {

	case P_COSSIN_P : {	// case of only even harmonics in phi

	    xtff.set_etat_zero() ; 
	    break ; 
	}
	    
	case P_COSSIN :	{	// general case: a Fourier transform must be
				//  done

	    xtff.set_etat_qcq() ; 
	    double* pxtff = xtff.t ; 
	    pxtgg = xtgg.t ; 

	    assert(np>=4) ;
	    int deg[3] ;
	    int dimc[3] ;
	    deg[0] = np ;
	    deg[1] = nt ;
	    deg[2] = 1; 
	    dimc[0] = np + 2 ;
	    dimc[1] = nt ;
	    dimc[2] = 1 ; 
	    double* cf = new double[(np+2)*nt] ; 
	    double* cf0 = new double[(np+2)*nt] ; 
	    double* ff0 = new double[np*nt] ; 

	    for (int i=0; i < np*nt; i++) {
		cf[i] = *pxtgg ;
		pxtgg++ ;
	    }
	    
	    cfpcossin(deg, dimc, cf) ;	    // Fourier transform 
	    
	    // Even harmonics
	    // --------------
	    double* pcf0 = cf0 ; 
	    double* pcf = cf ; 
	    for (int k=0; k<np-1; k += 4) {
		for(int j=0; j<2*nt; j++) {
		    *pcf0 = *pcf ; 
		    pcf0++ ; 
		    pcf++ ;
		}
		for(int j=0; j<2*nt; j++) {
		    *pcf0 = 0 ; 
		    pcf0++ ; 
		    pcf++ ;
		}
	    }	
	    if (np%4 == 0) {	    // particular case of np multiple of 4
		for(int j=0; j<2*nt; j++) {
		    *pcf0 = *pcf ; 
		    pcf0++ ; 
		    pcf++ ;
		}		
	    }

	    cipcossin(deg, dimc, deg, cf0, ff0) ; // Inverse Fourier transform 

	    pxtgg = xtgg.t ; 
	    for (int i=0; i < np*nt; i++) {
		*pxtgg = ff0[i] ;
		pxtgg++ ;
	    }
	    
	    // Odd harmonics
	    // -------------
	    pcf0 = cf0 ; 
	    pcf = cf ; 
	    for (int k=0; k<np-1; k += 4) {
		for(int j=0; j<2*nt; j++) {
		    *pcf0 = 0 ; 
		    pcf0++ ; 
		    pcf++ ;
		}
		for(int j=0; j<2*nt; j++) {
		    *pcf0 = *pcf ; 
		    pcf0++ ; 
		    pcf++ ;
		}
	    }	
	    if (np%4 == 0) {	    // particular case of np multiple of 4
		for(int j=0; j<2*nt; j++) {
		    *pcf0 = 0 ; 
		    pcf0++ ; 
		}		
	    }

	    cipcossin(deg, dimc, deg, cf0, ff0) ; // TF inverse

	    pxtff = xtff.t ; 
	    for (int i=0; i < np*nt; i++) {
		*pxtff = ff0[i] ;
		pxtff++ ;
	    }
	    
	    delete [] cf ; 
	    delete [] cf0 ; 
	    delete [] ff0 ; 
	    break ;
	}
	
	default : { 
	    cout << "Map_et::adapt: unknown phi basis !" << endl ; 
	    cout << "  base_p = " << base_p << endl ; 
	    abort() ; 
	}
    }
    
    // Computation of lambda and mu in the nucleus :
    // --------------------------------------------  
    
    double lambda = 0 ;	    // lambda is set to zero because F(theta, phi)
			    // must not have constant terms in the nucleus
     
    double mu = - fact_lamu * min(xtgg) ; 

    // Computation of alpha and beta in the nucleus : 
    // --------------------------------------------  

    nalpha[0] = r_eq - lambda - mu ; 
    nbeta[0] = 0 ; 
    
    // Computation of F_0, G_0 and {\tilde F_1} :
    // ------------------------------------------

    Mtbl nff(ff.get_mg()) ; 
    Mtbl ngg(gg.get_mg()) ; 
    nff.set_etat_qcq() ; 
    ngg.set_etat_qcq() ; 
    
    *(nff.t[0]) = ( xtff + lambda ) / nalpha[0] ; 
    *(ngg.t[0]) = ( xtgg + mu	  ) / nalpha[0] ; 
    xtff += xtgg ; 
    
    
    //----------------------------------------------------------------
    //   Adaptation of the shells
    //----------------------------------------------------------------
    
    double r_eqlm1 = r_eq ; 


    // Loop on the shells where the adaptation must be performed 
    
    for (int l=1; l<nzadapt; l++) {
    
	ent0 = ent_limit(l) ; 

	// Search for the equipotential surface ent = ent0 
	//    --->  l = l_ext(theta, phi) 
	//		xi = x_ext(theta, phi)
	// -----------------------------------------------

	(ent.va).equipot_outward(ent0, nz_search, precis, nitermax, niter0, 
				 l_ext, x_ext) ; 
		     
	niter = ( niter0 > niter ) ? niter0 : niter ;

	// The new mapping must be such that the found isosurface corresponds
	//  to xi = 1
    
	// Computation of r_ext(theta, phi) - r_eq ---> xtgg
	// -------------------------------------------------
    
	xtgg.set_etat_qcq() ; 
	assert(l_ext.get_etat() == ETATQCQ) ; 

	r_eq = val_r_jk(l, xi_max, j_bord, k_bord) ;

	pxtgg = xtgg.t ; 
	pl_ext = l_ext.t ;
	px_ext = x_ext.t ;
    
	for (int k=0; k<np; k++) {
	    for (int j=0; j<nt; j++) {

		*pxtgg = val_r_jk(*pl_ext, *px_ext, j, k) - r_eq  ;
	    			    
		// ... next one
		pl_ext++ ;
		px_ext++ ;
		pxtgg++ ;
	    }	   
	}	

	// Computation of lambda and mu in domain no. l :
	// --------------------------------------------  
    
	lambda = - fact_lamu * max(xtff) ; 
	mu = - fact_lamu * min(xtgg) ; 

	// Computation of alpha and beta in domain no. l : 
	// --------------------------------------------  

	nalpha[l] = .5 * ( r_eq - r_eqlm1 + lambda - mu ) ;
	nbeta[l] =  .5 * ( r_eq + r_eqlm1 - lambda - mu ) ; 
    
	// Computation of F_l, G_l and {\tilde F_(l+1)} :
	// ------------------------------------------
    
	*(nff.t[l]) = ( xtff + lambda ) / nalpha[l] ; 
	*(ngg.t[l]) = ( xtgg + mu     ) / nalpha[l] ; 
	xtff = xtgg ; 
    
	r_eqlm1 = r_eq ; 
	
    } // end of the loop on the shells where the adaptation must be performed
    
    //----------------------------------------------------------------
    //   Adaptation of the domain of index nzadapt
    //----------------------------------------------------------------
    
    if (mg->get_type_r(nzadapt) == UNSURR) {
	
	// Case of a compactified domain 
	// -----------------------------
	
	xtff = 1 / (xtff + r_eqlm1) - double(1) / r_eqlm1 ; 
	
	lambda = - fact_lamu * min(xtff) ;
	
	nalpha[nzadapt] = .5 * ( lambda - double(1) / r_eqlm1 ) ; 
	nbeta[nzadapt] = - nalpha[nzadapt] ;
	
	// Computation of F_nzadapt : 
	*(nff.t[nzadapt]) = ( xtff + lambda ) / nalpha[nzadapt] ; 
	
    }
    else {
	// Case of a non-compactified domain 
	// ----------------------------------
	
	r_eq = val_r_jk(nzadapt, xi_max, j_bord, k_bord) ;
	
	lambda = - fact_lamu * max(xtff) ;

	nalpha[nzadapt] = .5 * ( r_eq - r_eqlm1 + lambda ) ; 
	nbeta[nzadapt]  = .5 * ( r_eq + r_eqlm1 - lambda ) ; 

	// Computation of F_l : 

	*(nff.t[nzadapt]) = ( xtff + lambda ) / nalpha[nzadapt] ;	

    }

    // In all cases, G_nzadapt is set to zero : 
    ngg.t[nzadapt]->set_etat_zero() ; 

    //----------------------------------------------------------------
    //  Values of alpha, beta, F and G in the most external domains
    //   where the mapping is unchanged
    //----------------------------------------------------------------

    for (int l=nzadapt+1; l<nz; l++) {

	nalpha[l] = alpha[l] ; 
	nbeta[l]  = beta[l] ;

    }

    if (ff.get_etat() == ETATZERO) {
	for (int l=nzadapt+1; l<nz; l++) {
	    nff.t[l]->set_etat_zero() ; 
	}	
    }
    else {
	assert(ff.get_etat() == ETATQCQ) ; 
	assert(ff.c != 0x0) ; 
	assert(ff.c->get_etat() == ETATQCQ) ; 
	for (int l=nzadapt+1; l<nz; l++) {
	    *(nff.t[l]) = *(ff.c->t[l]) ; 
	}	
    }
    
    if (gg.get_etat() == ETATZERO) {
	for (int l=nzadapt+1; l<nz; l++) {
	    ngg.t[l]->set_etat_zero() ; 
	}	
    }
    else {
	assert(gg.get_etat() == ETATQCQ) ; 
	assert(gg.c != 0x0) ; 
	assert(gg.c->get_etat() == ETATQCQ) ; 
	for (int l=nzadapt+1; l<nz; l++) {
	    *(ngg.t[l]) = *(gg.c->t[l]) ; 
	}	
    }
    

    //=============================================================================
    //  Assignment of the mapping parameters alpha, beta, ff and gg to
    //   the newly computed values
    //=============================================================================

    for (int l=0; l<nz; l++) {
    
	if (mg->get_type_r(l) == UNSURR) {
	    alpha[l] = nalpha[l] / fact_echelle ; 
	    beta[l] = nbeta[l] / fact_echelle  ; 
	}
	else {
	    alpha[l] = fact_echelle * nalpha[l] ; 
	    beta[l] = fact_echelle * nbeta[l] ; 	    
	}
    
    }

    ff = nff ; 
    gg = ngg ; 
    
    ff.std_base_scal() ;    // Standard spectral bases for F 
    gg.std_base_scal() ;    // Standard spectral bases for G
        

    if (nbr_filtre !=0) {
      ff.coef() ;
      gg.coef() ;
      ff.set_etat_cf_qcq() ;
      gg.set_etat_cf_qcq() ;
      for (int l=0 ; l<nzadapt+1 ; l++)
	for (int k=0 ; k<np ; k++) 
	  for (int j=nt-nbr_filtre-1 ; j<nt ; j++) {
	    if (ff.c_cf->t[l]->get_etat() != ETATZERO)
	      ff.c_cf->set(l, k,j,0) = 0 ;
	  
	    if  (gg.c_cf->t[l]->get_etat() != ETATZERO)
	      gg.c_cf->set(l,k,j,0) = 0 ;
	  }
    }

    // The derived quantities must be reset
    // ------------------------------------
    
    reset_coord() ; 


    // Clean exit
    // ----------
    
    delete [] nalpha ;   
    delete [] nbeta ;  
        
}
