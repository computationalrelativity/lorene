/*
 * Member function of the Scalar class for initiating a Scalar from 
 * a Scalar defined on another mapping.
 */

/*
 *   Copyright (c) 2003 Eric Gourgoulhon & Jerome Novak
 *   Copyright (c) 1999-2001 Eric Gourgoulhon (Cmp version)
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


char scalar_import_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2003/10/10 15:57:29  j_novak
 * Added the state one (ETATUN) to the class Scalar
 *
 * Revision 1.2  2003/10/01 13:04:44  e_gourgoulhon
 * The method Tensor::get_mp() returns now a reference (and not
 * a pointer) onto a mapping.
 *
 * Revision 1.1  2003/09/25 09:07:05  j_novak
 * Added the functions for importing from another mapping (to be tested).
 *
 *
 * $Header$
 *
 */

// Headers C
#include <math.h>

// Headers Lorene
#include "tensor.h"
#include "param.h"
#include "nbr_spx.h"

			//-------------------------------//
			//  Importation in all domains   //
			//-------------------------------//

void Scalar::import(const Scalar& ci) {
    
    int nz = mp->get_mg()->get_nzone() ; 
    
    import(nz, ci) ; 
        
}

			//--------------------------------------//
			//  Importation in inner domains only   //
			//--------------------------------------//

void Scalar::import(int nzet, const Scalar& cm_d) {
    
    const Map* mp_d = &(cm_d.get_mp()) ; // Departure mapping

    // Trivial case : mappings identical !
    // -----------------------------------
    
    if (mp_d == mp) {
	*this = cm_d ; 
	return ; 
    }
    
    // Relative orientation of the two mappings
    // ----------------------------------------
 
    int align_rel = (mp->get_bvect_cart()).get_align()  
		    * (mp_d->get_bvect_cart()).get_align() ;
		    
    switch (align_rel) {

	case 1 : {  // the two mappings have aligned Cartesian axis
	    import_align(nzet, cm_d) ; 
	    break ; 
	}

	case -1 : {  // the two mappings have anti-aligned Cartesian axis
	    import_anti(nzet, cm_d) ; 
	    break ; 
	}

	case 0 : {  // general case 
	    import_gal(nzet, cm_d) ; 
	    break ; 
	}

	default : {  
	    cout << "Scalar::import : unexpected value of align_rel : " 
		 << align_rel << endl ; 
	    abort() ; 
	    break ; 
	}

    }
    
}    

			//--------------------------------------//
			//     General case (axis not aligned)  //
			//--------------------------------------//


void Scalar::import_gal(int nzet, const Scalar& cm_d) {
    
    const Map* mp_d = &(cm_d.get_mp()) ; // Departure mapping

    // Trivial case : mappings identical !
    // -----------------------------------
    
    if (mp_d == mp) {
	*this = cm_d ; 
	return ; 
    }
    
    // Another trivial case : null Scalar
    // -------------------------------
    
    if (cm_d.get_etat() == ETATZERO) {
	set_etat_zero() ; 
	return ; 
    }

    if (cm_d.get_etat() == ETATUN) {
	set_etat_one() ; 
	return ; 
    }

    // Protections
    // -----------
    
    assert(cm_d.get_etat() != ETATNONDEF) ; 

    if (cm_d.get_dzpuis() != 0) {
	cout << 
	"Scalar::import : the dzpuis of the Scalar to be imported must be zero !"
	     << endl ; 
	abort() ; 
    }


    const Mg3d* mg_a = mp->get_mg() ; 
    int nz_a = mg_a->get_nzone() ; 
    assert(nzet <= nz_a) ;     


    // General case :
    // -------------
    assert(cm_d.get_etat() == ETATQCQ) ;
    const Valeur& va_d = cm_d.get_spectral_va() ; 
    va_d.coef() ;		// The coefficients are required
    

    // Preparations for storing the result in *this 
    // --------------------------------------------
    del_t() ;	// delete all previously computed derived quantities
    
    set_etat_qcq() ;		// Set the state to ETATQCQ
    
    va.set_etat_c_qcq() ;	// Allocates the memory for the Mtbl va.c
				//  if it does not exist already
    va.c->set_etat_qcq() ;	// Allocates the memory for the Tbl's in each
				//  domain if they do not exist already


    // Absolute coordinates of the origin of the Departure mapping 
    double xo_d = mp_d->get_ori_x() ;	  
    double yo_d = mp_d->get_ori_y() ;	  
    double zo_d = mp_d->get_ori_z() ;	  

    // Orientation relative to the Absolute frame of the Departure mapping
    double rot_phi_d = mp_d->get_rot_phi() ; 

    // Orientation relative to the Absolute frame of the Arrival mapping
    double rot_phi_a = mp->get_rot_phi() ; 
    
    // r, theta, phi, X, Y and Z on the Arrival mapping 
    //  update of the corresponding Coord's if necessary
    
    if ( (mp->r).c == 0x0 ) (mp->r).fait() ; 
    if ( (mp->tet).c == 0x0 ) (mp->tet).fait() ; 
    if ( (mp->phi).c == 0x0 ) (mp->phi).fait() ; 
    if ( (mp->xa).c == 0x0 ) (mp->xa).fait() ; 
    if ( (mp->ya).c == 0x0 ) (mp->ya).fait() ; 
    if ( (mp->za).c == 0x0 ) (mp->za).fait() ; 

    const Mtbl* mr_a = (mp->r).c ; 
    const Mtbl* mtet_a = (mp->tet).c ; 
    const Mtbl* mphi_a = (mp->phi).c ; 
    const Mtbl* mxa_a = (mp->xa).c ;
    const Mtbl* mya_a = (mp->ya).c ;
    const Mtbl* mza_a = (mp->za).c ;

    Param par_precis ;	    // Required precision in the method Map::val_lx
    int nitermax = 100 ;    // Maximum number of iteration in the secant method
    int niter ; 
    double precis = 1e-15 ; // Absolute precision in the secant method
    par_precis.add_int(nitermax) ;	
    par_precis.add_int_mod(niter) ; 
    par_precis.add_double(precis) ; 


    // Loop of the Arrival domains where the computation is to be performed
    // --------------------------------------------------------------------
    
    for (int l=0; l < nzet; l++) {
	
	int nr = mg_a->get_nr(l) ;
	int nt = mg_a->get_nt(l) ;
	int np = mg_a->get_np(l) ;
	
	
	const double* pr_a = mr_a->t[l]->t ;	  // Pointer on the values of r
	const double* ptet_a = mtet_a->t[l]->t ;  // Pointer on the values of theta
	const double* pphi_a = mphi_a->t[l]->t ;  // Pointer on the values of phi
	const double* pxa_a = mxa_a->t[l]->t ;	  // Pointer on the values of X
	const double* pya_a = mya_a->t[l]->t ;	  // Pointer on the values of Y
	const double* pza_a = mza_a->t[l]->t ;	  // Pointer on the values of Z
	
	(va.c->t[l])->set_etat_qcq() ;	     // Allocates the array of double to
					     //  store the result 
	double* ptx = (va.c->t[l])->t ;	     // Pointer on the allocated array
	
	
	// Loop on all the grid points in the considered arrival domain:
		
	for (int k=0 ; k<np ; k++) {
	    for (int j=0 ; j<nt ; j++) {
		for (int i=0 ; i<nr ; i++) {

		    double r = *pr_a ; 
		    double rd, tetd, phid ;
		    if (r == __infinity) {
			rd = r ;
			tetd = *ptet_a ;
			phid = *pphi_a + rot_phi_a - rot_phi_d ;
			if (phid < 0) phid += 2*M_PI ;
		    }
		    else {
			// Coordinates in a Cartesian frame centered on
			//  the Departure mapping and whose axes are
			//  parallel to those of the Absolue Frame
			double xd = *pxa_a - xo_d ;
			double yd = *pya_a - yo_d ;
			double zd = *pza_a - zo_d ;
		    
			// Spherical coordinates on the Departure mapping
			double rhod2 = xd*xd + yd*yd ; 
			double rhod = sqrt( rhod2 ) ;
			rd = sqrt(rhod2 + zd*zd) ;
			tetd = atan2(rhod, zd) ;
			phid = atan2(yd, xd) - rot_phi_d ; // (rotation)
			if (phid < 0) phid += 2*M_PI ;
		    }		    


		    // NB: to increase the efficiency, the method Scalar::val_point
		    //  is not invoked; the method Mtbl_cf::val_point is 
		    //  called directly instead. 

		    // Value of the grid coordinates (l,xi) corresponding to
		    //  (rd,tetd,phid) :
		    
		    int ld ;	    // domain index
		    double xxd ;    // radial coordinate xi in [0,1] or [-1,1]
		    mp_d->val_lx(rd, tetd, phid, par_precis, ld, xxd) ;

		    // Value of the Departure Scalar at the obtained point:
		    *ptx = va_d.c_cf->val_point(ld, xxd, tetd, phid) ;

		    // Next point : 
		    ptx++ ;   
		    pr_a++ ;  
		    ptet_a++ ;  
		    pphi_a++ ; 
		    pxa_a++ ;  
		    pya_a++ ;  
		    pza_a++ ;  

		}
	    }
	}	
	
		
    }	// End of the loop on the Arrival domains
    
    // In the remaining domains, *this is set to zero:
    // ----------------------------------------------
    
    if (nzet < nz_a) {
	annule(nzet, nz_a - 1) ; 
    }
    
    // Treatment of dzpuis
    // -------------------
    
    
    set_dzpuis(0) ; 

}


		//-----------------------------------------//
		//   Case of Cartesian axis anti-aligned   //
		//-----------------------------------------//


void Scalar::import_anti(int nzet, const Scalar& cm_d) {
    
    // Trivial case : null Scalar
    // ------------------------
    
    if (cm_d.get_etat() == ETATZERO) {
	set_etat_zero() ; 
	return ; 
    }

    if (cm_d.get_etat() == ETATUN) {
	set_etat_one() ; 
	return ; 
    }

    const Map* mp_d = &(cm_d.get_mp()) ; // Departure mapping

    // Protections
    // -----------
    int align = (mp->get_bvect_cart()).get_align() ;

    assert( align  * (mp_d->get_bvect_cart()).get_align() == -1 ) ; 
    
    assert(cm_d.get_etat() == ETATQCQ) ;

    if (cm_d.get_dzpuis() != 0) {
	cout << 
	"Scalar::import : the dzpuis of the Scalar to be imported must be zero !"
	     << endl ; 
	abort() ; 
    }


    const Mg3d* mg_a = mp->get_mg() ; 
    int nz_a = mg_a->get_nzone() ; 
    assert(nzet <= nz_a) ;     

    const Valeur& va_d = cm_d.get_spectral_va() ; 
    va_d.coef() ;		// The coefficients are required
    

    // Preparations for storing the result in *this 
    // --------------------------------------------
    del_t() ;	// delete all previously computed derived quantities
    
    set_etat_qcq() ;		// Set the state to ETATQCQ
    
    va.set_etat_c_qcq() ;	// Allocates the memory for the Mtbl va.c
				//  if it does not exist already
    va.c->set_etat_qcq() ;	// Allocates the memory for the Tbl's in each
				//  domain if they do not exist already


    // Departure (x,y,z) coordinates of the origin of the Arrival mapping :

    double xx_a, yy_a, zz_a ; 
    if (align == 1) {
	xx_a = mp_d->get_ori_x() - mp->get_ori_x() ; 
	yy_a = mp_d->get_ori_y() - mp->get_ori_y() ; 
    }
    else {
	xx_a = mp->get_ori_x() - mp_d->get_ori_x() ; 
	yy_a = mp->get_ori_y() - mp_d->get_ori_y() ; 
    }
    zz_a = mp->get_ori_z() - mp_d->get_ori_z() ; 

    
    // r, theta, phi, x, y and z on the Arrival mapping 
    //  update of the corresponding Coord's if necessary
    
    if ( (mp->r).c == 0x0 ) (mp->r).fait() ; 
    if ( (mp->tet).c == 0x0 ) (mp->tet).fait() ; 
    if ( (mp->phi).c == 0x0 ) (mp->phi).fait() ; 
    if ( (mp->x).c == 0x0 ) (mp->x).fait() ; 
    if ( (mp->y).c == 0x0 ) (mp->y).fait() ; 
    if ( (mp->z).c == 0x0 ) (mp->z).fait() ; 

    const Mtbl* mr_a = (mp->r).c ; 
    const Mtbl* mtet_a = (mp->tet).c ; 
    const Mtbl* mphi_a = (mp->phi).c ; 
    const Mtbl* mx_a = (mp->x).c ;
    const Mtbl* my_a = (mp->y).c ;
    const Mtbl* mz_a = (mp->z).c ;

    Param par_precis ;	    // Required precision in the method Map::val_lx
    int nitermax = 100 ;    // Maximum number of iteration in the secant method
    int niter ; 
    double precis = 1e-15 ; // Absolute precision in the secant method
    par_precis.add_int(nitermax) ;	
    par_precis.add_int_mod(niter) ; 
    par_precis.add_double(precis) ; 


    // Loop of the Arrival domains where the computation is to be performed
    // --------------------------------------------------------------------
    
    for (int l=0; l < nzet; l++) {
	
	int nr = mg_a->get_nr(l) ;
	int nt = mg_a->get_nt(l) ;
	int np = mg_a->get_np(l) ;
	
	
	const double* pr_a = mr_a->t[l]->t ;	  // Pointer on the values of r
	const double* ptet_a = mtet_a->t[l]->t ;  // Pointer on the values of theta
	const double* pphi_a = mphi_a->t[l]->t ;  // Pointer on the values of phi
	const double* px_a = mx_a->t[l]->t ;	  // Pointer on the values of X
	const double* py_a = my_a->t[l]->t ;	  // Pointer on the values of Y
	const double* pz_a = mz_a->t[l]->t ;	  // Pointer on the values of Z
	
	(va.c->t[l])->set_etat_qcq() ;	     // Allocates the array of double to
					     //  store the result 
	double* ptx = (va.c->t[l])->t ;	     // Pointer on the allocated array
	
	
	// Loop on all the grid points in the considered arrival domain:
		
	for (int k=0 ; k<np ; k++) {
	    for (int j=0 ; j<nt ; j++) {
		for (int i=0 ; i<nr ; i++) {

		    double r = *pr_a ; 
		    double rd, tetd, phid ;
		    if (r == __infinity) {
			rd = r ;
			tetd = *ptet_a ;
			phid = *pphi_a + M_PI ;
			if (phid < 0) phid += 2*M_PI ;
		    }
		    else {

			// Cartesian coordinates on the Departure mapping
			double xd = - *px_a + xx_a ; 
			double yd = - *py_a + yy_a ; 
			double zd = *pz_a + zz_a ; 

			// Spherical coordinates on the Departure mapping
			double rhod2 = xd*xd + yd*yd ; 
			double rhod = sqrt( rhod2 ) ;
			rd = sqrt(rhod2 + zd*zd) ;
			tetd = atan2(rhod, zd) ;
			phid = atan2(yd, xd) ; 
			if (phid < 0) phid += 2*M_PI ;
		    }		    


		    // NB: to increase the efficiency, the method Scalar::val_point
		    //  is not invoked; the method Mtbl_cf::val_point is 
		    //  called directly instead. 

		    // Value of the grid coordinates (l,xi) corresponding to
		    //  (rd,tetd,phid) :
		    
		    int ld ;	    // domain index
		    double xxd ;    // radial coordinate xi in [0,1] or [-1,1]
		    mp_d->val_lx(rd, tetd, phid, par_precis, ld, xxd) ;

		    // Value of the Departure Scalar at the obtained point:
		    *ptx = va_d.c_cf->val_point(ld, xxd, tetd, phid) ;

		    // Next point : 
		    ptx++ ;   
		    pr_a++ ;  
		    ptet_a++ ;  
		    pphi_a++ ; 
		    px_a++ ;  
		    py_a++ ;  
		    pz_a++ ;  

		}
	    }
	}	
	
		
    }	// End of the loop on the Arrival domains
    
    // In the remaining domains, *this is set to zero:
    // ----------------------------------------------
    
    if (nzet < nz_a) {
	annule(nzet, nz_a - 1) ; 
    }
    
    // Treatment of dzpuis
    // -------------------
    
    
    set_dzpuis(0) ; 

}


		//-------------------------------------//
		//   Case of aligned Cartesian axis    //
		//-------------------------------------//


void Scalar::import_align(int nzet, const Scalar& cm_d) {
    
    // Trivial case : null Scalar
    // ------------------------
    
    if (cm_d.get_etat() == ETATZERO) {
	set_etat_zero() ; 
	return ; 
    }
    if (cm_d.get_etat() == ETATUN) {
	set_etat_one() ; 
	return ; 
    }

    const Map* mp_d = &(cm_d.get_mp()) ; // Departure mapping

    // Protections
    // -----------
    int align = (mp->get_bvect_cart()).get_align() ;

    assert( align  * (mp_d->get_bvect_cart()).get_align() == 1 ) ; 
    
    assert(cm_d.get_etat() == ETATQCQ) ;

    if (cm_d.get_dzpuis() != 0) {
	cout << 
	"Scalar::import : the dzpuis of the Scalar to be imported must be zero !"
	     << endl ; 
	abort() ; 
    }


    const Mg3d* mg_a = mp->get_mg() ; 
    int nz_a = mg_a->get_nzone() ; 
    assert(nzet <= nz_a) ;     

    const Valeur& va_d = cm_d.get_spectral_va() ; 
    va_d.coef() ;		// The coefficients are required
    

    // Preparations for storing the result in *this 
    // --------------------------------------------
    del_t() ;	// delete all previously computed derived quantities
    
    set_etat_qcq() ;		// Set the state to ETATQCQ
    
    va.set_etat_c_qcq() ;	// Allocates the memory for the Mtbl va.c
				//  if it does not exist already
    va.c->set_etat_qcq() ;	// Allocates the memory for the Tbl's in each
				//  domain if they do not exist already


    // Departure (x,y,z) coordinates of the origin of the Arrival mapping :

    double xx_a, yy_a, zz_a ; 
    if (align == 1) {
	xx_a = mp->get_ori_x() - mp_d->get_ori_x() ; 
	yy_a = mp->get_ori_y() - mp_d->get_ori_y() ; 
    }
    else {
	xx_a = mp_d->get_ori_x() - mp->get_ori_x() ; 
	yy_a = mp_d->get_ori_y() - mp->get_ori_y() ; 
    }
    zz_a = mp->get_ori_z() - mp_d->get_ori_z() ; 

    
    // r, theta, phi, x, y and z on the Arrival mapping 
    //  update of the corresponding Coord's if necessary
    
    if ( (mp->r).c == 0x0 ) (mp->r).fait() ; 
    if ( (mp->tet).c == 0x0 ) (mp->tet).fait() ; 
    if ( (mp->phi).c == 0x0 ) (mp->phi).fait() ; 
    if ( (mp->x).c == 0x0 ) (mp->x).fait() ; 
    if ( (mp->y).c == 0x0 ) (mp->y).fait() ; 
    if ( (mp->z).c == 0x0 ) (mp->z).fait() ; 

    const Mtbl* mr_a = (mp->r).c ; 
    const Mtbl* mtet_a = (mp->tet).c ; 
    const Mtbl* mphi_a = (mp->phi).c ; 
    const Mtbl* mx_a = (mp->x).c ;
    const Mtbl* my_a = (mp->y).c ;
    const Mtbl* mz_a = (mp->z).c ;

    Param par_precis ;	    // Required precision in the method Map::val_lx
    int nitermax = 100 ;    // Maximum number of iteration in the secant method
    int niter ; 
    double precis = 1e-15 ; // Absolute precision in the secant method
    par_precis.add_int(nitermax) ;	
    par_precis.add_int_mod(niter) ; 
    par_precis.add_double(precis) ; 


    // Loop of the Arrival domains where the computation is to be performed
    // --------------------------------------------------------------------
    
    for (int l=0; l < nzet; l++) {
	
	int nr = mg_a->get_nr(l) ;
	int nt = mg_a->get_nt(l) ;
	int np = mg_a->get_np(l) ;
	
	
	const double* pr_a = mr_a->t[l]->t ;	  // Pointer on the values of r
	const double* ptet_a = mtet_a->t[l]->t ;  // Pointer on the values of theta
	const double* pphi_a = mphi_a->t[l]->t ;  // Pointer on the values of phi
	const double* px_a = mx_a->t[l]->t ;	  // Pointer on the values of X
	const double* py_a = my_a->t[l]->t ;	  // Pointer on the values of Y
	const double* pz_a = mz_a->t[l]->t ;	  // Pointer on the values of Z
	
	(va.c->t[l])->set_etat_qcq() ;	     // Allocates the array of double to
					     //  store the result 
	double* ptx = (va.c->t[l])->t ;	     // Pointer on the allocated array
	
	
	// Loop on all the grid points in the considered arrival domain:
		
	for (int k=0 ; k<np ; k++) {
	    for (int j=0 ; j<nt ; j++) {
		for (int i=0 ; i<nr ; i++) {

		    double r = *pr_a ; 
		    double rd, tetd, phid ;
		    if (r == __infinity) {
			rd = r ;
			tetd = *ptet_a ;
			phid = *pphi_a ;
		    }
		    else {

			// Cartesian coordinates on the Departure mapping
			double xd = *px_a + xx_a ; 
			double yd = *py_a + yy_a ; 
			double zd = *pz_a + zz_a ; 

			// Spherical coordinates on the Departure mapping
			double rhod2 = xd*xd + yd*yd ; 
			double rhod = sqrt( rhod2 ) ;
			rd = sqrt(rhod2 + zd*zd) ;
			tetd = atan2(rhod, zd) ;
			phid = atan2(yd, xd) ; 
			if (phid < 0) phid += 2*M_PI ;
		    }		    


		    // NB: to increase the efficiency, the method Scalar::val_point
		    //  is not invoked; the method Mtbl_cf::val_point is 
		    //  called directly instead. 

		    // Value of the grid coordinates (l,xi) corresponding to
		    //  (rd,tetd,phid) :
		    
		    int ld ;	    // domain index
		    double xxd ;    // radial coordinate xi in [0,1] or [-1,1]
		    mp_d->val_lx(rd, tetd, phid, par_precis, ld, xxd) ;

		    // Value of the Departure Scalar at the obtained point:
		    *ptx = va_d.c_cf->val_point(ld, xxd, tetd, phid) ;

		    // Next point : 
		    ptx++ ;   
		    pr_a++ ;  
		    ptet_a++ ;  
		    pphi_a++ ; 
		    px_a++ ;  
		    py_a++ ;  
		    pz_a++ ;  

		}
	    }
	}	
	
		
    }	// End of the loop on the Arrival domains
    
    // In the remaining domains, *this is set to zero:
    // ----------------------------------------------
    
    if (nzet < nz_a) {
	annule(nzet, nz_a - 1) ; 
    }
    
    // Treatment of dzpuis
    // -------------------
    
    
    set_dzpuis(0) ; 

}
