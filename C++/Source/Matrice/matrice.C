/*
 *  Methods of class Matrice
 *
 *   (see file matrice.h for documentation)
 *
 */

/*
 *   Copyright (c) 1999-2001 Philippe Grandclement
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


char matrice_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.13  2005/10/24 09:22:24  p_grandclement
 * addition of annule_hard for matrices
 *
 * Revision 1.12  2005/09/16 12:29:02  j_novak
 * New method del_deriv() and reorganization of band, lu, permute handling.
 *
 * Revision 1.11  2005/01/25 12:47:34  j_novak
 * Added some member arithmetic and operator=(Tbl).
 *
 * Revision 1.10  2004/12/29 12:27:36  j_novak
 * permute is now a Itbl* which array is sent directly to the LAPACK routines.
 * It is now possible to solve a general system (i.e. even if the Matrice
 * is not in a banded form).
 *
 * Revision 1.9  2004/10/05 15:44:19  j_novak
 * Minor speed enhancements.
 *
 * Revision 1.8  2004/08/24 09:14:43  p_grandclement
 * Addition of some new operators, like Poisson in 2d... It now requieres the
 * GSL library to work.
 *
 * Also, the way a variable change is stored by a Param_elliptic is changed and
 * no longer uses Change_var but rather 2 Scalars. The codes using that feature
 * will requiere some modification. (It should concern only the ones about monopoles)
 *
 * Revision 1.7  2003/12/19 16:21:44  j_novak
 * Shadow hunt
 *
 * Revision 1.6  2002/10/16 14:36:42  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.5  2002/09/24 10:51:16  e_gourgoulhon
 *
 * The case of a 1D Tbl in the constructor from Tbl is now taken into account
 * (resulting in a single-column matrix).
 *
 * Revision 1.4  2002/09/24 08:36:44  e_gourgoulhon
 *
 * Corrected error in output (operator<<) : permutted number of rows and columns
 *
 * Added matrix multiplication
 * Added function transpose()
 *
 * Revision 1.3  2002/09/09 13:00:39  e_gourgoulhon
 * Modification of declaration of Fortran 77 prototypes for
 * a better portability (in particular on IBM AIX systems):
 * All Fortran subroutine names are now written F77_* and are
 * defined in the new file C++/Include/proto_f77.h.
 *
 * Revision 1.2  2002/01/03 13:18:41  j_novak
 * Optimization: the members set(i,j) and operator(i,j) of class Matrice are
 * now defined inline. Matrice is a friend class of Tbl.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 2.9  1999/12/24  10:19:16  eric
 * Suppression des definitions de nbl et nbc lignes 149 et 150.
 *
 * Revision 2.8  1999/11/30  17:45:16  phil
 * changerment prototypage
 *
 * Revision 2.7  1999/10/12  15:49:16  phil
 * apres set band, lu et permute ne sont plus a jour ...
 *
 * Revision 2.6  1999/10/12  09:42:17  phil
 * retrour versian anterieure
 *
 * Revision 2.5  1999/10/12  09:39:07  phil
 * passage en const
 *
 * Revision 2.4  1999/10/11  09:35:07  phil
 * ajout de determinant et val_propre  + modif de operator= (const Matrice&)
 *
 * Revision 2.3  1999/10/05  17:02:46  phil
 * ajout de determinant et val_propre
 *
 * Revision 2.2  1999/04/13  13:57:23  phil
 * ajout proto.h
 *
 * Revision 2.1  1999/04/07  14:18:51  phil
 * optimisation egalite
 *
 * Revision 2.0  1999/04/07  14:10:05  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */


//fichiers includes
#include <stdlib.h>
#include "matrice.h"
#include "proto_f77.h"

//Destructeur logique

void Matrice::del_t() {
    if (std != 0x0) delete std ;
    std = 0x0 ; 
    del_deriv() ;
}

//Destructeur des quantites derivees

void Matrice::del_deriv() {
    if (band != 0x0) delete band ;
    if (lu != 0x0) delete lu ;
    if (permute != 0x0) delete permute ;
    band = 0x0 ;
    lu = 0x0 ;
    permute = 0x0 ;
}

//Manipulation des etats

void Matrice::set_etat_qcq() {
    std->set_etat_qcq() ;
    del_deriv() ;
    etat = ETATQCQ ;
}

void Matrice::set_etat_zero() {
    std->set_etat_zero() ;
    del_deriv() ;
    etat = ETATZERO ;
}

void Matrice::set_etat_nondef() {
    if (std != 0x0) std->set_etat_nondef() ;
    del_deriv() ;
    etat = ETATNONDEF ;
}

void Matrice::annule_hard() {
    std->set_etat_qcq() ;
    del_deriv() ;
    etat = ETATQCQ ;

    for (int i=0 ; i<std->get_dim(0) ; i++)
    	 for (int j=0 ; j<std->get_dim(1) ; j++)
    		std->set(i,j) = 0 ;
}
	 
// Constructeurs
Matrice::Matrice (int i, int j) {
    etat = ETATNONDEF ;
    std = new Tbl(i, j) ;
    kl = 0 ;
    ku = 0 ;
    band = 0x0 ;
    lu = 0x0 ;
    permute = 0x0 ;
}


Matrice::Matrice (const Matrice & source) {
    etat = source.etat ;
    kl = source.kl ;
    ku = source.ku ;
    std = new Tbl(*source.std) ;
    if (source.band != 0x0) band = new Tbl(*source.band) ;
    else band = 0x0 ;
    if (source.lu != 0x0) lu = new Tbl(*source.lu) ;
    else lu = 0x0 ;
    if (source.permute != 0x0) permute = new Itbl(*source.permute) ;
    else permute = 0x0 ;
}


Matrice::Matrice (const Tbl & source) {
    etat = source.get_etat() ;
    kl = 0 ;
    ku = 0 ;
    if (source.get_ndim() == 1) {     // column vector
        int n = source.get_taille() ;
        std = new Tbl(n,1) ;
        if (source.get_etat() == ETATZERO) {
	        std->set_etat_zero() ;
        }
	else  {
                assert( source.get_etat() == ETATQCQ ) ;
                std->set_etat_qcq() ;
                for (int i=0; i<n; i++) {
                        std->t[i] = source.t[i] ;
                }
	}
    }
    else {                        // 2D Tbl
        std = new Tbl(source) ;
    }
    band = 0x0 ;
    lu = 0x0 ;
    permute = 0x0 ;
}

// destructeur
Matrice::~Matrice() {
    del_t() ;
}

// Extraction des dimensions
int Matrice::get_dim(int i) const {
    return std->get_dim(i) ;
}

// affectation
void Matrice::operator= (double x) {
    if (x == 0 ) set_etat_zero() ;
    else {
    set_etat_qcq() ;
    *std = x ;
    }
}

void Matrice::operator= (const Matrice &source) {
    
    assert (std->get_dim(0) == source.std->get_dim(0)) ;
    assert (std->get_dim(1) == source.std->get_dim(1)) ;
    
    switch (source.etat) {
	case ETATNONDEF :
	    set_etat_nondef() ;
	    break ;
	case ETATZERO :
	    set_etat_zero() ;
	    break ;
	case ETATQCQ :
	    set_etat_qcq() ;
	    del_t() ;
	    
	    if (source.std != 0x0)
		std = new Tbl(*source.std) ;
	    
	    if (source.band != 0x0) {
		band = new Tbl(*source.band) ;
		ku = source.ku ;
		kl = source.kl ;
	    }
	    
	    if (source.lu != 0x0) {
		lu = new Tbl(*source.lu) ;
		permute = new Itbl(*source.permute) ;
	    }
	    break ;
	}
}

void Matrice::operator= (const Tbl &source) {
    
    assert (std->get_dim(0) == source.get_dim(0)) ;
    assert (std->get_dim(1) == source.get_dim(1)) ;
    
    switch (source.etat) {
	case ETATNONDEF :
	    set_etat_nondef() ;
	    break ;
	case ETATZERO :
	    set_etat_zero() ;
	    break ;
	case ETATQCQ :
	    set_etat_qcq() ;
	    del_t() ;
	    
	    assert (source.t != 0x0) ;
	    std = new Tbl(source) ;
	    break ;
	}
}


//Impression
ostream& operator<< (ostream& flux, const Matrice & source) {
    switch (source.std->get_etat()) {
	case ETATZERO :
	    flux << "Null matrix. " << endl ;
	    break ;
	case ETATNONDEF :
	    flux << "Undefined matrix. " << endl ;
	    break ;
	case ETATQCQ :
	    int nbl = source.std->get_dim(1) ;
	    int nbc = source.std->get_dim(0) ;
	    flux << "Matrix " << nbl << " * " << nbc << endl ;
	    for (int i=0 ; i<nbl ; i++) {
		for (int j=0 ; j<nbc ; j++)
		    flux << (*source.std)(i, j) << "  " ;
		flux << endl ;		
	    }
	}
	
    flux << endl ;
    
    if ((source.band != 0x0) && (source.band->get_etat() != ETATNONDEF)) {
	flux << "Matrix : " << source.ku << " upper diags. and  "
		 << source.kl << " lower diags." << endl ;
    }
	//    else flux << "Diagonalisation non faite." << endl ;
    
    if ((source.lu != 0x0) && (source.lu->get_etat() != ETATNONDEF))
	flux << "LU factorization done." << endl ;

return flux ;
}

// Passage matrice a bande : stockage LAPACK
void Matrice::set_band (int u, int l) const {
    if (band != 0x0) return ;
    else {
	int n = std->get_dim(0) ;
	assert (n == std->get_dim(1)) ;
	
	ku = u ; kl = l ;
	int ldab = 2*l+u+1 ;
	band = new Tbl(ldab*n) ;
	
	band->set_etat_qcq() ;
	
	for (int i=0 ; i<u ; i++)
	    for (int j=u-i ; j<n ; j++)
		band->set(j*ldab+i+l) = (*this)(j-u+i, j) ;
	
	for (int j=0 ; j<n ; j++)
	    band->set(j*ldab+u+l) = (*this)(j, j) ;
	
	for (int i=u+1 ; i<u+l+1 ; i++)
	    for (int j=0 ; j<n-i+u ; j++)
		band->set(j*ldab+i+l) = (*this) (i+j-u, j) ;
    }
    return ; 
}

//Decomposition UL : stockage LAPACK
void Matrice::set_lu() const {
    if (lu != 0x0) {
	assert (permute != 0x0) ;
	return ;
    }
    else {
	// Decomposition LU
	int n = std->get_dim(0) ;
	int ldab, info ;
	permute = new Itbl(n) ;
	permute->set_etat_qcq() ;

	// Cas d'une matrice a bandes
	if (band != 0x0) {
	    assert (band->get_etat() == ETATQCQ) ;
	    ldab = 2*kl+ku+1 ;
	    lu = new Tbl(*band) ;
	    
	    F77_dgbtrf(&n, &n, &kl, &ku, lu->t, &ldab, permute->t, &info) ;
	}
	else { // matrice generale
	    assert (std->get_etat() == ETATQCQ) ;
	    ldab = n ;
	    lu = new Tbl(*std) ;
	    
	    F77_dgetrf(&n, &n, lu->t, &ldab, permute->t, &info) ;
	}
    }
    return ;
}

// Solution de Ax = B : utilisation de LAPACK et decomposition lu.
Tbl Matrice::inverse (const Tbl& source) const {
    
    assert(lu != 0x0) ;
    assert(lu->get_etat() == ETATQCQ) ;
    assert(permute != 0x0) ;
    assert(permute->get_etat() == ETATQCQ) ;
       
    int n = source.get_dim(0) ;
    assert (get_dim(1) == n) ;
    int ldab, info ;
    char* trans ;
    int nrhs = 1 ;
    int ldb = n ;
        
    Tbl res(source) ;
    
    if (band != 0x0) { //Cas d'une matrice a bandes
	ldab = 2*kl+ku+1 ;
	trans = "N" ;
	F77_dgbtrs(trans, &n, &kl, &ku, &nrhs, lu->t,
		   &ldab, permute->t, res.t, &ldb, &info);
    }
    else { // Cas general
	ldab = n ;
	trans = "T" ; // stockage different entre le C et le fortran
	F77_dgetrs(trans, &n, &nrhs, lu->t, &ldab, permute->t,
		   res.t, &ldb, &info) ;
    }
    
    return res ;
}

// Renvoit les valeurs propres de la matrice (appel de LAPACK) :
Tbl Matrice::val_propre() const {
    
    assert (etat != ETATNONDEF) ;
    assert (std != 0x0) ;
    
    char* jobvl = "N" ;
    char* jobvr = "N" ;
    
    int n = get_dim(0) ;
    assert (n == get_dim(1)) ;
    
    double* a = new double [n*n] ;
    for (int i=0 ; i<n*n ; i++)
	a[i] = std->t[i] ;
	
    int lda = n ;
    double* wr = new double[n] ;
    double* wi = new double[n] ;
    
    int ldvl = 1 ;
    double* vl = 0x0 ;
    int ldvr = 1 ;
    double* vr = 0x0 ;
    
    int ldwork = 3*n ;
    double* work = new double[ldwork] ;
    
    int info ;
    
    F77_dgeev(jobvl, jobvr, &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr,
	    work, &ldwork, &info) ;
    
    Tbl result(2, n) ;
    result.set_etat_qcq() ;
    
    for (int i=0 ; i<n ; i++) {
	result.set(0, i) = wr[n-i-1] ;
	result.set(1, i) = wi[n-i-1] ;
	}
    
    delete [] wr ;
    delete [] wi ;
    delete [] a ;
    delete [] work ;
    
    return result ; 
    
}

// les valeurs vecteurs propres de la matrice (appel de LAPACK) :
Matrice Matrice::vect_propre() const {
    
    assert (etat != ETATNONDEF) ;
    assert (std != 0x0) ;
    
    char* jobvl = "V" ;
    char* jobvr = "N" ;
    
    int n = get_dim(0) ;
    assert (n == get_dim(1)) ;
    
    double* a = new double [n*n] ;
    for (int i=0 ; i<n*n ; i++)
	a[i] = std->t[i] ;
	
    int lda = n ;
    double* wr = new double[n] ;
    double* wi = new double[n] ;
    
    int ldvl = n ;
    double* vl = new double[ldvl*ldvl] ;
    int ldvr = 1 ;
    double* vr = 0x0 ;
    
    int ldwork = 4*n ;
    double* work = new double[ldwork] ;
    
    int info ;
    
    F77_dgeev(jobvl, jobvr, &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr,
        work, &ldwork, &info) ;
    
   
    Matrice res (n,n) ;
    res.set_etat_qcq() ;

    int conte = 0 ;
    for (int i=0 ; i<n ; i++)
      for (int j=0 ; j<n ; j++) {
	res.set(j,n-i-1) = vl[conte] ;
	conte ++ ;
      }
    
    delete [] wr ;
    delete [] wi ;
    delete [] a ;
    delete [] work ;
    delete [] vl ;

    return res ;
}

// Calcul le determinant :
double Matrice::determinant() const {
    
    int n = get_dim(0) ;
    assert(n == get_dim(1)) ;
    
    Tbl valp(val_propre()) ;
    double result = 1 ;
    for (int i = 0 ; i<n ; i++)
	if (valp(1, i) == 0)
	    result *= valp(0, i) ;
	else {
	    result*= valp(0, i)*valp(0, i)+valp(1, i)*valp(1, i) ;
	    i++ ;
	}
    return result ;
}

// Transposee
Matrice Matrice::transpose() const {

	int nbl = std->get_dim(1) ;
	int nbc = std->get_dim(0) ;

	Matrice resu(nbc, nbl) ;
	
	if (etat == ETATZERO) {
		resu.set_etat_zero() ;
	}
	else{
		assert(etat == ETATQCQ) ;
		resu.set_etat_qcq() ;
		for (int i=0; i<nbc; i++) {
			for (int j=0; j<nbl; j++) {
				resu.set(i,j) = (*std)(j,i) ;
			}
		}
	}
	return resu ;
}


// Operateurs d'arithmetique
void Matrice::operator+=(const Matrice& a) {
    assert((std != 0x0)&&(a.std != 0x0)) ;
    std->operator+=(*a.std) ;
}

void Matrice::operator-=(const Matrice& a) {
    assert((std != 0x0)&&(a.std != 0x0)) ;
    std->operator-=(*a.std) ;
}

void Matrice::operator+=(double x) {
    assert(std != 0x0);
    std->operator+=(x) ;
}

void Matrice::operator-=(double x) {
    assert(std != 0x0);
    std->operator-=(x) ;
}

void Matrice::operator*=(double x) {
    assert(std != 0x0);
    std->operator*=(x) ;
}

void Matrice::operator/=(double x) {
    assert(std != 0x0);
    assert(x != 0) ;
    std->operator/=(x) ;
}

// Operateurs d'arithmetique non membres
Matrice operator+ (const Matrice& a, const Matrice& b) {
    assert((a.std != 0x0) && (b.std != 0x0)) ;
    Matrice res(*a.std+*b.std) ;
    return res ;
}

Matrice operator- (const Matrice& a, const Matrice& b) {
    assert((a.std != 0x0) && (b.std != 0x0)) ;
    Matrice res(*a.std-*b.std) ;
    return res ;
}

Matrice operator* (const Matrice& a, double x) {
    assert(a.std != 0x0) ;
    Matrice res(*a.std*x);
    return res ;
}

Matrice operator* (double x, const Matrice& a) {
    assert(a.std != 0x0) ;
    Matrice res(*a.std*x);
    return res ;
}

Matrice operator* (const Matrice& aa, const Matrice& bb) {

	int nbla = aa.std->get_dim(1) ;
	int nbca = aa.std->get_dim(0) ;
#ifndef NDEBUG
	int nblb = bb.std->get_dim(1) ;
#endif
	int nbcb = bb.std->get_dim(0) ;
	
	assert( nbca == nblb ) ;
	
	Matrice resu(nbla, nbcb) ;
	
	if ( (aa.get_etat() == ETATZERO) || (bb.get_etat() == ETATZERO) ) {
		resu.set_etat_zero() ;
	}
	else {
		assert( aa.get_etat() == ETATQCQ ) ;
		assert( bb.get_etat() == ETATQCQ ) ;
		resu.set_etat_qcq() ;
		for (int i=0; i<nbla; i++) {
			for (int j=0; j<nbcb; j++) {
				double sum = 0 ;
				for (int k=0; k<nbca; k++) {
					sum += aa(i,k) * bb(k, j) ;
				}
				resu.set(i,j) = sum ;
			}
			
		}
	}

    return resu ;
}

Matrice operator/ (const Matrice& a, double x) {
    assert (x != 0) ;
    assert(a.std != 0x0) ;
    Matrice res(*a.std/x);
    return res ;
}
