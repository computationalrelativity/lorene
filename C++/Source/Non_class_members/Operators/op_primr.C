/*
 *  Computation of primitive in a single domain 
 *
 *
 */

/*
 *   Copyright (c) 2004  Eric Gourgoulhon. 
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License version 2
 *   as published by the Free Software Foundation.
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

char op_primr_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.4  2004/11/23 15:16:01  m_forot
 *
 * Added the bases for the cases without any equatorial symmetry
 *  (T_COSSIN_C, T_COSSIN_S, T_LEG, R_CHEBPI_P, R_CHEBPI_I).
 *
 * Revision 1.3  2004/10/12 09:58:24  j_novak
 * Better memory management.
 *
 * Revision 1.2  2004/06/14 15:24:57  e_gourgoulhon
 * First operationnal version (tested).
 *
 * Revision 1.1  2004/06/13 21:33:13  e_gourgoulhon
 * First version.
 *
 *
 * $Header$
 *
 */


// C headers
#include <stdlib.h>

// Lorene headers
#include "tbl.h"

// Unexpected case
//----------------
void _primr_pas_prevu(const Tbl&, int bin, const Tbl&, Tbl&, int&, Tbl& ) {

    cout << "Unexpected basis in primr : basis = " << hex << bin << endl ; 
    abort() ;  

}

// case R_CHEB
//------------
void _primr_r_cheb(const Tbl& tin, int bin, const Tbl& valm1, Tbl& tout, 
        int& bout, Tbl& valp1) {

    assert(tin.dim == tout.dim) ;   

    // Output spectral basis
    bout = bin ; 

    // Number of coefficients
    int nr = tin.get_dim(0) ;	    
    int nt = tin.get_dim(1) ;	    
    int np = tin.get_dim(2) - 2 ;
    int borne_phi = np + 1 ; 
    if (np == 1) borne_phi = 1 ; 
    
    // Case of a zero input
    // --------------------
    if (tin.get_etat() == ETATZERO) {
        if (valm1.get_etat() == ETATZERO) {
            tout.set_etat_zero() ; 
            valp1.set_etat_zero() ;
            return ; 
        }
        else {
            assert(valm1.get_etat() == ETATQCQ) ; 
            tout.set_etat_qcq() ; 
            valp1.set_etat_qcq() ; 
            double* xco = tout.t ;	
            for (int k=0 ; k< borne_phi ; k++) {
	        if (k==1) {     // jump over the coefficient of sin(0*phi) 
	            xco += nr*nt ;
	        }
	        else {
	            for (int j=0 ; j<nt ; j++) {
                        xco[0] = valm1(k,j) ;  // constant value = boundary value
                        for (int i=1; i<nr; i++) xco[i] = 0 ; 
                        valp1.set(k,j) = xco[0]  ;
	                xco += nr ;
                    }
                }
            }
            return ; 
        }
    }

    // Case of a non-zero input
    // ------------------------

    assert(tin.get_etat() == ETATQCQ ) ; 
    tout.annule_hard() ; 
    valp1.annule_hard() ; 
    
    const double* xci = tin.t ;	
    double* xco = tout.t ;	

    for (int k=0 ; k< borne_phi ; k++) {
        if (k==1) {     // jump over the coefficient of sin(0*phi) 
            xci += nr*nt ;
            xco += nr*nt ;
        }
        else {
            for (int j=0 ; j<nt ; j++) {
            
                xco[1] = xci[0] - 0.5 * xci[2] ; // special case i = 1

                for (int i=2; i<nr-2; i++) {
                    xco[i] = (xci[i-1] - xci[i+1]) / double(2*i) ; 
                }
                
                xco[nr-2] = xci[nr-3] / double(2*nr - 4) ; 
                xco[nr-1] = xci[nr-2] / double(2*nr - 2) ; 

                // Determination of the T_0 coefficient by maching with
                // provided value at xi = - 1 : 
                double som = - xco[1] ; 
                for (int i=2; i<nr; i+=2) som += xco[i] ; 
                for (int i=3; i<nr; i+=2) som -= xco[i] ; 
                xco[0] = valm1(k,j) - som ;                 

                // Value of primitive at xi = + 1 : 
                som = xco[0] ; 
                for (int i=1; i<nr; i++) som += xco[i] ;
                valp1.set(k,j) = som ; 
                
                xci += nr ;
                xco += nr ;
            }   // end of theta loop
        }   
    }   // end of phi loop

}



// case R_CHEBP
//-------------
void _primr_r_chebp(const Tbl& tin, int bin, const Tbl&, Tbl& tout, 
        int& bout, Tbl& valp1) {

    assert(tin.dim == tout.dim) ;   

    // Output spectral basis
    int base_t = bin & MSQ_T ;
    int base_p = bin & MSQ_P ;
    bout = base_p | base_t | R_CHEBI ;

    // Number of coefficients
    int nr = tin.get_dim(0) ;	    
    int nt = tin.get_dim(1) ;	    
    int np = tin.get_dim(2) - 2 ;
    int borne_phi = np + 1 ; 
    if (np == 1) borne_phi = 1 ; 
        
    // Case of a zero input
    // --------------------
    if (tin.get_etat() == ETATZERO) {
        tout.set_etat_zero() ; 
        valp1.set_etat_zero() ; 
        return ; 
    }

    // Case of a non-zero input
    // ------------------------

    assert(tin.get_etat() == ETATQCQ ) ; 
    tout.annule_hard() ; 
    valp1.annule_hard() ; 
    
    const double* xci = tin.t ;	
    double* xco = tout.t ;	

    for (int k=0 ; k< borne_phi ; k++) {
        if (k==1) {     // jump over the coefficient of sin(0*phi) 
            xci += nr*nt ;
            xco += nr*nt ;
        }
        else {
            for (int j=0 ; j<nt ; j++) {
            
                xco[0] = xci[0] - 0.5*xci[1] ; // special case i = 0

                for (int i=1; i<nr-2; i++) {
                    xco[i] = (xci[i] - xci[i+1]) / double(4*i+2) ; 
                }
                
                xco[nr-2] = xci[nr-2] / double(4*nr - 6) ; 
                xco[nr-1] = 0 ; 

                // Value of primitive at xi = + 1 : 
                double som = xco[0] ; 
                for (int i=1; i<nr; i++) som += xco[i] ;
                valp1.set(k,j) = som ; 
                
                xci += nr ;
                xco += nr ;
            }   // end of theta loop
        }   
    }   // end of phi loop

}


// case R_CHEBI
//-------------
void _primr_r_chebi(const Tbl& tin, int bin, const Tbl& val0, Tbl& tout, 
        int& bout, Tbl& valp1) {

    assert(tin.dim == tout.dim) ;   

    // Output spectral basis
    int base_t = bin & MSQ_T ;
    int base_p = bin & MSQ_P ;
    bout = base_p | base_t | R_CHEBP ;

    // Number of coefficients
    int nr = tin.get_dim(0) ;	    
    int nt = tin.get_dim(1) ;	    
    int np = tin.get_dim(2) - 2 ;
    int borne_phi = np + 1 ; 
    if (np == 1) borne_phi = 1 ; 
    

    // Case of a zero input
    // --------------------
    if (tin.get_etat() == ETATZERO) {
        if (val0.get_etat() == ETATZERO) {
            tout.set_etat_zero() ; 
            valp1.set_etat_zero() ; 
            return ; 
        }
        else {
            assert(val0.get_etat() == ETATQCQ) ; 
            tout.annule_hard() ; 
            valp1.annule_hard() ; 
            double* xco = tout.t ;	
            for (int k=0 ; k< borne_phi ; k++) {
	        if (k==1) {     // jump over the coefficient of sin(0*phi) 
	            xco += nr*nt ;
	        }
	        else {
	            for (int j=0 ; j<nt ; j++) {
                        xco[0] = val0(k,j) ;  // constant value = boundary value
                        for (int i=1; i<nr; i++) xco[i] = 0 ; 
                        valp1.set(k,j) = xco[0]  ;
	                xco += nr ;
                    }
                }
            }
            return ; 
        }
    }

    // Case of a non-zero input
    // ------------------------

    assert(tin.get_etat() == ETATQCQ ) ; 
    tout.annule_hard() ; 
    valp1.annule_hard() ; 
    
    const double* xci = tin.t ;	
    double* xco = tout.t ;	

    for (int k=0 ; k< borne_phi ; k++) {
        if (k==1) {     // jump over the coefficient of sin(0*phi) 
            xci += nr*nt ;
            xco += nr*nt ;
        }
        else {
            for (int j=0 ; j<nt ; j++) {
            
                for (int i=1; i<nr-1; i++) {
                    xco[i] = (xci[i-1] - xci[i]) / double(4*i) ; 
                }
                
                xco[nr-1] = xci[nr-2] / double(4*nr - 4) ; 

                // Determination of the T_0 coefficient by maching with
                // provided value at xi = 0 : 
                double som = - xco[1] ; 
                for (int i=2; i<nr; i+=2) som += xco[i] ; 
                for (int i=3; i<nr; i+=2) som -= xco[i] ; 
                xco[0] = val0(k,j) - som ;                 

                // Value of primitive at xi = + 1 : 
                som = xco[0] ; 
                for (int i=1; i<nr; i++) som += xco[i] ;
                valp1.set(k,j) = som ; 

                xci += nr ;
                xco += nr ;
            }   // end of theta loop
        }   
    }   // end of phi loop

}



// case R_CHEBPIM_P
//-----------------
void _primr_r_chebpim_p(const Tbl& tin, int bin, const Tbl& val0, Tbl& tout, 
        int& bout, Tbl& valp1) {

    assert(tin.dim == tout.dim) ;   

    // Output spectral basis
    int base_t = bin & MSQ_T ;
    int base_p = bin & MSQ_P ;
    bout = base_p | base_t | R_CHEBPIM_I ;

    // Number of coefficients
    int nr = tin.get_dim(0) ;	    
    int nt = tin.get_dim(1) ;	    
    int np = tin.get_dim(2) - 2 ;
    

    // Case of a zero input
    // --------------------
    if (tin.get_etat() == ETATZERO) {
        if (val0.get_etat() == ETATZERO) {
            tout.set_etat_zero() ; 
            valp1.set_etat_zero() ; 
            return ; 
        }
        else {
            assert(val0.get_etat() == ETATQCQ) ; 
            tout.annule_hard() ; 
            valp1.annule_hard() ; 
            double* xco = tout.t ;	

            // m even part 
            for (int k=0 ; k<np+1 ; k += 4) {
        	int auxiliaire = (k==np) ? 1 : 2 ; // to avoid the last coef
	        for (int kmod=0 ; kmod<auxiliaire ; kmod++) {
	            if ((k==0) && (kmod == 1)) {     // jump over the coefficient of sin(0*phi) 
	                xco += nr*nt ;
	            }
	            else {
	                for (int j=0 ; j<nt ; j++) {
                            assert( val0(k+kmod,j) == double(0) ) ; 
                            for (int i=0; i<nr; i++) xco[i] = 0 ; 
                            valp1.set(k+kmod,j) = 0. ;
	                    xco += nr ;
                        }
                    }
                }
                xco += 2*nr*nt ;    // next even m
            }
            
            // m odd part
            xco = tout.t + 2*nr*nt ; 
            for (int k=2 ; k<np+1 ; k += 4) {
	        int auxiliaire = (k==np) ? 1 : 2 ; // to avoid the last coef
	        for (int kmod=0 ; kmod<auxiliaire ; kmod++) {
	            for (int j=0 ; j<nt ; j++) {
                        xco[0] = val0(k+kmod,j) ;  // constant value = boundary value
                        for (int i=1; i<nr; i++) xco[i] = 0 ; 
                        valp1.set(k+kmod,j) = xco[0] ;
	                xco += nr ;
                    }
                }
                xco += 2*nr*nt ;    // next odd m
            }
            return ; 
        }
    }

    // Case of a non-zero input
    // ------------------------

    assert(tin.get_etat() == ETATQCQ ) ; 
    tout.annule_hard() ; 
    valp1.annule_hard() ; 
    
    const double* xci = tin.t ;	
    double* xco = tout.t ;	

    // m even part 
    // -----------
    for (int k=0 ; k<np+1 ; k += 4) {
        int auxiliaire = (k==np) ? 1 : 2 ; // to avoid the last coef
        for (int kmod=0 ; kmod<auxiliaire ; kmod++) {
            if ((k==0) && (kmod == 1)) {     // jump over the coefficient of sin(0*phi) 
                xci += nr*nt ;
                xco += nr*nt ;
            }
            else {
                for (int j=0 ; j<nt ; j++) {
                    xco[0] = xci[0] - 0.5*xci[1] ; // special case i = 0

                    for (int i=1; i<nr-2; i++) {
                        xco[i] = (xci[i] - xci[i+1]) / double(4*i+2) ; 
                    }
                
                    xco[nr-2] = xci[nr-2] / double(4*nr - 6) ; 
                    xco[nr-1] = 0 ; 

                    // Value of primitive at xi = + 1 : 
                    double som = xco[0] ; 
                    for (int i=1; i<nr; i++) som += xco[i] ;
                    valp1.set(k+kmod,j) = som ; 
                
                    xci += nr ;
                    xco += nr ;
                }   // end of theta loop
                    
            }
        }
        xci += 2*nr*nt ;    // next even m
        xco += 2*nr*nt ;    // 
    }

    // m odd part 
    // ----------
    xci = tin.t + 2*nr*nt ;
    xco = tout.t + 2*nr*nt ;
    for (int k=2 ; k<np+1 ; k += 4) {
        int auxiliaire = (k==np) ? 1 : 2 ; // to avoid the last coef
        for (int kmod=0 ; kmod<auxiliaire ; kmod++) {
            for (int j=0 ; j<nt ; j++) {
            
                for (int i=1; i<nr-1; i++) {
                    xco[i] = (xci[i-1] - xci[i]) / double(4*i) ; 
                }
                
                xco[nr-1] = xci[nr-2] / double(4*nr - 4) ; 

                // Determination of the T_0 coefficient by maching with
                // provided value at xi = 0 : 
                double som = - xco[1] ; 
                for (int i=2; i<nr; i+=2) som += xco[i] ; 
                for (int i=3; i<nr; i+=2) som -= xco[i] ; 
                xco[0] = val0(k+kmod,j) - som ;                 

                // Value of primitive at xi = + 1 : 
                som = xco[0] ; 
                for (int i=1; i<nr; i++) som += xco[i] ;
                valp1.set(k+kmod,j) = som ; 
                
                xci += nr ;
                xco += nr ;
            }   // end of theta loop
        }
        xci += 2*nr*nt ;    // next odd m
        xco += 2*nr*nt ;    // 
    }


}



// case R_CHEBPIM_I
//-----------------
void _primr_r_chebpim_i(const Tbl& tin, int bin, const Tbl& val0, Tbl& tout, 
        int& bout, Tbl& valp1) {

    assert(tin.dim == tout.dim) ;   

    // Output spectral basis
    int base_t = bin & MSQ_T ;
    int base_p = bin & MSQ_P ;
    bout = base_p | base_t | R_CHEBPIM_P ;

    // Number of coefficients
    int nr = tin.get_dim(0) ;	    
    int nt = tin.get_dim(1) ;	    
    int np = tin.get_dim(2) - 2 ;
    

    // Case of a zero input
    // --------------------
    if (tin.get_etat() == ETATZERO) {
        if (val0.get_etat() == ETATZERO) {
            tout.set_etat_zero() ; 
            valp1.set_etat_zero() ; 
            return ; 
        }
        else {
            assert(val0.get_etat() == ETATQCQ) ; 
            tout.annule_hard() ; 
            valp1.annule_hard() ; 
            double* xco = tout.t ;	

            // m odd part 
            for (int k=0 ; k<np+1 ; k += 4) {
        	int auxiliaire = (k==np) ? 1 : 2 ; // to avoid the last coef
	        for (int kmod=0 ; kmod<auxiliaire ; kmod++) {
	            if ((k==0) && (kmod == 1)) {     // jump over the coefficient of sin(0*phi) 
	                xco += nr*nt ;
	            }
	            else {
	                for (int j=0 ; j<nt ; j++) {
                            xco[0] = val0(k+kmod,j) ;  // constant value = boundary value
                            for (int i=1; i<nr; i++) xco[i] = 0 ; 
                            valp1.set(k+kmod,j) = xco[0] ;
	                    xco += nr ;
                        }
                    }
                }
                xco += 2*nr*nt ;    // next odd m
            }
            
            // m even part
            xco = tout.t + 2*nr*nt ; 
            for (int k=2 ; k<np+1 ; k += 4) {
	        int auxiliaire = (k==np) ? 1 : 2 ; // to avoid the last coef
	        for (int kmod=0 ; kmod<auxiliaire ; kmod++) {
	                for (int j=0 ; j<nt ; j++) {
                            assert( val0(k+kmod,j) == double(0) ) ; 
                            for (int i=0; i<nr; i++) xco[i] = 0 ; 
                            valp1.set(k+kmod,j) = 0. ;
	                    xco += nr ;
                        }
                }
                xco += 2*nr*nt ;    // next even m
            }
            return ; 
        }
    }

    // Case of a non-zero input
    // ------------------------

    assert(tin.get_etat() == ETATQCQ ) ; 
    tout.annule_hard() ; 
    valp1.annule_hard() ; 
    
    const double* xci = tin.t ;	
    double* xco = tout.t ;	

    // m odd part 
    // ----------
    for (int k=0 ; k<np+1 ; k += 4) {
        int auxiliaire = (k==np) ? 1 : 2 ; // to avoid the last coef
        for (int kmod=0 ; kmod<auxiliaire ; kmod++) {
            if ((k==0) && (kmod == 1)) {     // jump over the coefficient of sin(0*phi) 
                xci += nr*nt ;
                xco += nr*nt ;
            }
            else {
                for (int j=0 ; j<nt ; j++) {
            
                    for (int i=1; i<nr-1; i++) {
                        xco[i] = (xci[i-1] - xci[i]) / double(4*i) ; 
                    }
                
                    xco[nr-1] = xci[nr-2] / double(4*nr - 4) ; 

                    // Determination of the T_0 coefficient by maching with
                    // provided value at xi = 0 : 
                    double som = - xco[1] ; 
                    for (int i=2; i<nr; i+=2) som += xco[i] ; 
                    for (int i=3; i<nr; i+=2) som -= xco[i] ; 
                    xco[0] = val0(k+kmod,j) - som ;                 

                    // Value of primitive at xi = + 1 : 
                    som = xco[0] ; 
                    for (int i=1; i<nr; i++) som += xco[i] ;
                    valp1.set(k+kmod,j) = som ; 
                
                    xci += nr ;
                    xco += nr ;
                }   // end of theta loop                    
            }
        }
        xci += 2*nr*nt ;    // next odd m
        xco += 2*nr*nt ;    // 
    }

    // m even part 
    // -----------
    xci = tin.t + 2*nr*nt ;
    xco = tout.t + 2*nr*nt ;
    for (int k=2 ; k<np+1 ; k += 4) {
        int auxiliaire = (k==np) ? 1 : 2 ; // to avoid the last coef
        for (int kmod=0 ; kmod<auxiliaire ; kmod++) {

            for (int j=0 ; j<nt ; j++) {
                xco[0] = xci[0] - 0.5*xci[1] ; // special case i = 0

                for (int i=1; i<nr-2; i++) {
                    xco[i] = (xci[i] - xci[i+1]) / double(4*i+2) ; 
                }
                
                xco[nr-2] = xci[nr-2] / double(4*nr - 6) ; 
                xco[nr-1] = 0 ; 

                // Value of primitive at xi = + 1 : 
                double som = xco[0] ; 
                for (int i=1; i<nr; i++) som += xco[i] ;
                valp1.set(k+kmod,j) = som ; 
                
                xci += nr ;
                xco += nr ;
            }   // end of theta loop

        }
        xci += 2*nr*nt ;    // next even m
        xco += 2*nr*nt ;    // 
    }


}

// case R_CHEBPI_P
//-------------
void _primr_r_chebpi_p(const Tbl& tin, int bin, const Tbl& val0, Tbl& tout, 
        int& bout, Tbl& valp1) {

    assert(tin.dim == tout.dim) ;   

    // Output spectral basis
    int base_t = bin & MSQ_T ;
    int base_p = bin & MSQ_P ;
    bout = base_p | base_t | R_CHEBPI_I ;

    // Number of coefficients
    int nr = tin.get_dim(0) ;	    
    int nt = tin.get_dim(1) ;	    
    int np = tin.get_dim(2) - 2 ;
    int borne_phi = np + 1 ; 
    if (np == 1) borne_phi = 1 ; 
      

     // Case of a zero input
    // --------------------
    if (tin.get_etat() == ETATZERO) {
        if (val0.get_etat() == ETATZERO) {
            tout.set_etat_zero() ; 
            valp1.set_etat_zero() ; 
            return ; 
        }
        else {
            assert(val0.get_etat() == ETATQCQ) ; 
            tout.annule_hard() ; 
            valp1.annule_hard() ; 
            double* xco = tout.t ;	
            for (int k=0 ; k< borne_phi ; k++) {
	        if (k==1) {     // jump over the coefficient of sin(0*phi) 
	            xco += nr*nt ;
	        }
	        else {
	            for (int j=0 ; j<nt ; j++) {
		      int l = j%2;
		      if(l==0){
			for (int i=0; i<nr; i++) xco[i] = 0 ; 
			valp1.set(k,j) = 0. ;
		      } else {
                        xco[0] = val0(k,j) ;  // constant value = boundary value
                        for (int i=1; i<nr; i++) xco[i] = 0 ; 
                        valp1.set(k,j) = xco[0]  ;
		      }
		      xco += nr ;
                    }
                }
            }
            return ; 
        }
    }
   
    // Case of a non-zero input
    // ------------------------

    assert(tin.get_etat() == ETATQCQ ) ; 
    tout.annule_hard() ; 
    valp1.annule_hard() ; 
    
    const double* xci = tin.t ;	
    double* xco = tout.t ;	

    for (int k=0 ; k< borne_phi ; k++) {
        if (k==1) {     // jump over the coefficient of sin(0*phi) 
            xci += nr*nt ;
            xco += nr*nt ;
        }
        else {
            for (int j=0 ; j<nt ; j++) {

	        int l = j%2;
		if(l==0){
		  xco[0] = xci[0] - 0.5*xci[1] ; // special case i = 0
		  
		  for (int i=1; i<nr-2; i++) {
                    xco[i] = (xci[i] - xci[i+1]) / double(4*i+2) ; 
		  }
		  
		  xco[nr-2] = xci[nr-2] / double(4*nr - 6) ; 
		  xco[nr-1] = 0 ; 
		  
		  // Value of primitive at xi = + 1 : 
		  double som = xco[0] ; 
		  for (int i=1; i<nr; i++) som += xco[i] ;
		  valp1.set(k,j) = som ;
		} else {
		  for (int i=1; i<nr-1; i++) {
                    xco[i] = (xci[i-1] - xci[i]) / double(4*i) ; 
		  }
		  
		  xco[nr-1] = xci[nr-2] / double(4*nr - 4) ; 
		  
		  // Determination of the T_0 coefficient by maching with
		  // provided value at xi = 0 : 
		  double som = - xco[1] ; 
		  for (int i=2; i<nr; i+=2) som += xco[i] ; 
		  for (int i=3; i<nr; i+=2) som -= xco[i] ; 
		  xco[0] = val0(k,j) - som ;                 
		  
		  // Value of primitive at xi = + 1 : 
		  som = xco[0] ; 
		  for (int i=1; i<nr; i++) som += xco[i] ;
		  valp1.set(k,j) = som ; 
                }
                xci += nr ;
                xco += nr ;
            }   // end of theta loop
        }   
    }   // end of phi loop

}
// case R_CHEBPI_I
//-------------
void _primr_r_chebpi_i(const Tbl& tin, int bin, const Tbl& val0, Tbl& tout, 
        int& bout, Tbl& valp1) {

    assert(tin.dim == tout.dim) ;   

    // Output spectral basis
    int base_t = bin & MSQ_T ;
    int base_p = bin & MSQ_P ;
    bout = base_p | base_t | R_CHEBPI_P ;

    // Number of coefficients
    int nr = tin.get_dim(0) ;	    
    int nt = tin.get_dim(1) ;	    
    int np = tin.get_dim(2) - 2 ;
    int borne_phi = np + 1 ; 
    if (np == 1) borne_phi = 1 ; 
      

     // Case of a zero input
    // --------------------
    if (tin.get_etat() == ETATZERO) {
        if (val0.get_etat() == ETATZERO) {
            tout.set_etat_zero() ; 
            valp1.set_etat_zero() ; 
            return ; 
        }
        else {
            assert(val0.get_etat() == ETATQCQ) ; 
            tout.annule_hard() ; 
            valp1.annule_hard() ; 
            double* xco = tout.t ;	
            for (int k=0 ; k< borne_phi ; k++) {
	        if (k==1) {     // jump over the coefficient of sin(0*phi) 
	            xco += nr*nt ;
	        }
	        else {
	            for (int j=0 ; j<nt ; j++) {
		      int l = j%2;
		      if(l==1){
			for (int i=0; i<nr; i++) xco[i] = 0 ; 
			valp1.set(k,j) = 0. ;
		      } else {
                        xco[0] = val0(k,j) ;  // constant value = boundary value
                        for (int i=1; i<nr; i++) xco[i] = 0 ; 
                        valp1.set(k,j) = xco[0]  ;
		      }
		      xco += nr ;
                    }
                }
            }
            return ; 
        }
    }
   
    // Case of a non-zero input
    // ------------------------

    assert(tin.get_etat() == ETATQCQ ) ; 
    tout.annule_hard() ; 
    valp1.annule_hard() ; 
    
    const double* xci = tin.t ;	
    double* xco = tout.t ;	

    for (int k=0 ; k< borne_phi ; k++) {
        if (k==1) {     // jump over the coefficient of sin(0*phi) 
            xci += nr*nt ;
            xco += nr*nt ;
        }
        else {
            for (int j=0 ; j<nt ; j++) {

	        int l = j%2;
		if(l==1){
		  xco[0] = xci[0] - 0.5*xci[1] ; // special case i = 0
		  
		  for (int i=1; i<nr-2; i++) {
                    xco[i] = (xci[i] - xci[i+1]) / double(4*i+2) ; 
		  }
		  
		  xco[nr-2] = xci[nr-2] / double(4*nr - 6) ; 
		  xco[nr-1] = 0 ; 
		  
		  // Value of primitive at xi = + 1 : 
		  double som = xco[0] ; 
		  for (int i=1; i<nr; i++) som += xco[i] ;
		  valp1.set(k,j) = som ;
		} else {
		  for (int i=1; i<nr-1; i++) {
                    xco[i] = (xci[i-1] - xci[i]) / double(4*i) ; 
		  }
		  
		  xco[nr-1] = xci[nr-2] / double(4*nr - 4) ; 
		  
		  // Determination of the T_0 coefficient by maching with
		  // provided value at xi = 0 : 
		  double som = - xco[1] ; 
		  for (int i=2; i<nr; i+=2) som += xco[i] ; 
		  for (int i=3; i<nr; i+=2) som -= xco[i] ; 
		  xco[0] = val0(k,j) - som ;                 
		  
		  // Value of primitive at xi = + 1 : 
		  som = xco[0] ; 
		  for (int i=1; i<nr; i++) som += xco[i] ;
		  valp1.set(k,j) = som ; 
                }
                xci += nr ;
                xco += nr ;
            }   // end of theta loop
        }   
    }   // end of phi loop

}





















 
