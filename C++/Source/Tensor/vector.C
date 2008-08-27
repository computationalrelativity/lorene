/*
 *  Methods of class Vector
 *
 *   (see file vector.h for documentation)
 *
 */

/*
 *   Copyright (c) 2003 Eric Gourgoulhon & Jerome Novak
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


char vector_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.27  2008/08/27 08:52:23  jl_cornou
 * Added fonctions for angular potential A
 *
 * Revision 1.26  2007/12/21 16:07:08  j_novak
 * Methods to filter Tensor, Vector and Sym_tensor objects.
 *
 * Revision 1.25  2005/02/14 13:01:50  j_novak
 * p_eta and p_mu are members of the class Vector. Most of associated functions
 * have been moved from the class Vector_divfree to the class Vector.
 *
 * Revision 1.24  2005/01/25 15:37:35  j_novak
 * Solved some dzpuis problem...
 *
 * Revision 1.23  2005/01/12 16:48:23  j_novak
 * Better treatment of the case where all vector components are null in
 * decompose_div .
 *
 * Revision 1.22  2004/10/12 09:58:25  j_novak
 * Better memory management.
 *
 * Revision 1.21  2004/10/11 09:46:31  j_novak
 * Speed improvements.
 *
 * Revision 1.20  2004/05/09 20:55:05  e_gourgoulhon
 * Added method flux.
 *
 * Revision 1.19  2004/03/29 11:57:45  e_gourgoulhon
 * Added methods ope_killing and ope_killing_conf.
 *
 * Revision 1.18  2004/02/26 22:48:50  e_gourgoulhon
 * -- Method divergence: call to Tensor::divergence and cast of the
 *    result.
 * -- Added method derive_lie.
 *
 * Revision 1.17  2004/02/24 09:46:20  j_novak
 * Correction to cope with SGI compiler's warnings.
 *
 * Revision 1.16  2004/02/20 10:53:04  j_novak
 * Added the matching of the potential adapted to the case where the
 * vector is the source of a Poisson equation (dzpuis = 4).
 *
 * Revision 1.15  2004/01/30 10:30:17  j_novak
 * Changed dzpuis handling in Vector::decompose_div (this may be temporary).
 *
 * Revision 1.14  2003/12/30 23:09:47  e_gourgoulhon
 * Change in methods derive_cov() and divergence() to take into account
 *  the change of name: Metric::get_connect() --> Metric::connect().
 *
 * Revision 1.13  2003/12/19 15:18:16  j_novak
 * Shadow variables hunt
 *
 * Revision 1.12  2003/10/29 11:04:34  e_gourgoulhon
 * dec2_dpzuis() replaced by dec_dzpuis(2).
 * inc2_dpzuis() replaced by inc_dzpuis(2).
 *
 * Revision 1.11  2003/10/22 14:24:19  j_novak
 * *** empty log message ***
 *
 * Revision 1.9  2003/10/20 13:00:38  j_novak
 * Memory error corrected
 *
 * Revision 1.8  2003/10/20 09:32:12  j_novak
 * Members p_potential and p_div_free of the Helmholtz decomposition
 * + the method decompose_div(Metric).
 *
 * Revision 1.7  2003/10/16 14:21:37  j_novak
 * The calculation of the divergence of a Tensor is now possible.
 *
 * Revision 1.6  2003/10/13 13:52:40  j_novak
 * Better managment of derived quantities.
 *
 * Revision 1.5  2003/10/06 13:58:48  j_novak
 * The memory management has been improved.
 * Implementation of the covariant derivative with respect to the exact Tensor
 * type.
 *
 * Revision 1.4  2003/10/05 21:14:20  e_gourgoulhon
 * Added method std_spectral_base().
 *
 * Revision 1.3  2003/10/03 14:10:32  e_gourgoulhon
 * Added constructor from Tensor.
 *
 * Revision 1.2  2003/10/03 14:08:46  j_novak
 * Removed old change_trid...
 *
 * Revision 1.1  2003/09/26 08:05:31  j_novak
 * New class Vector.
 *
 *
 * $Header$
 *
 */

// Headers C
#include <stdlib.h>
#include <assert.h>
#include <math.h>

// Headers Lorene
#include "metric.h"
#include "proto.h"
#include "matrice.h"
#include "nbr_spx.h"


			//--------------//
			// Constructors //
			//--------------//

// Standard constructor 
// --------------------
Vector::Vector(const Map& map, int tipe, const Base_vect& triad_i) 
		: Tensor(map, 1, tipe, triad_i) {
		
  set_der_0x0() ;

}

// Standard constructor with the triad passed as a pointer
// -------------------------------------------------------
Vector::Vector(const Map& map, int tipe, const Base_vect* triad_i) 
		: Tensor(map, 1, tipe, *triad_i) {
		
  set_der_0x0() ;
}
	
// Copy constructor
// ----------------
Vector::Vector (const Vector& source) : 
    Tensor(source) {
  
  assert(valence == 1) ;
  set_der_0x0() ;

}   


// Constructor from a {\tt Tensor}.
//--------------------------------
Vector::Vector(const Tensor& uu) : Tensor(uu) {

  assert(valence == 1) ;
  set_der_0x0() ;

}


// Constructor from a file
// -----------------------
Vector::Vector(const Map& mapping, const Base_vect& triad_i, FILE* fd) : 
  Tensor(mapping, triad_i, fd) {
   
  assert ( (valence == 1) && (n_comp == 3) ) ;
  set_der_0x0() ;

}


			//--------------//
			//  Destructor  //
			//--------------//


Vector::~Vector () {

  Vector::del_deriv() ;

}


			//-------------------//
			// Memory managment  //
			//-------------------//

void Vector::del_deriv() const {

  for (int i=0; i<N_MET_MAX; i++) 
    del_derive_met(i) ;
  
  if (p_A != 0x0) delete p_A ;
  if (p_eta != 0x0) delete p_eta ; 
  if (p_mu != 0x0) delete p_mu ; 
  set_der_0x0() ;
  Tensor::del_deriv() ;

}

void Vector::set_der_0x0() const {

  for (int i=0; i<N_MET_MAX; i++) 
    set_der_met_0x0(i) ;
  p_A = 0x0 ;
  p_eta = 0x0 ; 
  p_mu = 0x0 ; 

}

void Vector::del_derive_met(int j) const {

  assert( (j>=0) && (j<N_MET_MAX) ) ;

  if (met_depend[j] != 0x0) {
    if (p_potential[j] != 0x0)
      delete p_potential[j] ;
    if (p_div_free[j] != 0x0)
      delete p_div_free[j] ;
    
    set_der_met_0x0(j) ;
    
    Tensor::del_derive_met(j) ;
  }
}

void Vector::set_der_met_0x0(int i) const {

  assert( (i>=0) && (i<N_MET_MAX) ) ;

  p_potential[i] = 0x0 ;
  p_div_free[i] = 0x0 ;

}

void Vector::operator=(const Vector& t) {
    
    triad = t.triad ; 

    assert(t.type_indice(0) == type_indice(0)) ;

    for (int i=0 ; i<3 ; i++) {
      *cmp[i] = *t.cmp[i] ;
    }
    del_deriv() ;
}

void Vector::operator=(const Tensor& t) {
    
    assert (t.valence == 1) ;

    triad = t.triad ; 

    assert(t.type_indice(0) == type_indice(0)) ;

    for (int i=0 ; i<3 ; i++) {
      *cmp[i] = *t.cmp[i] ;
    }
    del_deriv() ;
}



// Affectation d'une composante :
Scalar& Vector::set(int index) {
    
  assert ( (index>=1) && (index<=3) ) ;

  del_deriv() ;
  
  return *cmp[index - 1] ;
}

const Scalar& Vector::operator()(int index) const {
    
  assert ( (index>=1) && (index<=3) ) ;

  return *cmp[index - 1] ;

}


// Sets the standard spectal bases of decomposition for each component

void Vector::std_spectral_base() {

	Base_val** bases = 0x0 ;

	if ( triad->identify() == (mp->get_bvect_cart()).identify() ) {

		// Cartesian case
		bases = mp->get_mg()->std_base_vect_cart() ;

	}
	else {
		// Spherical case
		assert( triad->identify() == (mp->get_bvect_spher()).identify()) ;
		bases = mp->get_mg()->std_base_vect_spher() ;
	}
	    
	for (int i=0 ; i<3 ; i++) {
		cmp[i]->set_spectral_base( *bases[i] ) ; 
	}
		
	for (int i=0 ; i<3 ; i++) {
		delete bases[i] ;
	}
	delete [] bases ;


}


                    //-------------------------------//
                    //    Computational methods      //
                    //-------------------------------//
                    

const Scalar& Vector::divergence(const Metric& metre) const {
  
    const Scalar* pscal = 
        dynamic_cast<const Scalar*>( &(Tensor::divergence(metre)) ) ;

    assert(pscal != 0x0) ;

    return *pscal ;
}


Vector Vector::derive_lie(const Vector& vv) const {

    Vector resu(*mp, type_indice(0), triad) ; 
    
    compute_derive_lie(vv, resu) ;
    
    return resu ; 
    
}


Sym_tensor Vector::ope_killing(const Metric& gam) const {

    Sym_tensor resu(*mp, type_indice(0), *triad) ; 

    if (type_indice(0) == CON ) {
        for (int i=1; i<=3; i++) {
            for (int j=i; j<=3; j++) {
                resu.set(i,j) = derive_con(gam)(i,j) + derive_con(gam)(j,i)  ;
            }
        }
    }
    else {
        for (int i=1; i<=3; i++) {
            for (int j=i; j<=3; j++) {
                resu.set(i,j) = derive_cov(gam)(i,j) + derive_cov(gam)(j,i)  ;
            }
        }
    }
    
    return resu ; 

} 


Sym_tensor Vector::ope_killing_conf(const Metric& gam) const {

    Sym_tensor resu(*mp, type_indice(0), *triad) ; 

    if (type_indice(0) == CON ) {
        for (int i=1; i<=3; i++) {
            for (int j=i; j<=3; j++) {
                resu.set(i,j) = derive_con(gam)(i,j) + derive_con(gam)(j,i)  
                                - 0.6666666666666666* divergence(gam) 
                                            * gam.con()(i,j) ;
            }
        }
    }
    else {
        for (int i=1; i<=3; i++) {
            for (int j=i; j<=3; j++) {
                resu.set(i,j) = derive_cov(gam)(i,j) + derive_cov(gam)(j,i)  
                                - 0.6666666666666666* derive_con(gam).trace() 
                                            * gam.cov()(i,j) ;
            }
        }
    }

    return resu ; 

} 




const Scalar& Vector::potential(const Metric& metre) const {

  set_dependance(metre) ;
  int j = get_place_met(metre) ;
  assert ((j>=0) && (j<N_MET_MAX)) ;
  if (p_potential[j] == 0x0) {
    decompose_div(metre) ;
  }

  return *p_potential[j] ;
}

const Vector_divfree& Vector::div_free(const Metric& metre) const {

  set_dependance(metre) ;
  int j = get_place_met(metre) ;
  assert ((j>=0) && (j<N_MET_MAX)) ;
  if (p_div_free[j] == 0x0) {
    decompose_div(metre) ;
  }

  return *p_div_free[j] ;
}

void Vector::decompose_div(const Metric& metre) const {

  assert( type_indice(0) == CON ) ; //Only for contravariant vectors...

  set_dependance(metre) ;
  int ind =  get_place_met(metre) ;
  assert ((ind>=0) && (ind<N_MET_MAX)) ;

  if ( (p_potential[ind] != 0x0) && (p_div_free[ind] != 0x0) ) 
    return ; // Nothing to do ...

  else {
    if (p_div_free[ind] != 0x0)
      delete p_div_free[ind] ;  
    if (p_potential[ind] != 0x0) 
      delete p_potential[ind] ;

    const Mg3d* mmg = mp->get_mg() ;
    int nz = mmg->get_nzone() ;

    int dzp = cmp[0]->get_dzpuis() ;
#ifndef NDEBUG
    bool dz_zero = cmp[0]->check_dzpuis(0) ;
    bool dz_four = cmp[0]->check_dzpuis(4) ;
#endif
    assert( dz_zero || dz_four) ;
    assert (cmp[1]->check_dzpuis(dzp)) ;
    assert (cmp[2]->check_dzpuis(dzp)) ;

    Scalar dive = divergence(metre) ;

    if (dive.get_etat() == ETATZERO) {
	p_potential[ind] = new Scalar(*mp) ;
	p_potential[ind]->set_etat_zero() ;
	p_potential[ind]->set_dzpuis(dzp) ;
    }
    else {
        //No matching of the solution at this point
	p_potential[ind] = new Scalar(dive.poisson()) ; 

	if (dzp == 4) {
	    assert (mmg->get_type_r(nz-1) == UNSURR) ;
	    // Let's now do the matching ... in the case dzpuis = 4
	    const Map_af* mapping = dynamic_cast<const Map_af*>(mp) ;
	    assert (mapping != 0x0) ;
	    Valeur& val = p_potential[ind]->set_spectral_va() ;
	    val.ylm() ;
	    Mtbl_cf& mtc = *val.c_cf ;
	    if (nz > 1) {
		int np = mmg->get_np(0) ;
		int nt = mmg->get_nt(0) ;
#ifndef NDEBUG
		for (int lz=0; lz<nz; lz++) {
		    assert (mmg->get_np(lz) == np) ;
		    assert (mmg->get_nt(lz) == nt) ;
		}
#endif
		int lmax = 0 ;
		for (int k=0; k<np+1; k++) 
		    for (int j=0; j<nt; j++) 
			if ( nullite_plm(j, nt, k, np, val.base)) {
			    int m_quant, l_quant, base_r ;
			    donne_lm (nz, 0, j, k, val.base, m_quant, 
				      l_quant, base_r) ; 
			    lmax = (l_quant > lmax ? l_quant : lmax) ;
			}
		Scalar** ri = new Scalar*[lmax+1] ;
		Scalar** rmi = new Scalar*[lmax+1] ;
		Scalar erre(*mp) ;
		erre = mp->r ;
		for (int l=0; l<=lmax; l++) {
		    ri[l] = new Scalar(*mp) ;
		    rmi[l] = new Scalar(*mp) ;
		    if (l == 0) *(ri[l]) = 1. ;
		    else *(ri[l]) = pow(erre, l) ;
		    ri[l]->annule_domain(nz - 1) ;
		    ri[l]->std_spectral_base() ; //exact base in r will be set later
		    if (l==2) *(rmi[l]) = 1 ;
		    else *(rmi[l]) = pow(erre, 2-l) ;
		    rmi[l]->annule(0,nz-2) ;
		    Scalar tmp = pow(erre, -l-1) ;
		    tmp.annule_domain(nz-1) ;
		    tmp.annule_domain(0) ;
		    *(rmi[l]) += tmp ;
		    rmi[l]->set_dzpuis(3) ;
		    rmi[l]->std_spectral_base() ;//exact base in r will be set later
		}

		for (int k=0; k<np+1; k++) {
		    for (int j=0; j<nt; j++) {
			if ( nullite_plm(j, nt, k, np, val.base)) {
			    int m_quant, l_quant, base_r, l, c ;
			    donne_lm (nz, 0, j, k, val.base, m_quant, l_quant, base_r) ; 
			    bool lzero = (l_quant == 0) || (l_quant == 1) ;
			    int nb_hom_ced  = (l_quant < 2 ? 0 : 1) ; 
	      
			    int taille = 2*(nz-1) - 1 + nb_hom_ced ;
			    Tbl deuz(taille) ;
			    deuz.set_etat_qcq() ;
			    Matrice systeme(taille,taille) ;
			    systeme.set_etat_qcq() ;
			    for (l=0; l<taille; l++) 
				for (c=0; c<taille; c++) systeme.set(l,c) = 0 ;
			    for (l=0; l<taille; l++) deuz.set(l) = 0 ;
	      
			    //---------
			    // Nucleus
			    //---------
			    assert(mmg->get_type_r(0) == RARE) ;
			    assert ((base_r == R_CHEBP)||(base_r == R_CHEBI)) ;
			    ri[l_quant]->set_spectral_va().set_base_r(0, base_r) ;
	  
			    double alpha = mapping->get_alpha()[0] ;
			    int nr = mmg->get_nr(0) ;
			    Tbl partn(nr) ;
			    partn.set_etat_qcq() ;
			    for (int i=0; i<nr; i++)
				partn.set(i) = mtc(0, k, j, i) ;
			    l=0 ; c=0 ;

			    systeme.set(l,c) += pow(alpha, double(l_quant)) ;

			    deuz.set(l) -= val1_dern_1d(0, partn, base_r) ;
			    l++ ;

			    if (!lzero) {
				systeme.set(l,c) += l_quant * pow(alpha, double(l_quant - 1)) ;
   
				deuz.set(l) -= val1_dern_1d(1, partn, base_r) / alpha ;
			    }

			    //----------
			    //  Shells
			    //----------

			    for (int lz=1; lz<nz-1; lz++) {
				alpha = mapping->get_alpha()[lz] ;
				double beta = mapping->get_beta()[lz] ;
				double xm1 = beta - alpha ;
				double xp1 = alpha + beta ;
				nr = mmg->get_nr(lz) ;
				Tbl parts(nr) ;
				parts.set_etat_qcq() ;
				for (int i=0; i<nr; i++)
				    parts.set(i) = mtc(lz, k, j, i) ;
  
				//Function at x = -1
				l-- ; 
				c = 2*lz - 1 ;
				systeme.set(l,c) -= pow(xm1, double(l_quant)) ;
				c++;
				systeme.set(l,c) -= pow(xm1, double(-l_quant-1)) ;
	      
				deuz.set(l) += valm1_dern_1d(0, parts, R_CHEB) ;
	    
				if ((lz>1) || (!lzero)) {
				    //First derivative at x=-1
				    l++ ;
				    c-- ;
				    systeme.set(l,c) -=  l_quant*pow(xm1, double(l_quant - 1)) ;
				    c++;
				    systeme.set(l,c) -= (-l_quant - 1)*
					pow(xm1, double(-l_quant - 2)) ;

				    deuz.set(l) += valm1_dern_1d(1, parts, R_CHEB) / alpha ;
				}

				//Function at x = 1
				l++ ; 
				c-- ;
				systeme.set(l,c) += pow(xp1, double(l_quant)) ;
				c++;
				systeme.set(l,c) += pow(xp1, double(-l_quant-1)) ;

				deuz.set(l) -= val1_dern_1d(0, parts, R_CHEB) ;
	    
				//First derivative at x = 1
				l++ ; 
				c-- ;
				systeme.set(l,c) +=  l_quant*pow(xp1, double(l_quant - 1)) ;
				c++;
				systeme.set(l,c) += (-l_quant - 1)*
				    pow(xp1, double(-l_quant - 2)) ;
	    
				deuz.set(l) -= val1_dern_1d(1, parts, R_CHEB) / alpha ;
	    
			    } //End of the loop on different shells

			    //-------------------------------
			    //  Compactified external domain
			    //-------------------------------
			    assert(mmg->get_type_r(nz-1) == UNSURR) ;
			    nr = mmg->get_nr(nz-1) ;
			    Tbl partc(nr) ;
			    partc.set_etat_qcq() ;
			    for (int i=0; i<nr; i++)
				partc.set(i) = mtc(nz-1, k, j, i) ;

			    alpha = mapping->get_alpha()[nz-1] ;
			    double beta =  mapping->get_beta()[nz-1] ;
			    double xm1 = beta - alpha ; // 1 / r_left
			    double part0, part1 ;
			    part0 = valm1_dern_1d(0, partc, R_CHEB) ;
			    part1 = xm1*valm1_dern_1d(1, partc, R_CHEB) / alpha ;
			    assert (p_potential[ind]->get_dzpuis() == 3) ;
		
			    //Function at x = -1
			    l--;
			    if (!lzero) {
				c++;
				systeme.set(l,c) -= pow(xm1, double(l_quant+1)) ;
			    }
			    deuz.set(l) += part0*xm1*xm1*xm1 ;

			    // First derivative at x = -1
			    l++ ;
			    if (!lzero) {
				systeme.set(l,c) -= (-l_quant - 1)*
				    pow(xm1, double(l_quant + 2)) ;
			    }
			    if ( (nz > 2) || (!lzero))
				deuz.set(l) += -xm1*xm1*xm1*xm1*(3*part0 + part1) ;

			    //--------------------------------------
			    //   Solution of the linear system
			    //--------------------------------------
	  
			    int inf = 1 + nb_hom_ced;
			    int sup = 3 - nb_hom_ced;
			    systeme.set_band(sup, inf) ;
			    systeme.set_lu() ;
			    Tbl facteur(systeme.inverse(deuz)) ;
			    ri[l_quant]->set_spectral_va().coef() ;
			    rmi[l_quant]->set_spectral_va().coef() ;

			    //Linear combination in the nucleus
			    nr = mmg->get_nr(0) ;
			    for (int i=0; i<nr; i++) 
				mtc.set(0, k, j, i) += 
	     facteur(0)*(*(ri[l_quant]->get_spectral_va().c_cf))(0, 0, 0, i) ;
	    
			    //Linear combination in the shells
			    for (int lz=1; lz<nz-1; lz++) {
				nr = mmg->get_nr(lz) ;
				for (int i=0; i<nr; i++) 
				    mtc.set(lz, k, j, i) += facteur(2*lz - 1)*
	       (*(ri[l_quant]->get_spectral_va().c_cf))(lz, 0, 0, i) ;
				for (int i=0; i<nr; i++) 
				    mtc.set(lz, k, j, i) += facteur(2*lz)*
	       (*(rmi[l_quant]->get_spectral_va().c_cf))(lz, 0, 0, i) ;
			    }

			    //Linear combination in the CED
			    nr = mmg->get_nr(nz-1) ;
			    if (!lzero) {
				for (int i=0; i<nr; i++) 
				    mtc.set(nz-1, k, j, i) += 
					facteur(taille - 1)*
	      (*(rmi[l_quant]->get_spectral_va().c_cf))(nz-1, 0, 0, i) ;
			    }
			} //End of nullite_plm ...

		    } //End of j/theta loop   
		} //End of k/phi loop 

		for (int l=0; l<=lmax; l++) {
		    delete ri[l] ;
		    delete rmi[l] ;
		}
		delete [] ri ;
		delete [] rmi ;

	    } //End of the case of more than one domain

	    val.ylm_i() ;

	} //End of the case dzp = 4
    }
    p_div_free[ind] = new Vector_divfree(*mp, *triad, metre) ;

    Vector gradient = p_potential[ind]->derive_con(metre) ;
    if (dzp != 4)  gradient.dec_dzpuis(2) ;

    *p_div_free[ind] = ( *this - gradient ) ;

  }
  
}



double Vector::flux(double radius, const Metric& met) const {

    assert(type_indice(0) == CON) ; 
    
    const Map_af* mp_af = dynamic_cast<const Map_af*>(mp) ; 
    if (mp_af == 0x0) {
        cerr << 
        "Vector::flux : the case of a mapping outside the class Map_af\n"
            << " is not implemented yet !" << endl ; 
        abort() ; 
    } 
    
    const Metric_flat* ff = dynamic_cast<const Metric_flat*>(&met) ;
    if (ff == 0x0) {
        cerr << 
        "Vector::flux : the case of a non flat metric is not implemented yet !"
             << endl ; 
        abort() ; 
    }
    
    const Base_vect_cart* bcart = dynamic_cast<const Base_vect_cart*>(triad) ; 
    Vector* vspher = 0x0 ; 
    if (bcart != 0x0) { // switch to spherical components:
        vspher = new Vector(*this) ; 
        vspher->change_triad(mp->get_bvect_spher()) ;
    }
    else assert( dynamic_cast<const Base_vect_spher*>(triad) != 0x0 ) ; 

    const Vector* vv = (bcart != 0x0) ? vspher : this ; 

    const Scalar& vr = vv->operator()(1) ; 
    
    double resu ; 
    if (radius == __infinity) {
        resu = mp_af->integrale_surface_infini(vr) ; 
    }
    else {
        resu = mp_af->integrale_surface(vr, radius) ;
    }
    
    return resu ; 
}

void Vector::exponential_filter_r(int lzmin, int lzmax, int p, 
			  double alpha) {
    if( triad->identify() == (mp->get_bvect_cart()).identify() ) 
	for (int i=0; i<n_comp; i++)
	    cmp[i]->exponential_filter_r(lzmin, lzmax, p, alpha) ;
    else {
	assert( triad->identify() == (mp->get_bvect_spher()).identify()) ;
	Scalar vr_tmp = operator()(1) ; 
	vr_tmp.exponential_filter_r(lzmin, lzmax, p, alpha) ;
	Scalar eta_tmp = eta() ;
	eta_tmp.exponential_filter_r(lzmin, lzmax, p, alpha) ;
	Scalar mu_tmp = mu() ;
	mu_tmp.exponential_filter_r(lzmin, lzmax, p, alpha) ;
	set_vr_eta_mu(vr_tmp, eta_tmp, mu_tmp) ;
    }
}

void Vector::exponential_filter_ylm(int lzmin, int lzmax, int p, 
			  double alpha) {
    if( triad->identify() == (mp->get_bvect_cart()).identify() ) 
	for (int i=0; i<n_comp; i++)
	    cmp[i]->exponential_filter_ylm(lzmin, lzmax, p, alpha) ;
    else {
	assert( triad->identify() == (mp->get_bvect_spher()).identify()) ;
	Scalar vr_tmp = operator()(1) ; 
	vr_tmp.exponential_filter_ylm(lzmin, lzmax, p, alpha) ;
	Scalar eta_tmp = eta() ;
	eta_tmp.exponential_filter_ylm(lzmin, lzmax, p, alpha) ;
	Scalar mu_tmp = mu() ;
	mu_tmp.exponential_filter_ylm(lzmin, lzmax, p, alpha) ;
	set_vr_eta_mu(vr_tmp, eta_tmp, mu_tmp) ;
    }
}
