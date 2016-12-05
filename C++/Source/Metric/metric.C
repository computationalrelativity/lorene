/*
 *  Definition of methods for the class Metric.
 *
 */

/*
 *   Copyright (c) 2003 Eric Gourgoulhon & Jerome Novak
 *
 *   Copyright (c) 1999-2001 Philippe Grandclement (for previous class Metrique)
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

 

/*
 * $Id$
 * $Log$
 * Revision 1.14  2016/12/05 16:17:59  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.13  2014/10/13 08:53:07  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.12  2014/10/06 15:13:14  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.11  2005/03/02 15:03:46  f_limousin
 * p_radial_vect is added in del_deriv() and set_der_0x0.
 *
 * Revision 1.10  2004/11/18 12:22:33  jl_jaramillo
 * Method to compute the unit radial vector field with respect
 * spherical surfaces
 *
 * Revision 1.9  2004/02/18 18:45:36  e_gourgoulhon
 * Computation of p_ricci_scal thanks to the new method
 * Tensor::trace(const Metric& ).
 *
 * Revision 1.8  2004/01/22 14:35:23  e_gourgoulhon
 * Corrected bug in operator=(const Sym_tensor& ): adresses of deleted
 * p_met_cov and p_met_con are now set to 0x0.
 *
 * Revision 1.7  2004/01/20 09:51:40  f_limousin
 * Correction in metric::determinant() : indices of tensors go now from 1 to
 * 3 !
 *
 * Revision 1.6  2003/12/30 23:06:30  e_gourgoulhon
 * Important reorganization of class Metric:
 *   -- suppression of virtual methods fait_* : the actual computations
 *      are now performed via the virtual methods con(), cov(), connect(),
 *      ricci(), ricci_scal(), determinant()
 *   -- the member p_connect is now treated as an ordinary derived data
 *      member
 *   -- the construction of the associated connection (member p_connect)
 *      is performed thanks to the new methods Map::flat_met_spher() and
 *      Map::flat_met_cart().
 *
 * Revision 1.5  2003/10/28 21:23:59  e_gourgoulhon
 * Method Tensor::contract(int, int) renamed Tensor::scontract(int, int).
 *
 * Revision 1.4  2003/10/06 16:17:30  j_novak
 * Calculation of contravariant derivative and Ricci scalar.
 *
 * Revision 1.3  2003/10/06 13:58:47  j_novak
 * The memory management has been improved.
 * Implementation of the covariant derivative with respect to the exact Tensor
 * type.
 *
 * Revision 1.2  2003/10/03 11:21:47  j_novak
 * More methods for the class Metric
 *
 * Revision 1.1  2003/10/02 15:45:50  j_novak
 * New class Metric
 *
 *
 *
 * $Header$
 *
 */

// C headers
#include <cstdlib>

// Lorene headers
#include "metric.h"
#include "utilitaires.h"

                    //-----------------//
                    //  Constructors   //
                    //-----------------//

namespace Lorene {
Metric::Metric(const Sym_tensor& symti) : mp(&symti.get_mp()),
					  p_met_cov(0x0), 
                                          p_met_con(0x0) {
  
  int type_index = symti.get_index_type(0) ;
  assert (symti.get_index_type(1) == type_index) ;

  if (type_index == COV) {
    p_met_cov = new Sym_tensor(symti) ;
  }
  else {
    assert(type_index == CON) ;
    p_met_con = new Sym_tensor(symti) ;
  }

  set_der_0x0() ;
  set_tensor_depend_0x0() ;

}


Metric::Metric(const Metric& meti) : mp(meti.mp),
                                     p_met_cov(0x0),
				     p_met_con(0x0) {

  if (meti.p_met_cov != 0x0) p_met_cov = new Sym_tensor(*meti.p_met_cov) ;

  if (meti.p_met_con != 0x0) p_met_con = new Sym_tensor(*meti.p_met_con) ;

  set_der_0x0() ;
  set_tensor_depend_0x0() ;

}

Metric::Metric(const Map& mpi, FILE* ) : mp(&mpi), 
					 p_met_cov(0x0), 
                                         p_met_con(0x0) {

  cout << "Metric::Metric(FILE*) : not implemented yet!" << endl ;

  abort() ;
}

Metric::Metric(const Map& mpi) : mp(&mpi),  
			         p_met_cov(0x0), 
                                 p_met_con(0x0) {
  set_der_0x0() ;
  set_tensor_depend_0x0() ;

}


                    //---------------//
                    //  Destructor   //
                    //---------------//

Metric::~Metric() {

  if (p_met_cov != 0x0) delete p_met_cov ;

  if (p_met_con != 0x0) delete p_met_con ;

  del_deriv() ;

  del_tensor_depend() ;

}

                //-------------------//
                // Memory management //
                //-------------------//
                
void Metric::del_deriv() const {
  
  if (p_connect != 0x0) delete p_connect ; 
  if (p_ricci_scal != 0x0) delete p_ricci_scal ;
  if (p_determinant != 0x0) delete p_determinant ;
  if (p_radial_vect != 0x0) delete p_radial_vect ;
  
  set_der_0x0() ;

    //## call to del_tensor_depend() ?
}

void Metric::set_der_0x0() const {

  p_connect = 0x0 ; 
  p_ricci_scal = 0x0 ;
  p_determinant = 0x0 ;
  p_radial_vect = 0x0 ;

}

void Metric::del_tensor_depend() const {

    for (int i=0 ; i<N_TENSOR_DEPEND ; i++)
	if (tensor_depend[i] != 0x0) {
	  int j = tensor_depend[i]->get_place_met(*this) ;
	  if (j!=-1) tensor_depend[i]->del_derive_met(j) ;
	}
    set_tensor_depend_0x0() ;
 
}

void Metric::set_tensor_depend_0x0() const {

  for (int i=0 ; i<N_TENSOR_DEPEND ; i++) {
    tensor_depend[i] = 0x0 ;
  }
}

  
                    //-----------------------//
                    // Mutators / assignment //
                    //-----------------------//

void Metric::operator=(const Metric& meti) {

  assert( mp == meti.mp) ;

  if (p_met_cov != 0x0) {
    delete p_met_cov ;
    p_met_cov = 0x0 ; 
  }
  
  if (p_met_con != 0x0) {
    delete p_met_con ;
    p_met_con = 0x0 ; 
  }
  
  if (meti.p_met_cov != 0x0) {
    p_met_cov = new Sym_tensor(*meti.p_met_cov) ;
  }

  if (meti.p_met_con != 0x0) {
    p_met_con = new Sym_tensor(*meti.p_met_con) ;
  }

  del_deriv() ;

}

void Metric::operator=(const Sym_tensor& symti) {
  
  assert(mp == &symti.get_mp()) ;

  int type_index = symti.get_index_type(0) ;
  assert (symti.get_index_type(1) == type_index) ;

  if (p_met_cov != 0x0) {
    delete p_met_cov ;
    p_met_cov = 0x0 ; 
  }
  
  if (p_met_con != 0x0) {
    delete p_met_con ;
    p_met_con = 0x0 ; 
  }
  
  if (type_index == COV) {
    p_met_cov = new Sym_tensor(symti) ;
  }
  else {
    assert(type_index == CON) ;
    p_met_con = new Sym_tensor(symti) ;
  }

  del_deriv() ;

}

            //----------------//
            //   Accessors    //
            //----------------//


const Sym_tensor& Metric::cov() const {
  
    if (p_met_cov == 0x0) {   // a new computation is necessary
        assert( p_met_con != 0x0 ) ;
        p_met_cov = p_met_con->inverse() ;
    }

    return *p_met_cov ; 
}

const Sym_tensor& Metric::con() const {
  
    if (p_met_con == 0x0) {   // a new computation is necessary
        assert( p_met_cov != 0x0 ) ;
        p_met_con = p_met_cov->inverse() ;
    }

    return *p_met_con ; 
}


const Connection& Metric::connect() const {

    if (p_connect == 0x0) {   // a new computation is necessary
    
        // The triad is obtained from the covariant or contravariant representation:
        const Base_vect_spher* triad_s ; 
        const Base_vect_cart* triad_c ; 
        if (p_met_cov != 0x0) {
            triad_s = 
              dynamic_cast<const Base_vect_spher*>(p_met_cov->get_triad()) ; 
            triad_c = 
              dynamic_cast<const Base_vect_cart*>(p_met_cov->get_triad()) ; 
        }
        else {
            assert(p_met_con != 0x0) ; 
            triad_s = 
              dynamic_cast<const Base_vect_spher*>(p_met_con->get_triad()) ; 
            triad_c = 
              dynamic_cast<const Base_vect_cart*>(p_met_con->get_triad()) ; 
        }
    
        // Background flat metric in spherical or Cartesian components
        if ( triad_s != 0x0 ) {
            p_connect = new Connection(*this, mp->flat_met_spher()) ;
        }
        else {
            assert( triad_c != 0x0 ) ;
            p_connect = new Connection(*this, mp->flat_met_cart()) ;
        }
    
    }

    return *p_connect ; 

}


const Sym_tensor& Metric::ricci() const {

    const Tensor& ricci_connect = connect().ricci() ; 
    
    // Check: the Ricci tensor of the connection associated with 
    //  the metric must be symmetric:
    assert( typeid(ricci_connect) == typeid(const Sym_tensor&) ) ; 

    return dynamic_cast<const Sym_tensor&>( ricci_connect ) ; 
}


const Scalar& Metric::ricci_scal() const {

    if (p_ricci_scal == 0x0) {   // a new computation is necessary

        p_ricci_scal = new Scalar( ricci().trace(*this) ) ; 

    }

    return *p_ricci_scal  ; 

}

const Vector& Metric::radial_vect() const {

  if (p_radial_vect == 0x0) { // a new computation is necessary


    p_radial_vect = new Vector ((*this).get_mp(), CON, *((*this).con().get_triad()) ) ;

    Scalar prov ( sqrt((*this).con()(1,1)) ) ;

    prov.std_spectral_base() ;

       
    p_radial_vect->set(1) = (*this).con()(1,1)/ prov ;
 
    p_radial_vect->set(2) = (*this).con()(1,2)/ prov ;

    p_radial_vect->set(3) = (*this).con()(1,3)/ prov ;



    //    p_radial_vect.std_spectral_base() ;


  }

  return *p_radial_vect   ;

}


const Scalar& Metric::determinant() const {

    if (p_determinant == 0x0) {   // a new computation is necessary

        p_determinant = new Scalar(*mp) ;
        *p_determinant = cov()(1, 1)*cov()(2, 2)*cov()(3, 3) 
	    + cov()(1, 2)*cov()(2, 3)*cov()(3, 1)
	    + cov()(1, 3)*cov()(2, 1)*cov()(3, 2) 
	    - cov()(3, 1)*cov()(2, 2)*cov()(1, 3)
	    - cov()(3, 2)*cov()(2, 3)*cov()(1, 1) 
	    - cov()(3, 3)*cov()(2, 1)*cov()(1, 2) ;
    }

    return *p_determinant ; 
}


 
                //---------//
                // Outputs //
                //---------//

void Metric::sauve(FILE* fd) const {

  // Which representation is to be saved
  int indic ;
  if (p_met_cov != 0x0)
    indic = COV ;
  else if (p_met_con != 0x0)
    indic = CON ;
  else indic = 0 ;
  fwrite_be(&indic, sizeof(int), 1, fd) ;
  switch (indic) {
  case COV : {
    p_met_cov->sauve(fd) ;
    break ;
  }
  case CON : {
    p_met_con->sauve(fd) ;
    break ;
  }
  default : {
    break ;
  }
  } 
}

ostream& operator<<(ostream& ost, const Metric& meti) {

  meti >> ost ;
  return ost ;
}


ostream& Metric::operator>>(ostream& ost) const {

  ost << '\n' ;

  ost << "General type metric" << '\n' ;
  ost << "-------------------" << '\n' ;
  ost << '\n' ;

  if (p_met_cov == 0x0) {
    ost << "Covariant representation unknown!" << '\n' ;
    assert( p_met_con != 0x0) ;
    ost << "CONTRA-variant representation: " << '\n' ;
    ost << *p_met_con ;
  }
  else {
    ost << "Covariant representation: " << '\n' ;
    ost << *p_met_cov ;
  }


  if (p_connect == 0x0)
    ost << "Associated connection not computed yet." << '\n' ;
  else
    ost << "Associated connection computed." << '\n' ;

  if (p_ricci_scal == 0x0)
    ost << "Ricci scalar not computed yet." << '\n' ;
  else
    ost << "Ricci scalar computed." << '\n' ;
  
  if (p_determinant == 0x0)
    ost << "determinant not computed yet." << '\n' ;
  else
    ost << "determinant computed." << '\n' ;

  ost << endl ;
  return ost ;
}

}
