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

char metric_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
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
#include<stdlib.h>

// Lorene headers
#include "metric.h"
#include "utilitaires.h"

Metric::Metric(const Sym_tensor& symti) : mp(&symti.get_mp()), p_connect(0x0),
					  p_met_cov(0x0), p_met_con(0x0) {
  
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


Metric::Metric(const Metric& meti) : mp(meti.mp), p_connect(0x0), p_met_cov(0x0),
				     p_met_con(0x0) {

  if (meti.p_met_cov != 0x0) p_met_cov = new Sym_tensor(*meti.p_met_cov) ;

  if (meti.p_met_con != 0x0) p_met_con = new Sym_tensor(*meti.p_met_con) ;

  set_der_0x0() ;
  set_tensor_depend_0x0() ;

}

Metric::Metric(const Map& mpi, FILE* ) : mp(&mpi), p_connect(0x0), 
					 p_met_cov(0x0), p_met_con(0x0) {

  cout << "Metric::Metric(FILE*) : not implemented yet!" << endl ;

  abort() ;
}

Metric::Metric(const Map& mpi) : mp(&mpi), p_connect(0x0), 
					 p_met_cov(0x0), p_met_con(0x0) {
  set_der_0x0() ;
  set_tensor_depend_0x0() ;

}

Metric::~Metric() {

  if (p_connect != 0x0) delete p_connect ;

  if (p_met_cov != 0x0) delete p_met_cov ;

  if (p_met_con != 0x0) delete p_met_con ;

  del_deriv() ;

  del_tensor_depend() ;

}

void Metric::del_deriv() const {
  
  if (p_ricci_scal != 0x0) delete p_ricci_scal ;

  if (p_determinant != 0x0) delete p_determinant ;

  set_der_0x0() ;

}

void Metric::set_der_0x0() const {

  p_ricci_scal = 0x0 ;

  p_determinant = 0x0 ;

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

  

void Metric::operator=(const Metric& meti) {

  assert( mp == meti.mp) ;

  if (p_connect != 0x0) delete p_connect ;

  if (meti.p_connect != 0x0) {
    p_connect = new Connection(*meti.p_connect) ;
  }
  else {
    p_connect = 0x0 ;
  }

  if (p_met_cov != 0x0) delete p_met_cov ;

  if (meti.p_met_cov != 0x0) {
    p_met_cov = new Sym_tensor(*meti.p_met_cov) ;
  }
  else {
    p_met_cov = 0x0 ;
  }

  if (p_met_con != 0x0) delete p_met_con ;

  if (meti.p_met_con != 0x0) {
    p_met_con = new Sym_tensor(*meti.p_met_con) ;
  }
  else {
    p_met_con = 0x0 ;
  }

  del_deriv() ;

}

void Metric::operator=(const Sym_tensor& symti) {
  
  assert(mp == &symti.get_mp()) ;

  if (p_connect != 0x0) delete p_connect ; //The connection will be done by 
                                           // fait_connection() later...
  int type_index = symti.get_index_type(0) ;
  assert (symti.get_index_type(1) == type_index) ;

  if (p_met_cov != 0x0) delete p_met_cov ;
  if (p_met_con != 0x0) delete p_met_con ;

  if (type_index == COV) {
    p_met_cov = new Sym_tensor(symti) ;
  }
  else {
    assert(type_index == CON) ;
    p_met_con = new Sym_tensor(symti) ;
  }

  del_deriv() ;

}
  
const Connection& Metric::get_connect() const {

  if (p_connect == 0x0) 
    fait_connection() ;

  return *p_connect ;

}

const Sym_tensor& Metric::ricci() const {

  if (p_connect == 0x0) 
    fait_connection() ;
  
  //A check of consistency (the Ricci tensor associated with the connection
  // (based upon a metric) should be symmetric.
  assert( typeid(p_connect->ricci()) == typeid(const Sym_tensor) ) ;

  return dynamic_cast<const Sym_tensor&>(p_connect->ricci()) ; 
}

const Scalar& Metric::ricci_scal() const {

  if (p_ricci_scal == 0x0)
    fait_ricci_scal() ;

  return *p_ricci_scal ;

}

const Scalar& Metric::determinant() const {

  if (p_determinant == 0x0)
    fait_determinant() ;

  return *p_determinant ;

}

void Metric::fait_cov() const {
  
  assert( p_met_cov == 0x0 ) ;
  assert( p_met_con != 0x0 ) ;

  p_met_cov = p_met_con->inverse() ;

}

void Metric::fait_con() const {

  assert( p_met_con == 0x0 ) ;
  assert( p_met_cov != 0x0 ) ;

  p_met_con = p_met_cov->inverse() ;

}

void Metric::fait_connection() const {

  assert( p_connect == 0x0 ) ;
  
  p_connect = new Connection(*this) ;

}

void Metric::fait_ricci_scal() const {

  assert( p_ricci_scal == 0x0 ) ;

  cout << "Metric::fait_ricci_scal : not implemented yet!" << endl ;
  abort() ;

}

void Metric::fait_determinant() const {

  assert( p_determinant == 0x0 ) ;

  p_determinant = new Scalar(*mp) ;
  *p_determinant = cov()(0, 0)*cov()(1, 1)*cov()(2, 2) 
	+ cov()(0, 1)*cov()(1, 2)*cov()(2, 0)
	+ cov()(0, 2)*cov()(1, 0)*cov()(2, 1) 
	- cov()(2, 0)*cov()(1, 1)*cov()(0, 2)
	- cov()(2, 1)*cov()(1, 2)*cov()(0, 0) 
	- cov()(2, 2)*cov()(1, 0)*cov()(0, 1) ;

}
 
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
  } //## what to do with the connection??
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

//##   if (p_connect == 0x0) 
//     ost << "Connection not defined!" << '\n' ;
//   else {
//     ost << "Associated Connection : " << '\n' ;
//     ost << *p_connect ;
//##   }

  if (p_ricci_scal == 0x0)
    ost << "Ricci scalar unknown." << '\n' ;
  else
    ost << "Ricci scalar known." << '\n' ;
  
  if (p_determinant == 0x0)
    ost << "determinant unknown." << '\n' ;
  else
    ost << "determinant known." << '\n' ;

  ost << endl ;
  return ost ;
}

