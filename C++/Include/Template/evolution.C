/*
 *  Methods of template class Evolution
 *
 *    (see file evolution.h for documentation).
 *
 */

/*
 *   Copyright (c) 2004  Eric Gourgoulhon & Jerome Novak
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
 * Revision 1.18  2014/10/13 08:52:38  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.17  2014/03/27 16:59:41  j_novak
 * Added methods next_position(int) and previous_position(int). Changed (corrected + simplified) the interpolation method.
 *
 * Revision 1.16  2013/07/19 15:50:24  j_novak
 * Implementation of the interpolation function for Evolution, with order=0, 1 or 2.
 *
 * Revision 1.15  2008/09/19 13:24:29  j_novak
 * Added the third-order scheme for the time derivativre computation.
 *
 * Revision 1.14  2004/05/14 08:51:01  p_grandclement
 * *** empty log message ***
 *
 * Revision 1.13  2004/05/14 08:41:05  e_gourgoulhon
 * Added declaration "class Tbl ;" before the declarations of
 * write_formatted.
 *
 * Revision 1.12  2004/05/13 21:30:32  e_gourgoulhon
 * Use of function write_formatted in method save( ).
 *
 * Revision 1.11  2004/05/11 20:12:49  e_gourgoulhon
 * Added methods j_min, j_max and save.
 *
 * Revision 1.10  2004/05/03 15:23:22  e_gourgoulhon
 * Method downdate: changed the order of conditions (pos_jtop>=0)
 * and (val[pos_jtop] == 0x0) in the while(...) test.
 *
 * Revision 1.9  2004/03/31 20:24:04  e_gourgoulhon
 * Method time_derive: result object created by the copy constructor of
 * class TyT, since the arithmetics may not return an object of exactly
 * class TyT.
 *
 * Revision 1.8  2004/03/26 13:31:09  j_novak
 * Definition of the macro UNDEF_STEP for non-defined time-steps.
 * Changes in the way the time derivative is calculated.
 *
 * Revision 1.7  2004/03/26 08:22:13  e_gourgoulhon
 * *** Full reorganization of class Evolution ***
 * Introduction of the notion of absoluteuniversal time steps,
 * stored in the new array 'step'.
 * The new function position(int j) makes a correspondence
 * between a universal time step j and the position in the
 * arrays step, the_time and val.
 * Only method update is now virtual.
 * Methods operator[], position, is_known, downdate belong to
 * the base class.
 *
 * Revision 1.6  2004/03/24 14:55:47  e_gourgoulhon
 * Added method last_value().
 *
 * Revision 1.5  2004/03/23 14:50:41  e_gourgoulhon
 * Added methods is_updated, downdate, get_jlast, get_size,
 * as well as constructors without any initial value.
 * Formatted documentation for Doxygen.
 *
 * Revision 1.4  2004/03/06 21:13:15  e_gourgoulhon
 * Added time derivation (method time_derive).
 *
 * Revision 1.3  2004/02/17 22:13:34  e_gourgoulhon
 * Suppressed declaration of global char[] evolution_C = ...
 *
 * Revision 1.2  2004/02/15 21:55:33  e_gourgoulhon
 * Introduced derived classes Evolution_full and Evolution_std.
 * Evolution is now an abstract base class.
 *
 * Revision 1.1  2004/02/13 15:53:20  e_gourgoulhon
 * New (template) class for time evolution.
 *
 *
 * $Header$
 *
 */

// C++ headers
#include "headcpp.h" 

// C headers
#include <cstdlib>
#include <cassert>
#include <cmath>

namespace Lorene {
class Tbl ;

void write_formatted(const double&, ostream& ) ; 
void write_formatted(const Tbl&, ostream& ) ; 


                    //-------------------------//
                    //      Constructors       //
                    //-------------------------//

                    
template<typename TyT> 
Evolution<TyT>::Evolution(const TyT& initial_value, int initial_step,
                          double initial_time, int size_i)
      : size(size_i),
        pos_jtop(0) {

    step = new int[size] ; 
    step[0] = initial_step ; 
    for (int j=1; j<size; j++) {
        step[j] = UNDEF_STEP ; 
    }
    
    the_time = new double[size] ; 
    the_time[0] = initial_time ; 
    for (int j=1; j<size; j++) {
        the_time[j] = -1e20 ; 
    }
    
    val = new TyT*[size] ; 
    val[0] = new TyT(initial_value) ; 
    for (int j=1; j<size; j++) {
        val[j] = 0x0 ; 
    }
        
}                    

                    
template<typename TyT> 
Evolution<TyT>::Evolution(int size_i)
      : size(size_i),
        pos_jtop(-1) {

    step = new int[size] ; 
    for (int j=0; j<size; j++) {
        step[j] = UNDEF_STEP ; 
    }
    
    the_time = new double[size] ; 
    for (int j=0; j<size; j++) {
        the_time[j] = -1e20 ; 
    }
    
    val = new TyT*[size] ; 
    for (int j=0; j<size; j++) {
        val[j] = 0x0 ; 
    }    
    
}                    
                    


template<typename TyT> 
Evolution<TyT>::Evolution(const Evolution<TyT>& evo)
      : size(evo.size),
        pos_jtop(evo.pos_jtop) {

    step = new int[size] ; 
    for (int j=0; j<size; j++) {
        step[j] = evo.step[j] ; 
    }
    
    the_time = new double[size] ; 
    for (int j=0; j<size; j++) {
        the_time[j] = evo.the_time[j] ; 
    }
    
    val = new TyT*[size] ; 
    for (int j=0; j<size; j++) {
        if (evo.val[j] != 0x0) {
            val[j] = new TyT( *(evo.val[j]) ) ; 
        }
        else {
            val[j] = 0x0 ; 
        }
    }
    
    
}                    
                    

                    
                    
                    //-----------------------//
                    //      Destructor       //
                    //-----------------------//

                    
template<typename TyT> 
Evolution<TyT>::~Evolution(){

    delete [] step ; 
    delete [] the_time ; 

    for (int j=0; j<size; j++) {
        if (val[j] != 0x0) delete val[j] ; 
    }
    
    delete [] val ;
    
}
                    
                    

                    //-----------------------//
                    //      Mutators         //
                    //-----------------------//

                    
template<typename TyT> 
void Evolution<TyT>::operator=(const Evolution<TyT>& ) {

    cerr << "void Evolution<TyT>::operator= : not implemented yet ! \n" ; 
    abort() ; 

}


template<typename TyT> 
void Evolution<TyT>::downdate(int j) {

    if ( !(is_known(j)) ) return ;  // a never updated step cannot
                                    // be downdated
    
    int pos = position(j) ; 
    
    assert( val[pos] != 0x0) ; 

    delete val[pos] ; 
    val[pos] = 0x0 ; 
    step[pos] = UNDEF_STEP ; 
    the_time[pos] = -1e20 ; 

    if (pos == pos_jtop) {  // pos_jtop must be decreased
        pos_jtop-- ; 
        while ( (pos_jtop>=0) && (val[pos_jtop] == 0x0) ) pos_jtop-- ;
    }
    
}


                    
                    
                        //------------//
                        // Accessors  //
                        //------------//

template<typename TyT> 
int Evolution<TyT>::position(int j) const {
    
    assert(pos_jtop >= 0) ; 
    int jmax = step[pos_jtop] ; 
    
    if (j == jmax) return pos_jtop ;   // for efficiency purpose
    
    int pos = - 1 ; 

    if ( (j>=step[0]) && (j<jmax) ) {

        for (int i=pos_jtop-1; i>=0; i--) {  // cas i=pos_jtop treated above
            if (step[i] == j) {
                pos = i ;
                break ; 
            }
        }
    }
    
    if (pos == -1) {
        cerr << "Evolution<TyT>::position: time step j = " <<
            j << " not found !" << endl ; 
        abort() ; 
    }
    
    return pos ; 
}
   
template<typename TyT> 
int Evolution<TyT>::next_position(int j) const {

  assert( (j>=0) && (j<=pos_jtop) ) ;
  assert( step[j] != UNDEF_STEP ) ;

  int n_pos = -1 ;

  while ( (n_pos == -1) && ( j < pos_jtop ) ) {
    
    j++ ;
    if (step[j] != UNDEF_STEP) n_pos = j ;
  }
  return n_pos ;
}

template<typename TyT> 
int Evolution<TyT>::previous_position(int j) const {              
                    
  assert( (j>=0) && (j<=pos_jtop) ) ;
  assert( step[j] != UNDEF_STEP ) ;

  int n_pos = -1 ;

  while ( (n_pos == -1) && ( j > 0 ) ) {
    
    j-- ;
    if (step[j] != UNDEF_STEP) n_pos = j ;
  }
  return n_pos ;
}


template<typename TyT> 
bool Evolution<TyT>::is_known(int j) const {

    if (pos_jtop == -1) return false ; 
    
    assert(pos_jtop >= 0) ; 
    
    int jmax = step[pos_jtop] ; 
    
    if (j == jmax) {
        return ( val[pos_jtop] != 0x0 ) ; 
    }

    if ( (j>=step[0]) && (j<jmax) ) {

        for (int i=pos_jtop-1; i>=0; i--) {  // cas i=pos_jtop treated above

            if (step[i] == j) return ( val[i] != 0x0 ) ; 
        }
    }
    
    return false ; 
}                  
                    

template<typename TyT> 
const TyT& Evolution<TyT>::operator[](int j) const {

    TyT* pval = val[position(j)] ; 
    assert(pval != 0x0) ; 
    return *pval ; 

}


template<typename TyT> 
TyT Evolution<TyT>::operator()(double t, int order) const {

  int imin = position( j_min() ) ;

  if ( ( t < the_time[ imin ] ) || ( t > the_time[pos_jtop] ) ) {
    cerr << "Evolution<TyT>::operator()(double t, int order) : \n"
	 << "Requested time outside stored range!" << endl ;
    abort() ;
  }
  assert ( order <= 2 ) ;
  assert ( pos_jtop > order ) ;

  int j0 = imin ;

  while ((the_time[j0] < t) && ( j0 < pos_jtop )) {
    j0 = next_position(j0) ;
    assert( j0 != -1 ) ;
  }

  switch (order) {

  case 0: {

    return (*val[j0]) ;

    break ;
  }
    
  case 1: {
    
    int j1 = ( (j0 > imin) ? previous_position(j0) : next_position(j0) ) ;
    assert( j1 != -1) ;

    double x0 = the_time[j1] - t ;
    double x1 = the_time[j0] - t ;
    double dx = the_time[j0] - the_time[j1] ;
    
    double a0 = x1 / dx ;
    double a1 = x0 / dx ;

    return (a0*(*val[j0]) - a1*(*val[j1])) ;

    break ;
  }

  case 2: {

    int j1 = ( (j0 == imin) ? next_position(j0) : previous_position(j0) ) ;
    assert( j1 != -1) ;

    int j2 = -1 ;
    if (j0 == imin) 
      j2 = next_position(j1) ;
    else {
      if (j1 == imin) 
	j2 = next_position(j0) ;
      else
	j2 = previous_position(j1) ;
    }
    assert (j2 != -1 ) ;

    double x0 = t - the_time[j0] ;
    double x1 = t - the_time[j1] ;
    double dx = the_time[j0] - the_time[j1] ;
    
    double x2 = t - the_time[j2] ;

    double dx0 = the_time[j1] - the_time[j2] ;
    double dx1 = the_time[j0] - the_time[j2] ;
    double dx2 = dx ;

    double a0 = ( x2*x1 ) / ( dx2*dx1 ) ;
    double a1 = ( x0*x2 ) / ( dx0*dx2 ) ;
    double a2 = ( x0*x1 ) / ( dx0*dx1 ) ;

    return ( a0*(*val[j0]) - a1*(*val[j1]) + a2*(*val[j2]) ) ;

    break ;
  }

  default: {
    cerr << " Evolution<TyT>::operator()(double t, int order) : \n" << endl ;
    cerr << " The case order = " << order << " is not implemented!" << endl ;
    abort() ;
    break ;
  }
  }

  return *val[j0] ;

}                  
 
template<typename TyT> 
int Evolution<TyT>::j_min() const {

    int resu = UNDEF_STEP ; 
    for (int i=0; i<=pos_jtop; i++) {
        if (step[i] != UNDEF_STEP) {
            resu = step[i] ; 
            break ; 
        }
    }
    
    if (resu == UNDEF_STEP) {
        cerr << "Evolution<TyT>::j_min() : no valid time step found !" << endl ; 
        abort() ; 
    }
    
    return resu ; 
}

template<typename TyT> 
int Evolution<TyT>::j_max() const {
    
    if (pos_jtop == -1) {
        cerr << "Evolution<TyT>::j_max() : no valid time step found !" << endl ; 
        abort() ; 
    }
    
    assert(pos_jtop >=0) ; 
    int jmax = step[pos_jtop] ; 
    assert(jmax != UNDEF_STEP) ; 
    return jmax ; 
}





                    //-----------------------//
                    //   Time derivative     //
                    //-----------------------//

template<typename TyT> 
TyT Evolution<TyT>::time_derive(int j, int n) const {

  if (n == 0) { 
    TyT resu ( operator[](j) ) ;
    resu = 0 * resu ;

    return resu ;

  }
  
  else { 
    
    int pos = position(j) ;
    assert ( pos > 0 ) ;
    assert ( step[pos-1] != UNDEF_STEP ) ;

    switch (n) {
    
        case 1 : {

            double dt = the_time[pos] - the_time[pos-1] ; 

            return ( (*val[pos]) - (*val[pos-1]) ) / dt  ;
            break ;
        } 
           
        case 2 : {
	  
	  assert ( pos > 1 ) ;
	  assert ( step[pos-2] != UNDEF_STEP ) ;
	  double dt = the_time[pos] - the_time[pos-1] ;
#ifndef NDEBUG
	  double dt2 = the_time[pos-1] - the_time[pos-2] ;
	  if (fabs(dt2 -dt) > 1.e-13) {
	      cerr << 
		  "Evolution<TyT>::time_derive: the current version is"
		   << " valid only for \n"
		   << " a constant time step !" << endl ; 
	      abort() ;
	  }
#endif
	  return ( 1.5* (*val[pos]) - 2.* (*val[pos-1]) + 0.5* (*val[pos-2]) ) / dt ;  
	  break ;
        } 
	    
        case 3 : {
	    
	    assert ( pos > 2 ) ;
	    assert ( step[pos-2] != UNDEF_STEP ) ;
	    assert ( step[pos-3] != UNDEF_STEP ) ;
	    double dt = the_time[pos] - the_time[pos-1] ;
#ifndef NDEBUG
	    double dt2 = the_time[pos-1] - the_time[pos-2] ;
	    double dt3 = the_time[pos-2] - the_time[pos-3] ;
            if ((fabs(dt2 -dt) > 1.e-13)||(fabs(dt3 -dt2) > 1.e-13)) {
                cerr << 
  "Evolution<TyT>::time_derive: the current version is  valid only for \n"
    << " a constant time step !" << endl ; 
                abort() ;
            }
#endif
	    return ( 11.*(*val[pos]) - 18.*(*val[pos-1]) + 9.*(*val[pos-2])
		     - 2.*(*val[pos-3])) / (6.*dt) ;
	    break ;
	}
        default : {
            cerr << "Evolution<TyT>::time_derive: the case n = " << n 
                 << "  is not implemented !" << endl ; 
            abort() ;
            break ;
        }    
    }

  }
  return operator[](j) ;
}                    
 
                   

                    //-----------------------//
                    //        Outputs        //
                    //-----------------------//

                   
template<typename TyT> 
void Evolution<TyT>::save(const char* filename) const {

    ofstream fich(filename) ; 

    time_t temps = time(0x0) ; 
    
    fich << "# " << filename << "    " << ctime(&temps) ; 
    fich << "# " << size << "  size" << '\n' ; 
    fich << "# " << pos_jtop << "  pos_jtop" << '\n' ; 
    fich << "#         t                  value...                   \n" ; 
    
    fich.precision(14) ; 
    fich.setf(ios::scientific) ; 
    fich.width(20) ; 
    for (int i=0; i<=pos_jtop; i++) {
        if (step[i] != UNDEF_STEP) {
            fich << the_time[i] ; fich.width(23) ; 
            assert(val[i] != 0x0) ;
            write_formatted(*(val[i]), fich) ;  
            fich << '\n' ; 
        }
    }
        
    fich.close() ; 
}


















}
