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
#include <stdlib.h>
#include <assert.h>
#include <math.h>



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
        while ( (val[pos_jtop] == 0x0) && (pos_jtop>=0) ) pos_jtop-- ;
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
TyT Evolution<TyT>::operator()(double ) const {

    cerr << 
    "Evolution<TyT>::operator()(double ) const : not implemented yet !"
        << endl ; 
    abort() ; 

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
            double dt2 = the_time[pos-1] - the_time[pos-2] ;
            if (fabs(dt2 -dt) > 1.e-13) {
                cerr << 
  "Evolution<TyT>::time_derive: the current version is  valid only for \n"
    << " a constant time step !" << endl ; 
                abort() ;
            }

            return ( 1.5 * (*val[pos]) - 2.* (*val[pos-1])
                        + 0.5 * (*val[pos-2]) ) / dt  ;  
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
 
                   
