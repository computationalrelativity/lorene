/*
 *  Methods of template class Evolution_std
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

char evolution_std_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2004/02/16 12:37:34  e_gourgoulhon
 * Added an assert in method update.
 *
 * Revision 1.1  2004/02/15 21:55:33  e_gourgoulhon
 * Introduced derived classes Evolution_full and Evolution_std.
 * Evolution is now an abstract base class.
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



                    //-------------------------//
                    //      Constructors       //
                    //-------------------------//

                    
template<typename TyT> 
Evolution_std<TyT>::Evolution_std(const TyT& initial_value, double initial_time,
    int nstored) 
      : Evolution<TyT>(initial_value, initial_time, nstored), 
        pos_jlast(0)
{ }                    
                                        

template<typename TyT> 
Evolution_std<TyT>::Evolution_std(const Evolution_std<TyT>& evo)
      : Evolution<TyT>(evo), 
        pos_jlast(evo.pos_jlast) 
{ }


                    
                    
                    //-----------------------//
                    //      Destructor       //
                    //-----------------------//

                    
template<typename TyT> 
Evolution_std<TyT>::~Evolution_std(){ }

                    
                    
                    //-----------------------//
                    //      Mutators         //
                    //-----------------------//

                    
template<typename TyT> 
void Evolution_std<TyT>::operator=(const Evolution_std<TyT>& ) {

    cerr << "void Evolution_std<TyT>::operator= : not implemented yet ! \n" ; 
    abort() ; 

}

template<typename TyT> 
void Evolution_std<TyT>::operator=(const Evolution<TyT>& ) {

    cerr << "void Evolution_std<TyT>::operator= : not implemented yet ! \n" ; 
    abort() ; 

}


                    
template<typename TyT> 
void Evolution_std<TyT>::update(const TyT& new_value, double new_time) {

    jlast++ ; 
    pos_jlast++ ; 
         
    if (pos_jlast == size) {  // re-organization of arrays val and the_time is necessary
    
        assert( val[0] != 0x0 ) ; 
        delete val[0] ; 
            
        for (int j=0; j<size-1; j++) {
            val[j] = val[j+1] ; 
        }
        
        for (int j=0; j<size-1; j++) {
            the_time[j] = the_time[j+1] ; 
        }
            
        pos_jlast-- ;
        
    }
    else {
        assert( pos_jlast < size ) ; 
    }
    
    val[pos_jlast] = new TyT( new_value ) ; 
    the_time[pos_jlast] = new_time ; 

}                   
                    
                    
                    //-----------------------//
                    //      Accessors        //
                    //-----------------------//

                 
template<typename TyT> 
const TyT& Evolution_std<TyT>::operator[](int j) const {

    int pos = j - jlast + pos_jlast ; 

    assert(pos >= 0) ;
    assert(pos < size) ; 
    
    TyT* pval = val[pos] ; 
    assert(pval != 0x0) ; 
    
    return *pval ; 

}                  
                    

template<typename TyT> 
double Evolution_std<TyT>::get_time(int j) const {

    int pos = j - jlast + pos_jlast ; 

    assert(pos >= 0) ;
    assert(pos < size) ; 
            
    return the_time[pos] ; 

}                  
                    
