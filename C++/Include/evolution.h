/*
 *  Definition of Lorene template classes Evolution, Evolution_full
 *  and Evolution_std
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


#ifndef __EVOLUTION_H_ 
#define __EVOLUTION_H_ 

/*
 * $Id$
 * $Log$
 * Revision 1.5  2004/03/06 21:13:13  e_gourgoulhon
 * Added time derivation (method time_derive).
 *
 * Revision 1.4  2004/02/16 17:37:17  j_novak
 * Arguments named for doc++.
 *
 * Revision 1.3  2004/02/16 10:36:03  e_gourgoulhon
 * Replaced " = 0x0" by " = 0" in the declaration of pure virtual functions.
 *
 * Revision 1.2  2004/02/15 21:55:32  e_gourgoulhon
 * Introduced derived classes Evolution_full and Evolution_std.
 * Evolution is now an abstract base class.
 *
 * Revision 1.1  2004/02/13 15:53:20  e_gourgoulhon
 * New (template) class for time evolution.
 *
 *
 *
 *
 * $Header$
 *
 */

                        //---------------------------//
                        //      Class Evolution      //
                        //---------------------------//


/** Time evolution (*** under development ***). 
 * 
 * The template class {\tt Evolution} has been devised to store and
 * manipulate evolving quantities of any type, for instance {\tt TyT = double}
 * or {\tt TyT = Scalar}.
 *
 * {\tt Evolution} is an abstract base class for classes
 * {\tt Evolution\_full} and {\tt Evolution\_std}. 
 *
 */
template<typename TyT> class Evolution {

        
    // Data:
    // -----
    
    protected: 
        /// Maximum number of stored time steps
        int size ; 
        
        /// Array of pointers onto the values (size {\tt size}) 
        TyT** val ; 
        
        /// Array of successive time steps
        double* the_time ; 
      
        /// Index of last updated time step for the object
        int jlast ; 
      
        
    // Constructors - Destructor
    // -------------------------
    protected:
        /** Constructor (to be used by derived classes)
         *
         */
        Evolution(const TyT& initial_value, double initial_time, int size_i) ;			
	
        Evolution(const Evolution<TyT>& t_in) ;		/// Copy constructor

    public: 

	virtual ~Evolution() ;			/// Destructor
 
    // Mutators 
    // --------
        /// Update
        virtual void update(const TyT& new_value, double new_time) = 0 ; 
    
        /// Assignement
        virtual void operator=(const Evolution<TyT>& t_in) ;

    
    // Accessors
    // ---------
        /// Returns value at time step j
        virtual const TyT& operator[](int j) const = 0 ;

        /// Returns value at time t
        TyT operator()(double t) const ;

        /// Returns the time t at time step j
        virtual double get_time(int j) const = 0 ;

    // Computational methods
    // ---------------------
        /** Computes the time derivative at time step {\tt j} by means of a 
         *  n-th order scheme, from the values at steps {\tt j}, 
         *  {\tt j-1}, ..., {\tt j-n}.
         * 
         * @param j [input] : value of the time step at which the time
         *      derivative is required
         * @param n [input] : order of the time scheme (default value=2)
         * @return time derivative at time step {\tt j} 
         *   
         */
        TyT time_derive(int j, int n = 2) const ; 

    // Outputs
    // -------
    
    

};


                        //---------------------------//
                        //   Class Evolution_full    //
                        //---------------------------//


/** Time evolution with full storage (*** under development ***). 
 * 
 * The template class {\tt Evolution\_full} has been devised to store and
 * manipulate evolving quantities of any type, for instance {\tt TyT = double}
 * or {\tt TyT = Scalar}.
 * The quantity is stored at all time steps since the beginning of the
 * time evolution. For large objects, this might result in some memory
 * problem. The class {\tt Evolution\_std}, which stores only a limited
 * number of time steps, is to be prefered then.  
 *
 */
template<typename TyT> class Evolution_full : public Evolution<TyT> {

        
    // Data:
    // -----
    
    private:
        /** Factor by which the size {\tt size} of the arrays 
         *  {\tt val} and {\tt the\_time} are to be multiplied when 
         *  their limits have been reached.
         */        
         int fact_resize ; 
        
    // Constructors - Destructor
    // -------------------------
    public:
        /** Standard constructor.
         *  
         * @param initial_value value to be stored at time step 0
         * @param initial_time  time t corresponding to time step 0
         * @param fact_resize_i factor by which the size {\tt size} of the arrays 
         *  {\tt val} and {\tt the\_time} are to be multiplied when 
         *  their limits have been reached.
         *
         */
        Evolution_full(const TyT& initial_value, double initial_time, 
                  int fact_resize_i = 2) ;			
	
	
        Evolution_full(const Evolution_full<TyT>& t_in) ;  /// Copy constructor

        virtual ~Evolution_full() ;			/// Destructor
 
    // Mutators 
    // --------
        /// Update
        virtual void update(const TyT& new_value, double new_time) ; 
    
        /// Assignement to another Evolution\_full
        virtual void operator=(const Evolution_full<TyT>& t_in) ;

        /// Assignement to a generic Evolution
        virtual void operator=(const Evolution<TyT>& t_in) ;

    
    // Accessors
    // ---------
        /// Returns value at time step j
        virtual const TyT& operator[](int j) const ;

        /// Returns the time t at time step j
        virtual double get_time(int j) const ;

    // Outputs
    // -------
    
    

};


                        //---------------------------//
                        //   Class Evolution_std     //
                        //---------------------------//


/** Time evolution with partial storage (*** under development ***). 
 * 
 * The template class {\tt Evolution\_std} has been devised to store and
 * manipulate evolving quantities of any type, for instance {\tt TyT = double}
 * or {\tt TyT = Scalar}.
 * The quantity is stored only for a limited number of time steps (the
 * n last ones).
 * For a full storage, use instead the class {\tt Evolution\_full}.
 *
 */
template<typename TyT> class Evolution_std : public Evolution<TyT> {

        
    // Data:
    // -----
    protected:
    
        /** Position of the time step jlast in the arrays 
         * {\tt val} and {\tt the\_time}
         */
        int pos_jlast ; 
    
        
    // Constructors - Destructor
    // -------------------------
    public:
        /** Standard constructor.
         *  
         * @param initial_value value to be stored at time step 0
         * @param initial_time  time t corresponding to time step 0
         * @param nstored total number of time steps to be stored
         *
         */
        Evolution_std(const TyT& initial_value, double initial_time, 
                  int nstored) ;			
	
	
        Evolution_std(const Evolution_std<TyT>& t_in) ;	/// Copy constructor

        virtual ~Evolution_std() ;			/// Destructor
 
    // Mutators 
    // --------
        /// Update
        virtual void update(const TyT& new_value, double new_time) ; 
    
        /// Assignement to another Evolution\_std
        virtual void operator=(const Evolution_std<TyT>& t_in) ;

        /// Assignement to a generic Evolution
        virtual void operator=(const Evolution<TyT>& t_in) ;

    // Accessors
    // ---------
        /// Returns value at time step j
        virtual const TyT& operator[](int j) const ;

        /// Returns the time t at time step j
        virtual double get_time(int j) const ;

    // Outputs
    // -------
    
    

};


#include "Template/evolution.C"
#include "Template/evolution_full.C"
#include "Template/evolution_std.C"

#endif

