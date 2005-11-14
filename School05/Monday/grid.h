/*
 * Declaration of class Grid and derived classes
 */

/*
 *   Copyright (c) 2005 Eric Gourgoulhon
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

/*
 * $Id$
 * $Log$
 * Revision 1.1  2005/11/14 01:56:58  e_gourgoulhon
 * First version
 *
 *
 * $Header$
 *
 */


#ifndef __GRID_H_ 
#define __GRID_H_ 

            //----------------------------------//
            //          class Grid              //
            //----------------------------------//

/** Set of nodes in the interval [-1,1]. \ingroup (grid)
 *
 */
class Grid {

    // Data : 
    // -----
    protected:
        const int nn ; ///< \e N = number of nodes - 1
        double* xx ;  ///< Values of the nodes (pointer to an array of size \c nn+1) 

    // Constructors - Destructor
    // -------------------------
    public:
        /** Standard constructor.
         *  @param nb_nodes number of nodes 
         *  @param xi array of size \c nb_nodes containing the abscissas of the nodes
         */
	Grid(int nb_nodes, double* xi) ;
        
	Grid(const Grid& ) ;		///< Copy constructor
        
    protected:
        /** Constructor to be used only by derived classes (hence protected)
         *  @param nb_nodes number of nodes 
         */
        Grid(int nb_nodes) ;
         
    public:
	virtual ~Grid() ;			///< Destructor
 
    // Mutators / assignment
    // ---------------------
    public:
	/// Assignment to another Grid
	void operator=(const Grid&) ;	
	
    // Accessors
    // ---------
    public:
        int n() const ; ///< returns \e N, i.e. the number of nodes - 1

        double operator()(int i) const ; ///< returns value of node no. \c i

        /** Graphical display of the node points.
        *
        *  @param color color of the points (drawn as small circles) : 
        *      \li 0 : black (background)
        *      \li 1 : white (default)
        *      \li 2 : red
        *      \li 3 : green
        *      \li 4 : blue
        *      \li 5 : cyan
        *      \li 6 : magenta
        *      \li 7 : yellow
        *      \li 8 : orange
        *      \li 9 : green + yellow
        *      \li 10 : green + cyan
        *      \li 11 : blue + cyan
        *      \li 12 : blue + magenta
        *      \li 13 : red + magenta
        *      \li 14 : dark gray 
        *      \li 15 : light gray
        *  @param nfig index of the figure (in the range [0,99])
        *  to be used for the plot: if this figure does not exist, 
        *    it will be created with the device name \c device  provided by the last
        *      argument. 
        *  @param ymin lower bound on y of the graphical window (used only if a new 
        *      figure must be created)
        *  @param ymax upper bound on y of the graphical window (used only if a new 
        *      figure must be created)
        *  @param title title of the figure (used only if a new figure must be created)
        *  @param label_y y legend of the figure (used only if a new 
        *      figure must be created)
        *  @param device type of graphical device (default value = 0x0, will result in
        *  interactive choice) (used only if a new 
        *      figure must be created)
        */
        void plot(int color = 1, int nfig = 0, double ymin = -1., 
                  double ymax = 1., const char* title = 0x0, 
                  const char* label_y = 0x0, const char* device = 0x0) const ;
                        
    // Computation
    // -----------
    public: 
        /** Lagrange polynomials (characteristic polynomials).
         *  The Lagrange polynomial no. \e i is the unique polynomial 
         *  \f$\ell_i(x)\f$ of degree \e N such that 
         *  \f$\ell_i(x_{j})=\delta_{ij}\f$. It is also called \e cardinal
         *  \e polynomial associated with the node \f$x_i\f$
         *  @param i number of Lagrange polynomial \f$\ell_i\f$
         *  @param x point \e x in [-1,1] where the Lagrange polynomial is to 
         *      be evaluated
         *  @return value of \f$\ell_i(x)\f$
         */
        double lagrange(int i, double x) const ; 
         
        /** Nodal polynomial. 
         * The nodal polynomoal is the unique polynomial of degree \e N + 1
         *  and leading coefficient 1, the roots of which are the \e N +1 nodes
         *  \f$x_i\f$.
         *
         *  @param x point \e x in [-1,1] where the Lagrange polynomial is to 
         *      be evaluated
         *  @return value of the nodal polynomial at \e x
         */
        double nodal_polynomial(double x) const ; 
          
        /** Interpolation of a given function at the nodes
         *
         *  @param (*f)(double) function to be interpolated         
         *  @param x point \e x in [-1,1] where the interpolating polynomial
         *    associated with the nodes \f$x_i\f$ is to evaluated
         *  @return value of the interpolating polynomial at \e x       
         */
        double interpole(double (*f)(double), double x) const ;  
        
        /// Computes (an approximate value of) the Lebesgue constant
        double lebesgue_constant() const ; 
    
};

/// Display of an object of class \c Grid
ostream& operator<<(ostream& , const Grid& ) ;	


            //----------------------------------//
            //          class Grid_uniform      //
            //----------------------------------//



/** Uniform sampling of the interval [-1,1]. \ingroup (grid)
 *
 */
class Grid_uniform : public Grid {

    // Constructors - Destructor
    // -------------------------
    public:
        /** Standard constructor.
         *  @param nb_nodes number of nodes 
         */
	Grid_uniform(int nb_nodes) ;
        
	Grid_uniform(const Grid_uniform& ) ;		///< Copy constructor
        
	virtual ~Grid_uniform() ;			///< Destructor
 
    // Mutators / assignment
    // ---------------------
    public:
	/// Assignment to another Grid_uniform
	void operator=(const Grid_uniform&) ;	
	
                       
};

/// Display of an object of class \c Grid_uniform
ostream& operator<<(ostream& , const Grid_uniform& ) ;	


            //----------------------------------//
            //    class Grid_Legendre_Gauss    //
            //----------------------------------//


/** Legendre Gauss nodes in the interval [-1,1]. \ingroup (grid)
 *
 */
class Grid_Legendre_Gauss : public Grid {

    // Constructors - Destructor
    // -------------------------
    public:
        /** Standard constructor.
         *  @param nb_nodes number of nodes 
         */
        Grid_Legendre_Gauss(int nb_nodes) ;
        
        Grid_Legendre_Gauss(const Grid_Legendre_Gauss& ) ;	///< Copy constructor
        
	    virtual ~Grid_Legendre_Gauss() ;			///< Destructor
 
    // Mutators / assignment
    // ---------------------
    public:
	    /// Assignment to another Grid_Legendre_Gauss
	    void operator=(const Grid_Legendre_Gauss&) ;	
	
                       
};

/// Display of an object of class \c Grid_Legendre_Gauss
ostream& operator<<(ostream& , const Grid_Legendre_Gauss& ) ;	


            //----------------------------------//
            //     class Grid_Legendre_GL      //
            //----------------------------------//


/** Legendre Gauss-Lobatto nodes in the interval [-1,1]. \ingroup (grid)
 *
 */
class Grid_Legendre_GL : public Grid {

    // Constructors - Destructor
    // -------------------------
    public:
        /** Standard constructor.
         *  @param nb_nodes number of nodes 
         */
        Grid_Legendre_GL(int nb_nodes) ;
        
        Grid_Legendre_GL(const Grid_Legendre_GL& ) ;	///< Copy constructor
        
	    virtual ~Grid_Legendre_GL() ;			///< Destructor
 
    // Mutators / assignment
    // ---------------------
    public:
	    /// Assignment to another Grid_Legendre_GL
	    void operator=(const Grid_Legendre_GL&) ;	
	
                       
};

/// Display of an object of class \c Grid_Legendre_GL
ostream& operator<<(ostream& , const Grid_Legendre_GL& ) ;	



            //----------------------------------//
            //    class Grid_Chebyshev_Gauss    //
            //----------------------------------//


/** Chebyshev Gauss nodes in the interval [-1,1]. \ingroup (grid)
 *
 */
class Grid_Chebyshev_Gauss : public Grid {

    // Constructors - Destructor
    // -------------------------
    public:
        /** Standard constructor.
         *  @param nb_nodes number of nodes 
         */
        Grid_Chebyshev_Gauss(int nb_nodes) ;
        
        Grid_Chebyshev_Gauss(const Grid_Chebyshev_Gauss& ) ;	///< Copy constructor
        
	    virtual ~Grid_Chebyshev_Gauss() ;			///< Destructor
 
    // Mutators / assignment
    // ---------------------
    public:
	    /// Assignment to another Grid_Chebyshev_Gauss
	    void operator=(const Grid_Chebyshev_Gauss&) ;	
	
                       
};

/// Display of an object of class \c Grid_Chebyshev_Gauss
ostream& operator<<(ostream& , const Grid_Chebyshev_Gauss& ) ;	


            //----------------------------------//
            //     class Grid_Chebyshev_GL      //
            //----------------------------------//


/** Chebyshev Gauss-Lobatto nodes in the interval [-1,1]. \ingroup (grid)
 *
 */
class Grid_Chebyshev_GL : public Grid {

    // Constructors - Destructor
    // -------------------------
    public:
        /** Standard constructor.
         *  @param nb_nodes number of nodes 
         */
        Grid_Chebyshev_GL(int nb_nodes) ;
        
        Grid_Chebyshev_GL(const Grid_Chebyshev_GL& ) ;	///< Copy constructor
        
	    virtual ~Grid_Chebyshev_GL() ;			///< Destructor
 
    // Mutators / assignment
    // ---------------------
    public:
	    /// Assignment to another Grid_Chebyshev_GL
	    void operator=(const Grid_Chebyshev_GL&) ;	
	
                       
};

/// Display of an object of class \c Grid_Chebyshev_GL
ostream& operator<<(ostream& , const Grid_Chebyshev_GL& ) ;	

#endif
