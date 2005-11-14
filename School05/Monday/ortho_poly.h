/*
 * Declaration of class Ortho_poly and derived classes
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
 * Revision 1.1  2005/11/14 01:57:00  e_gourgoulhon
 * First version
 *
 *
 * $Header$
 *
 */

#ifndef __ORTHO_POLY_H_ 
#define __ORTHO_POLY_H_ 

class Grid ; 

            //----------------------------------//
            //          class Ortho_poly        //
            //----------------------------------//

/** Base class for orthogonal polynomials on [-1,1]. \ingroup (poly)
 * Storage of polynomials of an orthogonal family \f$(p_i)\f$ up to a 
 * given degree \e N.
 * The class \c Ortho_poly is abstract: it cannot be instanciated. 
 *
 */
class Ortho_poly {

    // Data : 
    // -----
    protected:
        const int nn ; ///< \e N = maximum degree of the polynomials
        
        /** Pointer on the Gauss nodes
         */
        mutable const Grid* p_gauss_nodes ; 

        /** Pointer on the Gauss weights
         */
        mutable double* p_gauss_weights ; 

        /** Pointer on the gamma factors (squares of the polynomials with
         *  respect to the discrete scalar product associated with the
         *   Gauss nodes)
         */
        mutable double* p_gauss_gamma ; 

        /** Pointer on the Gauss-Lobatto nodes
         */
        mutable const Grid* p_gauss_lobatto_nodes ; 
        
        /** Pointer on the Gauss-Lobatto weights
         */
        mutable double* p_gauss_lobatto_weights ; 

        /** Pointer on the gamma factors (square of the polynomials with
         *  respect to the discrete scalar product associated with the
         *   Gauss-Lobatto nodes)
         */
        mutable double* p_gauss_lobatto_gamma ; 

        
    // Constructors - Destructor
    // -------------------------
    protected:
	    Ortho_poly(const Ortho_poly& ) ; ///< Copy constructor
        
        /** Constructor to be used only by derived classes (hence protected)
         *  @param ni maximum degree \e N of the polynomials 
         */
        Ortho_poly(int ni) ;
         
    public:
	    virtual ~Ortho_poly() ;			///< Destructor
 
    // Mutators / assignment
    // ---------------------
    public:
	    /// Assignment to another Ortho_poly
	    void operator=(const Ortho_poly&)  ;	
	
    // Accessors
    // ---------
    public:
        int n() const ; ///< returns \e N, i.e. the maximum degree of the polynomials
        
    // Computation
    // -----------
    public: 
        /** Weight function \f$w\f$.
         *
         * @param x point in [-1,1] where the weight function is to be evaluated
         * @return value of  \f$w(x)\f$
         */
        virtual double weight(double x) const = 0 ; 

        /** Value of the polynomial \f$p_i\f$
         *
         * @param i no. of the polynomial
         * @param x point in [-1,1] where the polynomial is to be evaluated
         * @return value of  \f$p_i(x)\f$
         *
         */
        virtual double operator()(int i, double x) const = 0 ; 
        
        /// Gauss nodes
        virtual const Grid& gauss_nodes() const = 0 ;
        
        /** Gauss weights
         * 
         * @param i index of the weight
         * @return weight \f$w_i\f$
         */
        virtual double gauss_weight(int i) const = 0 ;
        
        /** Gamma factor \f$\gamma_i\f$ (square of the polynomial 
         * number \e i with
         *  respect to the discrete scalar product associated with the
         *   Gauss nodes)
         */
        virtual double gauss_gamma(int i) const = 0 ;

        /// Gauss-Lobatto nodes
        virtual const Grid& gauss_lobatto_nodes() const = 0 ;
        
        /** Gauss-Lobatto weights
         * 
         * @param i index of the weight
         * @return weight \f$w_i\f$
         */
        virtual double gauss_lobatto_weight(int i) const = 0 ;
        
        /** Gamma factor \f$\gamma_i\f$ (square of the polynomial number
         *  \e i with
         *  respect to the discrete scalar product associated with the
         *   Gauss-Lobatto nodes)
         */
        virtual double gauss_lobatto_gamma(int i) const = 0 ;

        /** Coefficients of the interpolant polynomial through the 
         *  Gauss nodes.
         *
         *  @param (*f)(double) function to be interpolated         
         *  @param cf [output] array containing the \e N+1 coefficients
         *          of the expansion of the interpolant on the basis polynomials.
         *          This array must be of size at least \e N+1 and must have
         *          been allocated by the user prior to the call of this method.       
         */
        void coef_interpolant_Gauss(double (*f)(double), double* cf) const ;  

        /** Coefficients of the interpolant polynomial through the 
         *  Gauss-Lobatto nodes.
         *
         *  @param (*f)(double) function to be interpolated         
         *  @param cf [output] array containing the \e N+1 coefficients
         *          of the expansion of the interpolant on the basis polynomials.
         *          This array must be of size at least \e N+1 and must have
         *          been allocated by the user prior to the call of this method.       
         */
        void coef_interpolant_GL(double (*f)(double), double* cf) const ;  

        /** Coefficients of the expansion over the basis polynomials of the 
         * orthogonal projection on the space of polynomials of maximum 
         * degree \e N. 
         *
         *  @param (*f)(double) function to be projected      
         *  @param cf [output] array containing the \e N+1 coefficients of the 
         *   expansion over the basis polynomials of the 
         *      orthogonal projection on the space of polynomials of maximum 
         *      degree \e N.          
         *          This array must be of size at least \e N+1 and must have
         *          been allocated by the user prior to the call of this method.       
         */
        virtual void coef_projection(double (*f)(double), double* cf) const = 0 ;  

        /** Value of some expansion over the polynomials:
         * \f$\sum_{i=0}^N a_i\, p_i(x)\f$
         *
         *  @param a array of size at least \e N+1 and storing the \e N+1 
         *          coefficients \f$a_i\f$ of the expansion on the basis 
         *          polynomials \f$p_i(x)\f$.
         *          For an interpolant polynomial, \e a is the array \c cf
         *          computed by the method
         *          \c coef_interpolant_Gauss (Gauss nodes) or
         *          \c coef_interpolant_GL (Gauss-Lobatto nodes) 
         *  @param x point in [-1,1] where the series is to be evaluated
         *  @return value of \f$\sum_{i=0}^N a_i\, p_i(x)\f$
         */
        double series(const double* a, double x) const ;  

};
    
            //----------------------------------//
            //      class Chebyshev_poly        //
            //----------------------------------//

/** Chebyshev polynomials. \ingroup (poly)
 * Each instance of this class represents a family of Chebyshev polynomials 
 * \f$(T_i)\f$ up to a given degree \e N.
 *
 */
class Chebyshev_poly : public Ortho_poly {
        
    // Constructors - Destructor
    // -------------------------
    public:
        /** Standard constructor
         *  @param ni maximum degree of the polynomials 
         */
        Chebyshev_poly(int ni ) ; 
        
	    Chebyshev_poly(const Chebyshev_poly& ) ; ///< Copy constructor
        
    public:
	    virtual ~Chebyshev_poly() ;			///< Destructor
 
    // Mutators / assignment
    // ---------------------
    public:
	    /// Assignment to another Chebyshev_poly
	    void operator=(const Chebyshev_poly&) ;	
	
            
    // Computation
    // -----------
    public: 
        /** Weight function \f$w\f$.
         *
         * @param x point in [-1,1] where the weight function is to be evaluated
         * @return \f$w(x)=1/\sqrt{1-x^2}\f$
         */
        virtual double weight(double x) const ; 

        /** Value of the polynomial \f$p_i\f$
         *
         * @param i no. of the polynomial
         * @param x point in [-1,1] where the polynomial is to be evaluated
         * @return value of  \f$p_i(x)\f$
         *
         */
        virtual double operator()(int i, double x) const ; 

        /// Gauss nodes
        virtual const Grid& gauss_nodes() const ;
        
        /** Gauss weights
         * 
         * @param i index of the weight
         * @return weight \f$w_i\f$
         */
        virtual double gauss_weight(int i) const ;
        
        /** Gamma factor \f$\gamma_i\f$ (square of the polynomial no. \e i with
         *  respect to the discrete scalar product associated with the
         *   Gauss nodes)
         */
        virtual double gauss_gamma(int i) const ;

        /// Gauss-Lobatto nodes
        virtual const Grid& gauss_lobatto_nodes() const ;
        
        /** Gauss-Lobatto weights
         * 
         * @param i index of the weight
         * @return weight \f$w_i\f$
         */
        virtual double gauss_lobatto_weight(int i) const ;
        
        /** Gamma factor \f$\gamma_i\f$ (square of the polynomial no. \e i with
         *  respect to the discrete scalar product associated with the
         *   Gauss-Lobatto nodes)
         */
        virtual double gauss_lobatto_gamma(int i) const ;

        /** Coefficients of the interpolant polynomial through the 
         *  Gauss-Lobatto nodes via a FFT.
         *
         *  @param (*f)(double) function to be interpolated         
         *  @param cf [output] array of size containing the \e N+1 coefficients
         *          of the expansion of the interpolant on the Chebyshev 
         *          polynomials.
         *          This array must be of size at least \e N+1 and must have
         *          allocated by the user prior to the call of this method.       
         */
        void coef_interpolant_GL_FFT(double (*f)(double), double* cf) const ;  

        /** Coefficients of the expansion over the Chebyshev polynomials of the 
         * orthogonal projection on the space of polynomials of maximum 
         * degree \e N. 
         *
         *  @param (*f)(double) function to be projected      
         *  @param cf [output] array containing the \e N+1 coefficients of the 
         *      expansion over the Chebyshev polynomials of the 
         *      orthogonal projection on the space of polynomials of maximum 
         *      degree \e N.          
         *          This array must be of size at least \e N+1 and must have
         *          been allocated by the user prior to the call of this method.       
         */
        virtual void coef_projection(double (*f)(double), double* cf) const  ;  

};
    
/// Display of an object of class \c Chebyshev_poly
ostream& operator<<(ostream& , const Chebyshev_poly& ) ;	


            //----------------------------------//
            //      class Legendre_poly        //
            //----------------------------------//

/** Legendre polynomials. \ingroup (poly)
 * Each instance of this class represents a family of Legendre polynomials 
 * \f$(T_i)\f$ up to a given degree \e N.
 *
 */
class Legendre_poly : public Ortho_poly {
        
    // Constructors - Destructor
    // -------------------------
    public:
        /** Standard constructor
         *  @param ni maximum degree of the polynomials 
         */
        Legendre_poly(int ni ) ; 
        
	    Legendre_poly(const Legendre_poly& ) ; ///< Copy constructor
        
    public:
	    virtual ~Legendre_poly() ;			///< Destructor
 
    // Mutators / assignment
    // ---------------------
    public:
	    /// Assignment to another Legendre_poly
	    void operator=(const Legendre_poly&) ;	
	
            
    // Computation
    // -----------
    public: 
        /** Weight function \f$w\f$.
         *
         * @param x point in [-1,1] where the weight function is to be evaluated
         * @return \f$w(x)=1\f$ (= constant value 1 for Legendre polynomials)
         */
        virtual double weight(double x) const ; 

        /** Value of the polynomial \f$p_i\f$
         *
         * @param i no. of the polynomial
         * @param x point in [-1,1] where the polynomial is to be evaluated
         * @return value of  \f$p_i(x)\f$
         *
         */
        virtual double operator()(int i, double x) const ; 

        /// Gauss nodes
        virtual const Grid& gauss_nodes() const ;
        
        /** Gauss weights
         * 
         * @param i index of the weight
         * @return weight \f$w_i\f$
         */
        virtual double gauss_weight(int i) const ;
        
        /** Gamma factor \f$\gamma_i\f$ (square of the polynomial no. \e i with
         *  respect to the discrete scalar product associated with the
         *   Gauss nodes)
         */
        virtual double gauss_gamma(int i) const ;

        /// Gauss-Lobatto nodes
        virtual const Grid& gauss_lobatto_nodes() const ;
        
        /** Gauss-Lobatto weights
         * 
         * @param i index of the weight
         * @return weight \f$w_i\f$
         */
        virtual double gauss_lobatto_weight(int i) const ;
        
        /** Gamma factor \f$\gamma_i\f$ (square of the polynomial no. \e i with
         *  respect to the discrete scalar product associated with the
         *   Gauss-Lobatto nodes)
         */
        virtual double gauss_lobatto_gamma(int i) const ;

        /** Coefficients of the expansion over the Legendre polynomials of the 
         * orthogonal projection on the space of polynomials of maximum 
         * degree \e N. 
         *
         *  @param (*f)(double) function to be projected      
         *  @param cf [output] array containing the \e N+1 coefficients of the 
         *      expansion over the Legendre polynomials of the 
         *      orthogonal projection on the space of polynomials of maximum 
         *      degree \e N.          
         *          This array must be of size at least \e N+1 and must have
         *          been allocated by the user prior to the call of this method.       
         */
        virtual void coef_projection(double (*f)(double), double* cf) const  ;  

};
    
/// Display of an object of class \c Legendre_poly
ostream& operator<<(ostream& , const Legendre_poly& ) ;	


#endif

