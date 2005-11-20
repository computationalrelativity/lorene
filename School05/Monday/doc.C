/*
 * $Id$
 * $Log$
 * Revision 1.2  2005/11/20 17:50:18  e_gourgoulhon
 * Better documentation.
 *
 * Revision 1.1  2005/11/14 01:56:58  e_gourgoulhon
 * First version
 *
 *
 * $Header$
 *
 */


/*!\mainpage POLYNOMIAL APPROXIMATION
 *
 * This directory contains simple C++ classes, which are not part of
 * the main <A HREF="http://www.lorene.obspm.fr/"> LORENE </A> library.
 *
 *
 * These classes are intended to illustrate polynomial 
 * approximation to functions (mainly interpolation on various
 * grids). There are two families of classes:
 *  \li \c Grid : grids as sampling of [-1,1], declared in the file \c grid.h
 *  \li \c Ortho_poly : orthogonal polynomials, declared in the file 
 *          \c ortho_poly.h 
 *
 * There are also some simple graphical routines (declared in the file
 *  \c plot.h), which provides provides an interface to PGPLOT, in order
 * to perform easily some drawings (X11, EPS or PNG). 
 * See Modules -> <A HREF="group__graphbasic.html"> Basic graphical routines </A>  
 *
 * The main codes are
 *  \li \c es.C : test code for class \c Grid
 *  \li \c runge.C : exhibition of Runge's phenomenon
 *  \li \c cheby.C : Chebyshev expansions
 *  \li \c legendre.C : Legendre expansions
 *  \li \c interpole.C : comparison of various interpolations of a function
 *  \li \c approx_ortho.C : approximation of a function by orthogonal 
 *                          polynomials
 *  
 */
 
/**
 * \defgroup grid Grids
 *
 * A \e grid is a set of points (called \e nodes) in [-1,1].
 *
 */
 
/**
 * \defgroup poly Orthogonal polynomials
 *
 *
 */
 
 
