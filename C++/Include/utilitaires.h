/*
 *  Prototypes of various utilities for Lorene
 *
 */

/*
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
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


#ifndef	__UTILITAIRES_H_
#define	__UTILITAIRES_H_


/*
 * $Id$
 * $Log$
 * Revision 1.8  2003/12/17 23:12:30  r_prix
 * replaced use of C++ <string> by standard ANSI char* to be backwards compatible
 * with broken compilers like MIPSpro Compiler 7.2 on SGI Origin200. ;-)
 *
 * Revision 1.7  2003/12/05 15:05:53  r_prix
 * added read_variable() for (C++-type) strings.
 *
 * Revision 1.6  2003/12/04 12:33:21  r_prix
 * added prototypes and documentation for variable-reading functions
 * (read_variable, load_file, load_file_buffered)
 *
 * Revision 1.5  2002/05/02 15:16:22  j_novak
 * Added functions for more general bi-fluid EOS
 *
 * Revision 1.4  2002/04/16 08:06:44  j_novak
 * Addition of zerosec_borne
 *
 * Revision 1.3  2002/04/11 09:19:46  j_novak
 * Back to old version of zerosec
 *
 * Revision 1.2  2001/12/04 21:24:33  e_gourgoulhon
 *
 * New functions fwrite_be and fread_be for writing/reading in a
 * binary file according to the big endian convention.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 1.6  2001/09/14  14:23:53  eric
 * Ajout de la fonction zero_list.
 *
 * Revision 1.5  2001/05/29  16:11:21  eric
 * Modif commentaires (mise en conformite Doc++ 3.4.7).
 *
 * Revision 1.4  1999/12/24  12:59:54  eric
 * Ajout de la routine zero_premier.
 *
 * Revision 1.3  1999/12/15  15:39:42  eric
 * *** empty log message ***
 *
 * Revision 1.2  1999/12/15  15:17:03  eric
 * *** empty log message ***
 *
 * Revision 1.1  1999/12/15  09:41:47  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */
 
#include "stdio.h"
#include <string.h>

class Param ;
class Tbl ;

    /** @name Miscellaneous.
     */
    //@{     

/** Setting a stop point in a code.
 * 
 *  Stops the execution of a code, until the 'Enter' case is hit.
 *  @param   a	[input] stops the run if, and only if, a=0.
 *			Default value : 0.
 * 
 * 
 */
void arrete(int a = 0) ;

/** Locates the sub-interval containing the first zero of a function in 
 *  a given interval.
 * 
 *  @param (*f)(double, const Param\&) [input] Function the zero of which is 
 *		    to be searched: 
 *			f(x0, par) = 0 , where par are the parameters of the
 *		    function f, stored in an object of the Lorene class 
 *		    {\tt Param}. 
 *  @param a [input] Lower bound of the search interval [a, b]
 *  @param b [input] Higher bound of the search interval [a, b]
 *  @param n [input] Number of subdivisions of the interval [a, b]
 *  @param a0 [output] Lower bound of the first (i.e. closest to a) interval 
 *		      [a0, b0] which contains a zero of f 
 *  @param b0 [output] Higher bound of the first (i.e. closest to a) interval 
 *		      [a0, b0] which contains a zero of f 
 *  @return  true if the interval [a0, b0] containing a zero of f has been
 *	    found, false otherwise
 */
bool zero_premier(double (*f)(double, const Param&), const Param& par,
		  double a, double b, int n, double& a0, double& b0) ;



/** Finding the zero a function.
 * 
 *  This routine locates a zero by means of the secant method. 
 * 
 *  @param (*f)(double, const Param\&) [input] Function the zero of which is 
 *		    to be searched: the routine computes x0 in a given
 *		    interval [a, b] such that 
 *			f(x0, par) = 0 , where par are the parameters of the
 *		    function f, stored in an object of the Lorene class 
 *		    {\tt Param}. 
 *  @param par [input] Parameters of the function f.
 *  @param a [input] Lower bound of the search interval [a, b]
 *  @param b [input] Higher bound of the search interval [a, b]
 *  @param precis [input] Required precision in the determination of x0 : 
 *			the returned solution will be x0 +/- precis
 *  @param nitermax [input] Maximum number of iterations in the secant 
 *			    method to compute x0.
 *  @param niter [output] Number of iterations effectively used in computing x0				
 *  @return x0 (zero of function f)
 *
 */
double zerosec( double (*f)(double, const Param&), const Param& par, 
		double a, double b, double precis, int nitermax, 
		int& niter) ;

/** Finding the zero a function on a bounded domain.
 * 
 *  Same as {\tt zerosec} but insures that all calls to the function f
 *  are within [x1,x2]. Namely, it requires that f(x1)*f(x2)<0. 
 * 
 *  @param (*f)(double, const Param\&) [input] Function the zero of which is 
 *		    to be searched: the routine computes x0 in a given
 *		    interval [a, b] such that 
 *			f(x0, par) = 0 , where par are the parameters of the
 *		    function f, stored in an object of the Lorene class 
 *		    {\tt Param}. 
 *  @param par [input] Parameters of the function f.
 *  @param a [input] Lower bound of the search interval [a, b]
 *  @param b [input] Higher bound of the search interval [a, b]
 *  @param precis [input] Required precision in the determination of x0 : 
 *			the returned solution will be x0 +/- precis
 *  @param nitermax [input] Maximum number of iterations in the secant 
 *			    method to compute x0.
 *  @param niter [output] Number of iterations effectively used in computing x0				
 *  @return x0 (zero of function f)
 *
 */
double zerosec_b( double (*f)(double, const Param&),const Param& par, 
		  double a, double b, double precis, int nitermax, 
		  int& niter) ;

/** Locates approximatively all the zeros of a function in a given interval.
 *  The N zeros are located in N intervals [az(i), bz(i)] with
 *   $0\leq i \leq N-1$.  
 * 
 *  @param (*f)(double, const Param\&) [input] Function the zeros of which are 
 *		    to be located: a zero x0 is defined by 
 *			f(x0, par) = 0 , where par are the parameters of the
 *		    function f, stored in an object of the Lorene class 
 *		    {\tt Param}. 
 *  @param par [input] Parameters of the function f.
 *  @param xmin [input] Lower bound of the search interval 
 *  @param xmax [input] Higher bound of the search interval 
 *  @param nsub [input] Number of subdivision of the interval [xmin, xmax]
 *		        to locate the zeros
 *  @param az [output] 1-D array (Lorene {\tt Tbl}) contain the lower bounds
 *			of the intervals containing a zero. This {\tt Tbl}\ 
 *			is allocated by the routine via a {\tt new Tbl}\ 
 *			command (hence the pointer type). 
 *  @param bz [output] 1-D array (Lorene {\tt Tbl}) contain the higher bounds
 *			of the intervals containing a zero. This {\tt Tbl}\ 
 *			is allocated by the routine via a {\tt new Tbl}\ 
 *			command (hence the pointer type). 
 *  
 */
void zero_list( double (*f)(double, const Param&), const Param& par,
		double xmin, double xmax, int nsub, 
		Tbl*& az, Tbl*& bz ) ;  
		
/** Writes integer(s) into a binary file according to the
 *  big endian convention.
 *
 *  This function has the same prototype and return value than
 *  the {\tt fwrite} function of the {\tt stdio} C library.
 *  The difference is that it ensures that the binary file is
 *  written in the big endian format, whatever the system is
 *  using little endian or big endian.
 *	@param aa [input] integer array to be written (in case of one
 *		element, address of this integer)
 *	@param size [input] number of bytes of one {\tt int} (must
 *		be 4)
 *	@param nb [input] number of elements in the array {\tt aa}
 *	@param fich [input] binary file (must have been
 *		open by {\tt fopen})
 *	@return number of integers effectively written in the file
 */		
int fwrite_be(const int* aa, int size, int nb, FILE* fich) ;

/** Writes double precision number(s) into a binary file according to the
 *  big endian convention.
 *
 *  This function has the same prototype and return value than
 *  the {\tt fwrite} function of the {\tt stdio} C library.
 *  The difference is that it ensures that the binary file is
 *  written in the big endian format, whatever the system is
 *  using little endian or big endian.
 *	@param aa [input] array of {\tt double} to be written (in case of one
 *		element, address of this {\tt double})
 *	@param size [input] number of bytes of one {\tt double} (must
 *		be 8)
 *	@param nb [input] number of elements in the array {\tt aa}
 *	@param fich [input] binary file (must have been
 *		open by {\tt fopen})
 *	@return number of {\tt double} effectively written in the file
 */		
int fwrite_be(const double* aa, int size, int nb, FILE* fich) ;

/** Reads integer(s) from a binary file according to the
 *  big endian convention.
 *
 *  This function has the same prototype and return value than
 *  the {\tt fread} function of the {\tt stdio} C library.
 *  The difference is that it assumes that the binary file is
 *  written in the big endian format, whatever the system is
 *  using little endian or big endian.
 *	@param aa [output] integer array to be read (in case of one
 *		element, address of this integer)
 *	@param size [input] number of bytes of one {\tt int} (must
 *		be 4)
 *	@param nb [input] number of elements in the array {\tt aa}
 *	@param fich [input] binary file (must have been
 *		open by {\tt fopen})
 *	@return number of integers effectively read in the file
 */		
int fread_be(int* aa, int size, int nb, FILE* fich) ;

/** Reads double precision number(s) from a binary file according to the
 *  big endian convention.
 *
 *  This function has the same prototype and return value than
 *  the {\tt fread} function of the {\tt stdio} C library.
 *  The difference is that it assumes that the binary file is
 *  written in the big endian format, whatever the system is
 *  using little endian or big endian.
 *	@param aa [output] array of {\tt double} to be read (in case of one
 *		element, address of this {\tt double})
 *	@param size [input] number of bytes of one {\tt double} (must
 *		be 8)
 *	@param nb [input] number of elements in the array {\tt aa}
 *	@param fich [input] binary file (must have been
 *		open by {\tt fopen})
 *	@return number of {\tt double} effectively read in the file
 */		
int fread_be(double* aa, int size, int nb, FILE* fich) ;


/** Read file into memory and returns pointer to data.
 * 	@return: pointer to allocated memory or NULL on error.\\
 * 
 *  NOTE: don't forget to free the memory after use !
 *
 */
char *load_file(char *fname);

/** Returns pointer to data from a file using a buffer.
 *  This function only reads from disk if the file has not been buffered yet.
 *  If a new file is read, the buffer is free'ed and the new data allocated.
 *	@param fname [input] name of file to read in. You can use NULL for previous file.
 *  	@return : pointer to buffered data or NULL on error.\\
 *
 *   NOTE: do _NEVER_ free the (buffer-)data pointer, or the next call will crash!!
 */
char *load_file_buffered(char *fname);

/** Reads a variable from file. 
 * Variable definitions can be of the type 
 * "variable = value", all other lines are ignored as comments.
 * (alternatively, you can use ":" or whitespace instead of "=")
 * 
 * NOTE: the variable-definition has to be at the beginning of a 
 * line (modulo whitespace), or it will be considered a comment!
 *
 *	@param fname [input] Name of config-file to read from. Use NULL to use buffered file.
 *	@param var_name [input] Variable-name to read from config-file
 *	@param fmt [input]  C-style format string for reading (see {\tt sscanf}).
 *	@param varp [output] Pointer to C-variable to read value into (has to be big enough!!)
 *
 *	@return 0 on success, -1 on error. \\
 *
 *   NOTE: rather use one of the type-specific overloaded functions below whenever possible
 *   (safer due to type-checking)
 */
int read_variable(char *fname, char *var_name, char *fmt, void *varp);

/// Read an integer-variable from file (cf {\tt read_variable(char *, char *, char *, void *)}).
int read_variable(char *fname, char *var_name, int &var);
/// Read a bool variable from file (cf {\tt read_variable(char *, char *, char *, void *)}).
int read_variable(char *fname, char *var_name, bool &var);
/// Read a double variable from file (cf {\tt read_variable(char *, char *, char *, void *)}).
int read_variable(char *fname, char *var_name, double &var);
/// Read a (ANSI C) string variable from file.
int read_variable (char *fname, char *var_name, char *str);

/// 'Improved' malloc that sets memory to 0 and also auto-terminates on error.
void *MyMalloc (long bytes);

/// A portable routine to determine the length of a file
int FS_filelength (FILE *f);


    
    //@}
    
#endif
