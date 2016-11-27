/*
 *  Definition of special numbers (__infinity)
 *
 */

/*
 *   Copyright (c) 1999-2000 Jean-Alain Marck
 *   Copyright (c) 1999-2003 Eric Gourgoulhon
 *   Copyright (c) 2001 Joachim Frieben
 *   Copyright (c) 2003 Jerome Novak
 *   Copyright (c) 2003 Christian Klein
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


#ifndef __NBR_SPX_
#define __NBR_SPX_

/*
 * Definition des nombres particuliers (norme ieee)
 * 
 * Actuellement:
 *   infini code comme __infinity
 * 
 * On pourra mettre par la suite, au fur et a mesure des besoins:
 * 
 * extern double __libm_qnan_d ;
 * extern double __libm_inf_d ;
 * extern double __libm_neginf_d ;
 * extern double __infinity ;
 * 
 * Les machines reconnuees sont:
 *   SGI (par defaut)
 *   DEC Alpha (__alpha)
 *   HP-PA (__hppa)
 *   PC Linux  (__linux)
 * 
 * Pour info:
 * 
 *  --------------------------------------------------------
 *  |                                                      |
 *  | 		Double Precision Format                    |
 *  |                                                      |
 *  |  | bit 0 | bits 1 - 11 | bits 12 - 63                |
 *  |  -----------------------------------------------     |
 *  |  | Sign  | Exponent    | Mantissa              |     |
 *  |          |             |                       |     |
 *  |                                                      |
 *  |           \ msb   lsb / \ msb             lsb /      |
 *  |                                                      |
 *  --------------------------------------------------------
 * 
 *  ----------------------------------------------
 *  | Condition: plus or minus ZERO              |
 *  |                                            |
 *  | Sign    |  Exponent     |  Mantissa        |
 *  |         |  (All Zeros)  |  (All Zeros)     |
 *  |                                            |
 *  | Condition: plus or minus INFINITY          |
 *  |                                            |
 *  | Sign    |  Exponent     |  Mantissa        |
 *  |         |  (All Ones)   |  (All Zeros)     |
 *  |                                            |
 *  | Condition: NaN (Not a Number)              |
 *  |                                            |
 *  | Sign    |  Exponent     |  Mantissa        |
 *  |         |  (All Ones)   |  (Not All Zeros) |
 *  |                                            |
 *  ----------------------------------------------
 */


/*
 * $Id$
 * $Log$
 * Revision 1.7  2016/11/27 18:00:17  j_novak
 * Simplifiaction du nbr_spx.h, 'INFINITY' etant standard
 *
 * Revision 1.6  2014/10/06 15:09:40  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.5  2009/11/22 17:06:11  j_novak
 * Patched to build on a Macbook Pro (thanks E. Schnetter)
 *
 * Revision 1.4  2007/06/16 22:19:23  m_saijo
 * Add the case __i386__ for supporting Intel Mac.
 *
 * Revision 1.3  2003/02/05 13:51:40  e_gourgoulhon
 * Added the case __ppc__ for MacOS X.
 *
 * Revision 1.2  2002/09/09 12:57:22  e_gourgoulhon
 * Added the case of IBM AIX with xlC compiler.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.2  2001/05/29  15:59:27  eric
 * Ajout du cas des machines HP-PA (contribution de Joachim).
 *
 * Revision 2.1  2000/03/17  08:17:35  eric
 * Mise en conformite Linux : huge_val.h ---> math.h
 *
 * Revision 2.0  1999/02/15  10:41:51  hyc
 * *** empty log message ***
 *
 * Revision 2.1  1999/02/15  09:59:50  hyc
 * *** empty log message ***
 *
 * Revision 1.3  1998/06/10  07:25:37  eric
 * Ajout du cas PC Linux.
 *
 * Revision 1.2  1997/10/21  12:48:04  hyc
 * *** empty log message ***
 *
 * Revision 1.1  1997/10/17 15:40:01  hyc
 * Initial revision
 *
 *
 * $Header$
 *
 */

#include <cmath>
#define __infinity INFINITY

#endif
