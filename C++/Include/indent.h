/*
 *  Indentation
 *
 */

/*
 *   Copyright (c) 1999-2000 Jean-Alain Marck
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


#ifndef __INDENT_H_
#define __INDENT_H_

/*
 * $Id$
 * $Log$
 * Revision 1.2  2002/10/16 14:36:29  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.0  1999/02/15  10:41:51  hyc
 * *** empty log message ***
 *
 * Revision 2.0  1998/12/03  13:48:02  hyc
 * Version 2
 *
 * Revision 2.0  1998/12/03  13:48:02  hyc
 * Version 2
 *
 *
 * $Header$
 *
 */


#include "headcpp.h"

extern int iostream_index_indent ;

ostream& iendl(ostream& ) ;

inline ostream& incindent(ostream& o) {
    ++o.iword(iostream_index_indent) ;
    return o ;
}
inline ostream& decindent(ostream& o) {
    if (o.iword(iostream_index_indent)) {
	--o.iword(iostream_index_indent) ;
    }
    return o ;
}
inline ostream& resetindent(ostream& o) {
    o.iword(iostream_index_indent) = 0 ;
    return o ;
}
inline ostream& incendl(ostream& o) {
    return o << incindent << iendl ;
}
inline ostream& decendl(ostream& o) {
    return o << decindent << iendl ;
}
inline ostream& operator<<(ostream& o, ostream& ) {
    return o ;
}

#endif
