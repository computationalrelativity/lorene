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


#include "indent.h"

int iostream_index_indent = ios::xalloc() ;

ostream& iendl(ostream& o) {
    o << '\n' ;
    int mem = o.width(o.iword(iostream_index_indent)*2) ;
    o << "" << flush ;
    o.width(mem) ;
    return o ;
}

static struct CInitStdStream {
    CInitStdStream() {
	cout << resetindent ;
	cerr << resetindent ;
	clog << resetindent ;
    }
} InitStdStream ;

