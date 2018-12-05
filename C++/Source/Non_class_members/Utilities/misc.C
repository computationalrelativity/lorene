/*
 *  Various utilities...
 *
 */

/*
 *   Copyright (c) 2018  Jerome Novak
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

char misc_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2018/12/05 15:03:20  j_novak
 * New Mg3d constructor from a formatted file.
 *
 * Revision 1.4  2003/10/19 20:01:10  e_gourgoulhon
 * Template file
 *
 * $Header$
 *
 */

// Lorene headers
#include "headcpp.h"


namespace Lorene {

  // Searches the file 'infile' for a given 'pattern'. The file stream is
  // positionned just after the occurrence of the pattern.
  //=====================================================================
  bool search_file(ifstream& infile, const string& pattern) {
    string line ;
    while( getline(infile, line) ) {
      if (line.find(pattern, 0) != string::npos)
	return true ;
    }
    return false ;
  }
  
} // End of namespace Lorene
