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
 * Revision 1.2  2019/12/02 14:51:37  j_novak
 * Moved piecewise parabolic computation of dpdnb to a separate function.
 *
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
#include "tbl.h"


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

  // Computes the first derivative of the 2nd Tbl, with respect to the 1st one
  // using piecewise parabolic scheme.
  void compute_derivative(const Tbl& xx, const Tbl& yy, Tbl& dydx)
  {
    double y0, y1, y2, x0, x1, x2, der; 
    int np = xx.get_dim(0) ;
    assert ((yy.get_dim(0) == np) && (dydx.get_dim(0) == np)) ;
    
    // special case: i=0
    
    y0 = yy(0);
    y1 = yy(1);
    y2 = yy(2);
    
    x0 = xx(0);
    x1 = xx(1);
    x2 = xx(2);
    
    der = y0*(2*x0-x1-x2)/(x0-x1)/(x0-x2) +
      y1*(x0-x2)/(x1-x0)/(x1-x2) +
      y2*(x0-x1)/(x2-x0)/(x2-x1) ;
    
    dydx.set(0) = der ; 
    
    for(int i=1;i<np-1;i++) { 
      
      y0 = yy(i-1);
      y1 = yy(i);
      y2 = yy(i+1);
      
      x0 = xx(i-1);
      x1 = xx(i);
      x2 = xx(i+1);
      
      der = y0*(x1-x2)/(x0-x1)/(x0-x2) +
	y1*(2*x1-x0-x2)/(x1-x0)/(x1-x2) +
	y2*(x1-x0)/(x2-x0)/(x2-x1) ;
      
      dydx.set(i) = der ;
      
    } 
    
    // special case: i=np-1
    
    y0 = yy(np-3);
    y1 = yy(np-2);
    y2 = yy(np-1);
    
    x0 = xx(np-3);
    x1 = xx(np-2);
    x2 = xx(np-1);
    
    der = y0*(x2-x1)/(x0-x1)/(x0-x2) +
      y1*(x2-x0)/(x1-x0)/(x1-x2) +
      y2*(2*x2-x0-x1)/(x2-x0)/(x2-x1) ;
    
    dydx.set(np-1) = der ;

    return ;
  }
  
} // End of namespace Lorene
