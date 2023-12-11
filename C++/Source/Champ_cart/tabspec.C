/*
 * Methods of class TabSpec.
 *
 *  (see file tabspec.h for documentation)
 *
 */

/*
 *   Copyright (c) 2019 Jerome Novak
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


#include<exception>
#include<fstream>
#include<iomanip>
#include "tabspec.h"

namespace Lorene {

  TabSpec::TabSpec(int dim1, int dim2, int dim3) : sizex(dim1), sizey(dim2),
						   sizez(dim3), tableau(0x0) {
    if ( (dim1<=0) || (dim2<=0) || (dim3<=0) ) {
      throw(std::invalid_argument("Non-positive size in TabSpec constructor")) ;
    }
    int nelts = dim1*dim2*dim3 ;
    tableau = new double[nelts] ;
    //#pragma omp parallel for
    for (int i=0; i<nelts; i++)
      tableau[i] = 0 ;
  }
  
  TabSpec::TabSpec(const string& nom_fich) {
    ifstream fich ;
    fich.exceptions ( ifstream::badbit );
    try {
      fich.open(nom_fich.c_str()) ;
      if (!fich) throw std::ios::failure( "Error opening file!" ) ;
    }
    catch (const ifstream::failure& e) {
      cerr << "Error opening " << nom_fich ;
      cerr << " in TabSpec constructor..." << endl ;
      throw ;
    }
    fich >> sizex >> sizey >> sizez ;
    if ( (sizex<=0) || (sizey<=0) || (sizez<=0) ) {
      throw(std::invalid_argument("Non-positive size in TabSpec constructor")) ;
    }
    int t_tot = sizex*sizey*sizez ;
    tableau = new double[t_tot] ;
    try {
      for (int i=0; i<t_tot; i++) {
	fich >> tableau[i] ;
      }
    }
    catch (const ifstream::failure& e) {
      cerr << "Error reading " << nom_fich << endl ;
      cerr << "in TabSpec constructor" << endl ;
      throw ;
    }  
  }

  TabSpec::TabSpec(const TabSpec& tab_in) : sizex(tab_in.sizex), sizey(tab_in.sizey),
					    sizez(tab_in.sizez), tableau(0x0) {

    int t_tot = sizex*sizey*sizez ;
    tableau = new double[t_tot] ;
    //#pragma omp parallel for
    for (int i=0; i<t_tot; i++) 
      tableau[i] = tab_in.tableau[i] ;
  }

  TabSpec::~TabSpec() {
    delete [] tableau ;
  }

  double TabSpec::operator()(int i, int j, int k) const {

#ifndef NDEBUG
    if ( (i<0) || (i>=sizex) || (j<0) || (j>= sizey) || (k<0) || (k>=sizez) )
      throw(std::out_of_range("Invalid access to TabSpec")) ;
#endif
    return tableau[(i*sizey + j)*sizez + k] ;

  }

  double& TabSpec::set(int i, int j, int k) {

#ifndef NDEBUG
    if ( (i<0) || (i>=sizex) || (j<0) || (j>= sizey) || (k<0) || (k>=sizez) )
      throw(std::out_of_range("Invalid access to TabSpec")) ;
#endif
    return tableau[(i*sizey + j)*sizez + k] ;


  }; 



  bool TabSpec::check_sizes(const TabSpec& t) const {
    bool resu = true ;
    resu = resu && (sizex == t.sizex) ;
    resu = resu && (sizey == t.sizey) ;
    resu = resu && (sizez == t.sizez) ;

    return resu;
  }

  void TabSpec::resize(int newx, int newy, int newz) {

    if ( (newx == sizex) && (newy == sizey) && (newz == sizez) )
      return ;
    else {
      delete [] tableau ;
      if ( (newx<=0) || (newy<=0) || (newz<=0) ) {
	throw(std::invalid_argument("Non-positive size in TabSpec::resize")) ;
      }
      sizex = newx ;
      sizey = newy ;
      sizez = newz ;
      int nelts = newx*newy*newz ;
      tableau = new double[nelts] ;
      //#pragma omp parallel for
      for (int i=0; i<nelts; i++)
	tableau[i] = 0. ;
    }
  }

  int TabSpec::get_dim() const {
    int dim = 3 ;
    if (sizez==1) dim-- ;
    if (sizey==1) dim-- ;
    if (sizex==1) dim-- ;
    return dim ;
  }

  void TabSpec::operator=(const TabSpec& tab_in) {
    if (!check_sizes(tab_in)) throw(range_error("Different array sizes")) ;

    int t_tot = sizex*sizey*sizez ;
    //#pragma omp parallel for
    for (int i=0; i<t_tot; i++)
      tableau[i] = tab_in.tableau[i] ;
  }
 
  void TabSpec::operator=(double x) {

    int t_tot = sizex*sizey*sizez ;
    //#pragma omp parallel for
    for (int i=0; i<t_tot; i++)
      tableau[i] = x ;
  }
 
  void TabSpec::write_file(const string& nom_fich) const {
  
    ofstream fich(nom_fich.c_str()) ;
    int t_tot = sizex*sizey*sizez ;
    fich << sizex << '\t' << sizey << '\t' << sizez << endl ;
    fich << setprecision(16) ;
    for (int i=0; i<t_tot; i++) 
      fich << tableau[i] << '\t' ;
    fich << endl ;
  }

  void TabSpec::display(ostream& ost) const {

    int dimens = get_dim() ;
    ost << "TabSpec with " << dimens << " dimension" << (dimens>1 ? "s:\n" : ":\n") ;
    if (sizex > 1) {
      ost << sizex ;
      if (dimens>1) {
	ost << "x" ;
	dimens-- ;
      }
    }
    if (sizey>1) {
      ost << sizey ;
      if (dimens>1) ost << "x" ;
    }
    if (sizez>1) ost << sizez ;
    ost << " elements" << endl ;
    ost << setprecision(5) ;

    if (sizez == 1) {
      for (int i=0; i<sizex; i++) {
	for (int j=0; j<sizey; j++) {
	  ost << tableau[i*sizey + j] << '\t' ;
	}
	if (sizey >1) ost << endl ;
      }
      ost << endl ;
    }
    else {
      for (int k=0; k<sizez; k++) {
	ost << "k=" << k << '\n' ;
	for (int j=0; j<sizey; j++) {
	  for (int i=0; i<sizex; i++) {
	    ost << tableau[(i*sizey+j)*sizez + k] << '\t' ;
	  }
	  ost << endl ;
	}
	ost << endl ;
      }
      ost << endl ;
    }
    return ;
  }
  ostream& operator<<(ostream& ost, const TabSpec& tab_in ) {
    tab_in.display(ost) ;
    return ost ;
  }

}
