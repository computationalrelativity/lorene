 /*
 *  Definition of Lorene class TabSpec
 *
 */

/*
 *   Copyright (c) 2019 Jerome Novak
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

#ifndef __TABSPEC_H_ 
#define __TABSPEC_H_ 

#include<string>
#include<iostream>
#include<fftw3.h>

using namespace std ;

namespace Lorene {

  /** 
   * 3-indices array to be used with the representation of functions
   * with Cartesian coordinates. \ingroup (map_cart)
   */
  class TabSpec {
    
    // Donnees : 
    // -------
  protected:
    int sizex ;
    int sizey ;
    int sizez ;
    double* tableau ;
    
    // Constructeurs - Destructeur
    // ---------------------------
  public:
    explicit TabSpec(int dim1, int dim2=1, int dim3=1) ; //constructeur standard
    explicit TabSpec(const string& ) ; // constructeur depuis un fichier
    TabSpec(const TabSpec&) ; //constructeur par recopie
    
    virtual ~TabSpec() ; //destructeur
    
    // Affectation
    // -----------
    void operator=(const TabSpec&) ; //affectation depuis un autre TabSpec
    void operator=(double) ; //affectation depuis un double
    
    // Acces aux donnees
    // -----------------
    int get_sizex() const {return sizex ; };
    int get_sizey() const {return sizey ; };
    int get_sizez() const {return sizez ; };
    int get_nelt() const {return sizex*sizey*sizez ; } ;
    int get_dim() const ;
    
    double operator()(int i, int j=0, int k=0) const ;
    
    double& set(int i, int j=0, int k=0) ; 
    
    // Sauvegarde dans un fichier
    // --------------------------
    void write_file(const string&) const ; 
    
  protected:
    // Verif des tailles 
    bool check_sizes(const TabSpec&) const ;
    void resize(int, int=1, int=1) ;
    
    // Affichage
    virtual void display(ostream& ) const ; //pour utiliser le polymorphisme

  public:  //"public" inutile ici...
    friend ostream& operator<<(ostream&, const TabSpec& ) ; //la fonction externe
    //a appeler
    // operations mathematiques:
    friend TabSpec apply(const TabSpec& t, double (*p_fonc)(double));
    friend TabSpec apply(const TabSpec& t1, const TabSpec& t2,
			 double (*p_fonc)(double,double));
    friend TabSpec apply(const TabSpec& t, double r,
			 double (*p_fonc)(double,double));
    friend TabSpec apply(double r, const TabSpec& t,
			 double (*p_fonc)(double,double));
    friend double max(const TabSpec& t);
    friend TabSpec operator-(const TabSpec&) ;
    
    friend class FuncSpec ;
    friend fftw_plan prepare_fft(int, TabSpec*&) ;
    friend fftw_plan back_fft(int, TabSpec*&) ;
    
  };
  
  TabSpec operator+(const TabSpec&, const TabSpec&) ;
  TabSpec operator+(const TabSpec&, double) ;
  TabSpec operator+(double, const TabSpec&) ;
  TabSpec operator-(const TabSpec&, const TabSpec&) ;
  TabSpec operator-(const TabSpec&, double) ;
  TabSpec operator-(double, const TabSpec&) ;
  TabSpec operator-(const TabSpec&) ;
  TabSpec operator*(const TabSpec&, const TabSpec&) ;
  TabSpec operator*(const TabSpec&, double) ;
  TabSpec operator*(double, const TabSpec&) ;
  TabSpec operator/(const TabSpec&, const TabSpec&) ;
  TabSpec operator/(const TabSpec&, double) ;
  TabSpec operator/(double, const TabSpec&) ;
  
  TabSpec sin(const TabSpec&) ;
  TabSpec cos(const TabSpec&) ;
  TabSpec tan(const TabSpec&) ;
  TabSpec exp(const TabSpec&) ;
  TabSpec log(const TabSpec&) ;
  TabSpec sqrt(const TabSpec&) ;
  TabSpec pow(const TabSpec&, double) ;
  TabSpec abs(const TabSpec&) ;
  double max(const TabSpec&) ;

}
#endif
