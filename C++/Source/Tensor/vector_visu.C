/*
 *  3D visualization of a Vector via OpenDX 
 *
 *    (see file vector.h for documentation).
 *
 */

/*
 *   Copyright (c) 2003  Eric Gourgoulhon
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

 

/*
 * $Id$
 * $Log$
 * Revision 1.6  2016/12/05 16:18:18  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.5  2014/10/13 08:53:45  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.4  2014/10/06 15:13:21  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.3  2005/02/16 15:31:56  m_forot
 * Add the visu_streamline function
 *
 * Revision 1.2  2003/12/15 08:30:39  p_grandclement
 * Addition of #include <string.h>
 *
 * Revision 1.1  2003/12/14 21:48:26  e_gourgoulhon
 * First version
 *
 *
 * $Header$
 *
 */

// C headers
#include <cstdlib>
#include <cstring>

// Lorene headers
#include "tensor.h"

                    
namespace Lorene {
void Vector::visu_arrows(double xmin, double xmax, double ymin, double ymax,
    double zmin, double zmax, const char* title0, const char* filename0, 
    bool start_dx, int nx, int ny, int nz) const {

    const Vector* vect ; 
    Vector* vect_tmp = 0x0 ; 
    
    // The drawing is possible only in Cartesian coordinates:
    if (*triad == mp->get_bvect_cart()) {
        vect = this ; 
    }
    else {
        if (*triad == mp->get_bvect_spher()) {
            vect_tmp = new Vector(*this) ; 
            vect_tmp->change_triad( mp->get_bvect_cart() ) ; 
            vect = vect_tmp ;
        }
        else {
            cerr << "Vector::visu_arrows : unknown triad !" << endl ; 
            abort() ; 
        }
    }
    
    // The drawing is possible only if dzpuis = 0 
    bool dzpnonzero = false ; 
    for (int i=1; i<=3; i++) {
        dzpnonzero = dzpnonzero || !(operator()(i).check_dzpuis(0)) ; 
    }
    if (dzpnonzero) {
        if (vect_tmp == 0x0) {
            vect_tmp = new Vector(*this) ;                
        }
        for (int i=1; i<=3; i++) {
            Scalar& cvect = vect_tmp->set(i) ; 
            int dzpuis = cvect.get_dzpuis() ; 
            if (dzpuis != 0) {
                cvect.dec_dzpuis(dzpuis) ; 
            }
        }
        vect = vect_tmp ; 
    }
    
    char* title ;
    char* title_quotes ;
    if (title0 == 0x0) {
        title = new char[2] ; 
        strcpy(title, " ") ; 

        title_quotes = new char[4] ; 
        strcpy(title_quotes, "\" \"") ; 
    }
    else {
        title = new char[ strlen(title0)+1 ] ; 
        strcpy(title, title0) ;
         
        title_quotes = new char[ strlen(title0)+3 ] ; 
        strcpy(title_quotes, "\"") ; 
        strcat(title_quotes, title0) ; 
        strcat(title_quotes, "\"") ; 
    }
    
    // --------------------------------------------------------
    // Data file for OpenDX
    // --------------------------------------------------------

    char* filename ;
    if (filename0 == 0x0) {
        filename = new char[30] ; 
        strcpy(filename, "vector3d_arrows.dxdata") ; 
    }
    else {
        filename = new char[ strlen(filename0)+8 ] ; 
        strcpy(filename, filename0) ; 
        strcat(filename, ".dxdata") ; 
    }

    ofstream fdata(filename) ; // output file
    
    fdata << title << "\n" ; 
    fdata << "size : " << nx << " x " << ny << " x " << nz << "\n" ;
    fdata << "x_min = " << xmin << "  x_max = " << xmax << "\n" ; 
    fdata << "y_min = " << ymin << "  y_max = " << ymax << "\n" ; 
    fdata << "z_min = " << zmin << "  z_max = " << zmax << "\n" ; 
    
    // The spectral coefficients are required
    const Valeur& vax = (vect->operator()(1)).get_spectral_va() ; 
    const Valeur& vay = (vect->operator()(2)).get_spectral_va() ; 
    const Valeur& vaz = (vect->operator()(3)).get_spectral_va() ; 
    vax.coef() ; 
    vay.coef() ; 
    vaz.coef() ; 
    const Mtbl_cf& cvax = *(vax.c_cf) ; 
    const Mtbl_cf& cvay = *(vay.c_cf) ; 
    const Mtbl_cf& cvaz = *(vaz.c_cf) ; 
    
    // What follows assumes that the mapping is radial:
    assert( dynamic_cast<const Map_radial*>(mp) != 0x0 ) ; 
    
    fdata.precision(5) ; 
    fdata.setf(ios::scientific) ; 
    
    // Loop on the points of the drawing box
    // ---------------------------------------
    double dx = (xmax - xmin) / double(nx-1) ; 
    double dy = (ymax - ymin) / double(ny-1) ; 
    double dz = (zmax - zmin) / double(nz-1) ; 
    
    int npoint = 0 ;    // number of data points per line in the file

    for (int k=0; k<nz; k++) {

        double zz = zmin + dz * k ;

        for (int j=0; j<ny; j++) {
            
            double yy = ymin + dy * j ; 
            
            for (int i=0; i<nx; i++) {
                
                double xx = xmin + dx * i ; 
                    
                // Values of (r,theta,phi) corresponding to (xa,ya,za) :
                double rr, th, ph ;  // polar coordinates of the mapping associated
                                     // to *this
            
                mp->convert_absolute(xx, yy, zz, rr, th, ph) ; 

                // Values of (l,xi,theta',phi') corresponding to (r,theta,phi):
                double xi ; int l ; 
            
                mp->val_lx(rr, th, ph, l, xi) ;   // radial mapping assumption
            
                // Field value at this point:
            
                double vx = cvax.val_point(l, xi, th, ph) ;
                double vy = cvay.val_point(l, xi, th, ph) ;
                double vz = cvaz.val_point(l, xi, th, ph) ;

                fdata.width(14) ; fdata << vx ; 
                fdata.width(14) ; fdata << vy ; 
                fdata.width(14) ; fdata << vz ; 
                npoint++ ; 
             
                if (npoint == 3) {
                    fdata << "\n" ; 
                    npoint = 0 ; 
                }
            
            }
        }

    }
    
    if (npoint != 0) fdata << "\n" ; 

    fdata.close() ; 
    
    // --------------------------------------------------------
    // Header file for OpenDX
    // --------------------------------------------------------

    char* headername ;
    if (filename0 == 0x0) {
        headername = new char[30] ; 
        strcpy(headername, "vector3d_arrows.dxhead") ; 
    }
    else {
        headername = new char[ strlen(filename0)+9 ] ; 
        strcpy(headername, filename0) ; 
        strcat(headername, ".dxhead") ; 
    }

    ofstream fheader(headername) ;
    
    fheader << "file = " << filename << endl ; 
    fheader << "grid = " << nx << " x " << ny << " x " << nz << endl ; 
    fheader << "format = ascii" << endl ;  
    fheader << "interleaving = record-vector"  << endl ;
    fheader << "majority = column" << endl ; 
    fheader << "header = lines 5" << endl ; 
    fheader << "field = " << title_quotes << endl ; 
    fheader << "structure = 3-vector" << endl ; 
    fheader << "type = float" << endl ; 
    fheader << "dependency = positions" << endl ; 
    fheader << "positions = regular, regular, regular, " 
            << xmin << ", " << dx << ", " 
            << ymin << ", " << dy << ", " 
            << zmin << ", " << dz << endl ; 
    fheader << endl ; 
    fheader << "end" << endl ; 
    
    fheader.close() ; 
    

    if ( start_dx ) {       // Launch of OpenDX
        
        char* commande = new char[ strlen(headername) + 60 ] ;
        strcpy(commande, "ln -s ") ; 
        strcat(commande, headername) ; 
        strcat(commande, " visu_vector3d.dxhead") ; 
    
        system("rm -f visu_vector3d.dxhead") ; 
        system(commande) ;                      // ln -s headername visu_section.general
        system("dx -image visu_vector3d.net &") ; 
    
        delete [] commande ;    
    }


    // Final cleaning
    // --------------
    
    if (vect_tmp != 0x0) delete vect_tmp ;
    delete [] title ; 
    delete [] title_quotes ; 
    delete [] filename ; 
    delete [] headername ;    
   

}   

void Vector::visu_streamline(double xmin, double xmax, double ymin, double ymax,
			     double zmin, double zmax, const char* title0, const char* filename0, 
			     bool start_dx, int nx, int ny, int nz) const {

    const Vector* vect ; 
    Vector* vect_tmp = 0x0 ; 
    
    // The drawing is possible only in Cartesian coordinates:
    if (*triad == mp->get_bvect_cart()) {
        vect = this ; 
    }
    else {
        if (*triad == mp->get_bvect_spher()) {
            vect_tmp = new Vector(*this) ; 
            vect_tmp->change_triad( mp->get_bvect_cart() ) ; 
            vect = vect_tmp ;
        }
        else {
            cerr << "Vector::visu_streamline : unknown triad !" << endl ; 
            abort() ; 
        }
    }
    
    // The drawing is possible only if dzpuis = 0 
    bool dzpnonzero = false ; 
    for (int i=1; i<=3; i++) {
        dzpnonzero = dzpnonzero || !(operator()(i).check_dzpuis(0)) ; 
    }
    if (dzpnonzero) {
        if (vect_tmp == 0x0) {
            vect_tmp = new Vector(*this) ;                
        }
        for (int i=1; i<=3; i++) {
            Scalar& cvect = vect_tmp->set(i) ; 
            int dzpuis = cvect.get_dzpuis() ; 
            if (dzpuis != 0) {
                cvect.dec_dzpuis(dzpuis) ; 
            }
        }
        vect = vect_tmp ; 
    }
    
    char* title ;
    char* title_quotes ;
    if (title0 == 0x0) {
        title = new char[2] ; 
        strcpy(title, " ") ; 

        title_quotes = new char[4] ; 
        strcpy(title_quotes, "\" \"") ; 
    }
    else {
        title = new char[ strlen(title0)+1 ] ; 
        strcpy(title, title0) ;
         
        title_quotes = new char[ strlen(title0)+3 ] ; 
        strcpy(title_quotes, "\"") ; 
        strcat(title_quotes, title0) ; 
        strcat(title_quotes, "\"") ; 
    }
    
    // --------------------------------------------------------
    // Data file for OpenDX
    // --------------------------------------------------------

    char* filename ;
    if (filename0 == 0x0) {
        filename = new char[30] ; 
        strcpy(filename, "vector3d_streamline.dxdata") ; 
    }
    else {
        filename = new char[ strlen(filename0)+8 ] ; 
        strcpy(filename, filename0) ; 
        strcat(filename, ".dxdata") ; 
    }

    ofstream fdata(filename) ; // output file
    
    fdata << title << "\n" ; 
    fdata << "size : " << nx << " x " << ny << " x " << nz << "\n" ;
    fdata << "x_min = " << xmin << "  x_max = " << xmax << "\n" ; 
    fdata << "y_min = " << ymin << "  y_max = " << ymax << "\n" ; 
    fdata << "z_min = " << zmin << "  z_max = " << zmax << "\n" ; 
    
    // The spectral coefficients are required
    const Valeur& vax = (vect->operator()(1)).get_spectral_va() ; 
    const Valeur& vay = (vect->operator()(2)).get_spectral_va() ; 
    const Valeur& vaz = (vect->operator()(3)).get_spectral_va() ; 
    vax.coef() ; 
    vay.coef() ; 
    vaz.coef() ; 
    const Mtbl_cf& cvax = *(vax.c_cf) ; 
    const Mtbl_cf& cvay = *(vay.c_cf) ; 
    const Mtbl_cf& cvaz = *(vaz.c_cf) ; 
    
    // What follows assumes that the mapping is radial:
    assert( dynamic_cast<const Map_radial*>(mp) != 0x0 ) ; 
    
    fdata.precision(5) ; 
    fdata.setf(ios::scientific) ; 
    
    // Loop on the points of the drawing box
    // ---------------------------------------
    double dx = (xmax - xmin) / double(nx-1) ; 
    double dy = (ymax - ymin) / double(ny-1) ; 
    double dz = (zmax - zmin) / double(nz-1) ; 
    
    int npoint = 0 ;    // number of data points per line in the file

    for (int k=0; k<nz; k++) {

        double zz = zmin + dz * k ;

        for (int j=0; j<ny; j++) {
            
            double yy = ymin + dy * j ; 
            
            for (int i=0; i<nx; i++) {
                
                double xx = xmin + dx * i ; 
                    
                // Values of (r,theta,phi) corresponding to (xa,ya,za) :
                double rr, th, ph ;  // polar coordinates of the mapping associated
                                     // to *this
            
                mp->convert_absolute(xx, yy, zz, rr, th, ph) ; 

                // Values of (l,xi,theta',phi') corresponding to (r,theta,phi):
                double xi ; int l ; 
            
                mp->val_lx(rr, th, ph, l, xi) ;   // radial mapping assumption
            
                // Field value at this point:
            
                double vx = cvax.val_point(l, xi, th, ph) ;
                double vy = cvay.val_point(l, xi, th, ph) ;
                double vz = cvaz.val_point(l, xi, th, ph) ;

                fdata.width(14) ; fdata << vx ; 
                fdata.width(14) ; fdata << vy ; 
                fdata.width(14) ; fdata << vz ; 
                npoint++ ; 
             
                if (npoint == 3) {
                    fdata << "\n" ; 
                    npoint = 0 ; 
                }
            
            }
        }

    }
    
    if (npoint != 0) fdata << "\n" ; 

    fdata.close() ; 
    
    // --------------------------------------------------------
    // Header file for OpenDX
    // --------------------------------------------------------

    char* headername ;
    if (filename0 == 0x0) {
        headername = new char[30] ; 
        strcpy(headername, "vector3d_streamline.dxhead") ; 
    }
    else {
        headername = new char[ strlen(filename0)+9 ] ; 
        strcpy(headername, filename0) ; 
        strcat(headername, ".dxhead") ; 
    }

    ofstream fheader(headername) ;
    
    fheader << "file = " << filename << endl ; 
    fheader << "grid = " << nx << " x " << ny << " x " << nz << endl ; 
    fheader << "format = ascii" << endl ;  
    fheader << "interleaving = record-vector"  << endl ;
    fheader << "majority = column" << endl ; 
    fheader << "header = lines 5" << endl ; 
    fheader << "field = " << title_quotes << endl ; 
    fheader << "structure = 3-vector" << endl ; 
    fheader << "type = float" << endl ; 
    fheader << "dependency = positions" << endl ; 
    fheader << "positions = regular, regular, regular, " 
            << xmin << ", " << dx << ", " 
            << ymin << ", " << dy << ", " 
            << zmin << ", " << dz << endl ; 
    fheader << endl ; 
    fheader << "end" << endl ; 
    
    fheader.close() ; 
    

    if ( start_dx ) {       // Launch of OpenDX
        
        char* commande = new char[ strlen(headername) + 60 ] ;
        strcpy(commande, "ln -s ") ; 
        strcat(commande, headername) ; 
        strcat(commande, " visu_vector3d_SL.general") ; 
    
        system("rm -f visu_vector3d_SL.general") ; 
        system(commande) ;                      // ln -s headername visu_section.general
        system("dx -image visu_vector3d_SL.net &") ; 
    
        delete [] commande ;    
    }


    // Final cleaning
    // --------------
    
    if (vect_tmp != 0x0) delete vect_tmp ;
    delete [] title ; 
    delete [] title_quotes ; 
    delete [] filename ; 
    delete [] headername ;    
   

}   

}
