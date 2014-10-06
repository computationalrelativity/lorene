/*
 *  3D visualization of a Scalar via OpenDX 
 *
 *    (see file scalar.h for documentation).
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

char scalar_visu_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.8  2014/10/06 15:16:16  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.7  2005/02/18 13:14:19  j_novak
 * Changing of malloc/free to new/delete + suppression of some unused variables
 * (trying to avoid compilation warnings).
 *
 * Revision 1.6  2004/03/11 12:07:55  e_gourgoulhon
 * Added method visu_section_anim.
 *
 * Revision 1.5  2003/12/19 15:18:17  j_novak
 * Shadow variables hunt
 *
 * Revision 1.4  2003/12/16 06:32:57  e_gourgoulhon
 * Added method visu_box.
 *
 * Revision 1.3  2003/12/15 08:30:40  p_grandclement
 * Addition of #include <string.h>
 *
 * Revision 1.2  2003/12/14 21:49:14  e_gourgoulhon
 * Added argument start_dx (to launch OpenDX as a subprocess).
 *
 * Revision 1.1  2003/12/11 16:20:25  e_gourgoulhon
 * First version.
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

                    //-----------------------------------------//
                    //          visu_section : special cases   //
                    //-----------------------------------------//
                    
void Scalar::visu_section(const char section_type, double aa, double umin, 
        double umax, double vmin, double vmax, const char* title0, 
        const char* filename0, bool start_dx, int nu, int nv) const {

    Tbl plane(3,3) ;        
    plane.set_etat_qcq() ;  // Memory allocation for the Tbl

    switch (section_type) {
    
        case 'x' : {     
            plane.set(0,0) = aa ;   // Origin in the plane
            plane.set(0,1) = 0. ;   //  (absolute Cartesian coordinates)
            plane.set(0,2) = 0. ;   //
    
            plane.set(1,0) = 0. ;   // u-coordinate unit vector
            plane.set(1,1) = 1. ;   //  (absolute Cartesian components)
            plane.set(1,2) = 0. ;
    
            plane.set(2,0) = 0. ;   // v-coordinate unit vector
            plane.set(2,1) = 0. ;   //  (absolute Cartesian components)
            plane.set(2,2) = 1. ;

            visu_section(plane, umin, umax, vmin, vmax, title0, filename0, 
                         start_dx, nu, nv) ;
            break ; 
        }

        case 'y' : {     
            plane.set(0,0) = 0. ;   // Origin in the plane
            plane.set(0,1) = aa ;   //  (absolute Cartesian coordinates)
            plane.set(0,2) = 0. ;   //
    
            plane.set(1,0) = 1. ;   // u-coordinate unit vector
            plane.set(1,1) = 0. ;   //  (absolute Cartesian components)
            plane.set(1,2) = 0. ;
    
            plane.set(2,0) = 0. ;   // v-coordinate unit vector
            plane.set(2,1) = 0. ;   //  (absolute Cartesian components)
            plane.set(2,2) = 1. ;

            visu_section(plane, umin, umax, vmin, vmax, title0, filename0, 
                         start_dx, nu, nv) ;
            break ; 
        }

        case 'z' : {     
            plane.set(0,0) = 0. ;   // Origin in the plane
            plane.set(0,1) = 0. ;   //  (absolute Cartesian coordinates)
            plane.set(0,2) = aa ;   //
    
            plane.set(1,0) = 1. ;   // u-coordinate unit vector
            plane.set(1,1) = 0. ;   //  (absolute Cartesian components)
            plane.set(1,2) = 0. ;
    
            plane.set(2,0) = 0. ;   // v-coordinate unit vector
            plane.set(2,1) = 1. ;   //  (absolute Cartesian components)
            plane.set(2,2) = 0. ;

            visu_section(plane, umin, umax, vmin, vmax, title0, filename0, 
                         start_dx, nu, nv) ;
            break ; 
        }

        default : {
            cerr << "Scalar::visu_section : unknown type of section ! \n" ;
            cerr << "   section_type = " << section_type << endl ; 
            break ; 
        }
    }
        
}


                    //-----------------------------------------//
                    //          visu_section : general case    //
                    //-----------------------------------------//
                    

void Scalar::visu_section(const Tbl& plane, double umin, double umax, 
        double vmin, double vmax, const char* title0, const char* filename0,
        bool start_dx, int nu, int nv) const {
        
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
        strcpy(filename, "scalar_section.dxdata") ; 
    }
    else {
        filename = new char[ strlen(filename0)+8 ] ; 
        strcpy(filename, filename0) ; 
        strcat(filename, ".dxdata") ; 
    }

    ofstream fdata(filename) ; // output file
    
    fdata << title << "\n" ; 
    fdata << "size : " << nu << " x " << nv << "\n" ;
    fdata << "u_min = " << umin << "  u_max = " << umax << "\n" ; 
    fdata << "v_min = " << vmin << "  v_max = " << vmax << "\n" ; 

    // Plane characterization
    // ----------------------
    
    double xa0 = plane(0,0) ; 
    double ya0 = plane(0,1) ; 
    double za0 = plane(0,2) ;

    double eux = plane(1,0) ;  
    double euy = plane(1,1) ;  
    double euz = plane(1,2) ;  

    double evx = plane(2,0) ;  
    double evy = plane(2,1) ;  
    double evz = plane(2,2) ;  
    
    
    // The spectral coefficients are required
    va.coef() ; 
    const Mtbl_cf& cva = *(va.c_cf) ; 
    
    // What follows assumes that the mapping is radial:
    assert( dynamic_cast<const Map_radial*>(mp) != 0x0 ) ; 
    
    fdata.precision(5) ; 
    fdata.setf(ios::scientific) ; 
    
    // Loop on the points in the section plane
    // ---------------------------------------
    double du = (umax - umin) / double(nu-1) ; 
    double dv = (vmax - vmin) / double(nv-1) ; 
    
    int npoint = 0 ;    // number of data points per line in the file
    
    for (int j=0; j<nv; j++) {
        
        double v = vmin + dv * j ; 
        
        for (int i=0; i<nu; i++) {
        
            double u = umin + du * i ; 
            
            double xa = xa0 + u * eux + v * evx ;    
            double ya = ya0 + u * euy + v * evy ;    
            double za = za0 + u * euz + v * evz ;   
            

            // Values of (r,theta,phi) corresponding to (xa,ya,za) :
            double rr, th, ph ;  // polar coordinates of the mapping associated
                                 // to *this
            
            mp->convert_absolute(xa, ya, za, rr, th, ph) ; 

            // Values of (l,xi,theta',phi') corresponding to (r,theta,phi):
            double xi ; 
            int l ; 
            
            mp->val_lx(rr, th, ph, l, xi) ;   // radial mapping assumption
            
            // Field value at this point:
            
            double ff = cva.val_point(l, xi, th, ph) ;
            
            fdata.width(14) ; 
            fdata << ff ; 
            npoint++ ; 
             
            if (npoint == 6) {
                fdata << "\n" ; 
                npoint = 0 ; 
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
        strcpy(headername, "scalar_section.dxhead") ; 
    }
    else {
        headername = new char[ strlen(filename0)+9 ] ; 
        strcpy(headername, filename0) ; 
        strcat(headername, ".dxhead") ; 
    }

    ofstream fheader(headername) ;
    
    fheader << "file = " << filename << endl ; 
    fheader << "grid = " << nu << " x " << nv << endl ; 
    fheader << "format = ascii" << endl ;  
    fheader << "interleaving = record"  << endl ;
    fheader << "majority = column" << endl ; 
    fheader << "header = lines 4" << endl ; 
    fheader << "field = " << title_quotes << endl ; 
    fheader << "structure = scalar" << endl ; 
    fheader << "type = float" << endl ; 
    fheader << "dependency = positions" << endl ; 
    fheader << "positions = regular, regular, " << umin << ", " << du 
        << ", " << vmin << ", " << dv << endl ; 
    fheader << endl ; 
    fheader << "end" << endl ; 
    
    fheader.close() ; 
    

    if ( start_dx ) {       // Launch of OpenDX
        
        char* commande = new char[ strlen(headername) + 60 ] ;
        strcpy(commande, "ln -s ") ; 
        strcat(commande, headername) ; 
        strcat(commande, " visu_section.dxhead") ; 
    
        system("rm -f visu_section.dxhead") ; 
        system(commande) ;                      // ln -s headername visu_section.general
        system("dx -image visu_section.net &") ; 
    
        delete [] commande ;    
    }

    // Final cleaning
    // --------------
    delete [] title ; 
    delete [] title_quotes ; 
    delete [] filename ; 
    delete [] headername ;    
    
}   


                    //------------------------------//
                    //          visu_box            //
                    //------------------------------//

void Scalar::visu_box(double xmin, double xmax, double ymin, double ymax,
    double zmin, double zmax, const char* title0, const char* filename0, 
    bool start_dx, int nx, int ny, int nz) const {

    const Scalar* scal ; 
    Scalar* scal_tmp = 0x0 ; 
        
    // Decrease of dzpuis if dzpuis != 0 
    if ( !check_dzpuis(0) ) {
        scal_tmp = new Scalar(*this) ;                
        scal_tmp->dec_dzpuis(dzpuis) ; 
        scal = scal_tmp ; 
    }
    else{
        scal = this ; 
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
        strcpy(filename, "scalar_box.dxdata") ; 
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
    const Valeur& val = scal->va ; 
    val.coef() ; 
    const Mtbl_cf& cva = *(val.c_cf) ; 
    
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
            
                double ff = cva.val_point(l, xi, th, ph) ;

                fdata.width(14) ; 
                fdata << ff ; 
                npoint++ ; 
             
                if (npoint == 9) {
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
        strcpy(headername, "scalar_box.dxhead") ; 
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
    fheader << "interleaving = record"  << endl ;
    fheader << "majority = column" << endl ; 
    fheader << "header = lines 5" << endl ; 
    fheader << "field = " << title_quotes << endl ; 
    fheader << "structure = scalar" << endl ; 
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
        strcat(commande, " visu_scalar_box.dxhead") ; 
    
        system("rm -f visu_scalar_box.dxhead") ; 
        system(commande) ;                      // ln -s headername visu_section.general
        system("dx -image visu_scalar_box.net &") ; 
    
        delete [] commande ;    
    }


    // Final cleaning
    // --------------
    
    if (scal_tmp != 0x0) delete scal_tmp ;
    delete [] title ; 
    delete [] title_quotes ; 
    delete [] filename ; 
    delete [] headername ;    
   

}


                    //-------------------------------------//
                    //          visu_section_anim          //
                    //-------------------------------------//

                    
void Scalar::visu_section_anim(const char section_type, double aa, double umin, 
        double umax, double vmin, double vmax, int jtime, double , 
        int jgraph, const char* title, const char* filename_root, bool start_dx, 
        int nu, int nv) const {
        
    if ( jtime % jgraph != 0 ) return ;     
        
    // Preparation of the name of output file
    // --------------------------------------
    int k = jtime / jgraph ;
        
    char* filename ;
    if (filename_root == 0x0) {
        filename = new char[40] ; 
        strcpy(filename, "anim") ; 
    }
    else {
        filename = new char[ strlen(filename_root)+10 ] ; 
        strcpy(filename, filename_root) ; 
    }

    char nomk[5] ; 
    sprintf(nomk, "%04d", k) ; 
    strcat(filename, nomk) ; 
        
    // Call to visu_section to create the output file
    // ----------------------------------------------

    visu_section(section_type, aa, umin, umax, vmin, vmax, title, filename, 
                 false, nu, nv) ;   
         
    // Shall we start OpenDX ?
    // ---------------------

    if ( start_dx ) {       // Launch of OpenDX
        
        system("dx -edit anime.net &") ; 
    
    }

    // Final cleaning
    // --------------
        
    delete [] filename ;        
        
}           
