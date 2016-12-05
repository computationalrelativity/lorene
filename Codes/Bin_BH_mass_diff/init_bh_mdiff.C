/*
 * Main code for computing Misner-Lindquist solution 
 *
 */

/*
 *   Copyright (c) 2005 Philippe Grandclement
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
 * Revision 1.4  2016/12/05 16:18:22  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.3  2014/10/13 08:53:52  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2014/10/06 15:09:41  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.1  2005/08/29 15:10:19  p_grandclement
 * Addition of things needed :
 *   1) For BBH with different masses
 *   2) Provisory files for the mixted binaries (Bh and NS) : THIS IS NOT
 *   WORKING YET !!!
 *
 * 
 * $Header$
 *
 */

#include "headcpp.h"

//standard
#include <cstdlib>
#include <cmath>

// LORENE
#include "type_parite.h"
#include "nbr_spx.h"
#include "proto.h"
#include "coord.h"
#include "tenseur.h"
#include "bhole.h"
#include "utilitaires.h"
#include "graphique.h"


using namespace Lorene ;

int main(int argc, char** argv) {
    
    // Lecture du fichier de parametres :
     if (argc <2) {
	cout <<" Passer nom du ficher en arguments SVP !" << endl ;
	abort() ;
    }
    char blabla [119] ;
    char* name_fich = argv[1] ;
    ifstream param(name_fich) ;
    
    double  precis, relax, distance ;
    int nz_un, nz_deux, nbrer, nbret, nbrep ;
    
    param >> distance ; param.getline(blabla, 120) ;
    param >> nz_un ; param >> nz_deux ; param.getline(blabla, 120) ;
    param >> nbrer ; param >> nbret ; param >> nbrep ; param.getline(blabla, 120) ;
    
    double* bornes_un = new double[nz_un+1] ;
    for (int i=0 ; i<nz_un ; i++)
	param >> bornes_un[i] ;
    bornes_un[nz_un] = __infinity ;
    param.getline(blabla, 120) ;
    
    double* bornes_deux = new double[nz_deux+1] ;
    for (int i=0 ; i<nz_deux ; i++)
	param >> bornes_deux[i] ;
    bornes_deux[nz_deux] = __infinity ;
    param.getline(blabla, 120) ;
    
    param >> precis ; param.getline(blabla, 120) ;
    param >> relax ; param.getline(blabla, 120) ;
    
    param.close() ;
    
    int symetrie = NONSYM ; 
    // echantillonnage en phi :
    int* np_un = new int [nz_un] ;
    for (int l=0 ; l<nz_un ; l++)
	np_un[l] = nbrep ;
    int type_p = symetrie ;
    
    // echantillonnage en theta :
    int* nt_un = new int [nz_un] ;
    for (int l=0 ; l<nz_un ; l++)
	nt_un[l] = nbret ;
    int type_t = SYM ;
     
    // echantillonage en r :
    int* nr_un = new int [nz_un] ;
    for (int l=0 ; l<nz_un ; l++)
	nr_un[l] = nbrer ;

    int* type_r_un = new int[nz_un] ;
    type_r_un[0] = RARE ;
    for (int l=1 ; l<nz_un-1 ; l++)
	type_r_un[l] = FIN ;
    type_r_un[nz_un-1] = UNSURR ;
    
    Mg3d grille_un (nz_un, nr_un, type_r_un, nt_un, type_t, np_un, type_p) ;
    
     // echantillonnage en phi :
    int* np_deux = new int [nz_deux] ;
    for (int l=0 ; l<nz_deux ; l++)
	np_deux[l] = nbrep ;
    
    // echantillonnage en theta :
    int* nt_deux = new int [nz_deux] ;
    for (int l=0 ; l<nz_deux ; l++)
	nt_deux[l] = nbret ;
     
    // echantillonage en r :
    int* nr_deux = new int [nz_deux] ;
    for (int l=0 ; l<nz_deux ; l++)
	nr_deux[l] = nbrer ;

    int* type_r_deux = new int[nz_deux] ;
    type_r_deux[0] = RARE ;
    for (int l=1 ; l<nz_deux-1 ; l++)
	type_r_deux[l] = FIN ;
    type_r_deux[nz_deux-1] = UNSURR ;
    
    Mg3d grille_deux (nz_deux, nr_deux, type_r_deux, nt_deux, type_t, np_deux, type_p) ;   
    
    Map_af mapping_un (grille_un, bornes_un) ;
    Map_af mapping_deux (grille_deux, bornes_deux) ;
    
    mapping_un.set_ori (distance/2.,0,   0) ;
    mapping_deux.set_ori (-distance/2., 0, 0) ;
    mapping_deux.set_rot_phi (M_PI) ;
    
    Bhole_binaire bin (mapping_un, mapping_deux) ;
    bin.set_statiques(precis, relax) ;
    
    FILE* fich = fopen("statiques.dat", "w") ;
    grille_un.sauve(fich) ;
    grille_deux.sauve(fich) ;
    mapping_un.sauve(fich) ;
    mapping_deux.sauve(fich) ;
    bin(1).sauve(fich) ;
    bin(2).sauve(fich) ;
    fclose(fich) ;
  
    delete [] nr_un ;
    delete [] nt_un ;
    delete [] np_un ;
    delete [] type_r_un ;
    delete [] nr_deux ;
    delete [] nt_deux ;
    delete [] np_deux ;
    delete [] type_r_deux ; 
    delete [] bornes_un ;
    delete [] bornes_deux ;
    
    return 0 ;
}
