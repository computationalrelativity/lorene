/*
 * Main code for computing Misner-Lindquist solution 
 *
 */

/*
 *   Copyright (c) 2002 Philippe Grandclement
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
 * Revision 1.6  2016/12/05 16:18:22  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.5  2014/10/13 08:53:52  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.4  2014/10/06 15:09:41  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.3  2003/11/14 08:14:52  p_grandclement
 * Correction of the codes for binary black holes in circular orbit
 *
 * Revision 1.2  2003/04/09 15:05:09  e_gourgoulhon
 * Added copyright and log messages.
 *
 *
 * Revision 1.1
 * date: 2003/03/31 16:08:23;  author: p_grandclement;  state: Exp;
 * Add codes for the binary black holes computation :
 *       - init_bh.C : computes the Misner-Lindquist solution (used as initial guess for the coal code)
 *       - param_statiques.d : example of parameter file for init_bh
 *       - coal_bh.C : computation of the binary
 *       - param_coal.d : example of the paramater file for coal_bh
 *       - lit_holes_bin.C : reads a configuration.
 * 
 * $Header$
 *
 */

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
    char blabla [120] ;
    char* name_fich = argv[1] ;
    ifstream param(name_fich) ;
    
    double  precis, relax, rayon, beta ;
    int nz, nbrer, nbret, nbrep ;
    
    param >> beta ; param.getline(blabla, 120) ;
    param >> nz ; param.getline(blabla, 120) ;
    param >> nbrer ; param >> nbret ; param >> nbrep ; param.getline(blabla, 120) ;
    
    double* bornes = new double[nz+1] ;
    for (int i=0 ; i<nz ; i++)
	param >> bornes[i] ;
    bornes[nz] = __infinity ;
    rayon = bornes[1] ;
    param.getline(blabla, 120) ;
    
    param >> precis ; param.getline(blabla, 120) ;
    param >> relax ; param.getline(blabla, 120) ;
   
    double distance = rayon*beta ;
    
    param.close() ;
    
    int symetrie = NONSYM ; 
    // echantillonnage en phi :
    int* np = new int [nz] ;
    for (int l=0 ; l<nz ; l++)
	np[l] = nbrep ;
    int type_p = symetrie ;
    
    // echantillonnage en theta :
    int* nt = new int [nz] ;
    for (int l=0 ; l<nz ; l++)
	nt[l] = nbret ;
    int type_t = SYM ;
     
    // echantillonage en r :
    int* nr = new int [nz] ;
    for (int l=0 ; l<nz ; l++)
	nr[l] = nbrer ;

    int* type_r = new int[nz] ;
    type_r[0] = RARE ;
    for (int l=1 ; l<nz-1 ; l++)
	type_r[l] = FIN ;
    type_r[nz-1] = UNSURR ;
    
    Mg3d grille (nz, nr, type_r, nt, type_t, np, type_p) ;
    
    Map_af mapping_un (grille, bornes) ;
    Map_af mapping_deux (grille, bornes) ;
    
    mapping_un.set_ori (distance/2.,0,   0) ;
    mapping_deux.set_ori (-distance/2., 0, 0) ;
    mapping_deux.set_rot_phi (M_PI) ;

    Bhole_binaire bin (mapping_un, mapping_deux) ;
    bin.set_statiques(precis, relax) ;
    
    FILE* fich = fopen("statiques.dat", "w") ;
    grille.sauve(fich) ;
    mapping_un.sauve(fich) ;
    mapping_deux.sauve(fich) ;
    bin(1).sauve(fich) ;
    bin(2).sauve(fich) ;
    fclose(fich) ;
  
    delete [] nr ;
    delete [] nt ;
    delete [] np ;
    delete [] type_r ;
    delete [] bornes ;

    return 1 ;
}
