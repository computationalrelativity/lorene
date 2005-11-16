#include "type_parite.h"
#include "nbr_spx.h"
#include "grilles.h"
#include "map.h"
#include "graphique.h"
#include "monopole.h"
#include "utilitaires.h"

int main() {

  // THE GRID
  int nz = 4 ;
  int nt = 5 ;
  int np = 4 ;
  int nr = 45 ;
  int typet = SYM ;
  int typep = SYM ;
  
  bool compact = true ;
  Mg3d grid (nz, nr, nt, np, typet, typep, compact) ;

  // THE MAPPING AFFINE
  double r_lim = 1. ;
  double* bounds = new double[nz+1] ;
  
  bounds[0] = 0 ;
  bounds[1] = r_lim ;
  for (int i=2 ; i<nz ; i++)
    bounds[i] = 2*bounds[i-1] ;
  bounds[nz] = __infinity ;
  
  Map_af mp (grid, bounds) ;
  delete [] bounds ;
   
  double beta ;
  cout << "Beta ? " << endl ;
  cin >> beta ;
  while (cin.get() != '\n') ;

  // The Monopole :
  Monopole champs(mp, beta) ; 
  
  return EXIT_SUCCESS ;
}
