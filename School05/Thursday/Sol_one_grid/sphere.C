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
  int nt = 1 ;
  int np = 1 ;
  int nr = 33 ;
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

  // The Monopole :
  double precis = 1e-10 ;
  int itemax = 200 ;
  double relax = 0.5 ;

  Monopole champs(mp, beta) ; 
  champs.solve_config(precis, itemax, relax) ;

  cout << "Parameter a = " << champs.give_a() << endl ;
  cout << "Parameter b = " << champs.give_b() << endl ;
  
  // Some plots
  
  des_profile (champs.get_big_W(), 0, 5, 0, 0, "W") ;
  des_profile (champs.get_small_w(), 0, 5, 0, 0, "w") ;
  des_profile (champs.get_big_H(), 0, 5, 0, 0, "H") ;
  des_profile (champs.get_small_h(), 0, 5, 0, 0, "h") ;

  return 0 ;
}
