#include "type_parite.h"
#include "nbr_spx.h"
#include "grilles.h"
#include "map.h"
#include "graphique.h"
#include "monopole.h"
#include "utilitaires.h"

int main() {

  // THE GRID
  int nz_W = 6 ;
  int nz_H = 5 ;
  int nt = 5 ;
  int np = 4 ;
  int nr = 33 ;
  int typet = SYM ;
  int typep = SYM ;
  
  bool compact = true ;
  Mg3d grid_W (nz_W, nr, nt, np, typet, typep, compact) ;
  Mg3d grid_H (nz_H, nr, nt, np, typet, typep, compact) ;
  
  double beta ;
  cout << "Beta ? " << endl ;
  cin >> beta ;
  
  // THE MAPPING AFFINE
  double r_lim = 0.1 ;
  double* bounds = new double[nz_W+1] ;
  
  bounds[0] = 0 ;
  bounds[1] = r_lim ;
  for (int i=2 ; i<nz_W ; i++)
    bounds[i] = 2*bounds[i-1] ;
  bounds[nz_W] = __infinity ;
  
  Map_af mp_W (grid_W, bounds) ;
  delete [] bounds ;
  
  r_lim = (beta>5) ? 1./beta : 1. ;
  bounds = new double[nz_H+1] ;
  
  bounds[0] = 0 ;
  bounds[1] = r_lim ;
  for (int i=2 ; i<nz_H ; i++)
    bounds[i] = 2*bounds[i-1] ;
  bounds[nz_H] = __infinity ;
  
  Map_af mp_H (grid_H, bounds) ;
  delete [] bounds ;
   

  // The Monopole :
  double precis = 1e-10 ;
  int itemax = 1000 ;
  double relax = 0.2 ;

  Monopole champs(mp_W, mp_H, beta) ; 
  champs.solve_config(precis, itemax, relax) ;

  cout.precision (14) ;
  cout << "Parameter a = " << champs.give_a() << endl ;
  cout << "Parameter b = " << champs.give_b() << endl ;
  cout << "Energy      = " << champs.energy() << endl ;
  
  // Some plots
  //des_profile (champs.get_big_W(), 0, 1.1*mp_W.val_r (nz_W-1,-1, 0,0), 0, 0, "W") ;
  //des_profile (champs.get_big_H(), 0, 1.1*mp_H.val_r (nz_H-1,-1, 0,0) , 0, 0, "H") ;

  return 0 ;
}
