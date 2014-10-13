#ifndef __MONOPOLE_H_ 
#define __MONOPOLE_H_ 

#include "map.h"
#include "scalar.h"

using namespace Lorene ;

class Monopole {

 protected:
  const Map_af& mp ; /// Mapping 

  // radius
  Scalar radius ;
  int nz ; ///Number of domains
  
  // DATA : 
  double beta ;
  
  Scalar big_W ; /// Function W (spherical symetry ?)
  Scalar small_w ; /// Representation of W in terms of w
  
  Scalar big_H ; /// Function H
  Scalar small_h ; /// Auxiliary variable close to r=0 ;
  
  // Constructors destructor
 public:
  explicit Monopole (const Map_af&, double) ;
  /// Standard constructor 
  Monopole (const Monopole&) ; /// Constructor by copy
  Monopole(const Map_af&, FILE*) ; /// Constructor from a file
  
  virtual ~Monopole() ; /// Destructor

  // Assigment
 public:
  void operator= (const Monopole&) ; /// Assignment to another monopole

  // Accessor :
 public:
  const Map_af& get_mp() const {return mp ;} ; /// Return the mapping
  
  Scalar get_big_W() const {return big_W ;} ; /// Return the field
  /*
   * Return the w-representation of the field
   */
  Scalar get_small_w() const {return small_w ;} ;
  
  Scalar get_big_H() const {return big_H ;} ; /// Return the field
  /*
   * Return the w-representation of the field
   */
  Scalar get_small_h() const {return small_h ;} ;
    
  // Output :
  public :
    void sauve(FILE *) const ; /// Save to a file
  
  // The numerical methods (at last) :
 private:
  void init_big_W() ; /// Initialize big W to an anzats
  
  void init_big_H() ; /// Initialize big H to an anzats
  /*
   * Computes the auxiliary variable w. Outside the domain  zone, we keep 
   * big_W
   */
  void do_small_w() ;

  /*
   * Computes big_W from small_w
   */
  void do_big_W() ;

  /*
   * Computes the auxiliary variable h. Outside the domain  zone, we keep 
   * big_H
   */
  void do_small_h() ;

  /*
   * Computes big_H from small_h
   */
  void do_big_H() ;

  /*
   * Compute the source of the equation for W , using either the small_w version 
   * or the W one, depending on the region
   */

  Scalar compute_source_W() const ;
 
  /*
   * Compute the source of the equation for H , using either the small_h version 
   * or the H one, depending on the region
   */

  Scalar compute_source_H() const ;
 
 public:
  /**
   * Solve system
   *
   **/
  int solve_config (double, int, double) ;

  double give_a () const ;
  double give_b() const ;
 
} ;

# endif
