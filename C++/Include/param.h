/*
 *  Definition of Lorene class Param
 *
 */

/*
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
 *   Copyright (c) 2000-2001 Jerome Novak
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


#ifndef __PARAM_H_ 
#define __PARAM_H_ 


/*
 * $Id$
 * $Log$
 * Revision 1.3  2003/09/25 12:08:02  j_novak
 * Tensors can be stored in Param objects
 *
 * Revision 1.2  2002/09/19 09:52:42  j_novak
 * Added objects Qtenseur and Qmetrique for 4D tensor and metric handling.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 1.10  2001/10/27  09:26:24  novak
 * *** empty log message ***
 *
 * Revision 1.9  2001/10/11 07:44:12  eric
 * Ajout du stokage des Etoile's
 *
 * Revision 1.8  2000/10/24  14:54:49  novak
 * Added the function clean_all()
 *
 * Revision 1.7  2000/05/25 12:39:19  eric
 * MODIFICATION MAJEURE: pour les int et les double, ce sont desormais les
 * adresses qui sont stokees, et non plus les nombres eux-memes
 * (le traitement des int et des double est donc desormais completement
 * aligne sur celui des Tbl, Cmp, etc...)
 *
 * Revision 1.6  1999/12/29  13:10:39  eric
 * Ajout du stokage des Mtbl_cf.
 *
 * Revision 1.5  1999/12/27  12:16:43  eric
 * Ajout du stokage des mappings (class Map).
 *
 * Revision 1.4  1999/12/16  10:27:40  eric
 * Ajout des membres modifiables.
 * Par defaut, les objets listes sont const.
 *
 * Revision 1.3  1999/12/15  16:49:52  eric
 * *** empty log message ***
 *
 * Revision 1.2  1999/12/15  16:22:36  eric
 * Changement de l'ordre des arguments dans add_*
 * Argument par defaut: position = 0
 * Ajout du stokage des int et des double.
 *
 * Revision 1.1  1999/12/13  14:35:56  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */

class Tbl ; 
class Itbl ; 
class Mtbl_cf ; 
class Map ; 
class Cmp ; 
class Tenseur ;
class Qtenseur ;
class Tensor ;
class Etoile ;

/** Parameter storage.
 * 
 *  This class is intended to store addresses of various Lorene objects to
 *  pass them as parameters in some subroutines. 
 * 
 */
class Param {

    // Data : 
    // -----
    private:
	int n_int ;	/// Number of {\tt int}'s (integers). 
	/// Array (size {\tt n\_int}) of the {\tt int}'s addresses.
	const int** p_int ;	

	int n_int_mod ;	/// Number of modifiable {\tt int}'s (integers). 
	/// Array (size {\tt n\_int\_mod}) of the modifiable {\tt int}'s addresses.
	int** p_int_mod ;	

	int n_double ; /// Number of {\tt double}'s (double precis. numbers). 
	/// Array (size {\tt n\_double}) of the {\tt double}'s addresses.
	const double** p_double ; 

	int n_double_mod ; /// Number of modifiable {\tt double}'s (double precis. numbers). 
	/// Array (size {\tt n\_double\_mod}) of the {\tt double}'s addresses
	double** p_double_mod ; 


	int n_tbl ;	/// Number of {\tt Tbl}'s 
	/// Array (size {\tt n\_tbl}) of the {\tt Tbl}'s addresses
	const Tbl** p_tbl ;	
    
	int n_tbl_mod ;	/// Number of modifiable {\tt Tbl}'s 
	/// Array (size {\tt n\_tbl\_mod}) of the modifiable {\tt Tbl}'s addresses
	Tbl** p_tbl_mod ;	
    
	int n_itbl ;	/// Number of {\tt Itbl}'s 
	/// Array (size {\tt n\_itbl}) of the {\tt Itbl}'s addresses
	const Itbl** p_itbl ;	
    
	int n_itbl_mod ;    /// Number of modifiable {\tt Itbl}'s 
	/// Array (size {\tt n\_itbl\_mod}) of the modifiable {\tt Itbl}'s addresses
	Itbl** p_itbl_mod ;	
    
	int n_cmp ;	/// Number of {\tt Cmp}'s 
	/// Array (size {\tt n\_cmp}) of the {\tt Cmp}'s addresses
	const Cmp** p_cmp ;	
    
	int n_cmp_mod ;	/// Number of modifiable {\tt Cmp}'s 
	/// Array (size {\tt n\_cmp\_mod}) of the modifiable {\tt Cmp}'s addresses
	Cmp** p_cmp_mod ;	
    
	int n_tenseur ;	/// Number of {\tt Tenseur}'s 
	/// Array (size {\tt n\_tenseur}) of the {\tt Tenseur}'s addresses
	const Tenseur** p_tenseur ;	
    
	int n_tenseur_mod ;	/// Number of modifiable {\tt Tenseur}'s 
	/// Array (size {\tt n\_tenseur\_mod}) of the modifiable {\tt Tenseur}'s addresses
	Tenseur** p_tenseur_mod ;	

	int n_qtenseur ;	/// Number of {\tt Qtenseur}'s 
	/// Array (size {\tt n\_qtenseur}) of the {\tt Qtenseur}'s addresses
	const Qtenseur** p_qtenseur ;	
    
	int n_qtenseur_mod ;	/// Number of modifiable {\tt Qtenseur}'s 
	/// Array (size {\tt n\_qtenseur\_mod}) of the modifiable {\tt Qtenseur}'s addresses
	Qtenseur** p_qtenseur_mod ;	

	int n_map ;	/// Number of {\tt Map}'s 
	/// Array (size {\tt n\_map}) of the {\tt Map}'s addresses
	const Map** p_map ;	
    
	int n_mtbl_cf ;	/// Number of {\tt Mtbl\_cf}'s 
	/// Array (size {\tt n\_mtbl\_cf}) of the {\tt Mtbl\_cf}'s addresses
	const Mtbl_cf** p_mtbl_cf ;	
    
	int n_tensor ;	/// Number of {\tt Tensor}'s 
	/// Array (size {\tt n\_tensor}) of the {\tt Tensor}'s addresses
	const Tensor** p_tensor ;	
    
	int n_tensor_mod ;	/// Number of modifiable {\tt Tensor}'s 
	/// Array (size {\tt n\_tensor\_mod}) of the modifiable {\tt Tensor}'s addresses
	Tensor** p_tensor_mod ;
	
	int n_etoile ;	/// Number of {\tt Etoile}'s
	/// Array (size {\tt n\_etoile}) of the {\tt Etoile}'s addresses
	const Etoile** p_etoile ;	


    // Constructors - Destructor
    // -------------------------
	
    public:
	Param() ;	/// Default constructor is the only constructor

    private:
	/** Copy constructor (private and not implemented to make {\tt Param}
	 * a non-copyable class)
	 */ 
	Param(const Param& ) ;
	
    public: 
	~Param() ;	/// Destructor

	/** Deletes all the objects stored as modifiables,
	 *  i.e. all quantities with the suffix {\tt mod}.
	 */
	void clean_all() ; 	



    // Assignment
    // -----------
    private: 
	/** Assignment operator (private and not implemented to make 
	 *   {\tt Param} a non-copyable class)
	 */
	void operator=(const Param& ) ;
	 	

    // Addition/Extraction of one element
    // ----------------------------------
    public:
    
	///Returns the number of stored {\tt int}'s addresses.
	int get_n_int() const ; 
    
	/** Adds the address of a new {\tt int} to the list.
	 * 
	 *  @param n [input] {\tt int} the address of which is 
	 *                             to be stored
	 *  @param position [input] position of the {\tt int} in the list
	 *			    of stored {\tt int} addresses (default
	 *			    value = 0)
	 * 
	 */
	void add_int(const int& n, int position = 0) ;
	
	/** Returns the reference of a {\tt int} stored in the list.
	 * 
	 *  @param position [input] position of the {\tt int} in the list
	 *			    of stored {\tt int} addresses (default
	 *			    value = 0)
	 *  @return Reference to the {\tt int} the address of 
	 *           which is stored at the location  {\tt position} in the 
	 *           list
	 */
	const int& get_int(int position = 0) const; 

	///Returns the number of modifiable {\tt int}'s addresses in the list.
	int get_n_int_mod() const ; 
    
	/** Adds the address of a new modifiable {\tt int} to the list.
	 * 
	 *  @param n [input] modifiable {\tt int} the address of which is 
	 *                   to be stored
	 *  @param position [input] position of the modifiable {\tt int} 
	 *                in the list of stored modifiable {\tt int} addresses
	 *                (default value = 0)
	 * 
	 */
	void add_int_mod(int& n, int position = 0) ;
	
	/** Returns the reference of a modifiable {\tt int} stored in the list.
	 * 
	 *  @param position [input] position of the modifiable {\tt int} 
	 *              in the list of stored modifiable {\tt int} addresses 
	 *                    (default value = 0)
	 *  @return Reference to the modifiable {\tt int} the address of 
	 *           which is stored at the location  {\tt position} in the 
	 *           list
	 */
	int& get_int_mod(int position = 0) const; 
	

	///Returns the number of stored {\tt double}'s addresses.
	int get_n_double() const ; 
    
	/** Adds the the address of a new {\tt double} to the list.
	 * 
	 *  @param x [input] {\tt double} the address of which is 
	 *                             to be stored
	 *  @param position [input] position of the {\tt double} in the list
	 *			    of stored {\tt double} addresses (default
	 *			    value = 0)
	 * 
	 */
	void add_double(const double& x, int position = 0) ;
	
	/** Returns the reference of a {\tt double} stored in the list.
	 * 
	 *  @param position [input] position of the {\tt double} in the list
	 *			    of stored {\tt double} addresses (default
	 *			    value = 0)
	 *  @return Reference to the {\tt double} the address of 
	 *           which is stored at the location  {\tt position} in the 
	 *           list
	 */
	const double& get_double(int position = 0) const; 
	

	///Returns the number of stored modifiable {\tt double}'s addresses.
	int get_n_double_mod() const ; 
    
	/** Adds the address of a new modifiable {\tt double} to the list.
	 * 
	 *  @param x [input] modifiable {\tt double} the address of which is 
	 *                    to be stored
	 *  @param position [input] position of the modifiable {\tt double} 
	 *                          in the list of stored modifiable 
	 *                      {\tt double} addresses (default value = 0)
	 * 
	 */
	void add_double_mod(double& x, int position = 0) ;
	
	/** Returns the reference of a stored modifiable {\tt double}.
	 * 
	 *  @param position [input] position of the modifiable {\tt double} 
	 *                          in the list of stored modifiable 
	 *                      {\tt double} addresses (default value = 0)
	 * 
	 *  @return Reference to the modifiable {\tt double} the address of 
	 *                  which is stored at 
	 *		    the location  {\tt position} in the list
	 */
	double& get_double_mod(int position = 0) const; 
	

	///Returns the number of {\tt Tbl}'s addresses in the list.
	int get_n_tbl() const ; 
    
	/** Adds the address of a new {\tt Tbl} to the list.
	 * 
	 *  @param ti [input] {\tt Tbl} the address of which is to be stored
	 *  @param position [input] position of the {\tt Tbl} in the list
	 *			    of stored {\tt Tbl} addresses (default
	 *			    value = 0)
	 * 
	 */
	void add_tbl(const Tbl& ti, int position = 0) ;
	
	/** Returns the reference of a {\tt Tbl} stored in the list.
	 * 
	 *  @param position [input] position of the {\tt Tbl} in the list
	 *			    of stored {\tt Tbl} addresses (default
	 *			    value = 0)
	 *  @return Reference to the {\tt Tbl} the address of which is stored at 
	 *		    the location  {\tt position} in the list
	 */
	const Tbl& get_tbl(int position = 0) const; 
	

	///Returns the number of modifiable {\tt Tbl}'s addresses in the list.
	int get_n_tbl_mod() const ; 
    
	/** Adds the address of a new modifiable {\tt Tbl} to the list.
	 * 
	 *  @param ti [input] modifiable {\tt Tbl} the address of which is 
	 *                    to be stored
	 *  @param position [input] position of the {\tt Tbl} in the list
	 *			    of stored modifiable {\tt Tbl} addresses 
	 *                          (default value = 0)
	 */
	void add_tbl_mod(Tbl& ti, int position = 0) ;
	
	/** Returns the reference of a modifiable {\tt Tbl} stored in the list.
	 * 
	 *  @param position [input] position of the {\tt Tbl} in the list
	 *			    of stored modifiable {\tt Tbl} addresses 
	 *                          (default value = 0)
	 *  @return Reference to the modifiable {\tt Tbl} the address of 
	 *           which is stored at  the location  {\tt position} in the 
	 *           list
	 */
	 Tbl& get_tbl_mod(int position = 0) const; 
	

	///Returns the number of {\tt Itbl}'s addresses in the list.
	int get_n_itbl() const ; 
    
	/** Adds the address of a new {\tt Itbl} to the list.
	 * 
	 *  @param ti [input] {\tt Itbl} the address of which is to be stored
	 *  @param position [input] position of the {\tt Itbl} in the list
	 *			    of stored {\tt Itbl} addresses (default
	 *			    value = 0)
	 * 
	 */
	void add_itbl(const Itbl& ti, int position = 0) ;
	
	/** Returns the reference of a {\tt Itbl} stored in the list.
	 * 
	 *  @param position [input] position of the {\tt Itbl} in the list
	 *			    of stored {\tt Itbl} addresses (default
	 *			    value = 0)
	 *  @return Reference to the {\tt Itbl} the address of which is stored at 
	 *		    the location  {\tt position} in the list
	 */
	const Itbl& get_itbl(int position = 0) const; 
	

	///Returns the number of modifiable {\tt Itbl}'s addresses in the list.
	int get_n_itbl_mod() const ; 
    
	/** Adds the address of a new modifiable {\tt Itbl} to the list.
	 * 
	 *  @param ti [input] modifiable {\tt Itbl} the address of which is 
	 *         to be stored
	 *  @param position [input] position of the {\tt Itbl} in the list
	 *			    of stored modifiable {\tt Itbl} addresses 
	 *                          (default value = 0)
	 * 
	 */
	void add_itbl_mod(Itbl& ti, int position = 0) ;
	
	/** Returns the reference of a stored modifiable {\tt Itbl}.
	 * 
	 *  @param position [input] position of the {\tt Itbl} in the list
	 *			    of stored modifiable {\tt Itbl} addresses 
	 *                          (default value = 0)
	 *  @return Reference to the modifiable {\tt Itbl} the address of 
	 *            which is stored at the location  {\tt position} 
	 *            in the list
	 */
	Itbl& get_itbl_mod(int position = 0) const; 
	

	///Returns the number of {\tt Cmp}'s addresses in the list.
	int get_n_cmp() const ; 
    
	/** Adds the address of a new {\tt Cmp} to the list.
	 * 
	 *  @param ti [input] {\tt Cmp} the address of which is to be stored
	 *  @param position [input] position of the {\tt Cmp} in the list
	 *			    of stored {\tt Cmp} addresses (default
	 *			    value = 0)
	 * 
	 */
	void add_cmp(const Cmp& ti, int position = 0) ;
	
	/** Returns the reference of a {\tt Cmp} stored in the list.
	 * 
	 *  @param position [input] position of the {\tt Cmp} in the list
	 *			    of stored {\tt Cmp} addresses (default
	 *			    value = 0)
	 *  @return Reference to the {\tt Cmp} the address of which is stored at 
	 *		    the location  {\tt position} in the list
	 */
	const Cmp& get_cmp(int position = 0) const; 
	

	///Returns the number of modifiable {\tt Cmp}'s addresses in the list.
	int get_n_cmp_mod() const ; 
    
	/** Adds the address of a new modifiable {\tt Cmp} to the list.
	 * 
	 *  @param ti [input] modifiable {\tt Cmp} the address of which is 
	 *                    to be stored
	 *  @param position [input] position of the {\tt Cmp} in the list
	 *			    of stored modifiable {\tt Cmp} addresses 
	 *                          (default value = 0)
	 */
	void add_cmp_mod(Cmp& ti, int position = 0) ;
	
	/** Returns the reference of a modifiable {\tt Cmp} stored in the list.
	 * 
	 *  @param position [input] position of the {\tt Cmp} in the list
	 *			    of stored modifiable {\tt Cmp} addresses 
	 *                          (default value = 0)
	 *  @return Reference to the modifiable {\tt Cmp} the address of 
	 *           which is stored at  the location  {\tt position} in the 
	 *           list
	 */
	 Cmp& get_cmp_mod(int position = 0) const; 
	

	///Returns the number of {\tt Tenseur}'s addresses in the list.
	int get_n_tenseur() const ; 
    
	/** Adds the address of a new {\tt Tenseur} to the list.
	 * 
	 *  @param ti [input] {\tt Tenseur} the address of which is to be stored
	 *  @param position [input] position of the {\tt Tenseur} in the list
	 *			    of stored {\tt Tenseur} addresses (default
	 *			    value = 0)
	 * 
	 */
	void add_tenseur(const Tenseur& ti, int position = 0) ;
	
	/** Returns the reference of a {\tt Tenseur} stored in the list.
	 * 
	 *  @param position [input] position of the {\tt Tenseur} in the list
	 *			    of stored {\tt Tenseur} addresses (default
	 *			    value = 0)
	 *  @return Reference to the {\tt Tenseur} the address of which is stored at 
	 *		    the location  {\tt position} in the list
	 */
	const Tenseur& get_tenseur(int position = 0) const; 
	

	///Returns the number of modifiable {\tt Tenseur}'s addresses in the list.
	int get_n_tenseur_mod() const ; 
    
	/** Adds the address of a new modifiable {\tt Tenseur} to the list.
	 * 
	 *  @param ti [input] modifiable {\tt Tenseur} the address of which is 
	 *                    to be stored
	 *  @param position [input] position of the {\tt Tenseur} in the list
	 *			    of stored modifiable {\tt Tenseur} addresses 
	 *                          (default value = 0)
	 */
	void add_tenseur_mod(Tenseur& ti, int position = 0) ;
	
	/** Returns the reference of a modifiable {\tt Tenseur} stored in the list.
	 * 
	 *  @param position [input] position of the {\tt Tenseur} in the list
	 *			    of stored modifiable {\tt Tenseur} addresses 
	 *                          (default value = 0)
	 *  @return Reference to the modifiable {\tt Tenseur} the address of 
	 *           which is stored at  the location  {\tt position} in the 
	 *           list
	 */
	 Tenseur& get_tenseur_mod(int position = 0) const; 
	
	///Returns the number of {\tt Qtenseur}'s addresses in the list.
	int get_n_qtenseur() const ; 
    
	/** Adds the address of a new {\tt Qtenseur} to the list.
	 * 
	 *  @param ti [input] {\tt Qtenseur} the address of which is to be stored
	 *  @param position [input] position of the {\tt Qtenseur} in the list
	 *			    of stored {\tt Qtenseur} addresses (default
	 *			    value = 0)
	 * 
	 */
	void add_qtenseur(const Qtenseur& ti, int position = 0) ;
	
	/** Returns the reference of a {\tt Qtenseur} stored in the list.
	 * 
	 *  @param position [input] position of the {\tt Qtenseur} in the list
	 *			    of stored {\tt Qtenseur} addresses (default
	 *			    value = 0)
	 *  @return Reference to the {\tt Qtenseur} the address of which is stored at 
	 *		    the location  {\tt position} in the list
	 */
	const Qtenseur& get_qtenseur(int position = 0) const; 
	

	///Returns the number of modifiable {\tt Qtenseur}'s addresses in the list.
	int get_n_qtenseur_mod() const ; 
    
	/** Adds the address of a new modifiable {\tt Qtenseur} to the list.
	 * 
	 *  @param ti [input] modifiable {\tt Qtenseur} the address of which is 
	 *                    to be stored
	 *  @param position [input] position of the {\tt Qtenseur} in the list
	 *			    of stored modifiable {\tt Qtenseur} addresses 
	 *                          (default value = 0)
	 */
	void add_qtenseur_mod(Qtenseur& ti, int position = 0) ;
	
	/** Returns the reference of a modifiable {\tt Qtenseur} stored in the list.
	 * 
	 *  @param position [input] position of the {\tt Qtenseur} in the list
	 *			    of stored modifiable {\tt Qtenseur} addresses 
	 *                          (default value = 0)
	 *  @return Reference to the modifiable {\tt Qtenseur} the address of 
	 *           which is stored at  the location  {\tt position} in the 
	 *           list
	 */
	 Qtenseur& get_qtenseur_mod(int position = 0) const; 
	
	///Returns the number of {\tt Map}'s addresses in the list.
	int get_n_map() const ; 
    
	/** Adds the address of a new {\tt Map} to the list.
	 * 
	 *  @param mi [input] {\tt Map} the address of which is to be stored
	 *  @param position [input] position of the {\tt Map} in the list
	 *			    of stored {\tt Map} addresses (default
	 *			    value = 0)
	 * 
	 */
	void add_map(const Map& mi, int position = 0) ;
	
	/** Returns the reference of a {\tt Map} stored in the list.
	 * 
	 *  @param position [input] position of the {\tt Map} in the list
	 *			    of stored {\tt Map} addresses (default
	 *			    value = 0)
	 *  @return Reference to the {\tt Map} the address of which is stored at 
	 *		    the location  {\tt position} in the list
	 */
	const Map& get_map(int position = 0) const; 
	
	///Returns the number of {\tt Mtbl\_cf}'s addresses in the list.
	int get_n_mtbl_cf() const ; 
    
	/** Adds the address of a new {\tt Mtbl\_cf} to the list.
	 * 
	 *  @param mi [input] {\tt Mtbl\_cf} the address of which is to be stored
	 *  @param position [input] position of the {\tt Mtbl\_cf} in the list
	 *			    of stored {\tt Mtbl\_cf} addresses (default
	 *			    value = 0)
	 * 
	 */
	void add_mtbl_cf(const Mtbl_cf& mi, int position = 0) ;
	
	/** Returns the reference of a {\tt Mtbl\_cf} stored in the list.
	 * 
	 *  @param position [input] position of the {\tt Mtbl\_cf} in the list
	 *			    of stored {\tt Mtbl\_cf} addresses (default
	 *			    value = 0)
	 *  @return Reference to the {\tt Mtbl\_cf} the address of which is stored at 
	 *		    the location  {\tt position} in the list
	 */
	const Mtbl_cf& get_mtbl_cf(int position = 0) const; 
	
	///Returns the number of {\tt Tensor}'s addresses in the list.
	int get_n_tensor() const ; 
    
	/** Adds the address of a new {\tt Tensor} to the list.
	 * 
	 *  @param ti [input] {\tt Tensor} the address of which is to be stored
	 *  @param position [input] position of the {\tt Tensor} in the list
	 *			    of stored {\tt Tensor} addresses (default
	 *			    value = 0)
	 * 
	 */
	void add_tensor(const Tensor& ti, int position = 0) ;
	
	/** Returns the reference of a {\tt Tensor} stored in the list.
	 * 
	 *  @param position [input] position of the {\tt Tensor} in the list
	 *			    of stored {\tt Tensor} addresses (default
	 *			    value = 0)
	 *  @return Reference to the {\tt Tensor} the address of which is stored at 
	 *		    the location  {\tt position} in the list
	 */
	const Tensor& get_tensor(int position = 0) const; 
	

	///Returns the number of modifiable {\tt Tensor}'s addresses in the list.
	int get_n_tensor_mod() const ; 
    
	/** Adds the address of a new modifiable {\tt Tensor} to the list.
	 * 
	 *  @param ti [input] modifiable {\tt Tensor} the address of which is 
	 *                    to be stored
	 *  @param position [input] position of the {\tt Tensor} in the list
	 *			    of stored modifiable {\tt Tensor} addresses 
	 *                          (default value = 0)
	 */
	void add_tensor_mod(Tensor& ti, int position = 0) ;
	
	/** Returns the reference of a modifiable {\tt Tensor} stored in the list.
	 * 
	 *  @param position [input] position of the {\tt Tensor} in the list
	 *			    of stored modifiable {\tt Tensor} addresses 
	 *                          (default value = 0)
	 *  @return Reference to the modifiable {\tt Tensor} the address of 
	 *           which is stored at  the location  {\tt position} in the 
	 *           list
	 */
	 Tensor& get_tensor_mod(int position = 0) const; 
	///Returns the number of {\tt Etoile}'s addresses in the list.
	int get_n_etoile() const ;

	/** Adds the address of a new {\tt Etoile} to the list.
	 *
	 *  @param mi [input] {\tt Etoile} the address of which is to be stored
	 *  @param position [input] position of the {\tt Etoile} in the list
	 *			    of stored {\tt Etoile} addresses (default
	 *			    value = 0)
	 *
	 */
	void add_etoile(const Etoile& eti, int position = 0) ;
	
	/** Returns the reference of a {\tt Etoile} stored in the list.
	 *
	 *  @param position [input] position of the {\tt Etoile} in the list
	 *			    of stored {\tt Etoile} addresses (default
	 *			    value = 0)
	 *  @return Reference to the {\tt Etoile} the address of which is stored at
	 *		    the location  {\tt position} in the list
	 */
	const Etoile& get_etoile(int position = 0) const;
	

};

#endif
