			//------------------------------------//
			//         class Map_star             //
			//------------------------------------//

/*
 * Affine and starlike radial mapping to describe 3D star. \ingroup (map)
 * 
 * The affine radial mapping is the simplest one between the grid coordinates
 * \f$(\xi, \theta', \phi')\f$ and the physical coordinates \f$(r, \theta, \phi)\f$. It comprises a single domain (nucleus, shells to be added in the future)
 * It is defined by \f$\theta=\theta'\f$, \f$\phi=\phi'\f$ and 
 *  \li \f$r=\alpha \xi + \beta\f$, in non-compactified domains, 
 * where \f$\alpha\f$ and \f$\beta\f$ depend upon the angular direction. 
 * 
 *
 */

class Map_star : public Map_radial{


	// Data :
    // ----
    private:
	/// Array (size: \c mg->nzone*Nt*Np ) of the values of \f$\alpha\f$ in each domain
	Valeur alpha ;	 
	Valeur beta  ;


	// Constructors, destructor : 
    // ------------------------
    public:
	/**
	 * Standard Constructor
	 * @param mgrille  [input] Multi-domain grid on which the mapping is defined
	 * @param r_limits [input] Array (size: \c Nt*Np ) (only 1 domain at the moment) of the
	 *			   value of \e r  at the boundaries of the various 
	 *			   domains : 
	 *			   \li \c r_limits[l,k,j] : outer boundary of the 
	 *				 domain no. \c l in the angular direction \c j,k 
	 */
	Map_star(const Mg3d& mgrille, const double* r_limits) ;	
	/**
	 * Standard Constructor with Tbl
	 * @param mgrille  [input] Multi-domain grid on which the mapping is defined
	 * @param r_limits [input] Array (size: \c Nt*Np ) (only 1 domain at the moment) of the
	 *			   value of \e r  at the boundaries of the various 
	 *			   domains : 
	 *			   \li \c r_limits[l,k,j] : outer boundary of the 
	 *				 domain no. \c l in the angular direction \c j,k 
	 */
	Map_star(const Mg3d& mgrille, const Tbl& r_limits) ;	
	
	Map_star(const Map_star& ) ;      ///< Copy constructor
	Map_star(const Mg3d&, const string&) ;///< Constructor from a formatted file
	Map_star(const Mg3d&, FILE* ) ; ///< Constructor from a file (see \c sauve(FILE*) )


	virtual ~Map_star() ;	      ///< Destructor

	/** Returns the "angular" mapping for the outside of domain \c l_zone.
	 * Valid only for the class \c Map_af.
	 */
	virtual const Map_af& mp_angu(int) const ;

	// Assignment
    // ----------
    public: 
	/// Assignment to an affine mapping.
	virtual void operator=(const Map_af& ) ;
	/// Assignment to another Map_star. 
	virtual void operator=(const Map_star& ) ;
        
    // Memory management
    // -----------------
    private:
	/// Assignment of the building functions to the member \c Coords
	void set_coord() ;	
    //protected:
	//virtual void reset_coord() ;  ///< Resets all the member \c Coords	

    // Extraction of information
    // -------------------------
    public:
	/// Returns the reference on the Tbl \c alpha 
	const Valeur& get_alpha() const ; 
	const Valeur& get_beta() const ;
	
	
	/// Modifies the value of \f$\alpha\f$ in domain no. \e l 
	void set_alpha(const Tbl& alpha0, int l) ;
	void set_beta(const Tbl& beta0, int l) ;

	void set_alpha(const Valeur& alpha0) ;
	void set_beta(const Valeur& beta0) ;

	/**
	 *  Returns the value of the radial coordinate \e r  for a given
	 *  \f$(\xi, \theta', \phi')\f$ in a given domain. 
	 *	@param l [input] index of the domain
	 *	@param xi [input] value of \f$\xi\f$
	 *	@param theta [input] value of \f$\theta'\f$
	 *	@param pphi [input] value of \f$\phi'\f$
	 *	@return value of \f$r=R_l(\xi, \theta', \phi')\f$
	 */
	virtual double val_r(int l, double xi, double theta, double pphi) const ; 

	/**
	 * Computes the domain index \e l  and the value of \f$\xi\f$ corresponding
	 * to a point given by its physical coordinates \f$(r, \theta, \phi)\f$. 
	 *	@param rr [input] value of \e r 
	 *	@param theta [input] value of \f$\theta\f$
	 *	@param pphi [input] value of \f$\phi\f$
	 *	@param l [output] value of the domain index
	 *	@param xi [output] value of \f$\xi\f$
	 */
	virtual void val_lx(double rr, double theta, double pphi,
			    int& l, double& xi) const ; 
	
	/** Computes the domain index \e l  and the value of \f$\xi\f$ corresponding
	 * to a point given by its physical coordinates \f$(r, \theta, \phi)\f$. 
	 *	@param rr [input] value of \e r 
	 *	@param theta [input] value of \f$\theta\f$
	 *	@param pphi [input] value of \f$\phi\f$
	 *	@param par [] unused by the \c Map_star  version
	 *	@param l [output] value of the domain index
	 *	@param xi [output] value of \f$\xi\f$
	 */
	virtual void val_lx(double rr, double theta, double pphi, 
			    const Param& par, int& l, double& xi) const ; 
		
	/** Returns the value of the radial coordinate \e r  for a given
	 *  \f$\xi\f$ and a given collocation point in \f$(\theta', \phi')\f$ 
	 *   in a given domain. 
	 *	@param l [input] index of the domain
	 *	@param xi [input] value of \f$\xi\f$
	 *	@param j [input] index of the collocation point in \f$\theta'\f$
	 *	@param k [input] index of the collocation point in \f$\phi'\f$
	 *	@return value of \f$r=R_l(\xi, {\theta'}_j, {\phi'}_k)\f$
	 */
	virtual double val_r_jk(int l, double xi, int j, int k) const ; 
	
	/** Computes the domain index \e l  and the value of \f$\xi\f$ corresponding
	 * to a point of arbitrary \e r  but collocation values of \f$(\theta, \phi)\f$ 
	 *	@param rr [input] value of \e r 
	 *	@param j [input] index of the collocation point in \f$\theta\f$
	 *	@param k [input] index of the collocation point in \f$\phi\f$
	 *	@param par [] unused by the \c Map_star  version
	 *	@param l [output] value of the domain index
	 *	@param xi [output] value of \f$\xi\f$
	 */
	virtual void val_lx_jk(double rr, int j, int k, const Param& par,
			       int& l, double& xi) const ; 

	/// Comparison operator (egality)
	virtual bool operator==(const Map& ) const ;  

	/** Computes \f$\partial/ \partial r\f$ of a \c Scalar.
	 *  Note that in the  compactified external domain (CED), it computes
	 *  \f$-\partial/ \partial u = r^2 \partial/ \partial r\f$.
	 *  @param ci [input] field to consider
	 *  @param resu [output] derivative of \c ci 
	 */
	virtual void dsdr (const Scalar& ci, Scalar& resu) const ;

	/** Computes \f$1/(r\sin\theta) \partial/ \partial \phi\f$ of a \c Scalar.
	 *  Note that in the  compactified external domain (CED), the \c dzpuis 
	 *  flag of the output is 2 if the input has \c dzpuis  = 0, and 
	 *  is increased by 1 in other cases.
	 *  @param uu [input] field to consider
	 *  @param resu [output] derivative of \c uu 
	 */
	virtual void srstdsdp (const Scalar&, Scalar&) const ;

	/** Computes \f$1/r \partial/ \partial \theta\f$ of a \c Scalar.
	 *  Note that in the  compactified external domain (CED), the \c dzpuis 
	 *  flag of the output is 2 if the input has \c dzpuis  = 0, and 
	 *  is increased by 1 in other cases.
	 *  @param uu [input] field to consider
	 *  @param resu [output] derivative of \c uu 
	 */
	virtual void srdsdt (const Scalar&, Scalar&) const ;

	/** Computes \f$\partial/ \partial \theta\f$ of a \c Scalar .
	 *  @param uu [input] scalar field 
	 *  @param resu [output] derivative of \c uu 
	 */
	virtual void dsdt (const Scalar&, Scalar&) const ;


	/** Computes \f$1/\sin\theta \partial/ \partial \varphi\f$ of a \c Scalar .
	 *  @param uu [input] scalar field 
	 *  @param resu [output] derivative of \c uu 
	 */
	virtual void stdsdp (const Scalar&, Scalar&) const ;

    // Outputs
    // -------
    public:
	virtual void sauve(FILE* ) const ;	  ///< Save in a file
    
    private:
	virtual ostream& operator>>(ostream &) const ;    ///< Operator >>

	virtual void homothetie (double) ; /// < Not implemented
	virtual void resize (int, double) ;/// < Not implemented
	virtual void adapt (const Cmp&, const Param&, int) ;/// < Not implemented
	virtual void dsdr (const Cmp&, Cmp&) const ;/// < Not implemented
	virtual void dsdxi (const Cmp&, Cmp&) const ;/// < Not implemented
	virtual void dsdxi (const Scalar&, Scalar&) const ;/// < Not implemented
	virtual void dsdradial (const Scalar& uu, Scalar& resu) const ;/// < Not implemented
	virtual void srdsdt (const Cmp&, Cmp&) const ;/// < Not implemented
	virtual void srstdsdp (const Cmp&, Cmp&) const ;/// < Not implemented
	virtual void laplacien (const Scalar&, int, Scalar&) const ;/// < Not implemented
	virtual void laplacien (const Cmp&, int, Cmp&) const ;/// < Not implemented
	virtual void lapang (const Scalar&, Scalar&) const ;/// < Not implemented
	virtual void primr(const Scalar&, Scalar&, bool) const ;/// < Not implemented
	virtual Tbl* integrale (const Cmp&) const ;/// < Not implemented
	virtual void poisson (const Cmp&, Param&, Cmp&) const ;/// < Not implemented
	virtual void poisson_tau (const Cmp&, Param&, Cmp&) const ;/// < Not implemented
	virtual void poisson_falloff(const Cmp&, Param&, Cmp&, int) const ;/// < Not implemented
	virtual void poisson_ylm(const Cmp&, Param&, Cmp&, int, double*) const ;/// < Not implemented
	virtual void poisson_regular (const Cmp&, int, int, double, Param&, Cmp&, Cmp&, Cmp&, 
				      Tenseur&, Cmp&, Cmp&) const ;/// < Not implemented
	virtual void poisson_angu (const Scalar&, Param&, Scalar&, double=0) const ;/// < Not implemented
	virtual void poisson_angu (const Cmp&, Param&, Cmp&, double=0) const ;/// < Not implemented
	virtual Param* donne_para_poisson_vect (Param&, int) const ;/// < Not implemented
	virtual void poisson_frontiere (const Cmp&, const Valeur&, int, int, Cmp&, double = 0., double = 0.) const ;/// < Not implemented
	virtual void poisson_frontiere_double (const Cmp&, const Valeur&, const Valeur&, int, Cmp&) const ;/// < Not implemented
	virtual void poisson_interne (const Cmp&, const Valeur&, Param&, Cmp&) const ;/// < Not implemented
	virtual void poisson2d (const Cmp&, const Cmp&, Param&, Cmp&) const ;/// < Not implemented
	virtual void dalembert (Param&, Scalar&, const Scalar&, const Scalar&, const Scalar&) const ;/// < Not implemented
	// Building functions for the Coord's
	// ----------------------------------
	friend Mtbl* map_star_fait_r(const Map* ) ;
	friend Mtbl* map_star_fait_tet(const Map* ) ;
	friend Mtbl* map_star_fait_phi(const Map* ) ;
	friend Mtbl* map_star_fait_sint(const Map* ) ;
	friend Mtbl* map_star_fait_cost(const Map* ) ;
	friend Mtbl* map_star_fait_sinp(const Map* ) ;
	friend Mtbl* map_star_fait_cosp(const Map* ) ;
	
	friend Mtbl* map_star_fait_x(const Map* ) ;
	friend Mtbl* map_star_fait_y(const Map* ) ;
	friend Mtbl* map_star_fait_z(const Map* ) ;
	
	friend Mtbl* map_star_fait_xa(const Map* ) ;
	friend Mtbl* map_star_fait_ya(const Map* ) ;
	friend Mtbl* map_star_fait_za(const Map* ) ;
	
	friend Mtbl* map_star_fait_xsr(const Map* ) ;
	friend Mtbl* map_star_fait_dxdr(const Map* ) ;
	friend Mtbl* map_star_fait_drdt(const Map* ) ;
	friend Mtbl* map_star_fait_stdrdp(const Map* ) ;
	friend Mtbl* map_star_fait_srdrdt(const Map* ) ;
	friend Mtbl* map_star_fait_srstdrdp(const Map* ) ;
	friend Mtbl* map_star_fait_sr2drdt(const Map* ) ;
	friend Mtbl* map_star_fait_sr2stdrdp(const Map* ) ;
	friend Mtbl* map_star_fait_d2rdx2(const Map* ) ;
	friend Mtbl* map_star_fait_lapr_tp(const Map* ) ;
	friend Mtbl* map_star_fait_d2rdtdx(const Map* ) ;
	friend Mtbl* map_star_fait_sstd2rdpdx(const Map* ) ;
	friend Mtbl* map_star_fait_sr2d2rdt2(const Map* ) ;
};
		Mtbl* map_star_fait_r(const Map* ) ;
		Mtbl* map_star_fait_tet(const Map* ) ;
		Mtbl* map_star_fait_phi(const Map* ) ;
		Mtbl* map_star_fait_sint(const Map* ) ;
		Mtbl* map_star_fait_cost(const Map* ) ;
		Mtbl* map_star_fait_sinp(const Map* ) ;
		Mtbl* map_star_fait_cosp(const Map* ) ;
			
		Mtbl* map_star_fait_x(const Map* ) ;
		Mtbl* map_star_fait_y(const Map* ) ;
		Mtbl* map_star_fait_z(const Map* ) ;
			
		Mtbl* map_star_fait_xa(const Map* ) ;
		Mtbl* map_star_fait_ya(const Map* ) ;
		Mtbl* map_star_fait_za(const Map* ) ;
			
		Mtbl* map_star_fait_xsr(const Map* ) ;
		Mtbl* map_star_fait_dxdr(const Map* ) ;
		Mtbl* map_star_fait_drdt(const Map* ) ;
		Mtbl* map_star_fait_stdrdp(const Map* ) ;
		Mtbl* map_star_fait_srdrdt(const Map* ) ;
		Mtbl* map_star_fait_srstdrdp(const Map* ) ;
		Mtbl* map_star_fait_sr2drdt(const Map* ) ;
		Mtbl* map_star_fait_sr2stdrdp(const Map* ) ;
		Mtbl* map_star_fait_d2rdx2(const Map* ) ;
		Mtbl* map_star_fait_lapr_tp(const Map* ) ;
		Mtbl* map_star_fait_d2rdtdx(const Map* ) ;
		Mtbl* map_star_fait_sstd2rdpdx(const Map* ) ;
		Mtbl* map_star_fait_sr2d2rdt2(const Map* ) ;


			//------------------------------------//
			//         class Map_eps              //
			//------------------------------------//

/*
 * Affine and starlike radial mapping to describe 3D star. \ingroup (map)
 * 
 * The affine radial mapping is the simplest one between the grid coordinates
 * \f$(\xi, \theta', \phi')\f$ and the physical coordinates \f$(r, \theta, \phi)\f$. It comprises a single domain (nucleus, shells to be added in the future)
 * It is defined by \f$\theta=\theta'\f$, \f$\phi=\phi'\f$ and 
 *  \li \f$r=\alpha \xi + \beta\f$, in non-compactified domains, 
 * where \f$\alpha\f$ and \f$\beta\f$ depend upon the angular direction. 
 * 
 *
 */

class Map_eps : public Map_radial{


	// Data :
    // ----
    private:
	/// Array (size: \c mg->nzone*Nt*Np ) of the values of \f$\alpha\f$ in each domain
	double aa, bb, cc ;
	Valeur alpha ;	 
	Valeur beta  ;


	// Constructors, destructor : 
    // ------------------------
    public:
	/**
	 * Standard Constructor
	 * @param mgrille  [input] Multi-domain grid on which the mapping is defined
	 * @param r_limits [input] Array (size: \c Nt*Np ) (only 1 domain at the moment) of the
	 *			   value of \e r  at the boundaries of the various 
	 *			   domains : 
	 *			   \li \c r_limits[l,k,j] : outer boundary of the 
	 *				 domain no. \c l in the angular direction \c j,k 
	 */
	Map_eps(const Mg3d& mgrille, const double* r_limits) ;	
	/**
	 * Standard Constructor with Tbl
	 * @param mgrille  [input] Multi-domain grid on which the mapping is defined
	 * @param r_limits [input] Array (size: \c Nt*Np ) (only 1 domain at the moment) of the
	 *			   value of \e r  at the boundaries of the various 
	 *			   domains : 
	 *			   \li \c r_limits[l,k,j] : outer boundary of the 
	 *				 domain no. \c l in the angular direction \c j,k 
	 */
	Map_eps(const Mg3d& mgrille, double a, double b, double c) ;	
	
	Map_eps(const Map_eps& ) ;      ///< Copy constructor
	Map_eps(const Mg3d&, const string&) ;///< Constructor from a formatted file
	Map_eps(const Mg3d&, FILE* ) ; ///< Constructor from a file (see \c sauve(FILE*) )


	virtual ~Map_eps() ;	      ///< Destructor

	/** Returns the "angular" mapping for the outside of domain \c l_zone.
	 * Valid only for the class \c Map_af.
	 */
	virtual const Map_af& mp_angu(int) const ;

	// Assignment
    // ----------
    public: 
	/// Assignment to an affine mapping.
	virtual void operator=(const Map_af& ) ;
	/// Assignment to another Map_eps. 
	virtual void operator=(const Map_eps& ) ;
        
    // Memory management
    // -----------------
    private:
	/// Assignment of the building functions to the member \c Coords
	void set_coord() ;	
    //protected:
	//virtual void reset_coord() ;  ///< Resets all the member \c Coords	

    // Extraction of information
    // -------------------------
    public:
	/// Returns the reference on the Tbl \c alpha 
	const Valeur& get_alpha() const ; 
	const Valeur& get_beta() const ;
	double get_aa() const { return aa;} ;
	double get_bb() const { return bb;} ;
	double get_cc() const { return cc;} ;
	
	
	/// Modifies the value of \f$\alpha\f$ in domain no. \e l 
	void set_alpha(const Tbl& alpha0, int l) ;
	void set_beta(const Tbl& beta0, int l) ;

	void set_alpha(const Valeur& alpha0) ;
	void set_beta(const Valeur& beta0) ;

	/**
	 *  Returns the value of the radial coordinate \e r  for a given
	 *  \f$(\xi, \theta', \phi')\f$ in a given domain. 
	 *	@param l [input] index of the domain
	 *	@param xi [input] value of \f$\xi\f$
	 *	@param theta [input] value of \f$\theta'\f$
	 *	@param pphi [input] value of \f$\phi'\f$
	 *	@return value of \f$r=R_l(\xi, \theta', \phi')\f$
	 */
	virtual double val_r(int l, double xi, double theta, double pphi) const ; 

	/**
	 * Computes the domain index \e l  and the value of \f$\xi\f$ corresponding
	 * to a point given by its physical coordinates \f$(r, \theta, \phi)\f$. 
	 *	@param rr [input] value of \e r 
	 *	@param theta [input] value of \f$\theta\f$
	 *	@param pphi [input] value of \f$\phi\f$
	 *	@param l [output] value of the domain index
	 *	@param xi [output] value of \f$\xi\f$
	 */
	virtual void val_lx(double rr, double theta, double pphi,
			    int& l, double& xi) const ; 
	
	/** Computes the domain index \e l  and the value of \f$\xi\f$ corresponding
	 * to a point given by its physical coordinates \f$(r, \theta, \phi)\f$. 
	 *	@param rr [input] value of \e r 
	 *	@param theta [input] value of \f$\theta\f$
	 *	@param pphi [input] value of \f$\phi\f$
	 *	@param par [] unused by the \c Map_eps  version
	 *	@param l [output] value of the domain index
	 *	@param xi [output] value of \f$\xi\f$
	 */
	virtual void val_lx(double rr, double theta, double pphi, 
			    const Param& par, int& l, double& xi) const ; 
		
	/** Returns the value of the radial coordinate \e r  for a given
	 *  \f$\xi\f$ and a given collocation point in \f$(\theta', \phi')\f$ 
	 *   in a given domain. 
	 *	@param l [input] index of the domain
	 *	@param xi [input] value of \f$\xi\f$
	 *	@param j [input] index of the collocation point in \f$\theta'\f$
	 *	@param k [input] index of the collocation point in \f$\phi'\f$
	 *	@return value of \f$r=R_l(\xi, {\theta'}_j, {\phi'}_k)\f$
	 */
	virtual double val_r_jk(int l, double xi, int j, int k) const ; 
	
	/** Computes the domain index \e l  and the value of \f$\xi\f$ corresponding
	 * to a point of arbitrary \e r  but collocation values of \f$(\theta, \phi)\f$ 
	 *	@param rr [input] value of \e r 
	 *	@param j [input] index of the collocation point in \f$\theta\f$
	 *	@param k [input] index of the collocation point in \f$\phi\f$
	 *	@param par [] unused by the \c Map_eps  version
	 *	@param l [output] value of the domain index
	 *	@param xi [output] value of \f$\xi\f$
	 */
	virtual void val_lx_jk(double rr, int j, int k, const Param& par,
			       int& l, double& xi) const ; 

	/// Comparison operator (egality)
	virtual bool operator==(const Map& ) const ;  

	/** Computes \f$\partial/ \partial r\f$ of a \c Scalar.
	 *  Note that in the  compactified external domain (CED), it computes
	 *  \f$-\partial/ \partial u = r^2 \partial/ \partial r\f$.
	 *  @param ci [input] field to consider
	 *  @param resu [output] derivative of \c ci 
	 */
	virtual void dsdr (const Scalar& ci, Scalar& resu) const ;

	/** Computes \f$1/(r\sin\theta) \partial/ \partial \phi\f$ of a \c Scalar.
	 *  Note that in the  compactified external domain (CED), the \c dzpuis 
	 *  flag of the output is 2 if the input has \c dzpuis  = 0, and 
	 *  is increased by 1 in other cases.
	 *  @param uu [input] field to consider
	 *  @param resu [output] derivative of \c uu 
	 */
	virtual void srstdsdp (const Scalar&, Scalar&) const ;

	/** Computes \f$1/r \partial/ \partial \theta\f$ of a \c Scalar.
	 *  Note that in the  compactified external domain (CED), the \c dzpuis 
	 *  flag of the output is 2 if the input has \c dzpuis  = 0, and 
	 *  is increased by 1 in other cases.
	 *  @param uu [input] field to consider
	 *  @param resu [output] derivative of \c uu 
	 */
	virtual void srdsdt (const Scalar&, Scalar&) const ;

	/** Computes \f$\partial/ \partial \theta\f$ of a \c Scalar .
	 *  @param uu [input] scalar field 
	 *  @param resu [output] derivative of \c uu 
	 */
	virtual void dsdt (const Scalar&, Scalar&) const ;


	/** Computes \f$1/\sin\theta \partial/ \partial \varphi\f$ of a \c Scalar .
	 *  @param uu [input] scalar field 
	 *  @param resu [output] derivative of \c uu 
	 */
	virtual void stdsdp (const Scalar&, Scalar&) const ;

    // Outputs
    // -------
    public:
	virtual void sauve(FILE* ) const ;	  ///< Save in a file
    
    private:
	virtual ostream& operator>>(ostream &) const ;    ///< Operator >>

	virtual void homothetie (double) ; /// < Not implemented
	virtual void resize (int, double) ;/// < Not implemented
	virtual void adapt (const Cmp&, const Param&, int) ;/// < Not implemented
	virtual void dsdr (const Cmp&, Cmp&) const ;/// < Not implemented
	virtual void dsdxi (const Cmp&, Cmp&) const ;/// < Not implemented
	virtual void dsdxi (const Scalar&, Scalar&) const ;/// < Not implemented
	virtual void dsdradial (const Scalar& uu, Scalar& resu) const ;/// < Not implemented
	virtual void srdsdt (const Cmp&, Cmp&) const ;/// < Not implemented
	virtual void srstdsdp (const Cmp&, Cmp&) const ;/// < Not implemented
	virtual void laplacien (const Scalar&, int, Scalar&) const ;/// < Not implemented
	virtual void laplacien (const Cmp&, int, Cmp&) const ;/// < Not implemented
	virtual void lapang (const Scalar&, Scalar&) const ;/// < Not implemented
	virtual void primr(const Scalar&, Scalar&, bool) const ;/// < Not implemented
	virtual Tbl* integrale (const Cmp&) const ;/// < Not implemented
	virtual void poisson (const Cmp&, Param&, Cmp&) const ;/// < Not implemented
	virtual void poisson_tau (const Cmp&, Param&, Cmp&) const ;/// < Not implemented
	virtual void poisson_falloff(const Cmp&, Param&, Cmp&, int) const ;/// < Not implemented
	virtual void poisson_ylm(const Cmp&, Param&, Cmp&, int, double*) const ;/// < Not implemented
	virtual void poisson_regular (const Cmp&, int, int, double, Param&, Cmp&, Cmp&, Cmp&, 
				      Tenseur&, Cmp&, Cmp&) const ;/// < Not implemented
	virtual void poisson_angu (const Scalar&, Param&, Scalar&, double=0) const ;/// < Not implemented
	virtual void poisson_angu (const Cmp&, Param&, Cmp&, double=0) const ;/// < Not implemented
	virtual Param* donne_para_poisson_vect (Param&, int) const ;/// < Not implemented
	virtual void poisson_frontiere (const Cmp&, const Valeur&, int, int, Cmp&, double = 0., double = 0.) const ;/// < Not implemented
	virtual void poisson_frontiere_double (const Cmp&, const Valeur&, const Valeur&, int, Cmp&) const ;/// < Not implemented
	virtual void poisson_interne (const Cmp&, const Valeur&, Param&, Cmp&) const ;/// < Not implemented
	virtual void poisson2d (const Cmp&, const Cmp&, Param&, Cmp&) const ;/// < Not implemented
	virtual void dalembert (Param&, Scalar&, const Scalar&, const Scalar&, const Scalar&) const ;/// < Not implemented
	// Building functions for the Coord's
	// ----------------------------------
	friend Mtbl* map_eps_fait_r(const Map* ) ;
	friend Mtbl* map_eps_fait_tet(const Map* ) ;
	friend Mtbl* map_eps_fait_phi(const Map* ) ;
	friend Mtbl* map_eps_fait_sint(const Map* ) ;
	friend Mtbl* map_eps_fait_cost(const Map* ) ;
	friend Mtbl* map_eps_fait_sinp(const Map* ) ;
	friend Mtbl* map_eps_fait_cosp(const Map* ) ;
	
	friend Mtbl* map_eps_fait_x(const Map* ) ;
	friend Mtbl* map_eps_fait_y(const Map* ) ;
	friend Mtbl* map_eps_fait_z(const Map* ) ;
	
	friend Mtbl* map_eps_fait_xa(const Map* ) ;
	friend Mtbl* map_eps_fait_ya(const Map* ) ;
	friend Mtbl* map_eps_fait_za(const Map* ) ;
	
	friend Mtbl* map_eps_fait_xsr(const Map* ) ;
	friend Mtbl* map_eps_fait_dxdr(const Map* ) ;
	friend Mtbl* map_eps_fait_drdt(const Map* ) ;
	friend Mtbl* map_eps_fait_stdrdp(const Map* ) ;
	friend Mtbl* map_eps_fait_srdrdt(const Map* ) ;
	friend Mtbl* map_eps_fait_srstdrdp(const Map* ) ;
	friend Mtbl* map_eps_fait_sr2drdt(const Map* ) ;
	friend Mtbl* map_eps_fait_sr2stdrdp(const Map* ) ;
	friend Mtbl* map_eps_fait_d2rdx2(const Map* ) ;
	friend Mtbl* map_eps_fait_lapr_tp(const Map* ) ;
	friend Mtbl* map_eps_fait_d2rdtdx(const Map* ) ;
	friend Mtbl* map_eps_fait_sstd2rdpdx(const Map* ) ;
	friend Mtbl* map_eps_fait_sr2d2rdt2(const Map* ) ;
};
		Mtbl* map_eps_fait_r(const Map* ) ;
		Mtbl* map_eps_fait_tet(const Map* ) ;
		Mtbl* map_eps_fait_phi(const Map* ) ;
		Mtbl* map_eps_fait_sint(const Map* ) ;
		Mtbl* map_eps_fait_cost(const Map* ) ;
		Mtbl* map_eps_fait_sinp(const Map* ) ;
		Mtbl* map_eps_fait_cosp(const Map* ) ;
			
		Mtbl* map_eps_fait_x(const Map* ) ;
		Mtbl* map_eps_fait_y(const Map* ) ;
		Mtbl* map_eps_fait_z(const Map* ) ;
			
		Mtbl* map_eps_fait_xa(const Map* ) ;
		Mtbl* map_eps_fait_ya(const Map* ) ;
		Mtbl* map_eps_fait_za(const Map* ) ;
			
		Mtbl* map_eps_fait_xsr(const Map* ) ;
		Mtbl* map_eps_fait_dxdr(const Map* ) ;
		Mtbl* map_eps_fait_drdt(const Map* ) ;
		Mtbl* map_eps_fait_stdrdp(const Map* ) ;
		Mtbl* map_eps_fait_srdrdt(const Map* ) ;
		Mtbl* map_eps_fait_srstdrdp(const Map* ) ;
		Mtbl* map_eps_fait_sr2drdt(const Map* ) ;
		Mtbl* map_eps_fait_sr2stdrdp(const Map* ) ;
		Mtbl* map_eps_fait_d2rdx2(const Map* ) ;
		Mtbl* map_eps_fait_lapr_tp(const Map* ) ;
		Mtbl* map_eps_fait_d2rdtdx(const Map* ) ;
		Mtbl* map_eps_fait_sstd2rdpdx(const Map* ) ;
		Mtbl* map_eps_fait_sr2d2rdt2(const Map* ) ;

