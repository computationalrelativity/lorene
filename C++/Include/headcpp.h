#ifdef OBSOLETE_HEADERS

#include <iomanip.h>
#include <iostream.h>
#include <fstream.h>

#else

#include <iomanip>
using std::setw ;
using std::hex ;
#include <iostream>
using std::ostream ;
using std::istream ;
using std::cout ;
using std::cin ;
using std::cerr ;
using std::clog ;
using std::endl ;
using std::flush ;
using std::ios ;
#include <fstream>
using std::ifstream ;
using std::ofstream ;

#endif

