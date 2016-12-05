//======================================================================
// some general-purpose routines for config-file reading
//======================================================================

/*
 *   Copyright (c) 2003 Reinhard Prix
 *
 *   This file is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This file is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with LORENE; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */
 

#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <cctype>

#include "headcpp.h"
#include "utilitaires.h"

namespace Lorene {

/*----------------------------------------------------------------------
 * load_file: does just that plus takes care of memory allocation
 *
 *    !!!!! --> don't forget to free the data after use !!!!!
 *
 * returns NULL on error
 *----------------------------------------------------------------------*/
char *load_file (const char *fname)
{
  FILE *fp;
  char *data;
  size_t size, read_size;
  
  if( (fp = fopen (fname, "r")) == NULL)
    {
      cout << "ERROR: Could not open config-file: '" << fname << "'\n";
      return (NULL);
    }
  
  size = FS_filelength (fp);
  
  // Now read in the raw data
  data = static_cast<char*> (MyMalloc (size+10));
  read_size = fread ( data, 1, size, fp);
  data [read_size] = '\0';  // properly terminate as string!
  
  fclose (fp);
  
  return (data);

} // load_file()

/*----------------------------------------------------------------------
 *  char *read_file_buffered (fname)
 *	
 *  just a buffer for load_file: if fname is NULL or the same as last one, 
 * 	then return buffered copy, otherwise free old one and get a new one
 *  
 * ----------------------------------------------------------------------*/
char *
load_file_buffered (const char *fname)
{
  static char *prev_fname = NULL;
  static char *data = NULL;

  if (!fname)   // fname=NULL: return buffered copy
    return (data);

  if ( prev_fname && !strcmp(prev_fname, fname) )   // fname=prev_name: return buffered copy
    return (data);

  // we're dealing with a new file here, so read it in properly
  if (data)   // free old data
    free(data);
  
  data = load_file (fname);

  return (data);

} // load_file_buffered()


#define ERR -1
#define OK  0
#define FMT_STRING "string"    // reading in strings needs some special treatment

/*----------------------------------------------------------------------
 *  parser for config-files: can read config-variables of the form
 *	VARIABLE [=: \t] VALUE
 *  everything else is ignored as comments
 *
 * RETURN: 	-1 if VARIABLE was not found or could not be read, 
 * 		0 if found&read
 * ----------------------------------------------------------------------*/
int
read_variable (const char *fname, const char *var_name, char *fmt, void *varp)
{
  char *found = NULL;
  char *seek, *pos, *bol;
  int ret;
  int before, after;
  char *data;
  int len;

  if ( (data = load_file_buffered (fname)) == NULL )
    return ERR;

  seek = data;
  while (!found)
    {
      if ( (pos = strstr(seek, var_name)) == NULL)
	break;  // we didn't find anything
      seek = pos + strlen (var_name);       // move data-pointer beyond current position in case we continue

      // make sure we didn't just find a substring:
      if (pos > data)
	before = *(pos-1);
      else
	before = 0;

      after = *(pos+strlen(var_name));

      if ( isalnum(before) || isalnum(after) || (before == '_') || (after == '_') )
	continue;

      // find beginning of this line: bol
      bol = (pos > data) ? pos - 1 : data;
      while ( (bol > data) && (*bol != '\n') ) 
	bol --;
      if ( *bol == '\n' ) bol++;

      // don't allow anything but whitespace before variable-name
      if (pos > bol)
	if ( strspn(bol, " \t") != static_cast<size_t>(pos - bol) )
	  continue;  // must have been a commentary ...

      found = pos;  // ok, that's it
    }

  if (!found)
    {
      cout << "ERROR: variable " << var_name << " was not found in config-file!\n";
      return (ERR);
    }

  found += strlen (var_name);
  
  // skip {space,tab,=,:}
  found += strspn(found, " \t=:");

  // now read the value into the variable
  
  // reading a string needs some special treatment:
  if ( !strcmp(fmt, FMT_STRING) )
    {
      if ( *found == '"')  // skip quotes
	{
	  if ( (pos = strchr(found+1, '"')) == NULL )  // find delimiting quotes
	    {
	      cout << "ERROR: no closing quotes found \n";
	      return (ERR);
	    }
	  found ++;
	} /* if quoted string */
      else
	{
	  if ( (pos = strchr (found, '\n')) == NULL)  // end of line? 
	    pos = data + strlen(data);		// end of file
	} /* if not quoted */

      // NOTE: varp here is supposed to be a pointer to char* !!
      char **cstr = static_cast<char**>(varp);
      len = int(pos - found);  // length of string excluding \0
      (*cstr) = static_cast<char*>(MyMalloc(len+1)); 
      strncpy ((*cstr), found, len);
      (*cstr)[len] = '\0'; 
      ret = 1;  
    } /* if fmt == string */
  else  // the default case is just sscanf...
    ret = sscanf (found, fmt, varp);


  if ( (ret == 0) || (ret == EOF) )
    {
      cout << "WARNING: Variable " << var_name <<" was not readable using the format '"<<fmt<<"'\n";
      return (ERR);
    }

  return (OK);

} // read_variable

/* ----------------------------------------------------------------------
 * specialize to a few common types:
 *----------------------------------------------------------------------*/
int 
read_variable (const char *fname, const char *var_name, int &var)
{
    int ret = read_variable(fname, var_name, const_cast<char*>("%d"), &var);

  cout << "DEBUG: " << var_name << " = " << var <<endl;

  return (ret);
}

int 
read_variable (const char *fname, const char *var_name, bool &var)
{
  int buf;
  int ret = read_variable(fname, var_name, const_cast<char*>("%d"), &buf);

  var = static_cast<bool>(buf);

  cout << "DEBUG: " << var_name << " = " << var <<endl;

  return (ret);
}

int 
read_variable (const char *fname, const char *var_name, double &var)
{
    int ret = read_variable(fname, var_name, const_cast<char*>("%lf"), &var);

  cout << "DEBUG: " << var_name << " = " << var <<endl;

  return (ret);
}

int
read_variable (const char *fname, const char *var_name, char **str)
{
  char *cstr;

  if (*str != NULL)
    {
      cout << "ERROR: return-string needs to be NULL in read_variable()\n";
      return (ERR);
    }

  int ret = read_variable(fname, var_name, const_cast<char*>(FMT_STRING), &cstr);

  if ((ret == OK) && cstr)
    *str = cstr;

  cout << "DEBUG: " << var_name << " = " << *str <<endl;

  return (ret);

}



/*----------------------------------------------------------------------
 * FS_filelength().. (taken from quake2)
 * 		contrary to stat() this fct is nice and portable, 
 *----------------------------------------------------------------------*/
int
FS_filelength (FILE *f)
{
  int		pos;
  int		end;

  pos = int(ftell (f));
  fseek (f, 0, SEEK_END);
  end = int(ftell (f));
  fseek (f, pos, SEEK_SET);
  
  return end;
}

/*@Function============================================================
@Desc: This function works like malloc, except that it also checks for
       success and terminates in case of "out of memory", so we dont
       need to do this always in the code.

@Ret: 
* $Function----------------------------------------------------------*/
void *
MyMalloc (long bytes)
{
  void *Mptr = NULL;

  // make Gnu-compatible even if on a broken system:
  if (bytes == 0)
    bytes = 1;

  if ((Mptr = calloc (1, bytes)) == NULL)
    {
      cout << "MyMalloc("<< bytes << ") did not succeed! Terminating...\n";
      exit (-1);
    }

  return (Mptr);

} // void* MyMalloc

}
