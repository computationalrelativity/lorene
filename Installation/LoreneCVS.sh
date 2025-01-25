#!/bin/bash

usage () {
    echo "Checkout Lorene from CVS server and inspect CVS history from date"
    echo "Usage: $(basename "$0") [OPTIONS]"
    echo " -h  help"
    echo " -n  HOME_LORENE : new cvs checkout"
    echo " -u  HOME_LORENE : update existing Lorene repo"
    echo " -g  HOME_LORENE : check CVS history from given date"
    echo " -d  YYYY-MM-DD  : date (default is '2024-08-13', date of git import)"
}

cvs_checkout() {
    cd $1
    cvs -d :pserver:anonymous:anonymous@octane.obspm.fr:/cvsroot login
    cvs -z5 -d :pserver:anonymous@octane.obspm.fr:/cvsroot checkout Lorene 
}

cvs_update() {
    cd $1
    cvs update -d
}

# CVS tips
tip="
Useful CVS commands
-------------------
Show commits (-c) from all users (-a) since a specified date:
  cvs history -c <date> -a

Read the log message of a file version:
  cvs log -r <version> <path/to/file>

Diff two versions of a file:
  cvs diff -r <version> -r <version> <path/to/file>

Diff file version relative to two dates:
  cvs diff -D <date1> -D <date2> <path/to/file>
"

# Some info on CoRe git cvsimport 
# https://github.com/computationalrelativity/lorene/wiki
initial_git_import=2024-08-13
latest_git_import=$initial_git_import

# Defaults
NEW=
UPDATE=
HISTORY=
DATE=$initial_git_import

options='hn:u:g:d:'
while getopts $options opt; do
  case "$opt" in
      h)
	  usage && exit 1
	  ;;
      n)
	  NEW=$OPTARG
	  ;;
      u)
	  UPDATE=$OPTARG
	  ;;
      g)
	  HISTORY=$OPTARG
	  ;;
      d)
	  DATE=$OPTARG
	  ;;
      :)
	  echo "Option -$OPTARG requires an argument." >&2
	  usage >&2 && exit 1
	  ;;
      \?)
	  echo "Invalid option: -$OPTARG" >&2
	  usage >&2 && exit 1
	  ;;
  esac
done
shift $((OPTIND-1))

# if [ "$#" -le 0 ]; then
#     usage && exit 1
# fi

if [ ! -z "$NEW" ] ; then
    mkdir $NEW && cvs_checkout $NEW && exit 1
fi

if [ ! -z "$UPDATE" ] ; then
    cd $UPDATE && ./lorene_up && exit 1
    #cvs_update $UPDATE && exit 1
fi

if [ ! -z "$HISTORY" ] ; then
    # Show commits (-c) from all users (-a) since a specified date:
    #   cvs history -c -D 2012-04-01 -a
    echo "Lorene CVS history from $DATE :"
    cd $HISTORY && cvs history -c -D $DATE -a && printf $tip$ && exit 1
fi

