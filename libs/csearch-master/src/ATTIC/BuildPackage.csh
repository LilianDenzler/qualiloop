#!/bin/csh -f
#############################################################################
#      Name: BldPackTmplt.csh                                               #
#  Function: Template for building a package                                #
# Copyright: (C) OXFORD MOLECULAR LTD, 1992                                 #
#---------------------------------------------------------------------------#
#    Author: Dr. D. Walker, Oxford Molecular Ltd                            #
#      Date: 27/10/92                                                       #
#---------------------------------------------------------------------------#
# Arguments: NONE                                                           #
#    Locals:                                                                #
#   Outputs: NONE                                                           #
# Externals: NONE                                                           #
# Exit code: 0 for OK, 2 for ERROR                                          #
#---------------------------------------------------------------------------#
# MODIFICATION RECORD:                                                      #
# DD/MM/YY   Initials   Comments                                            #
# 01/03/93   DW         Disabled debuggable library builds                  #
# 03/03/93   DW         Build CONGEN with -cckr                             #
#############################################################################

# Run the fortran pre-processor
echo "Running the FORTRAN pre-processor ..."
$OML_TOPDIR/DEVTOOLS/RunFPP.csh

# Set up the package name
set package = $cwd:t

### OMUPD DW 03/03/93 Make sure -cckr is in the Irix 4 build flags
if (`uname` == "IRIX") then
   setenv OML_CFLAGS "$OML_CFLAGS -cckr"
endif

# Build the two make files

echo "Building NO DEBUG make file ..."
$OML_TOPDIR/DEVTOOLS/MakeMake.csh CENTRAL $package.mak OFF

### OMUPD DW 01/03/93 Disable the debug builds
# echo "Building DEBUG make file ..."
# $OML_TOPDIR/DEVTOOLS/MakeMake.csh CENTRAL $package.dbgmak ON


# Build using the two make files

make -i -k -f $package.mak
$OML_TOPDIR/DEVTOOLS/AddHistory.csh BUILDNODBG $package

### OMUPD DW 01/03/93 Disable the debug builds
# make -i -k -f $package.dbgmak
# $OML_TOPDIR/DEVTOOLS/AddHistory.csh BUILDDBG $package
