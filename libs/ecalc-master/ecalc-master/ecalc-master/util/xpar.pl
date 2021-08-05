#!/usr/bin/perl
#*************************************************************************
#
#   Program:    xpar
#   File:       xpar.perl
#   
#   Version:    V1.0
#   Date:       29.09.94
#   Function:   Convert XPLOR parameter file to ECalc format
#   
#   Copyright:  (c) Dr. Andrew C. R. Martin 1994
#   Author:     Dr. Andrew C. R. Martin
#   Address:    Biomolecular Structure & Modelling Unit,
#               Department of Biochemistry & Molecular Biology,
#               University College,
#               Gower Street,
#               London.
#               WC1E 6BT.
#   Phone:      (Home) +44 (0372) 275775
#   EMail:      INTERNET: amartin@scitec.adsp.sub.org
#                         martin@bsm.bioc.ucl.ac.uk
#               UUCP:     ...{uunet|rutgers}!cbmehq!cbmuk!scitec!amartin
#               JANET:    martin@uk.ac.ucl.bioc.bsm
#               
#*************************************************************************
#
#   This program is not in the public domain, but it may be copied
#   according to the conditions laid out in the accompanying file
#   COPYING.DOC
#
#   The code may be modified as required, but any modifications must be
#   documented so that the person responsible can be identified. If 
#   someone else breaks this code, I don't want to be blamed for code 
#   that does not work! 
#
#   The code may not be sold commercially or included as part of a 
#   commercial product except as described in the file COPYING.DOC.
#
#*************************************************************************
#
#   Description:
#   ============
#
#*************************************************************************
#
#   Usage:
#   ======
#
#   xpar [xplor.par] >ecalc.par
#
#   Converts an XPLOR parameter file to ECalc format
#
#   NOTE: You must edit in the atom masses and radii which are set to 
#   0.000 and 1.500 respectively. For ECalc, the masses are ignored, but
#   the radii are required. They are found in the XPLOR topology file.
#
#*************************************************************************
#
#   Revision History:
#   =================
#   V1.0  29.09.94 Original
#
#*************************************************************************
#  main()
#  ------
#  Main program for converting XPLOR parameter file to ECalc format
#
#  29.09.94 Original    By: ACRM
#
if($ARGV[0] eq "-h")
{
   print "\nxpar V1.0 (c) Andrew C.R. Martin, UCL\n\n";
   print "Usage: xpar [xplor.par] >ecalc.par\n";
   print "Converts an XPLOR parameter file to ECalc format\n\n";
   print "NOTE: You must edit in the atom masses and radii\n";
   print "ECalc actually ignores the masses, but the radii are required\n";
   print "and may be found in the XPLOR topology file.\n\n";

   exit;
}

$DoCommands = 1;
$NBondNum   = 0;

print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
print "!            YOU MUST EDIT IN THE ATOM MASSES & RADII            \n";
print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";

while(<>)
{
   ($key, $rest) = split(' ',$_,2);


   if($key =~ /NBONDS/i)
   {
      $DoCommands = 0;
   }
   
   if($DoCommands == 1)
   {
      if($key eq "REMARKS")
      {
         print "!", $_;
         
      }
      elsif(substr($key,0,1) eq "!")
      {
         print;
      }
      elsif($key =~ /NONB/i)
      {
         # Case must come *before* BOND
         # Split the record properly
         ($key, $atom, $eps, $sigma) = split(' ');
         printf "NONBOND %s %s 0.000 %s %s 1.500\n", 
                $atom, $NBondNum, $eps, $sigma;

         $NBondNum++;
      }
      elsif($key =~ /BOND/i)
      {
         print;
      }
      elsif($key =~ /ANGL/i)
      {
         print;
      }
      elsif($key =~ /DIHE/i)
      {
         print "TORSION ", $rest;
      }
      elsif($key =~ /IMPR/i)
      {
         print "IMPROPER ", $rest;
      }
   }
   else
   {
      if($key =~ /END/i)
      {
         $DoCommands = 1;
      }
   }
}


print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
print "!            YOU MUST EDIT IN THE ATOM MASSES & RADII            \n";
print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";



