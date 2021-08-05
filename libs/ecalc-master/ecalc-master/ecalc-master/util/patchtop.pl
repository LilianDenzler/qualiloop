#!/usr/bin/perl
#*************************************************************************
#
#   Program:    patchtop
#   File:       patchtop.perl
#   
#   Version:    V1.0
#   Date:       30.09.94
#   Function:   Patch exclusions into an ECalc topology file converted 
#               from XPLOR by xtop.perl. Also adds in the NTER topology
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
#   This program patches an ECalc topology which has been converted from
#   XPLOR format by xtop.perl. A full ECalc topology file is searched for
#   exclusions data which is patched into the new topology file.
#   Finally copies the NTER topology from the reference topology file.
#
#*************************************************************************
#
#   Usage:
#   ======
#
#*************************************************************************
#
#   Revision History:
#   =================
#   V1.0  30.09.94 Original
#
#*************************************************************************
#  main()
#  ------
#  Main program for patching exclusions on converted topology file
#
#  30.09.94 Original    By: ACRM
#
if($#ARGV < 0 || $ARGV[0] eq "-h")
{
   print "\nUsage: patchtop oldtop [newtop] >fixedtop\n\n";
   print "Patches in the atom exclusion data for a topology file (newtop) which\n";
   print "has been converted by xtop from an XPLOR topology file. The oldtop file\n";
   print "is a full ECalc topology file used to find the exclusions.\n\n";
   print "N.B. This is a rather slow and inefficient implementation!\n";
   print "Check the resulting output file for any exclusions not found\n";
   print "in the reference file\n\n";

   exit;
}

#  Open the reference topology file and shift arguments
$reffile = $ARGV[0];
open(REFTOP, $reffile) || die "Unable to open file: $reffile ($!)" ;
close(REFTOP);

shift;

# Flag to ignore NTER
$DoCopy = 1;

while(<>)
{
   ($key, $rest) = split(' ',$_,2);
   
   if($key =~ /RES/i)
   {
      # Store the residue type
      ($residue) = split(' ',$rest);

      # Switch flag if it's an NTER
      if($residue =~ /NTER/i)
      {
         $DoCopy = 0;
      }
      else
      {
         $DoCopy = 1;
      }
      if($DoCopy == 1)
      {
         print;
      }
   }
   elsif($key =~ /ATOM/i)
   {
      if($DoCopy == 1)
      {
         # Process the atom information
         &DoAtom($_, $residue, $reffile);
      }
   }
   else
   {
      if($DoCopy == 1)
      {
         print;
      }
   }
}

&PatchNTER($reffile);

#*************************************************************************
#  sub DoAtom(string buffer, string residue, string reffilename)
#  -------------------------------------------------------------
#  Searches the reference file for the atom specified in buffer in the
#  specified residue type. Outputs the atom specification from the buffer
#  with the exclusions found for that atom name in the reference file
#
#  NOTE: This is terribly slow as it re-reads the reference file for each
#  atom to be processed.
#
#  30.09.94 Original    By: ACRM
#
sub DoAtom
{
   local($buffer, $residue, $reffile) = @_;
   local($filbuff, $key, $rest, $filres, $InRes, $Found);

   $Found = 0;

   # Split the input buffer to get the input atom file spec.
   ($spec) = split(/"/,$buffer);

   # Get the atom name out of the buffer
   ($junk,$atom) = split(' ',$buffer);

   # Set the InRes flag to FALSE
   $InRes = 0;

   # Open the reference file
   open(REFTOP, $reffile) || die "Unable to open file: $reffile ($!)" ;

   # Search through the file for this residue type
   while($filbuff = <REFTOP>)
   {
      ($key, $rest) = split(' ',$filbuff,2);
      
      if($key =~ /RES/i)
      {
         # Get the residue type
         ($filres) = split(' ',$rest);
         if($filres eq $residue)
         {
            $InRes = 1;
         }
         else
         {
            $InRes = 0;
         }
      }
      elsif($InRes == 1)
      {
         if($key =~ /ATOM/i)
         {
            # Split this record into the atom spec and the exclusions
            ($filspec, $exclusions) = split(/"/, $filbuff);

            # Split the atom spec into its parts
            ($junk, $filatom, $filtype, $charge) = split(' ',$filspec);

            # See if we've got the right atom
            if($filatom eq $atom)
            {
               # Write out the spec from the input buffer and the 
               # exclusions from the file
               print $spec, " \"", $exclusions, "\"\n";
               $Found = 1;
               last;
            }
         }
      }
   }

   if($Found == 0)
   {
      print "! No exclusions found in reference for this atom\n";
      print $buffer;
   }

   close(REFTOP);
}

#*************************************************************************
#  sub PatchNTER
#  -------------
#  Adds the NTER information from the reference topology file into the
#  output
#
#  30.09.94 Original    By: ACRM
#
sub PatchNTER
{
   local($reffile) = @_;
   local($InRes);

   # Set the InRes flag to FALSE
   $InRes = 0;

   # Open the reference file
   open(REFTOP, $reffile) || die "Unable to open file: $reffile ($!)" ;

   # Search through the file for this residue type
   while($filbuff = <REFTOP>)
   {
      ($key, $rest) = split(' ',$filbuff,2);
      
      if($key =~ /RES/i)
      {
         # See if residue type is NTER
         if($rest =~ /NTER/i)
         {
            $InRes = 1;
         }
         else
         {
            $InRes = 0;
         }
      }

      if($InRes == 1)
      {
         print $filbuff;
      }
   }

   close(REFTOP);
}




   
