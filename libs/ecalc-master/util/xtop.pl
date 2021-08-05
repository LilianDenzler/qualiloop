#!/usr/bin/perl
#*************************************************************************
#
#   Program:    xtop
#   File:       xtop.perl
#   
#   Version:    V1.0
#   Date:       29.09.94
#   Function:   Convert XPLOR topology file to ECalc format
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
#   Usage: xtop [xplor.top] >ecalc.top
#
#   Converts an XPLOR topology file to ECalc format
#
#   NOTE: You must edit in the ATOM exclusions and must add to the
#   topology for the NTER residue as the XPLOR NTER specification is 
#   incomplete for the requirements of ECalc
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
#  Main program for conversion of XPLOR topology files to ECalc format
#
#  29.09.94 Original    By: ACRM
#
if($ARGV[0] eq "-h")
{
   print "\nxtop V1.0 (c) Andrew C.R. Martin, UCL\n\n";
   print "Usage: xtop [xplor.top] >ecalc.top\n";
   print "Converts an XPLOR topology file to ECalc format\n\n";
   print "NOTE: You must edit in the ATOM exclusions and must add to\n";
   print "the NTER topology as the XPLOR NTER specification is\n\n";
   print "incomplete for the requirements of ECalc\n\n";

   exit;
}

print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
print "!            YOU MUST EDIT IN THE ATOM EXCLUSIONS                \n";
print "!                  AND ADD TO NTER TOPOLOGY                      \n";
print "!                                                                \n";
print "!          ALSO CHECK ATOM ORDER IS N,H,CA,C,O, s/c              \n";
print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";


while(<>)
{
   ($key, $rest) = split(' ',$_,2);

   if($key eq "REMARKS")
   {
      print "!", $_;
      
   }
   elsif(substr($key,0,1) eq "!")
   {
      print;
   }
   elsif($key =~ /RESI/i)
   {
      # Flag for adding IMPROPER records
      $ImpDone = 0;
      $TorDone = 0;

      # Resplit to get residue name
      chop;
      ($key, $restype) = split(' ',$_);
      
      print "RESIDUE ", $restype, " 0.0\n";
   }
   elsif($key =~ /ATOM/i)
   {
      &DoAtom($_);
   }
   elsif($key =~ /BOND/i)
   {
      print;
   }
   elsif($key =~ /IMPR/i)
   {
      if($ImpDone == 0)
      {
         $ImpDone = 1;
         print " IMPROPER N    -C   CA   H\n";
         print " IMPROPER C    CA   +N   O\n";
      }
      print;
   }
   elsif($key =~ /DIHE/i)
   {
      if($TorDone == 0)
      {
         $TorDone = 1;

         print " TORSION -C   N    CA   C\n";
         print " TORSION N    CA   C    +N\n";
         print " TORSION CA   C    +N   +CA\n";
      }
      chop($rest);
      print " TORSION ", $rest, "\n";
   }
   elsif($key =~ /DONO/i)
   {
      chop;
      
      print $_, "    X    X\n";
   }
   elsif($key =~ /ACCE/i)
   {
      ($atom) = split(' ', $rest);
      
      print " ACCEPTOR ", $atom, "\n";
   }
   elsif($key =~ /END/i)
   {
      if($ImpDone == 0)
      {
         $ImpDone = 1;
         print " IMPROPER N    -C   CA   H\n";
         print " IMPROPER C    CA   +N   O\n";
      }

      if($TorDone == 0)
      {
         $TorDone = 1;

         print " TORSION -C   N    CA   C\n";
         print " TORSION N    CA   C    +N\n";
         print " TORSION CA   C    +N   +CA\n";
      }
   }
}

print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
print "!            YOU MUST EDIT IN THE ATOM EXCLUSIONS                \n";
print "!                  AND ADD TO NTER TOPOLOGY                      \n";
print "!                                                                \n";
print "!          ALSO CHECK ATOM ORDER IS N,H,CA,C,O, s/c              \n";
print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";


#*************************************************************************
#  sub DoAtom(string buffer)
#  -------------------------
#  Subroutine to extract and print required information for an ATOM recor
#
#  29.09.94 Original    By: ACRM
#
sub DoAtom
{
   local($buffer) = @_;
   local($junk1,$atom,$type,$charge);

   ($junk0,$junk1,$atom,$junk2,$type,$junk3,$charge) = split(/[ =\s]+/,$buffer);

   print " ATOM ", $atom, "   ", $type, "   ", $charge, "   \" \"\n";
   
}
