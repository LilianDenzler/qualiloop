
ProFit V3.3
===========

The ultimate protein least squares fitting program!
---------------------------------------------------


> NOTE: This document is rather out of date - it does not describe
> all the current features of ProFit. You should look at the
> ProFit manual at http://www.bioinf.org.uk/software/profit/ProFit.pdf

-----------------------------------------------------------------------

```
(c) 1992-2020 Prof. Andrew C.R. Martin, SciTech Software,
    abYinformatics, University College London.

    Some features added and bugs fixed while working for the
    University of Reading under consultancy to Inpharmatica, Ltd.

    Features in V3.0 written by Dr. Craig T. Porter, UCL, 2008-2009
    funded by the BBSRC.
```

-----------------------------------------------------------------------

   To access the ProFit download page, got to
http://www.bioinf.org.uk/software/profit/

   ProFit (pronounced Pro-Fit, not profit!) is designed to be the
ultimate least squares fitting program and is written to be as easily
portable between systems as possible. It performs the basic function
of fitting one protein structure to another, but allows as much
flexibility as possible in this procedure. Thus one can specify
subsets of atoms to be considered, specify zones to be fitted by
number, sequence, or by sequence alignment.

   The program will output an RMS deviation (RMSD) and optionally the
fitted coordinates. RMS deviations may also be calculated without
actually performing a fit. Zones for calculating the RMSD can be
different from those used for fitting.

   The interface is command driven, but may offer a graphical 
environment at a later date.

   This document gives an overview of the features and command
language, but isn't regularly kept up-to-date to reflect new features.
There is a full comprehensive manual describing all the command in the
'doc' subdirectory of the full distribution.

Features:
1. Portability
2. Specify atom subsets
3. Specify zones:
  * Numerically
  * By sequence
  * By auto sequence alignment
4. Output RMSD over
  * Fitted region
  * Any other region
5. Output fitted coordinates
6. Help facility

(1) Portability
---------------

The code is written in highly portable ANSI-C. A future feature may
be a Tcl/Tk-based graphical interface, but this will be separate
from the main program. Earlier version of the code supported
windowing, but this has been dropped from the release version in
favour of a separate GUI layer which will call ProFit with a 
control script.

(2) Atom Subsets
----------------

These may be specified for both fitting and RMSD calculation.
e.g.  `ATOMS N,CA,C,O`

(3a) Numeric zones
------------------

These are specified as follows:

* `ZONE *`            Fits all residues in both structures, clearing
                      other zone specifications.                     
* `ZONE *:1-100`      Fits all in Reference with 1-100 in Mobile     
* `ZONE 10-20`        Fits 10-20 in Reference with 10-20 in Mobile   
* `ZONE 10-20:30-40`  Fits 10-20 in Reference with 30-40 in Mobile   
* `ZONE -10:50-59`    Fits 1-10 in Reference with 50-59 in Mobile    

If you have multiple chains, the chain name may preceed the residue
number (e.g. `ZONE L24-L24`); if unspecified, the first chain is
assumed. If the PDB file has insertions, the zone may include an
insert code by placing it after the residue number (e.g. L25B or
50C).

Optionally, the chain name may be separated from the residue number by
a full stop. Using the full stop also makes the statement
case-sensitive.  In practice, the full stop separator is used with
either numeric or lowercase chain names.

(3b) Sequence zone specification
--------------------------------

Typing a sequence will search for this seq in both structures and fit
these regions. 

* `ZONE CAR:VNS`        Fits the first occurence of CAR in Reference with
                        the first occurence of VNS in Mobile                
* `ZONE CAR,10:VNS,10`  Fits 10 residues from the first occurence of CAR
                        in Reference with the first occurence of VNS
                        in Mobile
* `ZONE CAR,5/2`        Fits 5 residues from the second occurence of
                        CAR in both structures                              
* `ZONE 24-34:EIR,11`   Fits 24-34 in Reference with 11 residues from
                        the first occurence of EIR in Mobile                

(3c) Auto sequence alignment
----------------------------

Needleman & Wunsch sequence alignment is provided to perform 
fitting on sequence aligned residues. The keyword `ALIGN` sets up zones
from the matched parts of the sequence alignment. Alternatively, one
may read a sequence alignment from a file using the `READALIGN`
keyword.

(4a) RMSD over fitted region
----------------------------

The program will report the RMS deviation over the fitted region once 
fitting has been performed with the command `FIT`. The `RATOMS` keyword 
has the same syntax as the `ATOMS` keyword and will cause the specified 
atoms to be included in the RMSD calculation.

(4b) RMSD over other regions
----------------------------

The keyword `RZONE` has the same syntax as the keywords `ZONE`, but causes
the RMSD to be calculated over the specified residue zone.

(5) Output fitted coordinates
-----------------------------

* `WRITE <file>`  causes the fitted coordinates to be written to file.

(6) Help facility
-----------------

A comprehensive VAX/VMS style help facility is provided. Type `HELP`
once in ProFit to enter the help facility.

Running ProFit
==============

```
profit [-h] [reference.pdb[%n] mobile.pdb[%n]]
```

Optionally you may specify the reference and mobile PDB files on
the command line. Otherwise you may use the `REFERENCE` and `MOBILE`
commands once in the program.

The optional `%n` added to a filename specifies the MODEL number for
PDB files containing multiple models (from V3.3 onwards).

The `-h` flag causes the program to read HETATM records from the
PDB files. The default is to ignore them. This may be changed 
using the `(NO)HETATOMS` commands.

Keywords
========

Basic commands
--------------

* `REFERENCE file.pdb` Specify the reference structure
* `MOBILE file.pdb`    Specify the structure to be fitted
* `FIT`                Causes the structures to be fitted on the current
                       residue/atom selection. The RMS deviation will be
                       reported over the fitted selection.

Specifying atoms
----------------

* `ATOMS atm[,atm]...`      Specify atoms to be considered in fitting.
* `BVALUE cutoff [REF|MOB]` Ignore atoms (in fitting and RMSD calculation)
                            if B-value greater than specified cutoff.
                            Optional `REF` or `MOB` causes calculation to be
                            restricted to the specified structure. A
                            negative cutoff may be used to ignore atoms
                            with B-values less than the absolute
                            (i.e. unsigned) value of the cutoff.

Specifying zones
----------------

* `ZONE CLEAR|((*|(X...[,n][/m])|(j[-k]))[:(*|(X...[,n][/m])|(j[-k]))])`
                            Specify the zone to be fitted either
                            numerically, or as a sequence. Repeating the
                            specification after a colon refers to the
                            mobile structure.

*X* represents a residue type or *?* wildcard.

*n* represents a number of residues.

*m* represents the *m*'th occurence of the sequence.

*j* and *k* represent residue numbers.


* `DELZONE CLEAR|((*|(X...[,n][/m])|(j[-k]))[:(*|(X...[,n][/m])|(j[-k]))])`
                                    Delete fitting zone. This command uses
                                    the same syntax as `ZONE`.

* `SETCENTRE CLEAR|(*|j[:j])`       Specifies a single residue as the
                                    centre of fitting.

* `NUMBER (RESIDUE|SEQUENTIAL)`     Specify zones based on residue numbering
                                    or sequential numbering.

* `ALIGN`                           Perform Needleman & Wunsch sequence alignment.
                                    Set the fitting zones based on the
                                    sequence aligned residues.

* `GAPPEN penalty [extend_penalty]` Specify gap penalty and extension penalty
                                    for `ALIGN` command (default: 10 and 2
                                    respectively)

* `READALIGNMENT file.pir`          Read a PIR alignment file and set zones
                                    from that.

Calculating RMSD over different atoms/zones
-------------------------------------------

* `RATOMS  atm[,atm]...`            Specify atoms to be considered in RMSD calculation.
                                    Then reports the RMSD.
* `RZONE CLEAR|((*|(X...[,n][/m])|(j[-k]))[:(*|(X...[,n][/m])|(j[-k]))])`
                                    Specify zone to be considered in RMSD calculation.
                                    Then reports the RMSD.
* `DELRZONE CLEAR|((*|(X...[,n][/m])|(j[-k]))[:(*|(X...[,n][/m])|(j[-k]))])`
                                    Remove zone for RMSD calculation.

Obtaining output
----------------

* `WRITE file.pdb`       Writes the fitted mobile coordinates
* `MATRIX`               Displays the rotation and translation matrices 
                         for fitting.
* `STATUS`               Shows information on the current selections.
* `RMS`                  Redisplays the RMSD calculated over the current
                         zones.
* `RESIDUE [filename]`   Display a by-residue RMSD over current `RATOMS`
                         and `RZONES` (`ATOMS` and `ZONES` if `RATOMS` and `RZONES`
                         are not specified)
* `NFITTED`              Report number of atom pairs fitted

Reading non-protein atoms
-------------------------

* `HETATOMS`             Read HETATOMS (default if started with -h)
* `NOHETATOMS`           Do not read HETATOMS (default unless started with -h)

Weighting the fit
-----------------

* `WEIGHT`               Weight the fitting by the B-value
* `BWEIGHT`              Weight the fitting by the inverse of the B-value
* `NOWEIGHT`             No weighting (default)

Miscellaneous
-------------

* `HELP [keywd]`         Gives help on keywd if specified; otherwise
                         lists valid keywords.
* `IGNOREMISSING`        Ignore missing atoms during fitting
* `NOIGNOREMISSING`      Restore default of generating an error if there are
                         any missing atoms
* `QUIT`                 Exits the program.


Not yet implemented:
====================

* `GRAPHIC`              Start graphical alignment
