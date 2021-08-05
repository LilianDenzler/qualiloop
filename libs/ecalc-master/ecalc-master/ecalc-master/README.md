ECalc
=====

(c) 1994-2020 UCL, Andrew C.R. Martin
=====================================

ECalc is a program designed to calculate the energy of a protein. It
may also be used to calculate the energies of a number of conformations 
of a region within a protein. The calculations may be restricted to
the region of interest (though non-bonded contacts with other parts of
the protein will also be considered). In addition, other parts of the
protein may be specified which should be excluded from consideration
in the non-bonded contact calculations. Parts of the potential may be
switched on or off and scaled with respect to one another. This
allows, for example, the simple calculation of just the van der Waals
repulsion or just the bond energy.

The program may is run using a control file which can be created
automatically using a Tcl/Tk interface.

See the doc/ecalc.pdf file for more information.

