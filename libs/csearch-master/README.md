CSearch
=======

CSearch is based on Bob Bruccoleri's CONGEN code (available from
http://www.congenomics.com)

Around 1991, I extracted the conformational search code from CONGEN
(which was based on CHARMM) for Oxford Molecular Ltd., to use as part
of their AbM antibody modelling software (which was based on my DPhil
research). Bob Bruccoleri gave us permission to do this since CHARMM
was commercial while his conformational search additions were owned by
him. Consequently I rewrote all the energy calculation code from
scratch.

Currently this code is incomplete and doesn't compile. The original
was written on VAX/VMS. Around 1994 I ported this to work with f2c
under Linux, but it doesn't currently work with gfortran. I need
another day or so's work to get this up and running!
