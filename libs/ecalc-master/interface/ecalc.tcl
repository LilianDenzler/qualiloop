#!/usr/bin/wish -f
#*************************************************************************
#
#   Program:    ecalc interface
#   File:       ecalc.tcl
#   
#   Version:    V1.4
#   Date:       23.05.95
#   Function:   Write a control file for ECalc
#   
#   Copyright:  (c) Dr. Andrew C. R. Martin 1994-5
#   Author:     Dr. Andrew C. R. Martin
#   Address:    Biomolecular Structure & Modelling Unit,
#               Department of Biochemistry & Molecular Biology,
#               University College,
#               Gower Street,
#               London.
#               WC1E 6BT.
#   Phone:      (Home) +44 (0)1372 275775
#               (Work) +44 (0)171 387 7050 X 3284
#   EMail:      martin@biochem.ucl.ac.uk
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
#*************************************************************************
#
#   Revision History:
#   =================
#   V0.0  09.09.94 Original
#   V0.1  12.09.94 Added advanced options button
#   V0.2  13.09.94 Added load and save configuration
#   V0.3  14.09.94 Added disulphide button
#   V0.4  15.09.94 Added regrid option
#   V0.5  16.09.94 Added zones button 
#   V0.6  22.09.94 Added ignores button
#                  Added <return> binding for nkeep entry
#   V0.7  29.09.94 Calls local.tcl to get local data
#   V1.0  29.09.94 First release version
#   V1.1  24.10.94 Added residue component to potential
#   V1.2  Skipped
#   V1.3  Skipped
#   V1.4  23.05.95 Added RELAX options
#
#*************************************************************************

##########################################################################
#   proc InitWM 
#   -----------
#   Initialise the window manager, set up the name bar, etc.
#
#   22.09.94 Original extracted from main source code   By: ACRM
#   24.10.94 V1.1
#   23.05.94 V1.4
#
proc InitWM { } {
    wm title . "ECalc V1.4 (c) 1994-5, Dr. Andrew C.R. Martin, UCL"
    wm iconname . "ECalc"
}

##########################################################################
#   proc InitVariables
#   ------------------
#   Initialise all variables
#
#   22.09.94 Original extracted from main source code   By: ACRM
#   24.10.94 Added potresidue, dispresidue, residuescale
#   23.05.95 Added Relax, shakeTol and vdwTol
#
proc InitVariables { } {
    #  Reference global variables
    global conffile pdbfile outfile
    global potbond potangle pottor potimp pothbond potvdwa potvdwr \
           potelect potresidue
    global dispbond dispangle disptor dispimp disphbond dispvdwa   \
           dispvdwr dispelect dispresidue
    global bondscale anglescale torscale impscale hbondscale vdwascale \
           vdwrscale electscale residuescale
    global nkeep
    global paramfile rtopfile Cutgrid Cutnonbond
    global ShowTimings ShowParams ShowDebug ShowRTop
    global ConstDielectric eta 
    global Cutdiston Cutangon Cutdistoff Cutangoff
    global AutoSS Regrid
    global Relax shakeTol vdwTol

    set potbond     1
    set potangle    1
    set pottor      1
    set potimp      1
    set pothbond    1
    set potvdwa     1
    set potvdwr     1
    set potelect    1
    set potresidue  1
    
    set dispbond    0
    set dispangle   0
    set disptor     0
    set dispimp     0
    set disphbond   0
    set dispvdwa    0
    set dispvdwr    0
    set dispelect   0
    set dispresidue 0
    
    set bondscale    1.0
    set anglescale   1.0
    set torscale     1.0
    set impscale     1.0
    set hbondscale   1.0
    set vdwascale    1.0
    set vdwrscale    1.0
    set electscale   1.0
    set residuescale 1.0
    
    set nkeep      1
    
    # Variables used by AdvancedOptions procedure    
    set paramfile ""
    set rtopfile  ""
    set Cutgrid   30.0
    set Regrid       0
    set Cutnonbond 8.0
    
    set ShowTimings 0
    set ShowParams  0
    set ShowDebug   0
    set ShowRTop    0
    
    set ConstDielectric 0
    set eta 50.0
    
    set Cutdiston 4.0
    set Cutangon  90.0
    set Cutdistoff 5.0
    set Cutangoff 90.0

    set Relax 0
    set shakeTol 0.001
    set vdwTol 0.5
    
    # Variables used by Disulphides procedure    
    set AutoSS    1
}

##########################################################################
#  proc quit button
#  ----------------
#  The procedure is used to end the program.
#  It is invoked with either the quit or run parameter
#
#  09.09.94 Original    By: ACRM
#  14.09.94 Sets the focus back to one of the text entries, so that
#           a FocusOut from the run button doesn't occur after the
#           button has been destroyed
#
proc quit button {
   if {$button == "run"} { 
      WriteControl ecalc.dat
      RunECalc     ecalc.dat
   }

   focus .top.pdb.entry
   destroy .
}

##########################################################################
#  proc relax relbutton activecolour
#  ---------------------------------
#  The procedure is used when the relax button is pressed
#
#  23.05.95 Original    By: ACRM
#
proc relax { relbutton activecolour } {
    global Relax

    if {$Relax == 0} {
        set Relax 1
        $relbutton configure -relief sunken 
        $relbutton configure -background $activecolour
    } else {
        set Relax 0
        $relbutton configure -relief raised
        $relbutton configure -background bisque
    }
        
}

##########################################################################
#  proc WriteControl filename
#  --------------------------
#  Write the control file from the values of the global variables
#
#  12.09.94 Original    By: ACRM
#  14.09.94 Added writing of disulphides
#  15.09.94 Added regrid
#  16.09.94 Added zones
#  22.09.94 Added ignores
#  24.10.94 Added potresidue, dispresidue, residuescale
#  23.05.95 Added Relax, shakeTol and vdwTol
#
proc WriteControl filename {
#  Reference global variables
   global conffile pdbfile outfile
   global potbond  potangle  pottor  potimp  pothbond  potvdwa  potvdwr \
          potelect potresidue
   global dispbond dispangle disptor dispimp disphbond dispvdwa dispvdwr \
          dispelect dispresidue
   global bondscale anglescale torscale impscale hbondscale vdwascale \
          vdwrscale electscale residuescale
   global nkeep
   global paramfile rtopfile Cutgrid Cutnonbond
   global ShowTimings ShowParams ShowDebug ShowRTop
   global ConstDielectric eta 
   global Cutdiston Cutangon Cutdistoff Cutangoff
   global AutoSS Regrid
   global Relax shakeTol vdwTol

   set file [open $filename w]

#  Write the filenames
   if {$pdbfile  != ""} { puts $file [format "PDBFILE %s" $pdbfile]   }
   if {$conffile != ""} { puts $file [format "CONFFILE %s" $conffile] }
   if {$outfile != ""}  { puts $file [format "OUTFILE %s" $outfile] }

#  Now the potential
   puts $file "POTENTIAL"
   if {$potbond    == 1} { 
      puts $file [format "   BONDS     %f" $bondscale]
   }
   if {$potangle   == 1} { 
      puts $file [format "   ANGLES    %f" $anglescale]
   }
   if {$pottor     == 1} { 
      puts $file [format "   TORSIONS  %f" $torscale]  
   }
   if {$potimp     == 1} { 
      puts $file [format "   IMPROPERS %f" $impscale]  
   }
   if {$pothbond   == 1} { 
      puts $file [format "   HBONDS    %f" $hbondscale]
   }
   if {$potvdwa    == 1} { 
      puts $file [format "   VDWA      %f" $vdwascale] 
   }
   if {$potvdwr    == 1} { 
      puts $file [format "   VDWR      %f" $vdwrscale] 
   }
   if {$potelect   == 1} { 
      puts $file [format "   ELECTR    %f" $electscale]
   }
   if {$potresidue == 1} { 
      puts $file [format "   RESIDUE   %f" $residuescale]
   }
   puts $file "END"

#  Now the display
   set sum [expr $dispbond+$dispangle+$disptor+$dispimp+$dispresidue]
   set sum [expr $sum+$disphbond+$dispvdwa+$dispvdwr+$dispelect]
   if {$sum != 0} { puts $file "DISPLAY" }
   if {$dispbond    == 1} { puts $file "   BONDS"     }
   if {$dispangle   == 1} { puts $file "   ANGLES"    }
   if {$disptor     == 1} { puts $file "   TORSIONS"  }
   if {$dispimp     == 1} { puts $file "   IMPROPERS" }
   if {$disphbond   == 1} { puts $file "   HBONDS"    }
   if {$dispvdwa    == 1} { puts $file "   VDWA"      }
   if {$dispvdwr    == 1} { puts $file "   VDWR"      }
   if {$dispelect   == 1} { puts $file "   ELECTR"    }
   if {$dispresidue == 1} { puts $file "   RESIDUE"   }
   if {$sum != 0} { puts $file "END" }

#  And the number of confs to keep
   puts $file [format "CACHE %d" $nkeep]

#  Advanced options
   if {$paramfile != ""} { puts $file [format "PARAMFILE %s" $paramfile] }
   if {$rtopfile  != ""} { puts $file [format "RTOPFILE %s"  $rtopfile]  }

   puts $file [format "GRIDCUT %f" $Cutgrid]
   puts $file [format "REGRID  %d" $Regrid]
   puts $file [format "NONBONDCUT %f" $Cutnonbond]
   if {$ShowTimings == 1} { puts $file "SHOWTIMINGS" }
   if {$ShowParams  == 1} { puts $file "SHOWPARAMS" }
   if {$ShowRTop    == 1} { puts $file "SHOWRTOP" }

   if {$ShowDebug   == 1} { puts $file "DEBUG" }

   if {$ConstDielectric == 1} { 
       puts $file "CONSTDIELECTRIC" 
       puts $file [format "ETA %f" $eta]
   } else { 
       puts $file "DISTDIELECTRIC"  
   }

   puts $file [format "CUTONHB %f"   $Cutdiston]
   puts $file [format "CUTONANG %f"  $Cutangon]
   puts $file [format "CUTOFFHB %f"  $Cutdistoff]
   puts $file [format "CUTOFFANG %f" $Cutangoff]

   if {$Relax == 1} {
       puts $file "RELAX"
       puts $file [format "TOL SHAKE %f" $shakeTol]
       puts $file [format "TOL VDW %f"   $vdwTol]
   }

#  Disulphides
   if {$AutoSS == 0} {
      puts $file "DISULPHIDES OFF"
   } elseif {$AutoSS == 2} {
      WriteDisulphides $file
   }
   # $AutoSS==1 is the default auto disulphide search

#  Zones
   WriteZones $file

#  Zones
   WriteIgnores $file

   close $file
}

##########################################################################
#   proc RunECalc file
#   ------------------
#   Run the main ecalc program using the control file we have generated
#
#   12.09.94 Original   By: ACRM
#   29.09.94 Program specified as external variable
#
proc RunECalc file {
    global ecalc
    exec $ecalc $file &
}

##########################################################################
#   ConfigError filename
#   --------------------
#   This routine displays a message when the configuration file has 
#   been corrupted and then causes the program to die.
#
#   13.09.94 Original   By: ACRM
#   14.09.94 Now calls exit
#
proc ConfigError filename {
    puts "Error! Missing lines in configuration file $filename."
    puts "Delete $filename and start with the default values."
    destroy .
    exit
}

##########################################################################
#   proc LoadConfig file
#   --------------------
#   Load the configuration file into the global variables
#
#   13.09.94 Original   By: ACRM
#   14.09.94 Added disulphide code
#   15.09.94 Added regrid
#   16.09.94 Added zones
#   22.09.94 Added ignores
#   24.10.94 Added potresidue, dispresidue, residuescale
#   23.05.95 Added Relax, shakeTol, vdwTol
#
proc LoadConfig filename {
    global conffile pdbfile outfile
    global potbond  potangle  pottor  potimp  pothbond  potvdwa  potvdwr \
           potelect potresidue
    global dispbond dispangle disptor dispimp disphbond dispvdwa dispvdwr\
           dispelect dispresidue
    global bondscale anglescale torscale impscale hbondscale vdwascale \
           vdwrscale electscale residuescale
    global nkeep
    global paramfile rtopfile Cutgrid Cutnonbond
    global ShowTimings ShowParams ShowDebug ShowRTop
    global ConstDielectric eta 
    global Cutdiston Cutangon Cutdistoff Cutangoff
    global AutoSS maxss SSfrom SSto Regrid
    global maxzones Zonefrom Zoneto
    global maxignores Ignorefrom Ignoreto Ignoresc
    global Relax shakeTol vdwTol

    # If the configuration file doesn't exist, simply return
    if {[file exists $filename] == 0} return

    # Open the confirguration file
    set file [open $filename r]
    
    if {[gets $file line] < 0} {ConfigError $filename}
    set conffile $line

    if {[gets $file line] < 0} {ConfigError $filename}
    set pdbfile $line

    if {[gets $file line] < 0} {ConfigError $filename}
    set outfile $line

    if {[gets $file line] < 0} {ConfigError $filename}
    set potbond $line

    if {[gets $file line] < 0} {ConfigError $filename}
    set potangle $line

    if {[gets $file line] < 0} {ConfigError $filename}
    set pottor $line

    if {[gets $file line] < 0} {ConfigError $filename}
    set potimp $line

    if {[gets $file line] < 0} {ConfigError $filename}
    set pothbond $line

    if {[gets $file line] < 0} {ConfigError $filename}
    set potvdwa $line

    if {[gets $file line] < 0} {ConfigError $filename}
    set potvdwr $line

    if {[gets $file line] < 0} {ConfigError $filename}
    set potelect $line

    if {[gets $file line] < 0} {ConfigError $filename}
    set potresidue $line

    if {[gets $file line] < 0} {ConfigError $filename}
    set dispbond $line

    if {[gets $file line] < 0} {ConfigError $filename}
    set dispangle $line

    if {[gets $file line] < 0} {ConfigError $filename}
    set disptor $line

    if {[gets $file line] < 0} {ConfigError $filename}
    set dispimp $line

    if {[gets $file line] < 0} {ConfigError $filename}
    set disphbond $line

    if {[gets $file line] < 0} {ConfigError $filename}
    set dispvdwa $line

    if {[gets $file line] < 0} {ConfigError $filename}
    set dispvdwr $line

    if {[gets $file line] < 0} {ConfigError $filename}
    set dispelect $line

    if {[gets $file line] < 0} {ConfigError $filename}
    set dispresidue $line

    if {[gets $file line] < 0} {ConfigError $filename}
    set bondscale $line

    if {[gets $file line] < 0} {ConfigError $filename}
    set anglescale $line

    if {[gets $file line] < 0} {ConfigError $filename}
    set torscale $line

    if {[gets $file line] < 0} {ConfigError $filename}
    set impscale $line

    if {[gets $file line] < 0} {ConfigError $filename}
    set hbondscale $line

    if {[gets $file line] < 0} {ConfigError $filename}
    set vdwascale $line

    if {[gets $file line] < 0} {ConfigError $filename}
    set vdwrscale $line

    if {[gets $file line] < 0} {ConfigError $filename}
    set electscale $line

    if {[gets $file line] < 0} {ConfigError $filename}
    set residuescale $line

    if {[gets $file line] < 0} {ConfigError $filename}
    set nkeep $line

    if {[gets $file line] < 0} {ConfigError $filename}
    set paramfile $line

    if {[gets $file line] < 0} {ConfigError $filename}
    set rtopfile $line

    if {[gets $file line] < 0} {ConfigError $filename}
    set Cutgrid $line

    if {[gets $file line] < 0} {ConfigError $filename}
    set Regrid $line

    if {[gets $file line] < 0} {ConfigError $filename}
    set Cutnonbond $line

    if {[gets $file line] < 0} {ConfigError $filename}
    set ShowTimings $line

    if {[gets $file line] < 0} {ConfigError $filename}
    set ShowParams $line

    if {[gets $file line] < 0} {ConfigError $filename}
    set ShowDebug $line

    if {[gets $file line] < 0} {ConfigError $filename}
    set ShowRTop $line

    if {[gets $file line] < 0} {ConfigError $filename}
    set ConstDielectric $line

    if {[gets $file line] < 0} {ConfigError $filename}
    set eta $line

    if {[gets $file line] < 0} {ConfigError $filename}
    set Cutdiston $line

    if {[gets $file line] < 0} {ConfigError $filename}
    set Cutangon $line

    if {[gets $file line] < 0} {ConfigError $filename}
    set Cutdistoff $line

    if {[gets $file line] < 0} {ConfigError $filename}
    set Cutangoff $line

    if {[gets $file line] < 0} {ConfigError $filename}
    set Relax $line

    if {[gets $file line] < 0} {ConfigError $filename}
    set shakeTol $line

    if {[gets $file line] < 0} {ConfigError $filename}
    set vdwTol $line

    # Do disulphides

    if {[gets $file line] < 0} {ConfigError $filename}
    set AutoSS $line

    # Read variables
    for {set i 1} {$i <= $maxss} {incr i} {
        if {[gets $file line] < 0} {ConfigError $filename}
        set SSfrom($i)  $line

        if {[gets $file line] < 0} {ConfigError $filename}
        set SSto($i)  $line
    }

    # Do zones
    for {set i 1} {$i <= $maxzones} {incr i} {
        if {[gets $file line] < 0} {ConfigError $filename}
        set Zonefrom($i)  $line

        if {[gets $file line] < 0} {ConfigError $filename}
        set Zoneto($i)  $line
    }

    # Do ignores
    for {set i 1} {$i <= $maxignores} {incr i} {
        if {[gets $file line] < 0} {ConfigError $filename}
        set Ignorefrom($i)  $line

        if {[gets $file line] < 0} {ConfigError $filename}
        set Ignoreto($i)  $line

        if {[gets $file line] < 0} {ConfigError $filename}
        set Ignoresc($i)  $line
    }
}

##########################################################################
#   proc SaveConfig filename
#   ------------------------
#   Write a configuration file from global variables
#
#   13.09.94 Original   By: ACRM
#   14.09.94 Added disulphide code
#   15.09.94 Added regrid
#   16.09.94 Added zones
#   22.09.94 Added ignores
#   24.10.94 Added potresidue, dispresidue, residuescale
#   23.05.95 Added Relax, shakeTol and vdwTol
#
proc SaveConfig filename {
    global conffile pdbfile outfile
    global potbond  potangle  pottor  potimp  pothbond  potvdwa  potvdwr \
           potelect potresidue
    global dispbond dispangle disptor dispimp disphbond dispvdwa dispvdwr\
           dispelect dispresidue
    global bondscale anglescale torscale impscale hbondscale vdwascale \
           vdwrscale electscale residuescale
    global nkeep
    global paramfile rtopfile Cutgrid Cutnonbond
    global ShowTimings ShowParams ShowDebug ShowRTop
    global ConstDielectric eta 
    global Cutdiston Cutangon Cutdistoff Cutangoff
    global AutoSS maxss SSfrom SSto Regrid
    global maxzones Zonefrom Zoneto
    global maxignores Ignorefrom Ignoreto Ignoresc
    global Relax shakeTol vdwTol

    # Open the confirguration file
    set file [open $filename w]
    
    puts $file $conffile
    puts $file $pdbfile
    puts $file $outfile
    puts $file $potbond
    puts $file $potangle
    puts $file $pottor
    puts $file $potimp
    puts $file $pothbond
    puts $file $potvdwa
    puts $file $potvdwr
    puts $file $potelect
    puts $file $potresidue
    puts $file $dispbond
    puts $file $dispangle
    puts $file $disptor
    puts $file $dispimp
    puts $file $disphbond
    puts $file $dispvdwa
    puts $file $dispvdwr
    puts $file $dispelect
    puts $file $dispresidue
    puts $file $bondscale
    puts $file $anglescale
    puts $file $torscale
    puts $file $impscale
    puts $file $hbondscale
    puts $file $vdwascale
    puts $file $vdwrscale
    puts $file $electscale
    puts $file $residuescale
    puts $file $nkeep
    puts $file $paramfile
    puts $file $rtopfile
    puts $file $Cutgrid
    puts $file $Regrid
    puts $file $Cutnonbond
    puts $file $ShowTimings
    puts $file $ShowParams
    puts $file $ShowDebug
    puts $file $ShowRTop
    puts $file $ConstDielectric
    puts $file $eta
    puts $file $Cutdiston
    puts $file $Cutangon
    puts $file $Cutdistoff
    puts $file $Cutangoff
    puts $file $Relax
    puts $file $shakeTol
    puts $file $vdwTol

    # Do disulphides
    puts $file $AutoSS

    # Write variables
    for {set i 1} {$i <= $maxss} {incr i} {
        global SSfrom($i) SSto($i)

        puts $file $SSfrom($i)
        puts $file $SSto($i)
    }

    # Do zones
    for {set i 1} {$i <= $maxzones} {incr i} {
        global Zonefrom($i) Zoneto($i)

        puts $file $Zonefrom($i)
        puts $file $Zoneto($i)
    }

    # Do ignores
    for {set i 1} {$i <= $maxzones} {incr i} {
        global Ignorefrom($i) Ignoreto($i)

        puts $file $Ignorefrom($i)
        puts $file $Ignoreto($i)
        puts $file $Ignoresc($i)
    }

    close $file
}


##########################################################################
#                                                                        #
#                       START OF MAIN PROGRAM                            #
#                                                                        #
##########################################################################
# Read local.tcl which defines stuff for this local setup
#
source "$env(ECALCDATA)/local.tcl"

# Initialise various stuff and load configuration file if present
# ---------------------------------------------------------------
InitWM
InitVariables
ClearDisulphides
ClearZones
ClearIgnores
LoadConfig ecalc.config

##########################################################################
# Create 2 frames and attach to the main window
# ---------------------------------------------
frame .top -relief raised -border 1
frame .bot -relief raised -border 1
pack append . .top {top fill expand padx 20} .bot {top fill expand}

##########################################################################
# To the top window, we add the text entry gadgets
# ------------------------------------------------
# Create a text gadget 'entry' for PDB filename
frame .top.pdb -bd 1m
entry .top.pdb.entry -relief sunken -width 40 -textvariable pdbfile
bind  .top.pdb.entry <Return> "focus .top.conf.entry"
label .top.pdb.label -text "PDB File:"
# Attach the entry and label to the frame
pack append .top.pdb .top.pdb.entry right .top.pdb.label left
# Attach to the main window
pack append .top .top.pdb {top fillx}

# Create a text gadget 'entry' for CONF filename
frame .top.conf -bd 1m
entry .top.conf.entry -relief sunken -width 40 -textvariable conffile
bind  .top.conf.entry <Return> "focus .top.out.entry"
label .top.conf.label -text "Conformation File:"
# Attach the entry and label to the frame
pack append .top.conf .top.conf.entry right .top.conf.label left
# Attach to the main window
pack append .top .top.conf {top fillx}

# Create a text gadget 'entry' for OUTPUT filename
frame .top.out -bd 1m
entry .top.out.entry -relief sunken -width 40 -textvariable outfile
bind  .top.out.entry <Return> "focus .bot.r.bot.run"
label .top.out.label -text "Output File:"
# Attach the entry and label to the frame
pack append .top.out .top.out.entry right .top.out.label left
# Attach to the main window
pack append .top .top.out {top fillx}


##########################################################################
# Create frames for the bottom left buttons and bottom right gadgets
# and add to the bottom window
frame .bot.l
frame .bot.r
pack append .bot .bot.l {left fillx} .bot.r {right}

##########################################################################
# Create a frame to house the quit and run buttons and another empty frame
# Add the buttons
# ---------------
frame .bot.r.bot -relief sunken -border 4
frame .bot.r.top 
# -height 10
pack append .bot.r \
        .bot.r.bot {bottom padx 20 pady 20} \
        .bot.r.top {top expand fill}

# Create options, run and quit buttons and append to frame
button .bot.r.bot.save -text "Save Config" \
                       -command "SaveConfig ecalc.config"\
                       -width 13 -height 2 -background $AdvColour \
                       -activebackground $ActiveColour
button .bot.r.bot.ss   -text "Disulphides" -command Disulphides \
                       -width 13 -height 2 -background $AdvColour \
                       -activebackground $ActiveColour
button .bot.r.bot.adv  -text "Options" -command AdvancedOptions \
                       -width 13 -height 2 -background $AdvColour \
                       -activebackground $ActiveColour
button .bot.r.bot.zone -text "Zones" -command Zones \
                       -width 13 -height 2 -background $AdvColour \
                       -activebackground $ActiveColour
button .bot.r.bot.ignores -text "Ignores" -command Ignores \
                       -width 13 -height 2 -background $AdvColour \
                       -activebackground $ActiveColour
button .bot.r.bot.run  -text "Run"  -command "quit run"  \
                       -width 13 -height 2 -background $PassiveColour \
                       -activebackground $ActiveColour
button .bot.r.bot.quit -text "Quit" -command "quit quit" \
                       -width 13 -height 2 -background $PassiveColour \
                       -activebackground $ActiveColour
pack append .bot.r.bot .bot.r.bot.save    {top padx 10 pady 10}
pack append .bot.r.bot .bot.r.bot.adv     {top padx 10 pady 10}
pack append .bot.r.bot .bot.r.bot.ss      {top padx 10 pady 10}
pack append .bot.r.bot .bot.r.bot.ignores {top padx 10 pady 10}
pack append .bot.r.bot .bot.r.bot.zone    {top padx 10 pady 10}
pack append .bot.r.bot .bot.r.bot.quit    {top padx 10 pady 10}
pack append .bot.r.bot .bot.r.bot.run     {top padx 10 pady 10}

# Make a return on the run button the same as LMB
bind .bot.r.bot.run  <Return> "quit run"

# Set default colour of Run button to cyan, but make it go red when
# it gets the focus
bind .bot.r.bot.run <FocusIn>  \
        ".bot.r.bot.run configure -background $ActiveColour"
bind .bot.r.bot.run <FocusOut> \
        ".bot.r.bot.run configure -background $PassiveColour"

##########################################################################
# Build frames for the potential and display buttons
frame .bot.l.pot
frame .bot.l.disp
frame .bot.l.labels
frame .bot.l.scale

##########################################################################
# Build a set of labels for the potential etc
label .bot.l.labels.title
label .bot.l.labels.bonds     -text "Bonds:" 
label .bot.l.labels.angles    -text "Angles:"
label .bot.l.labels.torsions  -text "Torsions:"
label .bot.l.labels.impropers -text "Impropers:"
label .bot.l.labels.hbonds    -text "HBonds:"
label .bot.l.labels.vdwa      -text "VDW Attraction:"
label .bot.l.labels.vdwr      -text "VDW Repulsion:"
label .bot.l.labels.elect     -text "Electrostatic:"
label .bot.l.labels.resid     -text "Residue (empirical):"
label .bot.l.labels.nkeep     -text "Number to keep:"
# Add the labels to the frame
pack append .bot.l.labels \
        .bot.l.labels.title     {top pady 4 fillx} \
        .bot.l.labels.bonds     {top pady 4 fillx} \
        .bot.l.labels.angles    {top pady 4 fillx} \
        .bot.l.labels.torsions  {top pady 4 fillx} \
        .bot.l.labels.impropers {top pady 4 fillx} \
        .bot.l.labels.hbonds    {top pady 4 fillx} \
        .bot.l.labels.vdwa      {top pady 4 fillx} \
        .bot.l.labels.vdwr      {top pady 4 fillx} \
        .bot.l.labels.elect     {top pady 4 fillx} \
        .bot.l.labels.resid     {top pady 4 fillx} \
        .bot.l.labels.nkeep     {top pady 12 fillx} 


##########################################################################
# Build a set of check boxes for the potential
label       .bot.l.pot.title     -text "Potential"
checkbutton .bot.l.pot.bonds     -variable potbond    -relief flat
checkbutton .bot.l.pot.angles    -variable potangle   -relief flat
checkbutton .bot.l.pot.torsions  -variable pottor     -relief flat
checkbutton .bot.l.pot.impropers -variable potimp     -relief flat
checkbutton .bot.l.pot.hbonds    -variable pothbond   -relief flat
checkbutton .bot.l.pot.vdwa      -variable potvdwa    -relief flat
checkbutton .bot.l.pot.vdwr      -variable potvdwr    -relief flat
checkbutton .bot.l.pot.elect     -variable potelect   -relief flat
checkbutton .bot.l.pot.resid     -variable potresidue -relief flat
entry       .bot.l.pot.nkeep     -relief sunken       -width 10 \
                                 -textvariable nkeep
bind        .bot.l.pot.nkeep <Return> "focus .bot.r.bot.run"

# Add these to the frame
pack append .bot.l.pot \
        .bot.l.pot.title     {top pady 4} \
        .bot.l.pot.bonds     {top pady 4} \
        .bot.l.pot.angles    {top pady 4} \
        .bot.l.pot.torsions  {top pady 4} \
        .bot.l.pot.impropers {top pady 4} \
        .bot.l.pot.hbonds    {top pady 4} \
        .bot.l.pot.vdwa      {top pady 4} \
        .bot.l.pot.vdwr      {top pady 4} \
        .bot.l.pot.elect     {top pady 4} \
        .bot.l.pot.resid     {top pady 4} \
        .bot.l.pot.nkeep     {top pady 12} 

##########################################################################
# Build a set of check boxes for display
label       .bot.l.disp.title     -text "Display"
checkbutton .bot.l.disp.bonds     -variable dispbond    -relief flat
checkbutton .bot.l.disp.angles    -variable dispangle   -relief flat
checkbutton .bot.l.disp.torsions  -variable disptor     -relief flat
checkbutton .bot.l.disp.impropers -variable dispimp     -relief flat
checkbutton .bot.l.disp.hbonds    -variable disphbond   -relief flat
checkbutton .bot.l.disp.vdwa      -variable dispvdwa    -relief flat
checkbutton .bot.l.disp.vdwr      -variable dispvdwr    -relief flat
checkbutton .bot.l.disp.elect     -variable dispelect   -relief flat
checkbutton .bot.l.disp.resid     -variable dispresidue -relief flat
label       .bot.l.disp.nkeep 
# Add these to the frame
pack append .bot.l.disp \
        .bot.l.disp.title     {top pady 4} \
        .bot.l.disp.bonds     {top pady 4} \
        .bot.l.disp.angles    {top pady 4} \
        .bot.l.disp.torsions  {top pady 4} \
        .bot.l.disp.impropers {top pady 4} \
        .bot.l.disp.hbonds    {top pady 4} \
        .bot.l.disp.vdwa      {top pady 4} \
        .bot.l.disp.vdwr      {top pady 4} \
        .bot.l.disp.elect     {top pady 4} \
        .bot.l.disp.resid     {top pady 4} \
        .bot.l.disp.nkeep     {top pady 12} 


##########################################################################
# Build a set of text entry boxes for scale
label .bot.l.scale.title     -text "Scale"
entry .bot.l.scale.bonds     -relief sunken -width 10 \
                             -textvariable bondscale
entry .bot.l.scale.angles    -relief sunken -width 10 \
                             -textvariable anglescale
entry .bot.l.scale.torsions  -relief sunken -width 10 \
                             -textvariable torscale
entry .bot.l.scale.impropers -relief sunken -width 10 \
                             -textvariable impscale
entry .bot.l.scale.hbonds    -relief sunken -width 10 \
                             -textvariable hbondscale
entry .bot.l.scale.vdwa      -relief sunken -width 10 \
                             -textvariable vdwascale
entry .bot.l.scale.vdwr      -relief sunken -width 10 \
                             -textvariable vdwrscale
entry .bot.l.scale.elect     -relief sunken -width 10 \
                             -textvariable electscale
entry .bot.l.scale.resid     -relief sunken -width 10 \
                             -textvariable residuescale
button .bot.l.scale.relax \
       -text "Relax" \
       -command "relax .bot.l.scale.relax $ActiveColour" \
       -activebackground $ActiveColour \
       -width 10 -height 1 -relief raised

# Bind <return> in each of these to step to the next one
bind  .bot.l.scale.bonds     <Return> "focus .bot.l.scale.angles"
bind  .bot.l.scale.angles    <Return> "focus .bot.l.scale.torsions"
bind  .bot.l.scale.torsions  <Return> "focus .bot.l.scale.impropers"
bind  .bot.l.scale.impropers <Return> "focus .bot.l.scale.hbonds"
bind  .bot.l.scale.hbonds    <Return> "focus .bot.l.scale.vdwa"
bind  .bot.l.scale.vdwa      <Return> "focus .bot.l.scale.vdwr"
bind  .bot.l.scale.vdwr      <Return> "focus .bot.l.scale.elect"
bind  .bot.l.scale.elect     <Return> "focus .bot.l.scale.resid"
bind  .bot.l.scale.resid     <Return> "focus .bot.l.scale.bonds"


# Add these to the frame
pack append .bot.l.scale \
        .bot.l.scale.title     {top pady 6} \
        .bot.l.scale.bonds     {top pady 6} \
        .bot.l.scale.angles    {top pady 6} \
        .bot.l.scale.torsions  {top pady 6} \
        .bot.l.scale.impropers {top pady 6} \
        .bot.l.scale.hbonds    {top pady 6} \
        .bot.l.scale.vdwa      {top pady 6} \
        .bot.l.scale.vdwr      {top pady 6} \
        .bot.l.scale.elect     {top pady 6} \
        .bot.l.scale.resid     {top pady 6} \
        .bot.l.scale.relax     {top pady 18} 

##########################################################################
# Add the label, potential, scale and display frames to the main 
# bottom window
pack append .bot.l .bot.l.labels {left}         \
                   .bot.l.pot    {left padx 30} \
                   .bot.l.scale  {left padx 30} \
                   .bot.l.disp   {left padx 30}

##########################################################################
# Finally set the focus to the PDB filename
# -----------------------------------------
focus .top.pdb.entry

##########################################################################
#                                                                        #
#                         END OF MAIN PROGRAM                            #
#                                                                        #
##########################################################################
