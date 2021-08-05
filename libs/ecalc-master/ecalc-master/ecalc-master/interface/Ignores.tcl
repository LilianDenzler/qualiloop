#*************************************************************************
#
#   Program:    Ecalc interface
#   File:       Ignores.tcl
#   
#   Version:    V1.4
#   Date:       23.05.95
#   Function:   Handle ignores for ECalc control file
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
#   EMail:      INTERNET: martin@biochem.ucl.ac.uk
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
#
#   V0.1  22.09.94 Original    By: ACRM
#   V1.0  29.09.94 First release version
#   V1.1           Skipped
#   V1.2           Skipped
#   V1.3           Skipped
#   V1.4  23.05.95 Skipped
#
#*************************************************************************

##########################################################################
#   proc ClearIgnores
#   -----------------
#   Clear all ignore  specifications
#
#   16.09.94 Original   By: ACRM
#
proc ClearIgnores { } {
    # Reference the global maxss variable
    global maxignores Ignorefrom Ignoreto Ignoresc

    # Create a list of variable subnames
    set Ignorespecs {from to}

    # Set all variables to blank string
    for {set i 1} {$i <= $maxignores} {incr i} {
        foreach ignore $Ignorespecs {
            set Ignore${ignore}($i)  ""
        }
        set Ignoresc($i) 0
    }
}
    
##########################################################################
#   proc WriteIgnores
#   -----------------
#   Write all ignore specifications
#
#   16.09.94 Original   By: ACRM
#
proc WriteIgnores {file} {
    # Reference the global variables
    global maxignores Ignorefrom Ignoreto Ignoresc

    for {set i 1} {$i <= $maxignores} {incr i} {
        if { ($Ignorefrom($i) != "") && ($Ignoreto($i) != "") } {
            if { $Ignoresc($i) == 1 } {
                puts $file "IGNORE $Ignorefrom($i) $Ignoreto($i) SIDE"
            } else {
                puts $file "IGNORE $Ignorefrom($i) $Ignoreto($i)"
            }
        }
    }
}
    
##########################################################################
#   proc Ignores [window]
#   -------------------
#   Create and handle the ignore entry window
#
#   16.09.94 Original   By: ACRM
#
proc Ignores {{w .ignores}} {
    toplevel $w

    ######################################################################
    # Set up the window name bar
    # --------------------------
    wm title $w "ECalc Ignores (c) 1994, Dr. Andrew C.R. Martin, UCL"
    wm iconname $w "ECalc Ignores"

    ######################################################################
    # Reference global variables
    # --------------------------
    global maxignores

    ######################################################################
    # Create 3 frames for the sections
    # --------------------------------
    set frames {flags specs exit}
    foreach frame $frames {
        frame $w.$frame -relief raised -border 1
        pack append $w $w.$frame {top fill expand pady 2}
    }

    ######################################################################
    # Create a button for exit
    # ------------------------
    button $w.exit.button -relief raised -text "Exit" \
            -width 6 -height 1 -command "destroy $w"
    pack append $w.exit $w.exit.button {right padx 20 pady 20}

    ######################################################################
    # Create text entries for the Ignores
    # ---------------------------------
    set ignorespecs {from to}
    for {set i 1} {$i <= $maxignores} {incr i} {
        # Build a frame for the 2 parts of the ignore spec and side button
        frame $w.specs.$i

        # For each part of the ignore spec, build a text entry
        foreach ignore $ignorespecs {
            global Ignore${ignore}($i)

            frame $w.specs.$i.$ignore -bd 1m
            entry $w.specs.$i.$ignore.entry -relief sunken -width 17 \
                                        -textvariable Ignore${ignore}($i)
            # Add appropriate label
            if {$ignore == "from"} {
                label $w.specs.$i.$ignore.label -text "From: "
            } else {
                label $w.specs.$i.$ignore.label -text "To: "
            }

            # Join entry and label onto frame
            pack append $w.specs.$i.$ignore \
                        $w.specs.$i.$ignore.entry right \
                        $w.specs.$i.$ignore.label left

            # Join frame onto parent window
            pack append $w.specs.$i $w.specs.$i.$ignore {left filly}
        }

        # Build a button for the SIDE flag
        checkbutton $w.specs.$i.side -text "Sidechains" \
                -variable Ignoresc($i) \
                -relief flat

        # Pack the button onto the specification
        pack append $w.specs.$i $w.specs.$i.side {left filly}

        # Pack the ignore spec trio into the parent window
        pack append $w.specs $w.specs.$i {top}
    }

    ######################################################################
    # Set up focus bindings for the text enries
    # -----------------------------------------
    for {set i 1} {$i < $maxignores} {incr i} {
        set j [expr $i+1]
        bind $w.specs.$i.from.entry <Return> "focus $w.specs.$i.to.entry"
        bind $w.specs.$i.to.entry <Return> "focus $w.specs.$j.from.entry"
    }    
    bind $w.specs.$maxignores.from.entry <Return> \
        "focus $w.specs.$maxignores.to.entry"
    bind $w.specs.$maxignores.to.entry <Return> \
        "focus $w.specs.1.from.entry"

    ######################################################################
    # Create button for top part
    # --------------------------

    button $w.flags.clear -text "Clear Ignores" \
            -relief raised -command ClearIgnores

    pack append $w.flags $w.flags.clear {top fillx} 

    ######################################################################
    # Focus the text into the first text entry
    # ----------------------------------------

    focus $w.specs.1.from.entry
}

