#*************************************************************************
#
#   Program:    Ecalc interface
#   File:       Zones.tcl
#   
#   Version:    V1.4
#   Date:       23.05.95
#   Function:   Handle zones for ECalc control file
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
#
#   V0.1  16.09.94 Original    By: ACRM
#   V1.0  29.09.94 First release
#   V1.1           Skipped
#   V1.2           Skipped
#   V1.3           Skipped
#   V1.4  23.05.95 Skipped
#
#*************************************************************************


##########################################################################
#   proc ClearZones
#   ---------------
#   Clear all zone  specifications
#
#   16.09.94 Original   By: ACRM
#
proc ClearZones { } {
    # Reference the global maxss variable
    global maxzones Zonefrom Zoneto

    # Create a list of variable subnames
    set Zonespecs {from to}

    # Set all variables to blank string
    for {set i 1} {$i <= $maxzones} {incr i} {
        foreach zone $Zonespecs {
            set Zone${zone}($i)  ""
        }
    }
}
    
##########################################################################
#   proc WriteZones
#   ---------------
#   Write all zone specifications
#
#   16.09.94 Original   By: ACRM
#
proc WriteZones {file} {
    # Reference the global variables
    global maxzones
    global Zonefrom Zoneto

    for {set i 1} {$i <= $maxzones} {incr i} {
        if { ($Zonefrom($i) != "") && ($Zoneto($i) != "") } {
            puts $file "ZONE $Zonefrom($i) $Zoneto($i)"
        }
    }
}
    
##########################################################################
#   proc Zones [window]
#   -------------------
#   Create and handle the zone entry window
#
#   16.09.94 Original   By: ACRM
#
proc Zones {{w .zones}} {
    toplevel $w

    ######################################################################
    # Set up the window name bar
    # --------------------------
    wm title $w "ECalc Zones (c) 1994, Dr. Andrew C.R. Martin, UCL"
    wm iconname $w "ECalc Zones"

    ######################################################################
    # Reference global variables
    # --------------------------
    global maxzones

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
    # Create text entries for the Zones
    # ---------------------------------
    set zonespecs {from to}
    for {set i 1} {$i <= $maxzones} {incr i} {
        # Build a frame for the 2 parts of the zone spec
        frame $w.specs.$i

        # For each part of the zone spec, build a text entry
        foreach zone $zonespecs {
            global Zone${zone}($i)

            frame $w.specs.$i.$zone -bd 1m
            entry $w.specs.$i.$zone.entry -relief sunken -width 17 \
                                        -textvariable Zone${zone}($i)
            # Add appropriate label
            if {$zone == "from"} {
                label $w.specs.$i.$zone.label -text "From: "
            } else {
                label $w.specs.$i.$zone.label -text "To: "
            }

            # Join entry and label onto frame
            pack append $w.specs.$i.$zone $w.specs.$i.$zone.entry right \
                        $w.specs.$i.$zone.label left

            # Join frame onto parent window
            pack append $w.specs.$i $w.specs.$i.$zone {left filly}
        }

        # Pack the zone spec pair into the parent window
        pack append $w.specs $w.specs.$i {top}
    }

    ######################################################################
    # Set up focus bindings for the text enries
    # -----------------------------------------
    for {set i 1} {$i < $maxzones} {incr i} {
        set j [expr $i+1]
        bind $w.specs.$i.from.entry <Return> "focus $w.specs.$i.to.entry"
        bind $w.specs.$i.to.entry <Return> "focus $w.specs.$j.from.entry"
    }    
    bind $w.specs.$maxzones.from.entry <Return> \
        "focus $w.specs.$maxzones.to.entry"
    bind $w.specs.$maxzones.to.entry <Return> \
        "focus $w.specs.1.from.entry"

    ######################################################################
    # Create button for top part
    # --------------------------

    button $w.flags.clear -text "Clear Zones" \
            -relief raised -command ClearZones

    pack append $w.flags $w.flags.clear {top fillx} 

    ######################################################################
    # Focus the text into the first text entry
    # ----------------------------------------

    focus $w.specs.1.from.entry
}

