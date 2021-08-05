#*************************************************************************
#
#   Program:    Ecalc interface
#   File:       Disulphides.tcl
#   
#   Version:    V1.0
#   Date:       29.09.94
#   Function:   Handle disulphides for ECalc control file
#   
#   Copyright:  (c) Dr. Andrew C. R. Martin 1994
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
#   V0.1  14.09.94 Original    By: ACRM
#   V1.0  29.09.94 First release version
#   V1.1           Skipped
#   V1.2           Skipped
#   V1.3           Skipped
#   V1.4  23.05.95 Skipped
#
#*************************************************************************

##########################################################################
#   proc ClearDisulphides
#   ---------------------
#   Clear all disulphide specifications
#
#   14.09.94 Original   By: ACRM
#
proc ClearDisulphides { } {
    # Reference the global maxss variable
    global maxss SSfrom SSto

    # Create a list of variable subnames
    set disulphidespecs {from to}

    # Set all variables to blank string
    for {set i 1} {$i <= $maxss} {incr i} {
        foreach disul $disulphidespecs {
            set SS${disul}($i)  ""
        }
    }
}
    
##########################################################################
#   proc WriteDisulphides
#   ---------------------
#   Write all disulphide specifications
#
#   14.09.94 Original   By: ACRM
#
proc WriteDisulphides {file} {
    # Reference the global variables
    global maxss
    global SSfrom SSto

    for {set i 1} {$i <= $maxss} {incr i} {
        if { ($SSfrom($i) != "") && ($SSto($i) != "") } {
            puts $file "DISULPHIDE $SSfrom($i) $SSto($i)"
        }
    }
}
    
##########################################################################
#   proc Disulphides [window]
#   -------------------------
#   Create and handle the disulphide entry window
#
#   14.09.94 Original   By: ACRM
#
proc Disulphides {{w .ss}} {
    toplevel $w

    ######################################################################
    # Set up the window name bar
    # --------------------------
    wm title $w "ECalc Disulphides (c) 1994, Dr. Andrew C.R. Martin, UCL"
    wm iconname $w "ECalc Disulphides"

    ######################################################################
    # Reference global variables
    # --------------------------
    global maxss AutoSS

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
    # Create text entries for the disulphides
    # ---------------------------------------
    set disulphidespecs {from to}
    for {set i 1} {$i <= $maxss} {incr i} {
        # Build a frame for the 2 parts of the SS spec
        frame $w.specs.$i

        # For each part of the SS spec, build a text entry
        foreach disul $disulphidespecs {
            global SS${disul}($i)

            frame $w.specs.$i.$disul -bd 1m
            entry $w.specs.$i.$disul.entry -relief sunken -width 17 \
                                        -textvariable SS${disul}($i)
            # Add appropriate label
            if {$disul == "from"} {
                label $w.specs.$i.$disul.label -text "From: "
            } else {
                label $w.specs.$i.$disul.label -text "To: "
            }

            # Join entry and label onto frame
            pack append $w.specs.$i.$disul $w.specs.$i.$disul.entry right \
                        $w.specs.$i.$disul.label left

            # Join frame onto parent window
            pack append $w.specs.$i $w.specs.$i.$disul {left filly}
        }

        # Pack the SS spec pair into the parent window
        pack append $w.specs $w.specs.$i {top}
    }

    ######################################################################
    # Set up focus bindings for the text enries
    # -----------------------------------------
    for {set i 1} {$i < $maxss} {incr i} {
        set j [expr $i+1]
        bind $w.specs.$i.from.entry <Return> "focus $w.specs.$i.to.entry"
        bind $w.specs.$i.to.entry <Return> "focus $w.specs.$j.from.entry"
    }    
    bind $w.specs.$maxss.from.entry <Return> \
        "focus $w.specs.$maxss.to.entry"
    bind $w.specs.$maxss.to.entry <Return> \
        "focus $w.specs.1.from.entry"

    ######################################################################
    # Create radio buttons and button for top part
    # --------------------------------------------

    radiobutton $w.flags.auto -text "Automatic disulphide search" \
            -relief flat -anchor w \
            -variable AutoSS -value 1

    radiobutton $w.flags.noss -text "Ignore disulphides" \
            -relief flat -anchor w \
            -variable AutoSS -value 0

    radiobutton $w.flags.user -text "Manual disulphide entry" \
            -relief flat -anchor w \
            -variable AutoSS -value 2

    button $w.flags.clear -text "Clear disulphides" \
            -relief raised -command ClearDisulphides

    pack append $w.flags $w.flags.auto  {top fillx} \
                         $w.flags.noss  {top fillx} \
                         $w.flags.user  {top fillx} \
                         $w.flags.clear {top fillx} 

    ######################################################################
    # Focus the text into the first text entry
    # ----------------------------------------

    focus $w.specs.1.from.entry
}

