#*************************************************************************
#
#   Program:    Ecalc interface
#   File:       AdvancedOptions.tcl
#   
#   Version:    V1.4
#   Date:       23.05.95
#   Function:   Handle advanced options for ECalc control file
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
#   V0.1  12.09.94 Original    By: ACRM
#   V0.2  15.09.94 Added regrid option
#   V1.0  29.09.94 First release version
#   V1.1           Skipped
#   V1.2           Skipped
#   V1.3           Skipped
#   V1.4  23.05.95 Added Relax options
#
#*************************************************************************

##########################################################################
proc AdvancedOptions {{w .adv}} {
    toplevel $w

    ######################################################################
    # Set up the window name bar
    # --------------------------
    wm title $w "ECalc Options (c) 1994, Dr. Andrew C.R. Martin, UCL"
    wm iconname $w "ECalc Options"

    ######################################################################
    # Reference global variables
    # --------------------------
    global  paramfile rtopfile Cutgrid Cutnonbond
    global  ShowTimings ShowParams ShowDebug ShowRTop
    global  ConstDielectric eta 
    global  Cutdiston Cutangon Cutdistoff Cutangoff
    global  Regrid
    global  shakeTol vdwTol

    ######################################################################
    # Create 7 frames for the sections
    # --------------------------------
    set frames {cutoff hbond elect relax flags files exit}
    foreach frame $frames {
        frame $w.$frame -relief raised -border 1
        pack append $w $w.$frame {top fill expand}
    }

    ######################################################################
    # Create a button for exit
    # ------------------------
    button $w.exit.button -relief raised -text "Exit" \
            -width 6 -height 1 -command "destroy $w"
    pack append $w.exit $w.exit.button {right padx 20 pady 20}

    ######################################################################
    # Create two text entries for the files
    # -------------------------------------
    set files {param rtop}
    foreach file $files {
        frame $w.files.$file -bd 1m
        entry $w.files.$file.entry -relief sunken -width 40 \
                                    -textvariable ${file}file
        label $w.files.$file.label
        # Join entry and label onto frame
        pack append $w.files.$file $w.files.$file.entry right \
                    $w.files.$file.label left
        # Join frame onto main window
        pack append $w.files $w.files.$file {top fillx}
    }
        
    $w.files.param.label configure -text "Parameter file:"
    $w.files.rtop.label  configure -text "Residue topology file: "
    bind $w.files.param.entry <Return> "focus $w.files.rtop.entry"
    bind $w.files.rtop.entry  <Return> "focus $w.files.param.entry"

    ######################################################################
    # Create 4 buttons for the flags
    # ------------------------------
    # Build 2 sub-frames
    frame $w.flags.l
    frame $w.flags.r

    checkbutton $w.flags.l.timings -variable ShowTimings -relief flat \
            -text "Timings" -anchor w
    checkbutton $w.flags.l.params  -variable ShowParams  -relief flat \
            -text "Show parameters" -anchor w
    checkbutton $w.flags.r.debug   -variable ShowDebug   -relief flat \
            -text "Debug" -anchor w
    checkbutton $w.flags.r.rtop    -variable ShowRTop    -relief flat \
            -text "Show residue topology" -anchor w

    pack append $w.flags.l $w.flags.l.timings {top expand fillx} \
            $w.flags.l.params  {bottom expand fillx} 
    pack append $w.flags.r $w.flags.r.debug   {top expand fillx} \
            $w.flags.r.rtop    {bottom expand fillx} 

    pack append $w.flags $w.flags.l {left} $w.flags.r {right}

    ######################################################################
    # Create two text entries for the Relax tolerence values
    # ------------------------------------------------------

    set tols {shake vdw}
    foreach tol $tols {
        # Build the sub-frame, label and entry
        frame $w.relax.$tol -bd 1m
        entry $w.relax.$tol.entry -relief sunken -width 10 \
                                  -textvariable ${tol}Tol
        label $w.relax.$tol.label

        # Join entry and label onto sub-entry
        pack append $w.relax.$tol $w.relax.$tol.entry right \
                                  $w.relax.$tol.label left
    }

    # Join frames onto parent window
    pack append $w.relax $w.relax.shake {left} \
                         $w.relax.vdw   {right}

    # Set the labels
    $w.relax.shake.label configure -text "Shake Tolerence:           "
    $w.relax.vdw.label   configure -text "VDW Tolerence:              "

    ######################################################################
    # Create radio buttons and text entry for electrostatics
    # ------------------------------------------------------
    # Build 2 sub-frames
    frame $w.elect.l
    frame $w.elect.r

    radiobutton $w.elect.l.const -text "Constant dielectric" \
            -relief flat -anchor w\
            -variable ConstDielectric -value 1

    radiobutton $w.elect.l.dist -text "Distance dependent dielectric" \
            -relief flat -anchor w\
            -variable ConstDielectric -value 0

    frame $w.elect.r.eta -bd 1m
    entry $w.elect.r.eta.entry -relief sunken -width 10 \
            -textvariable eta
    label $w.elect.r.eta.label -text "Eta: "
    # Join entry and label onto frame
    pack append $w.elect.r.eta $w.elect.r.eta.entry right \
            $w.elect.r.eta.label left
    label $w.elect.r.dummy
    # Join frame onto main window
    pack append $w.elect.r $w.elect.r.eta {top fillx} \
            $w.elect.r.dummy {top fillx}

    pack append $w.elect.l $w.elect.l.const {top fillx} \
            $w.elect.l.dist {bottom fillx}
    pack append $w.elect $w.elect.l {left} $w.elect.r {right}

    bind $w.elect.r.eta.entry <Return> "set ConstDielectric 1"

    ######################################################################
    # Create text entries for the HBond cutoffs
    # -----------------------------------------
    # Build 2 sub-frames
    frame $w.hbond.l
    frame $w.hbond.r

    set hbonds {diston angon}
    foreach hbond $hbonds {
        frame $w.hbond.l.$hbond -bd 1m
        entry $w.hbond.l.$hbond.entry -relief sunken -width 10 \
                                    -textvariable Cut$hbond
        label $w.hbond.l.$hbond.label
        # Join entry and label onto frame
        pack append $w.hbond.l.$hbond $w.hbond.l.$hbond.entry right \
                    $w.hbond.l.$hbond.label left
        # Join frame onto parent window
        pack append $w.hbond.l $w.hbond.l.$hbond {top fillx}
    }
        
    $w.hbond.l.diston.label configure -text "HBond Distance Cut On: "
    $w.hbond.l.angon.label  configure -text "HBond Angle Cut On:"

    set hbonds {distoff angoff}
    foreach hbond $hbonds {
        frame $w.hbond.r.$hbond -bd 1m
        entry $w.hbond.r.$hbond.entry -relief sunken -width 10 \
                                    -textvariable Cut$hbond
        label $w.hbond.r.$hbond.label
        # Join entry and label onto frame
        pack append $w.hbond.r.$hbond $w.hbond.r.$hbond.entry right \
                    $w.hbond.r.$hbond.label left
        # Join frame onto parent window
        pack append $w.hbond.r $w.hbond.r.$hbond {top fillx}
    }
        
    $w.hbond.r.distoff.label configure -text "HBond Distance Cut Off: "
    $w.hbond.r.angoff.label  configure -text "HBond Angle Cut Off:"

    # Join the 2 sub-windows onto the parent
    pack append $w.hbond $w.hbond.l {left} $w.hbond.r {right}

    # Get return key to modify focus
    bind $w.hbond.l.diston.entry <Return> "focus $w.hbond.r.distoff.entry"
    bind $w.hbond.r.distoff.entry <Return> "focus $w.hbond.l.angon.entry"
    bind $w.hbond.l.angon.entry <Return> "focus $w.hbond.r.angoff.entry"
    bind $w.hbond.r.angoff.entry <Return> "focus $w.hbond.l.diston.entry"

    ######################################################################
    # Create text entries for the Grid and non-bonded cutoffs
    # -------------------------------------------------------
    set cutoffs {grid nonbond regrid}
    foreach cutoff $cutoffs {
        frame $w.cutoff.$cutoff -bd 1m
        entry $w.cutoff.$cutoff.entry -relief sunken -width 10 \
                                    -textvariable Cut$cutoff
        label $w.cutoff.$cutoff.label
        # Join entry and label onto frame
        pack append $w.cutoff.$cutoff $w.cutoff.$cutoff.entry right \
                    $w.cutoff.$cutoff.label left
    }
        
    $w.cutoff.grid.label    configure -text "Grid Cutoff: "
    $w.cutoff.nonbond.label configure -text "Non-bond Cutoff: "
    $w.cutoff.regrid.label  configure -text "Regrid: "
    $w.cutoff.regrid.entry  configure -textvariable Regrid

    pack append $w.cutoff $w.cutoff.grid {left}   \
                          $w.cutoff.regrid {left} \
                          $w.cutoff.nonbond {right}

    # Get return to change focus
    bind $w.cutoff.grid.entry    <Return> "focus $w.cutoff.regrid.entry"
    bind $w.cutoff.regrid.entry  <Return> "focus $w.cutoff.nonbond.entry"
    bind $w.cutoff.nonbond.entry <Return> "focus $w.cutoff.grid.entry"
}

