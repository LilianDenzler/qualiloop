/***********************************************************************
 *      Name: LayoutProto.h                                            *
 *  Function: Screen layout prototypes                                 *
 * Copyright: (C) Oxford Molecular Limited                             *
 *---------------------------------------------------------------------*
 *    Author: D. Walker, Tessella Support Services plc                 *
 *      Date: 12/06/89                                                 *
 *---------------------------------------------------------------------*
 * Modification Record                                                 *
 * Date     Inits   Comments                                           *
 * dd/mm/yy                                                            *
 ***********************************************************************/

/* Flag to indicate inclusion of this file */
#ifndef __LAYOUTPROTO__
#define __LAYOUTPROTO__

/* Define the Oxford Molecular dialog colour scheme */

#define OMRED_RED     255
#define OMRED_GREEN   0
#define OMRED_BLUE    0

#define OMBLUE_RED    0
#define OMBLUE_GREEN  0
#define OMBLUE_BLUE   255

#define OMWHITE_RED   255
#define OMWHITE_GREEN 255
#define OMWHITE_BLUE  255

#define OMGREY_RED    128
#define OMGREY_GREEN  128
#define OMGREY_BLUE   128

#define OMLBLU_RED    0
#define OMLBLU_GREEN  160
#define OMLBLU_BLUE   255

/* Layout parameter definitions */

#define PIMMS_TITLE_HEIGHT    0.0292 /* Title box pixel height */
#define PIMMS_TITLE_FORERED   0    /* Title text foreground colour (red gun) */
#define PIMMS_TITLE_FOREGREEN 0  /* Title text foreground colour (Green gun) */
#define PIMMS_TITLE_FOREBLUE  0   /* Title text foreground colour (Blue gun) */

#define PIMMS_BACKGND_RED    0   /* Main screen background colour (red gun) */
#define PIMMS_BACKGND_GREEN  0   /* Main screen background colour (green gun) */
#define PIMMS_BACKGND_BLUE   0   /* Main screen background colour (blue gun) */

#define PIMMS_GAP_SIZE 0.001 /* Size of gap between menu & title etc.. */

#define SMANIP_WIDTH     0.26     /* Fractional width of manipulator */
#define SMANIP_HEIGHT    0.34     /* Fractional height of manipulator */

#define SMANIP_BACKRED   OMBLUE_RED   /* Manipulator background red gun */
#define SMANIP_BACKGREEN OMBLUE_GREEN /* Manipulator background green gun */
#define SMANIP_BACKBLUE  OMBLUE_BLUE  /* Manipulator background blue gun */

#define SMANIP_XMARGIN   0.0125   /* Size of X-Margins in manipulator */
#define SMANIP_YMARGIN   0.0125   /* Size of Y-Margins in manipulator */
#define SMANIP_ACTIVE    0.2      /* Fractional length of Active lablel */
#define SMANIP_NUMINC    5        /* Number of incremental steps in buttons */
#define SMANIP_CONTIN    SMANIP_NUMINC+1 /* Main button */
#define SMANIP_DIRECTION 0        /* Direction button */

#define ACTVTOR_WIDTH    0.45
#define ACTVTOR_HEIGHT   0.07
#define ACTVTOR_LWIDTH   0.15
#define ACTVTOR_XMARGIN  0.0125
#define ACTVTOR_YMARGIN  0.05

/* drm 29/10/91 update to be different on HP and SG */
#if defined (__SILICONGRAPHICS__) || defined (__IBM__) 
#define ACTVTOR_XFUDGE   1.02
#define ACTVTOR_YFUDGE   1.45
#define ACTVTOR_X        0.01
#endif
#if defined (__HEWLETTPACKARD__) 
#define ACTVTOR_XFUDGE   1.02
#define ACTVTOR_YFUDGE   0.90
#define ACTVTOR_X        -0.01
#endif

#define ACTVTOR_PAGE     8        /* Number of Active/In view per "page" */

/* drm 29/10/91 update to be different on HP and SG */
#if defined (__SILICONGRAPHICS__) || defined (__IBM__) 
#define SMANIP_XFUDGE    1.02     /* X position fudge for border */
#define SMANIP_YFUDGE    1.08     /* Y position fudge for border */
#endif
#if defined (__HEWLETTPACKARD__) 
#define SMANIP_XFUDGE    1.02     /* X position fudge for border */
#define SMANIP_YFUDGE    0.98     /* Y position fudge for border */
#endif

#define SMONITOR_WIDTH	0.95
#define SMONITOR_HEIGHT	0.2
#define SMONITOR_XFUDGE	1.02
#define SMONITOR_YFUDGE	1.08

#define SMONITOR_BACKRED	255
#define SMONITOR_BACKGREEN	255
#define SMONITOR_BACKBLUE	255

#define SMONITOR_DASHED_RED	0
#define SMONITOR_DASHED_GREEN	255
#define SMONITOR_DASHED_BLUE	0

#define SMONITOR_GEOML 0.236
#define SMONITOR_ANGLEL 0.318
#define SMONITOR_TORSIONL 0.4

#define SMONITOR_GAP ((1.0 - SMONITOR_GEOML - SMONITOR_ANGLEL - SMONITOR_TORSIONL) / 4.0)

#define SMONITOR_GEOMX SMONITOR_GAP
#define SMONITOR_ANGLEX (2.0*SMONITOR_GAP + SMONITOR_GEOML)
#define SMONITOR_TORSIONX (3.0*SMONITOR_GAP + SMONITOR_GEOML + SMONITOR_ANGLEL)

#define SMONITOR_GEOMH 0.6
#define SMONITOR_ANGLEH SMONITOR_GEOMH
#define SMONITOR_TORSIONH SMONITOR_GEOMH

#define SMONITOR_GEOMY 0.2
#define SMONITOR_ANGLEY SMONITOR_GEOMY
#define SMONITOR_TORSIONY SMONITOR_GEOMY

#define SMONITOR_ENERGYX 0.35
#define SMONITOR_ENERGYY 0.04
#define SMONITOR_ENERGYL 0.25
#define SMONITOR_ENERGYH 0.12

#define FITTING_DASHED_RED	255
#define FITTING_DASHED_GREEN	0
#define FITTING_DASHED_BLUE	0

#define MOL_COL_RED   255  /* Default molecule colour red gun value */
#define MOL_COL_GREEN 0   /* Default molecule colour gren gun value */
#define MOL_COL_BLUE  0   /* Default molecule colour blue gun value */

#define MSGBOX_XPOS   0.25  /* X position of DoMessage window */
#define MSGBOX_YPOS   0.4  /* Y Position of DoMessage window */
#define MSGBOX_LENGTH 0.5  /* Fractional length of message box */
#define MSGBOX_HEIGHT 0.2  /* Fractional Height of message box */
#define MSGBOX_RED    255  /* Red gun for message box background */
#define MSGBOX_GREEN  0   /* Green gun for message box background */
#define MSGBOX_BLUE   0   /* Blue gun for message box background */

#define MSGBOX_BORDER_WIDTH 0.01 /* Border width for Message box */

#define MSGBOX_BORDER_RED   255  /* Red gun for message box border */
#define MSGBOX_BORDER_GREEN 255  /* Green gun for message box border */
#define MSGBOX_BORDER_BLUE  255  /* Blue gun for message box border */

#define MSGBOX_TEXTX_INDENT 0.1  /* Text indentation for X */
#define MSGBOX_TEXTY_INDENT 0.55 /* Text indentation for Y */

#define MSGBOX_OK_XINDENT   0.7  /* X Indentation for OK button */
#define MSGBOX_OK_YINDENT   0.1  /* Y Indentation for OK button */
#define MSGBOX_OK_LENGTH    0.2  /* Length of OK button */
#define MSGBOX_OK_HEIGHT    0.2  /* Height of OK button */

#define MSGBOX_CANCEL_XINDENT 0.1 /* X Indentation for cancel button */

#define DIRBOX_XPOS	0.05
#define DIRBOX_YPOS	0.05
#define DIRBOX_LENGTH	0.3
#define DIRBOX_HEIGHT	0.65
#define DIRMAX_LINES    24
#define DIRBOX_RED    255
#define DIRBOX_GREEN  0
#define DIRBOX_BLUE   0

#define DIRBOX_BORDER_WIDTH 0.01 /* Border width for Message box */

#define DIRBOX_BORDER_RED   OMWHITE_RED   /* Red gun for message box border */
#define DIRBOX_BORDER_GREEN OMWHITE_GREEN /* Green gun for message box border */
#define DIRBOX_BORDER_BLUE  OMWHITE_BLUE  /* Blue gun for message box border */

#define DIRDIR_NAME_RED   OMBLUE_RED
#define DIRDIR_NAME_GREEN OMBLUE_GREEN
#define DIRDIR_NAME_BLUE  OMBLUE_BLUE

#define DIRALTDIR_NAME_RED   OMWHITE_RED
#define DIRALTDIR_NAME_GREEN OMWHITE_GREEN
#define DIRALTDIR_NAME_BLUE  OMWHITE_BLUE

#define DIRBACKGROUND_RED   OMGREY_RED
#define DIRBACKGROUND_GREEN OMGREY_GREEN
#define DIRBACKGROUND_BLUE  OMGREY_BLUE

#define DIRBUTTON_RED   OMBLUE_RED
#define DIRBUTTON_GREEN OMBLUE_GREEN
#define DIRBUTTON_BLUE  OMBLUE_BLUE

#define DIRBOX_OK_XINDENT  0.05  /* X Indentation for OK button */
#define DIRBOX_OK_YINDENT   0.05  /* Y Indentation for OK button */
#define DIRBOX_OK_LENGTH    0.2  /* Length of OK button */
#define DIRBOX_OK_HEIGHT    0.05  /* Height of OK button */

#define DIRBOX_CANCEL_XINDENT 0.8 /* X Indentation for cancel button */
#define DIRBOX_PGUP_XINDENT 0.325 /* X Indentation for cancel button */
#define DIRBOX_PGDN_XINDENT 0.525 /* X Indentation for cancel button */

#define DIRBOX_TEXTX_INDENT 0.025

#define DIRTITLE_LINES 2
#define DIRALTDIR_TITLE_LINES 1
#define DIRALTDIR_LINES 1
#define DIRBUTTON_LINES 7

#define DIRRADIO_HEIGHT 0.035
#define DIRRADIO_LENGTH 0.035

#define CURS_COL_DEFRED   OMRED_RED   /* Default cursor colour */
#define CURS_COL_DEFGREEN OMRED_GREEN /* Default cursor colour */
#define CURS_COL_DEFBLUE  OMRED_BLUE  /* Default cursor colour */

#define CURS_COL_MSGRED   255    /* Cursor colour in message box (Red Gun) */
#define CURS_COL_MSGGREEN 255    /* Cursor colour in message box (Green Gun) */
#define CURS_COL_MSGBLUE  0      /* Cursor colour in message box (Blue Gun) */

#define CURS_COL_DIRRED   OMRED_RED   /* Cursor colour in message box (Red Gun) */
#define CURS_COL_DIRGREEN OMRED_GREEN /* Cursor colour in message box (Green Gun) */
#define CURS_COL_DIRBLUE  OMRED_BLUE  /* Cursor colour in message box (Blue Gun) */

#define CURS_COL_DLOG_RED   OMRED_RED   /* Cursor colour in dialog box (Red Gun) */
#define CURS_COL_DLOG_GREEN OMRED_GREEN /* Cursor colour in dialog box (Green Gun) */
#define CURS_COL_DLOG_BLUE  OMRED_BLUE  /* Cursor colour in dialog box (Blue Gun) */

#define PIMMS_TPORT_X      0.05  /* X-indentation for text port */
#define PIMMS_TPORT_Y      0.05  /* Y-indentation for text port */

#define PIMMS_TPORT_LENGTH 0.45  /* X-Length of text port */
#define PIMMS_TPORT_HEIGHT 0.15  /* Y-Height of text port */

/* Default molecule colours */
#define DEFMOLCOL_RED   OMRED_RED
#define DEFMOLCOL_GREEN OMRED_GREEN
#define DEFMOLCOL_BLUE  OMRED_BLUE

/* default atom colours */
#define DEFATMCOL_RED   180
#define DEFATMCOL_GREEN 180
#define DEFATMCOL_BLUE  180

/* Clipping plane default values */
#define PIMMS_CLIP_XLEFT   -50.0
#define PIMMS_CLIP_XRIGHT   50.0

#define PIMMS_CLIP_YBOTTOM -50.0
#define PIMMS_CLIP_YTOP     50.0

#define PIMMS_CLIP_ZNEAR    50.0
#define PIMMS_CLIP_ZFAR    -50.0

#define PIMMS_CLIP_ZMINSEP   0.1

/* Define the user's viewpoint */

#define PIMMS_VIEWPOINTX 0.0
#define PIMMS_VIEWPOINTY 0.0
#define PIMMS_VIEWANGLE  0.0

/* define atom label positioning */
#define HIGHLIGHTED_ATOM_RADIUS 0.003
#define ATOM_LABEL_XOFFSET      0.01
#define ATOM_LABEL_YOFFSET      0.0
#define DUMMY_ATOM_SIZE         0.006

/* Define a structure type to hold clipping planes */
typedef struct {
   double XLeft;
   double XRight;
   double YBottom;
   double YTop;
   double ZNear;
   double ZFar;
} ClipType;

/* Colour slider dialog definitions */

#define COLSLID_WIDTH     0.3      /* Fractional width of colour slider */
#define COLSLID_HEIGHT    0.2      /* Fractional height of colour slider */

#define COLSLID_BACKRED   OMBLUE_RED   /* Colour Slider background red gun */
#define COLSLID_BACKGREEN OMBLUE_GREEN /* Colour Slider background green gun */
#define COLSLID_BACKBLUE  OMBLUE_BLUE  /* Colour Slider background blue gun */

#define COLSLID_MARGIN    0.0125   /* Size of X-Margins in colour slider */
#define SLIDER_HEIGHT     0.2      /* Height of the sliders */
#define SLIDER_LENGTH     0.6      /* Length of the sliders */
#define COLBOX_HEIGHT     0.6      /* Colour box fractional height */
#define COLBOX_LENGTH     0.3      /* Colour box fractional length */

#define COLSLID_BORDER    2        /* Border pixel width */

#define COLSLID_OUTRED    OMRED_RED   /* Outer border colour red gun value */
#define COLSLID_OUTGREEN  OMRED_GREEN /* Outer border colour green gun value */
#define COLSLID_OUTBLUE   OMRED_BLUE  /* Outer border colour blue gun value */

#define COLSLID_INRED     OMWHITE_RED   /* Inner border colour red gun value */
#define COLSLID_INGREEN   OMWHITE_GREEN /* Inner border colour green gun value */
#define COLSLID_INBLUE    OMWHITE_BLUE  /* Inner border colour blue gun value */

#define COLSLID_BUTHEIGHT 0.15     /* OK/Cancel button heights */

#define COLSLID_XFUDGE    1.02     /* X position fudge for border */
#define COLSLID_YFUDGE    1.08     /* Y position fudge for border */

#endif	
