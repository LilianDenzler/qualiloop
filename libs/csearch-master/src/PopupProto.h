/***********************************************************************
 *      Name: PopupProtoTypes.h                                        *
 *  Function: Popup Menu Manager Prototypes                            *
 * Copyright: (C) Oxford Molecular Limited                             *
 *---------------------------------------------------------------------*
 *    Author: D. Walker, Tessella Support Services plc                 *
 *      Date: 12/06/89                                                 *
 *---------------------------------------------------------------------*
 * Modification Record                                                 *
 * Date     Inits   Comments                                           *
 * 11/12/92 DW      Allow for the different size of Motif menu bar     *
 * 30/03/93 JAW     Increase MAX_POPUP to 20 (CR1478)                  *
 ***********************************************************************/
#ifndef __POPUPPROTO__
#define __POPUPPROTO__

/* Flag to indicate this file has been included */
#define POPUP_DEF

/* Menu parameter definition */

/* OMUPD DW 11/12/92 Allow for the different size of the Motif menu bar */
#ifdef __MOTIF_GUI__
#define PIMMS_MENU_HEIGHT    0.080 /* Menu bar height */
#else
#define PIMMS_MENU_HEIGHT    0.039 /* Menu bar height */
#endif

#define PIMMS_MENU_NOPTS     6     /* Number of options in menu bar */

#define PIMMS_MAX_POPUP 20        /* Maximum number of items per pop up menu */

#define PIMMS_MENU_FORERED   0     /* Menu bar foreground colour (red gun) */
#define PIMMS_MENU_FOREGREEN 0     /* Menu bar foreground colour (green gun) */
#define PIMMS_MENU_FOREBLUE  0     /* Menu bar foreground colour (blue gun) */

/* Menu structure */

typedef struct MenuType {
   ButtonType Buttons[PIMMS_MAX_POPUP];	/* Menu buttons */
   short NumItems;			/* Number of items in menu */
   short XLeft;				/* Left X position of menu box */
   short XRight;			/* Right X position of menu box */
   short YTop;				/* Top Y position of menu box */
   short YBottom;			/* Bottom Y position of menu box */
} MenuType;

/* structure to hold details of a menu window contents */
typedef struct {
      short NButtons;			/* number of square buttons */
      short NRadioGroups;		/* number of radio groups */
      short NCrossGroups;		/* number of cross buttons */
      short BoxX, BoxY;			/* bottom left of window */
      short ButXPos, ButYPos;		/* coordinates of top left button */
      short ButXSpace, ButYSpace;	/* spacing of buttons */
      short ButPerLine;			/* number of buttons per line */
      short Border;			/* size of window border */
      short ButLength, ButHeight;	/* size of square buttons */
      short SquareButtonTop;		/* top of square button area */
      short LRButXPos;			/* x pos. of left col. of radio buts */
      short RRButXPos;			/* x pos. of right col. of radio buts */
      short RButYPos;			/* y pos. of first radio button */
      short RTitleYPos;			/* y pos. of first radio button title */
      short RadioButtonTop;		/* top of radio button area */
      short RButLength;			/* size of radio buttons */
      short LCButXPos;			/* x pos. of left col. of cross buts */
      short RCButXPos;			/* x pos. of right col. of cross buts */
      short CButYPos;			/* y pos. of first cross button */
      short CTitleYPos;			/* y pos. of first cross button title */
      short CrossButtonTop;		/* top of cross button area */
      short BoxHeight, BoxLength;	/* total height of menu box */
      short TitleXPos, TitleYPos;	/* position of title */
      short TitleLength, TitleHeight;	/* soze of title box */
      short CommentXPos, CommentYPos;	/* position of comment box */
      short CommentLength,CommentHeight;/* position of comment box */
      short CommentTop;			/* y position of top of comments area */
      short TextOffX, TextOffY;		/* offset of text entry box */
      short TextLength, TextHeight;	/* size of text entry box */
   } DlogMenuSize, *DlogMenuSizePtr;

/* structure to hold the details of a menu window */
typedef struct {
      long WindowID;			/* window identifier */
      ButtonTypePtr ButtonPtr;		/* pointer to button structure */
      ButtonTypePtr RadioButtonPtr;	/* pointer to radio button structure */
      ButtonTypePtr CrossButtonPtr;	/* pointer to cross button structure */
      ButtonTypePtr TextButtonPtr;	/* pointer to text button structure */
      ButtonType TitleButton;		/* title button */
      ButtonType CommentButton;		/* comment button */
      short *PixelArray;		/* pointer to pixel storage */
      short ActiveTextBox;		/* indicator of the active text box */
      long InWindowID;			/* indicator for cursor location */
      long BufferState;			/* state of the fore and back buffers */
   } DlogMenuWindow , * DlogMenuWindowPtr ;

/* New structure to hold the details of a menu window - see 'NewInitDlog' etc */
typedef struct {
      short         BoxX ;		/* bottom left of window */
      short         BoxY ;		/* bottom left of window */
      short         BoxHeight ;		/* total height, length of menu box */
      short         BoxLength ;		/* total height, length of menu box */
      short         BoxTop ;		/* top of box excluding title area */
      short         Border ;		/* size of window border */
      long          WindowID ;		/* window identifier */
      ButtonTypePtr TitlePtr ;		/* pointer to title button struct */
      ButtonTypePtr CommentPtr ;	/* pointer to comment button struct */
      ButtonTypePtr ButtonPtr ;		/* pointer to button structures */
      ButtonTypePtr GroupPtr ;		/* pointer to group title structures */
      ButtonTypePtr TextBoxPtr ;	/* pointer to text box structures */
      ButtonTypePtr SpecialPtr ;	/* pointer to special buttons */
      short       * PixelArray ;	/* pointer to pixel storage */
      short       * ActiveTextBox ;	/* indicator of the active text box */
      short         CursorPosition ;	/* cursor position within text box */
      short         WindowOpen ;	/* has window been opened ? */
      long          InWindowID ;	/* indicator for cursor location */
      long          BufferState ;	/* state of the fore and back buffers */
      short         CRUserSel ;		/* User selection to return for <CR> */
      short         TABtoNext ;		/* <TAB> goes to next text box ? */
   } NewWind , *NewWindPtr ;

#define DLOG_BOX_XPOS 0.05
#define DLOG_BOX_YPOS 0.05
#define DLOG_BOX_LENGTH 0.9
#define DLOG_BORDER_WIDTH 0.0025
#define DLOG_BUTTON_WIDTH 0.1
#define DLOG_BUTTON_HEIGHT 0.05
#define DLOG_RADIO_BUTTON_LENGTH 0.035

#define DLOG_BOX_RED    OMRED_RED	/* Red gun for message box background */
#define DLOG_BOX_GREEN  OMRED_GREEN	/* Green gun for msg box background */
#define DLOG_BOX_BLUE   OMRED_BLUE	/* Blue gun for msg box background */

#define DLOG_BOX_BORDER_RED   OMWHITE_RED	/* Red gun for msg box border */
#define DLOG_BOX_BORDER_GREEN OMWHITE_GREEN	/* Green gun for msg box bord */
#define DLOG_BOX_BORDER_BLUE  OMWHITE_BLUE	/* Blue gun for msg box bord */

#define DLOG_BACKGROUND_RED   OMBLUE_RED
#define DLOG_BACKGROUND_GREEN OMBLUE_GREEN
#define DLOG_BACKGROUND_BLUE  OMBLUE_BLUE

#define DLOG_BUTTON_TITLE_LENGTH 40
#define DLOG_MENU_TITLE_LENGTH 25
#define DLOG_MENU_COMMENT_LENGTH 80
#define DLOG_GROUP_TITLE_LENGTH 35

#define DLOG_RBUT_INDENT 0.05

#define DLOG_TEXT_BOX_HEIGHT 0.35
#define DLOG_TEXT_BOX_LENGTH 0.2
#define DLOG_TEXT_BOX_OFFX 0.2
#define DLOG_TEXT_ENTRY_LENGTH 8

enum DlogOptions {DLOG_MENU_INITIALISE, DLOG_MENU_USE, DLOG_MENU_DELETE};

#define CheckLength(String, Maxlength) if (strlen(String) > Maxlength) DoMessage("Internal error in string storage", PIMMS_FALSE)

#endif
