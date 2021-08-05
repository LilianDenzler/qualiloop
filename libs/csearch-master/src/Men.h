#ifndef __MEN__
#define __MEN__

#include "Colours.h"
#include "ButtonProto.h"

/* define main menu colours */
#define MEN_BACKGROUND_RED	OMRED_RED
#define MEN_BACKGROUND_GREEN	OMRED_GREEN
#define MEN_BACKGROUND_BLUE	OMRED_BLUE

#define MEN_FOREGROUND_RED	OMBLUE_RED
#define MEN_FOREGROUND_GREEN	OMBLUE_GREEN
#define MEN_FOREGROUND_BLUE	OMBLUE_BLUE

#define MEN_BORDER_RED		OMWHITE_RED
#define MEN_BORDER_GREEN	OMWHITE_GREEN
#define MEN_BORDER_BLUE		OMWHITE_BLUE

#define MEN_TEXT_RED		OMWHITE_RED
#define MEN_TEXT_GREEN		OMWHITE_GREEN
#define MEN_TEXT_BLUE		OMWHITE_BLUE

#define MEN_BORDER_WIDTH	0.0025

#define MEN_TITLE_HEIGHT	0.025

/* structure to hold a rectangle coordinates */
typedef struct MenRectangle
{
   short xleft;
   short ybottom;
   short xright;
   short ytop;
} MenRectangle, *MenRectanglePtr;

typedef struct MenButton
{
   struct MenButton *forlink;  /* pointer to next element */
   struct MenButton *backlink; /* pointer to previous element */
   ButtonType Button;                /* pointer to button structure */
   int ButtonCode;                   /* code associated with button */
} MenButton, *MenButtonPtr;

typedef struct MenButtonHeader
{
   struct MenButtonHeader *forlink;  /* pointer to next structure in list */
   struct MenButtonHeader *backlink; /* pointer to previous structure */
   MenButtonPtr ButtonPointer;       /* pointer to first button in list */
   int CurrentButton;                /* current button in the group */
   char *GroupTitle;                 /* title for the button group */
   short TitleX;                     /* position of group title */
   short TitleY;                     /* position of group title */
} MenButtonHeader, *MenButtonHeaderPtr;

typedef struct MenMenuType
{
   long windowID;                  /* window identifier */
   char *title;                    /* the window title */
   short titlex, titley;           /* the title position */
   short xpos;                     /* menu x position */
   short ypos;                     /* menu y position */
   short length;                   /* menu length */
   short height;                   /* menu height */
   RGBColour backgroundcolour;     /* menu background colour */
   RGBColour foregroundcolour;     /* menu foreground colour */
   RGBColour bordercolour;         /* menu border colour */
   RGBColour textcolour;           /* colour for text */
   MenRectangle borderposition;    /* position of border line */
   MenRectangle foregroundposition;/* position of foreground rectangle */
   MenRectangle titleborder;       /* border round title */
   MenRectangle titleunderline;    /* line below title */
   MenButtonHeaderPtr SquareButPtr;      /* pointer to list of square buttons */
   MenButtonPtr TextButPtr;        /* pointer to list of text entry buttons */
   MenButtonPtr ScrollButPtr;      /* pointer to list of scrolled buttons */
   MenButtonPtr LabelButPtr;       /* pointer to list of label buttons */
   MenButtonHeaderPtr RadioButPtr; /* pointer to list of radio buttons */
   MenButtonPtr CurrentTextButton; /* pointer to current text button */
} MenMenuType, *MenMenuTypePtr;

#endif
