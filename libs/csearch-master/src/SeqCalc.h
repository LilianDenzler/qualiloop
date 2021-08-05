#ifndef __SEQCAL_H__   
#define __SEQCAL_H__   

/***********************************************************************
 *      Name: SeqCalc.h                                                *
 *  Function: Declares constants for calculation routines              *
 * Copyright: (C) Oxford Molecular Limited                             *
 *---------------------------------------------------------------------*
 *    Author: KY Cockwell, Oxford Molecular Ltd.                       *
 *      Date: 01/09/1993                                               *
 *---------------------------------------------------------------------*
 * Modification Record                                                 *
 * Date     Inits   Comments                                           *
 * 29/11/93 DW      Motif review (CR 2773)                             *
 ***********************************************************************/

#include <stdlib.h>
#include <sys/types.h>
#include <malloc.h>
#include "SeqU.h"

/* Define some types for the sequence pointers etc. */

typedef struct motifstruct MOTType,    *MOTPtr;

struct motifstruct
{
   MOTPtr    next;
   MOTPtr    previous;
   short     length;
   short     start;
   char      accession[20];
};

#define   SEQCALC_OK        0 /* OK */
#define   SEQCALC_BRACKET   1 /* Unmatched bracket */
#define   SEQCALC_TOOLONG   2 /* Motif too long */
#define   SEQCALC_MIXED     3 /* Mixed brackets */

#define   MAXMATRIX 26
#define   ARRSIZE   30  /* max length of motif, allow extra for brackets */
#define   SSWINDOW  17  /* length of ss claculation window */

#define   SS_NOVAL          0 /* no secondary structure */
#define   SS_HELIX          1 /* helix secondary structure */
#define   SS_BETA           2 /* beta sheet secondary structure */
#define   SS_COIL           3 /* coil secondary structure */
#define   SS_TURN           4 /* beta turn secondary structure */

#define   HELIX_WEIGH     -40 /* weighting for helix scores */
#define   SHEET_WEIGH     -50 /* weighting for sheet scores */
#define   TURN_WEIGH      -30 /* weighting for turn scores */
#define   COIL_WEIGH        0 /* weighting for coil scores */

/* Motif search definitions */
#define SEQ_MOTCHR "(){}?"
#define SEQ_ORS    '('
#define SEQ_ORE    ')'
#define SEQ_NOTS   '{'
#define SEQ_NOTE   '}'
#define SEQ_ANY    '?'

#endif
