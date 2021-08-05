#ifndef __SEQIOSEQ_H__
#define __SEQIOSEQ_H__

/***********************************************************************
 *      Name: SeqIOseq.h                                               *
 *  Function: Declares constants for IO formats                        *
 * Copyright: (C) Oxford Molecular Limited                             *
 *---------------------------------------------------------------------*
 *    Author: KY Cockwell, Oxford Molecular Ltd.                       *
 *      Date: 05/07/19931                                              *
 *---------------------------------------------------------------------*
 * Modification Record                                                 *
 * Date     Inits   Comments                                           *
 * 26/11/93 DW      Sequence/alignment I/O review (CR 2757)            *
 ***********************************************************************/

#define SEQFT_PIR       0 /* PIR/OWL/NBRF format file */
#define SEQFT_EMBL      1 /* EMBL format file */
#define SEQFT_SWISS     2 /* SWISSPROT format file */
#define ALNFT_MULTAL    3 /* MULTAL format alignment file */
#define ALNFT_GCG       4 /* GCG format alignment file */
#define ALNFT_SOMAP     5 /* SOMAP format alignment file */
#define ALNFT_ASCII     6 /* ASCII format alignment file */
#define ALNFT_PS        7 /* PostScript format alignment file */
#define ALNFT_EPS       8 /* Encapsulated postscript format alignment file */
#define ALNFT_CLST      9 /* Clustal format alignment file */

#define SEQIO_NOTALN   -1 /* Sequence not aligned completely (Multal only)*/
#define SEQIO_OK        0 /* OK */
#define SEQIO_OPEN      1 /* Open error */
#define SEQIO_READ      2 /* Read Error */
#define SEQIO_WRITE     3 /* Write error */
#define SEQIO_CORRUPT   4 /* Corrupt Entry */
#define SEQIO_INVALFMT  5 /* Invalid format */
#define SEQIO_NULLINP   6 /* NULL input */
#define SEQIO_EXISTS    7 /* File Exists */
#define SEQIO_ALLOC     8 /* Problems allocating memory */

/* constants for postscript dump */
#define MARGIN    10
#define BOXWIDTH  11
#define BOXHEIGHT 11
#define OFFSET    9
#define ALIGNX    90

/* PostScript output flags */
#define PS_NORM    0
#define PS_SOLID   1
#define PS_MONO    2
#define PS_MONOBOX 3

#define MAXLINE      255

/* Sequence I/O parameters */
#define SEQ_PIREND   "*"
#define SEQ_SWISSEND "//"

/* Alignment I/O parameters: Accession code and sequence start positions */
#define SEQ_SPECIAL "-."

#define SOMAP_AP    0
#define SOMAP_AE    13
#define SOMAP_SS    16

#define CLUST_AP    0
#define CLUST_AE    15
#define CLUST_SS    16

#define GCG_AP      0
#define GCG_AE      11
#define GCG_SS      12

/* CGF I/O */
#define CGF_SYMBOLS "~!@#$%^&*:\\.+=<>_,;/'|"

#endif
