C
C ****: FILES
C
C -- Function:-
C      Set up parameters etc for file handling.
C      Define parameters to be passed to 'FILOPN','FILCLS'
C
C
C
C
C -- Parameters for min and max channel numbers
C
C -- Minimum channel number
      INTEGER FCHMIN
      PARAMETER ( FCHMIN = 0 )
C -- Maximum channel number
      INTEGER FCHMAX
      PARAMETER ( FCHMAX = 128 )
C
C -- Parameters for type of file
C
C -- Starting position
      INTEGER FFRMMN
      PARAMETER ( FFRMMN = 100 )
C
C -- List carriage control
      INTEGER FFRMLC
      PARAMETER ( FFRMLC = FFRMMN + 0 )
C -- Formatted FORTRAN
      INTEGER FFRMFC
      PARAMETER ( FFRMFC = FFRMMN + 1 )
C -- Unformatted FORTRAN
      INTEGER FFRMUF
      PARAMETER ( FFRMUF = FFRMMN + 2 )
C -- Formatted direct access
      INTEGER FFRMDA
      PARAMETER ( FFRMDA = FFRMMN + 3 )
C -- Unformatted direct access
      INTEGER FFRMDU
      PARAMETER ( FFRMDU = FFRMMN + 4 )
C
C -- Maximum value
      INTEGER FFRMMX
      PARAMETER ( FFRMMX = FFRMDU )
C
C -- Parameters for access
C
C -- Starting position
      INTEGER FACCMN
      PARAMETER ( FACCMN = 200 )
C
C -- Readonly
      INTEGER FACCRD
      PARAMETER ( FACCRD = FACCMN + 0 )
C -- Create
      INTEGER FACCWT
      PARAMETER ( FACCWT = FACCMN + 1 )
C -- Modify
      INTEGER FACCRW
      PARAMETER ( FACCRW = FACCMN + 2 )
C -- Append
      INTEGER FACCAP
      PARAMETER ( FACCAP = FACCMN + 3 )
C -- Scratch
      INTEGER FACCSC
      PARAMETER ( FACCSC = FACCMN + 4 )
C
C -- Maximum value
      INTEGER FACCMX
      PARAMETER ( FACCMX = FACCSC )
C
C
C -- Parameters for close
C
C -- Starting position
      INTEGER FCLSMN
      PARAMETER ( FCLSMN = 300 )
C
C -- Keep
      INTEGER FCLSKP
      PARAMETER ( FCLSKP = FCLSMN + 0 )
C -- Delete
      INTEGER FCLSDL
      PARAMETER ( FCLSDL = FCLSMN + 1 )
C
C -- Maximum value
      INTEGER FCLSMX
      PARAMETER ( FCLSMX = FCLSDL )
C
C -- Parameters passed to FILWRT, FILREA for (un)formatted
C
C -- Read/write formatted
      INTEGER FRWFRM
      PARAMETER ( FRWFRM = 500 )
C -- Read/write unformatted
      INTEGER FRWUNF
      PARAMETER ( FRWUNF = 501 )
C
C -- Parameters passed to FILWRT, FILREA for (no)verify
C
C -- Verify I/O
      INTEGER FVERIF
      PARAMETER ( FVERIF = 600 )
C -- Don't verify I/O
      INTEGER FNOVER
      PARAMETER ( FNOVER = 601 )
C
C
