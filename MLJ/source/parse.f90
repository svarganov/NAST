MODULE M_strings
!-----------------------------------------------------------------------------------------------------------------------------------
PRIVATE
PUBLIC split      ! subroutine parses a string using specified delimiter characters and store tokens into an array
PUBLIC to_lower   ! function converts string to lowercase
!-----------------------------------------------------------------------------------------------------------------------------------
CONTAINS
!-----------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE split(input_line,array,delimiters,order,nulls)
!-----------------------------------------------------------------------------------------------------------------------------------
!  @(#) parse a string using specified delimiter characters and store tokens into an array
!-----------------------------------------------------------------------------------------------------------------------------------
   IMPLICIT NONE
   INTRINSIC INDEX, MIN, PRESENT, LEN
!-----------------------------------------------------------------------------------------------------------------------------------
!  given a line of structure " par1 par2 par3 ... parn " store each par(n) into a separate variable in array.
!    o by default  adjacent delimiters in the input string do not create an empty string in the output array
!    o no quoting of delimiters is supported
   CHARACTER(LEN=*),INTENT(IN)              :: input_line  ! input string to tokenize
   CHARACTER(LEN=*),OPTIONAL,INTENT(IN)     :: delimiters  ! list of delimiter characters
   CHARACTER(LEN=*),OPTIONAL,INTENT(IN)     :: order       ! order of output array SEQUENTIAL|[REVERSE|RIGHT]
   CHARACTER(LEN=*),OPTIONAL,INTENT(IN)     :: nulls       ! return strings composed of delimiters or not IGNORE|RETURN|IGNOREEND
   CHARACTER(LEN=*),ALLOCATABLE,INTENT(OUT) :: array(:)    ! output array of tokens
!    CHARACTER(LEN=*),INTENT(INOUT) :: array(:)
!-----------------------------------------------------------------------------------------------------------------------------------
   INTEGER                       :: n                      ! max number of strings INPUT_LINE could split into if all delimiter
   INTEGER,ALLOCATABLE           :: ibegin(:)              ! positions in input string where tokens start
   INTEGER,ALLOCATABLE           :: iterm(:)               ! positions in input string where tokens end
   CHARACTER(LEN=:),ALLOCATABLE  :: dlim                   ! string containing delimiter characters
   CHARACTER(LEN=:),ALLOCATABLE  :: ordr                   ! string containing order keyword
   CHARACTER(LEN=:),ALLOCATABLE  :: nlls                   ! string containing order keyword
   INTEGER                       :: ii,iiii                ! loop parameters used to control print order
   INTEGER                       :: icount                 ! number of tokens found
   INTEGER                       :: ilen                   ! length of input string with trailing spaces trimmed
   INTEGER                       :: i10,i20,i30            ! loop counters
   INTEGER                       :: icol                   ! pointer into input string as it is being parsed
   INTEGER                       :: idlim                  ! number of delimiter characters
   INTEGER                       :: ifound                 ! where next delimiter character is found in remaining input string data
   INTEGER                       :: inotnull               ! count strings not composed of delimiters
   INTEGER                       :: ireturn                ! number of tokens returned
   INTEGER                       :: imax                   ! length of longest token
!-----------------------------------------------------------------------------------------------------------------------------------
   ! decide on value for optional DELIMITERS parameter
   IF (PRESENT(delimiters)) THEN                                   ! optional delimiter list was present
      IF(delimiters.NE.'')THEN                                     ! if DELIMITERS was specified and not null use it
         dlim=delimiters
      ELSE                                                         ! DELIMITERS was specified on call as empty string
         dlim=' '//char(9)//char(10)//char(11)//char(12)//char(13)//char(0) ! use default delimiter when not specified
      ENDIF
   ELSE                                                            ! no delimiter value was specified
      dlim=' '//char(9)//char(10)//char(11)//char(12)//char(13)//char(0)    ! use default delimiter when not specified
   ENDIF
   idlim=LEN(dlim)                                                 ! dlim a lot of blanks on some machines if dlim is a big string
!-----------------------------------------------------------------------------------------------------------------------------------
   ! decide on value for optional ORDER parameter
   IF (PRESENT(order)) THEN                                        ! allocate optional parameter value for specifying output order
      ordr=to_lower(order)
   ELSE                                                            ! no delimiter value was specified
      ordr='sequential'
   ENDIF
!-----------------------------------------------------------------------------------------------------------------------------------
   IF (PRESENT(nulls)) THEN                                        ! allocate optional parameter value for specifying output order
      nlls=to_lower(nulls)
   ELSE                                                            ! no delimiter value was specified
      nlls='ignore'
   ENDIF
!-----------------------------------------------------------------------------------------------------------------------------------
   n=LEN(input_line)+1                        ! max number of strings INPUT_LINE could split into if all delimiter
   ALLOCATE(ibegin(n))                        ! allocate enough space to hold starting location of tokens if string all tokens
   ALLOCATE(iterm(n))                         ! allocate enough space to hold ending  location of tokens if string all tokens
   ibegin(:)=1
   iterm(:)=1
!-----------------------------------------------------------------------------------------------------------------------------------
   ilen=LEN(input_line)                                           ! ILEN is the column position of the last non-blank character
   icount=0                                                       ! how many tokens found
   inotnull=0                                                     ! how many tokens found not composed of delimiters
   imax=0                                                         ! length of longest token found
!-----------------------------------------------------------------------------------------------------------------------------------
   SELECT CASE (ilen)
!-----------------------------------------------------------------------------------------------------------------------------------
   CASE (:0)                                                      ! command was totally blank
!-----------------------------------------------------------------------------------------------------------------------------------
   CASE DEFAULT                                                   ! there is at least one non-delimiter in INPUT_LINE if get here
      icol=1                                                      ! initialize pointer into input line
      INFINITE: DO i30=1,ilen,1                                   ! store into each array element
         ibegin(i30)=icol                                         ! assume start new token on the character
         IF(INDEX(dlim(1:idlim),input_line(icol:icol)).eq.0)THEN  ! if current character is not a delimiter
            iterm(i30)=ilen                                       ! initially assume no more tokens
            DO i10=1,idlim                                        ! search for next delimiter
               ifound=INDEX(input_line(ibegin(i30):ilen),dlim(i10:i10))
               IF(ifound.GT.0)THEN
                  iterm(i30)=MIN(iterm(i30),ifound+ibegin(i30)-2)
               ENDIF
            ENDDO
            icol=iterm(i30)+2                                     ! next place to look as found end of this token
            inotnull=inotnull+1                                   ! increment count of number of tokens not composed of delimiters
         ELSE                                                     ! character is a delimiter for a null string
            iterm(i30)=icol-1                                     ! record assumed end of string. Will be less than beginning
            icol=icol+1                                           ! advance pointer into input string
         ENDIF
         imax=max(imax,iterm(i30)-ibegin(i30)+1)
         icount=i30                                               ! increment count of number of tokens found
         IF(icol.GT.ilen)THEN                                     ! text left
            EXIT INFINITE
         ENDIF
      enddo INFINITE
!-----------------------------------------------------------------------------------------------------------------------------------
   END SELECT
!-----------------------------------------------------------------------------------------------------------------------------------
   SELECT CASE (trim(adjustl(nlls)))
   CASE ('ignore','','ignoreend')
      ireturn=inotnull
   CASE DEFAULT
      ireturn=icount
   END SELECT
   ALLOCATE(array(ireturn))                                       ! allocate the array to turn
!-----------------------------------------------------------------------------------------------------------------------------------
   SELECT CASE (trim(adjustl(ordr)))                              ! decide which order to store tokens
   CASE ('reverse','right') ; ii=ireturn; iiii=-1                 ! last to first
   CASE DEFAULT             ; ii=1       ; iiii=1                 ! first to last
   END SELECT
!-----------------------------------------------------------------------------------------------------------------------------------
   DO i20=1,icount                                                ! fill the array with the tokens that were found
!     write(*,*) i20,'@'//input_line(ibegin(i20):iterm(i20))//'@',ibegin(i20),iterm(i20)
      IF(iterm(i20).LT.ibegin(i20))then
         SELECT CASE (trim(adjustl(nlls)))
         CASE ('ignore','','ignoreend')
         CASE DEFAULT
            array(ii)=' '
            ii=ii+iiii
         END SELECT
      ELSE
         array(ii)=input_line(ibegin(i20):iterm(i20))
         ii=ii+iiii
      ENDIF
   ENDDO
!-----------------------------------------------------------------------------------------------------------------------------------
   END SUBROUTINE split
!===================================================================================================================================

PURE FUNCTION to_lower(instr) result(outstr) ! @(#) function converts ASCII instr to lowercase
   IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
   CHARACTER(LEN=*), INTENT(IN) :: instr                          ! mixed-case input string to change
   CHARACTER(LEN=LEN(instr))    :: outstr                         ! lowercase output string to generate
!-----------------------------------------------------------------------------------------------------------------------------------
   INTEGER                      :: i10                            ! loop counter for stepping thru string
   INTEGER                      :: ade                            ! ASCII Decimal Equivalent of current character
   INTEGER,PARAMETER            :: ade_a=IACHAR('A')
   INTEGER,PARAMETER            :: ade_z=IACHAR('Z')
!-----------------------------------------------------------------------------------------------------------------------------------
   outstr=instr                                                   ! initially assume output string equals input string
!-----------------------------------------------------------------------------------------------------------------------------------
   stepthru: DO i10=1,LEN(instr)
      ade=IACHAR(instr(i10:i10))                                  ! convert letter to its value in ASCII collating sequence
      IF(ade .GE. ade_a .AND. ade .LE. ade_z ) THEN               ! if current letter is uppercase change it
         outstr(i10:i10)=ACHAR(ade+32)                            ! change letter to lowercase
      ENDIF
   ENDDO stepthru
!-----------------------------------------------------------------------------------------------------------------------------------
END FUNCTION to_lower

END MODULE M_strings
