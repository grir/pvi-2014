#include <cmath>
#include <fstream>
#include <iostream>

      const double ZERO  = 0.0;
      const double ONE   = 1.0;
      const double TWO   = 2.0;
      const double THREE = 3.0;
      const double FOUR  = 4.0;
      const double FIVE  = 5.0;
      const double SIX   = 6.0;
      const double SEVEN = 7.0;
      const double EIGHT = 8.0;
      const double NINE  = 9.0;
      const double TEN   = 10.0;
      const double HALF  = ONE/TWO;
//      END MODULE CONSTANTS

//      MODULE TOOLS
//
//     This module contains various utilities.
//
//      USE CONSTANTS
//      IMPLICIT NONE

//      CONTAINS
//
//      SUBROUTINE KEY_SORT(LIST,PERMUTATION)
      void KEY_SORT(int N, double* LIST, int* PERMUTATION){
//
//     This subroutine performs a key sort on the real LIST by the
//     heap sort method.  The returned permutation has the property that
//     LIST(PERMUTATION(I))<=LIST(PERMUTATION(I+1)) for all relevant I.
//     See: Nijenhuis A and Wilf HS (1978) "Combinatorial Algorithms for
//     Computers and Calculators, 2nd Ed.", Academic Press.
//
         int I,J,K,L,N,PSTAR;

      for(int i=0;i<N;i++)
         PERMUTATION(i) = i;

      if (N<=1) return;
//
//     Carry out the heap sort on the permutation key.
//
      L = 1+N/2;
      K = N;
      do {
         if (L>1){
            L = L-1;
            PSTAR = PERMUTATION[ L ];
         }
         else{
            PSTAR = PERMUTATION[K];
            PERMUTATION[K] = PERMUTATION[1];
            K = K-1;
            if (K==1){
               PERMUTATION[1] = PSTAR;
               return;
            }
         }

         I = L;
         J = L + L;
         while (J<=K){
            if (J<K){
               if ( LIST[ PERMUTATION[J]] < LIST[PERMUTATION[J+1]] )
                  J = J+1;
            }
            if ( LIST[ PSTAR ] < LIST[ PERMUTATION[ J ]] ){
               PERMUTATION(I) = PERMUTATION(J)
               I = J;
               J = J+J;
            }
            else
               J = K+1;

         }
         PERMUTATION[I] = PSTAR;
      } while(true);
   }

///////////////////////////////////////////////////////////////////
      double POISSON_LOGLIKELIHOOD(double mean, int k){

           return -mean + k * log(mean) - gamlog(k + ONE);

      }
///////////////////////////////////////////////////////////////////
      double POISSON_TAIL(double mean,int k){
           return STANDARD_GAMMA((double) k, mean);
      }
///////////////////////////////////////////////////////////////////

      double LOG_POISSON_TAIL(double mean, int k){
//
//     This function returns the log (base 10) of the Poisson tail
//     probability.
//
       if ( (mean >= ONE) && ((k-mean) > SIX * sqrt(mean))) {
          double a = k + ONE;
          double b = k + TWO;
          double c = -mean + k * log(mean) - GAMLOG(a) +
                     log(ONE + mean * b / (a * (b - mean)));
          return  c / log(TEN);
       }
       else
          return log10( POISSON_TAIL(mean,k));
       }

      double STANDARD_GAMMA(double a, double k){
//
//     This routine returns the gamma distribution function with shape
//     parameter A and scale parameter 1 at the point X.
//
          int N;
          double  res, an, bn, cn, dn, pn, sum, term;
//
          if ((k <= ZERO) || ( a <= ZERO))
             return ZERO;
          else
             if ( k <= a + ONE ) {
//
//     Use the power series expansion.
//
                 term = exp(-k + a * log(k) - GAMLOG(a + ONE));
                 res = term;
                 for(int i=1;i<=100;i++){
                   term = term * k /(a + i);
                   sum = res + term;
                   if ((term / sum) < 1.0E-8 )
                      return res;
                   res = sum;
                 }
             }
             else{
//
//     Use the continued fraction expansion.
//
            bn = k + ONE - a;
            cn = ONE / 1.E-30;
            dn = ONE / bn;
            res = exp(-k + a * log(k) - GAMLOG(a)) * dn;
            for(int i = 1; i<=100;i++){
               an = -i*(i-a);
               bn = bn + TWO;
               cn = cn + an / cn;
               if ( abs(cn) < 1.0E-30)
                   cn = 1.0E-30;
               dn = bn + an * dn;
               if (abs(dn) < 1.0E-30 )
                   dn = 1.0E-30;
               dn = ONE / dn;
               pn = cn * dn;
               res = res * pn;
               if (abs(pn - ONE) < 1.0E-8){
                  res = ONE - res;
                  break;
               }
            }
          }
         return  max(res, ZERO);
      }
/////////////////////////////////////////////////////////////////////
      double  GAMLOG(X){
//
//     This routine computes log(gamma(X)) via a recurrence relation
//     and Stirling's formula.
//

      int i,n;
      double f,res,x,y;
//
//     Compute gamma as a factorial for small integer arguments.
//
      y = x - ONE;
      n = floor(y);
      if (( (y-n) <= ZERO) && (n <= TWO * TEN) ) {
         y = max(n,1);
         for(int i=n-1;i<=2;i--)
             y = y * i;
         return  log(y);
//    }
//     Compute log gamma via a recurrence relation and Stirling's formula.
//
      }
      else {
         f = sqrt(eight * atan (ONE));
         while (y<TEN){
            y = y + ONE;
            f = f / y;
         }
         return  log(f) + ( y + ONE/TWO) * log(y) - y + ONE / (THREE * FOUR * y) - ONE / (TEN * SIX * SIX * y * y *y);

      }


      void FIT_GRAPH_MODEL( INCOMING_PROPENSITY, OUTGOING_PROPENSITY, ARC_COUNT, &
         LIST,OUTPUT_UNIT)
//
      IMPLICIT NONE
      INTEGER :: OUTPUT_UNIT
      CHARACTER(LEN=*), DIMENSION(:) :: LIST
      INTEGER, DIMENSION(:,:) :: ARC_COUNT
      REAL(KIND=DBLE), DIMENSION(:) :: INCOMING_PROPENSITY,OUTGOING_PROPENSITY
//
      INTEGER :: I,ITERATION,J,K,N,VERTICES
      REAL(KIND=DBLE) :: LOGLIK,OLDLOGLIK,MEAN,S,T
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: PAIR
      INTEGER, ALLOCATABLE, DIMENSION(:) :: INCOMING,OUTGOING,PERMUTATION
      REAL(KIND=DBLE), ALLOCATABLE, DIMENSION(:) :: P,Q
      REAL(KIND=DBLE), ALLOCATABLE, DIMENSION(:) :: LOG_PVALUE
//
//     Compute the number of vertices and allocate working arrays.
//
      VERTICES = SIZE(INCOMING_PROPENSITY)
      ALLOCATE(P(VERTICES),Q(VERTICES))
      ALLOCATE(INCOMING(VERTICES),OUTGOING(VERTICES))
//
//     Determine the numbers of arcs coming into each vertex going out of
//     each vertex.
//
      INCOMING = 0
      OUTGOING = 0
      DO K = 1,VERTICES
         DO J = 1,VERTICES
            IF (K==J) CYCLE
            OUTGOING(J) = OUTGOING(J)+ARC_COUNT(J,K)
            INCOMING(K) = INCOMING(K)+ARC_COUNT(J,K)
         END DO
      END DO
//
//     Initialize the propensities and allocate their temporary copies.
//
      INCOMING_PROPENSITY = ONE
      OUTGOING_PROPENSITY = ONE
//
//     Enter the main iteration loop.
//
      OLDLOGLIK = -HUGE(ZERO)
      DO ITERATION = 1,50
         P = OUTGOING_PROPENSITY
         Q = INCOMING_PROPENSITY
//
//     Update the loglikelihood.
//
         LOGLIK = ZERO
         DO K = 1,VERTICES
            DO J = 1,VERTICES
               IF (K==J) CYCLE
               MEAN = P(J)*Q(K)
               LOGLIK = LOGLIK+POISSON_LOGLIKELIHOOD(MEAN,ARC_COUNT(J,K))
            END DO
         END DO
//
//     Update the propensities.
//
         T = SUM(Q)
         DO I = 1,VERTICES
            OUTGOING_PROPENSITY(I) = MAX(SQRT(P(I)*OUTGOING(I)/(T-Q(I))),TEN**(-12))
         END DO
         S = SUM(P)
         DO I = 1,VERTICES
            INCOMING_PROPENSITY(I) = MAX(SQRT(Q(I)*INCOMING(I)/(S-P(I))),TEN**(-12))
         END DO
//
//     Output the iteration number and loglikelihood.
//
         PRINT*," ITERATION = ",ITERATION," LOGLIKELIHOOD = ",LOGLIK
//
//     Check for a decrease in the loglikelihood.
//
         IF (OLDLOGLIK>LOGLIK+TEN**(-5)) THEN
            PRINT*," WARNING: DECREASE IN LOGLIKLIHOOD "
         END IF
//
//     Check for convergence.
//
         IF (ABS(OLDLOGLIK-LOGLIK)<TEN**(-10)*(ABS(OLDLOGLIK)+ONE)) EXIT
         OLDLOGLIK = LOGLIK
      END DO
//
//     Deallocate and allocate arrays.
//
      DEALLOCATE(P,Q)
      DEALLOCATE(INCOMING,OUTGOING)
      N = 0
      DO K = 1,VERTICES
         DO J = 1,VERTICES
            IF (K==J.OR.ARC_COUNT(J,K)==0) CYCLE
            N = N+1
         END DO
      END DO
      ALLOCATE(PAIR(2,N),LOG_PVALUE(N))
      ALLOCATE(PERMUTATION(N))
//
//     Load the nontrivial pairs.
//
      I = 0
      DO K = 1,VERTICES
         DO J = 1,VERTICES
            IF (K==J.OR.ARC_COUNT(J,K)==0) CYCLE
            MEAN = OUTGOING_PROPENSITY(J)*INCOMING_PROPENSITY(K)
            I = I+1
            LOG_PVALUE(I) = LOG_POISSON_TAIL(MEAN,ARC_COUNT(J,K))
            PAIR(1,I) = J
            PAIR(2,I) = K
         END DO
      END DO
//
//     Output the most significant pairs.
//
      WRITE(OUTPUT_UNIT,'(/,A)') " Most significant pairs:"
      WRITE(OUTPUT_UNIT,'(/,A,/)') "    RANK  LOGPVALUE  OBSERVED     EXPECTED    PAIR"
      CALL KEY_SORT(LOG_PVALUE,PERMUTATION)
      DO I = 1,400
         K = PERMUTATION(I)
         MEAN = OUTGOING_PROPENSITY(PAIR(1,K))*INCOMING_PROPENSITY(PAIR(2,K))
         WRITE(OUTPUT_UNIT,'(3X,I4,1X,F10.4,2X,I6,2X,F12.4,4X,A,1X,A)') I,LOG_PVALUE(K), &
            ARC_COUNT(PAIR(1,K),PAIR(2,K)),MEAN,TRIM(LIST(PAIR(1,K))),TRIM(LIST(PAIR(2,K)))
      END DO
      END SUBROUTINE FIT_GRAPH_MODEL
//
      END MODULE TOOLS

      MODULE STRING_MANIPULATION
//
//     This module contains routines for manipulating strings.
//
      CONTAINS
//
      FUNCTION POSITION_IN_ALPHABET(LETTER)
//
//     This subroutine computes the position in the English alphabet of
//     the given LETTER.
//
      IMPLICIT NONE
      CHARACTER(LEN=1) :: LETTER,APOSTROPHE = "'"
      INTEGER :: I,POSITION_IN_ALPHABET
//
      I = ICHAR(LETTER)
      IF (LETTER==APOSTROPHE) THEN
         POSITION_IN_ALPHABET = I
      ELSE IF (I>=ICHAR('a').AND.I<=ICHAR('z')) THEN
         POSITION_IN_ALPHABET = I
      ELSE IF (I>=ICHAR('A').AND.I<=ICHAR('Z')) THEN
         POSITION_IN_ALPHABET = I-ICHAR('A')+ICHAR('a')
      ELSE
         POSITION_IN_ALPHABET = 0
      END IF
      END FUNCTION POSITION_IN_ALPHABET

      SUBROUTINE REMOVE_FORBIDDEN_CHARACTERS(PERMITTED,STRING)
//
//     This subroutine blanks out all non-permitted characters.  For
//     example, one might want to remove all punctuation from STRING.
//
      IMPLICIT NONE
      CHARACTER(LEN=*) :: PERMITTED,STRING
      INTEGER :: I
//
      DO I = 1,LEN(STRING)
         IF (INDEX(PERMITTED,STRING(I:I))==0) THEN
            STRING(I:I) = " "
         END IF
      END DO
      END SUBROUTINE REMOVE_FORBIDDEN_CHARACTERS

      SUBROUTINE REMOVE_APOSTROPHES(STRING)
//
//     This subroutine deletes apostrophe signs preceded  by a blank
//     or followed by a blank or a lower case d.
//
      IMPLICIT NONE
      CHARACTER(LEN=*) :: STRING
      INTEGER :: I,J
//
      DO I = 1,LEN(STRING)
         SELECT CASE(STRING(I:I))
         CASE("'")
            IF (I==1) THEN
               STRING(I:I) = " "
            END IF
            IF (I==LEN(STRING)) THEN
               STRING(I:I) = " "
            END IF
            IF (I>1) THEN
               J = POSITION_IN_ALPHABET(STRING(I-1:I-1))
               IF (J<ICHAR('a').OR.J>ICHAR('z')) STRING(I:I) = " "
            END IF
            IF (I<LEN(STRING)) THEN
               J = POSITION_IN_ALPHABET(STRING(I+1:I+1))
               IF (J==ICHAR('d')) THEN
                  STRING(I:I) = "e"
               ELSE IF (J<ICHAR('a').OR.J>ICHAR('z')) THEN
                   STRING(I:I) = " "
               END IF
            END IF
         END SELECT
      END DO
      END SUBROUTINE REMOVE_APOSTROPHES

      SUBROUTINE REPLACE(STRING,SUBSTITUTE,TARGET)
//
//     This subroutine replaces TARGET with SUBSTITUTE in STRING.
//
      IMPLICIT NONE
      CHARACTER(LEN=*) :: STRING,SUBSTITUTE,TARGET
      CHARACTER(LEN=LEN(STRING)) :: RIGHT
      INTEGER :: I,J
//
      J = 1
      DO
         I = INDEX(STRING(J:),TARGET)
         IF (I==0) RETURN
         I = I+J-1
         J = I+LEN(SUBSTITUTE)
         RIGHT = STRING(I+LEN(TARGET):)
         STRING(I:I+LEN(SUBSTITUTE)-1) = SUBSTITUTE
         STRING(I+LEN(SUBSTITUTE):) = RIGHT
      END DO
      END SUBROUTINE REPLACE

      SUBROUTINE PROCESS_LINE(LINE)
//
//     This subroutine processes a line, deleting extraneous characters
//     and replacing abbreviations whenever possible.
//
      IMPLICIT NONE
      CHARACTER(LEN=*) :: LINE
      CHARACTER(LEN=60) :: PERMITTED
      INTEGER :: I,J
//
//     Remove all forbidden characters.
//
      PERMITTED = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz.//?:;-' "
      CALL REMOVE_FORBIDDEN_CHARACTERS(PERMITTED,LINE)
//
//     There are many unusual abbreviations in Shakespeare.
//
      CALL REPLACE(LINE,"est","'st")
      CALL REPLACE(LINE,"to it","to't")
      CALL REPLACE(LINE,"taken","ta'en")
      CALL REPLACE(LINE,"the ","th' ")
      CALL REPLACE(LINE," it","'t")
      CALL REPLACE(LINE," on "," o' ")
      CALL REPLACE(LINE," in "," i' ")
      CALL REPLACE(LINE,"tis","'tis")
      CALL REPLACE(LINE," in "," 'n")
      CALL REPLACE(LINE,"eve","e'e")
      CALL REPLACE(LINE,"er","'r")
      CALL REPLACE(LINE,"en","'n")
      CALL REPLACE(LINE,"over","o'er")
//
//     Remove the remaining extraneous apostrophes.
//
      CALL REMOVE_APOSTROPHES(LINE)
      END SUBROUTINE PROCESS_LINE

      SUBROUTINE SORT_STRINGS(LIST)
//
//     This subroutine performs a heap sort on a list of strings.  See:
//     Nijenhuis A and Wilf HS (1978) "Combinatorial Algorithms for
//     Computers and Calculators, 2nd ed", Chapter 15, Academic Press.
//
      IMPLICIT NONE
      CHARACTER(LEN=*), DIMENSION(:), INTENT(INOUT) :: LIST

      INTEGER :: I,J,K,L,N
      CHARACTER(LEN=LEN(LIST(1))) :: TEMP
//
      N = SIZE(LIST)
      IF (N<=1) RETURN
      L = 1+N/2
      K = N
      DO
         IF (L>1) THEN
            L = L-1
            TEMP = LIST(L)
         ELSE
            TEMP = LIST(K)
            LIST(K) = LIST(1)
            K = K-1
            IF (K<=1) THEN
               LIST(1) = TEMP
               RETURN
            END IF
         END IF
         I = L
         J = L+L
         DO WHILE (J<=K)
            IF (J<K) THEN
               IF (LIST(J)<LIST(J+1)) J = J+1
            END IF
            IF (TEMP<LIST(J)) THEN
               LIST(I) = LIST(J)
               I = J
               J = J+J
            ELSE
               J = K+1
            END IF
         END DO
         LIST(I) = TEMP
      END DO
      END SUBROUTINE SORT_STRINGS

      SUBROUTINE PURGE_STRINGS(LIST,UNIQUE)
//
//     This subroutine purges the ordered string array LIST of duplicate
//     entries.  The number of UNIQUE entries is calculated in the process.
//     Note that UNIQUE returns with the value 1 if the LIST is empty.
//
      IMPLICIT NONE
      CHARACTER(LEN=*), DIMENSION(:) :: LIST
      INTEGER :: I,J,UNIQUE
//
      J = 1
      DO I = 2,SIZE(LIST)
         IF (LIST(I)/=LIST(J)) THEN
            J = J+1
            LIST(J) = LIST(I)
         END IF
      END DO
      UNIQUE = J
      END SUBROUTINE PURGE_STRINGS

      FUNCTION BISECT_STRING_LIST(LIST,ITEM)
//
//     This function returns the position of the particular ITEM in the
//     sorted string list.  The search is conducted by bisection.  The
//     user should check that the proposed position actually contains
//     the item.
//
      IMPLICIT NONE
      CHARACTER(LEN=*) :: ITEM
      INTEGER :: BISECT_STRING_LIST,FIRST,LAST,MID
      CHARACTER(LEN=*), DIMENSION(:) :: LIST
//
      FIRST = 1
      LAST = SIZE(LIST)
      DO
         IF (FIRST==LAST) THEN
            IF (ITEM==LIST(FIRST)) THEN
               BISECT_STRING_LIST = FIRST
            ELSE
               BISECT_STRING_LIST = 0
            END IF
            RETURN
         END IF
         MID = (FIRST+LAST)/2
         IF (ITEM<=LIST(MID)) THEN
            LAST = MID
         ELSE
            FIRST = MID+1
         END IF
      END DO
      END FUNCTION BISECT_STRING_LIST

      SUBROUTINE INPUT_DATA(TEXT_FILE,LINE,INPUT_UNIT,LINES)
//
//     This subroutine opens the text file and determines the number of
//     lines in it.
//
      IMPLICIT NONE
      CHARACTER(LEN=*) :: TEXT_FILE,LINE
      INTEGER :: I,INPUT_UNIT,IOERROR,LINES
//
//     Open the text file.
//
      OPEN(UNIT=INPUT_UNIT,FILE=TEXT_FILE)
//
//     Continue reading lines until the end of file is encountered.
//
      LINES = 0
      DO
         READ(INPUT_UNIT,'(1000A)',IOSTAT=IOERROR) (LINE(I:I),I=1,LEN(LINE))
         IF (IOERROR<0) THEN
            RETURN
         ELSE
            LINES = LINES+1
         END IF
      END DO
      END SUBROUTINE INPUT_DATA

      SUBROUTINE EXTEND_WORD_LIST(LINE,LIST,WORDS,ERROR,HYPHEN_ACTIVE)
//
//     This subroutine extracts the words from the current LINE and
//     adds them to the current word LIST.
//
      IMPLICIT NONE
      CHARACTER(LEN=*) :: LINE
      CHARACTER(LEN=100) :: WORD
      INTEGER :: I,J,K,WORDS
      LOGICAL :: ERROR,HYPHEN_ACTIVE
      CHARACTER(LEN=*), DIMENSION(:) :: LIST
//
//     Process the line character by character.  K is the length of
//     the current word.
//
      ERROR = .FALSE.
      WORD = " "
      K = 0
      DO I = 1,LEN(LINE)
         J = POSITION_IN_ALPHABET(LINE(I:I))
//
//     Extend the current word.
//
         IF (J>0) THEN
            K = K+1
            WORD(K:K) = CHAR(J)
         END IF
//
//     Add the current completed word to the list.
//
         IF (J==0.OR.I==LEN(LINE)) THEN
            IF (HYPHEN_ACTIVE) THEN
               LIST(WORDS) = TRIM(LIST(WORDS))//WORD(1:K)
               HYPHEN_ACTIVE = .FALSE.
            ELSE IF (K>0) THEN
               WORDS = WORDS+1
               IF (WORDS>SIZE(LIST)) THEN
                  ERROR = .TRUE.
                  RETURN
               END IF
               LIST(WORDS) = WORD(1:K)   //adds the word to the list
            END IF
            K = 0
         END IF
      END DO
//
//     Check if the last nonblank character is a hyphen.
//
      K = LEN_TRIM(LINE)
      IF (K>0) THEN
         HYPHEN_ACTIVE = LINE(K:K)=="-"
      ELSE
         HYPHEN_ACTIVE = .FALSE.
      END IF
      END SUBROUTINE EXTEND_WORD_LIST

      SUBROUTINE UPDATE_WORD_PAIR_COUNTS(LIST,LINE,WORD1,WORD2,WORDPAIR_COUNT, &
         WORDS,ERROR)
//
//     This subroutine extracts word pairs from the current LINE and
//     updates the word pair count for each pair encountered.
//
      IMPLICIT NONE
      CHARACTER(LEN=800) :: LINE
      CHARACTER(LEN=100) :: WORD1,WORD2
      INTEGER :: A,B,ERROR_POSITION,I,J,K,LAST,WORDS
      LOGICAL :: ERROR
      CHARACTER(LEN=*), DIMENSION(WORDS) :: LIST
      INTEGER, DIMENSION(WORDS,WORDS)::WORDPAIR_COUNT
//
//     Process the line character by character.  K is the length of
//     the current second word.
//
      K = LEN_TRIM(WORD2)
      ERROR = .FALSE.
      LAST = LEN_TRIM(LINE)
      DO I = 1,LAST
//
//     Extend the current second word.
//
         J = POSITION_IN_ALPHABET(LINE(I:I))
         IF (J>0) THEN
            K = K+1
            WORD2(K:K) = CHAR(J)
         END IF
//
//     If the last nonblank character of the line is a hyphen, then exit.
//
         IF (I==LAST.AND.LINE(I:I)=="-") RETURN
//
//     Otherwise, if the current character is not a letter, then the second
//     word has ended.
//
         IF (J==0.OR.I==LAST) THEN
//
//     Update the count for the current pair of words.
//
            IF (WORD2/=" ") THEN
               IF (WORD1/=" ") THEN
                  A = BISECT_STRING_LIST(LIST(1:WORDS),WORD1)
                  B = BISECT_STRING_LIST(LIST(1:WORDS),WORD2)
                  IF (A*B>0) THEN
                     WORDPAIR_COUNT(A,B) = WORDPAIR_COUNT(A,B)+1
                  ELSE
                     ERROR = .TRUE.
                     PRINT*," MATCH ERROR"," WORD1 =",WORD1," WORD2 =",WORD2
                     STOP
                  END IF
               END IF
//
//     Copy the second word into the first word and reset the position in the
//     second word.
//
               WORD1 = " "
               WORD1 = TRIM(WORD2)
               WORD2 = " "
               K = 0
            END IF
//
//     Check if the current character is a punctuation mark.  If so, reinitialize
//     both words.
//
            SELECT CASE(LINE(I:I))
            CASE(".","?","//",":",";",",")
               WORD1 = " "
               WORD2 = " "
               K = 0
            END SELECT
         END IF
      END DO
      END SUBROUTINE UPDATE_WORD_PAIR_COUNTS

      SUBROUTINE COUNT_LETTERPAIRS(LINE,WORD,LETTERPAIR_COUNT)
//
//     Count the letter pairs within the words in LINE.  WORD is the current
//     partial word.
//
      IMPLICIT NONE
      CHARACTER(LEN=*) :: LINE,WORD
      INTEGER :: I,J,K,L,LAST,N
      INTEGER, DIMENSION(:,:) :: LETTERPAIR_COUNT
//
//     Process the line character by character.  N is the length of
//     the current word.
//
      N = LEN_TRIM(WORD)
      LAST = LEN_TRIM(LINE)
      DO I = 1,LAST
//
//     Extend the current word.
//
         J = POSITION_IN_ALPHABET(LINE(I:I))
         IF (J>0) THEN
            N = N+1
            WORD(N:N) = CHAR(J)
         END IF
//
//     If the last nonblank character of the line is a hyphen, then exit with
//     a partial word.
//
         IF (I==LAST.AND.LINE(I:I)=="-") RETURN
//
//     Otherwise, if the current character is not a letter, then the current
//     word has ended.
//
         IF (J==0.OR.I==LAST) THEN
            IF (N>1) THEN
               IF (IACHAR(WORD(2:2))>IACHAR('a')) THEN
//
//     Update the letter pair counts and reinitialize the current word.
//
                  DO L = 1,N-1
                     J = POSITION_IN_ALPHABET(WORD(L:L))
                     IF (J==ICHAR("'")) THEN
                        J = 27
                     ELSE
                        J = J-ICHAR('a')+1
                     END IF
                     K = POSITION_IN_ALPHABET(WORD(L+1:L+1))
                     IF (K==ICHAR("'")) THEN
                        K = 27
                     ELSE
                        K = K-ICHAR('a')+1
                     END IF
                     LETTERPAIR_COUNT(J,K) = LETTERPAIR_COUNT(J,K)+1
                  END DO
               END IF
            END IF
            N = 0
            WORD = " "
         END IF
      END DO
      END SUBROUTINE COUNT_LETTERPAIRS
//
      END MODULE STRING_MANIPULATION
//
      PROGRAM WORD_PAIRS
//
//     This program implements one or more of the computational tools
//     for literary analysis.
//
      USE CONSTANTS
      USE TOOLS
      USE STRING_MANIPULATION
//
      IMPLICIT NONE
      CHARACTER(LEN=100) :: WORD1,WORD2
      CHARACTER(LEN=800) :: TEXT_FILE,LINE
      INTEGER :: I,IOERROR,INPUT_UNIT = 1,J,LINES,OUTPUT_UNIT = 2,WORDS
      LOGICAL :: ERROR,HYPHEN_ACTIVE
      CHARACTER(LEN=1), DIMENSION(27) :: ALPHABET
      CHARACTER(LEN=25), DIMENSION(100000) :: LIST
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: LETTERPAIR_COUNT,WORDPAIR_COUNT
      REAL(KIND=DBLE), ALLOCATABLE, DIMENSION(:) :: INCOMING_PROPENSITY,OUTGOING_PROPENSITY
//
//     Set the number of words, the hyphen alert, and the first line.
//
      WORDS = 0
      HYPHEN_ACTIVE = .FALSE.
//
//     Find the number of lines in the text.
//
      TEXT_FILE = "JaneEyre.txt"
      CALL INPUT_DATA(TEXT_FILE,LINE,INPUT_UNIT,LINES)
//
//     Read the text line by line and form the list of words.
//
      REWIND(INPUT_UNIT)
      DO J = 1,LINES
         LINE = " "
         READ(INPUT_UNIT,'(1000A)',IOSTAT=IOERROR) (LINE(I:I),I=1,LEN(LINE))
         CALL PROCESS_LINE(LINE)
         CALL EXTEND_WORD_LIST(LINE,LIST,WORDS,ERROR,HYPHEN_ACTIVE)
//
//     Every 100 lines sort and purge the list of words.
//
         IF (.NOT. HYPHEN_ACTIVE) THEN
            IF (MOD(J,100)==0.OR.J==LINES) THEN
               CALL SORT_STRINGS(LIST(1:WORDS))
               CALL PURGE_STRINGS(LIST(1:WORDS),WORDS)
            END IF
         END IF
      END DO
//
//     Print the alphabetized list.
//
      OPEN(UNIT=OUTPUT_UNIT,FILE="OUTPUT.txt")
      PRINT*," Unique words in file = ",WORDS
      WRITE(OUTPUT_UNIT,'(A,I6)') " Unique words in file = ",WORDS
      WRITE(OUTPUT_UNIT,'(/,A,/)') " File dictionary:"
      DO I = 1,WORDS
         WRITE(OUTPUT_UNIT,'(1X,A)') LIST(I)
      END DO
//
//     Set the dimensions of the word-pair count matrix.
//
      ALLOCATE(WORDPAIR_COUNT(WORDS,WORDS))
      WORDPAIR_COUNT = 0
      WORD1 = " "
      WORD2 = " "
//
//    Read in text line by line and update the word-pair count matrix.
//
      REWIND(INPUT_UNIT)
      DO J = 1,LINES
         LINE = " "
         READ(INPUT_UNIT,'(1000A)',IOSTAT=IOERROR) (LINE(I:I),I=1,LEN(LINE))
         CALL PROCESS_LINE(LINE)
         CALL UPDATE_WORD_PAIR_COUNTS(LIST,LINE,WORD1,WORD2,WORDPAIR_COUNT, &
         WORDS,ERROR)
      END DO
//
//     Output the total word pairs.
//
      J = SUM(WORDPAIR_COUNT)
      WRITE(OUTPUT_UNIT,'(/,A,I8)') " Word pair count in file = ",J
      PRINT*," Word pair count in file = ",J
//
      ALLOCATE(INCOMING_PROPENSITY(WORDS),OUTGOING_PROPENSITY(WORDS))
      CALL FIT_GRAPH_MODEL(INCOMING_PROPENSITY,OUTGOING_PROPENSITY,WORDPAIR_COUNT, &
         LIST,OUTPUT_UNIT)
//
//     Prepare for letter pair analysis.
//
      DEALLOCATE(WORDPAIR_COUNT)
      DEALLOCATE(INCOMING_PROPENSITY,OUTGOING_PROPENSITY)
      ALLOCATE(LETTERPAIR_COUNT(27,27))
      LETTERPAIR_COUNT = 0
//
      DO I = 1,26
         ALPHABET(I) = CHAR(ICHAR('a')+I-1)
      END DO
      ALPHABET(27) = "'"
      WORD1 = " "
      REWIND(INPUT_UNIT)
      DO J = 1,LINES
         LINE = " "
         READ(INPUT_UNIT,'(1000A)',IOSTAT=IOERROR) (LINE(I:I),I=1,LEN(LINE))
         CALL PROCESS_LINE(LINE)
         CALL COUNT_LETTERPAIRS(LINE,WORD1,LETTERPAIR_COUNT)
      END DO
//
//     Output the total letter pairs.
//
      J = SUM(LETTERPAIR_COUNT)
      WRITE(OUTPUT_UNIT,'(/,A,I8)') " Letter pair count in file = ",J
      PRINT*," Letter pair count in file = ",J
//
//     Estimate and output the graph parameters.
//
      ALLOCATE(INCOMING_PROPENSITY(27),OUTGOING_PROPENSITY(27))
      CALL FIT_GRAPH_MODEL(INCOMING_PROPENSITY,OUTGOING_PROPENSITY,LETTERPAIR_COUNT, &
         ALPHABET,OUTPUT_UNIT)
      PAUSE
      END PROGRAM WORD_PAIRS
