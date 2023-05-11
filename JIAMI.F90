!===============================================================================
!===============================================================================
    MODULE LDCSEM_DATA
!===============================================================================
      TYPE :: TYPE_ISINF
	   REAL(8) :: AX
	   REAL(8) :: AY
	   REAL(8) :: AZ
	   REAL(8) :: BX
	   REAL(8) :: BY
	   REAL(8) :: BZ
	   REAL(8) :: IS
	  END TYPE TYPE_ISINF
!===============================================================================
      TYPE :: TYPE_CDINF
	   INTEGER :: N
	   REAL(8) :: X
	   REAL(8) :: Y
	   REAL(8) :: Z
	  END TYPE TYPE_CDINF
!===============================================================================
      TYPE :: TYPE_MODINF
	   REAL(8) :: H
	   REAL(8) :: R
	  END TYPE TYPE_MODINF
!===============================================================================
      TYPE :: TYPE_EHINF
	   REAL(8) :: F
	   REAL(8) :: EX
	   REAL(8) :: EY
	   REAL(8) :: EZ
	   REAL(8) :: HX
	   REAL(8) :: HY
	   REAL(8) :: HZ
	   DOUBLE COMPLEX :: CEX
	   DOUBLE COMPLEX :: CEY
	   DOUBLE COMPLEX :: CEZ
	   DOUBLE COMPLEX :: CHX
	   DOUBLE COMPLEX :: CHY
	   DOUBLE COMPLEX :: CHZ
	  END TYPE TYPE_EHINF
!===============================================================================
    CONTAINS
	END MODULE
!===============================================================================
!===============================================================================


!===============================================================================
      SUBROUTINE CX_YZJC()
	  USE DFPORT
	  IMPLICIT NONE
      INTEGER    :: IsSucc
!===============================================================================
      OPEN(23,FILE="反演过程\1.TXT",ERR=1123)
	  CLOSE(23,STATUS = 'DELETE' )
	  GOTO 1112
1123  IsSucc=system("md 反演过程" )
1112  CONTINUE
      OPEN(23,FILE="反演模型\1.TXT",ERR=1125)
	  CLOSE(23,STATUS = 'DELETE' )
	  GOTO 1114
1125  IsSucc=system("md 反演模型" )
1114  CONTINUE
      OPEN(23,FILE="数据拟合\1.TXT",ERR=1127)
	  CLOSE(23,STATUS = 'DELETE' )
	  GOTO 1116
1127  IsSucc=system("md 数据拟合" )
1116  CONTINUE
!==========================================================================
      RETURN
!===========================================================================	   
    END SUBROUTINE 