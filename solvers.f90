module solvers
use mpi
use parameters

implicit none


contains
!***********************************************************************
!************ Matrix Diagonalization/Inversion Subroutine **************
!***********************************************************************
subroutine MatrixDiagonalization(Eigenv,Pmat,PmatInv,Matrice,NN)
!Perform Matrix Diagonalization: Matrice=Pmat*Lambda*PmatInv
!Inputs: Matrice, NN(dimension of Matrice)
!Outputs are: Eigenv,Pmat,PmatInv
!Eigenv is the diagonal of Lambda
implicit none 

real(kind=prec)::Eigenv(1:NN)
real(kind=prec)::Pmat(1:NN,1:NN)
real(kind=prec)::PmatInv(1:NN,1:NN)
real(kind=prec)::Matrice(1:NN,1:NN)
real(kind=prec),allocatable,dimension(:,:)::AA
integer(kind=prec)::NN
integer::LDA, LDVL, LDVR
INTEGER::LWMAX
!PARAMETER        ( LWMAX = 4*NN )	
!....Local Scalars ..	
INTEGER          INFO, LWORK
!....Local Arrays ..
DOUBLE PRECISION:: VL(NN,NN), VR(NN,NN), WR(NN), WI(NN)
double precision,allocatable,dimension(:)::WORK
!....External Subroutines ..
EXTERNAL         DGEEV
!....Intrinsic Functions ..
INTRINSIC        INT, MIN
!....Query the optimal workspace.
	
	LWMAX = 4*NN 
	
allocate(AA(1:NN,1:NN))
allocate(WORK(LWMAX))	
AA=Matrice
	
LDA=NN
LDVL=NN
LDVR=NN
	
LWORK = -1
CALL DGEEV( 'Vectors', 'Vectors', NN, AA, LDA, WR, WI, VL, LDVL,VR, LDVR, WORK, LWORK, INFO )
LWORK = MIN( LWMAX, INT(WORK( 1 )))	
	
!....Solve eigenproblem.

CALL DGEEV( 'Vectors', 'Vectors', NN, AA, LDA, WR, WI, VL, LDVL,VR, LDVR, WORK, LWORK, INFO )

!....Check for convergence.

IF( INFO.GT.0 ) THEN
WRITE(*,*)'The algorithm failed to compute eigenvalues.'
STOP
END IF

!....Right-Eigenvector matrix
Pmat=VR

!....Left-Eigenvector matrix
call InvertiMatrice(PmatInv,Pmat,NN)

!....Eigenvalues
Eigenv=WR

deallocate(AA)	 
deallocate(work)	

return
end subroutine
!***********************************************************************
!***********************************************************************
subroutine InvertiMatrice(MatriceInv,Matrice,NN)
implicit none
!Performs  calculation of the inverse of Matrice (1:NN)x(1:NN)
real(kind=prec)::Matrice(1:NN,1:NN),MatriceInv(1:NN,1:NN)
real(kind=prec), allocatable,dimension(:)::work! work array for LAPACK
integer,dimension(1:NN)::ipiv   ! pivot indices
integer(kind=prec)::NN
integer::LWORK
integer::info
!....External procedures defined in LAPACK
external DGETRF
external DGETRI
LWORK=4*NN

allocate(work(LWORK))
!NN=int(Ny-1)
!....Store Qmat in QmatInv to prevent it from being overwritten by LAPACK
MatriceInv=Matrice


!....DGETRF computes an LU factorization of a general M-by-N matrix A
!using partial pivoting with row interchanges.
call DGETRF(NN,NN, MatriceInv, NN, ipiv, info)

if (info /= 0) then
stop 'Matrix is numerically singular!'
end if

!....DGETRI computes the inverse of a matrix using the LU factorization
!computed by DGETRF.
call DGETRI(NN, MatriceInv, NN, ipiv, work, LWORK, info)

if (info /= 0) then
stop 'Matrix inversion failed!'
end if


deallocate(work)
return
end subroutine
!***********************************************************************
!***********************************************************************
subroutine matinv2(B,A)
!Performs a direct calculation of the inverse of a 2×2 matrix.
real(kind=prec):: A(2,2)   !Matrix
real(kind=prec):: B(2,2)   !Inverse matrix
real(kind=prec):: detinv
	
!....Calculate the inverse determinant of the matrix
detinv = 1/(A(1,1)*A(2,2) - A(1,2)*A(2,1))

!....Calculate the inverse of the matrix
B(1,1) = +detinv * A(2,2)
B(2,1) = -detinv * A(2,1)
B(1,2) = -detinv * A(1,2)
B(2,2) = +detinv * A(1,1)

return
end subroutine
!***********************************************************************
!***********************************************************************
end module solvers

