module transformspost
  use, intrinsic :: iso_c_binding
use mpi
use parameterspost
include 'fftw3-mpi.f03'

! use, intrinsic :: iso_c_binding
! use mpi
! use parameters
! include 'fftw3-mpi.f03'

integer(kind=8) :: fftplan_fwd_y
integer(kind=8) :: fftplan_bwd_y
integer(kind=8) :: fftplan_fwd_2d
integer(kind=8) :: fftplan_bwd_2d


!....Variables for Dicrete Chebyshev Transforms (DCT_1D_FWD and DCT_1D_BWD)
!....FWD stands fo Forward
!....BWD stands for Backward
real(kind=prec):: corrFWD
integer(C_SIZE_T) :: SizeArr
type(C_PTR):: PLAN_CHEB_FWD, PLAN_CHEB_BWD
type(C_PTR) :: rinFWD,rinBWD
real(C_DOUBLE),pointer::datarFWD(:,:,:),datarBWD(:,:,:)
integer(C_INT):: DimIn(1), howmany
!....Variables for 2D Discrete Fourier Transforms 
integer(C_FFTW_R2R_KIND):: kindT(1)
real(kind=prec):: corrFFT2D
type(C_PTR)::PLAN_FFT2D_FWD, PLAN_FFT2D_BWD
type(C_PTR) :: coutF_FWD, coutF_BWD
complex(C_DOUBLE_COMPLEX),pointer :: datacF_FWD(:,:,:), datacF_BWD(:,:,:)
integer(C_INTPTR_T)::alloc_local, local_M, local_j_offset
integer(C_PTRDIFF_T):: DimInFFT(2), howmanyFFT



contains

!***********************************************************************
subroutine transforms_init()
!....Initialize Discrete Fourier Transforms and Discrete Chebyshev Transforms
implicit none

  fftplan_fwd_y = 0
  fftplan_bwd_y = 0
  fftplan_fwd_2d = 0
  fftplan_bwd_2d = 0


!.... Create plans once for all...
!....DCTs
DimIn=(/int(Ny+1)/)
howmany=int(Nx*Nzloc)
SizeArr=int(Nx*Nzloc*(Ny+1))
kindT(1) = FFTW_REDFT00

!....FWD
corrFWD=1.0d0/(Ny)

rinFWD = fftw_alloc_real(SizeArr)
call c_f_pointer(rinFWD, datarFWD, [int(Ny+1),int(Nx),int(Nzloc)])

PLAN_CHEB_FWD = fftw_plan_many_r2r(1, DimIn, howmany,&
                             datarFWD, DimIn,1,&
                             int(Ny+1), datarFWD, DimIn,1,&
                             int(Ny+1), kindT, FFTW_ESTIMATE)


!....BWD
rinBWD = fftw_alloc_real(SizeArr)
call c_f_pointer(rinBWD, datarBWD, [int(Ny+1),int(Nx),int(Nzloc)])

PLAN_CHEB_BWD = fftw_plan_many_r2r(1, DimIn, howmany,&
                             datarBWD, DimIn,1,&
                             int(Ny+1), datarBWD, DimIn,1,&
                             int(Ny+1), kindT, FFTW_ESTIMATE)



!....FFTs
!....Correction
 corrFFT2D=1.0d0/(Nx*Nz)
!....Dimensions and howmany 
 DimInFFT  = (/int(Nz),int(Nx)/) 	
 howmanyFFT =  int(Ny+1)	



!....FWD and BWD
!....get local data size and allocate (note dimension reversal)
 alloc_local = fftw_mpi_local_size_many(int(2), DimInFFT, howmanyFFT,&
 FFTW_MPI_DEFAULT_BLOCK , MPI_COMM_WORLD, local_M, local_j_offset)
 
 
 coutF_FWD = fftw_alloc_complex(alloc_local)
 coutF_BWD = fftw_alloc_complex(alloc_local)
 
 
 call c_f_pointer(coutF_FWD, datacF_FWD, [int(Ny+1),int(Nx),int(local_M)])
 !watch out: the lower bound of datacF_FWD in the three dimensions is 1 (not 0)
 
  
!....create MPI plan for in-place forward DFT (note dimension reversal)   
 PLAN_FFT2D_FWD=fftw_mpi_plan_many_dft(int(2),DimInFFT,howmanyFFT,&
                                     FFTW_MPI_DEFAULT_BLOCK ,&
                                     FFTW_MPI_DEFAULT_BLOCK ,&
                                     datacF_FWD,datacF_FWD,&
                                     MPI_COMM_WORLD,&
                                     FFTW_FORWARD,&
                                     FFTW_ESTIMATE)

!....create MPI plan for in-place backward DFT (note dimension reversal)   
 call c_f_pointer(coutF_BWD, datacF_BWD, [int(Ny+1),int(Nx),int(local_M)])
 !watch out: the lower bound of datacF_BWD in the three dimensions is 1 (not 0)
     
     
 PLAN_FFT2D_BWD=fftw_mpi_plan_many_dft(int(2),DimInFFT,howmanyFFT,&
                                     FFTW_MPI_DEFAULT_BLOCK ,&
                                     FFTW_MPI_DEFAULT_BLOCK ,&
                                     datacF_BWD,datacF_BWD,&
                                     MPI_COMM_WORLD,&
                                     FFTW_BACKWARD,&
                                     FFTW_ESTIMATE)










return
end subroutine
!***********************************************************************
subroutine transforms_finalization()
!....Initialize Discrete Fourier Transforms and Discrete Chebyshev Transforms
implicit none

!....DCTs

!.....FWD
 call fftw_destroy_plan(PLAN_CHEB_FWD)
 call fftw_free(rinFWD)

!.....BWD
 call fftw_destroy_plan(PLAN_CHEB_BWD)
 call fftw_free(rinBWD)




!....FFTs

!.....FWD
 call fftw_destroy_plan(PLAN_FFT2D_FWD)
 call fftw_free(coutF_FWD)

!.....BWD
 call fftw_destroy_plan(PLAN_FFT2D_BWD)
 call fftw_free(coutF_BWD)




return
end subroutine
!***********************************************************************
subroutine trasformChebFourier(uhatcomp,ucomp)
implicit none
complex(kind=prec)::uhatcomp(0:Ny,0:Nx-1,0:Nzloc-1)
real(kind=prec)::ucomp(0:Ny,0:Nx-1,0:Nzloc-1)
real(kind=prec),allocatable,dimension(:,:,:)::ucheb
integer(kind=prec)::ii,jj
allocate(ucheb(0:Ny,0:Nx-1,0:Nzloc-1))


!...Chebyshev Transform in y
call DCT_1D_FWD(ucheb(0:Ny,0:Nx-1,0:Nzloc-1),ucomp(0:Ny,0:Nx-1,0:Nzloc-1))



!...Fourier Transform in x and z
call transform_rows_batched(uhatcomp(0:Ny,0:Nx-1,0:Nzloc-1),ucheb(0:Ny,0:Nx-1,0:Nzloc-1))





deallocate(ucheb)

return
end subroutine
!***********************************************************************
subroutine trasformChebFourier_inv(ucomp,uhatcomp)
implicit none
complex(kind=prec)::uhatcomp(0:Ny,0:Nx-1,0:Nzloc-1)
real(kind=prec)::ucomp(0:Ny,0:Nx-1,0:Nzloc-1)
real(kind=prec),allocatable,dimension(:,:,:)::ucheb
!complex(kind=prec),allocatable,dimension(:,:,:)::uchebT
!complex(kind=prec),allocatable,dimension(:,:,:)::TempC
integer(kind=prec)::ii,jj




allocate(ucheb(0:Ny,0:Nx-1,0:Nzloc-1))
!allocate(uchebT(0:Nx-1,0:Nz-1,0:Nyloc-1))
!allocate(TempC(0:Nx-1,0:Nz-1,0:Nyloc-1))


ucheb=0.0d0



!...Inverse Fourier Transform in y and z
call transformInv_rows_batched(ucheb(0:Ny,0:Nx-1,0:Nzloc-1),uhatcomp(0:Ny,0:Nx-1,0:Nzloc-1))



!...Inverse Chebyshev Transform in y
call DCT_1D_BWD(ucomp(0:Ny,0:Nx-1,0:Nzloc-1),ucheb(0:Ny,0:Nx-1,0:Nzloc-1))



!deallocate(temp,temp2)
!deallocate(TempC)
!deallocate(uchebT)
deallocate(ucheb)
return
end subroutine
!***********************************************************************
subroutine DCT_1D_FWD(output,ucomp)
!....Discrete Chebyshev Transform 1D forward
!....It uses FFTW for Cosine trasforms (real to real transform)
implicit none
real(kind=prec)::ucomp(0:Ny,0:Nx-1,0:Nzloc-1)
real(kind=prec)::output(0:Ny,0:Nx-1,0:Nzloc-1)


 datarFWD=ucomp

!call MPI_BARRIER(MPI_Comm_World,ierror) 

!....compute transform Real to Real (Cosine Transform)
 call fftw_execute_r2r(PLAN_CHEB_FWD, datarFWD, datarFWD)
  
 output=datarFWD

 call MPI_BARRIER(MPI_Comm_World,ierror) 
 
 
output=corrFWD*output
output(0,:,:)=output(0,:,:)/2.0d0
output(Ny,:,:)=output(Ny,:,:)/2.0d0




return
end subroutine
!***********************************************************************
subroutine DCT_1D_BWD(output,ucomp)
!....Discrete Chebyshev Transform 1D forward
!....It uses FFTW for Cosine trasforms (real to real transform)
implicit none
real(kind=prec)::ucomp(0:Ny,0:Nx-1,0:Nzloc-1)
real(kind=prec)::output(0:Ny,0:Nx-1,0:Nzloc-1)



 datarBWD(1:Ny+1,1:Nx,1:Nzloc)=ucomp(0:Ny,0:Nx-1,0:Nzloc-1)
 datarBWD(2:Ny,:,:)=datarBWD(2:Ny,:,:)/2.0d0

!call MPI_BARRIER(MPI_Comm_World,ierror) 

!....compute transform Real to Real (Cosine Transform)
 call fftw_execute_r2r(PLAN_CHEB_BWD, datarBWD, datarBWD)
    
  
 output(0:Ny,0:Nx-1,0:Nzloc-1)=datarBWD(1:Ny+1,1:Nx,1:Nzloc)
  
                             
return
end subroutine

!***********************************************************************
!***********************************************************************
subroutine transform_rows_batched(output,ucomp)
!....Two-Dimensional Discrete Fourier Transforms in x and z
!....Parallel (MPI) FFTW Complex-to-Complex batched in y 
implicit none
complex(kind=prec)::output(0:Ny,0:Nx-1,0:Nzloc-1)
real(kind=prec)::ucomp(0:Ny,0:Nx-1,0:Nzloc-1)


 datacF_FWD=cpx0
 datacF_FWD=DCMPLX(ucomp)

!....compute transform (as many times as desired)
 call fftw_mpi_execute_dft(PLAN_FFT2D_FWD, datacF_FWD, datacF_FWD)
  
 output=datacF_FWD


 output=output*corrFFT2D

return
end subroutine
!***********************************************************************
subroutine transformInv_rows_batched(output,uhatcomp)
!....Inverse Two-Dimensional Discrete Fourier Transforms in x and z
!....Parallel (MPI) FFTW Complex-to-Complex batched in y 
implicit none
complex(kind=prec)::uhatcomp(0:Ny,0:Nx-1,0:Nzloc-1)
real(kind=prec)::output(0:Ny,0:Nx-1,0:Nzloc-1)




 datacF_BWD=uhatcomp

!....compute transform (as many times as desired)
 call fftw_mpi_execute_dft(PLAN_FFT2D_BWD, datacF_BWD, datacF_BWD)
  
  
 output=DREAL(datacF_BWD)
  

return
end subroutine

!***********************************************************************

end module transformspost
!***********************************************************************
!***********************************************************************
!***********************************************************************

!OLD SUBROUTINES

! module transformspost
! use, intrinsic :: iso_c_binding
! use mpi
! use parameterspost
! !use flowvariables
! include 'fftw3-mpi.f03'



! contains

! !***********************************************************************
! !***********************************************************************
! subroutine transform_rows(uhatcomp,ucomp)
! implicit none
! complex(kind=prec)::uhatcomp(0:Ny,0:Nx-1,0:Nzloc-1)
! real(kind=prec)::ucomp(0:Ny,0:Nx-1,0:Nzloc-1)
! real(kind=prec):: corr
! complex(kind=prec), allocatable::uaus(:,:,:)
! type(C_PTR) :: plan, cdata
! complex(C_DOUBLE_COMPLEX),pointer::data1(:,:)
! integer(C_INTPTR_T)::alloc_local, local_M, local_j_offset
! integer::ii

! !....Correction
! corr=1.0d0/(Nx*Nz)


! allocate(uaus(0:Ny,0:Nx-1,0:Nzloc-1))

! uhatcomp=DCMPLX(0.0d0,0.0d0)
! uaus=DCMPLX(ucomp) 

! !....get local data size and allocate (note dimension reversal)
!   alloc_local = fftw_mpi_local_size_2d(Nz, Nx, MPI_COMM_WORLD, &
!                                        local_M, local_j_offset)
!   cdata = fftw_alloc_complex(alloc_local)
!   call c_f_pointer(cdata, data1, [Nx,local_M])

! !....create MPI plan for in-place forward DFT (note dimension reversal)
!   plan = fftw_mpi_plan_dft_2d(Nz, Nx, data1, data1, MPI_COMM_WORLD, &
!                               FFTW_FORWARD, FFTW_MEASURE)

! do ii=0,Ny
! !....initialize data 
! data1(1:Nx,1:Nzloc)=uaus(ii,0:Nx-1,0:Nzloc-1)

! !....compute transform (as many times as desired)
!   call fftw_mpi_execute_dft(plan, data1, data1)
!   uhatcomp(ii,0:Nx-1,0:Nzloc-1)=data1(1:Nx,1:Nzloc)
  
! end do


!   call fftw_destroy_plan(plan)
!   call fftw_free(cdata)




! uhatcomp=uhatcomp*corr

! deallocate(uaus)
! return
! end subroutine
! !***********************************************************************
! !***********************************************************************
! subroutine transforminv_rows(ucomp,uhatcomp)
! implicit none
! complex(kind=prec)::uhatcomp(0:Ny,0:Nx-1,0:Nzloc-1)
! real(kind=prec)::ucomp(0:Ny,0:Nx-1,0:Nzloc-1)
! real(kind=prec):: corr
! complex(kind=prec), allocatable::uaus(:,:,:)
! type(C_PTR) :: plan, cdata
! complex(C_DOUBLE_COMPLEX),pointer::data1(:,:)
! integer(C_INTPTR_T)::alloc_local, local_M, local_j_offset
! integer::ii

! ucomp=0.0d0

! allocate(uaus(0:Ny,0:Nx-1,0:Nzloc-1))


! !....get local data size and allocate (note dimension reversal)
!   alloc_local = fftw_mpi_local_size_2d(Nz, Nx, MPI_COMM_WORLD, &
!                                        local_M, local_j_offset)
!   cdata = fftw_alloc_complex(alloc_local)
!   call c_f_pointer(cdata, data1, [Nx,local_M])

! !....create MPI plan for in-place backward DFT (note dimension reversal)
!   plan = fftw_mpi_plan_dft_2d(Nz, Nx, data1, data1, MPI_COMM_WORLD, &
!                               FFTW_BACKWARD, FFTW_MEASURE)
! do ii=0,Ny
! !....initialize data 
! data1(1:Nx,1:Nzloc)=uhatcomp(ii,0:Nx-1,0:Nzloc-1)

! !....compute transform (as many times as desired)
!   call fftw_mpi_execute_dft(plan, data1, data1)
  
!   uaus(ii,0:Nx-1,0:Nzloc-1)=data1(1:Nx,1:Nzloc)
! end do


!   call fftw_destroy_plan(plan)
!   call fftw_free(cdata)




! ucomp=DREAL(uaus)


! deallocate(uaus)
! return
! end subroutine
! !***********************************************************************
! !***********************************************************************
! !***********************************************************************
! !***********************************************************************
! subroutine ChebTransform1Dy(udmytr,udmy)
! implicit none
! real(kind=prec)::udmytr(0:Ny)
! real(kind=prec)::udmy(0:Ny)
! real(kind=prec), allocatable,dimension(:)::temp1,temp2
! integer::ll

! allocate(temp1(0:2*Ny-1),temp2(0:2*Ny-1))

! temp1=0.0d0
! temp2=0.0d0
! temp1(0:Ny)=udmy(0:Ny)
! do ll=1,Ny-1
! temp1(ll+Ny)=udmy(Ny-ll)
! end do

! call FFT1D_inv(temp2(0:2*Ny-1),temp1(0:2*Ny-1),(2*Ny-1))


! udmytr(0)=temp2(0)
! udmytr(Ny)=temp2(Ny)
! udmytr(1:Ny-1)=2.0d0*temp2(1:Ny-1)



! deallocate(temp1,temp2)

! return
! end subroutine
! !***********************************************************************
! !***********************************************************************
! subroutine ChebTransform1Dy_inv(udmy,udmytr)
! implicit none
! real(kind=prec)::udmy(0:Ny)
! real(kind=prec)::udmytr(0:Ny)
! real(kind=prec), allocatable,dimension(:)::temp1,temp2
! integer::ll

! allocate(temp1(0:2*Ny-1),temp2(0:2*Ny-1))

! temp1=0.0d0
! temp2=0.0d0
! temp1(0)=udmytr(0)
! temp1(1:Ny)=udmytr(1:Ny)/2.0d0
! do ll=1,Ny-1
! temp1(ll+Ny)=udmytr(Ny-ll)/2.0d0
! end do

! call FFT1D(temp2(0:2*Ny-1),temp1(0:2*Ny-1),(2*Ny-1))


! udmy(0:Ny)=temp2(0:Ny)


! deallocate(temp1,temp2)

! return
! end subroutine
! !***********************************************************************
! !***********************************************************************
! subroutine FFT1D(ucomptr,ucomp,nny)
! implicit none
! integer(kind=prec)::nny
! real(kind=prec)::ucomptr(0:nny)
! real(kind=prec)::ucomp(0:nny)
! complex(kind=prec), allocatable::uaus(:),uaus2(:)
! integer(kind=8) :: plan
! integer::FFTW_DIR

! allocate(uaus(0:nny))
! allocate(uaus2(0:nny))

! uaus=DCMPLX(ucomp)

! FFTW_DIR=-1 !direction of fft -1 for forward, 1 for backward (inverse transform)

! call dfftw_plan_dft_1d(plan,nny+1,uaus(0:nny),uaus2(0:nny),FFTW_DIR,FFTW_ESTIMATE)
! call dfftw_execute_dft(plan,uaus(0:nny),uaus2(0:nny))
! call dfftw_destroy_plan(plan)


! ucomptr=DREAL(uaus2)!/(nny+1)

! deallocate(uaus)
! deallocate(uaus2)
! return
! end subroutine

! !***********************************************************************
! !***********************************************************************
! subroutine FFT1D_inv(ucomp,ucomptr,nny)
! implicit none
! integer(kind=prec)::nny
! real(kind=prec)::ucomp(0:nny)
! real(kind=prec)::ucomptr(0:nny)
! complex(kind=prec), allocatable::uaus(:),uaus2(:)
! integer(kind=8) :: plan
! integer::FFTW_DIR

! allocate(uaus(0:nny))
! allocate(uaus2(0:nny))

! uaus=DCMPLX(ucomptr)

! FFTW_DIR=1 !direction of fft -1 for forward, 1 for backward (inverse transform)

! call dfftw_plan_dft_1d(plan,nny+1,uaus(0:nny),uaus2(0:nny),FFTW_DIR,FFTW_ESTIMATE)
! call dfftw_execute_dft(plan,uaus(0:nny),uaus2(0:nny))
! call dfftw_destroy_plan(plan)


! ucomp=DREAL(uaus2)/(nny+1)

! deallocate(uaus)
! deallocate(uaus2)
! return
! end subroutine
! !***********************************************************************
! !***********************************************************************
! subroutine transforminv_alongx(ucomp,uhatcomp)
! implicit none
! complex(kind=prec)::uhatcomp(0:Ny,0:Nx-1,0:Nzloc-1)
! complex(kind=prec)::ucomp(0:Ny,0:Nx-1,0:Nzloc-1)
! complex(kind=prec), allocatable::uaus(:),uaus2(:)
! integer(kind=8) :: plan
! integer::FFTW_DIR
! integer::ii,jj


! allocate(uaus(0:Nx-1))
! allocate(uaus2(0:Nx-1))
! do ii=0,Nzloc-1
! do jj=0,Ny

! uaus(0:Nx-1)=uhatcomp(jj,0:Nx-1,ii)
! uaus2=cpx0

! FFTW_DIR=1 !direction of fft -1 for forward, 1 for backward (inverse transform)

! call dfftw_plan_dft_1d(plan,Nx,uaus(0:Nx-1),uaus2(0:Nx-1),FFTW_DIR,FFTW_ESTIMATE)
! call dfftw_execute_dft(plan,uaus(0:Nx-1),uaus2(0:Nx-1))
! call dfftw_destroy_plan(plan)

! uaus2=uaus2/Nx

! ucomp(jj,0:Nx-1,ii)=uaus2(0:Nx-1)

! end do
! end do




! deallocate(uaus,uaus2)
! return
! end subroutine
! !***********************************************************************
! !***********************************************************************


! end module transformspost
