program convflow
use, intrinsic :: iso_c_binding
use mpi
use parameters
use physics
use flowvariables
use in_out
use setup
use transforms

implicit none

integer:: ii,jj,ll

!.....Variables for parallelization:
integer status(MPI_Status_size)

!.....Parallel initialization:
call MPI_Init(ierror)
call MPI_Comm_size(MPI_COMM_WORLD,noprocs,ierror)!noprocs=number of processors
call MPI_Comm_rank(MPI_COMM_WORLD,nid,ierror)    !nid=rank of each processor
call fftw_mpi_init()

!.....End of initialization...
print *,'nid',nid
call MPI_BARRIER(MPI_Comm_World,ierror)    

!....Allocate namefiles
allocate(nomefile(1:3))
allocate(nomefile1(1:3))
allocate(resfile(1:3))
! allocate(gacc(1:2))


!.....Parameter are read by processor 0 from file in_dns_tc.txt
!and spread to all processors
call read_inputs()

call MPI_BARRIER(MPI_Comm_World,ierror) 

!...Allocate Variables
call allocate_variabs()

call MPI_BARRIER(MPI_Comm_World,ierror) 

!...Initialize calculation
call InitializeFlow()

! !...Initialize calculation
! call initialize()

call MPI_BARRIER(MPI_Comm_World,ierror) 

call transforms_init()

call MPI_BARRIER(MPI_Comm_World,ierror) 

!....Zeroing the mean fields
uhatMT=cpx0
vhatMT=cpx0
whatMT=cpx0
! ThetahatMT=cpx0
! ThetaSpzMT=0.0d0

!Zeroing variables for statistical stationary
Ubulk=0.0d0
uumax=0.0d0



if(nid.eq.0)then
WRITE(*,*)'Evaluating Mean Profiles.'
WRITE(*,*)'Processing step: ',IndSave(1)
end if

!....Generate names of the files to be loaded
write(snapshot,'(i8)') (IndSave(1))
nomefile(1)='u_'//trim(adjustl(snapshot))//'.bin'
nomefile(2)='v_'//trim(adjustl(snapshot))//'.bin'
nomefile(3)='w_'//trim(adjustl(snapshot))//'.bin'



!....Read solutions at the time-step 1
call LoadScalarField(u0(0:Ny,0:Nx-1,0:Nzloc-1),nomefile(1))
call LoadScalarField(v0(0:Ny,0:Nx-1,0:Nzloc-1),nomefile(2))
call LoadScalarField(w0(0:Ny,0:Nx-1,0:Nzloc-1),nomefile(3))

!....Transorm to Fourier Space
call transform_rows_batched(uhat,u0(0:Ny,0:Nx-1,0:Nzloc-1))
call transform_rows_batched(vhat,v0(0:Ny,0:Nx-1,0:Nzloc-1))
call transform_rows_batched(what,w0(0:Ny,0:Nx-1,0:Nzloc-1))




call EvalUbulk(Ubulk(1),uhat(0:Ny,0,0))

uhatMT=uhat/2.0d0
vhatMT=vhat/2.0d0
whatMT=what/2.0d0

call MPI_BARRIER(MPI_Comm_World,ierror) 


!....Loop to calculate means
do Iext=2,Nsavings-1

if(nid.eq.0)then
WRITE(*,*)'Processing step: ',IndSave(Iext)
end if

!....Generate names of the files to be loaded
write(snapshot,'(i8)') (IndSave(Iext))
nomefile(1)='u_'//trim(adjustl(snapshot))//'.bin'
nomefile(2)='v_'//trim(adjustl(snapshot))//'.bin'
nomefile(3)='w_'//trim(adjustl(snapshot))//'.bin'




!....Read solutions at the time-step n
call LoadScalarField(u0(0:Ny,0:Nx-1,0:Nzloc-1),nomefile(1))
call LoadScalarField(v0(0:Ny,0:Nx-1,0:Nzloc-1),nomefile(2))
call LoadScalarField(w0(0:Ny,0:Nx-1,0:Nzloc-1),nomefile(3))



!....Transorm to Fourier Space
call transform_rows_batched(uhat,u0(0:Ny,0:Nx-1,0:Nzloc-1))
call transform_rows_batched(vhat,v0(0:Ny,0:Nx-1,0:Nzloc-1))
call transform_rows_batched(what,w0(0:Ny,0:Nx-1,0:Nzloc-1))


uhatMT=uhatMT+uhat
vhatMT=vhatMT+vhat
whatMT=whatMT+what

call EvalUbulk(Ubulk(Iext),uhat(0:Ny,0,0))

call MPI_BARRIER(MPI_Comm_World,ierror) 


end do
!End of loop for central savings

call MPI_BARRIER(MPI_Comm_World,ierror) 

if(nid.eq.0)then
WRITE(*,*)'Processing step: ',IndSave(Nsavings)
end if

!....Generate names of the files to be loaded
write(snapshot,'(i8)') (IndSave(Nsavings))
nomefile(1)='u_'//trim(adjustl(snapshot))//'.bin'
nomefile(2)='v_'//trim(adjustl(snapshot))//'.bin'
nomefile(3)='w_'//trim(adjustl(snapshot))//'.bin'

!....Read solutions at the time-step n
call LoadScalarField(u0(0:Ny,0:Nx-1,0:Nzloc-1),nomefile(1))
call LoadScalarField(v0(0:Ny,0:Nx-1,0:Nzloc-1),nomefile(2))
call LoadScalarField(w0(0:Ny,0:Nx-1,0:Nzloc-1),nomefile(3))


!....Transorm to Fourier Space
call transform_rows_batched(uhat,u0(0:Ny,0:Nx-1,0:Nzloc-1))
call transform_rows_batched(vhat,v0(0:Ny,0:Nx-1,0:Nzloc-1))
call transform_rows_batched(what,w0(0:Ny,0:Nx-1,0:Nzloc-1))


uhatMT=uhatMT+uhat/2.0d0
vhatMT=vhatMT+vhat/2.0d0
whatMT=whatMT+what/2.0d0


uhatMT=uhatMT/Ninterv
vhatMT=vhatMT/Ninterv
whatMT=whatMT/Ninterv


call EvalUbulk(Ubulk(Nsavings),uhat(0:Ny,0,0))

call MPI_BARRIER(MPI_Comm_World,ierror) 
!.....End of means evaluation


!....Saving Mean Fields
if(nid.eq.0)then
WRITE(*,*)'**********************************************************'
WRITE(*,*)'**********************************************************' 
end if


call transformInv_rows_batched(uMT(0:Ny,0:Nx-1,0:Nzloc-1),uhatMT)
call transformInv_rows_batched(vMT(0:Ny,0:Nx-1,0:Nzloc-1),vhatMT)
call transformInv_rows_batched(wMT(0:Ny,0:Nx-1,0:Nzloc-1),whatMT)


call MPI_BARRIER(MPI_Comm_World,ierror) 


 
!....Reynolds Stresses
ReStress=cpx0



if(nid.eq.0)then
WRITE(*,*)'**********************************************************'
WRITE(*,*)'**********************************************************' 
WRITE(*,*)'Evaluating Turbulence Statistics and Turbulent Kinetic Energy Budget' 
WRITE(*,*)
WRITE(*,*)'Processing step: ',IndSave(1)
end if
call MPI_BARRIER(MPI_Comm_World,ierror)  


write(snapshot,'(i8)') (IndSave(1))
nomefile(1)='u_'//trim(adjustl(snapshot))//'.bin'
nomefile(2)='v_'//trim(adjustl(snapshot))//'.bin'
nomefile(3)='w_'//trim(adjustl(snapshot))//'.bin'

call MPI_BARRIER(MPI_Comm_World,ierror) 


!....Read solutions at the time-step 1
call LoadScalarField(u0(0:Ny,0:Nx-1,0:Nzloc-1),nomefile(1))
call LoadScalarField(v0(0:Ny,0:Nx-1,0:Nzloc-1),nomefile(2))
call LoadScalarField(w0(0:Ny,0:Nx-1,0:Nzloc-1),nomefile(3))

!...Removing the time-averaged fields in Physical space
u0(0:Ny,0:Nx-1,0:Nzloc-1)=u0(0:Ny,0:Nx-1,0:Nzloc-1)-uMT(0:Ny,0:Nx-1,0:Nzloc-1)
v0(0:Ny,0:Nx-1,0:Nzloc-1)=v0(0:Ny,0:Nx-1,0:Nzloc-1)-vMT(0:Ny,0:Nx-1,0:Nzloc-1)
w0(0:Ny,0:Nx-1,0:Nzloc-1)=w0(0:Ny,0:Nx-1,0:Nzloc-1)-wMT(0:Ny,0:Nx-1,0:Nzloc-1)

!....Transorm to Fourier Space
call transform_rows_batched(uhat,u0(0:Ny,0:Nx-1,0:Nzloc-1))
call transform_rows_batched(vhat,v0(0:Ny,0:Nx-1,0:Nzloc-1))
call transform_rows_batched(what,w0(0:Ny,0:Nx-1,0:Nzloc-1))

!....Removing also the mean field in space to consider the fluctuations only
!.....Mean fields in x and z
if(nid.eq.0)then
uhat(0:Ny,0,0)=cpx0
vhat(0:Ny,0,0)=cpx0
what(0:Ny,0,0)=cpx0
end if

!....Transform back to Physical Space

call transformInv_rows_batched(u0(0:Ny,0:Nx-1,0:Nzloc-1),uhat)
call transformInv_rows_batched(v0(0:Ny,0:Nx-1,0:Nzloc-1),vhat)
call transformInv_rows_batched(w0(0:Ny,0:Nx-1,0:Nzloc-1),what)

call MPI_BARRIER(MPI_Comm_World,ierror) 

!....Evaluating Statistics (Re Stress)
call ConvectiveStress(ConvStr,u0(0:Ny,0:Nx-1,0:Nzloc-1),v0(0:Ny,0:Nx-1,0:Nzloc-1)&
,w0(0:Ny,0:Nx-1,0:Nzloc-1))


call EvaluateMaxuu(uumax(1),ConvStr(0:Ny,0:Nx-1,0:Nzloc-1,1))


ReStress=ConvStr/2.0d0


!....Loop to calculate means
do Iext=2,Nsavings-1

if(nid.eq.0)then
WRITE(*,*)'Processing step: ',IndSave(Iext)
end if

!....Generate names of the files to be loaded
write(snapshot,'(i8)') (IndSave(Iext))
nomefile(1)='u_'//trim(adjustl(snapshot))//'.bin'
nomefile(2)='v_'//trim(adjustl(snapshot))//'.bin'
nomefile(3)='w_'//trim(adjustl(snapshot))//'.bin'


!....Read solutions at the time-step n
call LoadScalarField(u0(0:Ny,0:Nx-1,0:Nzloc-1),nomefile(1))
call LoadScalarField(v0(0:Ny,0:Nx-1,0:Nzloc-1),nomefile(2))
call LoadScalarField(w0(0:Ny,0:Nx-1,0:Nzloc-1),nomefile(3))




!...Removing the time-averaged fields in Physical space
u0(0:Ny,0:Nx-1,0:Nzloc-1)=u0(0:Ny,0:Nx-1,0:Nzloc-1)-uMT(0:Ny,0:Nx-1,0:Nzloc-1)
v0(0:Ny,0:Nx-1,0:Nzloc-1)=v0(0:Ny,0:Nx-1,0:Nzloc-1)-vMT(0:Ny,0:Nx-1,0:Nzloc-1)
w0(0:Ny,0:Nx-1,0:Nzloc-1)=w0(0:Ny,0:Nx-1,0:Nzloc-1)-wMT(0:Ny,0:Nx-1,0:Nzloc-1)



!....Transorm to Fourier Space
call transform_rows_batched(uhat,u0(0:Ny,0:Nx-1,0:Nzloc-1))
call transform_rows_batched(vhat,v0(0:Ny,0:Nx-1,0:Nzloc-1))
call transform_rows_batched(what,w0(0:Ny,0:Nx-1,0:Nzloc-1))



!....Removing also the mean field in space to consider the fluctuations only
!.....Mean fields in x and z
if(nid.eq.0)then
uhat(0:Ny,0,0)=cpx0
vhat(0:Ny,0,0)=cpx0
what(0:Ny,0,0)=cpx0
end if

!....Transform back to Physical Space

call transformInv_rows_batched(u0(0:Ny,0:Nx-1,0:Nzloc-1),uhat)
call transformInv_rows_batched(v0(0:Ny,0:Nx-1,0:Nzloc-1),vhat)
call transformInv_rows_batched(w0(0:Ny,0:Nx-1,0:Nzloc-1),what)


call MPI_BARRIER(MPI_Comm_World,ierror) 

!....Evaluating Statistics (Re Stress)
call ConvectiveStress(ConvStr,u0(0:Ny,0:Nx-1,0:Nzloc-1),v0(0:Ny,0:Nx-1,0:Nzloc-1)&
,w0(0:Ny,0:Nx-1,0:Nzloc-1))



call EvaluateMaxuu(uumax(Iext),ConvStr(0:Ny,0:Nx-1,0:Nzloc-1,1))


ReStress=ReStress+ConvStr

call MPI_BARRIER(MPI_Comm_World,ierror) 


end do
!End of loop for central savings


call MPI_BARRIER(MPI_Comm_World,ierror) 

if(nid.eq.0)then
WRITE(*,*)'Processing step: ',IndSave(Nsavings)
end if

!....Generate names of the files to be loaded
write(snapshot,'(i8)') (IndSave(Nsavings))
nomefile(1)='u_'//trim(adjustl(snapshot))//'.bin'
nomefile(2)='v_'//trim(adjustl(snapshot))//'.bin'
nomefile(3)='w_'//trim(adjustl(snapshot))//'.bin'



!....Read solutions at the time-step n
call LoadScalarField(u0(0:Ny,0:Nx-1,0:Nzloc-1),nomefile(1))
call LoadScalarField(v0(0:Ny,0:Nx-1,0:Nzloc-1),nomefile(2))
call LoadScalarField(w0(0:Ny,0:Nx-1,0:Nzloc-1),nomefile(3))


!...Removing the time-averaged fields in Physical space
u0(0:Ny,0:Nx-1,0:Nzloc-1)=u0(0:Ny,0:Nx-1,0:Nzloc-1)-uMT(0:Ny,0:Nx-1,0:Nzloc-1)
v0(0:Ny,0:Nx-1,0:Nzloc-1)=v0(0:Ny,0:Nx-1,0:Nzloc-1)-vMT(0:Ny,0:Nx-1,0:Nzloc-1)
w0(0:Ny,0:Nx-1,0:Nzloc-1)=w0(0:Ny,0:Nx-1,0:Nzloc-1)-wMT(0:Ny,0:Nx-1,0:Nzloc-1)



!....Transorm to Fourier Space
call transform_rows_batched(uhat,u0(0:Ny,0:Nx-1,0:Nzloc-1))
call transform_rows_batched(vhat,v0(0:Ny,0:Nx-1,0:Nzloc-1))
call transform_rows_batched(what,w0(0:Ny,0:Nx-1,0:Nzloc-1))



!....Removing also the mean field in space to consider the fluctuations only
!.....Mean fields in x and z


if(nid.eq.0)then
uhat(0:Ny,0,0)=cpx0
vhat(0:Ny,0,0)=cpx0
what(0:Ny,0,0)=cpx0
end if

!....Transform back to Physical Space

call transformInv_rows_batched(u0(0:Ny,0:Nx-1,0:Nzloc-1),uhat)
call transformInv_rows_batched(v0(0:Ny,0:Nx-1,0:Nzloc-1),vhat)
call transformInv_rows_batched(w0(0:Ny,0:Nx-1,0:Nzloc-1),what)



call MPI_BARRIER(MPI_Comm_World,ierror) 
!....Evaluating Statistics (Re Stress)
call ConvectiveStress(ConvStr,u0(0:Ny,0:Nx-1,0:Nzloc-1),v0(0:Ny,0:Nx-1,0:Nzloc-1)&
,w0(0:Ny,0:Nx-1,0:Nzloc-1))



call EvaluateMaxuu(uumax(Nsavings),ConvStr(0:Ny,0:Nx-1,0:Nzloc-1,1))



ReStress=ReStress+ConvStr/2.0d0

ReStress=ReStress/Ninterv

call MPI_BARRIER(MPI_Comm_World,ierror) 


if(nid.eq.0)then
WRITE(*,*)'Turbulence Statistics  Evaluated.'
WRITE(*,*)'**********************************************************'
WRITE(*,*)'**********************************************************'  
end if
call MPI_BARRIER(MPI_Comm_World,ierror)  

if(nid.eq.0)then
WRITE(*,*)
WRITE(*,*)'**********************************************************'
WRITE(*,*)'**********************************************************' 
WRITE(*,*)'Saving mean profiles'
 
end if


resfile(1)='StationaryVarbs.txt'
FMT1="(2X,E12.5,2X,E12.5,2X,E12.5,)"

if(nid.eq.0)then
open(2,file=resfile(1))
write(2,*) '%       step         Ubulk       max(uu) '
do ii=1,Nsavings
write(2,FMT1) dfloat(ii),Ubulk(ii),uumax(ii)
end do

close(2)

end if

call MPI_BARRIER(MPI_Comm_World,ierror) 






! resfile(1)='PScalstats.txt'
! FMT1="(2X,E12.5,2X,E12.5,2X,E12.5,2X,E12.5,2X,E12.5,2X,E12.5,2X,E12.5,2X,E12.5)"

! if(nid.eq.0)then
! open(4,file=resfile(1))

! write(4,*) '%       y            y+         <Theta>            d<Theta>/dy         <theta theta>         &
! <u theta>                <v theta>               <w theta> '!          k_theta'
! do ii=0,Ny
! write(4,FMT1) y(Ny-ii),yp(ii),Theta(Ny-ii,0,0),dThetady(Ny-ii) ,&
! TheStressPh(Ny-ii,0,0,4),TheStressPh(Ny-ii,0,0,1),&
! TheStressPh(Ny-ii,0,0,2),TheStressPh(Ny-ii,0,0,3)

! end do


! close(4)

! end if

call transforms_finalization()


call MPI_BARRIER(MPI_Comm_World,ierror) 


!....Deallocate Variables  
call deallocate_variabs()

! deallocate(gacc)
deallocate(nomefile,nomefile1)
deallocate(resfile)


!.....Closing
write(*,*) 'Stop!!! Processore:',nid
call MPI_Finalize(ierror)
      
stop

end program convflow

