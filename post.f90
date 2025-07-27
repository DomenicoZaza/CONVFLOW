program post
use, intrinsic :: iso_c_binding
use mpi
use parameterspost
use physicspost
use flowvariablespost
use in_outpost
use setuppost
use transformspost


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


! fileScal='PScal_'//trim(adjustl(snapshot))//'.bin'


!....Read solutions at the time-step 1
call LoadScalarField(u0(0:Ny,0:Nx-1,0:Nzloc-1),nomefile(1))
call LoadScalarField(v0(0:Ny,0:Nx-1,0:Nzloc-1),nomefile(2))
call LoadScalarField(w0(0:Ny,0:Nx-1,0:Nzloc-1),nomefile(3))
! call LoadScalarField(Theta(0:Ny,0:Nx-1,0:Nzloc-1),fileScal)

! ThetaSpz(0:Ny,0:Nx-1,0:Nzloc-1)=Theta(0:Ny,0:Nx-1,0:Nzloc-1)


!....Transorm to Fourier Space
call transform_rows_batched(uhat,u0(0:Ny,0:Nx-1,0:Nzloc-1))
call transform_rows_batched(vhat,v0(0:Ny,0:Nx-1,0:Nzloc-1))
call transform_rows_batched(what,w0(0:Ny,0:Nx-1,0:Nzloc-1))

! call transform_rows(Thetahat,Theta(0:Ny,0:Nx-1,0:Nzloc-1))

! call transforminv_alongx(ThetahatSpz,Thetahat)

! ThetaSpz(0:Ny,0:Nx-1,0:Nzloc-1)=ABS(ThetahatSpz(0:Ny,0:Nx-1,0:Nzloc-1))

! call MeanX(ThetaSpz2(0:Ny,0:Nzloc-1),ThetaSpz(0:Ny,0:Nx-1,0:Nzloc-1))


call EvalUbulk(Ubulk(1),uhat(0:Ny,0,0))

uhatMT=uhat/2.0d0
vhatMT=vhat/2.0d0
whatMT=what/2.0d0
! ThetahatMT=Thetahat/2.0d0
! ThetaSpzMT=ThetaSpz2/2.0d0

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

! fileScal='PScal_'//trim(adjustl(snapshot))//'.bin'


!....Read solutions at the time-step n
call LoadScalarField(u0(0:Ny,0:Nx-1,0:Nzloc-1),nomefile(1))
call LoadScalarField(v0(0:Ny,0:Nx-1,0:Nzloc-1),nomefile(2))
call LoadScalarField(w0(0:Ny,0:Nx-1,0:Nzloc-1),nomefile(3))

! call LoadScalarField(Theta(0:Ny,0:Nx-1,0:Nzloc-1),fileScal)


!....Transorm to Fourier Space
call transform_rows_batched(uhat,u0(0:Ny,0:Nx-1,0:Nzloc-1))
call transform_rows_batched(vhat,v0(0:Ny,0:Nx-1,0:Nzloc-1))
call transform_rows_batched(what,w0(0:Ny,0:Nx-1,0:Nzloc-1))

! call transform_rows(Thetahat(0:Ny,0:Nx-1,0:Nzloc-1),Theta(0:Ny,0:Nx-1,0:Nzloc-1))

! call transforminv_alongx(ThetahatSpz,Thetahat)

! ThetaSpz(0:Ny,0:Nx-1,0:Nzloc-1)=ABS(ThetahatSpz(0:Ny,0:Nx-1,0:Nzloc-1))

! call MeanX(ThetaSpz2(0:Ny,0:Nzloc-1),ThetaSpz(0:Ny,0:Nx-1,0:Nzloc-1))


uhatMT=uhatMT+uhat
vhatMT=vhatMT+vhat
whatMT=whatMT+what
! ThetahatMT=ThetahatMT+Thetahat
! ThetaSpzMT=ThetaSpzMT+ThetaSpz2

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


! fileScal='PScal_'//trim(adjustl(snapshot))//'.bin'


!....Read solutions at the time-step n
call LoadScalarField(u0(0:Ny,0:Nx-1,0:Nzloc-1),nomefile(1))
call LoadScalarField(v0(0:Ny,0:Nx-1,0:Nzloc-1),nomefile(2))
call LoadScalarField(w0(0:Ny,0:Nx-1,0:Nzloc-1),nomefile(3))

! call LoadScalarField(Theta(0:Ny,0:Nx-1,0:Nzloc-1),fileScal)


!....Transorm to Fourier Space
call transform_rows_batched(uhat,u0(0:Ny,0:Nx-1,0:Nzloc-1))
call transform_rows_batched(vhat,v0(0:Ny,0:Nx-1,0:Nzloc-1))
call transform_rows_batched(what,w0(0:Ny,0:Nx-1,0:Nzloc-1))

! call transform_rows(Thetahat(0:Ny,0:Nx-1,0:Nzloc-1),Theta(0:Ny,0:Nx-1,0:Nzloc-1))

! call transforminv_alongx(ThetahatSpz,Thetahat)

! ThetaSpz(0:Ny,0:Nx-1,0:Nzloc-1)=ABS(ThetahatSpz(0:Ny,0:Nx-1,0:Nzloc-1))

! call MeanX(ThetaSpz2(0:Ny,0:Nzloc-1),ThetaSpz(0:Ny,0:Nx-1,0:Nzloc-1))


uhatMT=uhatMT+uhat/2.0d0
vhatMT=vhatMT+vhat/2.0d0
whatMT=whatMT+what/2.0d0

! ThetahatMT=ThetahatMT+Thetahat/2.0d0
! ThetaSpzMT=ThetaSpzMT+ThetaSpz2/2.0d0



uhatMT=uhatMT/Ninterv
vhatMT=vhatMT/Ninterv
whatMT=whatMT/Ninterv
! ThetahatMT=ThetahatMT/Ninterv
! ThetaSpzMT=ThetaSpzMT/Ninterv


call EvalUbulk(Ubulk(Nsavings),uhat(0:Ny,0,0))

call MPI_BARRIER(MPI_Comm_World,ierror) 
!.....End of means evaluation


!....Saving Mean Fields
if(nid.eq.0)then
WRITE(*,*)'**********************************************************'
WRITE(*,*)'**********************************************************' 
WRITE(*,*)'Saving the mean fields...'
WRITE(*,*)
end if

resfile(1)='u_mean.bin'
resfile(2)='v_mean.bin'
resfile(3)='w_mean.bin'
! fileScal1='PScal_mean.bin'


call transformInv_rows_batched(uMT(0:Ny,0:Nx-1,0:Nzloc-1),uhatMT)
call transformInv_rows_batched(vMT(0:Ny,0:Nx-1,0:Nzloc-1),vhatMT)
call transformInv_rows_batched(wMT(0:Ny,0:Nx-1,0:Nzloc-1),whatMT)

! call transforminv_rows(ThetaMT(0:Ny,0:Nx-1,0:Nzloc-1),ThetahatMT(0:Ny,0:Nx-1,0:Nzloc-1))



call MPI_BARRIER(MPI_Comm_World,ierror) 

call SaveScalarField(uMT(0:Ny,0:Nx-1,0:Nzloc-1),resfile(1))
call MPI_BARRIER(MPI_Comm_World,ierror) 
call SaveScalarField(vMT(0:Ny,0:Nx-1,0:Nzloc-1),resfile(2))
call MPI_BARRIER(MPI_Comm_World,ierror) 
call SaveScalarField(wMT(0:Ny,0:Nx-1,0:Nzloc-1),resfile(3))
call MPI_BARRIER(MPI_Comm_World,ierror)  
! call SaveScalarField(ThetaMT(0:Ny,0:Nx-1,0:Nzloc-1),fileScal1)
! call MPI_BARRIER(MPI_Comm_World,ierror) 
! resfile(1)='SpectrumZ.bin'
! call SaveSpectrumZ(ThetaSpzMT(0:Ny,0:Nzloc-1),resfile(1))
! call MPI_BARRIER(MPI_Comm_World,ierror) 



if(nid.eq.0)then
WRITE(*,*)
WRITE(*,*)'Mean Fields Saved.'
end if

call MPI_BARRIER(MPI_Comm_World,ierror)  
!....End of saving mean fields.

 
!....Reynolds Stresses
ReStress=cpx0
! TheStress=cpx0
DissipMT=0.0d0
TurbTrMT=0.0d0
PressTrMT=0.0d0


!...Da modificare 
DiffPress=0.0d0



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


! fileScal='PScal_'//trim(adjustl(snapshot))//'.bin'

call MPI_BARRIER(MPI_Comm_World,ierror) 


!....Read solutions at the time-step 1
call LoadScalarField(u0(0:Ny,0:Nx-1,0:Nzloc-1),nomefile(1))
call LoadScalarField(v0(0:Ny,0:Nx-1,0:Nzloc-1),nomefile(2))
call LoadScalarField(w0(0:Ny,0:Nx-1,0:Nzloc-1),nomefile(3))
! call LoadScalarField(Theta(0:Ny,0:Nx-1,0:Nzloc-1),fileScal)


!...Removing the time-averaged fields in Physical space
u0(0:Ny,0:Nx-1,0:Nzloc-1)=u0(0:Ny,0:Nx-1,0:Nzloc-1)-uMT(0:Ny,0:Nx-1,0:Nzloc-1)
v0(0:Ny,0:Nx-1,0:Nzloc-1)=v0(0:Ny,0:Nx-1,0:Nzloc-1)-vMT(0:Ny,0:Nx-1,0:Nzloc-1)
w0(0:Ny,0:Nx-1,0:Nzloc-1)=w0(0:Ny,0:Nx-1,0:Nzloc-1)-wMT(0:Ny,0:Nx-1,0:Nzloc-1)
! Theta(0:Ny,0:Nx-1,0:Nzloc-1)=Theta(0:Ny,0:Nx-1,0:Nzloc-1)-ThetaMT(0:Ny,0:Nx-1,0:Nzloc-1)



!....Transorm to Fourier Space
call transform_rows_batched(uhat,u0(0:Ny,0:Nx-1,0:Nzloc-1))
call transform_rows_batched(vhat,v0(0:Ny,0:Nx-1,0:Nzloc-1))
call transform_rows_batched(what,w0(0:Ny,0:Nx-1,0:Nzloc-1))
! call transform_rows(Thetahat(0:Ny,0:Nx-1,0:Nzloc-1),Theta(0:Ny,0:Nx-1,0:Nzloc-1))


!....Removing also the mean field in space to consider the fluctuations only
!.....Mean fields in x and z


if(nid.eq.0)then
uhat(0:Ny,0,0)=cpx0
vhat(0:Ny,0,0)=cpx0
what(0:Ny,0,0)=cpx0
! Thetahat(0:Ny,0,0)=cpx0
end if

!....Transform back to Physical Space

call transformInv_rows_batched(u0(0:Ny,0:Nx-1,0:Nzloc-1),uhat)
call transformInv_rows_batched(v0(0:Ny,0:Nx-1,0:Nzloc-1),vhat)
call transformInv_rows_batched(w0(0:Ny,0:Nx-1,0:Nzloc-1),what)
! call transforminv_rows(Theta(0:Ny,0:Nx-1,0:Nzloc-1),Thetahat(0:Ny,0:Nx-1,0:Nzloc-1))


call MPI_BARRIER(MPI_Comm_World,ierror) 

!....Evaluating Statistics (Re Stress)
call ConvectiveStress(ConvStr,u0(0:Ny,0:Nx-1,0:Nzloc-1),v0(0:Ny,0:Nx-1,0:Nzloc-1)&
,w0(0:Ny,0:Nx-1,0:Nzloc-1))

! !....Evaluating Statistics (Re Stress and Passive scalar stress)
! call ConvectiveStressSc(ConvStr,ConvStrSc,u0(0:Ny,0:Nx-1,0:Nzloc-1),v0(0:Ny,0:Nx-1,0:Nzloc-1)&
! ,w0(0:Ny,0:Nx-1,0:Nzloc-1),Theta(0:Ny,0:Nx-1,0:Nzloc-1))


call EvaluateMaxuu(uumax(1),ConvStr(0:Ny,0:Nx-1,0:Nzloc-1,1))

!...Evaluate Turbulent Dissipation epsilon=1/Re_tau*<dui/dxj*dui/dxj>
call DissipationTerm(Dissip,uhat,vhat,what)
!....Evaluate Turbulent Transport
call TurbulentTransportTerm(TurbTr,u0(0:Ny,0:Nx-1,0:Nzloc-1),v0(0:Ny,0:Nx-1,0:Nzloc-1),&
w0(0:Ny,0:Nx-1,0:Nzloc-1))

!!....Evaluate the pressure
!call EvalPressure()

!!....Evaluate Pressure Transport Term
!call PressureTransportTerm(PressTr,uhat,vhat,what)



ReStress=ConvStr/2.0d0
! TheStress=ConvStrSc/2.0d0


!...Mean Dissipation
DissipMT=Dissip/2.0d0
!...Mean Turbulent Transport
TurbTrMT=TurbTr/2.0d0

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


! fileScal='PScal_'//trim(adjustl(snapshot))//'.bin'


!....Read solutions at the time-step n
call LoadScalarField(u0(0:Ny,0:Nx-1,0:Nzloc-1),nomefile(1))
call LoadScalarField(v0(0:Ny,0:Nx-1,0:Nzloc-1),nomefile(2))
call LoadScalarField(w0(0:Ny,0:Nx-1,0:Nzloc-1),nomefile(3))
! call LoadScalarField(Theta(0:Ny,0:Nx-1,0:Nzloc-1),fileScal)




!...Removing the time-averaged fields in Physical space
u0(0:Ny,0:Nx-1,0:Nzloc-1)=u0(0:Ny,0:Nx-1,0:Nzloc-1)-uMT(0:Ny,0:Nx-1,0:Nzloc-1)
v0(0:Ny,0:Nx-1,0:Nzloc-1)=v0(0:Ny,0:Nx-1,0:Nzloc-1)-vMT(0:Ny,0:Nx-1,0:Nzloc-1)
w0(0:Ny,0:Nx-1,0:Nzloc-1)=w0(0:Ny,0:Nx-1,0:Nzloc-1)-wMT(0:Ny,0:Nx-1,0:Nzloc-1)
! Theta(0:Ny,0:Nx-1,0:Nzloc-1)=Theta(0:Ny,0:Nx-1,0:Nzloc-1)-ThetaMT(0:Ny,0:Nx-1,0:Nzloc-1)



!....Transorm to Fourier Space
call transform_rows_batched(uhat,u0(0:Ny,0:Nx-1,0:Nzloc-1))
call transform_rows_batched(vhat,v0(0:Ny,0:Nx-1,0:Nzloc-1))
call transform_rows_batched(what,w0(0:Ny,0:Nx-1,0:Nzloc-1))
! call transform_rows(Thetahat(0:Ny,0:Nx-1,0:Nzloc-1),Theta(0:Ny,0:Nx-1,0:Nzloc-1))


!....Removing also the mean field in space to consider the fluctuations only
!.....Mean fields in x and z


if(nid.eq.0)then
uhat(0:Ny,0,0)=cpx0
vhat(0:Ny,0,0)=cpx0
what(0:Ny,0,0)=cpx0
! Thetahat(0:Ny,0,0)=cpx0


end if

!....Transform back to Physical Space

call transformInv_rows_batched(u0(0:Ny,0:Nx-1,0:Nzloc-1),uhat)
call transformInv_rows_batched(v0(0:Ny,0:Nx-1,0:Nzloc-1),vhat)
call transformInv_rows_batched(w0(0:Ny,0:Nx-1,0:Nzloc-1),what)
! call transforminv_rows(Theta(0:Ny,0:Nx-1,0:Nzloc-1),Thetahat(0:Ny,0:Nx-1,0:Nzloc-1))


call MPI_BARRIER(MPI_Comm_World,ierror) 

!....Evaluating Statistics (Re Stress)
call ConvectiveStress(ConvStr,u0(0:Ny,0:Nx-1,0:Nzloc-1),v0(0:Ny,0:Nx-1,0:Nzloc-1)&
,w0(0:Ny,0:Nx-1,0:Nzloc-1))

! !....Evaluating Statistics (Re Stress and Passive scalar stress)
! call ConvectiveStressSc(ConvStr,ConvStrSc,u0(0:Ny,0:Nx-1,0:Nzloc-1),&
! v0(0:Ny,0:Nx-1,0:Nzloc-1),w0(0:Ny,0:Nx-1,0:Nzloc-1),Theta(0:Ny,0:Nx-1,0:Nzloc-1))


call EvaluateMaxuu(uumax(Iext),ConvStr(0:Ny,0:Nx-1,0:Nzloc-1,1))

!...Evaluate Turbulent Dissipation epsilon=1/Re_tau*<dui/dxj*dui/dxj>
call DissipationTerm(Dissip,uhat,vhat,what)
!....Evaluate Turbulent Transport
call TurbulentTransportTerm(TurbTr,u0(0:Ny,0:Nx-1,0:Nzloc-1),v0(0:Ny,0:Nx-1,0:Nzloc-1)&
,w0(0:Ny,0:Nx-1,0:Nzloc-1))


ReStress=ReStress+ConvStr
! TheStress=TheStress+ConvStrSc

!...Mean Dissipation
DissipMT=DissipMT+Dissip
!...Mean Turbulent Transport
TurbTrMT=TurbTrMT+TurbTr

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

! fileScal='PScal_'//trim(adjustl(snapshot))//'.bin'


!....Read solutions at the time-step n
call LoadScalarField(u0(0:Ny,0:Nx-1,0:Nzloc-1),nomefile(1))
call LoadScalarField(v0(0:Ny,0:Nx-1,0:Nzloc-1),nomefile(2))
call LoadScalarField(w0(0:Ny,0:Nx-1,0:Nzloc-1),nomefile(3))
! call LoadScalarField(Theta(0:Ny,0:Nx-1,0:Nzloc-1),fileScal)



!...Removing the time-averaged fields in Physical space
u0(0:Ny,0:Nx-1,0:Nzloc-1)=u0(0:Ny,0:Nx-1,0:Nzloc-1)-uMT(0:Ny,0:Nx-1,0:Nzloc-1)
v0(0:Ny,0:Nx-1,0:Nzloc-1)=v0(0:Ny,0:Nx-1,0:Nzloc-1)-vMT(0:Ny,0:Nx-1,0:Nzloc-1)
w0(0:Ny,0:Nx-1,0:Nzloc-1)=w0(0:Ny,0:Nx-1,0:Nzloc-1)-wMT(0:Ny,0:Nx-1,0:Nzloc-1)
! Theta(0:Ny,0:Nx-1,0:Nzloc-1)=Theta(0:Ny,0:Nx-1,0:Nzloc-1)-ThetaMT(0:Ny,0:Nx-1,0:Nzloc-1)



!....Transorm to Fourier Space
call transform_rows_batched(uhat,u0(0:Ny,0:Nx-1,0:Nzloc-1))
call transform_rows_batched(vhat,v0(0:Ny,0:Nx-1,0:Nzloc-1))
call transform_rows_batched(what,w0(0:Ny,0:Nx-1,0:Nzloc-1))
! call transform_rows(Thetahat(0:Ny,0:Nx-1,0:Nzloc-1),Theta(0:Ny,0:Nx-1,0:Nzloc-1))


!....Removing also the mean field in space to consider the fluctuations only
!.....Mean fields in x and z


if(nid.eq.0)then
uhat(0:Ny,0,0)=cpx0
vhat(0:Ny,0,0)=cpx0
what(0:Ny,0,0)=cpx0
! Thetahat(0:Ny,0,0)=cpx0


end if

!....Transform back to Physical Space

call transformInv_rows_batched(u0(0:Ny,0:Nx-1,0:Nzloc-1),uhat)
call transformInv_rows_batched(v0(0:Ny,0:Nx-1,0:Nzloc-1),vhat)
call transformInv_rows_batched(w0(0:Ny,0:Nx-1,0:Nzloc-1),what)
! call transforminv_rows(Theta(0:Ny,0:Nx-1,0:Nzloc-1),Thetahat(0:Ny,0:Nx-1,0:Nzloc-1))


call MPI_BARRIER(MPI_Comm_World,ierror) 
!....Evaluating Statistics (Re Stress)
call ConvectiveStress(ConvStr,u0(0:Ny,0:Nx-1,0:Nzloc-1),v0(0:Ny,0:Nx-1,0:Nzloc-1)&
,w0(0:Ny,0:Nx-1,0:Nzloc-1))


! !....Evaluating Statistics (Re Stress and Passive scalar stress)
! call ConvectiveStressSc(ConvStr,ConvStrSc,u0(0:Ny,0:Nx-1,0:Nzloc-1),v0(0:Ny,0:Nx-1,0:Nzloc-1)&
! ,w0(0:Ny,0:Nx-1,0:Nzloc-1),Theta(0:Ny,0:Nx-1,0:Nzloc-1))



call EvaluateMaxuu(uumax(Nsavings),ConvStr(0:Ny,0:Nx-1,0:Nzloc-1,1))

!...Evaluate Turbulent Dissipation epsilon=1/Re_tau*<dui/dxj*dui/dxj>
call DissipationTerm(Dissip,uhat,vhat,what)
!....Evaluate Turbulent Transport
call TurbulentTransportTerm(TurbTr,u0(0:Ny,0:Nx-1,0:Nzloc-1),v0(0:Ny,0:Nx-1,0:Nzloc-1),&
w0(0:Ny,0:Nx-1,0:Nzloc-1))


ReStress=ReStress+ConvStr/2.0d0
! TheStress=TheStress+ConvStrSc/2.0d0


!...Mean Dissipation
DissipMT=DissipMT+Dissip/2.0d0
!...Mean Turbulent Transport
TurbTrMT=TurbTrMT+TurbTr/2.0d0


ReStress=ReStress/Ninterv
! TheStress=TheStress/Ninterv

!..Divide the dissipation by Re_tau
DissipMT=nu*DissipMT/Ninterv
TurbTrMT=TurbTrMT/Ninterv

!....Mean Dissipation in x and z
call MeanXZ(epsilonM,DissipMT)

call MPI_BARRIER(MPI_Comm_World,ierror) 


if(nid.eq.0)then
WRITE(*,*)'Turbulence Statistics and Turbulent Kinetic Energy Budget Evaluated.'
WRITE(*,*)'**********************************************************'
WRITE(*,*)'**********************************************************'  
end if
call MPI_BARRIER(MPI_Comm_World,ierror)  


ReStressPh=0.0d0
! TheStressPh=0.0d0




call MPI_BARRIER(MPI_Comm_World,ierror) 
if(nid.eq.0)then
WRITE(*,*)
WRITE(*,*)'**********************************************************'
WRITE(*,*)'**********************************************************' 
WRITE(*,*)'Saving mean profiles'
 
end if

uhat=cpx0
vhat=cpx0
what=cpx0
ConvStr=cpx0


! Thetahat=cpx0
! ConvStrSc=cpx0



if(nid.eq.0)then
uhat(0:Ny,0,0)=uhatMT(0:Ny,0,0)
vhat(0:Ny,0,0)=matmul(D(0:Ny,0:Ny),uhat(0:Ny,0,0))
what(0:Ny,0,0)=whatMT(0:Ny,0,0)
ConvStr(0:Ny,0,0,1)=ReStress(0:Ny,0,0,1)
ConvStr(0:Ny,0,0,2)=ReStress(0:Ny,0,0,2)
ConvStr(0:Ny,0,0,3)=ReStress(0:Ny,0,0,3)
ConvStr(0:Ny,0,0,4)=ReStress(0:Ny,0,0,4)
ConvStr(0:Ny,0,0,5)=ReStress(0:Ny,0,0,5)
ConvStr(0:Ny,0,0,6)=ReStress(0:Ny,0,0,6)

! Thetahat(0:Ny,0,0)=ThetahatMT(0:Ny,0,0)
! ConvStrSc(0:Ny,0,0,1)=TheStress(0:Ny,0,0,1)
! ConvStrSc(0:Ny,0,0,2)=TheStress(0:Ny,0,0,2)
! ConvStrSc(0:Ny,0,0,3)=TheStress(0:Ny,0,0,3)
! ConvStrSc(0:Ny,0,0,4)=TheStress(0:Ny,0,0,4)

end if

call MPI_BARRIER(MPI_Comm_World,ierror) 

call transformInv_rows_batched(u0(0:Ny,0:Nx-1,0:Nzloc-1),uhat)!<U>
call transformInv_rows_batched(v0(0:Ny,0:Nx-1,0:Nzloc-1),vhat)!d<U>/dy
call transformInv_rows_batched(w0(0:Ny,0:Nx-1,0:Nzloc-1),what)!<W>
call transformInv_rows_batched(ReStressPh(0:Ny,0:Nx-1,0:Nzloc-1,1),ConvStr(0:Ny,0:Nx-1,0:Nzloc-1,1))!<uu>
call transformInv_rows_batched(ReStressPh(0:Ny,0:Nx-1,0:Nzloc-1,2),ConvStr(0:Ny,0:Nx-1,0:Nzloc-1,2))!<vv>
call transformInv_rows_batched(ReStressPh(0:Ny,0:Nx-1,0:Nzloc-1,3),ConvStr(0:Ny,0:Nx-1,0:Nzloc-1,3))!<ww>
call transformInv_rows_batched(ReStressPh(0:Ny,0:Nx-1,0:Nzloc-1,4),ConvStr(0:Ny,0:Nx-1,0:Nzloc-1,4))!<uv>
call transformInv_rows_batched(ReStressPh(0:Ny,0:Nx-1,0:Nzloc-1,5),ConvStr(0:Ny,0:Nx-1,0:Nzloc-1,5))!<uw>
call transformInv_rows_batched(ReStressPh(0:Ny,0:Nx-1,0:Nzloc-1,6),ConvStr(0:Ny,0:Nx-1,0:Nzloc-1,6))!<vw>


! call transforminv_rows(Theta(0:Ny,0:Nx-1,0:Nzloc-1),Thetahat(0:Ny,0:Nx-1,0:Nzloc-1))

! call transforminv_rows(TheStressPh(0:Ny,0:Nx-1,0:Nzloc-1,1),&
! ConvStrSc(0:Ny,0:Nx-1,0:Nzloc-1,1))

! call transforminv_rows(TheStressPh(0:Ny,0:Nx-1,0:Nzloc-1,2),&
! ConvStrSc(0:Ny,0:Nx-1,0:Nzloc-1,2))

! call transforminv_rows(TheStressPh(0:Ny,0:Nx-1,0:Nzloc-1,3),&
! ConvStrSc(0:Ny,0:Nx-1,0:Nzloc-1,3))

! call transforminv_rows(TheStressPh(0:Ny,0:Nx-1,0:Nzloc-1,4),&
! ConvStrSc(0:Ny,0:Nx-1,0:Nzloc-1,4))

! dThetady(0:Ny)=matmul(D(0:Ny,0:Ny),Theta(0:Ny,0,0))



!....Evaluate Turbulent Kinetic Energy
Kturb=0.0d0
Kturb(0:Ny)=0.5d0*(ReStressPh(0:Ny,0,0,1)+ReStressPh(0:Ny,0,0,2)&
+ReStressPh(0:Ny,0,0,3))! k_turb=1/2*(<uu>+<vv>+<ww>)

!....Evaluate the Production Term of the Turb. Kinetic En.
ProdTurb(0:Ny)=0.0d0
ProdTurb(0:Ny)=-ReStressPh(0:Ny,0,0,4)*v0(0:Ny,0,0)! Production_turb=-<uv>*d<U>/dy

!....Molecular Diffusion Term
DiffMolec=0.0d0
DiffMolec(0:Ny)=matmul(D2(0:Ny,0:Ny),Kturb(0:Ny))
DiffMolec=nu*DiffMolec

!....Turbulent Transport (or Diffusion) term
DiffTurb=0.0d0
call EvalMeanTurbTranspor(DiffTurb,TurbTrMT)







call MPI_BARRIER(MPI_Comm_World,ierror) 
resfile(1)='flowstats.txt'
FMT1="(2X,E12.5,2X,E12.5,2X,E12.5,2X,E12.5,2X,E12.5,2X,E12.5,&
2X,E12.5,2X,E12.5,2X,E12.5,2X,E12.5,2X,E12.5)"

if(nid.eq.0)then
open(2,file=resfile(1))
write(2,*) '%       y        y+       <U>       d<U>/dy     <W>      <uu>      <vv>      <ww>      <uv>    &
<uw>      <vw>    '
do ii=0,Ny
write(2,FMT1) y(Ny-ii),yp(ii),u0(Ny-ii,0,0),v0(Ny-ii,0,0),w0(Ny-ii,0,0),ReStressPh(Ny-ii,0,0,1),&
ReStressPh(Ny-ii,0,0,2),ReStressPh(Ny-ii,0,0,3),ReStressPh(Ny-ii,0,0,4),ReStressPh(Ny-ii,0,0,5),ReStressPh(Ny-ii,0,0,6)

end do

close(2)

end if


resfile(2)='kbudget.txt'
FMT1="(2X,E12.5,2X,E12.5,2X,E12.5,2X,E12.5,2X,E12.5,2X,E12.5,&
2X,E12.5,2X,E12.5,2X,E12.5)"

if(nid.eq.0)then
open(3,file=resfile(2))
write(3,*) '%       y            y+       k       Produc.      Dissp.   &
  P_transp.      Turb_diff.      Molec_diff.     Residual'
do ii=0,Ny
write(3,FMT1) y(Ny-ii),yp(ii),Kturb(Ny-ii),ProdTurb(Ny-ii),epsilonM(Ny-ii),&
DiffPress(Ny-ii),DiffTurb(Ny-ii),DiffMolec(Ny-ii),DiffPress(Ny-ii)


end do

close(3)

end if


call MPI_BARRIER(MPI_Comm_World,ierror) 


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

end program post

