module in_outpost
use mpi
use parameterspost
use flowvariablespost

implicit none


contains

!***********************************************************************
!***********************************************************************
subroutine read_inputs()
integer::ii

!.....parameter are read by processor 0 from file in_dns_tc.txt...
!.....parameter are read by all processors from file in_dns_tc.txt...     
if(nid.eq.0)then
    
open(1,file='in_dns_tc.txt')
read(1,*)Lx_over_hpi
read(1,*)Lz_over_hpi
read(1,*)Nx
read(1,*)Ny
read(1,*)Nz
read(1,*)deltat
read(1,*)ntot
read(1,*)nsalva
read(1,*)Re
read(1,*) 
read(1,*)
read(1,*)
read(1,*)
read(1,*)
read(1,*)
close(1)

! open(1,file='in_dns_tc.txt')
! read(1,*)Nx
! read(1,*)Ny
! read(1,*)Nz
! read(1,*)deltat
! read(1,*)ntot
! read(1,*)nsalva
! read(1,*)Re
! read(1,*) 
! read(1,*)
! read(1,*)
! read(1,*)
! read(1,*)Sch
! read(1,*)
! read(1,*)Rayl
! read(1,*)gacc(1:2)
! read(1,*)
! close(1)

!open(1,file='in_dns_tc.txt')
!read(1,*)Nx
!read(1,*)Ny
!read(1,*)Nz
!read(1,*)deltat
!read(1,*)ntot
!read(1,*)nsalva
!read(1,*)Re
!close(1)



open(2,file='in_post.txt')
read(2,*)Nmin
read(2,*)Nmax
read(2,*)Ipasso
close(2)
      
      
!open(3,file='in_particles.txt')
!read(3,*)NParticles
!read(3,*)radiusG1
!read(3,*)radiusG2
!read(3,*)radiusG3
!read(3,*)rhorat
!read(3,*)cprat
!read(3,*)Froude
!read(3,*)
!read(3,*)
!read(3,*)
!read(3,*)
!read(3,*)
!read(3,*)
!close(3)

!open(2,file='in_particles.txt')
!read(2,*)NParticles
!read(2,*)radiusG1
!read(2,*)radiusG2
!read(2,*)radiusG3
!read(2,*)rhorat
!read(2,*)cprat
!read(2,*)Froude
!read(2,*)GenPart
!read(2,*)filePartG1
!read(2,*)fileTopolG1
!read(2,*)filePartG2
!read(2,*)fileTopolG2
!read(2,*)filePartG3
!read(2,*)fileTopolG3
!close(2)


      
sim_time=deltat*ntot
 



WRITE(*,*)'                POST-DNSTC MPI PROGRAM                    '
WRITE(*,*)
WRITE(*,*)
WRITE(*,*)'******************** INPUT DATA **************************'
WRITE(*,*)'**********************************************************'
WRITE(*,*)'Domain size in x divided by h*pi:', Lx_over_hpi
WRITE(*,*)'Domain size in z divided by h*pi:', Lz_over_hpi
WRITE(*,*)'Number of points in x :', Nx
WRITE(*,*)'Number of points in y :', Ny
WRITE(*,*)'Number of points in z :', Nz
WRITE(*,*)'Reynolds_tau Number:',Re
WRITE(*,*)'Time-Step: ',deltat
WRITE(*,*)'Number of steps: ',ntot
WRITE(*,*)'Simulation Time:  ',sim_time
WRITE(*,*)'Steps for saving: ',nsalva
WRITE(*,*)
WRITE(*,*)
WRITE(*,*)'First step: ',Nmin
WRITE(*,*)'Last step: ',Nmax
WRITE(*,*)'Processing Frequency:  ',Ipasso

WRITE(*,*)'**********************************************************'
WRITE(*,*)'**********************************************************'  



end if


!.....spread to all processors...
call MPI_Bcast(Lx_over_hpi,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
call MPI_Bcast(Lz_over_hpi,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
call MPI_Bcast(Nx,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
call MPI_Bcast(Ny,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
call MPI_Bcast(Nz,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
call MPI_Bcast(ntot,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
call MPI_Bcast(nsalva,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
      
call MPI_Bcast(Re,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
call MPI_Bcast(deltat,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
call MPI_Bcast(sim_time,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror) 
      
      
call MPI_Bcast(Nmin,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
call MPI_Bcast(Nmax,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
call MPI_Bcast(Ipasso,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)



ncicli=CEILING(ntot*1.0d0/nsalva)
nu=1.0d0/Re
Nzloc=Nz/noprocs
deltat2=nsalva*deltat*Ipasso


!Count Nsavings
Nsavings=0
do ii=Nmin,Nmax,Ipasso
Nsavings=Nsavings+1
end do
Ninterv=Nsavings-1


!Npassi=ceiling((Nmax-Nmin+1)/dfloat(Ipasso))


if(nid.eq.0)then
WRITE(*,*)'Actual deltat : ',deltat2
WRITE(*,*)'**********************************************************'
WRITE(*,*)'**********************************************************'  
end if



end subroutine
!***********************************************************************
!***********************************************************************
subroutine SaveScalarField(ucomp,resultfile)
implicit none
real(kind=prec):: ucomp(0:Ny,0:Nx-1,0:Nzloc-1)
character(len=clen):: resultfile
integer::ierror
integer::idproc,j1,j2,j3

write(*,*) 'Ready to save data, proc ',nid
call MPI_BARRIER(MPI_Comm_World,ierror) 

do idproc=0,noprocs-1

if(idproc.eq.nid)then 
open(2,file=resultfile,form='unformatted',access='direct',recl=prec*Nx*(Ny+1)*Nzloc)
write(2,rec=idproc+1) (((ucomp(j1,j2,j3),j1=0,Ny),j2=0,Nx-1),j3=0,Nzloc-1)
close(2)
!write(*,*) 'Data saved by proc',nid
end if 
call MPI_BARRIER(MPI_Comm_World,ierror) 
end do

write(*,*) 'End of storage ',resultfile

return
end subroutine
!***********************************************************************
!***********************************************************************
subroutine LoadScalarField(ucomp,namefilein)
implicit none
real(kind=prec):: ucomp(0:Ny,0:Nx-1,0:Nzloc-1)
character(len=clen):: namefilein
integer::ierror
integer::idproc,j1,j2,j3


!write(*,*) 'Ready to read data, proc',nid
call MPI_BARRIER(MPI_Comm_World,ierror) 

do idproc=0,noprocs-1

if(idproc.eq.nid)then 
open(3,file=namefilein,form='unformatted',access='direct',recl=prec*Nx*(Ny+1)*Nzloc)
read(3,rec=idproc+1) (((ucomp(j1,j2,j3),j1=0,Ny),j2=0,Nx-1),j3=0,Nzloc-1)
close(3)
!write(*,*) 'Data read by proc',nid
end if 
call MPI_BARRIER(MPI_Comm_World,ierror) 

end do

!write(*,*) 'End of reading ',namefilein
call MPI_BARRIER(MPI_Comm_World,ierror) 
return
end subroutine
!***********************************************************************
!***********************************************************************
!***********************************************************************
!***********************************************************************
subroutine SaveSpectrumZ(ucomp,resultfile)
implicit none
real(kind=prec):: ucomp(0:Ny,0:Nzloc-1)
character(len=clen):: resultfile
integer::ierror
integer::idproc,j1,j3

!write(*,*) 'Ready to save data, proc ',nid
call MPI_BARRIER(MPI_Comm_World,ierror) 

do idproc=0,noprocs-1

if(idproc.eq.nid)then 
open(2,file=resultfile,form='unformatted',access='direct',recl=prec*(Ny+1)*Nzloc)
write(2,rec=idproc+1) ((ucomp(j1,j3),j1=0,Ny),j3=0,Nzloc-1)
close(2)
!write(*,*) 'Data saved by proc',nid
end if 
call MPI_BARRIER(MPI_Comm_World,ierror) 
end do

!write(*,*) 'End of storage ',resultfile

call MPI_BARRIER(MPI_Comm_World,ierror) 

return
end subroutine
!***********************************************************************
!***********************************************************************


end module in_outpost
