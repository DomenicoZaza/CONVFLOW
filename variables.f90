!***********************************************************************
!***********************************************************************
module parameters
use mpi

implicit none

!precision definition
integer, parameter:: prec=8

!char length
integer, parameter:: clen=32

! !...Maximum number of passive scalars
! integer,parameter:: NPSmax=9


!....wavenumbers
real(kind=prec),allocatable, dimension(:)::kx,kz,kzloc
complex(kind=prec),allocatable, dimension(:,:,:):: ikkx,ikkz
real(kind=prec),allocatable,dimension(:,:):: kkquad
real(kind=prec),allocatable,dimension(:,:)::kkquadno0

!....dealiasing array
logical,allocatable,dimension(:,:,:)::dealiasing

!!....sigma array
!real(kind=prec),allocatable,dimension(:,:)::sigma,sigma2


!....Chebyshev Differentiation Matrices
real(kind=prec),allocatable, dimension(:,:)::D,D2,D2tilde
real(kind=prec),allocatable,dimension(:,:)::D2even,D2odd

!....Diagonalization of D2tilde
real(kind=prec),allocatable, dimension(:,:)::Qmat,QmatInv
real(kind=prec),allocatable, dimension(:)::Lambda,Lambda2

!....space coordinates
real(kind=prec),allocatable, dimension(:):: x,y,z,zloc
real(kind=prec),allocatable, dimension(:):: yp!y_plus

!....useful complex numbers
complex(kind=prec),parameter:: imu=DCMPLX(0.0d0,1.0d0)
complex(kind=prec),parameter:: cpx0=DCMPLX(0.0d0,0.0d0)

!....PI definition
real(kind=prec),parameter::pi=4.0d0*datan(1.0d0)


!....Number of grid points
integer(kind=prec)::Ny,Nx,Nz,Nzloc
real(kind=prec)::Lx_over_hpi,Lz_over_hpi
real(kind=prec)::Lx,Lz,Ly,Lzloc
real(kind=prec)::deltax,deltaz 

! !....Number of grid points
! integer(kind=prec)::Ny,Nx,Nz,Nzloc
! real(kind=prec)::Lx,Lz,Ly,Lzloc
! real(kind=prec)::deltax,deltaz 
! !!....Number of Passive Scalars
! !integer(kind=prec)::NPS

!....time-steps
real(kind=prec)::deltat,deltat2



!....viscosity
real(kind=prec)::Re,nu
! !....Parameters for Passive Scalars
! real(kind=prec)::Sch, nuPS
! !....Parameters for Mixed Convection
! real(kind=prec)::Rayl, nuMC
! real(kind=prec),allocatable,dimension(:)::gacc


!....time-integration utility
integer(kind=prec)::Nmin,Nmax,Ipasso,Nsavings,Ninterv
integer(kind=prec),allocatable,dimension(:)::IndSave
integer(kind=prec)::ntot
integer(kind=prec):: Npassi
integer(kind=prec)::nsalva
real(kind=prec)::sim_time
integer::ll1,Iext,iintern
integer(kind=prec)::ncicli

!....MPI size and rank
integer::nid,noprocs,nidp1,nidm1
integer::ierror


! !...Number of particles
! integer(kind=prec)::NParticles
! integer(kind=prec)::Nploc
! integer(kind=prec)::Nplocmax


! !....Parameters for Particles
! real(kind=prec)::radius,rhorat,taup
! real(kind=prec)::cprat,tauth
! integer(kind=prec)::PartInfo,PartInfoSave


! !....Indices for particle info
! integer(kind=prec)::Itag
! integer(kind=prec)::Ixp,Iyp,Izp
! integer(kind=prec)::Iup,Ivp,Iwp
! integer(kind=prec)::Ixp1,Iyp1,Izp1
! integer(kind=prec)::Iup1,Ivp1,Iwp1
! integer(kind=prec)::Iaxp1,Iayp1,Iazp1
! integer(kind=prec)::Ixp2,Iyp2,Izp2
! integer(kind=prec)::Iup2,Ivp2,Iwp2
! integer(kind=prec)::Iaxp2,Iayp2,Iazp2

! integer(kind=prec)::Ithetp,Ithetp1,Iqp1,Ithetp2,Iqp2

!....File Names
character(len=clen), allocatable,dimension(:):: nomefile,nomefile1
! character(len=clen):: fileScal,fileScal1
character(len=clen), allocatable,dimension(:):: resfile
character(len=clen):: meanfile,meanfileSc
character(len=clen)::snapshot,ScalInd
character(len=200):: FMT1


!....Flags (Save Reynolds stress or no)
integer::SaveReStr

complex(kind=prec),allocatable, dimension(:)::uvhmean

integer(kind=prec),allocatable,dimension(:)::nneven,nnodd

end module parameters
!***********************************************************************
!***********************************************************************
module flowvariables
use mpi 
use parameters
implicit none

!....Fields in Physical Space
real(kind=prec),allocatable, dimension(:,:,:)::u0,v0,w0
!....Time-averaged fields in Physical space
real(kind=prec),allocatable, dimension(:,:,:)::uMT,vMT,wMT

real(kind=prec),allocatable, dimension(:,:,:,:)::ReStressPh
!....Tutbulent Kinetic Enrergy
real(kind=prec),allocatable,dimension(:)::Kturb
real(kind=prec),allocatable,dimension(:)::ProdTurb
real(kind=prec),allocatable,dimension(:,:,:)::Dissip
real(kind=prec),allocatable,dimension(:,:,:)::DissipMT
real(kind=prec),allocatable,dimension(:,:,:)::TurbTr
real(kind=prec),allocatable,dimension(:,:,:)::TurbTrMT
real(kind=prec),allocatable,dimension(:,:,:)::PressTr
real(kind=prec),allocatable,dimension(:,:,:)::PressTrMT
real(kind=prec),allocatable,dimension(:)::epsilonM
real(kind=prec),allocatable,dimension(:)::DiffTurb,DiffMolec,DiffPress



!....Fields in Fourier Space
complex(kind=prec),allocatable, dimension(:,:,:)::uhat,vhat,what
complex(kind=prec),allocatable, dimension(:,:,:)::uhatp1,vhatp1,whatp1
complex(kind=prec),allocatable, dimension(:,:,:)::uhatMS,vhatMS,whatMS

!....Time-averaged fields in Fourier Space
complex(kind=prec),allocatable, dimension(:,:,:)::uhatMT,vhatMT,whatMT
complex(kind=prec),allocatable, dimension(:,:,:,:)::ReStress
complex(kind=prec),allocatable, dimension(:,:,:,:)::ConvStr
complex(kind=prec),allocatable, dimension(:,:,:,:)::ConvStrp1


! !....Variables for Passive Scalars
! real(kind=prec),allocatable, dimension(:,:,:)::Theta
! real(kind=prec),allocatable, dimension(:,:,:)::ThetaMT
! real(kind=prec),allocatable, dimension(:,:,:)::ThetaSpz
! real(kind=prec),allocatable, dimension(:,:)::ThetaSpz2
! real(kind=prec),allocatable, dimension(:,:)::ThetaSpzMT
! complex(kind=prec),allocatable, dimension(:,:,:)::Thetahat
! complex(kind=prec),allocatable, dimension(:,:,:)::ThetahatMT
! complex(kind=prec),allocatable, dimension(:,:,:)::ThetahatSpz

! complex(kind=prec),allocatable, dimension(:,:,:)::Thetahatp1
! complex(kind=prec),allocatable, dimension(:,:,:,:)::ConvStrSc
! complex(kind=prec),allocatable, dimension(:,:,:,:)::ConvStrScp1
! complex(kind=prec),allocatable, dimension(:,:,:,:)::TheStress
! real(kind=prec),allocatable, dimension(:,:,:,:)::TheStressPh
! real(kind=prec),allocatable, dimension(:)::dThetady


real(kind=prec),allocatable,dimension(:)::Ubulk,uumax



end module flowvariables
!***********************************************************************
!***********************************************************************