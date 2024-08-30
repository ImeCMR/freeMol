! This program will calculate number of free drug molecules
!         (Molecules which are not belong to a clump)
!
!
 program freeMol
 implicit none
!
! Defining needed variables
!
 integer::ndrug,natom,nfree,i,j,k,ii,n,clscut,ncount
 real,dimension(100)::xmean,ymean,zmean
 real, dimension(99,99)::rmean
 real::x,y,z,xx,yy,zz,xm,ym,zm,dx,dy,dz,rcut,rij
 character(len=20)::run,gout
!
  write(6,*)'==================================='
  write(6,*)''
  write(6,*)'  Calculating free drug molecules  '
  write(6,*)''
  write(6,*)'==================================='
!
! Taking needed inputs from freeMol.dat
!
 open(unit=10,file='freeMol.dat',status='old')
!
 read(10,*) ndrug                ! Number of drug molecules inside the box
 read(10,*) natom                ! Number of atoms in a drug molecule
 read(10,*) run                  ! Input pdb file
 read(10,*) gout                 ! General output file
 read(10,*) clscut               ! Minimum number of molecules needed for a cluster
 read(10,*) rcut                 ! Cutoff value for distance between two molecules
!
! Writing data into a ouput file
!
 open(unit=11,file=gout,status='unknown')
!
 write(11,500) ndrug,natom,run,gout,clscut,rcut
! 
 500 format(/,2x,'Number of drug molecules                           = ',i5,/,&
              2x,'Number of atoms in a molecule                      = ',i5,/,&
              2x,'Input pdb file                                     = ',a20,/,&
              2x,'General output file                                = ',a20,/,&
              2x,'Minimum number of molecules per a cluster          = ',i5,/,&
              2x,'Cutoff value for distance between two molecules(A) = ',f7.3,/)
!
! Staring calculations
!
! Taking the mean points for molecules
!
 xx=0
 yy=0
 zz=0
 open(unit=12,file=run,status='old') 
 read(12,*)
 read(12,*)
 read(12,*)
 read(12,*)
 do i = 1,(ndrug*natom)
   read(12,'(30x,3f8.3)') x,y,z
!   write(6,*) x,y,z
   xx = xx+x
   yy = yy+y
   zz = zz+z
   !Calculating mean values for x,y,z for ndrug molecules
   if (mod(i,natom) == 0) then
     xm =xx/natom
     ym =yy/natom
     zm =zz/natom
!     write(6,*) xm,ym,zm
     !
     ! Writing mean values in to the mean arrays (1D)
     !
     k = i/natom
     xmean(k) = xm
     ymean(k) = ym
     zmean(k) = zm
     xx=0
     yy=0
     zz=0
   end if  
 end do
!
! Calculating mean distance between two molecules
!
 do i = 1,ndrug-1
   do j = (i+1),ndrug
     dx = xmean(i) - xmean(j)
     dy = ymean(i) - ymean(j)
     dz = zmean(i) - zmean(j)
     !
     rmean((j-1),i) = sqrt(dx**2+dy**2+dz**2) ! Saving mean distance into a 2D array
   end do
   !
 end do
!
! Calculating free molecules
!
 nfree=0
 do i = 1,ndrug-1
 n=0
   do j = 1,ndrug-1
     rij = rmean(i,j)
     if (rij<rcut .and. rij/=0) then   !Check whether the molecule is inside the sphere and ignoring null values in the array
         n = n+1   ! Counting the number of molecules
     end if
   end do
   if (n>1) then
     nfree=nfree+1
   end if
!
 end do
  write(6,*)'===================================='
  write(6,*)''
  write(6,*)' Number of free drug-like molecules '
  write(6,*)                nfree
  write(6,*)''
  write(6,*)'===================================='
end program freeMol
