!=========================================================================
!This program reads in polarizability of a 1x1 layer system from chimat1x1
!and expand it to general super cell with a regular general k mesh
!
! Unit cell lattice: a1 a2;  and reciprocal lattice: b1 b2
! Super cell lattice: A1 = n1 * a1 + n2 * a2   
!                     A2 = m1 * a1 + m2 * a2 
! and reciprocal lattice: B1 = ( m2 * b1 - m1 * b2 )/(n1m2-n2m1)
!                         B2 = (-n2 * b1 + n1 * b2 )/(n1m2-n2m1)
! b1 = n1 * B1 + m1 * B2
! b2 = n2 * B1 + m2 * B2
! Any q point: q = X1 * B1 + X2 * B2 = x1 * b1 + x2 * b2
! X1 = n1 * x1 + n2 * x2;  x1 = ( m2 * X1 - n2 * X2 )/(n1m2-n2m1)
! X2 = m1 * x1 + m2 * x2;  x2 = (-m1 * X1 + n1 * X2 )/(n1m2-n2m1)
! 
! Input: n1 n2 m1 m2; supercel q=mesh = p x q 
! 
! Output: polarizability of the expanded system is written into file chi0mat_expandmxm
! and compared with chi0mat_exactmxm
!=========================================================================

program chi0expand

  implicit none

  integer, parameter :: DP = kind(1.0d0)
  integer, parameter :: DPC = kind((1.0d0,1.0d0))
  real(DP), parameter :: BOHR = 0.52917721092_dp
  real(DP), parameter :: RYD = 13.60569253_dp
  real(DP), parameter :: PI_D = 3.1415926535897932384626433832795_dp
  real(DP), parameter :: E_D = 2.7182818284590452353602874713526_dp
  real(DP), parameter :: TOL_ZERO = 1.0d-12
  real(DP), parameter :: TOL_Small = 1.0d-6
!-------------------
  type symmetry
    integer :: ntran         !< number of operations in full group
    integer :: ntranq        !< number of operations in small group of q
    real(DP) :: rq(3)        !< The q-point this ntranq belongs to
    integer :: mtrx(3,3,48)  !< symmetry matrix
    real(DP) :: tnp(3,48)    !< fractional translations
    integer :: indsub(48)    !< symmetry operations in subgroup of q
    integer :: kgzero(3,48)  !< Umklapp vectors for subgroup symmetry operations
    integer :: cell_symmetry !< 0 = cubic, 1 = hexagonal
  end type symmetry
!------------------
  type gspace
    integer :: ng       !< number of G-vectors
    integer :: nFFTgridpts !< number in FFT grid = product(FFTgrid(1:3))
    real(DP) :: ecutrho !< charge-density cutoff, in Ry
    integer, pointer :: components(:,:) !< the G-vectors, in units of 2 pi / a
    integer :: FFTgrid(3)  !< gsm: FFTgrid is the size of the FFT grid, not the maximum G-vector   
    integer, pointer :: index_vec(:) ! mapping to FFT grid
    real(DP), pointer :: ekin(:) !< kinetic energy of G-vectors
  end type gspace
!------------------
  type (symmetry) :: syms
  type (gspace) :: gvec1, gvec2, gvec0
  logical :: use_q0, slab, initial_read, initial_readp,findqx,findqy
  logical :: scismetal,ucismetal
  character :: dummy,calcname*12
  character :: aheadinput*60, ajname*6, adate*11, atime*14
  integer :: p,pp,q,n1,n2,m1,m2,nspin
  integer :: nq,nqs,nq1,nq2,nq0,ngq,n,k,i,j,ii,jj,nmtx1,ig,igp,iq,iqs,ioerr,allocerr,np
  integer :: freq_dep,nFreq,nband,nrk,icount,icount1,icount2, ndum, ntemp,igs,igps,m,mm
  integer :: kmax(3), nmax(3), nfft(6), maxnfft(6), qgrid(3), qkv(3), qkvp(3), qkv2(3), qkv2p(3)
  integer, allocatable :: isrtx1(:), isrtx2(:,:),irow(:),nmtx2(:),indp(:)
  real(DP) :: ecuts,temp1,temp2,temp3,temp4
  real(DP) :: bdot(3,3)
  real(DP), allocatable :: qpt1(:,:), qpt2(:,:), q_fzeps(:,:)
  complex(DPC), allocatable :: chimat1_ggp(:,:),chimat2_ggp(:,:,:),chimat_ggp(:,:),chimat3_ggp(:,:)
  complex(DPC), allocatable :: spchimat1_ggp(:,:),spchimat2_ggp(:,:,:),spchimat_ggp(:,:),spchimat3_ggp(:,:)

!------------program start-----------------
  write(*,*) '*******************Program Start*******************'
  write(*,*) '--------------reading input---------------'
  open(14,file='chimat_expandmxm',form='unformatted',status='replace')
  open(15,file='chi0mat_expandmxm',form='unformatted',status='replace')

!read kmesh in sc pxqx1
  p = 6
  q = 6
  n1 = 3 
  n2 = 1
  m = 4
  m1 =-1
  m2 = 2
read(*,*) p, q, n1,n2,m1,m2
read(*,*) scismetal,ucismetal
read(*,*) nspin
write(*,*) 'scismetal,ucismetal',scismetal,ucismetal

  mm = n1*m2 - n2*m1
if(mm .lt. 0) then
  write(*,*)'Be careful that the determinant is negative'
  mm = -mm 
endif
  nq = mm*p*q
  nq0 = 1
 write(*,*) 'total number of q points in unit cell'
 write(*,*) nq
 write(*,*) 'total number of q points in super cell'
  nqs = p*q
 write(*,*) nqs
      allocate (indp(mm))
!
!!!--------------read chi0mat1x1 to chimat2_ggp(:,:,1)
!once for all
      allocate (nmtx2(nq))
  write(*,*) '-------start reading UnitCell chi matrix-------'
        open(12,file='chi0mat1x1',form='unformatted',status='old',iostat=ioerr)
        read(12)
        read(12)
        read(12)
        read(12)
        read(12)
        read(12)
        read(12)
        read(12)
        read(12) gvec0%ng
       backspace(12)
        allocate ( gvec0%components(3,gvec0%ng) )
        read(12) gvec0%ng,gvec0%nFFTgridpts,gvec0%FFTgrid(1:3), &
              ((gvec0%components(jj,ig),jj=1,3),ig=1,gvec0%ng)
!allocate once for all
        allocate (qpt2(3,nq))
        read(12) ((qpt2(jj,ii),jj=1,3),ii=1,nq0)
      read(12)  ! syms%
!once for all
      allocate (isrtx2(gvec0%ng,nq))
      allocate (gvec0%ekin(gvec0%ng))
      read(12) nmtx2(1),np,(isrtx2(i,1),gvec0%ekin(i),i=1,gvec0%ng)
      ntemp = 100+nmtx2(1)
!once for all
      allocate (chimat2_ggp (ntemp,ntemp,mm), stat=allocerr)
      if(nspin .eq. 2) then
      allocate (spchimat2_ggp (ntemp,ntemp,mm), stat=allocerr)
      endif
indp(1) = 1
      do igp = 1,nmtx2(1)
        read(12) (chimat2_ggp(ig,igp,1),ig=1,nmtx2(1))
      enddo
      if(nspin .eq. 2) then
      do igp = 1,nmtx2(1)
        read(12) (spchimat2_ggp(ig,igp,1),ig=1,nmtx2(1))
      enddo
      endif
      write(*,*) 'reading iq = 1'
!!!-------read chimat1x1 to chimat2_ggp(:,:,2-mm)
  initial_read = .true.
pp=1
  do iq = 2,nq !
      if (initial_read) then
        open(17,file='chimat1x1',form='unformatted',status='old',iostat=ioerr)
        read(17)
        read(17)
        read(17)
        read(17)
        read(17)
        read(17)
        read(17)
        read(17)
        read(17) gvec2%ng
       backspace(17)
        allocate ( gvec2%components(3,gvec2%ng) ) 
        read(17) gvec2%ng,gvec2%nFFTgridpts,gvec2%FFTgrid(1:3), &
              ((gvec2%components(jj,ig),jj=1,3),ig=1,gvec2%ng)
!check
      if(gvec0%ng .ne. gvec2%ng) then
        write(*,*) 'STOP: gvec0%ng .ne. gvec2%ng'
        stop
      endif
      do ig = 1 , gvec2%ng
       if(all( abs(gvec0%components(1:3,ig)-gvec2%components(1:3,ig)) .lt. TOL_Small )) then   
       else
        write(*,*) 'very likely to STOP: gvec0 .ne. gvec2'
        stop
       endif
      enddo
!endcheck
        read(17) nq2
       backspace(17)
!        allocate (qpt2(3,nq2))
      if (ucismetal) then
        read(17) nq2,nq0,((qpt2(jj,ii),jj=1,3),ii=2,nq)
        nq0 = 1
        write(*,*) 'nq2',nq2
       do ii = 1,nq2+1
        write(6,'(3f10.6)') qpt2(:,ii)
       enddo
      else
        read(17) nq2,nq0,((qpt2(jj,ii),jj=1,3),ii=1,nq2)
        write(*,*) 'nq2',nq2
       do ii = 1,nq2
        write(6,'(3f10.6)') qpt2(:,ii)
       enddo
      endif
        initial_read = .false.
      endif
! Read q-dependent info
      read(17)
!      allocate (isrtx2(gvec2%ng))
      allocate (gvec2%ekin(gvec2%ng))
      read(17) nmtx2(iq),np,(isrtx2(i,iq),gvec2%ekin(i),i=1,gvec2%ng)

        if (nmtx2(iq) .gt. ntemp) then
        write(*,*) 'STOP: nmtx2 .gt. ntemp'
        stop
       endif
      write(*,*) "reading iq =",iq

!check this for larger SC than 2x2
temp1 = MOD(abs(n1*qpt2(1,iq)+n2*qpt2(2,iq))+0.0001,1.0)
temp2 = MOD(abs(m1*qpt2(1,iq)+m2*qpt2(2,iq))+0.0001,1.0)
 if(temp1.lt.0.01 .and. temp2.lt.0.01)then
! if(MOD(m*qpt2(1,iq)+0.0001,1.0).lt.0.01 .and. MOD(m*qpt2(2,iq)+0.0001,1.0).lt.0.01)then
      pp=pp+1
 indp(pp) = iq
        write(6,'(3f10.6)') qpt2(:,iq)
      do igp = 1,nmtx2(iq)
        read(17) (chimat2_ggp(ig,igp,pp),ig=1,nmtx2(iq))
      enddo
      if(nspin .eq. 2) then
      do igp = 1,nmtx2(iq)
        read(17) (spchimat2_ggp(ig,igp,pp),ig=1,nmtx2(iq))
      enddo
      endif
 else
      do igp = 1,nmtx2(iq)
        read(17) 
      enddo
      if(nspin .eq. 2) then
      do igp = 1,nmtx2(iq)
        read(17)
      enddo
      endif
 endif
!
    deallocate (gvec2%ekin)
  enddo !iq
  write(*,*) '--------done reading UnitCell chi matrix-------'

!!!----------------done reading UnitCell chi matrix
!!!!!!!!!!!!!=======start to expand q=0 ===============
 write(*,*) '-------------start to read chi0mat_exactmxm----- '
!--------read chi0mat_exactmxm q=0 to chimat1_ggp----------------------
        open(11,file='chi0mat_exactmxm',form='unformatted',status='old',iostat=ioerr)
        read(11) aheadinput,ajname,adate
        read(11) freq_dep,nFreq
        read(11) qgrid(1:3)
        read(11)
        read(11)
        read(11)
        read(11) ecuts,nband
        read(11) nrk
        read(11) gvec1%ng,gvec1%nFFTgridpts

       backspace(11)
        allocate ( gvec1%components(3,gvec1%ng) )
        allocate ( gvec1%index_vec(gvec1%nFFTgridpts) )
        read(11) gvec1%ng,gvec1%nFFTgridpts,gvec1%FFTgrid(1:3), &
              ((gvec1%components(jj,ig),jj=1,3),ig=1,gvec1%ng), &
              ((bdot(ii,jj),jj=1,3),ii=1,3),(gvec1%index_vec(ig),ig=1,gvec1%nFFTgridpts)
!once for all
        allocate (qpt1(3,nqs))
        read(11) ((qpt1(jj,ii),jj=1,3),ii=1,nq0)
      ! Read q-dependent info
      read(11) syms%ntranq
      backspace(11)
      read(11) syms%ntranq,(((syms%mtrx(i,j,syms%indsub(n)),i=1,3),j=1,3), &
            (syms%tnp(k,syms%indsub(n)),syms%kgzero(k,n),k=1,3),n=1,syms%ntranq)
!          np=pol%nmtx*(pol%nmtx+1)/2
      allocate (isrtx1(gvec1%ng))
      allocate (gvec1%ekin(gvec1%ng))
      read(11) nmtx1
      allocate (irow(nmtx1))
      backspace(11)
      read(11) nmtx1,np,(isrtx1(i),gvec1%ekin(i),i=1,gvec1%ng),(irow(i),i=1,nmtx1)
!--------chi0mat_expandmxm q=0--------------------------------
        write(15) aheadinput,ajname,adate
        write(15) freq_dep,nFreq
        write(15) qgrid(1:3)
        write(15)
        write(15)
        write(15)
        write(15) ecuts,nband
        write(15) nrk, 1
        write(15) gvec1%ng,gvec1%nFFTgridpts,gvec1%FFTgrid(1:3),&
              ((gvec1%components(jj,ig),jj=1,3),ig=1,gvec1%ng), &
              ((bdot(ii,jj),jj=1,3),ii=1,3),(gvec1%index_vec(ig),ig=1,gvec1%nFFTgridpts)
        write(15) ((qpt1(jj,ii),jj=1,3),ii=1,nq0)

      write(15) syms%ntranq
      write(15) nmtx1
      allocate (chimat1_ggp (nmtx1,nmtx1), stat=allocerr)
      allocate (chimat_ggp (nmtx1,nmtx1), stat=allocerr)
      if(nspin .eq. 2) then
      allocate (spchimat1_ggp (nmtx1,nmtx1), stat=allocerr)
      allocate (spchimat_ggp (nmtx1,nmtx1), stat=allocerr)
      endif
      do igp = 1,nmtx1
        read(11) (chimat1_ggp(ig,igp),ig=1,nmtx1)
        do ig = 1,nmtx1
           chimat_ggp(ig,igp)=CMPLX(0.0)
        enddo
      enddo
      if(nspin .eq. 2) then
      do igp = 1,nmtx1
        read(11) (spchimat1_ggp(ig,igp),ig=1,nmtx1)
        do ig = 1,nmtx1
           spchimat_ggp(ig,igp)=CMPLX(0.0)
        enddo
      enddo
      endif
      write(*,*)'complete reading chi0mat_exactmxm to chimat1_ggp and initiating chimat_ggp'
      write(*,*)'------start expanding chi0mat_expandmxm-----' 
!nmtx2 number of G vectors in 1x1
      write(*,*) 'nmtx1',nmtx1
icount = 0 
     do iq = 1,mm
        write(*,*) iq
        write(6,'(3f10.6)') qpt2(:,indp(iq))
        write(6,'(3f10.6)') qpt1(:,1)
icount = icount + nmtx2(indp(iq))*nmtx2(indp(iq))
! Nic enable omp
!$omp parallel do default(private) shared(chimat_ggp,spchimat_ggp,chimat2_ggp,spchimat2_ggp,nspin,iq,nmtx1,nmtx2,indp,gvec1,gvec2,isrtx1,isrtx2,n1,n2,m1,m2,qpt2)
      do igp = 1,nmtx2(indp(iq))
        do ig = 1,nmtx2(indp(iq)) !nmtx2
           qkv(:) = gvec2%components(:,isrtx2(ig,indp(iq)))
           qkvp(:)= gvec2%components(:,isrtx2(igp,indp(iq)))
!find the gvec coordinate in sc bz: Gs = G + G_I
qkv2(1) = n1*qkv(1) + n2*qkv(2) + NINT(n1*qpt2(1,indp(iq))+n2*qpt2(2,indp(iq)))
qkv2(2) = m1*qkv(1) + m2*qkv(2) + NINT(m1*qpt2(1,indp(iq))+m2*qpt2(2,indp(iq)))
qkv2(3) = qkv(3)
qkv2p(1) = n1*qkvp(1) + n2*qkvp(2) + NINT(n1*qpt2(1,indp(iq))+n2*qpt2(2,indp(iq)))
qkv2p(2) = m1*qkvp(1) + m2*qkvp(2) + NINT(m1*qpt2(1,indp(iq))+m2*qpt2(2,indp(iq)))
qkv2p(3) = qkvp(3)
!           qkv(1:2) = m*qkv(1:2) + NINT(m*qpt2(1:2,indp(iq)))
!           qkvp(1:2)= m*qkvp(1:2)+ NINT(m*qpt2(1:2,indp(iq)))
!find index igs
           igs = findigs(gvec1,qkv2,isrtx1)
           igps = findigs(gvec1,qkv2p,isrtx1)         
if(igs .gt. nmtx1 .or. igps .gt. nmtx1) then
 write(*,*) 'igs or igps gt nmtx1'
 stop
else
endif
           chimat_ggp(igs,igps)=chimat2_ggp(ig,igp,iq)
      if(nspin .eq. 2) then
           spchimat_ggp(igs,igps)=spchimat2_ggp(ig,igp,iq)
      endif
        enddo ! ig
      enddo   !igp
!$omp end parallel do
     enddo    ! iq
!!!!!write expand q =0  results
      write(*,*)'------end expanding q=0 and start writing results-----'
      do igp = 1,nmtx1
        write(15) (chimat_ggp(ig,igp),ig=1,nmtx1)
      enddo   !igp
      if(nspin .eq. 2) then
      do igp = 1,nmtx1
        write(15) (spchimat_ggp(ig,igp),ig=1,nmtx1)
      enddo   !igp
      endif
!check
      write(*,*) '----------begin checking----------'
      write(*,*) 'nmtx1^2',nmtx1*nmtx1
      write(*,*) 'total number of nonzero elements',icount
      write(*,*) 'total number of zero elements',nmtx1*nmtx1-icount
icount = 0
do iq = 1,nq
 icount = icount+ nmtx2(iq)*nmtx2(iq)
 write(*,*) 'nmtx2^2',iq,nmtx2(iq),nmtx2(iq)*nmtx2(iq)
enddo
write(*,*) 'nmtx1^2 - sum nmtx2^2',nmtx1*nmtx1 - icount
          icount = 0
          icount1 = 0
          icount2 = 0
      do igp = 1,nmtx1
       do ig = 1,nmtx1
         temp1 = REAL(chimat_ggp(ig,igp))
         temp2 = REAL(chimat1_ggp(ig,igp))
         temp3 = abs((temp1 - temp2 )/temp2)
          if (abs(temp2) .lt. 0.000003)then
          icount = icount +1
 if(MOD(gvec1%components(1,isrtx1(ig))-gvec1%components(1,isrtx1(igp)),m) .eq. 0)then
 if(MOD(gvec1%components(2,isrtx1(ig))-gvec1%components(2,isrtx1(igp)),m) .eq. 0)then
          icount1 = icount1 +1
 endif
 endif
          else
!         if (abs(temp2) .lt. 0.001 .and. temp3 .gt. 0.1)then
         if (temp3 .gt. 0.05)then
          icount2 = icount2+1
         endif
          endif
       enddo
      enddo
        write(*,*) '==>checking results:',icount,icount1,icount2

!endcheck
   deallocate (isrtx1)
    deallocate (gvec1%ekin)
    deallocate(gvec1%components)
    deallocate ( gvec1%index_vec)
    deallocate (chimat_ggp)
    deallocate (chimat1_ggp)
      if(nspin .eq. 2) then
    deallocate (spchimat_ggp)
    deallocate (spchimat1_ggp)
      endif
    deallocate (irow)
    close(12)
    close(17)
    deallocate (isrtx2)
    deallocate(gvec0%components)
    deallocate (gvec0%ekin)
    deallocate (nmtx2)
    deallocate(chimat2_ggp)
      if(nspin .eq. 2) then
    deallocate(spchimat2_ggp)
      endif
    deallocate(gvec2%components)
    deallocate (qpt2)
!!!!!!!!!!!!!=======done expand q=0 ===============
!!!!!!=====expand to all q .ne. 0
  initial_read = .true.
  do iqs = 2, nqs
        write(*,*) '----------start to read chimat_exactmxm-------',iqs
      if (initial_read) then
      open(16,file='chimat_exactmxm',form='unformatted',status='old',iostat=ioerr)
        read(16) aheadinput,ajname,adate
       write(14) aheadinput,ajname,adate
        read(16) freq_dep,nFreq
       write(14) freq_dep,nFreq
        read(16) qgrid(1:3)
       write(14) qgrid(1:3)
        read(16)
        read(16)
        read(16)
       write(14)
       write(14)
       write(14)
        read(16) ecuts,nband
        read(16) nrk
        read(16) gvec1%ng,gvec1%nFFTgridpts
       write(14) ecuts,nband
       write(14) nrk, 1
       backspace(16)
        allocate ( gvec1%components(3,gvec1%ng) )
        allocate ( gvec1%index_vec(gvec1%nFFTgridpts) )
        read(16) gvec1%ng,gvec1%nFFTgridpts,gvec1%FFTgrid(1:3), &
              ((gvec1%components(jj,ig),jj=1,3),ig=1,gvec1%ng), &
              ((bdot(ii,jj),jj=1,3),ii=1,3),(gvec1%index_vec(ig),ig=1,gvec1%nFFTgridpts)
        read(16) nq1
       backspace(16)
!        allocate (qpt1(3,nq1))
      if (scismetal) then
        read(16) nq1,nq0,((qpt1(jj,ii),jj=1,3),ii=2,nqs)
       write(14) gvec1%ng,gvec1%nFFTgridpts,gvec1%FFTgrid(1:3), &
              ((gvec1%components(jj,ig),jj=1,3),ig=1,gvec1%ng), &
              ((bdot(ii,jj),jj=1,3),ii=1,3),(gvec1%index_vec(ig),ig=1,gvec1%nFFTgridpts)
       write(14) nq1,nq0,((qpt1(jj,ii),jj=1,3),ii=2,nqs)
        write(*,*) 'number of q point is SuperCell',nq1
       do ii = 1,nq1+1
        write(6,'(3f10.6)') qpt1(:,ii)
       enddo
       nq0 = 1
      else
        read(16) nq1,nq0,((qpt1(jj,ii),jj=1,3),ii=1,nq1)
       write(14) gvec1%ng,gvec1%nFFTgridpts,gvec1%FFTgrid(1:3), &
              ((gvec1%components(jj,ig),jj=1,3),ig=1,gvec1%ng), &
              ((bdot(ii,jj),jj=1,3),ii=1,3),(gvec1%index_vec(ig),ig=1,gvec1%nFFTgridpts)
       write(14) nq1,nq0,((qpt1(jj,ii),jj=1,3),ii=1,nq1)
        write(*,*) 'number of q point is SuperCell',nq1
       do ii = 1,nq1
        write(6,'(3f10.6)') qpt1(:,ii)
       enddo
       endif
!write into
        initial_read = .false.
      endif
! Read q-dependent info
      read(16) syms%ntranq
      backspace(16)
      read(16) syms%ntranq,(((syms%mtrx(i,j,syms%indsub(n)),i=1,3),j=1,3), &
            (syms%tnp(k,syms%indsub(n)),syms%kgzero(k,n),k=1,3),n=1,syms%ntranq)
     write(14) syms%ntranq
      allocate (isrtx1(gvec1%ng), stat=allocerr)
      allocate (gvec1%ekin(gvec1%ng), stat=allocerr)
      read(16) nmtx1
      allocate (irow(nmtx1), stat=allocerr)
      backspace(16)
      read(16) nmtx1,np,(isrtx1(i),gvec1%ekin(i),i=1,gvec1%ng),(irow(i),i=1,nmtx1)
     write(14) nmtx1
      write(*,*) "nmtx1(iqs)",iqs,nmtx1
      allocate (chimat1_ggp (nmtx1,nmtx1), stat=allocerr)
      allocate (chimat_ggp (nmtx1,nmtx1), stat=allocerr)
      if(nspin .eq. 2) then
      allocate (spchimat1_ggp (nmtx1,nmtx1), stat=allocerr)
      allocate (spchimat_ggp (nmtx1,nmtx1), stat=allocerr)
      endif
      do igp = 1,nmtx1
        read(16) (chimat1_ggp(ig,igp),ig=1,nmtx1)
        do ig = 1,nmtx1
           chimat_ggp(ig,igp)=CMPLX(0.0)
        enddo
      enddo
      if(nspin .eq. 2) then
      do igp = 1,nmtx1
        read(16) (spchimat1_ggp(ig,igp),ig=1,nmtx1)
        do ig = 1,nmtx1
           spchimat_ggp(ig,igp)=CMPLX(0.0)
        enddo
      enddo
      endif
      write(*,*)'complete reading chimat_exactmxm to chimat1_ggp and initiating chimat_ggp'
      write(*,*) 'we are in the loop for q in SC'
        write(6,'(3f10.6)') qpt1(:,iqs)
!===================start reading UnitCell chi matrix=================
!once for all
      allocate (nmtx2(nq))
! read chi0mat1x1 to chimat2_ggp(:,:,1)
  write(*,*) '-------start reading UnitCell chi matrix-------'
        open(12,file='chi0mat1x1',form='unformatted',status='old',iostat=ioerr)
        read(12)
        read(12)
        read(12)
        read(12)
        read(12)
        read(12)
        read(12)
        read(12)
        read(12) gvec0%ng
       backspace(12)
        allocate ( gvec0%components(3,gvec0%ng) )
        read(12) gvec0%ng,gvec0%nFFTgridpts,gvec0%FFTgrid(1:3), &
              ((gvec0%components(jj,ig),jj=1,3),ig=1,gvec0%ng)

!allocate once for all
        allocate (qpt2(3,nq))
        read(12) ((qpt2(jj,ii),jj=1,3),ii=1,nq0)
      read(12)  ! syms%
!once for all
      allocate (isrtx2(gvec0%ng,nq))
      allocate (gvec0%ekin(gvec0%ng))
      read(12) nmtx2(1),np,(isrtx2(i,1),gvec0%ekin(i),i=1,gvec0%ng)
      ntemp = 100+nmtx2(1)
!once for all
      allocate (chimat2_ggp (ntemp,ntemp,mm), stat=allocerr)
      if(nspin .eq. 2) then
      allocate (spchimat2_ggp (ntemp,ntemp,mm), stat=allocerr)
      endif
! read chimat1x1 to chimat2_ggp(:,:,2-mm)
  initial_readp = .true.
pp=0
  do iq = 2,nq !
      if (initial_readp) then
        open(17,file='chimat1x1',form='unformatted',status='old',iostat=ioerr)
        read(17)
        read(17)
        read(17)
        read(17)
        read(17)
        read(17)
        read(17)
        read(17)
        read(17) gvec2%ng
       backspace(17)
        allocate ( gvec2%components(3,gvec2%ng) )
        read(17) gvec2%ng,gvec2%nFFTgridpts,gvec2%FFTgrid(1:3), &
              ((gvec2%components(jj,ig),jj=1,3),ig=1,gvec2%ng)
!check
!endcheck
        read(17) nq2
       backspace(17)
!        allocate (qpt2(3,nq2))
      if (ucismetal) then
        read(17) nq2,nq0,((qpt2(jj,ii),jj=1,3),ii=2,nq)
        nq0=1
        write(*,*) 'nq2',nq2
       do ii = 1,nq2+1
        write(6,'(3f10.6)') qpt2(:,ii)
       enddo
      else
        read(17) nq2,nq0,((qpt2(jj,ii),jj=1,3),ii=1,nq2)
        write(*,*) 'nq2',nq2
       do ii = 1,nq2
        write(6,'(3f10.6)') qpt2(:,ii)
       enddo
      endif
        initial_readp = .false.
      endif
! Read q-dependent info
      read(17)
!      allocate (isrtx2(gvec2%ng))
      allocate (gvec2%ekin(gvec2%ng))
      read(17) nmtx2(iq),np,(isrtx2(i,iq),gvec2%ekin(i),i=1,gvec2%ng)

        if (nmtx2(iq) .gt. ntemp) then
        write(*,*) 'STOP: nmtx2 .gt. ntemp'
        stop
       endif
      write(*,*) "reading iq",iq
!only for 2x2 SuperCell
!temp1 = MOD(abs(m*qpt2(1,iq)-qpt1(1,iqs))+0.0001,1.0)
!temp2 = MOD(abs(m*qpt2(2,iq)-qpt1(2,iqs))+0.0001,1.0)

temp1 = MOD(abs(n1*qpt2(1,iq)+n2*qpt2(2,iq)-qpt1(1,iqs))+0.0001,1.0)
temp2 = MOD(abs(m1*qpt2(1,iq)+m2*qpt2(2,iq)-qpt1(2,iqs))+0.0001,1.0)
 if(temp1.lt.0.01 .and. temp2.lt.0.01)then
!if(temp1.lt.0.01 .and. temp2.lt.0.01)then

      pp=pp+1
 indp(pp) = iq
        write(6,'(3f10.6)') qpt2(:,iq)
        write(*,*) indp(pp),pp
   if(pp .gt. mm) then
    stop
   endif
      do igp = 1,nmtx2(iq)
        read(17) (chimat2_ggp(ig,igp,pp),ig=1,nmtx2(iq))
      enddo
      if(nspin .eq. 2) then
      do igp = 1,nmtx2(iq)
        read(17) (spchimat2_ggp(ig,igp,pp),ig=1,nmtx2(iq))
      enddo
      endif
 else
      do igp = 1,nmtx2(iq)
        read(17)
      enddo
      if(nspin .eq. 2) then
      do igp = 1,nmtx2(iq)
        read(17)
      enddo
      endif
 endif
!
    deallocate (gvec2%ekin)
  enddo !iq
  write(*,*) '--------done reading UnitCell chi matrix-------'

!=======================done reading 1x1 chi matrix=======================

      write(*,*)'------start expanding chimat_expandmxm-----'
!expand procedure for q ne 0
  write(*,*) '--------------expanding q .ne. 0--------------'
icount = 0 
     do iq = 1,mm
        write(*,*) iq
        write(6,'(3f10.6)') qpt2(:,indp(iq))
        write(6,'(3f10.6)') qpt1(:,iqs)
icount = icount + nmtx2(indp(iq))*nmtx2(indp(iq))
! Nic enable omp
!$omp parallel do default(private) shared(chimat_ggp,spchimat_ggp,chimat2_ggp,spchimat2_ggp,nspin,iq,iqs,nmtx1,nmtx2,indp,gvec1,gvec2,isrtx1,isrtx2,n1,n2,m1,m2,qpt2,qpt1)
      do igp = 1,nmtx2(indp(iq))
        do ig = 1,nmtx2(indp(iq)) !nmtx2
           qkv(:) = gvec2%components(:,isrtx2(ig,indp(iq)))
           qkvp(:)= gvec2%components(:,isrtx2(igp,indp(iq)))
! problem: absolute value ???
!find the gvec coordinate in sc bz
qkv2(1) = n1*qkv(1) + n2*qkv(2) + NINT(n1*qpt2(1,indp(iq))+n2*qpt2(2,indp(iq))-qpt1(1,iqs))
qkv2(2) = m1*qkv(1) + m2*qkv(2) + NINT(m1*qpt2(1,indp(iq))+m2*qpt2(2,indp(iq))-qpt1(2,iqs))
qkv2(3) = qkv(3)
qkv2p(1) = n1*qkvp(1) + n2*qkvp(2) + NINT(n1*qpt2(1,indp(iq))+n2*qpt2(2,indp(iq))-qpt1(1,iqs))
qkv2p(2) = m1*qkvp(1) + m2*qkvp(2) + NINT(m1*qpt2(1,indp(iq))+m2*qpt2(2,indp(iq))-qpt1(2,iqs))
qkv2p(3) = qkvp(3)
!           qkv(1:2) = m*qkv(1:2) + NINT(m*qpt2(1:2,indp(iq))-qpt1(1:2,iqs))
!           qkvp(1:2)= m*qkvp(1:2)+ NINT(m*qpt2(1:2,indp(iq))-qpt1(1:2,iqs))
           igs = findigs(gvec1,qkv2,isrtx1)
           igps = findigs(gvec1,qkv2p,isrtx1)
if(igs .gt. nmtx1 .or. igps .gt. nmtx1) then
 write(*,*) 'igs or igps gt nmtx1'
 write(*,*) igs, igps, nmtx1
 stop
else
           chimat_ggp(igs,igps)=chimat2_ggp(ig,igp,iq)
      if(nspin .eq. 2) then
           spchimat_ggp(igs,igps)=spchimat2_ggp(ig,igp,iq)
      endif
endif
!           chimat_ggp(igs,igps)=chimat2_ggp(ig,igp,iq)
        enddo ! ig
      enddo   !igp
!$omp end parallel do
     enddo    ! iq
!
      write(*,*)'------end expanding q ne 0 and start writing results-----'
!write to file
      do igp = 1,nmtx1
        write(14) (chimat_ggp(ig,igp),ig=1,nmtx1)
      enddo
      if(nspin .eq. 2) then
      do igp = 1,nmtx1
        write(14) (spchimat_ggp(ig,igp),ig=1,nmtx1)
      enddo
      endif
!begincheck
      write(*,*) '----------begin checking q ne 0----------'
      write(*,*) 'nmtx1^2,iqs,nmtx1',nmtx1*nmtx1,iqs,nmtx1
      write(*,*) 'total number of nonzero elements',icount
      write(*,*) 'total number of zero elements',nmtx1*nmtx1-icount
          icount = 0
          icount1 = 0
          icount2 = 0
      do igp = 1,nmtx1
       do ig = 1,nmtx1
         temp1 = REAL(chimat_ggp(ig,igp))
         temp2 = REAL(chimat1_ggp(ig,igp))
         temp3 = abs((temp1 - temp2 )/temp2)
          if (abs(temp2) .lt. 0.000003)then
          icount = icount +1
 if(MOD(gvec1%components(1,isrtx1(ig))-gvec1%components(1,isrtx1(igp)),m) .eq. 0)then
 if(MOD(gvec1%components(2,isrtx1(ig))-gvec1%components(2,isrtx1(igp)),m) .eq. 0)then
          icount1 = icount1 +1
 endif
 endif
          else
!         if (abs(temp2) .lt. 0.001 .and. temp3 .gt. 0.1)then
         if (temp3 .gt. 0.05)then
          icount2 = icount2+1
         endif
          endif
       enddo
      enddo
        write(*,*) '==> checking results:',icount,icount1,icount2
!endcheck
    deallocate (isrtx1)
    deallocate (gvec1%ekin)
    deallocate (chimat_ggp)
    deallocate (chimat1_ggp)
      if(nspin .eq. 2) then
    deallocate (spchimat_ggp)
    deallocate (spchimat1_ggp)
      endif
    deallocate (irow)
   close(12)
   close(17)
    deallocate(gvec0%components)
    deallocate (isrtx2)
    deallocate (gvec0%ekin)
    deallocate (nmtx2)
    deallocate(chimat2_ggp)
      if(nspin .eq. 2) then
    deallocate(spchimat2_ggp)
      endif
    deallocate(gvec2%components)
    deallocate (qpt2)

  enddo !iqs 2-nqs

!-------Cleaning----------
    deallocate (gvec1%components)
    deallocate (gvec1%index_vec)
    deallocate (qpt1)
    deallocate (indp)

  close(11)
  close(14)
  close(15)
  close(16)
!  write(*,*) "Done reading and writing"
!----------------------------------
  write(*,*) "*******************Finish!!!********************"
!----------------define function-----------------------------------
contains 
integer function findigs(gvecf,qkvf,isortx)
  implicit none
  integer, intent(in) :: qkvf(3)
  integer, intent(in) :: isortx(:)
  type (gspace), intent(in) :: gvecf
  integer :: igf
!Problem: without exactmxm chimat, how to find isrtx1???
   do igf = 1, gvecf%ng  
    if (all( abs(gvecf%components(1:3,isortx(igf))-qkvf(1:3)) .lt. TOL_Small )) then
     findigs = igf
     exit
    endif
   enddo
end function

end program chi0expand
