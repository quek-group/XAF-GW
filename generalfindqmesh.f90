program findqmesh

implicit none
  integer, parameter :: DP = kind(1.0d0)
  integer, parameter :: DPC = kind((1.0d0,1.0d0))
integer :: n1,n2,i,m1,m2,mm,mden,p,q,nx1,nx2,mx1,mx2,my1,my2
integer :: components(3),qkv(3)
real(DP) :: xm1,xm2,scqx,scqy,ucqx,ucqy
real(DP) :: qpt(3)
real(DP), parameter :: TOL_Small = 1.0d-6

  open(14,file='qlist',status='replace')
  open(15,file='qlistshift',status='replace')
read(*,*) p, q, n1,n2,m1,m2,scqx,scqy
  mm = n1*m2 - n2*m1
if(mm .lt. 0) then
  mm = -mm
endif
mden = mm*(p)-1
ucqx = (m2*scqx-n2*scqy)/mm 
ucqy = (-m1*scqx+n1*scqy)/mm
write(*,*)ucqx,ucqy 
i=0
mx1=MAX(i,n1 ,n2 ,n1+n2)
mx2=MAX(i,m1 ,m2 ,m1+m2)
my1=MIN(i,n1 ,n2 ,n1+n2)
my2=MIN(i,m1 ,m2 ,m1+m2)
qpt(3) = 0.0d0
do nx1 = my1*p, mx1*p-1
 do nx2 = my2*q,mx2*q-1
 xm1 = (m2*(1.0d0*nx1/p) - n2*(1.0d0*nx2/q))/mm
 xm2 =(-m1*(1.0d0*nx1/p) + n1*(1.0d0*nx2/q))/mm
  if(xm1 .ge. 0.0 .and. xm2 .ge. 0.0)then
   if(xm1 .lt. 1.0 .and. xm2 .lt. 1.0)then
    qpt(1) = xm1
    qpt(2) = xm2
    write(14,'(3(f13.9))') qpt(:)
    qpt(1) = xm1+ucqx
    qpt(2) = xm2+ucqy
    write(15,'(3(f13.9))') qpt(:)
   endif
  endif
 enddo
enddo
  close(14)
  close(15)
end program findqmesh
