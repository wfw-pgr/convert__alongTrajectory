program main
  implicit none
  integer         , parameter   :: lun  = 50
  integer         , parameter   :: cLen = 300
  character(cLen) , parameter   :: inpFile = "dat/trajectory_input.dat"
  character(cLen) , parameter   :: outFile = "dat/output.dat"
  integer                       :: i, j, k, ik, LI, LJ, LK, ntrajectory, ncmp
  integer                       :: i1, i2, i3, i4, icls
  double precision              :: closest, closest_p, closest_n, dist, total_length
  double precision, allocatable :: trajectory(:,:)
  double precision, allocatable :: ds(:,:), slength(:), sbar(:)
  character(cLen)               :: cmt
  double precision, parameter   :: xt =   0.6d0
  double precision, parameter   :: yt = - 0.7d0
  integer         , parameter   :: x_=1, y_=2, z_=3, a_=4
  integer         , parameter   :: nSize = 4
  double precision              :: Amat(4,4), xvec(4), yvec(4), xcoef(4), ycoef(4)
  
  
  ! ------------------------------------------------------ !
  ! --- [1] load trajectory                            --- !
  ! ------------------------------------------------------ !
  open(lun,file=trim(inpFile),status="old")
  read(lun,*)
  read(lun,*) cmt, ntrajectory, ncmp
  read(lun,*) 
  allocate( trajectory(ncmp,ntrajectory) )
  do ik=1, ntrajectory
     read(lun,*) trajectory(:,ik)
  enddo
  close(lun)
  
  ! ------------------------------------------------------ !
  ! --- [2] length calculation                         --- !
  ! ------------------------------------------------------ !
  allocate( ds(ntrajectory,4), slength(ntrajectory), sbar(ntrajectory)  )
  slength(:)   = 0.d0
  do ik=1, ntrajectory-1
     ds(ik,x_)     = trajectory(ik+1,x_) - trajectory(ik,x_)
     ds(ik,y_)     = trajectory(ik+1,y_) - trajectory(ik,y_)
     ds(ik,z_)     = trajectory(ik+1,z_) - trajectory(ik,z_)
     ds(ik,a_)     = sqrt( ds(ik,x_)**2 + ds(ik,y_)**2 )
     slength(ik+1) = slength(ik) + ds(ik,a_)
  enddo
  ds(ntrajectory,:) = ds(ntrajectory-1,:)
  total_length = slength(ntrajectory)
  do ik=1, ntrajectory
     sbar(ik) = dble(ik)
  enddo
  
  ! ------------------------------------------------------ !
  ! --- [3] closest point                              --- !
  ! ------------------------------------------------------ !
  closest = sqrt( ( trajectory(1,x_)-xt )**2 + ( trajectory(1,y_)-yt )**2 )
  icls    = 1
  do ik=2, ntrajectory
     dist = sqrt( ( trajectory(ik,x_)-xt )**2 + ( trajectory(ik,y_)-yt )**2 )
     if ( dist.lt.closest ) then
        dist = closest
        icls = ik
     endif
  enddo
  if ( ( icls.lt.2 ).or.( icls.gt.ntrajectory-1 ) ) then
     stop "out of range!"
  endif
  
  closest_p = sqrt( ( trajectory(icls+1,x_)-xt )**2 + ( trajectory(icls+1,y_)-yt )**2 )
  closest_n = sqrt( ( trajectory(icls-1,x_)-xt )**2 + ( trajectory(icls-1,y_)-yt )**2 )

  if ( closest_p.lt.closest_n ) then
     i1 = icls-1
     i2 = icls
     i3 = icls+1
     i4 = icls+2
  else
     i1 = icls-2
     i2 = icls-1
     i3 = icls
     i4 = icls+1
  endif

  
  ! ------------------------------------------------------ !
  ! --- [3] interpolation                              --- !
  ! ------------------------------------------------------ !
  do j=1, nSize
     do i=1, nSize
        Amat(i,j) = sbar( i1+(i-1) )**dble( nSize-(j-1) )
     enddo
  enddo
  do i=1, nSize
     xvec(i) = trajectory( i1+(i-1), x_ )
  enddo
  do i=1, nSize
     yvec(i) = trajectory( i1+(i-1), y_ )
  enddo
  call gaussElimin( Amat, xcoef, xvec, nSize )
  call gaussElimin( Amat, ycoef, yvec, nSize )
  

contains


  ! ========================================================== !
  ! === Gauss Elimination Solver                           === !
  ! ========================================================== !
  subroutine gaussElimin( Amat, xvec, bvec, nSize )
    implicit none
    integer         , intent(in)  :: nSize
    double precision, intent(in)  :: Amat(nSize,nSize)
    double precision, intent(in)  :: bvec(nSize)
    double precision, intent(out) :: xvec(nSize)
    integer                       :: i, j, k, ipivot
    double precision              :: Dinv, buff, vpivot
    double precision, parameter   :: eps = 1.d-10
    double precision, allocatable :: Umat(:,:), vvec(:,:)

    ! ----------------------------------------- !
    ! --- [1] Preparation                   --- !
    ! ----------------------------------------- !
    !  -- [1-1] allocate Umat               --  !
    allocate( Umat(nSize,nSize) )
    !  -- [1-2] Copy AMat & bvec            --  !
    Umat(:,:) = Amat(:,:)
    xvec(:)   = bvec(:)
    
    ! ----------------------------------------- !
    ! --- [2] Forward Ellimination          --- !
    ! ----------------------------------------- !
    do k=1, nSize

       !  -- [2-1] Pivoting                 --  !
       vpivot = abs( Umat(k,k) )
       ipivot = k
       do j=k+1, nSize
          if ( abs( Umat(j,k) ).gt.vpivot ) then
             vpivot = abs( Umat(j,k) )
             ipivot = j
          endif
       end do
       if ( ipivot.ne.k ) then
          do j=k, nSize
             buff           = Umat(ipivot,j)
             Umat(ipivot,j) = Umat(k     ,j)
             Umat(k     ,j) = buff
          enddo
          buff         = xvec(ipivot)
          xvec(ipivot) = xvec(k)
          xvec(k)      = buff
       end if
       if ( abs( Umat(k,k) ).lt.eps ) then
          write(6,*) '[gaussElimin] Amat :: Singular Matrix :: No Solution End :: @ k= ', k
          stop
       endif
       !  -- [2-2] Diagonal Component       --  !
       Dinv      = 1.d0 / Umat(k,k)
       Umat(k,k) = 1.d0

       !  -- [2-3] Non-Diagonal Component   --  !
       if ( k.eq.nSize ) then
          ! -- [    Last Row :: k == nSize ] -- !
          xvec(k) = Dinv * xvec(k)
       else
          ! -- [Not Last Row :: k != nSize ] -- !
          !  - Division    -  !
          Umat(k,k+1:nSize) = Dinv * Umat(k,k+1:nSize)
          xvec(k)           = Dinv * xvec(k)
          !  - subtraction -  !
          do j=k+1,nSize
             Umat(j,k+1:nSize) = Umat(j,k+1:nSize) - Umat(j,k) * Umat(k,k+1:nSize)
             xvec(j)           = xvec(j)           - Umat(j,k) * xvec(k)
             Umat(j,k)         = 0.d0
          enddo
       endif
       
    end do

    ! ----------------------------------------- !
    ! --- [3] Backward Substituition        --- !
    ! ----------------------------------------- !
    do k=nSize-1, 1, -1
       do i=nSize, k+1, -1
          xvec(k) = xvec(k) - Umat(k,i)*xvec(i)
       enddo
    enddo

    return
  end subroutine gaussElimin

  
end program main
