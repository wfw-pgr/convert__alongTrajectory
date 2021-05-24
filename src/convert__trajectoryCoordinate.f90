! ====================================================== !
! === convert__trajectoryCoordinate                  === !
! ====================================================== !
subroutine convert__trajectoryCoordinate( trajectory, coord_xy, coord_sn, ntrajectory, ncoord )
  implicit none
  integer         , intent(in)  :: ntrajectory, ncoord
  double precision, intent(in)  :: coord_xy(3,ncoord), trajectory(3,ntrajectory)
  double precision, intent(out) :: coord_sn(3,ncoord)
  integer                       :: ic
  double precision              :: sval, nval, xtgt(3), xcrs(3)
  integer         , parameter   :: s_=1, n_=2

  xcrs(:) = 0.d0
  xtgt(:) = 0.d0
  do ic=1, ncoord
     xtgt(:)      = coord_xy(:,ic)
     call find__perpendicularLine( xtgt, trajectory, sval, nval, xcrs, ntrajectory )
     coord_sn(s_,ic) = sval
     coord_sn(n_,ic) = nval
  enddo

  return
  
contains

  ! ====================================================== !
  ! === find root of the perpendicular line            === !
  ! ====================================================== !
  subroutine find__perpendicularLine( xtgt, trajectory, sval, nval, xcrs, ntrajectory )
    implicit none
    integer         , intent(in)  :: ntrajectory
    double precision, intent(in)  :: trajectory(3,ntrajectory), xtgt(3)
    double precision, intent(out) :: sval, nval, xcrs(3)
    integer                       :: i, j, ik, ic, ic_old, iter, ipts(4)
    double precision              :: closest, closest_p, closest_n, dist
    double precision              :: Amat(4,4), xvec(4), yvec(4), xcoef(4), ycoef(4)
    double precision              :: xfit, yfit, xprime, yprime
    double precision              :: sa, sx, sy, nx, ny, rx, ry, ra, dx, dy, da
    double precision              :: sdotd, rdotn, del_x, del_y, del_s, costh
    double precision, allocatable :: ds(:,:), slength(:)
    integer         , parameter   :: x_=1, y_=2, z_=3, a_=4
    integer         , parameter   :: nSize        = 4
    integer         , parameter   :: iterMax      = 10001
    double precision, parameter   :: alpha_picard = 0.5d0
    double precision, parameter   :: eps          = 1.d-8
    logical         , parameter   :: Flag__check  = .false.


    ! ------------------------------------------------------ !
    ! --- [1] length calculation                         --- !
    ! ------------------------------------------------------ !
    allocate( ds(4,ntrajectory), slength(ntrajectory) )
    slength(:)   = 0.d0
    do ik=1, ntrajectory-1
       ds(x_,ik)      = trajectory(x_,ik+1) - trajectory(x_,ik)
       ds(y_,ik)      = trajectory(y_,ik+1) - trajectory(y_,ik)
       ds(z_,ik)      = trajectory(z_,ik+1) - trajectory(z_,ik)
       ds(a_,ik)      = sqrt( ds(x_,ik)**2 + ds(y_,ik)**2 )
       slength(ik+1)  = slength(ik) + ds(a_,ik)
    enddo
    ds(ntrajectory,:) = ds(ntrajectory-1,:)

    ! ------------------------------------------------------ !
    ! --- [2] closest point                              --- !
    ! ------------------------------------------------------ !
    !  -- [2-1] closest point :: all search              --  !
    closest = sqrt( ( trajectory(x_,1)-xtgt(x_) )**2 + ( trajectory(y_,1)-xtgt(y_) )**2 )
    ic      = 1
    ic_old  = 0
    do ik=2, ntrajectory
       dist = sqrt( ( trajectory(x_,ik)-xtgt(x_) )**2 + ( trajectory(y_,ik)-xtgt(y_) )**2 )
       if ( dist.lt.closest ) then
          closest = dist
          ic      = ik
       endif
    enddo
    if ( ( ic.lt.2 ).or.( ic.gt.ntrajectory-1 ) ) then
       stop "out of range!"
    endif
    !  -- [2-2] find next point                          --  !
    closest_p = sqrt( ( trajectory(x_,ic+1)-xtgt(x_) )**2 + ( trajectory(y_,ic+1)-xtgt(y_) )**2 )
    closest_n = sqrt( ( trajectory(x_,ic-1)-xtgt(x_) )**2 + ( trajectory(y_,ic-1)-xtgt(y_) )**2 )
    !  -- [2-3] save ipts to be used in interpolation    --  !
    if ( closest_p.lt.closest_n ) then
       do i=1, 4
          ipts(i) = ic + ( i-2 )
       enddo
    else
       do i=1, 4
          ipts(i) = ic + ( i-3 )
       enddo
    endif
    sval = slength( ic )

    ! ------------------------------------------------------ !
    ! --- [3] Main Loop of finding root                  --- !
    ! ------------------------------------------------------ !
    do iter=1, iterMax

       ! ------------------------------------------------------ !
       ! --- [4] fitting / interpolation                    --- !
       ! ------------------------------------------------------ !
       !  -- [4-1] fitting                                  --  !
       if ( ic.ne.ic_old ) then
          do j=1, nSize
             do i=1, nSize
                Amat(i,j) = dble( slength(ipts(i)) )**dble( nSize-j )
             enddo
          enddo
          do i=1, nSize
             xvec(i) = trajectory( x_, ipts(i) )
             yvec(i) = trajectory( y_, ipts(i) )
          enddo
          call gaussElimin( Amat, xcoef, xvec, nSize )
          call gaussElimin( Amat, ycoef, yvec, nSize )
       endif
       !  -- [4-2] interpolation at sval                    --  !
       xfit   = xcoef(1)*sval**3 + xcoef(2)*sval**2 + xcoef(3)*sval + xcoef(4)
       yfit   = ycoef(1)*sval**3 + ycoef(2)*sval**2 + ycoef(3)*sval + ycoef(4)
       xprime = 3.d0*xcoef(1)*sval**2 + 2.d0*xcoef(2)*sval + xcoef(3)
       yprime = 3.d0*ycoef(1)*sval**2 + 2.d0*ycoef(2)*sval + ycoef(3)

       ! ------------------------------------------------------ !
       ! --- [5] Gram-Schmidt Orthogonalization             --- !
       ! ------------------------------------------------------ !
       !  -- [5-1] tangent vector                           --  !
       sa     = sqrt( xprime**2 + yprime**2 )
       sx     = xprime / sa
       sy     = yprime / sa
       !  -- [5-2] normal vector                            --  !
       nx     = + sy
       ny     = - sx
       !  -- [5-3] perpendicular line :: xtgt(x_) - xfit          --  !
       rx     = xtgt(x_) - xfit
       ry     = xtgt(y_) - yfit
       ra     = sqrt( rx**2 + ry**2 )
       !  -- [5-4] delta vector ( translate vector )        --  !
       rdotn  = rx*nx + ry*ny
       dx     = rx - rdotn * nx
       dy     = ry - rdotn * ny
       da     = sqrt( dx**2 + dy**2 )
       !  -- [5-5] projection on tangent vector :: delta    --  !
       sdotd  = sx*dx + sy*dy
       del_x  = sdotd * sx
       del_y  = sdotd * sy
       del_s  = sign( 1.d0, sdotd ) * sqrt( del_x**2 + del_y**2 )

       ! ------------------------------------------------------ !
       ! --- [6] update sval & ipts                         --- !
       ! ------------------------------------------------------ !
       !  -- [6-1] update sval & nval                       --  !
       sval     = sval + alpha_picard * del_s
       nval     = rdotn
       costh    = rdotn / ra
       xcrs(x_) = xfit
       xcrs(y_) = yfit
       xcrs(z_) = 0.d0
       !  -- [6-2] update ic & ipts                         --  !
       ic_old = ic
       do i=1, ntrajectory-1
          if ( ( sval.ge.slength(i) ).and.( sval.lt.slength(i+1) ) ) then
             ic = i
             exit
          endif
       enddo
       do j=1, 4
          ipts(j) = ic + (j-2)
       enddo
       !  -- [6-3] check & exit from loop                    --  !
       if ( Flag__check ) then
          if ( iter.eq.1 ) write(6,"(a)") "# iter  ic  DeltaStep(da)  sval  nval  cos(th)"
          write(6,"(2(i8,1x),4(f15.8,1x))") iter, ic, da, sval, nval, costh
       endif
       if ( da.lt.eps ) then
          if ( Flag__check ) write(6,*) "END"
          exit
       endif

    enddo

    return
  end subroutine find__perpendicularLine



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

end subroutine convert__trajectoryCoordinate
