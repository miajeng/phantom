!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! Glue grains of a settled DEM packing into boulders.
!
! Optional stage after the crop in the settle-crop-settle packing method.
! The body is already a dense, settled packing of equal grains, so a boulder
! is just a set of grains given the same clump ID (part%iclump); dem.f90
! bonds grains that share an ID. Grains left with ID 0 are the dust. Grains
! keep their radius and mass, so a boulder's mass is its grain count and the
! bulk density is unchanged.
!
! Boulders come from two sources, made in this order:
!
!  1. PLACED boulders (iplaced), made first so they get the grains they ask
!     for, and numbered 1, 2, ... in the order made:
!       0 = none
!       1 = contact binary: two touching lobes along the body's long axis.
!           The lobes are spheres with radii in the ratio lobe_ratio, sized
!           so that together they span the body end to end and touch at one
!           point; a lobe is whatever part of the body lies inside its
!           sphere, so the body surface trims it to the real shape.
!       2 = read from a file, one boulder per line (km, degrees):
!             x y z  a b c  [rx ry rz]
!           centre relative to the body's centre of mass, semi-axes, and
!           optional rotations about x, then y, then z. '#' starts a comment.
!
!  2. RANDOM boulders, until f_boulder of the whole body is in them. Each
!     one is seeded on a random free grain with a size R drawn from a
!     truncated power law dN/dR ~ R^-q between rb_min and rb_max (in grain
!     radii), and shaped by ishape:
!       0 = round: grains within R of the seed
!       1 = ellipsoid: semi-axes R, and R times two ratios drawn uniformly
!           in [axis_min,1], with a random orientation
!       2 = grown: grain by grain from the seed, each new grain a random
!           free neighbour of the boulder, up to the grain count a round
!           boulder of radius R would have. With probability jagged the
!           next grain grows off the newest few, which pushes out lobes and
!           arms; jagged=0 gives compact lumps.
!
! A boulder only ever takes FREE grains, so one made next to an earlier one
! is carved by it. Carving can leave a grain in a boulder with no bond to
! the rest, or split it, so every boulder is walked from its seed along
! bonds (surface gap < bond_reach*(R_i+R_j), as in dem.f90) and whatever the
! walk does not reach goes back to the dust.
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: io, part, physcon, prompting, units
!
 use part,      only:iclump
 use prompting, only:prompt
 use io,        only:id,master,fatal
 use units,     only:udist
 use physcon,   only:km,pi

 implicit none
 private
 character(len=*), parameter, public :: moddump_flags = ''
 public :: modify_dump

 ! random boulders
 real    :: f_boulder  = 0.3   ! fraction of the body to put in random boulders
 real    :: rb_min     = 2.5   ! smallest boulder radius, in grain radii
 real    :: rb_max     = 6.0   ! largest boulder radius, in grain radii
 real    :: q_slope    = 3.0   ! power-law slope, dN/dR ~ R^-q
 integer :: ishape     = 0     ! 0 = round, 1 = ellipsoid, 2 = grown
 real    :: axis_min   = 0.5   ! ellipsoids: smallest axis ratio
 real    :: jagged     = 0.5   ! grown: 0 = compact, 1 = lobes and arms
 ! placed boulders
 integer :: iplaced    = 0     ! 0 = none, 1 = contact binary, 2 = from file
 real    :: lobe_ratio = 1.0   ! contact binary: small lobe radius / large lobe radius
 character(len=256) :: placed_file = 'boulders.txt'
 ! common
 real    :: reach      = 0.1   ! bond_reach, must match the .in file
 integer :: nmin_clump = 4     ! smaller random boulders go back to the dust
 integer :: iseed_clump = 4321

 ! neighbour grid and bond list
 real    :: xmin(3),cellsize
 integer :: ncell(3)
 integer, allocatable :: ihead(:,:,:),inext(:),nbstart(:),nblist(:)
 ! work arrays
 integer, allocatable :: members(:),queue(:),cand(:)
 logical, allocatable :: flag(:)

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer, allocatable :: nper(:),seed(:)
 real    :: r_grain,rb,u,a1q,b1q,ax(3),rot(3,3),xc(3)
 integer :: k,iseed,nclump,nplaced,ntarget,nrandom,nmem,ntail,nattempt,maxattempt
 integer :: ndropped,nsmall,ntot_seed,nplaced_grains

 if (npart < 1) call fatal('moddump_clumps','no particles in dump')

 call prompt('Fraction of the body to glue into RANDOM boulders',f_boulder,0.,1.)
 call prompt('Smallest random boulder radius, in grain radii',rb_min,1.)
 call prompt('Largest random boulder radius, in grain radii',rb_max,rb_min)
 call prompt('Power-law slope q of the size distribution, dN/dR ~ R^-q',q_slope,0.)
 call prompt('Random boulder shape (0=round, 1=ellipsoid, 2=grown)',ishape,0,2)
 if (ishape == 1) call prompt('Smallest axis ratio of the ellipsoids',axis_min,0.1,1.)
 if (ishape == 2) call prompt('Jaggedness of grown boulders (0=compact, 1=lobes and arms)',jagged,0.,1.)
 call prompt('Placed boulders (0=none, 1=contact binary, 2=read from file)',iplaced,0,2)
 if (iplaced == 1) call prompt('Contact binary: small lobe radius / large lobe radius',lobe_ratio,0.1,1.)
 if (iplaced == 2) call prompt('File listing the placed boulders (x y z a b c [rx ry rz], km and deg)',placed_file)
 call prompt('bond_reach (must match the .in file)',reach,0.)
 call prompt('Smallest random boulder kept, in grains',nmin_clump,2)
 call prompt('Random seed',iseed_clump)

 call random_seed(size=ntot_seed)
 allocate(seed(ntot_seed))
 do k=1,ntot_seed
    seed(k) = iseed_clump + 37*(k-1)
 enddo
 call random_seed(put=seed)
 deallocate(seed)
 !
 ! every DEM grain carries its radius in h, and all grains are equal
 !
 r_grain = maxval(xyzh(4,1:npart))
 call build_bonds(npart,xyzh,r_grain)
 !
 ! start from all dust: any clumps already in the dump are replaced
 !
 iclump(1:npart) = 0
 allocate(members(npart),queue(npart),flag(npart),nper(npart))
 flag(:) = .false.
 nclump = 0

 if (id==master) print "(/,a)",' --- gluing grains into boulders ---'
 !
 ! 1. placed boulders
 !
 select case(iplaced)
 case(1)
    call place_contact_binary(npart,xyzh,nclump,nper)
 case(2)
    call place_from_file(npart,xyzh,nclump,nper)
 end select
 nplaced = nclump
 nplaced_grains = 0
 if (nplaced > 0) nplaced_grains = sum(nper(1:nplaced))
 !
 ! 2. random boulders
 !
 ntarget    = nint(f_boulder*real(npart))
 nrandom    = 0
 nsmall     = 0
 ndropped   = 0
 nattempt   = 0
 maxattempt = 100*npart
 a1q = rb_min**(1.-q_slope)
 b1q = rb_max**(1.-q_slope)

 do while (nrandom < ntarget .and. nattempt < maxattempt)
    nattempt = nattempt + 1
    !
    ! boulder radius from the truncated power law (inverse transform)
    !
    call random_number(u)
    if (abs(q_slope - 1.) < 1.e-6) then
       rb = rb_min*(rb_max/rb_min)**u
    else
       rb = (a1q + u*(b1q - a1q))**(1./(1.-q_slope))
    endif
    !
    ! seed on a random free grain
    !
    call random_number(u)
    iseed = min(npart,1 + int(u*npart))
    if (iclump(iseed) /= 0) cycle

    select case(ishape)
    case(2)
       call grow_clump(iseed,max(nmin_clump,nint(0.64*rb**3)),ntail)
    case default
       xc = xyzh(1:3,iseed)
       ax = rb*r_grain
       rot = reshape((/1.,0.,0.,0.,1.,0.,0.,0.,1./),(/3,3/))
       if (ishape == 1) then
          call random_number(u)
          ax(2) = ax(1)*(axis_min + (1. - axis_min)*u)
          call random_number(u)
          ax(3) = ax(1)*(axis_min + (1. - axis_min)*u)
          call random_rotation(rot)
       endif
       call gather_ellipsoid(npart,xyzh,xc,ax,rot,nmem)
       call walk_bonds(iseed,nmem,ntail)
       ndropped = ndropped + (nmem - ntail)
    end select
    if (ntail < nmin_clump) then
       nsmall = nsmall + 1
       cycle
    endif

    nclump = nclump + 1
    iclump(queue(1:ntail)) = nclump
    nper(nclump) = ntail
    nrandom = nrandom + ntail
 enddo

 if (id==master) then
    print "(a,i9)",        ' grains               = ',npart
    if (nplaced > 0) print "(a,i9,a,i9,a,f6.3,a)",' placed boulders      = ',nplaced,'  holding ',&
                            nplaced_grains,' grains (',real(nplaced_grains)/real(npart),' of the body)'
    print "(a,i9,a,a)",    ' random boulders      = ',nclump-nplaced,'  shape: ', &
                            trim(merge('round    ',merge('ellipsoid','grown    ',ishape==1),ishape==0))
    print "(a,i9,a,f6.3,a,f6.3,a)",' grains in random     = ',nrandom,'  (',real(nrandom)/real(npart),&
                                   ' of the body, target ',f_boulder,')'
    if (nclump > nplaced) print "(a,i6,a,i6,a,f8.1)",' grains per random boulder = ', &
       minval(nper(nplaced+1:nclump)),' to ',maxval(nper(nplaced+1:nclump)),', mean ', &
       real(nrandom)/real(nclump-nplaced)
    print "(a,f6.3)",      ' dust (free grains)   = ',real(count(iclump(1:npart)==0))/real(npart)
    print "(a,i9)",        ' grains cut off by carving, back to dust = ',ndropped
    print "(a,i9)",        ' random boulders under the minimum size, skipped = ',nsmall
    if (nrandom < ntarget) print "(a)",' *** WARNING: ran out of free grains or attempts before the random target ***'
    print "(a)",' now set kb_cgs > 0 (and the same bond_reach) in the .in file and settle again'
 endif

 deallocate(members,queue,flag,nper,cand)
 deallocate(ihead,inext,nbstart,nblist)

end subroutine modify_dump

!----------------------------------------------------------------
!+
!  contact binary: two touching lobes along the body's long axis
!+
!----------------------------------------------------------------
subroutine place_contact_binary(npart,xyzh,nclump,nper)
 integer, intent(in)    :: npart
 real,    intent(in)    :: xyzh(:,:)
 integer, intent(inout) :: nclump,nper(:)
 real    :: com(3),cov(3,3),e(3),enew(3),p,pmin,pmax,r1,r2,c(3,2),rad(2),ax(3),rot(3,3),dx(3)
 integer :: i,it,ilobe,nmem,ntail,iclose

 com = sum(xyzh(1:3,1:npart),dim=2)/real(npart)
 cov = 0.
 do i=1,npart
    dx = xyzh(1:3,i) - com
    cov = cov + spread(dx,2,3)*spread(dx,1,3)
 enddo
 !
 ! long axis = largest eigenvector of the position covariance (power iteration)
 !
 e = (/1.,0.3,0.1/)
 e = e/norm2(e)
 do it=1,200
    enew = matmul(cov,e)
    enew = enew/norm2(enew)
    if (norm2(enew - e) < 1.e-12) exit
    e = enew
 enddo
 e = enew
 pmin = huge(pmin)
 pmax = -huge(pmax)
 do i=1,npart
    p = dot_product(xyzh(1:3,i) - com,e)
    pmin = min(pmin,p)
    pmax = max(pmax,p)
 enddo
 !
 ! lobe spheres span the body end to end and touch at one point
 !
 r1 = (pmax - pmin)/(2.*(1. + lobe_ratio))
 r2 = lobe_ratio*r1
 c(:,1) = com + (pmin + r1)*e
 c(:,2) = com + (pmax - r2)*e
 rad = (/r1,r2/)
 rot = reshape((/1.,0.,0.,0.,1.,0.,0.,0.,1./),(/3,3/))

 if (id==master) then
    print "(a)",' contact binary along the long axis'
    print "(a,3(1pg11.3))",' long axis direction  = ',e
    print "(a)",' same lobes as a boulders file (x y z a b c, km, relative to the centre of mass):'
 endif
 do ilobe=1,2
    ax = rad(ilobe)
    call gather_ellipsoid(npart,xyzh,c(:,ilobe),ax,rot,nmem)
    if (nmem < 1) call fatal('moddump_clumps','contact binary lobe holds no grains')
    iclose = closest(npart,xyzh,c(:,ilobe),nmem)
    call walk_bonds(iclose,nmem,ntail)
    nclump = nclump + 1
    iclump(queue(1:ntail)) = nclump
    nper(nclump) = ntail
    if (id==master) then
       print "(a,6(1pg11.3))",'   ',(c(:,ilobe) - com)*udist/km,ax*udist/km
       print "(a,i1,a,i8,a,1pg10.3,a)",' lobe ',ilobe,': ',ntail,' grains, radius ',rad(ilobe)*udist/km,' km'
    endif
 enddo

end subroutine place_contact_binary

!----------------------------------------------------------------
!+
!  placed boulders read from a file: x y z a b c [rx ry rz]
!+
!----------------------------------------------------------------
subroutine place_from_file(npart,xyzh,nclump,nper)
 integer, intent(in)    :: npart
 real,    intent(in)    :: xyzh(:,:)
 integer, intent(inout) :: nclump,nper(:)
 character(len=512) :: line
 real    :: v(9),com(3),ang(3),rot(3,3),rx(3,3),ry(3,3),rz(3,3),cs(3),sn(3)
 integer :: iu,ios,nmem,ntail,iclose,nline,icomment

 com = sum(xyzh(1:3,1:npart),dim=2)/real(npart)
 open(newunit=iu,file=trim(placed_file),status='old',action='read',iostat=ios)
 if (ios /= 0) call fatal('moddump_clumps','could not open '//trim(placed_file))
 nline = 0
 do
    read(iu,'(a)',iostat=ios) line
    if (ios /= 0) exit
    nline = nline + 1
    icomment = index(line,'#')
    if (icomment > 0) line = line(1:icomment-1)
    if (len_trim(line) == 0) cycle
    v = 0.
    read(line,*,iostat=ios) v(1:9)
    if (ios /= 0) then
       v(7:9) = 0.
       read(line,*,iostat=ios) v(1:6)
       if (ios /= 0) call fatal('moddump_clumps','cannot read line of '//trim(placed_file)//': '//trim(line))
    endif
    if (any(v(4:6) <= 0.)) call fatal('moddump_clumps','semi-axes must be > 0 in '//trim(placed_file))
    !
    ! rotate about x, then y, then z: R = Rz Ry Rx
    !
    ang = v(7:9)*pi/180.
    cs = cos(ang)
    sn = sin(ang)
    rx = reshape((/1.,0.,0., 0.,cs(1),sn(1), 0.,-sn(1),cs(1)/),(/3,3/))
    ry = reshape((/cs(2),0.,-sn(2), 0.,1.,0., sn(2),0.,cs(2)/),(/3,3/))
    rz = reshape((/cs(3),sn(3),0., -sn(3),cs(3),0., 0.,0.,1./),(/3,3/))
    rot = matmul(rz,matmul(ry,rx))

    call gather_ellipsoid(npart,xyzh,com + v(1:3)*km/udist,v(4:6)*km/udist,rot,nmem)
    if (nmem < 1) then
       if (id==master) print "(a,i4,a)",' *** WARNING: placed boulder on line ',nline,' holds no free grains: skipped ***'
       cycle
    endif
    iclose = closest(npart,xyzh,com + v(1:3)*km/udist,nmem)
    call walk_bonds(iclose,nmem,ntail)
    nclump = nclump + 1
    iclump(queue(1:ntail)) = nclump
    nper(nclump) = ntail
    if (id==master) print "(a,i4,a,i8,a,i8,a)",' placed boulder ',nclump,': ',ntail,' grains (',nmem-ntail,&
                          ' unconnected, left as dust)'
 enddo
 close(iu)
 if (nclump == 0 .and. id==master) print "(a)",' *** WARNING: no placed boulders read from '//trim(placed_file)//' ***'

end subroutine place_from_file

!----------------------------------------------------------------
!+
!  grow a boulder grain by grain from iseed, up to ngoal grains.
!  Result in queue(1:ntail).
!+
!----------------------------------------------------------------
subroutine grow_clump(iseed,ngoal,ntail)
 integer, intent(in)  :: iseed,ngoal
 integer, intent(out) :: ntail
 integer :: k,m,j,g,nfree,nfail
 real    :: u

 ntail = 1
 queue(1) = iseed
 flag(iseed) = .true.
 nfail = 0
 do while (ntail < ngoal .and. nfail < 50)
    !
    ! grow off the newest few grains (arms) or anywhere (compact)
    !
    call random_number(u)
    if (u < jagged) then
       call random_number(u)
       k = ntail - int(u*min(4,ntail))
    else
       call random_number(u)
       k = 1 + int(u*ntail)
    endif
    g = queue(max(1,min(k,ntail)))
    nfree = 0
    do m=nbstart(g),nbstart(g+1)-1
       j = nblist(m)
       if (iclump(j) == 0 .and. .not.flag(j)) then
          nfree = nfree + 1
          cand(nfree) = j
       endif
    enddo
    if (nfree == 0) then
       nfail = nfail + 1
       cycle
    endif
    call random_number(u)
    j = cand(min(nfree,1 + int(u*nfree)))
    ntail = ntail + 1
    queue(ntail) = j
    flag(j) = .true.
    nfail = 0
 enddo
 flag(queue(1:ntail)) = .false.

end subroutine grow_clump

!----------------------------------------------------------------
!+
!  free grains with centres inside an ellipsoid (centre xc, semi-axes
!  ax, rotation rot from body frame to space). Result in members(1:nmem).
!+
!----------------------------------------------------------------
subroutine gather_ellipsoid(npart,xyzh,xc,ax,rot,nmem)
 integer, intent(in)  :: npart
 real,    intent(in)  :: xyzh(:,:),xc(3),ax(3),rot(3,3)
 integer, intent(out) :: nmem
 real    :: dx(3),db(3),rmax
 integer :: lo(3),hi(3),ix,iy,iz,j

 rmax = maxval(ax)
 lo = max(1,    int((xc - rmax - xmin)/cellsize) + 1)
 hi = min(ncell,int((xc + rmax - xmin)/cellsize) + 1)
 nmem = 0
 do iz=lo(3),hi(3)
    do iy=lo(2),hi(2)
       do ix=lo(1),hi(1)
          j = ihead(ix,iy,iz)
          do while (j > 0)
             if (iclump(j) == 0) then
                dx = xyzh(1:3,j) - xc
                db = matmul(transpose(rot),dx)
                if (sum((db/ax)**2) <= 1.) then
                   nmem = nmem + 1
                   members(nmem) = j
                endif
             endif
             j = inext(j)
          enddo
       enddo
    enddo
 enddo
 if (nmem > npart) call fatal('moddump_clumps','gather overflow')

end subroutine gather_ellipsoid

!----------------------------------------------------------------
!+
!  keep only the members reachable from istart along bonds.
!  Result in queue(1:ntail).
!+
!----------------------------------------------------------------
subroutine walk_bonds(istart,nmem,ntail)
 integer, intent(in)  :: istart,nmem
 integer, intent(out) :: ntail
 integer :: nhead,j,k,m

 ! flag marks candidates not yet reached
 flag(members(1:nmem)) = .true.
 flag(istart) = .false.
 queue(1) = istart
 nhead = 1
 ntail = 1
 do while (nhead <= ntail)
    j = queue(nhead)
    nhead = nhead + 1
    do m=nbstart(j),nbstart(j+1)-1
       k = nblist(m)
       if (flag(k)) then
          flag(k) = .false.
          ntail = ntail + 1
          queue(ntail) = k
       endif
    enddo
 enddo
 flag(members(1:nmem)) = .false.

end subroutine walk_bonds

!----------------------------------------------------------------
!+
!  member grain closest to a point
!+
!----------------------------------------------------------------
integer function closest(npart,xyzh,xc,nmem)
 integer, intent(in) :: npart,nmem
 real,    intent(in) :: xyzh(:,:),xc(3)
 real    :: d2,d2min
 integer :: k

 closest = members(1)
 d2min = huge(d2min)
 do k=1,nmem
    d2 = sum((xyzh(1:3,members(k)) - xc)**2)
    if (d2 < d2min) then
       d2min = d2
       closest = members(k)
    endif
 enddo

end function closest

!----------------------------------------------------------------
!+
!  uniformly random rotation matrix (Shoemake 1992 quaternion)
!+
!----------------------------------------------------------------
subroutine random_rotation(rot)
 real, intent(out) :: rot(3,3)
 real :: u(3),x,y,z,w

 call random_number(u)
 x = sqrt(1.-u(1))*sin(2.*pi*u(2))
 y = sqrt(1.-u(1))*cos(2.*pi*u(2))
 z = sqrt(u(1))*sin(2.*pi*u(3))
 w = sqrt(u(1))*cos(2.*pi*u(3))
 rot(1,:) = (/1.-2.*(y*y+z*z), 2.*(x*y-z*w),    2.*(x*z+y*w)/)
 rot(2,:) = (/2.*(x*y+z*w),    1.-2.*(x*x+z*z), 2.*(y*z-x*w)/)
 rot(3,:) = (/2.*(x*z-y*w),    2.*(y*z+x*w),    1.-2.*(x*x+y*y)/)

end subroutine random_rotation

!----------------------------------------------------------------
!+
!  cell grid, and the list of bonded neighbours of every grain
!  (compressed rows: grain i's bonds are nblist(nbstart(i):nbstart(i+1)-1))
!+
!----------------------------------------------------------------
subroutine build_bonds(npart,xyzh,r_grain)
 integer, intent(in) :: npart
 real,    intent(in) :: xyzh(:,:),r_grain
 integer :: i,j,ic(3),ipass,nb,ix,iy,iz,lo(3),hi(3),maxdeg

 ! any bonded pair is within one cell of each other
 cellsize = 2.*r_grain*(1. + reach)*1.0001
 xmin  = minval(xyzh(1:3,1:npart),dim=2)
 ncell = int((maxval(xyzh(1:3,1:npart),dim=2) - xmin)/cellsize) + 1
 allocate(ihead(ncell(1),ncell(2),ncell(3)),inext(npart),nbstart(npart+1))
 ihead = 0
 do i=1,npart
    ic = min(ncell,int((xyzh(1:3,i) - xmin)/cellsize) + 1)
    inext(i) = ihead(ic(1),ic(2),ic(3))
    ihead(ic(1),ic(2),ic(3)) = i
 enddo

 do ipass=1,2
    nb = 0
    do i=1,npart
       if (ipass==1) then
          nbstart(i) = nb + 1
       endif
       ic = min(ncell,int((xyzh(1:3,i) - xmin)/cellsize) + 1)
       lo = max(1,ic-1)
       hi = min(ncell,ic+1)
       do iz=lo(3),hi(3)
          do iy=lo(2),hi(2)
             do ix=lo(1),hi(1)
                j = ihead(ix,iy,iz)
                do while (j > 0)
                   if (j /= i) then
                      if (bonded(xyzh(:,i),xyzh(:,j))) then
                         nb = nb + 1
                         if (ipass==2) nblist(nb) = j
                      endif
                   endif
                   j = inext(j)
                enddo
             enddo
          enddo
       enddo
    enddo
    if (ipass==1) then
       nbstart(npart+1) = nb + 1
       allocate(nblist(nb))
    endif
 enddo
 maxdeg = 0
 do i=1,npart
    maxdeg = max(maxdeg,nbstart(i+1) - nbstart(i))
 enddo
 allocate(cand(max(1,maxdeg)))
 if (id==master) print "(a,f6.2,a,i3)",' bonded neighbours per grain: mean ', &
    real(nbstart(npart+1)-1)/real(npart),', max ',maxdeg

end subroutine build_bonds

!----------------------------------------------------------------
!+
!  same bond test as get_ssdem_force: surface gap within reach
!+
!----------------------------------------------------------------
logical function bonded(xyzhi,xyzhj)
 real, intent(in) :: xyzhi(4),xyzhj(4)
 real :: dx(3),r,sumr

 dx   = xyzhi(1:3) - xyzhj(1:3)
 r    = sqrt(dot_product(dx,dx))
 sumr = xyzhi(4) + xyzhj(4)
 bonded = (r - sumr < reach*sumr)

end function bonded

end module moddump
