!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module analysis
!
! Shape and packing diagnostics for a DEM rubble pile.
!
! Reports the principal semi-axes from the inertia tensor, the axis ratios,
! and the packing fraction, so that the slump of a body away from its
! intended shape can be measured directly: a cohesionless, frictionless
! pile has no angle of repose, so its only equilibrium is a sphere and any
! elongated shape relaxes towards b/a = c/a = 1.
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: centreofmass, part, physcon, units
!
 implicit none
 character(len=20), parameter, public :: analysistype = 'demshape'
 logical, private :: firstcall = .true.

 public :: do_analysis

 private

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use part,         only:nptmass,xyzmh_ptmass,vxyz_ptmass,isdead_or_accreted
 use centreofmass, only:get_centreofmass
 use units,        only:udist,utime
 use physcon,      only:pi,km
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 real    :: xpos(3),vpos(3),dx(3)
 real    :: inert(3,3),eig(3),semi(3),rgrain,rmax,vgrain,phi,vbody
 real    :: funbound
 integer :: i,j,k,nlive
 logical :: iexist
 character(len=200) :: fileout

 !
 ! Centre on the body before taking moments, or the inertia tensor picks up
 ! a parallel-axis term from the offset and the axes come out wrong.
 !
 ! Particles ONLY: the sink arguments are deliberately not passed. In a flyby
 ! the perturber outweighs the body by ~1e14, so including it puts the centre
 ! of mass on the planet, every grain shares one huge offset, the tensor goes
 ! rank-1 and two of the three axes collapse to zero.
 !
 call get_centreofmass(xpos,vpos,npart,xyzh,vxyzu)

 inert  = 0.
 rgrain = 0.
 rmax   = 0.
 nlive  = 0
 do i=1,npart
    if (isdead_or_accreted(xyzh(4,i))) cycle
    nlive = nlive + 1
    dx = xyzh(1:3,i) - xpos
    rmax = max(rmax,sqrt(dot_product(dx,dx)))
    rgrain = max(rgrain,xyzh(4,i))
    do j=1,3
       do k=1,3
          inert(j,k) = inert(j,k) + dx(j)*dx(k)
       enddo
    enddo
 enddo
 if (nlive < 2) return
 inert = inert/real(nlive)
 !
 ! eigenvalues of the 3x3 symmetric second-moment tensor. For a uniform
 ! ellipsoid the second moment about each axis is a^2/5, so the semi-axes
 ! follow as sqrt(5*lambda).
 !
 call eigenvalues_sym3(inert,eig)
 do j=1,3
    semi(j) = sqrt(max(5.*eig(j),0.))
 enddo
 !
 ! packing fraction against the ellipsoid the axes imply
 !
 vgrain = real(nlive)*(4./3.*pi*rgrain**3)
 vbody  = 4./3.*pi*semi(1)*semi(2)*semi(3)
 if (vbody > 0.) then
    phi = vgrain/vbody
 else
    phi = 0.
 endif

 !
 ! energy-based unbound fraction, iterated onto the largest remnant
 !
 call get_unbound_fraction(npart,xyzh,vxyzu,particlemass,funbound)

 fileout = trim(dumpfile(1:index(dumpfile,'_')-1))//'_shape.dat'
 inquire(file=fileout,exist=iexist)
 if (.not.iexist .or. firstcall) then
    firstcall = .false.
    open(iunit,file=fileout,status='replace')
    write(iunit,"('#',11(1x,'[',i2.2,1x,a11,']',2x))") &
          1,'time_s', 2,'npart', 3,'a_m', 4,'b_m', 5,'c_m', &
          6,'b_on_a', 7,'c_on_a', 8,'rmax_m', 9,'rgrain_m', 10,'phi', 11,'f_unbound'
 else
    open(iunit,file=fileout,position='append')
 endif
 write(iunit,'(es18.10,1x,i10,1x,9(es18.10,1x))') time*utime,nlive, &
       semi(1)*udist/km*1000., semi(2)*udist/km*1000., semi(3)*udist/km*1000., &
       semi(2)/semi(1), semi(3)/semi(1), rmax*udist/km*1000., &
       rgrain*udist/km*1000., phi, funbound
 close(iunit)

 print "(/,a,es10.3,a)",' --- DEM shape at t = ',time*utime,' s ---'
 print "(a,i10)",       ' live grains     = ',nlive
 print "(a,3(1x,f9.2))",' semi-axes (m)   = ',semi(1)*udist/km*1000., &
                                              semi(2)*udist/km*1000., &
                                              semi(3)*udist/km*1000.
 print "(a,2(1x,f7.4))",' b/a, c/a        = ',semi(2)/semi(1),semi(3)/semi(1)
 print "(a,f9.2,a)",    ' furthest grain  = ',rmax*udist/km*1000.,' m'
 print "(a,f7.4)",      ' packing frac    = ',phi
 print "(a,f7.4)",      ' unbound frac    = ',funbound

end subroutine do_analysis

!----------------------------------------------------------------
!+
!  Fraction of mass not gravitationally bound to the largest remnant.
!
!  Energy based, not friends-of-friends: the FoF metric is prone to
!  flipping between 0 and exactly 50% on an intact body, which makes it
!  useless near a disruption threshold. Here a grain is unbound if its
!  kinetic energy in the remnant's frame exceeds the potential well of
!  the remnant. Iterated a few times, since removing unbound mass moves
!  the centre of mass and shallows the well.
!
!  Self-gravity of the grains only: any sink (Earth) is excluded, so this
!  measures escape from the body rather than from the encounter. Evaluate
!  it once the perturber is far away or it will not mean what it says.
!+
!----------------------------------------------------------------
subroutine get_unbound_fraction(npart,xyzh,vxyzu,pmass,funbound)
 use part,    only:isdead_or_accreted
 use physcon, only:pi
 integer, intent(in)  :: npart
 real,    intent(in)  :: xyzh(:,:),vxyzu(:,:),pmass
 real,    intent(out) :: funbound
 logical, allocatable :: bound(:)
 real    :: xcm(3),vcm(3),dx(3),dv(3),phii,rij,etot
 integer :: i,j,iter,nb,nb_prev,nlive

 funbound = 0.
 if (npart < 2 .or. pmass <= 0.) return

 allocate(bound(npart))
 nlive = 0
 do i=1,npart
    bound(i) = .not.isdead_or_accreted(xyzh(4,i))
    if (bound(i)) nlive = nlive + 1
 enddo
 if (nlive < 2) then
    deallocate(bound)
    return
 endif
 nb_prev = -1

 do iter=1,5
    nb = count(bound(1:npart))
    if (nb < 2 .or. nb == nb_prev) exit
    nb_prev = nb
    !
    ! centre of mass and mean velocity of the current remnant
    !
    xcm = 0.; vcm = 0.
    do i=1,npart
       if (.not.bound(i)) cycle
       xcm = xcm + xyzh(1:3,i)
       vcm = vcm + vxyzu(1:3,i)
    enddo
    xcm = xcm/real(nb)
    vcm = vcm/real(nb)
    !
    ! energy of every live grain against the remnant's potential
    !
    !$omp parallel do default(none) schedule(guided) &
    !$omp shared(npart,xyzh,vxyzu,bound,xcm,vcm,pmass) &
    !$omp private(i,j,dx,dv,rij,phii,etot)
    do i=1,npart
       if (isdead_or_accreted(xyzh(4,i))) cycle
       phii = 0.
       do j=1,npart
          if (j == i) cycle
          if (.not.bound(j)) cycle
          dx = xyzh(1:3,i) - xyzh(1:3,j)
          rij = sqrt(dot_product(dx,dx))
          if (rij > tiny(rij)) phii = phii - pmass/rij   ! G = 1 in code units
       enddo
       dv = vxyzu(1:3,i) - vcm
       etot = 0.5*dot_product(dv,dv) + phii
       bound(i) = (etot < 0.)
    enddo
    !$omp end parallel do
 enddo

 nb = count(bound(1:npart))
 funbound = 1. - real(nb)/real(nlive)
 deallocate(bound)

end subroutine get_unbound_fraction

!----------------------------------------------------------------
!+
!  eigenvalues of a real symmetric 3x3 matrix, largest first,
!  by the closed-form trigonometric solution of the characteristic
!  cubic (Smith 1961). No external solver needed for 3x3.
!+
!----------------------------------------------------------------
subroutine eigenvalues_sym3(a,eig)
 use physcon, only:pi
 real, intent(in)  :: a(3,3)
 real, intent(out) :: eig(3)
 real :: p1,p2,p,q,r,phi,b(3,3),tmp
 integer :: i,j

 p1 = a(1,2)**2 + a(1,3)**2 + a(2,3)**2
 q  = (a(1,1) + a(2,2) + a(3,3))/3.

 if (p1 <= tiny(p1)) then          ! already diagonal
    eig(1) = a(1,1); eig(2) = a(2,2); eig(3) = a(3,3)
 else
    p2 = (a(1,1)-q)**2 + (a(2,2)-q)**2 + (a(3,3)-q)**2 + 2.*p1
    p  = sqrt(p2/6.)
    b  = a
    do i=1,3
       b(i,i) = b(i,i) - q
    enddo
    b = b/p
    r = ( b(1,1)*(b(2,2)*b(3,3) - b(2,3)*b(3,2)) &
        - b(1,2)*(b(2,1)*b(3,3) - b(2,3)*b(3,1)) &
        + b(1,3)*(b(2,1)*b(3,2) - b(2,2)*b(3,1)) )/2.
    !
    ! r is cos(3*phi) up to rounding; clamp so acos stays in range
    !
    r = max(-1.0,min(1.0,r))
    phi = acos(r)/3.
    eig(1) = q + 2.*p*cos(phi)
    eig(3) = q + 2.*p*cos(phi + 2.*pi/3.)
    eig(2) = 3.*q - eig(1) - eig(3)   ! trace is invariant
 endif
 !
 ! sort descending
 !
 do i=1,2
    do j=i+1,3
       if (eig(j) > eig(i)) then
          tmp = eig(i); eig(i) = eig(j); eig(j) = tmp
       endif
    enddo
 enddo

end subroutine eigenvalues_sym3

end module analysis
