!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! Put a sinkless particle body onto a hyperbolic flyby past a planet.
!
! Written for the packing comparison: the settled and cropped bodies exist
! as dumps of DEM particles with no sinks at all, so moddump_addflyby is
! not usable (it requires an existing sink to act as the primary and calls
! fatal otherwise). Here the particle body itself is the primary.
!
! The encounter is specified by pericentre distance and velocity at
! infinity rather than by semi-major axis, because those are the two
! numbers a tidal-disruption study actually cares about: rp sets the peak
! tidal stress and v_inf sets how long it is applied for. Eccentricity
! follows as e = 1 + rp*v_inf^2/mu, so holding v_inf fixed while varying
! rp gives a sweep at constant encounter speed.
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: centreofmass, io, orbits, part, physcon, prompting,
!   setbinary, units
!
 use part,         only:idem,igas,isdead_or_accreted,nptmass,xyzmh_ptmass,vxyz_ptmass,ihacc,ihsoft
 use prompting,    only:prompt
 use centreofmass, only:reset_centreofmass
 use units,        only:udist,umass,utime
 use physcon,      only:km,earthm,earthr
 use io,           only:id,master,fatal

 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

 real :: rp_km      = 9600.
 real :: vinf_kms   = 5.9
 real :: d_km       = 4.0e5
 real :: m2_earth   = 1.0

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use setbinary, only:set_binary
 use orbits,    only:get_true_anomaly_from_separation
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real    :: m1,m2,rp,vinf,d,ecc,semia,f,mu
 real    :: xyzmh_tmp(5,2),vxyz_tmp(3,2),x1(3),v1(3),x2(3),v2(3)
 integer :: i,itype,nlive,ntmp,ierr

 if (nptmass > 0) call fatal('moddump_earthflyby','dump already has sinks; this expects a bare body')

 itype = idem
 if (npartoftype(idem) < 1) itype = igas
 nlive = 0
 do i=1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) nlive = nlive + 1
 enddo
 if (nlive < 2) call fatal('moddump_earthflyby','no live particles')
 m1 = massoftype(itype)*real(nlive)

 call prompt('Pericentre distance in km',rp_km,0.)
 call prompt('Velocity at infinity in km/s',vinf_kms,0.)
 call prompt('Initial separation in km',d_km,0.)
 call prompt('Perturber mass in Earth masses',m2_earth,0.)
 !
 ! to code units
 !
 rp   = rp_km*km/udist
 d    = d_km*km/udist
 vinf = vinf_kms*km/udist*utime
 m2   = m2_earth*earthm/umass
 mu   = m1 + m2                      ! G = 1 in code units
 !
 ! hyperbolic orbit from the two numbers that matter. e > 1 strictly, and
 ! a is negative, which is what set_binary expects for an unbound orbit.
 !
 ecc   = 1. + rp*vinf**2/mu
 semia = rp/(1.-ecc)
 !
 ! start incoming: negative true anomaly at the requested separation
 !
 f = -abs(get_true_anomaly_from_separation(semia,ecc,d))
 !
 ! set_binary solves the Kepler problem for us. Give it a scratch array:
 ! we want the relative orbit, not its two sinks.
 !
 ntmp = 0
 xyzmh_tmp = 0.
 vxyz_tmp  = 0.
 call set_binary(m1,m2,semia,ecc,0.,earthr/udist, &
                 xyzmh_tmp,vxyz_tmp,ntmp,ierr,f=f,verbose=(id==master))
 if (ierr /= 0) call fatal('moddump_earthflyby','set_binary failed',var='ierr',ival=ierr)

 x1 = xyzmh_tmp(1:3,1); v1 = vxyz_tmp(1:3,1)
 x2 = xyzmh_tmp(1:3,2); v2 = vxyz_tmp(1:3,2)
 !
 ! the body is assumed centred on its own centre of mass (the crop leaves
 ! it that way); move it bodily onto the primary's slot in the orbit
 !
 call reset_centreofmass(npart,xyzh,vxyzu)
 do i=1,npart
    if (isdead_or_accreted(xyzh(4,i))) cycle
    xyzh(1:3,i)  = xyzh(1:3,i)  + x1
    vxyzu(1:3,i) = vxyzu(1:3,i) + v1
 enddo
 !
 ! and add the perturber as the one and only sink
 !
 nptmass = 1
 xyzmh_ptmass(:,1)     = 0.
 xyzmh_ptmass(1:3,1)   = x2
 xyzmh_ptmass(4,1)     = m2
 xyzmh_ptmass(ihacc,1) = earthr/udist
 xyzmh_ptmass(ihsoft,1)= earthr/udist
 vxyz_ptmass(1:3,1)    = v2

 if (id==master) then
    print "(/,a)",       ' --- hyperbolic flyby set up ---'
    print "(a,1pg12.4,a)",' body mass        = ',m1*umass,' g'
    print "(a,i10)",      ' live particles   = ',nlive
    print "(a,1pg12.4,a)",' pericentre       = ',rp_km,' km'
    print "(a,1pg12.4,a)",' v at infinity    = ',vinf_kms,' km/s'
    print "(a,1pg12.4)",  ' eccentricity     = ',ecc
    print "(a,1pg12.4,a)",' initial sep      = ',d_km,' km'
    print "(a,1pg12.4,a)",' perturber mass   = ',m2*umass,' g'
 endif

end subroutine modify_dump

end module moddump
