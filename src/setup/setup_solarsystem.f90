!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Setup asteroid orbits using data from the IAU Minor Planet Center
!
! :References: https://minorplanetcenter.net/data
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - asteroids  : *add distant minor bodies as km-sized dust particles*
!   - dtmax_in   : *time between dumps (e.g. 1 hr)*
!   - epoch      : *epoch to query ephemeris, YYYY-MMM-DD HH:MM:SS.fff, blank = today*
!   - np_apophis : *number of particles used to represent apophis (0=none; 1=sink; n=gas)*
!   - tmax_in    : *end time of simulation (e.g. 3 days)*
!   - scale_pos  : *scaling factor for apophis initial position (heliocentric)*
!   - scale_earth_sep : *scale geocentric Earth–Apophis separation (1=ephemeris; requires apophis_only=F)*
!   - apophis_spin_axis_x/y/z : *Apophis spin axis direction (normalized)*
!
! :Dependencies: centreofmass, eos_tillotson, infile_utils, io, kernel,
!   options, part, physcon, setbinary, setsolarsystem, setup_params,
!   spherical, timestep, units
!
 implicit none
 public :: setpart

 logical :: add_mars_moons ! new runtime parameter, default .true or .false

 integer :: np_apophis
 logical :: asteroids
 character(len=20) :: epoch,tmax_in,dtmax_in
 logical :: use_dem,apophis_only,use_dem_as_sinks
 character(len=256) :: apophis_shape_file

 real :: scale_vel
 real :: scale_pos
 real :: scale_earth_sep
 real :: scale_r_apophis
 real :: scale_rho
 real :: mass_apophis
 real :: apophis_spin_period
 real :: apophis_spin_axis(3)

 logical :: pack_settle
 real :: pack_expand
 real :: pack_phi

 private

contains
!----------------------------------------------------------------
!+
!  setup for solar system orbits
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use part,         only:nptmass,xyzmh_ptmass,vxyz_ptmass,idust,set_particle_type,&
                        grainsize,graindens,ndustlarge,ndusttypes,ndustsmall,ihacc,igas,idem
 use setbinary,     only:set_binary
 use units,         only:set_units,umass,udist,unit_density,unit_velocity,utime,in_code_units,in_units
 use physcon,       only:solarm,pi,au,km,solarr,ceresm,earthm,earthr,days,gg
 use io,            only:master,fatal,warning
 use timestep,      only:tmax,dtmax
 use damping,     only:idamp,tdyn_s
 use centreofmass,  only:reset_centreofmass
 use setsolarsystem,only:set_minor_planets,add_sun_and_planets,add_body
 use kernel,        only:hfact_default
 use eos_tillotson, only:rho_0,A
 use shape,         only:set_shape,get_mesh_geometry
 use options,       only:ieos
 use setup_params,  only:npart_total
 use orbits,        only:get_pericentre_distance,get_eccentricity
 use infile_utils,  only:get_options
 use ptmass,        only:isink_potential
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(inout) :: massoftype(:)
 real,              intent(inout) :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 integer :: ierr,i,nerr,n_apophis_part,i_apophis_first,i_apophis_last
 integer, parameter :: iearth = 4  ! Earth sink index when apophis_only=F
 !integer :: values(8),year,month,day
 real    :: period,semia,mtot,dx
 real    :: r_apophis,m_apophis,vol_apophis,rtidal,spsoundmin,r_grain,r_circ,r_cloud
 integer :: ierr_mesh,n_settle
 real    :: dr(3),sep_km,sep_re,rperi,rperi_km,rperi_re,ecc,vrel_kms
 real    :: dv(3)
!
! default runtime parameters
!
 tmax_in = '1000 yr'
 dtmax_in = '1 yr'
 asteroids = .true.
 np_apophis = 0
 use_dem_as_sinks = .false.
 use_dem = .false.
 apophis_only = .false.
 add_mars_moons = .false.
 !call date_and_time(values=values)
 !year = values(1); month = values(2); day = values(3)
 !write(epoch,"(i4.4,'-',i2.2,'-',i2.2)") year,month,day
 epoch='2029-04-10'   ! encounter is on Friday 13th April
 scale_vel=1.
 scale_pos=1.
 scale_earth_sep=1.
 scale_r_apophis=1.
 scale_rho=1.
 mass_apophis=0.
 apophis_shape_file='apophis.shape'
 apophis_spin_period = 0.
 apophis_spin_axis   = (/ 0., 0., 1. /)
 pack_settle = .false.
 pack_expand = 1.8
 pack_phi    = 0.64
 r_grain     = 0.
 r_circ      = 0.
 r_cloud     = 0.
 n_settle    = 0
!
! read runtime parameters from setup file
!
 if (id==master) print "(/,65('-'),1(/,a),/,65('-'),/)",&
   ' Welcome to the Superb Solar System Setup'

 call get_options(trim(fileprefix)//'.setup',id==master,ierr,&
                  read_setupfile,write_setupfile)
 if (ierr /= 0) stop 'rerun phantomsetup after editing .setup file'
!
! set units
!
 call set_units(mass=solarm,dist=km,G=1.d0)
!
! general parameters
!
 time  = 0.
 polyk = 0.
 gamma = 1.
 hfact = hfact_default
!
!--space available for injected gas particles
!
 npart = 0
 npart_total = 0
 npartoftype(:) = 0
 xyzh(:,:)  = 0.
 vxyzu(:,:) = 0.
 nptmass = 0

 semia  = 1.*au/udist  !  Earth
 mtot   = solarm/umass !  mass around which all bodies should orbit

 period = 2.*pi*sqrt(semia**3/mtot)
 tmax   = in_code_units(tmax_in,ierr,unit_type='time')
 if (ierr /= 0) call fatal('setup_solarsystem',' could not parse tmax')
 dtmax  = in_code_units(dtmax_in,ierr,unit_type='time')
 if (ierr /= 0) call fatal('setup_solarsystem',' could not parse dtmax')

 if (asteroids) then
    call set_minor_planets(npart,npartoftype,massoftype,xyzh,vxyzu,&
                           mtot,itype=idust,sample_orbits=.false.)
    print*,'npart = ',npart,' npartoftype = ',npartoftype(idust)
    !
    ! treat minor bodies as km-sized dust particles
    !
    ndustlarge = 1
    ndustsmall = 0
    ndusttypes = 1
    grainsize(ndustlarge) = km/udist         ! assume km-sized bodies
    graindens(ndustlarge) = 2./unit_density  ! 2 g/cm^3
 endif
 ! 
 ! add the planets
 !
 ierr = 0
 call add_sun_and_planets(nptmass,xyzmh_ptmass,vxyz_ptmass,mtot,nerr,epoch)
 if (nerr > 0) ierr = ierr + nerr

 if (apophis_only) nptmass = 0
 !
 ! add mars moons
 !
 if (add_mars_moons) then
    call add_body('phobos',nptmass,xyzmh_ptmass,vxyz_ptmass,mtot,nerr,epoch)
    xyzmh_ptmass(4,nptmass) = 1.08e16/umass ! mass in code units
    xyzmh_ptmass(5,nptmass) = 11.1 ! radius in code units (km)
    !override mass & radius
    call add_body('deimos',nptmass,xyzmh_ptmass,vxyz_ptmass,mtot,nerr,epoch)
    xyzmh_ptmass(4,nptmass) = 1.8e15/umass 
    xyzmh_ptmass(5,nptmass) = 6.0
    !override mass & radius
 end if
 !
 ! add the bringer of death
 !
 if (np_apophis > 0) then
    call add_body('apophis',nptmass,xyzmh_ptmass,vxyz_ptmass,mtot,nerr,epoch)
    if (nerr > 0) call warning('apophis','missing some information')

    r_apophis = xyzmh_ptmass(5,nptmass) * scale_r_apophis
    xyzmh_ptmass(5,nptmass) = r_apophis
    print "(a,1pg10.3)",' apophis radius scaled by ',scale_r_apophis

    vxyz_ptmass(1:3,nptmass) = vxyz_ptmass(1:3,nptmass)*scale_vel
    print "(a,1pg10.3)",' velocity of apophis scaled by ',scale_vel

    xyzmh_ptmass(1:3,nptmass) = xyzmh_ptmass(1:3,nptmass)*scale_pos
    print "(a,1pg10.3)",' initial position of apophis scaled by ',scale_pos

    if (.not.apophis_only .and. nptmass >= iearth) then
       if (abs(scale_earth_sep - 1.0) > 1.0e-6) then
          dr = xyzmh_ptmass(1:3,nptmass) - xyzmh_ptmass(1:3,iearth)
          xyzmh_ptmass(1:3,nptmass) = xyzmh_ptmass(1:3,iearth) + scale_earth_sep*dr
          print "(a,1pg10.3)",' Earth-Apophis separation scaled by ',scale_earth_sep
       endif
       dr = xyzmh_ptmass(1:3,nptmass) - xyzmh_ptmass(1:3,iearth)
       dv = vxyz_ptmass(1:3,nptmass) - vxyz_ptmass(1:3,iearth)
       rperi = get_pericentre_distance(xyzmh_ptmass(4,iearth),dr,dv)
       ecc = get_eccentricity(xyzmh_ptmass(4,iearth),dr,dv)
       rperi_km = rperi*udist/km
       rperi_re = rperi_km/(earthr/km)
       sep_km = sqrt(dot_product(dr,dr))*udist/km
       sep_re = sep_km/(earthr/km)
       vrel_kms = in_units(sqrt(dot_product(dv,dv)),'km/s')
       print "(a,1pg10.3,a,1pg10.3,a)",' pericentre distance = ',rperi_km,' km (',rperi_re,' R_Earth)'
       print "(a,1pg10.3)",' geocentric eccentricity = ',ecc
       print "(a,1pg10.3,a,1pg10.3,a)",' geocentric separation (initial) = ',sep_km,' km (',sep_re,' R_Earth)'
       print "(a,1pg10.3,a)",' Apophis-Earth relative velocity = ',vrel_kms,' km/s'
    elseif (abs(scale_earth_sep - 1.0) > 1.0e-6) then
       call warning('setup_solarsystem','scale_earth_sep ignored when apophis_only=T (Earth absent)')
    endif

    !
    ! volume of apophis: a sphere of r_apophis for a sink, or the body actually built by set_shape
    !
    vol_apophis = 4./3.*pi*r_apophis**3
    if (np_apophis > 1) then
       if (pack_settle) then
          !
          ! The body's volume is the SHAPE's volume, not a sphere's and not the
          ! cloud's: it fixes the mass and hence the bulk density. Taken exactly
          ! from the mesh here (divergence theorem), with no lattice probe.
          !
          call get_mesh_geometry(apophis_shape_file,r_apophis,vol_apophis,r_circ,&
                                 ierr_mesh,id,master)
          if (ierr_mesh /= 0) call warning('apophis',&
             'could not read shape for settling: falling back to a sphere of r_apophis')
          !
          ! The settled body must ENCLOSE the shape or the later crop truncates
          ! it, so the settle target is the circumradius, not r_apophis. Holding
          ! the grain radius fixed, that needs more grains than the shape alone
          ! will keep: n_settle/n_kept = V_sphere(r_circ)/V_shape.
          !
          r_grain = (3.*pack_phi*vol_apophis/(4.*pi*real(np_apophis)))**(1./3.)
          !
          ! 10% more grains than the circumscribing sphere needs exactly. Two
          ! things make the settled body come out smaller than the arithmetic
          ! says, and both shrink the radius: set_shape only matches the
          ! requested count to within a few per cent, and a settled random
          ! packing does not land exactly on pack_phi (0.66 measured against
          ! 0.64 assumed, which is a denser and so smaller body). 10% on the
          ! count is ~3% on the radius, which covers both. Costs a few per
          ! cent of settling time; the alternative is clipping the tips off
          ! the shape, which moddump_cropshape can only warn about.
          !
          n_settle = nint(1.10*np_apophis*(4./3.*pi*r_circ**3)/vol_apophis)
          r_cloud = pack_expand*r_circ
          call set_shape('random',id,master,n_settle,xyzmh_ptmass(1:3,nptmass),&
                         r_cloud,hfact,npart,xyzh,npart_total,&
                         objfile=apophis_shape_file,sphere_radius=r_cloud)
       else
          call set_shape('closepacked',id,master,np_apophis,xyzmh_ptmass(1:3,nptmass),r_apophis,&
                         hfact,npart,xyzh,npart_total,objfile=apophis_shape_file,vol=vol_apophis)
       endif
    endif
    !
    ! fix either the mass (mass_apophis > 0) or the density (scale_rho); the other follows from the volume
    !
    if (mass_apophis > 0.) then
       m_apophis = mass_apophis/umass
    else
       m_apophis = (rho_0*scale_rho/unit_density)*vol_apophis
    endif
    xyzmh_ptmass(4,nptmass) = m_apophis
    print "(a,2(es10.3,a))",' mass of apophis is ',m_apophis*umass,&
                            ' g or ',m_apophis*umass/ceresm,' ceres masses'
    print "(a,1pg10.3,a)",' density is ',m_apophis/vol_apophis*unit_density,' g/cm^3'
    if (np_apophis > 1 .and. .not.(use_dem .or. use_dem_as_sinks) .and. &
        abs(m_apophis/vol_apophis*unit_density - rho_0*scale_rho) > 0.02*rho_0*scale_rho) &
       call warning('apophis','SPH body density differs from rho_0*scale_rho used by the EOS')
    rtidal = (2.)**(1./3.)*(3.*vol_apophis/(4.*pi))**(1./3.)*(earthm/umass/m_apophis)**(1./3.)
    print "(3(a,1pg10.3),a)",' fluid Roche limit r_tidal is ',rtidal,' au, ',&
         rtidal*udist/km,' km, or ',rtidal*udist/earthr,' earth radii'

    if (np_apophis > 1) then
       !
       ! replace the sink particle with a ball of stuff
       !
       !call set_sphere('closepacked',id,master,0.,r_apophis,dx,hfact,npart,xyzh,npart_total,&
       !                xyz_origin=xyzmh_ptmass(1:3,nptmass),exactN=.true.,np_requested=np_apophis)

       n_apophis_part = npart !we are saving count here because will rewrite other variable later
       do i=1,npart
          vxyzu(1:3,i) = vxyz_ptmass(1:3,nptmass)
       enddo
       if (pack_settle) then
          !
          ! In settling mode the grain, not the body, is the physical object:
          ! its mass follows from its size and the grain material density, and
          ! must not change when the body is later cropped. So divide the body
          ! mass by the count the SHAPE will keep, not by the larger settling
          ! count. Dividing by npart instead would give every grain 1/2.4 of
          ! its proper density during the settle and then jump it back at the
          ! crop. A side effect worth having: the settling sphere then carries
          ! more than the body mass, but at the correct bulk density, so its
          ! dynamical time matches the final body's and tdyn_s is the same
          ! number for both settles.
          !
          massoftype(igas) = m_apophis / np_apophis
       else
          massoftype(igas) = m_apophis / npart
       endif
       npartoftype(igas) = npart
       nptmass = nptmass - 1
       i_apophis_first = 0 !initialize apophis sink index range to 0 before know if DEM used
       i_apophis_last  = 0 !initialize apophis sink index range to 0 before know if DEM used

       if (use_dem_as_sinks) then
          call replace_gas_with_dem(id,npart,npartoftype(igas),massoftype(igas),&
                                    xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass,hfact)
          isink_potential = 2
          i_apophis_first = nptmass - n_apophis_part + 1
          i_apophis_last  = nptmass
       elseif (use_dem) then
          massoftype(idem) = massoftype(igas)
          npartoftype(idem) = npartoftype(igas)
          massoftype(igas) = 0.
          npartoftype(igas) = 0
          if (pack_settle) then
             !
             ! r_grain was fixed above from the SHAPE volume and the requested
             ! kept-particle count, so that np_apophis grains fill the shape at
             ! pack_phi once the body has been cropped. It must NOT be re-derived
             ! from npart here: npart is the larger settling count, and it must
             ! not come from the cloud spacing, which is pack_expand times too
             ! large. h is both the grain radius and the neighbour search radius
             ! and is held fixed through settling (step_leapfrog skips h
             ! evolution when DEM particles are present).
             !
             do i=1,npart
                call set_particle_type(i,idem)
                xyzh(4,i) = r_grain
             enddo
             !
             ! Settling needs damping, so set it here rather than making the
             ! user know to add it by hand: write_options_damping returns
             ! early when idamp==0, so phantomsetup would otherwise never
             ! write the block and a fresh .in would silently run undamped.
             ! idamp=2 takes a dynamical time in seconds and scales itself,
             ! which is what we want: the constant-damping route (idamp=1)
             ! expects a rate in inverse code time, and the value carried
             ! over from relax_star is ~8 orders of magnitude too strong here.
             ! t_dyn = 1/sqrt(G*rho) from the body's own bulk density.
             !
             idamp  = 2
             tdyn_s = 1./sqrt(gg*rho_0*scale_rho)
             if (id==master) then
                print "(a)",' --- settling mode: loose cloud, shape cut after packing ---'
                print "(a,1pg10.3,a)",' damping t_dyn        = ',tdyn_s,' s (idamp=2)'
                print "(a,1pg10.3,a)",' shape volume         = ',vol_apophis*(udist/km)**3,' km^3'
                print "(a,1pg10.3,a)",' shape circumradius   = ',r_circ*udist/km,' km'
                print "(a,1pg10.3,a)",' DEM grain radius     = ',r_grain*udist/km,' km'
                print "(a,1pg10.3,a)",' initial cloud radius = ',r_cloud*udist/km,' km'
                print "(a,i9,a,i9)",   ' grains settling      = ',npart,'  -> after crop ~',np_apophis
                print "(a,1pg10.3)",   ' target packing frac  = ',pack_phi
             endif
          else
             do i=1,npart
                call set_particle_type(i,idem)
                xyzh(4,i) = xyzh(4,i) /hfact * 0.5 !set radius = to original particle spacing
             enddo
          endif
       endif

       if (apophis_spin_period > 0.) then
          call set_apophis_spin(id,apophis_spin_period,apophis_spin_axis,m_apophis,r_apophis,&
                                use_dem_as_sinks,n_apophis_part,i_apophis_first,i_apophis_last,&
                                xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass)
       endif
       !
       ! print quantities from the equation of state to give an idea of the timestep
       !
       if (ieos==23) then
          spsoundmin = sqrt(A/rho_0)/unit_velocity
          print "(a,1pg11.4,a)",'     sound speed min = ',spsoundmin*unit_velocity/km,' km/s'
          print "(a,1pg10.3,a)",' sound crossing time = ',(r_apophis/spsoundmin)*utime,' seconds'
       endif
    elseif (apophis_spin_period > 0.) then
       call set_apophis_spin(id,apophis_spin_period,apophis_spin_axis,m_apophis,r_apophis,&
                             .false.,0,nptmass,nptmass,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass)
    endif
 endif
 !
 ! set centre of mass as the origin
 !
 call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

 if (ierr /= 0) call fatal('setup','ERRORS during setup')

end subroutine setpart

!----------------------------------------------------------------
!+
!  impose solid-body spin on Apophis (SPH or DEM particles)
!+
!----------------------------------------------------------------
subroutine set_apophis_spin(id,period_s,spin_axis_in,m_body,r_body,dem_as_sinks,n_gas,&
                            i_sink_first,i_sink_last,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass)
 use part,        only:ispinx,ispiny,ispinz,iReff
 use units,       only:utime
 use physcon,     only:twopi,hours
 use io,          only:master,fatal
 use vectorutils, only:cross_product3D,unitvec,mag
 integer, intent(in) :: id,i_sink_first,i_sink_last,n_gas
 real,    intent(in) :: period_s,spin_axis_in(3),m_body,r_body
 logical, intent(in) :: dem_as_sinks
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 real :: omega,omega_crit,period_crit_s,spin_axis(3),spin_vec(3),r_cm(3),r_rel(3),dv(3)
 real :: pmass,reff
 integer :: i,n

 if (period_s <= 0.) return

 if (mag(spin_axis_in) <= 0.) call fatal('set_apophis_spin','zero spin axis vector')
 spin_axis = unitvec(spin_axis_in)

 omega = twopi/(period_s/utime)
 spin_vec = omega*spin_axis

 r_cm = 0.
 pmass = 0.
 if (dem_as_sinks) then
    n = i_sink_last - i_sink_first + 1
    do i=i_sink_first,i_sink_last
       pmass = pmass + xyzmh_ptmass(4,i)
       r_cm  = r_cm + xyzmh_ptmass(4,i)*xyzmh_ptmass(1:3,i)
    enddo
 else
    n = n_gas
    do i=1,n_gas
       r_cm = r_cm + xyzh(1:3,i)
    enddo
 endif
 if (n > 0) then
    if (dem_as_sinks) then
       r_cm = r_cm/pmass
    else
       r_cm = r_cm/real(n)
    endif
 endif

 if (dem_as_sinks) then
    do i=i_sink_first,i_sink_last
       r_rel = xyzmh_ptmass(1:3,i) - r_cm
       call cross_product3D(spin_vec,r_rel,dv)
       vxyz_ptmass(1:3,i) = vxyz_ptmass(1:3,i) + dv
       reff = xyzmh_ptmass(iReff,i)
       xyzmh_ptmass(ispinx,i) = omega*reff**2*spin_axis(1)
       xyzmh_ptmass(ispiny,i) = omega*reff**2*spin_axis(2)
       xyzmh_ptmass(ispinz,i) = omega*reff**2*spin_axis(3)
    enddo
 else
    do i=1,n_gas
       r_rel = xyzh(1:3,i) - r_cm
       call cross_product3D(spin_vec,r_rel,dv)
       vxyzu(1:3,i) = vxyzu(1:3,i) + dv
    enddo
    if (i_sink_first > 0) then
       reff = r_body
       xyzmh_ptmass(ispinx,i_sink_first) = omega*reff**2*spin_axis(1)
       xyzmh_ptmass(ispiny,i_sink_first) = omega*reff**2*spin_axis(2)
       xyzmh_ptmass(ispinz,i_sink_first) = omega*reff**2*spin_axis(3)
    endif
 endif

 if (r_body > 0.) then
    omega_crit = sqrt(m_body/r_body**3)
    period_crit_s = twopi/omega_crit*utime
    if (id == master) then
       print "(a,1pg10.3,a)",' Apophis spin period = ',period_s,' s (',period_s/hours,' hrs)'
       print "(a,3(1x,1pg10.3))",' Apophis spin axis  = ',spin_axis
       print "(a,1pg10.3,a)",' Apophis spin rate   = ',omega*utime,' rad/s'
       print "(a,1pg10.3)",' omega/omega_crit (fluid spin limit) = ',omega/omega_crit
       print "(a,1pg10.3,a)",' fluid spin-limit period = ',period_crit_s,' s (',period_crit_s/hours,' hrs)'
    endif
 endif

end subroutine set_apophis_spin

!----------------------------------------------------------------
!+
!  replace gas with discrete element method particles
!+
!----------------------------------------------------------------
subroutine replace_gas_with_dem(id,npart,ngas,pmass,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass,hfact)
 use part, only:iReff,ihacc
 use units, only:udist
 use physcon, only:km
 use io, only:master
 integer, intent(in)    :: id
 integer, intent(inout) :: npart,ngas,nptmass
 real, intent(inout) :: pmass,xyzh(:,:),vxyzu(:,:)
 real, intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 real, intent(in)    :: hfact
 integer :: i,j
 real :: dxij,dyij,dzij,dmin,reff,sep(npart)
 real :: sep_med,sep_min,sep_max

 if (npart < 1) return

 !
 ! DEM sphere radius from cropped lattice geometry (not r_apophis/40).
 ! Use per-particle Reff = 0.5*nearest-neighbour distance so ellipsoid/mesh
 ! surfaces do not share one radius (median) that overlaps close pairs on step 1.
 !
 sep_med = 0.
 sep_min = huge(sep_min)
 sep_max = 0.
 if (npart == 1) then
    sep(1) = xyzh(4,1)/hfact
 else
    do i=1,npart
       dmin = huge(dmin)
       do j=1,npart
          if (j == i) cycle
          dxij = xyzh(1,i) - xyzh(1,j)
          dyij = xyzh(2,i) - xyzh(2,j)
          dzij = xyzh(3,i) - xyzh(3,j)
          dmin = min(dmin,sqrt(dxij*dxij + dyij*dyij + dzij*dzij))
       enddo
       sep(i) = dmin
       sep_min = min(sep_min,dmin)
       sep_max = max(sep_max,dmin)
    enddo
    do i=2,npart
       dmin = sep(i)
       j = i - 1
       do while (j >= 1 .and. sep(j) > dmin)
          sep(j+1) = sep(j)
          j = j - 1
       enddo
       sep(j+1) = dmin
    enddo
    if (mod(npart,2) == 1) then
       sep_med = sep((npart+1)/2)
    else
       sep_med = 0.5*(sep(npart/2) + sep(npart/2+1))
    endif
 endif

 if (id == master) then
    print "(a,1pg12.4,a)",' DEM NN spacing (min/median/max) = ',sep_min*udist/km,&
          ' / ',sep_med*udist/km,' / ',sep_max*udist/km,' km'
 endif

 xyzmh_ptmass(:,nptmass+1:) = 0.
 do i=1,npart
    nptmass = nptmass + 1
    vxyz_ptmass(1:3,nptmass) = vxyzu(1:3,i)
    xyzmh_ptmass(1:3,nptmass) = xyzh(1:3,i)
    xyzmh_ptmass(4,nptmass) = pmass
    reff = 0.5*sep(i)
    xyzmh_ptmass(iReff,nptmass) = reff
    xyzmh_ptmass(ihacc,nptmass) = reff
 enddo
 npart = 0
 ngas = 0
 pmass = 0.

end subroutine replace_gas_with_dem

!----------------------------------------------------------------
!+
!  write setup parameters to file
!+
!----------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')

 write(iunit,"(a)") '# input file for solar system setup routines'
 call write_inopt(tmax_in,'tmax_in','end time of simulation (e.g. 3 days)',iunit)
 call write_inopt(dtmax_in,'dtmax_in','time between dumps (e.g. 1 hr)',iunit)
 call write_inopt(asteroids,'asteroids','add distant minor bodies as km-sized dust particles',iunit)
 call write_inopt(np_apophis,'np_apophis','number of particles used to represent apophis (0=none; 1=sink; n=gas)',iunit)
 call write_inopt(epoch,'epoch','epoch to query ephemeris, YYYY-MMM-DD HH:MM:SS.fff, blank = today',iunit)

 call write_inopt(use_dem,'use_dem','represent apophis with DEM particles (contact forces via the neighbour tree)',iunit)
 call write_inopt(use_dem_as_sinks,'use_dem_as_sinks','legacy: represent apophis with DEM sink particles (all-pairs, slow)',iunit)
 call write_inopt(apophis_only,'apophis_only','only add apophis',iunit)
 call write_inopt(add_mars_moons,'add_mars_moons','add Phobos and Deimos as point masses',iunit)

 call write_inopt(scale_vel,'scale_vel','scaling factor for apophis velocity',iunit)
 call write_inopt(scale_pos,'scale_pos','scaling factor for apophis initial position (heliocentric)',iunit)
 call write_inopt(scale_earth_sep,'scale_earth_sep',&
   'scale geocentric Earth-Apophis distance (1=ephemeris; apophis_only=F)',iunit)
 call write_inopt(scale_r_apophis,'scale_r_apophis','scaling factor for apophis radius',iunit)
 call write_inopt(scale_rho,'scale_rho','scaling factor for apophis bulk density',iunit)
 call write_inopt(mass_apophis,'mass_apophis','apophis mass in g (0 = from scale_rho and actual body volume)',iunit)
 call write_inopt(apophis_shape_file,'apophis_shape_file','shape config file for lattice cropping',iunit)
 call write_inopt(pack_settle,'pack_settle','start from a loose cloud and settle under gravity (shape cut afterwards)',iunit)
 call write_inopt(pack_expand,'pack_expand','initial cloud radius / target radius (settling mode)',iunit)
 call write_inopt(pack_phi,'pack_phi','packing fraction used to size DEM grains (0.64 = random close packing)',iunit)
 call write_inopt(apophis_spin_period,'apophis_spin_period','Apophis spin period in seconds (0=no spin)',iunit)
 call write_inopt(apophis_spin_axis(1),'apophis_spin_axis_x','Apophis spin axis, x component',iunit)
 call write_inopt(apophis_spin_axis(2),'apophis_spin_axis_y','Apophis spin axis, y component',iunit)
 call write_inopt(apophis_spin_axis(3),'apophis_spin_axis_z','Apophis spin axis, z component',iunit)

 close(iunit)

end subroutine write_setupfile

!----------------------------------------------------------------
!+
!  read setup parameters from file
!+
!----------------------------------------------------------------
subroutine read_setupfile(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 use io,           only:error
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 21
 integer :: nerr
 type(inopts), allocatable :: db(:)

 nerr = 0
 ierr = 0
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(tmax_in, 'tmax_in',db,errcount=nerr)
 call read_inopt(dtmax_in,'dtmax_in',db,errcount=nerr)
 call read_inopt(asteroids,'asteroids',db,errcount=nerr)
 call read_inopt(np_apophis,'np_apophis',db,min=0,errcount=nerr)
 call read_inopt(epoch,'epoch',db,errcount=nerr)

 call read_inopt(use_dem,'use_dem',db,errcount=nerr)
 call read_inopt(use_dem_as_sinks,'use_dem_as_sinks',db,default=.false.,errcount=nerr)
 call read_inopt(apophis_only,'apophis_only',db,errcount=nerr)
 call read_inopt(add_mars_moons,'add_mars_moons',db,errcount=nerr)

 call read_inopt(scale_vel,'scale_vel',db,default=1.0,errcount=nerr)
 call read_inopt(scale_pos,'scale_pos',db,default=1.0,errcount=nerr)
 call read_inopt(scale_earth_sep,'scale_earth_sep',db,default=1.0,errcount=nerr)
 call read_inopt(scale_r_apophis,'scale_r_apophis',db,default=1.0,errcount=nerr)
 call read_inopt(scale_rho,'scale_rho',db,default=1.0,errcount=nerr)
 call read_inopt(mass_apophis,'mass_apophis',db,min=0.,default=0.,errcount=nerr)
 call read_inopt(apophis_shape_file,'apophis_shape_file',db,default='apophis.shape',errcount=nerr)
 call read_inopt(pack_settle,'pack_settle',db,default=.false.,errcount=nerr)
 call read_inopt(pack_expand,'pack_expand',db,min=1.0,default=1.8,errcount=nerr)
 call read_inopt(pack_phi,'pack_phi',db,min=0.01,max=0.74,default=0.64,errcount=nerr)
 call read_inopt(apophis_spin_period,'apophis_spin_period',db,min=0.,errcount=nerr)
 call read_inopt(apophis_spin_axis(1),'apophis_spin_axis_x',db,errcount=nerr)
 call read_inopt(apophis_spin_axis(2),'apophis_spin_axis_y',db,errcount=nerr)
 call read_inopt(apophis_spin_axis(3),'apophis_spin_axis_z',db,errcount=nerr)

 call close_db(db)

 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

end subroutine read_setupfile

end module setup
