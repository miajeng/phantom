!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module dem
!
! This module implements the soft sphere discrete element method (DEM)
! for sink-sink interactions.
!
! epsilon is a parameter [0,1] that controls the restitution of the collision.
! epsilon = 1 is perfectly elastic, epsilon = 0 is perfectly inelastic.
! epsilon = 0.5 is a reasonable default.
!
! Rolling and twisting resistance (mu_r, mu_t > 0, with mu_s > 0) follow
! the pkdgrav spring-dashpot-plastic model (Zhang+2017; Hu+2021 Table 1),
! which stands in for the interlocking of irregular grains. The relative
! spin w_i - w_j of a contact splits into rolling (in the tangent plane)
! and twisting (about the normal); each winds up an angular displacement
! that drives a restoring couple, capped plastically:
!    M_R = -k_R d_R - C_R w_R,  |M_R| <= mu_r beta Rbar |F_n|
!    M_T = -k_T d_T - C_T w_T,  |M_T| <= mu_t beta Rbar mu_s |F_n|
! with k_R = k_n (beta Rbar)^2, C_R = C_n (beta Rbar)^2,
!      k_T = 2 k_s (beta Rbar)^2, C_T = 2 C_s (beta Rbar)^2.
! The couples act on spins only (no force), equal and opposite on i and j.
! Rigid co-rotation of the pair (w_i = w_j) feels nothing.
!
! :References: Schwartz+2012, Granular Matter 14, 363-380;
!              Zhang+2017, Icarus 294, 98-123 (rolling/twisting resistance);
!              Cundall & Strack 1979, Geotechnique 29, 47-65;
!              Luding 2008, Granular Matter 10, 235-246 (sliding projection);
!              Zhang+2018, ApJ 857, 15; Hu+2021, MNRAS 502, 5277 (cohesion)
!
! Two cohesion models (use one):
!
!  cohesion_pa > 0: the pkdgrav granular-bridge cohesion (Zhang+2018,
!  Hu+2021, after Sanchez & Scheeres 2014): touching grains attract with a
!  constant force F_C = c * A_eff, A_eff = (2*beta*Rbar)^2, Rbar =
!  RiRj/(Ri+Rj), where c (Pa) is the interparticle cohesion and beta the
!  contact-size (shape) parameter. It acts only while the grains overlap,
!  so a pair parts when pulled harder than F_C. Because F_C presses the
!  grains together, the contact normal force, and so the Coulomb friction
!  limit mu_s*|F_n|, grows with cohesion: a cohesive contact resists
!  sliding even with no external load.
!
!  kt_cgs > 0 (legacy, used in the thesis): a linear tensile spring that
!  pulls slightly separated spheres back together. It acts only across a
!  gap, so it never adds to the contact force or to friction.
!
! Tangential friction (mu_s > 0) is the linear spring-dashpot with a
! Coulomb limit of Cundall & Strack (1979), as used in pkdgrav (Schwartz+2012,
! Zhang+2017): each touching pair carries a tangential spring displacement
! xi, so a contact can stick (static friction) until |F_t| reaches
! mu_s*|F_n|, then slides. The history lives in part (icontact, xicontact).
!
! Optional clump bonds glue grains into boulders: two grains with the same
! nonzero clump ID (part%iclump) feel a two-sided spring toward touching,
! r = R_i + R_j, so a clump holds its shape in tension and compression.
! Grains are equal spheres, so the rest length is the same for every bond
! and no bond list is stored. The price is that there is no bond memory: a
! bond stretched past bond_reach stops pulling (the boulder breaks there),
! but re-forms if the same two grains come back within reach.
!
! :Owner: Daniel Price
!
 implicit none
 private

 public :: get_ssdem_force,dem_cohesion_summary,get_dem_dt,dem_friction_on,dem_friction_summary
 public :: get_fcoh,dem_contact_check

 real, public :: C_dem = 0.1          ! Safety factor on the contact timestep
 real, public :: epsilon_n_dem = 0.5  ! Normal coefficient of restitution (user-settable)
 real, public :: epsilon_t_dem = 0.5  ! Tangential coefficient of restitution, sets the tangential damping
 real, public :: mu_s = 0.5           ! Static (sliding) friction coefficient; 0 = frictionless
 real, public :: ks_cgs = 0.          ! Tangential spring constant (g/s^2); 0 = (2/7) kn_cgs
 real, public :: cohesion_pa = 0.     ! Interparticle cohesion c (Pa), pkdgrav granular bridge; 0 = off
 real, public :: beta_dem = 0.5       ! Shape parameter: contact radius = beta*Rbar (pkdgrav)
 real, public :: mu_r = 0.            ! Rolling friction coefficient (pkdgrav "gravel": 1.05); 0 = off
 real, public :: mu_t = 0.            ! Twisting friction coefficient (pkdgrav "gravel": 1.3); 0 = off
 real, public :: kn_cgs = 1e7         ! Spring constant (e.g. 10^4 kg/s^2 = 10^7 g/s^2)
 real, public :: kt_cgs = 0.          ! Tensile spring constant (g/s^2 per cm gap); 0 = no cohesion
 real, public :: coh_gap_max_cgs = 0. ! Max surface gap (cm) for cohesive bond; 0 = use 1% of mean radius
 real, public :: kb_cgs = 0.          ! Clump bond spring constant (g/s^2 per cm stretch); 0 = no clumps
 real, public :: bond_reach = 0.1     ! Max surface gap for a clump bond, as a fraction of R_i + R_j

contains

!----------------------------------------------------------------
!+
!  Soft-sphere DEM normal force (Hooke's law)
!  Implements Eq. (3) from Schwartz+2012 for overlapping spheres
!+
!----------------------------------------------------------------
subroutine get_ssdem_force(Rsinki,Rsinkj,mi,mj,ddr,dx,dy,dz,fx,fy,fz,veli,velj,wi,wj,dtmin,bonded,dwi,&
                           i,jorig,dt)
 use physcon,     only:pi
 use vectorutils, only:cross_product
 use units,       only:umass,utime,udist
 use part,        only:maxcontact,icontact,xicontact,icontact_new,xicontact_new,rotcontact,rotcontact_new
 use io,          only:fatal
 real, intent(in)    :: Rsinki,Rsinkj,mi,mj,ddr,dx,dy,dz,veli(3),velj(3),wi(3),wj(3)
 real, intent(inout) :: fx,fy,fz,dtmin
 logical, intent(in), optional :: bonded  ! true if i and j are in the same clump
 real, intent(inout), optional :: dwi(3)  ! spin rate of change of i (torque / moment of inertia)
 integer,         intent(in), optional :: i      ! local index of i, for its contact history
 integer(kind=8), intent(in), optional :: jorig  ! permanent ID (iorig) of j
 real,            intent(in), optional :: dt     ! step over which the contact history advances
 real :: ks,cs,log_epsilon_t_dem,fn,ftmax,ftmag,xi0mag,xi(3),fcoh
 real :: rbar,wrel(3),wrol(3),wtw,rol(3),rol0mag,tw,kr,cr,ktw,ctw,mr(3),mrmax,mrmag,mt,mtmax
 integer :: k,kfree,kold
 real :: r,overlap,gap,kn,kn_dem,kt_dem,kb_dem,coh_gap_max
 logical :: is_bond
 real :: cn,reduced_mass,log_epsilon_n_dem,li,lj
 real :: nvec(3),vrel(3),n_cross_wi(3),n_cross_wj(3),u_dot_n,u_n(3),u_t(3),ft(3)

 !----------------------------------------------------------------
 ! Normal force
 !----------------------------------------------------------------
 r = 1.0 / ddr
  ! Normal unit vector
 nvec(1) = dx * ddr
 nvec(2) = dy * ddr
 nvec(3) = dz * ddr

 overlap = Rsinki + Rsinkj - r
 gap = r - Rsinki - Rsinkj
 kn = 0.
 kn_dem = kn_cgs / (umass/utime**2)  ! convert to code units
 kt_dem = kt_cgs / (umass/utime**2)
 kb_dem = kb_cgs / (umass/utime**2)
 is_bond = .false.
 if (present(bonded)) is_bond = bonded .and. kb_dem > 0. .and. gap < bond_reach*(Rsinki + Rsinkj)
 if (coh_gap_max_cgs > 0.) then
    coh_gap_max = coh_gap_max_cgs / udist
 else
    coh_gap_max = 0.01 * (Rsinki + Rsinkj)
 endif
 if (is_bond) then
    ! clump bond: spring toward touching, pulling when stretched (overlap < 0)
    ! and pushing, on top of the contact spring, when squeezed
    kn = kb_dem
    if (overlap > 0.) kn = kn_dem + kb_dem
    fx = fx + kn * overlap * nvec(1) / mi
    fy = fy + kn * overlap * nvec(2) / mi
    fz = fz + kn * overlap * nvec(3) / mi
 elseif (overlap > 0.0) then
    kn = kn_dem
    fx = fx + kn * overlap * nvec(1) / mi
    fy = fy + kn * overlap * nvec(2) / mi
    fz = fz + kn * overlap * nvec(3) / mi
    if (cohesion_pa > 0.) then
       ! granular-bridge cohesion, pulling i towards j (nvec points j to i)
       fcoh = get_fcoh(Rsinki,Rsinkj)
       fx = fx - fcoh * nvec(1) / mi
       fy = fy - fcoh * nvec(2) / mi
       fz = fz - fcoh * nvec(3) / mi
    endif
 elseif (kt_dem > 0. .and. gap > 0. .and. gap < coh_gap_max) then
    ! tensile spring when spheres are bonded but not overlapping
    kn = kt_dem
    fx = fx - kn * gap * nvec(1) / mi
    fy = fy - kn * gap * nvec(2) / mi
    fz = fz - kn * gap * nvec(3) / mi
 endif
 !
 ! All of fx,fy,fz is the acceleration of i, so forces are divided by mi.
 ! (They were divided by mj, which only agreed for equal-mass grains.)
 !

 !----------------------------------------------------------------
 ! Damping force
 !----------------------------------------------------------------

 ! Cross products: n x omega, where omega is the spin vector of the sphere
 n_cross_wi = cross_product(nvec,wi)
 n_cross_wj = cross_product(nvec,wj)

 ! Eqs. (9) and (10)from Schwartz+2012
 li = (Rsinki**2 - Rsinkj**2 + r**2) / (2.0 * r)
 lj = (Rsinkj**2 - Rsinki**2 + r**2) / (2.0 * r)

 ! Relative velocity at contact point (Eq. 8 from Schwartz+2012).
 ! nvec points from j to i, so the contact point is at x_i - li*nvec on i
 ! and x_j + lj*nvec on j, moving at v_i + li*(n x w_i) and v_j - lj*(n x w_j).
 ! Both spin terms therefore enter with a plus sign: with a minus on the
 ! j term, two touching grains rotating rigidly together would show a
 ! spurious sliding velocity. (Harmless until now, as spins were zero.)
 vrel = veli - velj + li * n_cross_wi + lj * n_cross_wj

 ! Normal and tangential components
 u_dot_n = dot_product(vrel, nvec)
 u_n = u_dot_n * nvec
 u_t = vrel - u_n

 ! Eqn (15) from Schwartz+2012
 reduced_mass = mj *  mi / (mj + mi)
 log_epsilon_n_dem = log(epsilon_n_dem)
 cn = -2.0 * sqrt(reduced_mass * kn) * log_epsilon_n_dem / sqrt(pi**2 + log_epsilon_n_dem**2)

 ! Normal damping force
 fx = fx - cn * u_n(1) / mi
 fy = fy - cn * u_n(2) / mi
 fz = fz - cn * u_n(3) / mi

 !----------------------------------------------------------------
 ! Tangential (friction) force, with contact history
 !----------------------------------------------------------------
 ! Only for touching grains whose caller passes the history (the DEM
 ! particle path in force.F90); the sink-sink path stays frictionless.
 !
 ! Linear spring-dashpot with a Coulomb limit (Cundall & Strack 1979):
 !    ft = -ks*xi - cs*u_t,   |ft| <= mu_s*|F_n|
 ! xi is the tangential spring stretch, accumulated as u_t*dt over the
 ! life of the contact. Each step it is first rotated into the current
 ! tangent plane, keeping its length, as the pair rolls round each other.
 ! If |ft| would exceed the Coulomb limit the contact slides: ft is cut
 ! back to the limit and xi reset to match it (Luding 2008), so on
 ! re-sticking the spring starts from the sliding force, not beyond it.
 ! F_n is the contact normal force, spring plus dashpot, as in pkdgrav
 ! (Hu+2021, Table 1); a separating contact (F_n < 0) has no friction.
 ! Granular-bridge cohesion is not added here, but it squeezes the pair
 ! until the spring carries it, so the limit becomes mu_s*(load + F_C).
 !
 ! i keeps its own copy of xi for the pair and j keeps the opposite one;
 ! the two are built from equal and opposite u_t, so they stay equal and
 ! opposite to roundoff and the pair forces obey Newton's third law.
 !
 ft = 0.
 if (present(i) .and. present(jorig) .and. present(dt) .and. mu_s > 0. .and. overlap > 0.) then
    if (ks_cgs > 0.) then
       ks = ks_cgs / (umass/utime**2)
    else
       ks = (2./7.) * kn_dem   ! equal normal and tangential frequencies (Schwartz+2012)
    endif
    log_epsilon_t_dem = log(epsilon_t_dem)
    cs = -2.0 * sqrt(reduced_mass * ks) * log_epsilon_t_dem / sqrt(pi**2 + log_epsilon_t_dem**2)

    ! xi at the start of the step, from the committed list (0 for a new contact)
    xi  = 0.
    rol = 0.
    tw  = 0.
    kold = 0
    do k=1,maxcontact
       if (icontact(k,i) == jorig) then
          kold = k
          exit
       endif
    enddo
    if (kold > 0) then
       xi  = xicontact(:,kold,i)
       rol = rotcontact(1:3,kold,i)
       tw  = rotcontact(4,kold,i)
    endif
    xi0mag = sqrt(dot_product(xi,xi))
    if (xi0mag > 0.) then
       xi = xi - dot_product(xi,nvec) * nvec
       ftmag = sqrt(dot_product(xi,xi))
       if (ftmag > 0.) xi = xi * (xi0mag / ftmag)
    endif
    xi = xi + u_t * dt

    ! trial force, then the Coulomb limit
    ft    = -ks * xi - cs * u_t
    fn    = kn * overlap - cn * u_dot_n
    ftmax = mu_s * max(fn,0.)
    ftmag = sqrt(dot_product(ft,ft))
    if (ftmag > ftmax) then
       if (ftmag > 0.) then
          ft = ft * (ftmax / ftmag)
       else
          ft = 0.
       endif
       xi = -(ft + cs * u_t) / ks
    endif

    !
    ! rolling and twisting resistance: couples from the relative spin
    !
    mr = 0.
    mt = 0.
    if (mu_r > 0. .or. mu_t > 0.) then
       rbar = Rsinki*Rsinkj/(Rsinki + Rsinkj)
       wrel = wi - wj
       wtw  = dot_product(wrel,nvec)       ! twisting rate, about nvec
       wrol = wrel - wtw * nvec            ! rolling rate, in the tangent plane
       if (mu_r > 0.) then
          ! rotate the rolling displacement into the current tangent plane
          rol0mag = sqrt(dot_product(rol,rol))
          if (rol0mag > 0.) then
             rol   = rol - dot_product(rol,nvec) * nvec
             mrmag = sqrt(dot_product(rol,rol))
             if (mrmag > 0.) rol = rol * (rol0mag / mrmag)
          endif
          rol   = rol + wrol * dt
          kr    = kn * (beta_dem*rbar)**2
          cr    = cn * (beta_dem*rbar)**2
          mr    = -kr * rol - cr * wrol
          mrmax = mu_r * beta_dem * rbar * max(fn,0.)
          mrmag = sqrt(dot_product(mr,mr))
          if (mrmag > mrmax) then
             if (mrmag > 0.) then
                mr = mr * (mrmax / mrmag)
             else
                mr = 0.
             endif
             rol = -(mr + cr * wrol) / kr
          endif
       else
          rol = 0.
       endif
       if (mu_t > 0.) then
          ! the twisting angle is a scalar about nvec, the same for i and j
          tw    = tw + wtw * dt
          ktw   = 2. * ks * (beta_dem*rbar)**2
          ctw   = 2. * cs * (beta_dem*rbar)**2
          mt    = -ktw * tw - ctw * wtw
          mtmax = mu_t * beta_dem * rbar * mu_s * max(fn,0.)
          if (abs(mt) > mtmax) then
             mt = sign(mtmax,mt)
             tw = -(mt + ctw * wtw) / ktw
          endif
       else
          tw = 0.
       endif
       if (present(dwi)) dwi = dwi + (mr + mt * nvec) / (0.4 * mi * Rsinki**2)
    endif

    ! record the trial history of this contact for i
    kfree = 0
    do k=1,maxcontact
       if (icontact_new(k,i) == 0) then
          kfree = k
          exit
       endif
    enddo
    if (kfree == 0) call fatal('dem','more than maxcontact contacts on one grain',var='maxcontact',ival=maxcontact)
    icontact_new(kfree,i)       = jorig
    xicontact_new(:,kfree,i)    = xi
    rotcontact_new(1:3,kfree,i) = rol
    rotcontact_new(4,kfree,i)   = tw
 endif
 !
 ! ft acts at the contact point, so it also spins i up:
 ! torque = (-li*nvec) x ft, over I = 2/5 mi Ri^2
 !
 fx = fx + ft(1) / mi
 fy = fy + ft(2) / mi
 fz = fz + ft(3) / mi
 if (present(dwi)) dwi = dwi - li * cross_product(nvec,ft) / (0.4 * mi * Rsinki**2)

 if (kn > 0.) dtmin = min(dtmin,sqrt(reduced_mass/kn))

end subroutine get_ssdem_force

!----------------------------------------------------------------
!+
!  Global timestep constraint from the DEM contact springs.
!
!  The pairwise constraint applied in get_ssdem_force is
!  sqrt(reduced_mass/kn), with reduced_mass = mi*mj/(mi+mj). That is
!  smallest when both particles are the lightest in the run, giving
!  reduced_mass = m_min/2. The stiffest pair is therefore bounded by
!  the smallest particle mass alone, so this constraint is global: it
!  needs no per-pair reduction and is identical on every MPI rank.
!
!  Returns huge() when DEM is inactive so it never limits the step.
!+
!----------------------------------------------------------------
real function get_dem_dt(mass_dem)
 use units, only:umass,utime
 real, intent(in) :: mass_dem
 real :: kn_dem
 !
 ! The stiffest spring in the problem sets the step, and there are two:
 ! the normal contact spring kn and the tensile cohesion spring kt, both
 ! applied in get_ssdem_force. Using kn alone silently under-resolves any
 ! run with kt > kn, which is exactly what a cohesion sweep reaches for.
 ! A squeezed clump bond is kn + kb, so that is the stiffest when kb > 0.
 ! The tangential spring ks defaults to (2/7) kn, so it only matters if set.
 !
 kn_dem = max(kn_cgs + kb_cgs,kt_cgs,ks_cgs) / (umass/utime**2)
 if (kn_dem > 0. .and. mass_dem > 0.) then
    get_dem_dt = C_dem*sqrt(0.5*mass_dem/kn_dem)
 else
    get_dem_dt = huge(0.)
 endif

end function get_dem_dt

!----------------------------------------------------------------
!+
!  Granular-bridge cohesive force between touching grains of radii
!  Ri, Rj (code units): F_C = c * (2*beta*Rbar)^2, Rbar = RiRj/(Ri+Rj)
!+
!----------------------------------------------------------------
real function get_fcoh(Ri,Rj)
 use units, only:umass,udist,utime
 real, intent(in) :: Ri,Rj
 real :: rbar,c_code

 rbar   = Ri*Rj/(Ri + Rj)
 c_code = 10.*cohesion_pa / (umass/(udist*utime**2))  ! Pa -> dyn/cm^2 -> code units
 get_fcoh = c_code * (2.*beta_dem*rbar)**2

end function get_fcoh

!----------------------------------------------------------------
!+
!  True if DEM grains feel tangential friction (and so need the
!  contact history and spins)
!+
!----------------------------------------------------------------
logical function dem_friction_on()

 dem_friction_on = (mu_s > 0.)

end function dem_friction_on

!----------------------------------------------------------------
!+
!  Print the DEM friction settings (once, when friction is switched on)
!+
!----------------------------------------------------------------
subroutine dem_friction_summary
 use io, only:iprint

 write(iprint,"(/,a)") ' DEM tangential friction enabled'
 write(iprint,"(a,1pg12.4)") '   mu_s = ',mu_s
 if (ks_cgs > 0.) then
    write(iprint,"(a,1pg12.4,a)") '   ks_cgs = ',ks_cgs,' g/s^2'
 else
    write(iprint,"(a,1pg12.4,a)") '   ks_cgs = (2/7) kn_cgs = ',2./7.*kn_cgs,' g/s^2'
 endif
 write(iprint,"(a,1pg12.4)") '   epsilon_t = ',epsilon_t_dem
 if (mu_r > 0. .or. mu_t > 0.) then
    write(iprint,"(a,1pg12.4,a,1pg12.4,a,1pg12.4)") '   rolling mu_r = ',mu_r,', twisting mu_t = ',mu_t,&
                                                   ', beta = ',beta_dem
 endif

end subroutine dem_friction_summary

!----------------------------------------------------------------
!+
!  Once per run, for grains of radius R (code units): report the static
!  overlap the cohesive force alone causes, F_C/kn, and warn if it is
!  more than 1% of R. pkdgrav sets kn so overlaps stay below 1% of the
!  smallest radius (Hu+2021); a soft spring against strong cohesion
!  squashes grains into each other and makes the packing unphysical.
!+
!----------------------------------------------------------------
subroutine dem_contact_check(R)
 use io,    only:iprint,warning
 use units, only:umass,utime,udist
 real, intent(in) :: R
 real :: fcoh,kn_dem,ovfrac

 if (cohesion_pa <= 0.) return
 fcoh   = get_fcoh(R,R)
 kn_dem = kn_cgs / (umass/utime**2)
 ovfrac = fcoh/kn_dem/R
 write(iprint,"(/,a)") ' DEM granular-bridge cohesion enabled'
 write(iprint,"(a,1pg12.4,a,1pg12.4)") '   c = ',cohesion_pa,' Pa, beta = ',beta_dem
 write(iprint,"(a,1pg12.4,a)") '   F_C (equal grains) = ',fcoh*umass*udist/utime**2,' dyn'
 write(iprint,"(a,1pg12.4,a)") '   cohesive overlap F_C/kn = ',ovfrac,' grain radii'
 if (ovfrac > 0.01) call warning('dem','cohesive overlap F_C/kn exceeds 1% of the grain radius: raise kn_cgs',&
                                 var='F_C/(kn R)',val=ovfrac)

end subroutine dem_contact_check

!----------------------------------------------------------------
!+
!  Print DEM cohesion and clump bond settings after reading the .in file
!+
!----------------------------------------------------------------
subroutine dem_cohesion_summary
 use io,      only:iprint
 use units,   only:umass,utime
 real :: kt_dem

 if (kb_cgs > 0.) then
    write(iprint,"(/,a)") ' DEM clump bonds enabled'
    write(iprint,"(a,1pg12.4,a)") '   kb_cgs = ',kb_cgs,' g/s^2 per cm stretch'
    write(iprint,"(a,1pg12.4,a)") '   bond_reach = ',bond_reach,' x (R_i + R_j) surface gap'
 endif
 if (kt_cgs <= 0.) return
 kt_dem = kt_cgs / (umass/utime**2)
 write(iprint,"(/,a)") ' DEM tensile cohesion enabled'
 write(iprint,"(a,1pg12.4,a)") '   kt_cgs = ',kt_cgs,' g/s^2 per cm surface gap'
 if (coh_gap_max_cgs > 0.) then
    write(iprint,"(a,1pg12.4,a)") '   coh_gap_max = ',coh_gap_max_cgs,' cm'
 else
    write(iprint,"(a)") '   coh_gap_max = 1% of (R_i + R_j) per pair (default)'
 endif
 write(iprint,"(a,1pg12.4)") '   kt (code units) = ',kt_dem
end subroutine dem_cohesion_summary

end module dem