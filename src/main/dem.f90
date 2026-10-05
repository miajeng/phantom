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
! :References: Schwartz+2012, Granular Matter 14, 363-380
!
! Optional tensile cohesion (linear spring when spheres are slightly separated)
! mimics van der Waals / regolith bond strength at sub-contact scale.
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

 public :: get_ssdem_force,dem_cohesion_summary,get_dem_dt

 real, public :: C_dem = 0.1          ! Safety factor on the contact timestep
 real, public :: ct_dem = 0.1         ! Tangential damping coefficient
 real, public :: epsilon_n_dem = 0.5  ! Normal coefficient of restitution (user-settable)
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
subroutine get_ssdem_force(Rsinki,Rsinkj,mi,mj,ddr,dx,dy,dz,fx,fy,fz,veli,velj,wi,wj,dtmin,bonded,dwi)
 use physcon,     only:pi
 use vectorutils, only:cross_product
 use units,       only:umass,utime,udist
 real, intent(in)    :: Rsinki,Rsinkj,mi,mj,ddr,dx,dy,dz,veli(3),velj(3),wi(3),wj(3)
 real, intent(inout) :: fx,fy,fz,dtmin
 logical, intent(in), optional :: bonded  ! true if i and j are in the same clump
 real, intent(inout), optional :: dwi(3)  ! spin rate of change of i (torque / moment of inertia)
 real :: r,overlap,gap,kn,kn_dem,kt_dem,kb_dem,coh_gap_max
 logical :: is_bond
 real :: cn,ct,reduced_mass,log_epsilon_n_dem,li,lj
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
    !print*,' fx = ',fx,' fy = ',fy,' fz = ',fz
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
 ct = 0.
 !print*,' cn = ',cn

 ! Normal damping force
 fx = fx - cn * u_n(1) / mi
 fy = fy - cn * u_n(2) / mi
 fz = fz - cn * u_n(3) / mi

 !----------------------------------------------------------------
 ! Tangential force and torque
 !----------------------------------------------------------------
 ! The tangential (friction) force ft acts at the contact point, so it
 ! also spins i up: torque = (-li*nvec) x ft, over I = 2/5 mi Ri^2.
 ! ct is still zero here, so ft is zero until friction is added.
 !
 ft = -ct * u_t
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
 !
 kn_dem = max(kn_cgs + kb_cgs,kt_cgs) / (umass/utime**2)
 if (kn_dem > 0. .and. mass_dem > 0.) then
    get_dem_dt = C_dem*sqrt(0.5*mass_dem/kn_dem)
 else
    get_dem_dt = huge(0.)
 endif

end function get_dem_dt

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