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
!              Zhang+2018, ApJ 857, 15; Hu+2021, MNRAS 502, 5277 (cohesion);
!              Potyondy & Cundall 2004, IJRMMS 41, 1329-1364 (parallel bonds)
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
! Optional clump bonds glue grains into boulders, as parallel bonds
! (Potyondy & Cundall 2004). Bonds are made once, by the first force
! evaluation, between grains with the same nonzero clump ID (part%iclump)
! whose surface gap is below bond_reach*(R_i+R_j); each grain keeps a list
! of its bonded partners (part%ibond), written to dumps. A bond is a disc
! of radius r_b = bond_lambda*min(R_i,R_j) (area A, second moments
! I = pi r_b^4/4, J = 2I) joining the two grains, and carries
!   tension/compression  F_n = kb*(R_i+R_j-r), rest length touching
!   shear                F_s = -k_bs*xi_b,    k_bs = (2/7) kb
!   bending              M_b = -kb*(I/A)*th_b
!   twisting             M_t = -k_bs*(J/A)*th_t
! each with a dashpot, from the shear displacement and the bending and
! twisting angles accumulated since the bond was made. It breaks, for good,
! when the tensile stress  sigma = -F_n/A + |M_b| r_b/I  exceeds
! bond_sigma_pa, or the shear stress  tau = |F_s|/A + |M_t| r_b/J  exceeds
! bond_tau_pa (Potyondy & Cundall 2004, Eq. 6), or when stretched past
! bond_reach. Strengths of 0 switch that criterion off, leaving the stretch
! limit alone. A broken bond leaves an ordinary contact: friction and
! cohesion then apply. A bonded pair feels no contact friction.
!
! :Owner: Daniel Price
!
 implicit none
 private

 public :: get_ssdem_force,dem_cohesion_summary,get_dem_dt,dem_friction_on,dem_friction_summary
 public :: get_fcoh,dem_contact_check,dem_commit_history,dem_bonds_on

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
 real, public :: bond_lambda = 1.     ! Bond radius as a fraction of the smaller grain radius
 real, public :: bond_sigma_pa = 0.   ! Bond tensile strength (Pa); 0 = no stress criterion
 real, public :: bond_tau_pa = 0.     ! Bond shear strength (Pa); 0 = no stress criterion
 logical, public :: dem_bonds_forming = .false.  ! true only for the force evaluation that makes the bonds

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
 use part,        only:maxcontact,icontact,xicontact,icontact_new,xicontact_new,rotcontact,rotcontact_new,&
                       maxbond,ibond,ibond_new,xibond,xibond_new,rotbond,rotbond_new
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
 integer :: k,kfree,kold,kbond
 real :: xib(3),rob(3),twb,rb,ab,ib,jb,kbs,cbn,cbs,mrb(3),mtb,sig,tau,ftb(3)
 logical :: bond_hist
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

 reduced_mass = mj *  mi / (mj + mi)

 !----------------------------------------------------------------
 ! Clump bond: is there one, and does it hold?
 !----------------------------------------------------------------
 is_bond   = .false.
 bond_hist = present(i) .and. present(jorig) .and. present(dt)
 if (present(bonded)) then
    if (bonded .and. kb_dem > 0.) then
       if (bond_hist .and. allocated(ibond)) then
          kbond = 0
          do k=1,maxbond
             if (ibond(k,i) == jorig) then
                kbond = k
                exit
             endif
          enddo
          if (kbond > 0 .or. (dem_bonds_forming .and. gap < bond_reach*(Rsinki + Rsinkj))) then
             ! bond springs at the start of the step (zero for a bond being made),
             ! turned with the pair into the current tangent plane, then advanced
             xib = 0.; rob = 0.; twb = 0.
             if (kbond > 0) then
                xib = xibond(:,kbond,i)
                rob = rotbond(1:3,kbond,i)
                twb = rotbond(4,kbond,i)
             endif
             call turn_into_plane(xib,nvec)
             call turn_into_plane(rob,nvec)
             wrel = wi - wj
             wtw  = dot_product(wrel,nvec)
             wrol = wrel - wtw * nvec
             xib  = xib + u_t * dt
             rob  = rob + wrol * dt
             twb  = twb + wtw * dt

             rb  = bond_lambda * min(Rsinki,Rsinkj)
             ab  = pi * rb**2
             ib  = 0.25 * pi * rb**4
             jb  = 2. * ib
             kbs = (2./7.) * kb_dem
             log_epsilon_n_dem = log(epsilon_n_dem)
             log_epsilon_t_dem = log(epsilon_t_dem)
             cbn = -2.0 * sqrt(reduced_mass * kb_dem) * log_epsilon_n_dem / sqrt(pi**2 + log_epsilon_n_dem**2)
             cbs = -2.0 * sqrt(reduced_mass * kbs) * log_epsilon_t_dem / sqrt(pi**2 + log_epsilon_t_dem**2)
             ftb = -kbs * xib - cbs * u_t
             mrb = -kb_dem * (ib/ab) * rob - cbn * (ib/ab) * wrol
             mtb = -kbs * (jb/ab) * twb - cbs * (jb/ab) * wtw

             ! peak stresses in the bond from its springs (tension positive)
             sig = -kb_dem * overlap / ab + kb_dem * sqrt(dot_product(rob,rob)) * rb / ab
             tau = kbs * (sqrt(dot_product(xib,xib)) + abs(twb) * rb) / ab
             is_bond = (gap < bond_reach*(Rsinki + Rsinkj))
             if (bond_sigma_pa > 0.) is_bond = is_bond .and. sig <= stress_code(bond_sigma_pa)
             if (bond_tau_pa > 0.)   is_bond = is_bond .and. tau <= stress_code(bond_tau_pa)
          endif
       else
          ! no bond memory available (e.g. no history passed): stretch limit only
          is_bond = gap < bond_reach*(Rsinki + Rsinkj)
       endif
    endif
 endif
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

 ! Eqn (15) from Schwartz+2012
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
 if (bond_hist .and. .not.is_bond .and. mu_s > 0. .and. overlap > 0.) then
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
 ! an intact bond: its shear force goes on as ft, its bending and twisting
 ! couples straight onto the spin, and it is recorded in the trial list
 !
 if (is_bond .and. bond_hist .and. allocated(ibond)) then
    ft = ftb
    if (present(dwi)) dwi = dwi + (mrb + mtb * nvec) / (0.4 * mi * Rsinki**2)
    kfree = 0
    do k=1,maxbond
       if (ibond_new(k,i) == 0) then
          kfree = k
          exit
       endif
    enddo
    if (kfree == 0) call fatal('dem','more than maxbond bonds on one grain',var='maxbond',ival=maxbond)
    ibond_new(kfree,i)      = jorig
    xibond_new(:,kfree,i)   = xib
    rotbond_new(1:3,kfree,i) = rob
    rotbond_new(4,kfree,i)   = twb
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
!  turn a vector into the plane normal to nvec, keeping its length
!  (how a spring displacement follows a pair that rolls round itself)
!+
!----------------------------------------------------------------
pure subroutine turn_into_plane(v,nvec)
 real, intent(inout) :: v(3)
 real, intent(in)    :: nvec(3)
 real :: v0,v1

 v0 = sqrt(dot_product(v,v))
 if (v0 > 0.) then
    v  = v - dot_product(v,nvec) * nvec
    v1 = sqrt(dot_product(v,v))
    if (v1 > 0.) v = v * (v0 / v1)
 endif

end subroutine turn_into_plane

!----------------------------------------------------------------
!+
!  a stress in Pa, in code units
!+
!----------------------------------------------------------------
real function stress_code(pa)
 use units, only:umass,udist,utime
 real, intent(in) :: pa

 stress_code = 10.*pa / (umass/(udist*utime**2))

end function stress_code

!----------------------------------------------------------------
!+
!  True if boulders are glued: clump bonds on and a clump in the body
!+
!----------------------------------------------------------------
logical function dem_bonds_on(npart)
 use part, only:iclump
 integer, intent(in) :: npart

 dem_bonds_on = (kb_cgs > 0.)
 if (dem_bonds_on) dem_bonds_on = any(iclump(1:npart) /= 0)

end function dem_bonds_on

!----------------------------------------------------------------
!+
!  End of a converged step (or of the force evaluation that made the
!  bonds): commit the trial contact and bond lists built by the force.
!  A bond is kept only if both grains still list each other, so the
!  rare bond judged broken by one grain and not the other (roundoff at
!  the threshold) breaks for both, as do bonds to accreted grains.
!+
!----------------------------------------------------------------
subroutine dem_commit_history(npart,xyzh)
 use part, only:icontact,icontact_new,xicontact,xicontact_new,rotcontact,rotcontact_new,&
                ibond,ibond_new,xibond,xibond_new,rotbond,rotbond_new,iorig,isdead_or_accreted,maxbond
 integer, intent(in) :: npart
 real,    intent(in) :: xyzh(:,:)
 integer, allocatable, save :: imap(:)
 integer(kind=8) :: jo,maxid
 integer :: i,j,k

 if (allocated(icontact)) then
    !$omp parallel do default(none) shared(npart,icontact,icontact_new,xicontact,xicontact_new) &
    !$omp shared(rotcontact,rotcontact_new) private(i)
    do i=1,npart
       icontact(:,i)     = icontact_new(:,i)
       xicontact(:,:,i)  = xicontact_new(:,:,i)
       rotcontact(:,:,i) = rotcontact_new(:,:,i)
    enddo
    !$omp end parallel do
 endif

 if (allocated(ibond)) then
    !$omp parallel do default(none) shared(npart,ibond,ibond_new,xibond,xibond_new) &
    !$omp shared(rotbond,rotbond_new) private(i)
    do i=1,npart
       ibond(:,i)     = ibond_new(:,i)
       xibond(:,:,i)  = xibond_new(:,:,i)
       rotbond(:,:,i) = rotbond_new(:,:,i)
    enddo
    !$omp end parallel do
    !
    ! keep only mutual bonds
    !
    maxid = maxval(iorig(1:npart))
    if (allocated(imap)) then
       if (size(imap,kind=8) < maxid) deallocate(imap)
    endif
    if (.not.allocated(imap)) allocate(imap(maxid))
    imap(1:maxid) = 0
    do i=1,npart
       if (.not.isdead_or_accreted(xyzh(4,i))) imap(iorig(i)) = i
    enddo
    do i=1,npart
       do k=1,maxbond
          jo = ibond(k,i)
          if (jo == 0) cycle
          j = 0
          if (jo <= maxid) j = imap(jo)
          if (j == 0) then
             ibond(k,i) = 0
          elseif (.not.any(ibond(:,j) == iorig(i))) then
             ibond(k,i) = 0
          endif
       enddo
    enddo
 endif

end subroutine dem_commit_history

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
    write(iprint,"(a,1pg12.4,a)") '   bond radius = ',bond_lambda,' x min(R_i,R_j)'
    if (bond_sigma_pa > 0.) write(iprint,"(a,1pg12.4,a)") '   tensile strength = ',bond_sigma_pa,' Pa'
    if (bond_tau_pa > 0.)   write(iprint,"(a,1pg12.4,a)") '   shear strength   = ',bond_tau_pa,' Pa'
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