!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! Cut a body shape out of a settled DEM packing.
!
! Second stage of the settle-crop-settle packing method: a loose cloud is
! condensed under self-gravity into a sphere (phantomsetup with pack_settle),
! that sphere is cropped to the shape here, and the cropped body is settled
! again to relieve the cut surface.
!
! The grains keep the radius they settled with. Only the particle mass is
! reset, so that the kept grains carry the whole body mass at the requested
! bulk density: the settling sphere holds more grains than the shape keeps,
! so leaving massoftype alone would leave the body far too light.
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: centreofmass, io, part, physcon, prompting, shape, units
!
 use part,         only:idem,igas,kill_particle,shuffle_part,iphase,iamtype,maxphase,maxp
 use prompting,    only:prompt
 use shape,        only:inside_shape_file,get_mesh_geometry
 use centreofmass, only:reset_centreofmass
 use units,        only:udist,umass,unit_density
 use physcon,      only:pi,km
 use io,           only:id,master,fatal

 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

 character(len=256) :: shapefile = 'apophis.shape'
 real :: rho_bulk_cgs = 2.2

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 logical, allocatable :: keep(:)
 real    :: vol_shape,r_circ,r_grain,rmax_guess,m_body,phi_settled,rbody
 integer :: i,nkeep,nkill,ierr,itype

 if (npart < 1) call fatal('moddump_cropshape','no particles in dump')

 call prompt('Shape file to cut the body out of',shapefile)
 call prompt('Bulk density of the cropped body in g/cm^3',rho_bulk_cgs,0.)
 !
 ! The mesh is centred on the origin by the shape loader, so the body has to
 ! be centred too or the cut lands in the wrong place. The settle is isolated
 ! (no sinks) so the centre of mass should not have moved, but a drift of even
 ! a grain radius would bias which grains survive, so do not assume it.
 !
 call reset_centreofmass(npart,xyzh,vxyzu)
 !
 ! grain radius and the settled radius, measured from the dump rather than
 ! assumed: every DEM grain carries its radius in h.
 !
 r_grain = 0.
 rbody   = 0.
 do i=1,npart
    if (xyzh(4,i) > r_grain) r_grain = xyzh(4,i)
    rbody = max(rbody,sqrt(dot_product(xyzh(1:3,i),xyzh(1:3,i))))
 enddo
 !
 ! exact shape volume, for the mass, and circumradius, to check the body
 ! actually encloses the shape before cutting
 !
 rmax_guess = rbody
 call get_mesh_geometry(shapefile,rmax_guess,vol_shape,r_circ,ierr,id,master)
 if (ierr /= 0) call fatal('moddump_cropshape','could not read shape file '//trim(shapefile))

 if (id==master) then
    print "(/,a)",' --- cropping a settled packing to shape ---'
    print "(a,1pg10.3,a)",' settled body radius  = ',rbody*udist/km,' km'
    print "(a,1pg10.3,a)",' shape circumradius   = ',r_circ*udist/km,' km'
    print "(a,1pg10.3,a)",' shape volume         = ',vol_shape*(udist/km)**3,' km^3'
    print "(a,1pg10.3,a)",' grain radius         = ',r_grain*udist/km,' km'
 endif

 if (r_circ > rbody) then
    print "(a)",' *** WARNING: the shape sticks out of the settled body ***'
    print "(a)",' *** the crop will truncate it. Settle a larger sphere. ***'
 endif
 !
 ! mark the grains inside the shape
 !
 allocate(keep(npart))
 call inside_shape_file(shapefile,rmax_guess,npart,xyzh,keep,nkeep,ierr,id,master)
 if (nkeep < 1) call fatal('moddump_cropshape','shape kept no particles')
 !
 ! kill the rest. kill_particle is not thread safe, so this stays serial,
 ! and shuffle_part then closes the gaps it leaves.
 !
 nkill = 0
 do i=1,npart
    if (.not.keep(i)) then
       call kill_particle(i,npartoftype)
       nkill = nkill + 1
    endif
 enddo
 call shuffle_part(npart)
 deallocate(keep)
 !
 ! the kept grains must carry the whole body mass, or the body ends up light
 ! by the ratio of the settling sphere to the shape
 !
 m_body = (rho_bulk_cgs/unit_density)*vol_shape
 itype  = idem
 if (npartoftype(idem) < 1) itype = igas   ! dump predates the DEM type
 massoftype(itype) = m_body/real(npart)
 !
 ! Re-centre. load_obj_mesh centres the mesh on its BOUNDING BOX, not on its
 ! centre of volume, and for the Apophis mesh those differ by ~20 m (3.5 grain
 ! radii at 10k). The cropped body therefore sits off-origin by that much, so
 ! put its own centre of mass back at the origin before anything downstream
 ! places it on an orbit or settles it again.
 !
 call reset_centreofmass(npart,xyzh,vxyzu)
 !
 ! packing fraction actually achieved, as a check on the settle
 !
 phi_settled = real(npart)*(4./3.*pi*r_grain**3)/vol_shape

 if (id==master) then
    print "(a,i9,a,i9)",   ' grains kept          = ',nkeep,'  of ',nkeep+nkill
    print "(a,2(es10.3,a))",' body mass            = ',m_body*umass,' g'
    print "(a,1pg10.3,a)", ' bulk density         = ',m_body/vol_shape*unit_density,' g/cm^3'
    print "(a,1pg10.3)",   ' packing fraction     = ',phi_settled
    print "(a)",' now settle again: the cut surface is not at force balance'
 endif

end subroutine modify_dump

end module moddump
