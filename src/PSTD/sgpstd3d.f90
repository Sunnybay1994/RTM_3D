!  sgpstd3d.f90 
!
!  FUNCTIONS:
!  sgpstd3d - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: sgpstd3d
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

program SGPSTD3D
    use mkl_dfti
    use omp_lib
    implicit none
    interface
        subroutine FFTDiffr(A, B, N, coef,My_Desc_Handle) 
            use MKL_DFTI
            type(DFTI_DESCRIPTOR), POINTER,intent(inout) :: My_Desc_Handle
	        integer, intent(in) ::  N
            real(4), intent(in)    ::  A(N)
            real(4), intent(out)   ::  B(N-1)
            real(4), intent(in)    ::  coef(N-1)
        end subroutine

        subroutine FFTDiffl(A, B, N,coef,My_Desc_Handle)
            use MKL_DFTI
            type(DFTI_DESCRIPTOR), POINTER,intent(inout) :: My_Desc_Handle
	        integer, intent(in) ::  N
            real(4), intent(in)    ::  A(N)
            real(4), intent(out)   ::  B(N-1)
            real(4), intent(in)    ::  coef(N-1)
	    end subroutine 
    end interface 
    
    interface data_output
        subroutine data_output_vector(field, field_name, Nx, Ny, Nz, t)
	        integer, intent(in) :: Nx, Ny, Nz
	        real(4), intent(in)    :: field(Nx, Ny, Nz)
	        integer, intent(in) :: t
	        character*11        :: field_name
        end subroutine data_output_vector
        
        subroutine data_output_scalar(field, field_name, t)
	        real(4), intent(in)    :: field
	        integer, intent(in) :: t
	        character*10        :: field_name
        end subroutine data_output_scalar
    end interface data_output

    
    real(4), parameter   :: pi       = 3.14159265358979
    real(4), parameter   :: epsl0    = 8.85*1e-12
    real(4), parameter   :: mu0      = 4*pi*1e-7
    
    integer     :: Nx, Ny, Nz
    integer     :: PMLx, PMLy, PMLz
    real(4)        :: PMLx_factor1, PMLy_factor1, PMLz_factor1
    real(4)        :: PMLx_factor2, PMLy_factor2, PMLz_factor2
    
    real(4)        :: dx,dy,dz
    
    !!================for mt_wash model
    !real(4),allocatable,dimension(:,:) :: tp
    !real(4),allocatable,dimension(:)   :: z
    !!=================
    
    real(4)     :: dt
    integer     :: Nt
    real(4)     :: travelTime
    
    real(4) :: fm, pulseWidth
	integer ::  nt_src, nsrc, nrec
	integer ::  nslicex, nslicey, nslicez
    integer,allocatable,dimension(:)      :: islicex, islicey, islicez
	integer ::  xstep, tstep
    integer, allocatable, dimension(:)	::  sourcex_idx, sourcey_idx, sourcez_idx, compnt
	integer, allocatable, dimension(:)  ::  recx,recy,recz,compnt_rec
	CHARACTER(2)	:: component
	real(4), allocatable, dimension(:,:)   :: srcpulse
	real(4), allocatable, dimension(:,:)   :: gather
    
    real(4), allocatable, dimension(:,:,:)  :: Ex, Ey, Ez
    real(4), allocatable, dimension(:,:,:)  :: Hx, Hy, Hz
    real(4), allocatable, dimension(:,:,:)  :: epsl, sigma
    real(4), allocatable, dimension(:,:,:)  :: Dxy, Dxz, Dyx, Dyz, Dzx, Dzy
    real(4) :: epslx, epsly, epslz
    real(4) :: sigmax, sigmay, sigmaz
    
    real(4), allocatable, dimension(:)      :: coefx, coefy, coefz
    
    real(4), allocatable, dimension(:,:)    :: tranx, trany, tranz
    real(4), allocatable, dimension(:,:)    :: dtranx, dtrany, dtranz
    real(4), allocatable, dimension(:,:,:)  :: field1, field2, field3, field4
    type df_pointer
        type(dfti_descriptor),pointer       :: desc_handle
    end type
    
    type(df_pointer), allocatable, dimension(:) :: df_pointer_x, df_pointer_y, df_pointer_z
    
    integer     :: nthreads
    integer     :: pid
    
    
    !The PML boundary
    real(4),allocatable ::Feyx(:,:,:)
    real(4),allocatable ::Feyz(:,:,:)
    real(4),allocatable ::Fezx(:,:,:)
    real(4),allocatable ::Fezy(:,:,:)
    real(4),allocatable ::Fexy(:,:,:)
    real(4),allocatable ::Fexz(:,:,:)

    real(4),allocatable ::Fhyx(:,:,:)
    real(4),allocatable ::Fhyz(:,:,:)
    real(4),allocatable ::Fhzx(:,:,:)
    real(4),allocatable ::Fhzy(:,:,:)
    real(4),allocatable ::Fhxy(:,:,:)
    real(4),allocatable ::Fhxz(:,:,:)

    real(4),allocatable ::iv(:)
    real(4),allocatable ::jv(:)
    real(4),allocatable ::kv(:)

    real(4),allocatable ::wx(:)
    real(4),allocatable ::wy(:)
    real(4),allocatable ::wz(:)
    real(4),allocatable ::wx2(:)
    real(4),allocatable ::wy2(:)
    real(4),allocatable ::wz2(:)
    
    real(4),allocatable ::ewx(:)
    real(4),allocatable ::ewy(:)
    real(4),allocatable ::ewz(:)
    real(4),allocatable ::ewx2(:)
    real(4),allocatable ::ewy2(:)
    real(4),allocatable ::ewz2(:)
    
    integer     :: i,j,k,m,n
    integer     :: status
    
    character(len=10)   :: field_name
    integer,allocatable,dimension(:,:,:)      :: RecordSurfaceEz, RecordSurfaceHx, RecordSurfaceHy
    real(4),allocatable,dimension(:,:,:)      :: EzSurface, HxSurface, HySurface
    real(4),allocatable,dimension(:,:,:)      :: slicex, slicey, slicez
    
    integer         :: time1, time2, time3, time0, time4
	real(4)         :: cal_time, io_time
	character(100)  :: string100
	character(15)   :: s1,s2,s3,s4,s5
    
    integer         :: isrc
	CHARACTER(4)    :: isrc_s
    CHARACTER(6)    :: it_s
    
	CHARACTER(9)    :: cmd_arg
    
    
    call system_clock(time0)

    
    if (command_argument_count() > 0) then
        CALL GET_COMMAND_ARGUMENT(1,cmd_arg)
        read(cmd_arg,*) nthreads
	else
		nthreads = 8
	endif
    print *, "nthreads:",nthreads

	if (command_argument_count() > 1) then
        CALL GET_COMMAND_ARGUMENT(2,cmd_arg)
        read(cmd_arg,*) isrc
        write(isrc_s,'(i4.4)') isrc
	else
		isrc_s = "0000"
    endif
    print *, "isrc:",isrc
    
    
    open(1,file = 'input/par.in',status = 'old',action = 'read')   
    read(1,*)
    read(1,*) dx,dy,dz,dt
    print *, "dx=",dx,", dy=",dy,", dz=",dz,", dt=",dt,'\n'
    read(1,*)
    read(1,*) nx, ny, nz, nt
    print *, "nx=",nx,", ny=",ny,", nz=",nz,", nt=",nt,'\n'
    read(1,*)
    read(1,*)
    read(1,*)
    read(1,*) tstep, xstep
    print *, "tstep=",tstep,", xstep=",xstep,'\n'
    read(1,*)
    read(1,*)
    read(1,*)
    read(1,*) PMLx, PMLy, PMLz
    print *, "PMLx=",PMLx,", PMLy=",PMLy,", PMLz=",PMLz,'\n'
    ! read(1,*)
    ! read(1,*)
    ! read(1,*)
    ! read(1,*)
    ! read(1,*)
    ! read(1,*)
    ! read(1,*)
    ! read(1,*)
    ! read(1,*)
    ! read(1,*)
    ! read(1,*)
    ! read(1,*)
    ! read(1,*)
    ! read(1,*)
    close(1)
	
    open(6,file = 'input/slice.in',status = 'old',action = 'read')
    read(6,*) nslicex,nslicey,nslicez
    print *, "nslicex=",nslicex,", nslicey=",nslicey,", nslicez=",nslicez,'\n'
    allocate(islicex(nslicex),islicey(nslicey),islicez(nslicez))
    do i = 1,nslicex
        read(6,*) islicex(i),cmd_arg
    enddo
    do i = 1,nslicey
        read(6,*) islicey(i),cmd_arg
    enddo
    do i = 1,nslicez
        read(6,*) islicez(i),cmd_arg
    enddo
    print *, "islicex=",islicex
    print *, "islicey=",islicey
    print *, "islicez=",islicez
    close(6)
    
	
    !! for mt_wash model
    !allocate(tp(Nx,Ny))
    !allocate(z(Nz))
    !open(1,file = 'MWashDEM1024.txt',status = 'old')
    !do j = Ny,1,-1
    !    read(1,*) tp(:,j)
    !enddo
    !do k = 1,Nz
    !    z(k) = k*dz
    !enddo
    !z = z-z(Nz/2)+(minval(tp)+maxval(tp))/2.
    
!============================
 !   dt = 2./pi/sqrt(3.)*min(dx,dy,dz)*sqrt(mu0*epsl0)
 !   pulseWidth = 100
 !   fm = 2/pulseWidth/dt
 !   travelTime = -1./fm
	!print *, 'frq_max = ', fm, 'dt = ', dt
	
	allocate(Ex(Nx,Ny,Nz), Ey(Nx,Ny,Nz), Ez(Nx,Ny,Nz))
    allocate(Hx(Nx,Ny,Nz), Hy(Nx,Ny,Nz), Hz(Nx,Ny,Nz))
    allocate(epsl(Nx,Ny,Nz), sigma(Nx,Ny,Nz))
    
    !epsl    = epsl0
    !sigma   = 0.
	
    ! load 3D model
	open(2,file = 'input/eps.in',status = 'old')
	open(3,file = 'input/sig.in',status = 'old')
	do i = 1,nx
		do j = 1,ny
			!print *,i*100+j
			read(2,*) epsl(i,j,:)
			read(3,*) sigma(i,j,:)
		enddo
	enddo
	close(2)
	close(3)
	epsl = epsl * epsl0
	
	! load source
	open(4,file = 'input/src.in_'//isrc_s,status = 'old',action = 'read')   
    read(4,*), nsrc, nt_src
	allocate(sourcex_idx(nsrc), sourcey_idx(nsrc), sourcez_idx(nsrc),compnt(nsrc))
	allocate(srcpulse(nsrc,nt_src))
	do i = 1,nsrc
		read(4,*), sourcex_idx(i), sourcey_idx(i), sourcez_idx(i), component
		if (component == 'Ex') then
			compnt(i) = 1
		else if (component == 'Ey') then
			compnt(i) = 2
		else if (component == 'Ez') then
			compnt(i) = 3
		else if (component == 'Hx') then
			compnt(i) = 4
		else if (component == 'Hy') then
			compnt(i) = 5
		else if (component == 'Hz') then
			compnt(i) = 6
		else
			print *, 'unknown component:'//component
		endif
	enddo
	do i = 1,nsrc
		read(4,*), srcpulse(i,:)
	    !do j = 1,nt_src
	    !    srcpulse(i,j) = 2.*pi*fm*(j*dt-1./fm)*exp(-(pi*fm*(j*dt-1./fm))**2)
	    !enddo
	enddo
    close(4)
	
	
	open(5,file = 'input/rec.in',status = 'old',action = 'read')
	read(5,*), nrec
	allocate(recx(nrec), recy(nrec), recz(nrec),compnt_rec(nrec))
	allocate(gather(nt,nrec))
	do i = 1,nrec
		!print *,i
		read(5,*), recx(i), recy(i), recz(i), component
		if (component == 'Ex') then
			compnt_rec(i) = 1
		else if (component == 'Ey') then
			compnt_rec(i) = 2
		else if (component == 'Ez') then
			compnt_rec(i) = 3
		else if (component == 'Hx') then
			compnt_rec(i) = 4
		else if (component == 'Hy') then
			compnt_rec(i) = 5
		else if (component == 'Hz') then
			compnt_rec(i) = 6
		else
			print *, 'unknown component:'//component
		endif
	enddo
	
	
	
    !!for mt_wash model
    !do k = 1,Nz-1
    !    do j = 1,Ny-1
    !        do i = 1,Nx-1
    !            if((z(k+1)-dz/2.)<(tp(i,j)+tp(i,j+1)+tp(i+1,j)+tp(i+1,j+1))/4.) then
    !                epsl(i,j,k) = 4.*epsl0
    !                sigma(i,j,k) = 1e-3
    !            endif
    !        enddo
    !    enddo
    !enddo
    
	!! record surface 
    !allocate(slicex(Ny,Nz),slicey(Ny,Nz),slicez(Nx,Ny))
 !   allocate(RecordSurfaceEz(Nx,Ny),RecordSurfaceHx(Nx,Ny),RecordSurfaceHy(Nx,Ny))
 !   allocate(EzSurface(Nx,Ny),HxSurface(Nx,Ny),HySurface(Nx,Ny))
 !   RecordSurfaceEz = Nz/2
 !   RecordSurfaceHx = Nz/2
 !   RecordSurfaceHy = Nz/2
 !   EzSurface = 0.
 !   HxSurface = 0.
 !   HySurface = 0.
 !   
 !
 !   do j = 2,Ny-1
 !       do i = 2,Nx-1
 !           do k = 1,Nz-1
 !               if(max(epsl(i,j,k),epsl(i-1,j,k),epsl(i,j-1,k),epsl(i-1,j-1,k))<3.*epsl0) then
 !                   RecordSurfaceEz(i,j) = k
 !                   exit
 !               endif
 !           enddo
 !       enddo
 !   enddo
 !   
 !   do j = 1,Ny-1
 !       do i = 2,Nx-1
 !           do k = 1,Nz-1
 !               if(max(epsl(i,j,k),epsl(i-1,j,k))<3.*epsl0) then
 !                   RecordSurfaceHx(i,j) = k
 !                   exit
 !               endif
 !           enddo
 !       enddo
 !   enddo
 !
 !   do j = 2,Ny-1
 !       do i = 1,Nx-1
 !           do k = 1,Nz-1
 !               if(max(epsl(i,j,k),epsl(i,j-1,k))<3.*epsl0) then
 !                   RecordSurfaceHy(i,j) = k
 !                   exit
 !               endif
 !           enddo
 !       enddo
 !   enddo
    
	!! source 
 !   sourcex_idx = 62
 !   sourcey_idx = 1024-143+1
 !   i = sourcex_idx
 !   j = sourcey_idx
 !   do k = 1,Nz-1
 !       if(max(epsl(i,j,k),epsl(i-1,j,k),epsl(i,j-1,k),epsl(i-1,j-1,k))<3.*epsl0) then
 !           sourcez_idx = k+1
 !           exit
 !       endif
 !   enddo
 !   
 !   print*,'sourcez_idx = ',sourcez_idx
       
    
    Ex      = 0.
    Ey      = 0.
    Ez      = 0.
    Hx      = 0.
    Hy      = 0.
    Hz      = 0.
    
    allocate(Dxy(Nx,Ny,Nz),Dxz(Nx,Ny,Nz))
    allocate(Dyx(Nx,Ny,Nz),Dyz(Nx,Ny,Nz))
    allocate(Dzx(Nx,Ny,Nz),Dzy(Nx,Ny,Nz))
    
    allocate(Feyx(1:2*PMLx,1:Ny,1:Nz))
    allocate(Fezx(1:2*PMLx,1:Ny,1:Nz))
    allocate(Fexy(1:Nx,1:2*PMLy,1:Nz))
    allocate(Fezy(1:Nx,1:2*PMLy,1:Nz))
    allocate(Fexz(1:Nx,1:Ny,1:2*PMLz))
    allocate(Feyz(1:Nx,1:Ny,1:2*PMLz))

    allocate(Fhyx(1:2*PMLx,1:Ny,1:Nz))
    allocate(Fhzx(1:2*PMLx,1:Ny,1:Nz))
    allocate(Fhxy(1:Nx,1:2*PMLy,1:Nz))
    allocate(Fhzy(1:Nx,1:2*PMLy,1:Nz))
    allocate(Fhxz(1:Nx,1:Ny,1:2*PMLz))
    allocate(Fhyz(1:Nx,1:Ny,1:2*PMLz))

    allocate(wx(1:Nx))
    allocate(wy(1:Ny))
    allocate(wz(1:Nz))
    allocate(wx2(1:Nx))
    allocate(wy2(1:Ny))
    allocate(wz2(1:Nz))

    allocate(ewx(1:Nx))
    allocate(ewy(1:Ny))
    allocate(ewz(1:Nz))
    allocate(ewx2(1:Nx))
    allocate(ewy2(1:Ny))
    allocate(ewz2(1:Nz))
    
    allocate(iv(1:Nx))
    allocate(jv(1:Ny))
    allocate(kv(1:Nz))
    
    wx  = 0.
    wy  = 0.
    wz  = 0.
    wx2 = 0.
    wy2 = 0.
    wz2 = 0.
    PMLx_factor1 = 10./PMLx**3./dx/(epsl0*mu0)**0.5
    PMLx_factor2 = 10./PMLx**3./dx/(epsl0*mu0)**0.5
    PMLy_factor1 = 10./PMLy**3./dy/(epsl0*mu0)**0.5
    PMLy_factor2 = 10./PMLy**3./dy/(epsl0*mu0)**0.5
    PMLz_factor1 = 10./PMLz**3./dz/(epsl0*mu0)**0.5
    PMLz_factor2 = 10./PMLz**3./dz/(4*epsl0*mu0)**0.5
    do i = 1,PMLx
        wx2(i)      = PMLx_factor1*(PMLx-i)**2.
        wx(i)       = PMLx_factor1*(PMLx-i+0.5)**2.
        wx2(Nx-i+1) = PMLx_factor2*(PMLx-i+1)**2.
        wx(Nx-i+1)  = PMLx_factor2*(PMLx-i+0.5)**2.
    enddo
    do i = 1,PMLy
        wy2(i)      = PMLy_factor1*(PMLy-i)**2.
        wy(i)       = PMLy_factor1*(PMLy-i+0.5)**2.
        wy2(Ny-i+1) = PMLy_factor2*(PMLy-i+1)**2.
        wy(Ny-i+1)  = PMLy_factor2*(PMLy-i+0.5)**2.
    enddo
    do i = 1,PMLz
        wz2(i)      = PMLz_factor1*(PMLz-i)**2.
        wz(i)       = PMLz_factor1*(PMLz-i+0.5)**2.
        wz2(Nz-i+1) = PMLz_factor2*(PMLz-i+1)**2.
        wz(Nz-i+1)  = PMLz_factor2*(PMLz-i+0.5)**2.
    end do

    do i = 1,Nx
         ewx(i)     = exp((-wx(i)*dt))
         ewx2(i)    = exp((-wx2(i)*dt))
    end do

    do j = 1,Ny
         ewy(j)     = exp((-wy(j)*dt))
         ewy2(j)    = exp((-wy2(j)*dt))
    end do

    do k = 1,Nz
         ewz(k)     = exp((-wz(k)*dt))
         ewz2(k)    = exp((-wz2(k)*dt))
    end do
    
    do i = 1,PMLx
        iv(i) = i
        iv(Nx-PMLx+i) = PMLx+i
    enddo
    do i = 1,PMLy
        jv(i) = i
        jv(Ny-PMLy+i) = PMLy+i
    enddo
    do i = 1,PMLz
        kv(i) = i;
        kv(Nz-PMLz+i) = PMLz+i;
    end do
    
    
    allocate(coefx(Nx-1),coefy(Ny-1),coefz(Nz-1))
    do i = 1,Nx/2-1
        coefx(i) = 2*PI/Nx/dx*(i)*cos(PI*(i)/Nx)
        coefx(Nx/2-1+i) = 2*PI/Nx/dx*(i)*sin(PI*(i)/Nx)
    end do
    coefx(Nx-1) = PI/Nx/dx*Nx

    do i = 1,Ny/2-1
        coefy(i) = 2*PI/Ny/dy*(i)*cos(PI*(i)/Ny)
        coefy(Ny/2-1+i) = 2*PI/Ny/dy*(i)*sin(PI*(i)/Ny)
    end do
    coefy(Ny-1) = PI/Ny/dy*Ny

    do i = 1,Nz/2-1
        coefz(i) = 2*PI/Nz/dz*(i)*cos(PI*(i)/Nz)
        coefz(Nz/2-1+i) = 2*PI/Nz/dz*(i)*sin(PI*(i)/Nz)
    end do
    coefz(Nz-1) = PI/Nz/dz*Nz
    
    allocate(df_pointer_x(nthreads), df_pointer_y(nthreads), df_pointer_z(nthreads))
    
    do i = 1,nthreads
        Status = DftiCreateDescriptor(df_pointer_x(i)%desc_handle, dfti_single, dfti_real, 1, Nx)
        Status = DftiCommitDescriptor(df_pointer_x(i)%desc_handle)
        Status = DftiCreateDescriptor(df_pointer_y(i)%desc_handle, dfti_single, dfti_real, 1, Ny)
        Status = DftiCommitDescriptor(df_pointer_y(i)%desc_handle)
        Status = DftiCreateDescriptor(df_pointer_z(i)%desc_handle, dfti_single, dfti_real, 1, Nz)
        Status = DftiCommitDescriptor(df_pointer_z(i)%desc_handle)
    enddo
    
    allocate(tranx(Nx,nthreads), dtranx(Nx-1,nthreads))
    allocate(trany(Ny,nthreads), dtrany(Ny-1,nthreads))
    allocate(tranz(Nz,nthreads), dtranz(Nz-1,nthreads))
    allocate(field1(Nz,Nx,nthreads), field2(Nz,Nx,nthreads))
    allocate(field3(Ny,Nx,nthreads), field4(Ny,Nx,nthreads))
    
    !open(997,file = 'HxSurface.bin',access = 'stream', form = 'unformatted',status = 'replace')
    !open(998,file = 'HySurface.bin',access = 'stream', form = 'unformatted',status = 'replace')
    !open(999,file = 'EzSurface.bin',access = 'stream', form = 'unformatted',status = 'replace')
    !open(996,file = 'travelTimeE.txt',status = 'replace')
    !!open(995,file = 'travelTimeH.txt',status = 'replace')
	
    open(11,file = 'output/time_usage.txt',status = 'replace')

    print*, 'nthreads = ',nthreads
	
	cal_time = 0
	io_time = 0
    
    do n = 1,Nt
		call system_clock(time1)
        !print*, 'Time step: ', n
            
        !$omp parallel private(pid) num_threads(nthreads) 
        pid = omp_get_thread_num()
        !$omp do
        do j = 1, Ny
	        field1(:,:,pid+1) = transpose(Hx(:,j,:))
	        field2(:,:,pid+1) = transpose(Hy(:,j,:))
		    do i = 1, Nx
                if(i>1) then
                    tranz(:,pid+1) = field1(:,i,pid+1)
			        call fftdiffr(tranz(:,pid+1), dtranz(:,pid+1), Nz, coefz, df_pointer_z(pid+1)%desc_handle) 
                    field1(1:Nz-1,i,pid+1) = dtranz(:,pid+1)
                endif
                if(j>1) then
			        tranz(:,pid+1) = field2(:,i,pid+1)
			        call fftdiffr(tranz(:,pid+1), dtranz(:,pid+1), Nz, coefz, df_pointer_z(pid+1)%desc_handle) 
                    field2(1:Nz-1,i,pid+1) = dtranz(:,pid+1)
                endif
		    end do
		    Dxz(:,j,:) = transpose(field1(:,:,pid+1))
		    Dyz(:,j,:) = transpose(field2(:,:,pid+1))
        end do
        !$omp end do

        !$omp do
	    Do k = 1, Nz
		    Do j= 1, Ny
                if(k>1) then
			        call fftdiffr(Hz(:,j,k), Dzx(1:Nx-1,j,k), Nx, coefx, df_pointer_x(pid+1)%desc_handle) 
                endif
                
			    if(j>1) then
			        call fftdiffr(Hy(:,j,k), Dyx(1:Nx-1,j,k), Nx, coefx, df_pointer_x(pid+1)%desc_handle)
                endif
                
		    end do
		
        end do
        !$omp end do

        
        !$omp do
	    do k = 1, Nz
            field3(:,:,pid+1) = transpose(Hz(:,:,k))
            field4(:,:,pid+1) = transpose(Hx(:,:,k))
		    Do i= 1, Nx
                if(k>1) then
			        trany(:,pid+1) = field3(:,i,pid+1)
			        call fftdiffl(trany(:,pid+1), dtrany(:,pid+1), Ny, coefy, df_pointer_y(pid+1)%desc_handle)
                    field3(1:Ny-1,i,pid+1) = dtrany(:,pid+1)
                endif
                
			    if(i>1) then
			        trany(:,pid+1) = field4(:,i,pid+1)
			        call fftdiffl(trany(:,pid+1), dtrany(:,pid+1), Ny, coefy, df_pointer_y(pid+1)%desc_handle)
                    field4(1:Ny-1,i,pid+1) = dtrany(:,pid+1)
                endif
            end do
            Dzy(:,:,k) = transpose(field3(:,:,pid+1))
            Dxy(:,:,k) = transpose(field4(:,:,pid+1))
        end do
        !$omp end do
        !$omp end parallel

        !$omp parallel do private(epslx, epsly, epslz, sigmax, sigmay, sigmaz) num_threads(nthreads)
        do k = 1,Nz-1
            do j = 1,Ny-1
                do i = 1,Nx-1
            
                    if(i<=PMLx .or. i>Nx-PMLx) then
                        if(i>1 .and. j>1) then
                            Fhyx(iv(i),j,k) = Fhyx(iv(i),j,k)+ewx(i)*(Dyx(i-1,j,k)/2)-Dyx(i-1,j,k)/2
                        endif
                        if(i>1 .and. k>1) then
                            Fhzx(iv(i),j,k) = Fhzx(iv(i),j,k)+ewx(i)*(Dzx(i-1,j,k)/2)-Dzx(i-1,j,k)/2
                        endif
                    end if
                    if(j<=PMLy .or. j>Ny-PMLy) then
                        if(i>1 .and. j>1) then
                            Fhxy(i,jv(j),k) = Fhxy(i,jv(j),k)+ewy(j)*(Dxy(i,j-1,k)/2)-Dxy(i,j-1,k)/2
                        endif
                        if(j>1 .and. k>1) then
                            Fhzy(i,jv(j),k) = Fhzy(i,jv(j),k)+ewy(j)*(Dzy(i,j-1,k)/2)-Dzy(i,j-1,k)/2
                        endif
                    end if
                    if(k<=PMLz .or. k>Nz-PMLz) then
                        if(i>1 .and. k>1) then
                            Fhxz(i,j,kv(k)) = Fhxz(i,j,kv(k))+ewz(k)*(Dxz(i,j,k-1)/2)-Dxz(i,j,k-1)/2
                        endif
                        if(j>1 .and. k>1) then
                            Fhyz(i,j,kv(k)) = Fhyz(i,j,kv(k))+ewz(k)*(Dyz(i,j,k-1)/2)-Dyz(i,j,k-1)/2
                        endif
                    end if
                
                    if(j>1 .and. k>1) then
                        epslx   = epsl(i,j,k)+epsl(i,j-1,k)+epsl(i,j,k-1)+epsl(i,j-1,k-1)
                        sigmax  = sigma(i,j,k)+sigma(i,j-1,k)+sigma(i,j,k-1)+sigma(i,j-1,k-1)
                        Ex(i,j,k) = (dt*(Dzy(i,j-1,k)-Dyz(i,j,k-1))+(epslx-sigmax*dt/2)*Ex(i,j,k))/(epslx+sigmax*dt/2)
                        
                        if(j<=PMLy .or. j>Ny-PMLy) then                 
                            Ex(i,j,k) = Ex(i,j,k) + dt*(Fhzy(i,jv(j),k))/(epslx+sigmax*dt/2)
                            Fhzy(i,jv(j),k) = ewy(j)*(Fhzy(i,jv(j),k)+Dzy(i,j-1,k)/2)-Dzy(i,j-1,k)/2             
                        end if
                        if(k<=PMLz .or. k>Nz-PMLz) then                 
                            Ex(i,j,k) = Ex(i,j,k) + dt*(-Fhyz(i,j,kv(k)))/(epslx+sigmax*dt/2)  
                            Fhyz(i,j,kv(k)) = ewz(k)*(Fhyz(i,j,kv(k))+Dyz(i,j,k-1)/2)-Dyz(i,j,k-1)/2            
                        end if
                    endif
                    
                    if(i>1 .and. k>1) then
                        epsly   = epsl(i,j,k)+epsl(i-1,j,k)+epsl(i,j,k-1)+epsl(i-1,j,k-1)
                        sigmay  = sigma(i,j,k)+sigma(i-1,j,k)+sigma(i,j,k-1)+sigma(i-1,j,k-1)
                        Ey(i,j,k) = (dt*(Dxz(i,j,k-1)-Dzx(i-1,j,k))+(epsly-sigmay*dt/2)*Ey(i,j,k))/(epsly+sigmay*dt/2)
                        
                        if(i<=PMLx .or. i>Nx-PMLx) then       
                            Ey(i,j,k) = Ey(i,j,k) + dt*(-Fhzx(iv(i),j,k))/(epsly+sigmay*dt/2)
                            Fhzx(iv(i),j,k) = ewx(i)*(Fhzx(iv(i),j,k)+Dzx(i-1,j,k)/2)-Dzx(i-1,j,k)/2
                        end if
                        if(k<=PMLz .or. k>Nz-PMLz) then 
                            Ey(i,j,k) = Ey(i,j,k) + dt*(Fhxz(i,j,kv(k)))/(epsly+sigmay*dt/2)
                            Fhxz(i,j,kv(k)) = ewz(k)*(Fhxz(i,j,kv(k))+Dxz(i,j,k-1)/2)-Dxz(i,j,k-1)/2             
                        end if
                    endif
                    
                    if(i>1 .and. j>1) then
                        epslz   = epsl(i,j,k)+epsl(i-1,j,k)+epsl(i,j-1,k)+epsl(i-1,j-1,k)
                        sigmaz  = sigma(i,j,k)+sigma(i-1,j,k)+sigma(i,j-1,k)+sigma(i-1,j-1,k)
                        Ez(i,j,k) = (dt*(Dyx(i-1,j,k)-Dxy(i,j-1,k))+(epslz-sigmaz*dt/2)*Ez(i,j,k))/(epslz+sigmaz*dt/2)
                        
                        if(i<=PMLx .or. i>Nx-PMLx) then
                            Ez(i,j,k) = Ez(i,j,k) + dt*(Fhyx(iv(i),j,k))/(epslz+sigmaz*dt/2)
                            Fhyx(iv(i),j,k) = ewx(i)*(Fhyx(iv(i),j,k)+Dyx(i-1,j,k)/2)-Dyx(i-1,j,k)/2
                        end if
                
                        if(j<=PMLy .or. j>Ny-PMLy) then
                            Ez(i,j,k) = Ez(i,j,k) + dt*(-Fhxy(i,jv(j),k))/(epslz+sigmaz*dt/2)
                            Fhxy(i,jv(j),k) = ewy(j)*(Fhxy(i,jv(j),k)+Dxy(i,j-1,k)/2)-Dxy(i,j-1,k)/2              
                        end if
                    endif
                end do
            end do
        end do
        !$omp end parallel do
        
        travelTime = travelTime+dt/2.
		!Ez(Nx/2:Nx/2+1,Ny/2:Ny/2+1,Nz/2:Nz/2+1) = Ez(Nx/2:Nx/2+1,Ny/2:Ny/2+1,Nz/2:Nz/2+1) + (1-2*(pi*fm*(n*dt-1/fm))**2)*exp(-(pi*fm*(n*dt-1/fm))**2)
        
        !$omp parallel private(pid) num_threads(nthreads)
        pid = omp_get_thread_num()
        !$omp do
        do j = 1, Ny
	        field1(:,:,pid+1) = transpose(Ex(:,j,:))
	        field2(:,:,pid+1) = transpose(Ey(:,j,:))
		    do i = 1, Nx
                if(j>1) then
			        tranz(:,pid+1) = field1(:,i,pid+1)
			        call fftdiffr(tranz(:,pid+1), dtranz(:,pid+1), Nz, coefz, df_pointer_z(pid+1)%desc_handle) 
                    field1(1:Nz-1,i,pid+1) = dtranz(:,pid+1)
                endif
                if(i>1) then
			        tranz(:,pid+1) = field2(:,i,pid+1)
			        call fftdiffr(tranz(:,pid+1), dtranz(:,pid+1), Nz, coefz, df_pointer_z(pid+1)%desc_handle) 
                    field2(1:Nz-1,i,pid+1) = dtranz(:,pid+1)
                endif
		    end do
		    Dxz(:,j,:) = transpose(field1(:,:,pid+1))
		    Dyz(:,j,:) = transpose(field2(:,:,pid+1))
        end do
        !$omp end do
        

        !$omp do
	    Do k = 1, Nz
		    Do j= 1, Ny
                if(j>1) then
			        call fftdiffl(Ez(:,j,k), Dzx(1:Nx-1,j,k), Nx, coefx, df_pointer_x(pid+1)%desc_handle)
                endif
                
			    if(k>1) then
			        call fftdiffl(Ey(:,j,k), Dyx(1:Nx-1,j,k), Nx, coefx, df_pointer_x(pid+1)%desc_handle)
                endif
                
		    end do
		
	    end do
        !$omp end do

        !$omp do
	    do k = 1, Nz
            field3(:,:,pid+1) = transpose(Ez(:,:,k))
            field4(:,:,pid+1) = transpose(Ex(:,:,k))
		    Do i= 1, Nx
                if(i>1) then
                    trany(:,pid+1) = field3(:,i,pid+1)
			        call fftdiffl(trany(:,pid+1), dtrany(:,pid+1), Ny, coefy, df_pointer_y(pid+1)%desc_handle)
                    field3(1:Ny-1,i,pid+1) = dtrany(:,pid+1)
                endif
                
			    if(k>1) then
			        trany(:,pid+1) = field4(:,i,pid+1)
			        call fftdiffl(trany(:,pid+1), dtrany(:,pid+1), Ny, coefy, df_pointer_y(pid+1)%desc_handle)
                    field4(1:Ny-1,i,pid+1) = dtrany(:,pid+1)
                endif
            end do
            Dzy(:,:,k) = transpose(field3(:,:,pid+1))
            Dxy(:,:,k) = transpose(field4(:,:,pid+1))
        end do
        !$omp end do
        !$omp end parallel
        
        !$omp parallel do num_threads(nthreads)
        do k = 1,Nz-1
            do j = 1,Ny-1
                do i = 1,Nx-1
                    
                    if(i<=PMLx .or. i>Nx-PMLx) then    
                        if(k>1) then
                            Feyx(iv(i),j,k) = Feyx(iv(i),j,k)+ewx2(i)*(Dyx(i,j,k)/2)-Dyx(i,j,k)/2
                        endif
                        if(j>1) then
                            Fezx(iv(i),j,k) = Fezx(iv(i),j,k)+ewx2(i)*(Dzx(i,j,k)/2)-Dzx(i,j,k)/2
                        endif
                    end if
                
                    if(j<=PMLy .or. j>Ny-PMLy) then 
                        if(k>1) then
                            Fexy(i,jv(j),k) = Fexy(i,jv(j),k)+ewy2(j)*(Dxy(i,j,k)/2)-Dxy(i,j,k)/2
                        endif
                        if(i>1) then
                            Fezy(i,jv(j),k) = Fezy(i,jv(j),k)+ewy2(j)*(Dzy(i,j,k)/2)-Dzy(i,j,k)/2
                        endif
                    end if
                
                    if(k<=PMLz .or. k>Nz-PMLz) then      
                        if(i>1) then
                            Feyz(i,j,kv(k)) = Feyz(i,j,kv(k))+ewz2(k)*(Dyz(i,j,k)/2)-Dyz(i,j,k)/2
                        endif
                        if(j>1) then
                            Fexz(i,j,kv(k)) = Fexz(i,j,kv(k))+ewz2(k)*(Dxz(i,j,k)/2)-Dxz(i,j,k)/2
                        endif
                    end if   
                    
                    if(i>1) then
                        Hx(i,j,k) = (-dt*(Dzy(i,j,k)-Dyz(i,j,k))/mu0+Hx(i,j,k))
                    endif
                    if(j>1) then
                        Hy(i,j,k) = (-dt*(Dxz(i,j,k)-Dzx(i,j,k))/mu0+Hy(i,j,k))
                    endif
                    if(k>1) then
                        Hz(i,j,k) = (-dt*(Dyx(i,j,k)-Dxy(i,j,k))/mu0+Hz(i,j,k))
                    endif
                    
                
                    if(i<=PMLx .or. i>Nx-PMLx) then   
                        if(j>1) then
                            Hy(i,j,k) = (-dt*(-Fezx(iv(i),j,k))/mu0+Hy(i,j,k))
                            Fezx(iv(i),j,k) = ewx2(i)*(Fezx(iv(i),j,k)+Dzx(i,j,k)/2)-Dzx(i,j,k)/2
                        endif
                        if(k>1) then
                            Hz(i,j,k) = (-dt*(Feyx(iv(i),j,k))/mu0+Hz(i,j,k))
                            Feyx(iv(i),j,k) = ewx2(i)*(Feyx(iv(i),j,k)+Dyx(i,j,k)/2)-Dyx(i,j,k)/2
                        endif
                    end if
                
                    if(j<=PMLy .or. j>Ny-PMLy) then  
                        if(i>1) then
                            Hx(i,j,k) = (-dt*(Fezy(i,jv(j),k))/mu0+Hx(i,j,k))
                            Fezy(i,jv(j),k) = ewy2(j)*(Fezy(i,jv(j),k)+Dzy(i,j,k)/2)-Dzy(i,j,k)/2
                        endif
                        if(k>1) then
                            Hz(i,j,k) = (-dt*(-Fexy(i,jv(j),k))/mu0+Hz(i,j,k))
                            Fexy(i,jv(j),k) = ewy2(j)*(Fexy(i,jv(j),k)+Dxy(i,j,k)/2)-Dxy(i,j,k)/2
                        endif
                    end if
                
                    if(k<=PMLz .or. k>Nz-PMLz) then   
                        if(i>1) then
                            Hx(i,j,k) = (-dt*(-Feyz(i,j,kv(k)))/mu0+Hx(i,j,k))
                            Feyz(i,j,kv(k)) = ewz2(k)*(Feyz(i,j,kv(k))+Dyz(i,j,k)/2)-Dyz(i,j,k)/2
                        endif
                        if(j>1) then
                            Hy(i,j,k) = (-dt*(Fexz(i,j,kv(k)))/mu0+Hy(i,j,k))
                            Fexz(i,j,kv(k)) = ewz2(k)*(Fexz(i,j,kv(k))+Dxz(i,j,k)/2)-Dxz(i,j,k)/2
                        endif
                    endif
                   
                end do
            end do
		end do
        !$omp end parallel do
        
		! record gather
		do i = 1,nrec
		    select case (compnt_rec(i))
			    case (1)
				    gather(n,i) = Ex(recx(i),recy(i),recz(i))
			    case (2)
				    gather(n,i) = Ey(recx(i),recy(i),recz(i))
			    case (3)
				    gather(n,i) = Ez(recx(i),recy(i),recz(i))
			    case (4)
				    gather(n,i) = Hx(recx(i),recy(i),recz(i))
			    case (5)
				    gather(n,i) = Hy(recx(i),recy(i),recz(i))
			    case (6)
				    gather(n,i) = Hz(recx(i),recy(i),recz(i))
			endselect
		enddo
		
		! input source
        !Hx(sourcex_idx,sourcey_idx-1,sourcez_idx:sourcez_idx+1)    = Hx(sourcex_idx,sourcey_idx-1,sourcez_idx:sourcez_idx+1) + 2.*pi*fm*travelTime*exp(-(pi*fm*travelTime)**2)
        !Hx(sourcex_idx,sourcey_idx,sourcez_idx:sourcez_idx+1)      = Hx(sourcex_idx,sourcey_idx,sourcez_idx:sourcez_idx+1) - 2.*pi*fm*travelTime*exp(-(pi*fm*travelTime)**2)
        !Hy(sourcex_idx-1,sourcey_idx,sourcez_idx:sourcez_idx+1)    = Hy(sourcex_idx-1,sourcey_idx,sourcez_idx:sourcez_idx+1) - 2.*pi*fm*travelTime*exp(-(pi*fm*travelTime)**2)
        !Hy(sourcex_idx,sourcey_idx,sourcez_idx:sourcez_idx+1)      = Hy(sourcex_idx,sourcey_idx,sourcez_idx:sourcez_idx+1) + 2.*pi*fm*travelTime*exp(-(pi*fm*travelTime)**2)
		if (n <= nt_src) then
			do i = 1,nsrc
		        select case (compnt(i))
			        case (1)
				        Ex(sourcex_idx(i),sourcey_idx(i),sourcez_idx(i))    = Ex(sourcex_idx(i),sourcey_idx(i),sourcez_idx(i)) + srcpulse(i,n)
					case (2)
						Ey(sourcex_idx(i),sourcey_idx(i),sourcez_idx(i))    = Ey(sourcex_idx(i),sourcey_idx(i),sourcez_idx(i)) + srcpulse(i,n)
				  !      Ey(sourcex_idx(i)-1,sourcey_idx(i),sourcez_idx(i):sourcez_idx(i)+1)    = Ey(sourcex_idx(i)-1,sourcey_idx(i),sourcez_idx(i):sourcez_idx(i)+1) + srcpulse(i,n)
						!Ey(sourcex_idx(i),sourcey_idx(i),sourcez_idx(i):sourcez_idx(i)+1)    = Ey(sourcex_idx(i),sourcey_idx(i),sourcez_idx(i):sourcez_idx(i)+1) - srcpulse(i,n)
			        case (3)
				        Ez(sourcex_idx(i),sourcey_idx(i),sourcez_idx(i))    = Ez(sourcex_idx(i),sourcey_idx(i),sourcez_idx(i)) + srcpulse(i,n)
			        case (4)
				        Hx(sourcex_idx(i),sourcey_idx(i),sourcez_idx(i))    = Hx(sourcex_idx(i),sourcey_idx(i),sourcez_idx(i)) + srcpulse(i,n)
			        case (5)
				        Hy(sourcex_idx(i),sourcey_idx(i),sourcez_idx(i))    = Hy(sourcex_idx(i),sourcey_idx(i),sourcez_idx(i)) + srcpulse(i,n)
			        case (6)
				        Hz(sourcex_idx(i),sourcey_idx(i),sourcez_idx(i))    = Hz(sourcex_idx(i),sourcey_idx(i),sourcez_idx(i)) + srcpulse(i,n)
					endselect
			enddo
		endif
		
        
        travelTime = travelTime+dt/2.
        
		call system_clock(time2)
		
		write(s1,'(i5)') n
		write(s2,'(i5)') nt
		write(s3,'(f15.3)') (time2-time1)/10000.
		
		write (string100,*) 'time step: '//trim(adjustl(s1))//'/'//trim(adjustl(s2))//', calculate time: '//trim(adjustl(s3))//'s'
		
		cal_time = cal_time + (time2-time1)/10000.
		
		! output wavefield and slice 
        if(mod(n,tstep) == 1) then
            call data_output_vector(Ey(1:nx:xstep,1:ny:xstep,1:nz:xstep),'wvf_Ey_'//isrc_s, Nx/xstep,Ny/xstep,Nz/xstep, n)
			call data_output_vector(Ey(islicex,:,:),'slx_Ey_'//isrc_s, nslicex,Ny,Nz, n)
			call data_output_vector(Ey(:,islicey,:),'sly_Ey_'//isrc_s, nx,nslicey,Nz, n)
			call data_output_vector(Ey(:,:,islicez),'slz_Ey_'//isrc_s, nx,Ny,nslicez, n)
            !field_name = 'travelTime'
            !call data_output_scalar(travelTime-dt/2., field_name, n)
			
		    call system_clock(time3)
		    write(s4,'(f15.3)') (time3-time2)/10000.
		    write (string100,*) trim(adjustl(string100))//', I/O time: '//trim(adjustl(s4))//'s'
		    io_time = io_time + (time3-time2)/10000.
		endif
        
		print *,trim(adjustl(string100))
		write(11,*) trim(adjustl(string100))
		
        !if(mod(n,tstep) == 1) then
        !    do j = 2,Ny-1
        !        do i = 2,Nx-1
        !            EzSurface(i,j) = Ez(i,j,RecordSurfaceEz(i,j))
        !        enddo
        !    enddo
        !    
        !    do j = 1,Ny-1
        !        do i = 2,Nx-1
        !            HxSurface(i,j) = Hx(i,j,RecordSurfaceHx(i,j))
        !        enddo
        !    enddo
        !    
        !    do j = 2,Ny-1
        !        do i = 1,Nx-1
        !            HySurface(i,j) = Hy(i,j,RecordSurfaceHy(i,j))
        !        enddo
        !    enddo
        !    write(997) Hxsurface
        !    write(998) Hysurface
        !    write(999) EzSurface
        !    write(995,*) travelTime
        !    !write(996,*) travelTime-dt/2.
        !endif
		
        
    end do
    
	
	! output gather
    call system_clock(time3)
	call data_output_vector(gather,'gather_'//isrc_s, nt, nrec, 1, 0)
    
    do i = 1,nthreads
        status = DftiFreeDescriptor(df_pointer_x(i)%desc_handle)
        status = DftiFreeDescriptor(df_pointer_y(i)%desc_handle)
        status = DftiFreeDescriptor(df_pointer_z(i)%desc_handle)
    enddo

    call system_clock(time4)
	write(s5,'(f15.3)') (time4-time0)/10000.
	
	write(s3,'(f15.3)') cal_time
	write(s4,'(f15.3)') io_time + (time4-time3)/10000.
	write (string100,*) 'Total time: '//trim(adjustl(s5))//'s, Total calculate time: '//trim(adjustl(s3))//'s, Total I/O time: '//trim(adjustl(s4))//'s'
	print *,trim(adjustl(string100))
	write(11,*) trim(adjustl(string100))
	close(11)
    
    
end program SGPSTD3D
    
subroutine FFTDiffl(A, B, N,coef,My_Desc_Handle)
!=============================
! FFTDIFFL
!=============================
    use MKL_DFTI
    implicit none
    integer, intent(in) :: N
    type(dfti_descriptor),pointer,intent(inout) :: My_Desc_Handle
    real(4), intent(in)    :: A(N), coef(N-1)
    real(4), intent(out)   :: B(N-1)
    real(4)                :: work(N+2)
    real(4)                :: temp1, temp2
	integer             :: N2,STATUS
    integer             :: i


	N2 = N/2
	WORK(1:N) = A(:)

	Status = DftiComputeForward(My_Desc_Handle, Work)
	WORK(1) = 0
	WORK(N+1) = coef(N-1)*WORK(N+1)
	do i = 2, N2
		temp1 = WORK(2*i-1)
		temp2 = WORK(2*i)
		WORK(2*i-1) = -WORK(2*i)*coef(i-1)+coef(N2+i-2)*Temp1
		WORK(2*i) = temp1*coef(i-1)+coef(N2+i-2)*Temp2
    end do
    
	Status = DftiComputeBackward(My_Desc_Handle, Work)
	B(1:N-1) = WORk(2:N)/N

End subroutine FFTDiffl

subroutine FFTDiffr(A, B, N, coef,My_Desc_Handle) 
!=============================
! FFTDIFFR
!=============================
    use MKL_DFTI
    implicit none
    integer, intent(in) :: N
    type(dfti_descriptor),pointer,intent(inout) :: My_Desc_Handle
    real(4), intent(in)    :: A(N), coef(N-1)
    real(4), intent(out)   :: B(N-1)
    real(4)                :: work(N+2)
    real(4)                :: temp1, temp2
	integer             :: N2,STATUS
    integer             :: i
    
	N2 = N/2
	work(1:N)   = A(1:N)

	Status = DftiComputeForward(My_Desc_Handle, Work)
	WORK(1) = 0
	WORK(N+1) = -coef(N-1)*WORK(N+1)
	do i = 2, N2
		temp1 = WORK(2*i-1)
		temp2 = WORK(2*i)
		WORK(2*i-1) = -WORK(2*i)*coef(i-1)-coef(N2+i-2)*Temp1
		WORK(2*i) = temp1*coef(i-1)-coef(N2+i-2)*Temp2
    end do
    
	Status = DftiComputeBackward(My_Desc_Handle, Work)
	B(1:N-1) = WORk(1:N-1)/N
    
end subroutine FFTDiffr
    
    
subroutine data_output_vector(field, field_name, Nx, Ny, Nz, t)
	integer, intent(in) :: Nx, Ny, Nz
	real(4), intent(in)    :: field(Nx, Ny, Nz)
	integer, intent(in) :: t
	character(len=11)   :: field_name
	!character(len=4)    :: s2 = '.bin'
	character*50        :: s
    character*10        :: a
    
    write(a,'(i5.5)') t
	
    s = 'output/'//trim(adjustl(field_name))//'_'//trim(adjustl(a))//'.bin'
	open(1, file = trim(s), access = 'stream', form = 'unformatted', status = 'replace')
    ! NOTICE:
    ! Matlab & Fortran: store element in the first dimension first 
    ! C/CPP & python: store element in the last dimension first 
    write(1) field
	close(1)

end subroutine data_output_vector

subroutine data_output_scalar(field,field_name, t)
    implicit none
	real(4), intent(in)    :: field
	integer, intent(in) :: t
	character*10        :: field_name
	character(len=50)   :: s
    character*10        :: a

    
    write(a,'(i10)') t
	
        
    s = 'output/'//trim(adjustl(field_name))//trim(adjustl(a))//'.txt'
    open(1,file = trim(s), status = 'replace')
    write(1,*) field
    close(1)

end subroutine data_output_scalar
    

        

