! compile: gfortran lib/fftw2d.f90 lib/integrator.f90 examples/run.f90 -lfftw3 -o run
program run
	use integrator
	implicit none

	integer :: i, j, k, nnx, nny
	real ( kind = 8 ) :: systemsize(2)
	real (kind = 8) :: fps
	real (kind = 8) :: tempx, tempy
	real (kind = 8), allocatable :: tmp_field(:,:,:)
	real (kind = 8), allocatable :: boundary00(:,:,:), boundary0_val0(:,:,:), &
		absreal_data(:,:)
	CHARACTER(LEN=60) :: buffer1
	CHARACTER(LEN=60) :: buffer2

	nnx = 800
	nny = 800
	systemsize(1) = 16.0d0*real(nnx)/real(nny)
	systemsize(2) = 16.0d0

	call initintegrator( &
		1.0d-5, systemsize(1), systemsize(2), nnx, nny, .true., .true., &
		[ &
			625.0d0, & ! gamma
			30.0d0, & ! d
			0.9d0, & ! pa
			0.8d0, & ! pb
			6.5d0, & ! pe
			1.0d0, & ! ph
			1.0d0, & ! Tr
			0.006d0, & ! R
			30.0d0 & ! Amax
		] &
	)

	call roughen(0.01d0, 123)
	do i = -nmargin, nx - 1 + nmargin
		do j = -nmargin, ny - 1 + nmargin
			tempx = real(i)/real(nx) * systemsize(1)
			tempy = real(j)/real(ny) * systemsize(2)

			if ((norm2([tempx, tempy] - 0.5d0 * systemsize))**2 < &
					0.4d0**2*1.0d0*(1.0d0 + 0.25d0 * cos(3*atan2(tempy, tempx)))) then
				real_data(3,i, j) = 1.0d0
			else
				real_data(3,i, j) = 0.001d0
			end if
		end do
	end do

	real_data(1, 0:nx - 1, 0:ny - 1) = antialias(real_data(1, 0:nx - 1, 0:ny - 1))
	real_data(2, 0:nx - 1, 0:ny - 1) = antialias(real_data(2, 0:nx - 1, 0:ny - 1))
	! real_data(3, 0:nx - 1, 0:ny - 1) = antialias(real_data(3, 0:nx - 1, 0:ny - 1))

	real_data(1,:,:) = getcellmask(real_data(3,:,:)) * (real_data(1,:,:) + 0.5d0)
	real_data(2,:,:) = getcellmask(real_data(3,:,:)) * (real_data(2,:,:) + 0.1d0)
	! real_data(3,:,:) = getcellmask(abs(real_data(3,:,:))) * abs(real_data(3,:,:))

	call snapshot( real_data(3,:,:), 0, '', 'ic.dat' )
	call snapshot( getcellmask(real_data(3,:,:)), 0, '', 'ic_mask.dat' )

	allocate(tmp_field(-nmargin:nx - 1 + nmargin, -nmargin:ny - 1 + nmargin, 1))

	do i = 0, 480*5
		call step_RK4_2d( max(1, floor(100.0d0*1.25d-5/dtime + 0.5d0)))
		write(*,*) i, time

		call snapshot( real_data(1,0:nx-1, 0:ny-1), i, 'u' )
		call snapshot( real_data(2,0:nx-1, 0:ny-1), i, 'v' )
		call snapshot( real_data(3,0:nx-1, 0:ny-1), i, 'c' )
		tmp_field(:,:,1) = getborder(real_data(3,:,:))
		call snapshot( tmp_field(0:nx-1, 0:ny-1, 1), i, 'c_border' )
		tmp_field(:,:,1) = getcellmask(real_data(3,:,:))
		call snapshot( tmp_field(0:nx-1, 0:ny-1, 1), i, 'omega' )
		tmp_field(:,:,1) = getborder(getcellmask(real_data(3,:,:)))
		call snapshot( tmp_field(0:nx-1, 0:ny-1, 1), i, 'omega_border' )

		if (mod(i, 100) == 0) then
			call snapshot( real_data(1,0:nx-1, 0:ny-1), i, 'backupu' )
			call snapshot( real_data(2,0:nx-1, 0:ny-1), i, 'backupv' )
			call snapshot( real_data(3,0:nx-1, 0:ny-1), i, 'backupc' )
		end if
	end do
contains
	function makefilename( ichars, inumber, prefix, suffix, out_format )
		! use as trim(makefilename( 99, 4, 'asd', '.daat' ))
		integer :: ichars, inumber
		character( ichars ) :: makefilename, dummy
		character(*) :: prefix, suffix
		character(*), optional :: out_format
		character( len = 8 ) :: out_format0 = '(I5.5)'

		if (.not. present(out_format)) then
			out_format0 = '(I5.5)'
		else
			out_format0 = out_format
		end if

		! convert int to str
		write (dummy, out_format0) inumber
		makefilename = prefix//trim(adjustl(dummy))//suffix
	end function
end program
