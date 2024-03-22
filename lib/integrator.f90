module integrator
	use FFTW2d
	implicit none

	integer, parameter :: neq = 3, highestderivative0 = 2
	integer :: nx, ny, nmargin
	real( kind = 8 ) :: dtime, Lx, Ly, dx, dy
	real(kind = 8), allocatable :: parameters(:) ! gamma, d, pa, pb, pe, ph, Tr, R, Amax
	logical :: periodicx, periodicy
	real(kind = 8) :: time
	real(kind = 8), allocatable :: real_data(:, :, :) ! u, v, c

contains
	subroutine initintegrator( dtime0, Lx0, Ly0, nx0, ny0, periodicx0, periodicy0, &
			parameters0 )
		real (kind = 8) :: dtime0, Lx0, Ly0, parameters0(:)
		integer :: nx0, ny0, highestderivative0
		logical :: periodicx0, periodicy0

		dtime = dtime0
		nx = nx0
		ny = ny0
		Lx = Lx0
		Ly = Ly0
		dx = Lx/real(nx)
		dy = Ly/real(ny)
		periodicx = periodicx0
		periodicy = periodicy0

		allocate(parameters(size(parameters0)))
		parameters = parameters0

		select case(highestderivative0)
		case (1)
			nmargin = 2
		case (2)
			nmargin = 2
		case (3)
			nmargin = 3
		case (4)
			nmargin = 3
		case (5)
			nmargin = 4
		case (6)
			nmargin = 4
		end select
		! adjust for certain parameter ranges
		nmargin = max(nmargin, max(ceiling(sqrt(parameters(1)) * parameters(8) / dx), &
			ceiling(sqrt(parameters(1)) * parameters(8) / dy)))

		allocate(real_data(neq, -nmargin:nx-1+nmargin, -nmargin:ny-1+nmargin))
		real_data = 0.0d0
	end subroutine
	function dynamicalequation(sometime, someinput, somemask)
		integer :: di, dj
		real (kind = 8) :: laplaceu(-nmargin:nx - 1 + nmargin, -nmargin:ny - 1 + nmargin), &
			laplacev(-nmargin:nx - 1 + nmargin, -nmargin:ny - 1 + nmargin)
		real (kind = 8) :: dynamicalequation(neq, -nmargin:nx - 1 + nmargin, -nmargin:ny - 1 + nmargin), &
			someinput(neq, -nmargin:nx - 1 + nmargin, -nmargin:ny - 1 + nmargin), &
			somemask(-nmargin:nx - 1 + nmargin, -nmargin:ny - 1 + nmargin)
		real (kind = 8) :: sometime

		laplaceu = getderivative(2, [1,0], dx, someinput(1, :, :)) + getderivative(2, [0,1], dy, someinput(1, :, :))
		laplacev = getderivative(2, [1,0], dx, someinput(2, :, :)) + getderivative(2, [0,1], dy, someinput(2, :, :))

		! equation u
		dynamicalequation(1, :, :) = somemask * (laplaceu + parameters(1) * f(someinput(1, :, :), someinput(2, :, :)))

		! equation v
		dynamicalequation(2, :, :) = parameters(2) * laplacev + parameters(1) * g(someinput(1, :, :), someinput(2, :, :))

		! equation c
		dynamicalequation(3, :, :) = dcdt(someinput(1, :, :), someinput(3, :, :))

		dynamicalequation(:, -nmargin:-1, :) = 0.0d0
		dynamicalequation(:, :, -nmargin:-1) = 0.0d0
		dynamicalequation(:, :, ny:ny-1+nmargin) = 0.0d0
		dynamicalequation(:, nx:nx-1+nmargin, :) = 0.0d0
	end function
	! helper functions
	function f(u, v)
		real (kind = 8) :: u(-nmargin:nx - 1 + nmargin, -nmargin:ny - 1 + nmargin), &
			v(-nmargin:nx - 1 + nmargin, -nmargin:ny - 1 + nmargin)
		real (kind = 8) :: f(-nmargin:nx - 1 + nmargin, -nmargin:ny - 1 + nmargin)

		f = parameters(3) * u + u * u - parameters(4) * u * v

	end function
	function g(u, v)
		real (kind = 8) :: u(-nmargin:nx - 1 + nmargin, -nmargin:ny - 1 + nmargin), &
			v(-nmargin:nx - 1 + nmargin, -nmargin:ny - 1 + nmargin)
		real (kind = 8) :: g(-nmargin:nx - 1 + nmargin, -nmargin:ny - 1 + nmargin)

		g = parameters(5) * u * u * u - v

	end function
	function dcdt(u, c)
		real (kind = 8) :: u(-nmargin:nx - 1 + nmargin, -nmargin:ny - 1 + nmargin), &
			c(-nmargin:nx - 1 + nmargin, -nmargin:ny - 1 + nmargin)
		real (kind = 8) :: dcdt(-nmargin:nx - 1 + nmargin, -nmargin:ny - 1 + nmargin)

		dcdt = parameters(1) * parameters(6) * c * (a(u) - c) * (c - 1.0d0)
	end function
	function a(u)
		real (kind = 8) :: u(-nmargin:nx - 1 + nmargin, -nmargin:ny - 1 + nmargin)
		real (kind = 8) :: a(-nmargin:nx - 1 + nmargin, -nmargin:ny - 1 + nmargin)
		integer :: somei, somej

		do somei = 0, nx - 1
			do somej = 0, ny - 1
				if (u(somei, somej) <= parameters(7)) then
					a(somei, somej) = 0.49d0
				else
					a(somei, somej) = 0.49d0 - 2.5d0 * (u(somei, somej) - parameters(7))
				end if
			end do
		end do

	end function
	! helper functions end
	function getderivative(nderiv, direction, delta, someinput)
		integer :: gi, gj, gk
		integer :: nderiv, direction(2), & ! x/[1,0] or y/[0,1]
				currentindex(2)
		real (kind = 8) :: delta
		real (kind = 8) :: someinput(-nmargin:nx - 1 + nmargin, -nmargin:ny - 1 + nmargin), &
			getderivative(-nmargin:nx - 1 + nmargin, -nmargin:ny - 1 + nmargin)
		real (kind = 8), parameter :: coeff1(-2:2) = [1.0d0/12.0d0, -2.0d0/3.0d0, 0.0d0, 2.0d0/3.0d0, -1.0d0/12.0d0], &
			coeff2(-2:2) = [-1.0d0/12.0d0, 4.0d0/3.0d0, -5.0d0/2.0d0, 4.0d0/3.0d0, -1.0d0/12.0d0]

		! https://en.wikipedia.org/wiki/Finite_difference_coefficient
		getderivative = 0.0d0
		select case(nderiv)
		case (1)
			do gi = 0, nx - 1
				do gj = 0, ny - 1
					do gk = -2, 2
						currentindex = [gi, gj] + gk * direction
						getderivative(gi, gj) = getderivative(gi, gj) + coeff1(gk) * someinput(currentindex(1), currentindex(2))
					end do
					getderivative(gi, gj) = getderivative(gi, gj) / delta
				end do
			end do
		case (2)
			do gi = 0, nx - 1
				do gj = 0, ny - 1
					do gk = -2, 2
						currentindex = [gi, gj] + gk * direction
						getderivative(gi, gj) = getderivative(gi, gj) + coeff2(gk) * someinput(currentindex(1), currentindex(2))
					end do
					getderivative(gi, gj) = getderivative(gi, gj) / delta**2
				end do
			end do
		case (3)
		case (4)
		case (5)
		case (6)
		end select

	end function
	subroutine step_RK4_2d( nstep )
		real (kind = 8) :: time0
		real (kind = 8), dimension(neq, -nmargin:nx - 1 + nmargin, -nmargin:ny - 1 + nmargin) :: k1, k2, k3, k4
		real (kind = 8), dimension(neq, -nmargin:nx - 1 + nmargin, -nmargin:ny - 1 + nmargin) :: tmp
		integer :: nstep, istep, jstep

		time0 = time
		do istep = 1, nstep
			tmp = real_data
			call handleperiodicity(tmp)
			k1 = dynamicalequation(time, tmp, getcellmask(tmp(3, :, :)))
			tmp(1, :, :) = max(0.0d0, tmp(1, :, :))
			tmp(2, :, :) = max(0.0d0, tmp(2, :, :))

			tmp = real_data + 0.5d0 * dtime * k1
			call handleperiodicity(tmp)
			k2 = dynamicalequation(time + 0.5d0 * dtime, tmp, getcellmask(tmp(3, :, :)))
			tmp(1, :, :) = max(0.0d0, tmp(1, :, :))
			tmp(2, :, :) = max(0.0d0, tmp(2, :, :))

			tmp = real_data + 0.5d0 * dtime * k2
			call handleperiodicity(tmp)
			k3 = dynamicalequation(time + 0.5d0 * dtime, tmp, getcellmask(tmp(3, :, :)))
			tmp(1, :, :) = max(0.0d0, tmp(1, :, :))
			tmp(2, :, :) = max(0.0d0, tmp(2, :, :))

			tmp = real_data + dtime * k3
			call handleperiodicity(tmp)
			k4 = dynamicalequation(time + dtime, tmp, getcellmask(tmp(3, :, :)))

			real_data = real_data + dtime / 6.0d0 * (k1 + 2.0d0 * k2 + 2.0d0 * k3 + k4)


			real_data(1, :, :) = max(0.0d0, min(parameters(9), &
				real_data(1, :, :) * getcellmask(real_data(3, :, :))))

			real_data(3, 0:nx - 1, 0:ny - 1) = antialias(real_data(3, 0:nx - 1, 0:ny - 1))

			real_data(3, :, :) = real_data(3, :, :) * getcellmask(real_data(3, :, :))

			real_data(2, 0:nx - 1, 0:ny - 1) = max(0.0d0, real_data(2, 0:nx - 1, 0:ny - 1))

			real_data(1, 0:nx - 1, 0:ny - 1) = antialias(real_data(1, 0:nx - 1, 0:ny - 1))
			real_data(2, 0:nx - 1, 0:ny - 1) = antialias(real_data(2, 0:nx - 1, 0:ny - 1))
			real_data(3, 0:nx - 1, 0:ny - 1) = antialias(real_data(3, 0:nx - 1, 0:ny - 1))

			time = time + dtime
		end do
		call handleperiodicity(real_data)
		time = time0 + nstep * dtime

	end subroutine
	function getcellmask(c)
		real (kind = 8), dimension(-nmargin:nx - 1 + nmargin, -nmargin:ny - 1 + nmargin) :: c, c_clipped, getcellmask
		integer :: somei, somej, somek, somel
		real (kind = 8) :: searchrange

		c_clipped = 0.0d0
		do somei = -nmargin, nx - 1 + nmargin
			do somej = -nmargin, ny - 1 + nmargin
				if (c(somei, somej) > 0.5d0) then
					c_clipped(somei, somej) = 1.0d0
				else
					c_clipped(somei, somej) = 0.0d0
				end if
			end do
		end do

		searchrange = sqrt(parameters(1)) * parameters(8)
		getcellmask = 0.0d0
		do somei = 0, nx - 1
			do somej = 0, ny - 1
				if (c_clipped(somei, somej) > 0.5d0) then
					outer : do somek = somei - ceiling(searchrange/dx), somei + ceiling(searchrange/dx)
						do somel = somej - ceiling(searchrange/dy), somej + ceiling(searchrange/dy)
							if (sqrt(real(somek - somei)**2*dx**2 + real(somel - somej)**2*dy**2) <= searchrange) then
								getcellmask(somek, somel) = 1.0d0
							end if
						end do
					end do outer
				end if
			end do
		end do
	end function
	function getborder(c)
		real (kind = 8), dimension(-nmargin:nx - 1 + nmargin, -nmargin:ny - 1 + nmargin) :: c, c_clipped, getborder
		integer :: somei, somej
		real (kind = 8) :: searchrange

		c_clipped = 0.0d0
		do somei = -nmargin, nx - 1 + nmargin
			do somej = -nmargin, ny - 1 + nmargin
				if (c(somei, somej) > 0.5d0) then
					c_clipped(somei, somej) = 1.0d0
				else
					c_clipped(somei, somej) = 0.0d0
				end if
			end do
		end do

		searchrange = sqrt(parameters(1)) * parameters(8)
		getborder = 0.0d0
		do somei = 0, nx - 1
			do somej = 0, ny - 1
				if (c_clipped(somei, somej) > 0.5d0) then
					if (c_clipped(somei + 1, somej + 0) < 0.5d0 .or. &
							c_clipped(somei + 0, somej + 1) < 0.5d0 .or. &
							c_clipped(somei + 0, somej - 1) < 0.5d0 .or. &
							c_clipped(somei - 1, somej + 0) < 0.5d0) then
						getborder(somei, somej) = 1.0d0
					end if
				end if
			end do
		end do
		c_clipped = getborder
		do somei = 0, nx - 1
			do somej = 0, ny - 1
				if (c_clipped(somei + 1, somej + 0) > 0.5d0 .or. &
						c_clipped(somei + 0, somej + 1) > 0.5d0 .or. &
						c_clipped(somei + 0, somej - 1) > 0.5d0 .or. &
						c_clipped(somei - 1, somej + 0) > 0.5d0) then
					getborder(somei, somej) = 1.0d0
				end if
			end do
		end do

	end function
	subroutine handleperiodicity(someinput)
		integer :: hi
		real (kind = 8), intent(inout) :: someinput(neq, -nmargin:nx - 1 + nmargin, -nmargin:ny - 1 + nmargin)

		if (periodicx .and. periodicy) then
			do hi = 1, neq
				call makeperiodic(someinput(hi, :, :))
			end do
		else if (periodicx) then
			do hi = 1, neq
				call makexperiodic(someinput(hi, :, :))
			end do
		else if (periodicy) then
			do hi = 1, neq
				call makeyperiodic(someinput(hi, :, :))
			end do
		end if
	end subroutine
	subroutine makeperiodic(someinput)
		integer :: mi, mj
		real (kind = 8), intent(inout) :: someinput(-nmargin:nx - 1 + nmargin, -nmargin:ny - 1 + nmargin)

		do mi = 0, ny - 1
			do mj = -nmargin, -1
				someinput(mj, mi) = someinput(mj + nx, mi)
			end do
			do mj = nx, nx -1 + nmargin
				someinput(mj, mi) = someinput(mj - nx, mi)
			end do
		end do
		do mi = 0, nx - 1
			do mj = -nmargin, -1
				someinput(mi, mj) = someinput(mi, mj + ny)
			end do
			do mj = ny, ny -1 + nmargin
				someinput(mi, mj) = someinput(mi, mj - ny)
			end do
		end do
	end subroutine
	subroutine makexperiodic(someinput)
		integer :: mi, mj
		real (kind = 8), intent(inout) :: someinput(-nmargin:nx - 1 + nmargin, -nmargin:ny - 1 + nmargin)

		do mi = 0, ny - 1
			do mj = -nmargin, -1
				someinput(mj, mi) = someinput(mj + nx, mi)
			end do
			do mj = nx, nx -1 + nmargin
				someinput(mj, mi) = someinput(mj - nx, mi)
			end do
		end do
	end subroutine
	subroutine makeyperiodic(someinput)
		integer :: mi, mj
		real (kind = 8), intent(inout) :: someinput(-nmargin:nx - 1 + nmargin, -nmargin:ny - 1 + nmargin)

		do mi = 0, nx - 1
			do mj = -nmargin, -1
				someinput(mi, mj) = someinput(mi, mj + ny)
			end do
			do mj = ny, ny -1 + nmargin
				someinput(mi, mj) = someinput(mi, mj - ny)
			end do
		end do
	end subroutine
	function antialias(someinput)
		integer :: mi, mj
		real (kind = 8) :: someinput(1:nx, 1:ny), antialias(1:nx, 1:ny)
		complex ( kind = 8 ) :: somecomplex(nx/2 + 1, ny)
		real (kind = 8), allocatable, save :: karrayx(:), karrayy(:)
		logical, save :: initialized = .false.

		if (.not. initialized) then
			allocate(karrayx(1:nx/2 + 1), karrayy(1:ny))
			do mi = 1, ny/2 + 1
				karrayy(mi) = mi - 1
			end do
			do mi = ny/2 + 2, ny
				karrayy(mi) = -ny + mi - 1
			end do
			do mi = 1, nx/2 + 1
				karrayx(mi) = mi - 1
			end do
			initialized = .true.
		end if

		somecomplex = dfft( nx, ny, someinput )

		do mi = 1, nx/2 + 1
			do mj = 1, ny
				if (sqrt((real(karrayx(mi)) / real(nx/2))**2 + (real(karrayy(mj)) / real(ny/2))**2) > &
						0.5d0) then
					somecomplex(mi, mj) = 0.0d0
				end if
			end do
		end do

		antialias = idfft( nx, ny, somecomplex )
	end function
	subroutine roughen(ramplitude, rseed)
		real(kind = 8) :: ramplitude
		integer :: rseed, irx, iry

		call srand(rseed)
		do irx = 0, nx - 1
			do iry = 0, ny - 1
				real_data(1, irx, iry) = real_data(1, irx, iry) &
							+ 2.0d0 * ramplitude * (rand() - 0.5d0)
				real_data(2, irx, iry) = real_data(2, irx, iry) &
							+ 2.0d0 * ramplitude * (rand() - 0.5d0)
				real_data(3, irx, iry) = real_data(3, irx, iry) &
							+ 2.0d0 * ramplitude * (rand() - 0.5d0)
			end do
		end do
	end subroutine
	subroutine snapshot( rinput, ksnapshots, cidentifier, filename0 )
		integer :: ksnapshots
		integer :: isnapshot, is
		character( len = 1000 ) :: filename, dummy
		character( len = 8 ) :: out_format = '(I5.5)'
		real ( kind = 8 ) :: rinput( :, : )
		character( * ) :: cidentifier
		character( * ), optional :: filename0

		if (present(filename0)) then
			open(unit = 110, file = filename0, status = 'replace')
		else
			! convert int to str
			write (dummy, out_format) ksnapshots
			! combine parts
			filename = 'data/out'//cidentifier//trim(dummy)//'.dat'
			open(unit = 110, file = trim(filename), status = 'replace')
		end if

		do is = 1, size(rinput(1,:))
			write(110, *) rinput(:, is)
		end do

		close(110)

	end subroutine
end module
