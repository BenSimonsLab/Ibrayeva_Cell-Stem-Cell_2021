! Stochastic simulation: developmental-like programme

! To compile with gfortran, run compile.sh.

! How to call the simulation, provide the parameters as command line arguments:
!   nsc_devlike [lambda] [lambda_1] [m_n] [omega_n] [q] [R0] [R1] [total simulation time]
!       [output time step] [no. of realisations] [seed for RNG] [output path (relative)]

program devlike
  implicit none

	integer,          parameter                           :: FATE_RR = 0, FATE_RN = 1, FATE_Q = 2

	character(255),   parameter                           :: data_directory = './results/'

  real(kind=8),     parameter                           :: pi2 = 6.283185307179586
	integer,          parameter                           :: max_clone_size = 50

	! parameters
	real(kind=8)                                          :: lambda, omega_n, m_n, lambda_1, q
	integer                                               :: R0, R1
	real(kind=8)                                          :: sim_time, dt_out
	integer                                               :: runs, seed

	! derived parameters
	integer                                               :: readout_points, r_tot

  ! fields
	integer,          dimension(:),           allocatable :: cells_r
	integer                                               :: last_cells_r
  integer                                               :: cells_n, last_cells_n

	! computation variables
	real(kind=8)                                          :: reaction, time_step
	real(kind=8),     dimension(:),           allocatable :: prop
	real(kind=8)                                          :: total_prop, p_sum, die, die2
	integer                                               :: gen
	logical                                               :: searching, first_div
	integer                                               :: lrp, rp

	real(kind=8)                                          :: nrm

  ! observables
  integer,          dimension(:,:,:,:),     allocatable :: clone_sizes
	integer,          dimension(:,:),         allocatable :: avg_clone_content, unrecorded_clones, avg_fate
	integer,          dimension(:),           allocatable :: surv_clones, r_clones, avg_rrcontent
	integer,          parameter                           :: CLONES_RN = 1, CLONES_RA = 2, CLONES_NA = 3
	integer                                               :: clone_type
	character(2)                         	                :: suffix

  ! counters
  integer                                               :: i, j, k, r
	real(kind=8)                                          :: t

	character(255)                                        :: directory, prefix

	! command line arguments
	integer,          parameter                           :: req_num_args = 12
  integer                                               :: num_args
  character(255),   dimension(1:req_num_args)           :: cl_args

  !write(*,'(a)',advance='no') achar(27)//'[36mnsc devlike simulation' // achar(27)// '[0m'

  num_args=iargc()
	if (num_args.ne.req_num_args) then
		write(*,'(a)') 'Wrong number of arguments. Terminated.'
		call exit()
	end if
  do i=1,num_args
    call getarg(i, cl_args(i))
  end do

	read(cl_args(1),*) lambda
	read(cl_args(2),*) lambda_1
	read(cl_args(3),*) m_n
  read(cl_args(4),*) omega_n
	read(cl_args(5),*) q
	read(cl_args(6),*) R0
	read(cl_args(7),*) R1
	read(cl_args(8),*) sim_time
	read(cl_args(9),*) dt_out
	read(cl_args(10),*) runs
	read(cl_args(11),*) seed
	prefix = trim(cl_args(12))

	readout_points = floor(sim_time / dt_out)
	r_tot = R0 + R1 + 1

	! allocate arrays
	allocate(cells_r(0:r_tot))	! index 0 corresponds to quiescent cells
	allocate(prop(0:r_tot))

	allocate(clone_sizes(1:3,0:readout_points,0:max_clone_size,0:max_clone_size))
	allocate(avg_clone_content(1:3,0:readout_points))
	allocate(surv_clones(0:readout_points), r_clones(0:readout_points))
	allocate(unrecorded_clones(1:3,0:readout_points))
	allocate(avg_fate(0:2,0:readout_points))
	allocate(avg_rrcontent(0:readout_points))

  ! create output directory
	directory = trim(data_directory) // trim(prefix) // '/'
  call system('mkdir -p ' // trim(directory))

	! initialize random number generator
  call init_random_seed(seed)

  ! initial conditions for observables
	r_clones = 0
	surv_clones = 0
	avg_clone_content = 0
	avg_fate = 0
	avg_rrcontent = 0
  clone_sizes = 0

	! realizations
	do r = 1, runs
		write (*,'(a,f5.1,a)',advance='no') '\b\b\b\b\b\b', (r/runs)*100.0, '%'

		! initial conditions
		cells_r = 0
		cells_n = 0

		call random_number(die)
		if (die < q) then
			call random_number(die2)
			cells_r(ceiling(die2 * r_tot)) = 1
		else
			cells_n = 1
		end if

		last_cells_r = sum(cells_r)
		last_cells_n = cells_n

		! initiate counters
		first_div = .true.
		rp = -1
		lrp = -1

		! time evolution
		t = 0.
		do while (t < sim_time)
			! ------------------------------------------------------------------------
			! partial propensities

			! probability of N death
			prop(0) = omega_n * real(cells_n,8)

			! probability of R division
			if (first_div) then
				do i = 1, r_tot
					prop(i) = lambda_1 * real(cells_r(i), 8)
				end do
				first_div = .false.
			else
				do i = 1, r_tot
					prop(i) = lambda * real(cells_r(i), 8)
				end do
			end if

			! ------------------------------------------------------------------------

			! total propensity
			total_prop = sum(prop)

			if (total_prop > 0.) then
				! advance time
				call random_number(time_step)
				t = t - log(time_step) / total_prop

				! compute reaction
				call random_number(reaction)
				reaction = total_prop * reaction

				! ----------------------------------------------------------------------
				! reactions

				p_sum = 0.
				searching = .true.

				! reaction 0: N death
				p_sum = p_sum + prop(0)
				if (reaction < p_sum .and. searching) then
					cells_n = cells_n - 1
					searching = .false.
				end if

				! R cell divisions
				do i = 1, r_tot
					p_sum = p_sum + prop(i)
					if (reaction < p_sum .and. searching) then
						if (i <= R0) then
							cells_r(i) = cells_r(i) - 1
							cells_r(i+1) = cells_r(i+1) + 2
						else if (i <= R0 + R1) then
							cells_r(i) = cells_r(i) - 1
							cells_r(i+1) = cells_r(i+1) + 1
							call random_number(die)
							cells_n = cells_n + nint(-log(1. - die) * m_n)
						else
							call random_number(die)
							cells_n = cells_n + 2 * nint(-log(1. - die) * m_n)
							cells_r(i) = cells_r(i) - 1
						end if

						searching = .false.
						exit
					end if
				end do

			else
				t = sim_time
			end if

			! ------------------------------------------------------------------------
			! update surviving clone size distribution

			rp = floor(t / dt_out)
			if (rp > lrp) then
				rp = min(rp, readout_points)

				! only record surviving clones
				if ((last_cells_r + last_cells_n) > 0) then

					! surviving clones and average clone sizes
					do k = lrp + 1, rp
						surv_clones(k) = surv_clones(k) + 1
						avg_clone_content(1,k) = avg_clone_content(1,k) + last_cells_r
						avg_clone_content(2,k) = avg_clone_content(2,k) + last_cells_n
					end do

					if (last_cells_r > 0) then
						do k = lrp + 1, rp
							r_clones(k) = r_clones(k) + 1
							avg_rrcontent(k) = avg_rrcontent(k) + last_cells_r
						end do
					end if

					! fraction of clones by 'fate'
					if (last_cells_r >= 2 .and. last_cells_n == 0) then
						do k = lrp + 1, rp
							avg_fate(FATE_RR,k) = avg_fate(FATE_RR,k) + 1
						end do
					end if
					if (last_cells_r >= 1 .and. last_cells_n > 0) then
						do k = lrp + 1, rp
							avg_fate(FATE_RN,k) = avg_fate(FATE_RN,k) + 1
						end do
					end if
					if (last_cells_r == 1 .and. last_cells_n == 0) then
						do k = lrp + 1, rp
							avg_fate(FATE_Q,k) = avg_fate(FATE_Q,k) + 1
						end do
					end if

					! joint clone size distributions
					if (last_cells_r <= max_clone_size .and. last_cells_n <= max_clone_size) then
						do k = lrp + 1, rp
							clone_sizes(CLONES_RN,k,last_cells_r,last_cells_n) = clone_sizes(CLONES_RN,k,last_cells_r,last_cells_n) + 1
						end do
					else
						do k = lrp + 1, rp
							unrecorded_clones(CLONES_RN,k) = unrecorded_clones(CLONES_RN,k) + 1
						end do
					end if

				end if

				lrp = rp
			end if

			last_cells_r = sum(cells_r)
			last_cells_n = cells_n

			if (total_prop == 0.) exit
		end do
	end do

	! --------------------------------------------------------------------------------------
	! write fraction of surviving clones
	open(unit=1,file=trim(directory) // 'surviving_fraction.csv')
		do k=0,readout_points
			write(1,'(f17.6,a1,f17.6,a1,f17.6)') real(k,8) * dt_out, ',', real(surv_clones(k),8) / real(runs,8), ',', real(r_clones(k),8) / real(surv_clones(k),8)
		end do
	close(1)

	! write average clone content
	open(unit=1,file=trim(directory) // 'avg_clone_content.csv')
		do k=0,readout_points
			write(1,'(f17.6,a1,f17.6,a1,f17.6,a1,f17.6,a1,f17.6)') real(k,8) * dt_out, ',', real(avg_clone_content(1,k),8) / real(surv_clones(k),8), ',', real(avg_clone_content(2,k),8) / real(surv_clones(k),8), ',', real(avg_clone_content(1,k),8) / real(runs,8), ',', real(avg_clone_content(2,k),8) / real(runs,8)
		end do
	close(1)

	! write 'fate' composition
	open(unit=1,file=trim(directory) // 'avg_fate.csv')
		do k=0,readout_points
			nrm = real(sum(avg_fate(:,k)), 8)
			if (nrm == 0.) then
				write(1,'(f17.6,a1,f17.6,a1,f17.6,a1,f17.6,a1)') real(k,8) * dt_out, ',', 0., ',', 0., ',', 0.
			else
				write(1,'(f17.6,a1,f17.6,a1,f17.6,a1,f17.6,a1)') real(k,8) * dt_out, ',', real(avg_fate(FATE_Q,k),8) / nrm, ',', real(avg_fate(FATE_RR,k),8) / nrm, ',', real(avg_fate(FATE_RN,k),8) / nrm
			end if
		end do
	close(1)

	! write R content of R clones
	open(unit=1,file=trim(directory) // 'avg_rrcontent.csv')
		do k=0,readout_points
			if (r_clones(k) == 0) then
				write(1,'(f17.6,a1,f17.6)') real(k,8) * dt_out, ',', 0.
			else
				write(1,'(f17.6,a1,f17.6)') real(k,8) * dt_out, ',', real(avg_rrcontent(k),8) / real(r_clones(k),8)
			end if
		end do
	close(1)

	! write clone size distributions
	do clone_type=1,1
		select case (clone_type)
			case (CLONES_RN)
				suffix = 'rn'
			case (CLONES_RA)
				suffix = 'ra'
			case (CLONES_NA)
				suffix = 'na'
		end select

		open(unit=1,file=trim(directory) // 'clone_sizes_' // suffix // '.csv')
			do k=0,readout_points
				do i=0,max_clone_size
					do j=0,max_clone_size-1
						write(1,'(i10,a1)',advance='no') clone_sizes(clone_type,k,i,j), ','
					end do
					write(1,'(i10)') clone_sizes(clone_type,k,i,max_clone_size)
				end do
			end do
		close(1)
	end do

	open(unit=1,file=trim(directory) // 'unrecorded_fraction.csv')
		do k=0,readout_points
			write(1,'(f17.6,a1,f17.6,a1,f17.6,a1,f17.6)') real(k,8) * dt_out, ',', real(unrecorded_clones(CLONES_RN,k),8) / real(runs,8), ',', real(unrecorded_clones(CLONES_RA,k),8) / real(runs,8), ',', real(unrecorded_clones(CLONES_NA,k),8) / real(runs,8)
		end do
	close(1)

  open(unit=1,file=trim(directory) // 'parameters.csv')
  	call write_parameters(1, lambda, lambda_1, R0, R1, omega_n, m_n, q, sim_time, dt_out, runs, seed)
	close(1)
end program devlike


subroutine write_parameters(output_unit, lambda, lambda_1, R0, R1, omega_n, m_n, q, sim_time, dt_out, runs, seed)
	implicit none

	integer,															   intent(in)   :: output_unit
	real(kind=8),                            intent(in)   :: lambda, lambda_1, omega_n, m_n, sim_time, dt_out, q
	integer,															   intent(in)   :: R0, R1, seed, runs

  write(output_unit,'(a20,f17.6)') 'lambda,',	    lambda
	write(output_unit,'(a20,f17.6)') 'lambda_1,',   lambda_1
	write(output_unit,'(a20,i7)')    'R0,',     		R0
	write(output_unit,'(a20,i7)')    'R1,',     		R1
	write(output_unit,'(a20,f17.6)') 'm_n,',				m_n
	write(output_unit,'(a20,f17.6)') 'omega_n,',	  omega_n
	write(output_unit,'(a20,f17.6)') 'q,',	    		q

	write(output_unit,'(a20,f17.6)') 'sim_time,',   sim_time
	write(output_unit,'(a20,f17.6)') 'dt_out,',     dt_out
  write(output_unit,'(a20,i10)')   'num_runs,',		runs
	write(output_unit,'(a20,i7)')    'seed,',    		seed
end subroutine write_parameters


subroutine init_random_seed(m)
  implicit none
  integer, intent(in) :: m
  integer :: i, n
  integer, dimension(:), allocatable :: seed

  call random_seed(size = n)
  allocate(seed(n))

  seed = m + 37 * (/ (i - 1, i = 1, n) /)
  call random_seed(PUT = seed)

  deallocate(seed)
end subroutine
