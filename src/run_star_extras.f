! ***********************************************************************
!
!   Copyright (C) 2011  Bill Paxton
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************
     

    module run_star_extras

    use star_lib
    use star_def
    use const_def

    implicit none

    integer :: time0, time1, clock_rate
    double precision, parameter :: expected_runtime = 28.2 ! minutes
    integer :: numTacc,j                          !track location of accretion times read

    contains

    subroutine extras_controls(id, ierr)
        integer, intent(in) :: id
        integer, intent(out) :: ierr
        integer :: i
        type (star_info), pointer :: s
        ierr = 0
        call star_ptr(id, s, ierr)
        if (ierr /= 0) return

        s% extras_startup => extras_startup
        s% extras_check_model => extras_check_model
        s% extras_finish_step => extras_finish_step
        s% extras_after_evolve => extras_after_evolve
        s% how_many_extra_history_columns => how_many_extra_history_columns
        s% data_for_extra_history_columns => data_for_extra_history_columns
        s% how_many_extra_profile_columns => how_many_extra_profile_columns
        s% data_for_extra_profile_columns => data_for_extra_profile_columns 

        s% other_energy => energy_routine
        s% other_adjust_mdot=> mass_loss
    end subroutine extras_controls


    subroutine mass_loss(id, ierr)
        use star_def
        type(star_info), pointer :: s
        integer, intent(in) :: id
        integer, intent(out) :: ierr
        	 
        !x controls (n, a, M, f)
        double precision :: frac_absorbed_euv, frac_absorbing_radius, host_star_mass, right_side
        double precision :: escape_rate_reduction_factor, orbital_distance, eddy_coeff, homopause_pressure

        double precision :: planet_radius_cgs, planet_age_cgs, planet_mass_cgs, r_constant, molar_mass, height
        double precision :: mass_fractionation_effect, homopause_temp, homopause_radius
        double precision :: h1_number_frac, he4_number_frac,atomic_mass_he4,atomic_mass_h1
        double precision :: he3_num_frac, c12_num_frac, n14_num_frac, o16_num_frac, ne20_num_frac, mg24_num_frac


        double precision :: h1_atoms, he3_atoms, he4_atoms, c12_atoms, mg_mass
        double precision :: n14_atoms,o16_atoms, ne20_atoms,mg24_atoms, total_atoms

        !Dependent Variables
        double precision :: luminosity_euv, epsilon, K_tides, LOG_LEUV, total_loss, binary_temp_constant

        !Diffusion and energy limited escape rate
        double precision :: escape_dl, escape_el, envelope_mass, initial_hydrogen_fraction, temp_star_mass
        double precision :: escape_rate_h1, escape_rate_he4, total_loss_rate, initial_helium_fraction

        !Final Values for helium and hydrogen escape
        double precision :: escape_h1, escape_he4, frac_change_h1, frac_change_he4, total_frac_change, frac_other_isotopes
        double precision :: scale_height, bi_diff_co, surface_pressure, radius_above_surface

        integer :: k,i,j

        call star_ptr(id, s, ierr)
        if (ierr /= 0) return
        ierr = 0

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!CONSTANTS, AND SIMPLE CONVERSIONS!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !Converting all the planetary paremeters to cgs
        planet_radius_cgs= (10 ** (s% log_surface_radius)) * Rsun
        planet_age_cgs= (s% star_age) * 365 * 24 * 3600 + (1.892 * (10 ** 14))
        planet_mass_cgs = (s% star_mass) * msun

        !Setting Parameters defined by the x_controls
        frac_absorbed_euv= (s% x_ctrl(50))
        frac_absorbing_radius = s% x_ctrl(51)
        host_star_mass= s% x_ctrl(52) * msun
        escape_rate_reduction_factor= s% x_ctrl(53)
        orbital_distance= s% x_ctrl(54)* au
        eddy_coeff = s% x_ctrl(55) ! Eddy coefficient, should be ~1e9
        homopause_pressure = s% x_ctrl(56) !Homopause Pressure in pascals

        !Heigh parameters
        homopause_temp = 10 ** (s% log_surface_temperature)
        r_constant = 8.314 * (10 ** 7) ! Newton cm / mol K

        !From Hu and Seager
        mass_fractionation_effect = 8.0 * (10 ** 20)

        !The number of atoms of each species
        h1_atoms = (s% star_mass_h1 * msun)/(amu)
        he3_atoms = (s% star_mass_he3 * msun)/(3 * amu)
        he4_atoms = (s% star_mass_he4 * msun)/(4 * amu)
        c12_atoms = (s% star_mass_c12 * msun)/(12 * amu)
        n14_atoms = (s% star_mass_n14 * msun)/(14 * amu)
        o16_atoms = (s% star_mass_o16 * msun)/(16 * amu)
        ne20_atoms = (s% star_mass_ne20 * msun)/(20 * amu)

        !Calculating the mg atoms
        mg_mass = 0
        mg_mass = s% xa(8,1) * s% star_mass * msun
        mg24_atoms = (mg_mass)/(24 * amu)

        !Calculating Abundances
        total_atoms = h1_atoms + he3_atoms + he4_atoms + c12_atoms &
        + n14_atoms + n14_atoms + o16_atoms + ne20_atoms + mg24_atoms

        !Calculating Mole Fraction
        atomic_mass_h1 = 1 * amu
        atomic_mass_he4 = 4 * amu
        h1_number_frac  = h1_atoms / total_atoms
        he4_number_frac = he4_atoms / total_atoms

        he3_num_frac = he3_atoms / total_atoms
        c12_num_frac = c12_atoms / total_atoms
        n14_num_frac = n14_atoms / total_atoms
        o16_num_frac = o16_atoms / total_atoms
        ne20_num_frac = ne20_atoms / total_atoms
        mg24_num_frac = mg24_atoms / total_atoms

        !Calculating the molar mass from abundances
        molar_mass = (h1_number_frac + (he3_num_frac * 3) + (he4_number_frac * 4) + (c12_num_frac * 12) + &
        (n14_num_frac * 14) + (o16_num_frac * 16) + (ne20_num_frac * 20) + (mg24_num_frac * 24))! grams/mol


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!BEGINING OF MORE COMPLES EQUATIONS!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        frac_absorbing_radius = 1
        scale_height = (r_constant * homopause_temp) / (molar_mass * (10 ** (s% log_surface_gravity)))

        !The 3.036 is from eddy dif equations. I need to put this in the paper
        homopause_pressure = (3.036d-10 * ((10 ** s% log_surface_temperature) ** 3.5) * 10) ! 10 for barye
        radius_above_surface  = -1 * scale_height * LOG(homopause_pressure / (10 ** (s% log_surface_pressure)))

        !Homopause Radius
        homopause_radius = planet_radius_cgs + radius_above_surface

        !Calculating the luminosity
        !The 29.12 instead of 22.12 is to convert to ergs
        LOG_LEUV = 29.12 - 1.24 * (LOG10((planet_age_cgs / (3.154 * (10 ** 16)))))
        luminosity_euv = 10 ** LOG_LEUV

        !K_tides doesn't matter too much. It's usually .99, whereas other factors vary by orders of magnitude
        epsilon = (((planet_mass_cgs / (host_star_mass)) / 3) ** (1 / 3)) * ((orbital_distance) / planet_radius_cgs)
        K_tides = (1 - (3 / (2 * epsilon)) + (1 / (2 * (epsilon ** 3))))

        envelope_mass = ((s% star_mass_h1 * msun)/ s% xa(1,1))
        temp_star_mass = s% star_mass * msun

        !The beginning of the interesting flux rates
        !This is in terms of the number of atoms
        escape_dl = (standard_cgrav * planet_mass_cgs * (atomic_mass_he4 - atomic_mass_h1) &
        * mass_fractionation_effect) / ((homopause_radius ** 2) * kerg * homopause_temp)

        !Energy limited escape rate in the subsonic regime
        !This is in terms of the mass per second
        escape_el = (luminosity_euv * frac_absorbed_euv * (frac_absorbing_radius ** 2) * &
        (homopause_radius ** 3)) / (4 * K_tides * (orbital_distance ** 2) * standard_cgrav * planet_mass_cgs)

        right_side = escape_dl * h1_number_frac * atomic_mass_h1 * 4 * pi * (homopause_radius ** 2)
        IF (escape_el < right_side) THEN
            escape_rate_h1 = escape_el
            escape_rate_he4 = 0

            total_loss_rate = (escape_rate_h1 + escape_rate_he4)
            s% mstar_dot = -total_loss_rate

            initial_hydrogen_fraction = s% xa(1,1)
            !do i = 1, s% nz
            !    s% xa(1,i) = ((envelope_mass * initial_hydrogen_fraction) - (escape_rate_h1 * s% dt)) &
            !    / (envelope_mass - (escape_rate_h1 * s% dt))

            !    s% xa(2,i) = (s% xa(2,i) * envelope_mass) / (envelope_mass - (total_loss_rate * s% dt))
            !    s% xa(3,i) = (s% xa(3,i) * envelope_mass) / (envelope_mass - (total_loss_rate * s% dt))
            !    s% xa(4,i) = (s% xa(4,i) * envelope_mass) / (envelope_mass - (total_loss_rate * s% dt))
            !    s% xa(5,i) = (s% xa(5,i) * envelope_mass) / (envelope_mass - (total_loss_rate * s% dt))
            !    s% xa(6,i) = (s% xa(6,i) * envelope_mass) / (envelope_mass - (total_loss_rate * s% dt))
            !    s% xa(7,i) = (s% xa(7,i) * envelope_mass) / (envelope_mass - (total_loss_rate * s% dt))
            !    s% xa(8,i) = (s% xa(8,i) * envelope_mass) / (envelope_mass - (total_loss_rate * s% dt))
            !end do
        END IF

        IF (escape_el > right_side) THEN
            escape_rate_h1 = ((escape_el * atomic_mass_h1 * h1_number_frac) + (escape_dl * atomic_mass_h1 * atomic_mass_he4 &
            * h1_number_frac * he4_number_frac * 4 * pi * (homopause_radius ** 2))) &
            / ((atomic_mass_h1 * h1_number_frac) + (atomic_mass_he4 * he4_number_frac))

            escape_rate_he4 = ((escape_el * atomic_mass_he4 * he4_number_frac) - (escape_dl * atomic_mass_h1 * atomic_mass_he4 &
            * h1_number_frac * he4_number_frac * 4 * pi * (homopause_radius ** 2))) &
            / ((atomic_mass_h1 * h1_number_frac) + (atomic_mass_he4 * he4_number_frac))

            total_loss_rate = (escape_rate_h1 + escape_rate_he4)
            s% mstar_dot = -total_loss_rate

            initial_hydrogen_fraction = s% xa(1,1)
            initial_helium_fraction = s% xa(3,1)

            !do i = 1, s% nz
            !    s% xa(1,i) = ((envelope_mass * initial_hydrogen_fraction) - (escape_rate_h1 * s% dt)) &
            !    / (envelope_mass - (total_loss_rate * s% dt))

            !    s% xa(3,i) = ((envelope_mass * initial_helium_fraction) - (escape_rate_he4 * s% dt)) &
            !    / (envelope_mass - (total_loss_rate * s% dt))

            !    s% xa(2,i) = (s% xa(2,i) * envelope_mass) / (envelope_mass - (total_loss_rate * s% dt))
            !    s% xa(4,i) = (s% xa(4,i) * envelope_mass) / (envelope_mass - (total_loss_rate * s% dt))
            !    s% xa(5,i) = (s% xa(5,i) * envelope_mass) / (envelope_mass - (total_loss_rate * s% dt))
            !    s% xa(6,i) = (s% xa(6,i) * envelope_mass) / (envelope_mass - (total_loss_rate * s% dt))
            !    s% xa(7,i) = (s% xa(7,i) * envelope_mass) / (envelope_mass - (total_loss_rate * s% dt))
            !    s% xa(8,i) = (s% xa(8,i) * envelope_mass) / (envelope_mass - (total_loss_rate * s% dt))
            !end do
        END IF

    end subroutine mass_loss


    subroutine energy_routine(id, ierr)
        type (star_info), pointer :: s
        !use const_def, only: Rsun
        integer, intent(in) :: id
        integer, intent(out) :: ierr

        !use star_lib, only: star_ptr
        double precision :: extraheat,junk,diag_heat
        integer :: k,i,n,counter,z,p,numOfAcc,zeta,jafter,jbefore,indexI
        double precision :: core_epsrad,Rpl,pressureDep,pressureHere,random_dp
        real(dp) :: tauHere,Vesc, KE,massTot,Cd,phi,Mpl,Rhopl,H,g,mH,Tacc
        real(dp), DIMENSION(700) :: arrayKE
        real(dp), DIMENSION(5000) :: arrayI

        call star_ptr(id, s, ierr)
        if (ierr /= 0) return
        ierr = 0

        ! INITALIZE
        s% extra_heat = 0.d0
        s% d_extra_heat_dlndm1(:) = 0d0
        s% d_extra_heat_dlnTm1(:) = 0d0

        k = s% nz
        extraheat = -s% x_ctrl(3) * s% M_center * s% T(s% nz) * s% dlnT_dt(s% nz) / s% dm(s% nz) ! erg/g/sec
        !assuming dlnT_dt is in s^-1

        ! EXTRA HEATING CONSTANT IN SPACE AND TIME
        ! Heat produced by radioactive decay due to 40K, 235U, 238U, 232Th, respectively

        k = s% nz
        core_epsrad = 36.7d-8 * exp(-5.543d-10 * s% star_age) ! erg/g/sec   40K
        core_epsrad = core_epsrad + 4.6d-8 * exp(-9.8485d-10 * s% star_age) ! erg/g/sec  235U
        core_epsrad = core_epsrad + 2.3d-8 * exp( -1.5513d-10 * s% star_age)! erg/g/sec  238U
        core_epsrad = core_epsrad + 1.3d-8 * exp( -0.4948d-10 * s% star_age)! erg/g/sec  232Th

        s% extra_heat(k) = (extraheat+ s% x_ctrl(4) * s% M_center * core_epsrad / s% dm(k)) ! erg/g/sec, core heat flux density
    end subroutine energy_routine


    integer function extras_startup(id, restart, ierr)
        integer, intent(in) :: id
        logical, intent(in) :: restart
        integer, intent(out) :: ierr
        type (star_info), pointer :: s
        ierr = 0
        call star_ptr(id, s, ierr)
        if (ierr /= 0) return
        extras_startup = 0
        call system_clock(time0,clock_rate)
        if (.not. restart) then
        call alloc_extra_info(s)
        else ! it is a restart
        call unpack_extra_info(s)
        end if
    end function extras_startup

          
    subroutine extras_after_evolve(id, id_extra, ierr)
        integer, intent(in) :: id, id_extra
        integer, intent(out) :: ierr
        type (star_info), pointer :: s
        real(dp) :: dt
        character (len=strlen) :: test
        
        ierr = 0
        call star_ptr(id, s, ierr)
        if (ierr /= 0) return
        call system_clock(time1,clock_rate)
        dt = dble(time1 - time0) / clock_rate / 60
        call GET_ENVIRONMENT_VARIABLE( &
           "MESA_TEST_SUITE_CHECK_RUNTIME", test, status=ierr, trim_name=.true.)
        if (ierr == 0 .and. trim(test) == 'true' .and. dt > 1.5*expected_runtime) then
           write(*,'(/,a70,2f12.1,99i10/)') &
              'failed: EXCESSIVE runtime, prev time, retries, backups, steps', &
              dt, expected_runtime, s% num_retries, s% num_backups, s% model_number
        else
           write(*,'(/,a50,2f12.1,99i10/)') 'runtime, prev time, retries, backups, steps', &
              dt, expected_runtime, s% num_retries, s% num_backups, s% model_number
        end if
        ierr = 0
    end subroutine extras_after_evolve

    

    !returns either keep_going, retry, backup, or terminate.
    integer function extras_check_model(id, id_extra)
        integer, intent(in) :: id, id_extra
        integer :: ierr
        type (star_info), pointer :: s
        ierr = 0
        call star_ptr(id, s, ierr)
        if (ierr /= 0) return
        extras_check_model = keep_going         
    end function extras_check_model


    integer function how_many_extra_history_columns(id, id_extra)
        integer, intent(in) :: id, id_extra
        integer :: ierr
        type (star_info), pointer :: s
        ierr = 0
        call star_ptr(id, s, ierr)
        if (ierr /= 0) return
        how_many_extra_history_columns = 19
    end function how_many_extra_history_columns
          
          
    subroutine data_for_extra_history_columns(id, id_extra, n, names, vals, ierr)
        integer, intent(in) :: id, id_extra, n
        character (len=maxlen_history_column_name) :: names(n)
        real(dp) :: vals(n)


        double precision :: luminosity_euv, frac_absorbed_euv,frac_absorbing_radius, K_tides, amu
        double precision :: orbital_distance, escape_el, escape_dl, atomic_mass_h1,atomic_mass_he4,planet_age_cgs

        integer, intent(out) :: ierr
        type (star_info), pointer :: s
        ierr = 0
        call star_ptr(id, s, ierr)
        if (ierr /= 0) return

        names(1)= "planet_radius_cgs"
        vals(1)= (10**(s% log_surface_radius))*Rsun

        names(2) = "planet_age_cgs"
        vals(2) = s% star_age*365*24*3600

        names(3) = "planet_mass_cgs"
        vals(3) = (s% star_mass)*msun

        names(4) = "Hydrogen_Mass"
        vals(4) = (s% star_mass_h1)*msun

        names(5) = "He4_Mass"
        vals(5) = (s% star_mass_he4)*msun

        names(6) = "Time Step dt"
        vals(6) = s% dt

        names(7) = "mdot"
        vals(7) = s% mstar_dot

        names(8) = "T(1)"
        vals(8) = s% T(1)

        names(9) = "xa(1,1)"
        vals(9) = s% xa(1,1)

        names(10) = "xa(2,1)"
        vals(10) = s% xa(2,1)

        names(11) = "xa(3,1)"
        vals(11) = s% xa(3,1)

        names(12) = "xa(4,1)"
        vals(12) = s% xa(4,1)

        names(13) = "xa(5,1)"
        vals(13) = s% xa(5,1)

        names(14) = "xa(6,1)"
        vals(14) = s% xa(6,1)

        names(15) = "xa(7,1)"
        vals(15) = s% xa(7,1)

        names(16) = "xa(8,1)"
        vals(16) = s% xa(8,1)

        names(17) = "xa(8,1)"
        vals(17) = s% xa(8,1)

        names(18) = "Temperature of base"
        vals(18) = s% T(s% nz)

        names(19) = "Pressure of Base"
        vals(19) = s% Pgas(s% nz)

    end subroutine data_for_extra_history_columns


    integer function how_many_extra_profile_columns(id, id_extra)
        use star_def, only: star_info
        integer, intent(in) :: id, id_extra
        integer :: ierr
        type (star_info), pointer :: s
        ierr = 0
        call star_ptr(id, s, ierr)
        if (ierr /= 0) return
        how_many_extra_profile_columns = 0
    end function how_many_extra_profile_columns
          
          
    subroutine data_for_extra_profile_columns(id, id_extra, n, nz, names, vals, ierr)
        use star_def, only: star_info, maxlen_profile_column_name
        use const_def, only: dp
        integer, intent(in) :: id, id_extra, n, nz
        character (len=maxlen_profile_column_name) :: names(n)
        real(dp) :: vals(nz,n)
        integer, intent(out) :: ierr
        type (star_info), pointer :: s
        integer :: k
        ierr = 0
        call star_ptr(id, s, ierr)
        if (ierr /= 0) return
    end subroutine data_for_extra_profile_columns
          







    ! returns either keep_going or terminate.
    integer function extras_finish_step(id, id_extra)
        use star_def
        integer, intent(in) :: id, id_extra
        integer :: ierr

        double precision :: frac_absorbed_euv, frac_absorbing_radius, host_star_mass, right_side
        double precision :: escape_rate_reduction_factor, orbital_distance, eddy_coeff, homopause_pressure

        double precision :: planet_radius_cgs, planet_age_cgs, planet_mass_cgs, r_constant, molar_mass, height
        double precision :: mass_fractionation_effect, homopause_temp, homopause_radius
        double precision :: h1_number_frac, he4_number_frac,atomic_mass_he4,atomic_mass_h1
        double precision :: he3_num_frac, c12_num_frac, n14_num_frac, o16_num_frac, ne20_num_frac, mg24_num_frac


        double precision :: h1_atoms, he3_atoms, he4_atoms, c12_atoms, mg_mass
        double precision :: n14_atoms,o16_atoms, ne20_atoms,mg24_atoms, total_atoms

        !Dependent Variables
        double precision :: luminosity_euv, epsilon, K_tides, LOG_LEUV, total_loss, binary_temp_constant

        !Diffusion and energy limited escape rate
        double precision :: escape_dl, escape_el, envelope_mass, initial_hydrogen_fraction, temp_star_mass
        double precision :: escape_rate_h1, escape_rate_he4, total_loss_rate, initial_helium_fraction

        !Final Values for helium and hydrogen escape
        double precision :: escape_h1, escape_he4, frac_change_h1, frac_change_he4, total_frac_change, frac_other_isotopes
        double precision :: scale_height, bi_diff_co, surface_pressure, radius_above_surface

        integer :: k,i,j

        type (star_info), pointer :: s
        ierr = 0
        call star_ptr(id, s, ierr)
        if (ierr /= 0) return
        extras_finish_step = keep_going
        call store_extra_info(s)


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!CONSTANTS, AND SIMPLE CONVERSIONS!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !Converting all the planetary paremeters to cgs
        planet_radius_cgs= (10 ** (s% log_surface_radius)) * Rsun
        planet_age_cgs= (s% star_age) * 365 * 24 * 3600 + (1.892 * (10 ** 14))
        planet_mass_cgs = (s% star_mass) * msun

        !Setting Parameters defined by the x_controls
        frac_absorbed_euv= (s% x_ctrl(50))
        frac_absorbing_radius = s% x_ctrl(51)
        host_star_mass= s% x_ctrl(52) * msun
        escape_rate_reduction_factor= s% x_ctrl(53)
        orbital_distance= s% x_ctrl(54)* au
        eddy_coeff = s% x_ctrl(55) ! Eddy coefficient, should be ~1e9
        homopause_pressure = s% x_ctrl(56) !Homopause Pressure in pascals


        !Heigh parameters
        homopause_temp = 10 ** (s% log_surface_temperature)
        r_constant = 8.314 * (10 ** 7) ! Newton cm / mol K


        !From Hu and Seager
        mass_fractionation_effect = 8.0 * (10 ** 20)

        !The number of atoms of each species
        h1_atoms = (s% star_mass_h1 * msun)/(amu)
        he3_atoms = (s% star_mass_he3 * msun)/(3 * amu)
        he4_atoms = (s% star_mass_he4 * msun)/(4 * amu)
        c12_atoms = (s% star_mass_c12 * msun)/(12 * amu)
        n14_atoms = (s% star_mass_n14 * msun)/(14 * amu)
        o16_atoms = (s% star_mass_o16 * msun)/(16 * amu)
        ne20_atoms = (s% star_mass_ne20 * msun)/(20 * amu)

        !Calculating the mg atoms
        mg_mass = 0
        mg_mass = s% xa(8,1) * s% star_mass * msun
        mg24_atoms = (mg_mass)/(24 * amu)

        !Calculating Abundances
        total_atoms = h1_atoms + he3_atoms + he4_atoms + c12_atoms &
        + n14_atoms + n14_atoms + o16_atoms + ne20_atoms + mg24_atoms

        !Calculating Mole Fraction
        atomic_mass_h1 = 1 * amu
        atomic_mass_he4 = 4 * amu
        h1_number_frac  = h1_atoms / total_atoms
        he4_number_frac = he4_atoms / total_atoms

        he3_num_frac = he3_atoms / total_atoms
        c12_num_frac = c12_atoms / total_atoms
        n14_num_frac = n14_atoms / total_atoms
        o16_num_frac = o16_atoms / total_atoms
        ne20_num_frac = ne20_atoms / total_atoms
        mg24_num_frac = mg24_atoms / total_atoms

        !Calculating the molar mass from abundances
        molar_mass = (h1_number_frac + (he3_num_frac * 3) + (he4_number_frac * 4) + (c12_num_frac * 12) + &
        (n14_num_frac * 14) + (o16_num_frac * 16) + (ne20_num_frac * 20) + (mg24_num_frac * 24))! grams/mol


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!BEGINING OF MORE COMPLES EQUATIONS!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        frac_absorbing_radius = 1
        scale_height = (r_constant * homopause_temp) / (molar_mass * (10 ** (s% log_surface_gravity)))

        !The 3.036 is from eddy dif equations. I need to put this in the paper
        homopause_pressure = (3.036d-10 * ((10 ** s% log_surface_temperature) ** 3.5) * 10) ! 10 for barye
        radius_above_surface  = -1 * scale_height * LOG(homopause_pressure / (10 ** (s% log_surface_pressure)))

        !Homopause Radius
        homopause_radius = planet_radius_cgs + radius_above_surface

        !Calculating the luminosity
        !The 29.12 instead of 22.12 is to convert to ergs
        LOG_LEUV = 29.12 - 1.24 * (LOG10((planet_age_cgs / (3.154 * (10 ** 16)))))
        luminosity_euv = 10 ** LOG_LEUV

        !K_tides doesn't matter too much. It's usually .99, whereas other factors vary by orders of magnitude
        epsilon = (((planet_mass_cgs / (host_star_mass)) / 3) ** (1 / 3)) * ((orbital_distance) / planet_radius_cgs)
        K_tides = (1 - (3 / (2 * epsilon)) + (1 / (2 * (epsilon ** 3))))

        envelope_mass = ((s% star_mass_h1 * msun)/ s% xa(1,1))
        temp_star_mass = s% star_mass * msun

        !The beginning of the interesting flux rates
        !This is in terms of the number of atoms
        escape_dl = (standard_cgrav * planet_mass_cgs * (atomic_mass_he4 - atomic_mass_h1) &
        * mass_fractionation_effect) / ((homopause_radius ** 2) * kerg * homopause_temp)

        !Energy limited escape rate in the subsonic regime
        !This is in terms of the mass per second
        escape_el = (luminosity_euv * frac_absorbed_euv * (frac_absorbing_radius ** 2) * &
        (homopause_radius ** 3)) / (4 * K_tides * (orbital_distance ** 2) * standard_cgrav * planet_mass_cgs)

        write(*,*) 'Running a change after finishing the step'
        right_side = escape_dl * h1_number_frac * atomic_mass_h1 * 4 * pi * (homopause_radius ** 2)
        IF (escape_el < right_side) THEN
            escape_rate_h1 = escape_el
            escape_rate_he4 = 0

            total_loss_rate = (escape_rate_h1 + escape_rate_he4)
            !s% mstar_dot = -total_loss_rate

            initial_hydrogen_fraction = s% xa(1,1)
            do i = 1, s% nz
                s% xa(1,i) = ((envelope_mass * initial_hydrogen_fraction) - (escape_rate_h1 * s% dt)) &
                / (envelope_mass - (escape_rate_h1 * s% dt))

                s% xa(2,i) = (s% xa(2,i) * envelope_mass) / (envelope_mass - (total_loss_rate * s% dt))
                s% xa(3,i) = (s% xa(3,i) * envelope_mass) / (envelope_mass - (total_loss_rate * s% dt))
                s% xa(4,i) = (s% xa(4,i) * envelope_mass) / (envelope_mass - (total_loss_rate * s% dt))
                s% xa(5,i) = (s% xa(5,i) * envelope_mass) / (envelope_mass - (total_loss_rate * s% dt))
                s% xa(6,i) = (s% xa(6,i) * envelope_mass) / (envelope_mass - (total_loss_rate * s% dt))
                s% xa(7,i) = (s% xa(7,i) * envelope_mass) / (envelope_mass - (total_loss_rate * s% dt))
                s% xa(8,i) = (s% xa(8,i) * envelope_mass) / (envelope_mass - (total_loss_rate * s% dt))
            end do
        END IF

        IF (escape_el > right_side) THEN
            escape_rate_h1 = ((escape_el * atomic_mass_h1 * h1_number_frac) + (escape_dl * atomic_mass_h1 * atomic_mass_he4 &
            * h1_number_frac * he4_number_frac * 4 * pi * (homopause_radius ** 2))) &
            / ((atomic_mass_h1 * h1_number_frac) + (atomic_mass_he4 * he4_number_frac))

            escape_rate_he4 = ((escape_el * atomic_mass_he4 * he4_number_frac) - (escape_dl * atomic_mass_h1 * atomic_mass_he4 &
            * h1_number_frac * he4_number_frac * 4 * pi * (homopause_radius ** 2))) &
            / ((atomic_mass_h1 * h1_number_frac) + (atomic_mass_he4 * he4_number_frac))

            total_loss_rate = (escape_rate_h1 + escape_rate_he4)
            !s% mstar_dot = -total_loss_rate

            initial_hydrogen_fraction = s% xa(1,1)
            initial_helium_fraction = s% xa(3,1)

            do i = 1, s% nz
                s% xa(1,i) = ((envelope_mass * initial_hydrogen_fraction) - (escape_rate_h1 * s% dt)) &
                / (envelope_mass - (total_loss_rate * s% dt))

                s% xa(3,i) = ((envelope_mass * initial_helium_fraction) - (escape_rate_he4 * s% dt)) &
                / (envelope_mass - (total_loss_rate * s% dt))

                s% xa(2,i) = (s% xa(2,i) * envelope_mass) / (envelope_mass - (total_loss_rate * s% dt))
                s% xa(4,i) = (s% xa(4,i) * envelope_mass) / (envelope_mass - (total_loss_rate * s% dt))
                s% xa(5,i) = (s% xa(5,i) * envelope_mass) / (envelope_mass - (total_loss_rate * s% dt))
                s% xa(6,i) = (s% xa(6,i) * envelope_mass) / (envelope_mass - (total_loss_rate * s% dt))
                s% xa(7,i) = (s% xa(7,i) * envelope_mass) / (envelope_mass - (total_loss_rate * s% dt))
                s% xa(8,i) = (s% xa(8,i) * envelope_mass) / (envelope_mass - (total_loss_rate * s% dt))
            end do
        END IF
    end function extras_finish_step






















          
          
    subroutine alloc_extra_info(s)
        integer, parameter :: extra_info_alloc = 1
        type (star_info), pointer :: s
        call move_extra_info(s,extra_info_alloc)
    end subroutine alloc_extra_info
          
          
    subroutine unpack_extra_info(s)
        integer, parameter :: extra_info_get = 2
        type (star_info), pointer :: s
        call move_extra_info(s,extra_info_get)
    end subroutine unpack_extra_info
          

    subroutine store_extra_info(s)
        integer, parameter :: extra_info_put = 3
        type (star_info), pointer :: s
        call move_extra_info(s,extra_info_put)
    end subroutine store_extra_info
      


    subroutine move_extra_info(s,op)
        integer, parameter :: extra_info_alloc = 1
        integer, parameter :: extra_info_get = 2
        integer, parameter :: extra_info_put = 3
        type (star_info), pointer :: s
        integer, intent(in) :: op
     
        integer :: i, j, num_ints, num_dbls, ierr
     
        i = 0
        ! call move_int or move_flg    
        num_ints = i
     
        i = 0
        ! call move_dbl       
     
        num_dbls = i
     
        if (op /= extra_info_alloc) return
        if (num_ints == 0 .and. num_dbls == 0) return
     
        ierr = 0
        call star_alloc_extras(s% id, num_ints, num_dbls, ierr)
        if (ierr /= 0) then
            write(*,*) 'failed in star_alloc_extras'
            write(*,*) 'alloc_extras num_ints', num_ints
            write(*,*) 'alloc_extras num_dbls', num_dbls
            stop 1
        end if
        contains

    subroutine move_dbl(dbl)
        real(dp) :: dbl
        i = i+1
        select case (op)
        case (extra_info_get)
           dbl = s% extra_work(i)
        case (extra_info_put)
           s% extra_work(i) = dbl
        end select
    end subroutine move_dbl
     
    subroutine move_int(int)
        integer :: int
        i = i+1
        select case (op)
        case (extra_info_get)
           int = s% extra_iwork(i)
        case (extra_info_put)
           s% extra_iwork(i) = int
        end select
    end subroutine move_int
     
    subroutine move_flg(flg)
        logical :: flg
        i = i+1
        select case (op)
        case (extra_info_get)
           flg = (s% extra_iwork(i) /= 0)
        case (extra_info_put)
           if (flg) then
              s% extra_iwork(i) = 1
           else
              s% extra_iwork(i) = 0
           end if
        end select
    end subroutine move_flg
  
    end subroutine move_extra_info

    subroutine enxa ( n, x, en )

        implicit none
        
        integer ( kind = 4 ) n

        real ( kind = 8 ) e1
        real ( kind = 8 ) ek
        real ( kind = 8 ) en(0:n)
        integer ( kind = 4 ) k
        real ( kind = 8 ) x
        
        en(0) = exp ( - x ) / x 
        call e1xb ( x, e1 )
        
        en(1) = e1
        do k = 2, n
           ek = ( exp ( - x ) - x * e1 ) / ( k - 1.0D+00 )
           en(k) = ek
           e1 = ek
        end do

        return
    end subroutine enxa

      subroutine e1xb ( x, e1 )
!*****************************************************************************80
!
!! E1XB computes the exponential integral E1(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    06 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) E1, the function value.
!
        implicit none

        real ( kind = 8 ) e1
        real ( kind = 8 ) ga
        integer ( kind = 4 ) k
        integer ( kind = 4 ) m
        real ( kind = 8 ) r
        real ( kind = 8 ) t
        real ( kind = 8 ) t0
        real ( kind = 8 ) x
        
        if ( x == 0.0D+00 ) then
           
           e1 = 1.0D+300

        else if ( x <= 1.0D+00 ) then

           e1 = 1.0D+00
           r = 1.0D+00

           do k = 1, 25
              r = -r * k * x / ( k + 1.0D+00 )**2
              e1 = e1 + r
              if ( abs ( r ) <= abs ( e1 ) * 1.0D-15 ) then
                 exit
              end if
           end do
    
           ga = 0.5772156649015328D+00
           e1 = - ga - log ( x ) + x * e1

        else

           m = 20 + int ( 80.0D+00 / x )
           t0 = 0.0D+00
           do k = m, 1, -1
              t0 = k / ( 1.0D+00 + k / ( x + t0 ) )
           end do
           t = 1.0D+00 / ( x + t0 )
           e1 = exp ( -x ) * t
    
        end if

        return
      end subroutine e1xb

      end module run_star_extras
      
