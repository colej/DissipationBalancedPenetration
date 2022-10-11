! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
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
      use math_lib

      implicit none

      real (dp) :: m_core, dm_core, dr_core, dr_core_div_h, r_core, rho_core_top

      ! these routines are called by the standard run_star check_model
      contains

      ! include 'standard_run_star_extras.inc'

      subroutine extras_controls(id, ierr)
          integer, intent(in) :: id
          integer, intent(out) :: ierr
          type (star_info), pointer :: s

          ierr = 0
          call star_ptr(id, s, ierr)
          if (ierr /= 0) return

          s% extras_startup => extras_startup
          s% extras_start_step => extras_start_step
          s% extras_check_model => extras_check_model
          s% extras_finish_step => extras_finish_step
          s% extras_after_evolve => extras_after_evolve
          s% how_many_extra_history_columns => how_many_extra_history_columns
          s% data_for_extra_history_columns => data_for_extra_history_columns
          s% how_many_extra_profile_columns => how_many_extra_profile_columns
          s% data_for_extra_profile_columns => data_for_extra_profile_columns

          s% how_many_extra_history_header_items => how_many_extra_history_header_items
          s% data_for_extra_history_header_items => data_for_extra_history_header_items
          s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
          s% data_for_extra_profile_header_items => data_for_extra_profile_header_items

          s% other_adjust_mlt_gradT_fraction => other_adjust_mlt_gradT_fraction_Peclet
          s% other_overshooting_scheme => extended_convective_penetration

          ! Add extra meshing
          ! s% use_other_mesh_delta_coeff_factor = .true.
          ! call read_inlist_xtra_coeff_core_boundary(ierr) ! Read inlist
          ! if (ierr /= 0) return
          ! s% other_mesh_delta_coeff_factor => mesh_delta_coeff_core_boundary


          ! s% kap_rq% Zbase = Z_ini ! set Z for opacity tables


      end subroutine extras_controls


      subroutine extras_startup(id, restart, ierr)
          integer, intent(in) :: id
          logical, intent(in) :: restart
          integer, intent(out) :: ierr
          type (star_info), pointer :: s
          ierr = 0
          call star_ptr(id, s, ierr)
          if (ierr /= 0) return
      end subroutine extras_startup


      integer function extras_start_step(id)
          integer, intent(in) :: id
          integer :: ierr
          type (star_info), pointer :: s
          ierr = 0
          call star_ptr(id, s, ierr)
          if (ierr /= 0) return
          extras_start_step = 0
      end function extras_start_step


      ! returns either keep_going, retry, backup, or terminate.
      integer function extras_check_model(id)
          integer, intent(in) :: id
          integer :: ierr
          type (star_info), pointer :: s
          logical :: do_retry
          ierr = 0
          call star_ptr(id, s, ierr)
          if (ierr /= 0) return
          extras_check_model = keep_going


          do_retry = .false.
          ! Terminate and save the pre-main sequence model when the convective core appears.
          ! The central hydrogen fraction needs to have decreased by a small amount,
          ! to make sure that core H-burning has started, and the star is near the ZAMS.
          ! If the model is not on the pre-main sequence, check if it needs to be saved,
          ! or if a retry needs to be made to save within precision at a desired Xc value.
          ! if (s%job% create_pre_main_sequence_model) then
          !     if ((s%mixing_type(s%nz) .eq. convective_mixing) .and. (X_ini-s% center_h1 > 1d-6)) then
          !         extras_check_model = terminate
          !     endif
          ! else
          !     call save_at_Xc(id, Xc_save, Xc_precision, Xc_save_step, need_to_save, do_retry, ierr)
          !     if (do_retry) extras_check_model = retry
          ! endif

          ! by default, indicate where (in the code) MESA terminated
          if (extras_check_model == terminate) s% termination_code = t_extras_check_model
      end function extras_check_model


      integer function how_many_extra_history_columns(id)
          integer, intent(in) :: id
          integer :: ierr
          type (star_info), pointer :: s
          ierr = 0
          call star_ptr(id, s, ierr)
          if (ierr /= 0) return
          how_many_extra_history_columns = 6
      end function how_many_extra_history_columns


      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
          integer, intent(in) :: id, n
          character (len=maxlen_history_column_name) :: names(n)
          real(dp) :: vals(n)
          integer, intent(out) :: ierr
          type (star_info), pointer :: s
          ierr = 0
          call star_ptr(id, s, ierr)
          if (ierr /= 0) return


          names(1) = 'm_core_pen'
          names(2) = 'dm_core_pen'
          names(3) = 'dr_core_pen'
          names(4) = 'dr_core_div_h_pen'
          names(5) = 'r_core_pen'
          names(6) = 'rho_core_top_pen'

          vals(1) = m_core
          vals(2) = dm_core
          vals(3) = dr_core
          vals(4) = dr_core_div_h
          vals(5) = r_core
          vals(6) = rho_core_top


      end subroutine data_for_extra_history_columns


      integer function how_many_extra_profile_columns(id)
          integer, intent(in) :: id
          integer :: ierr
          type (star_info), pointer :: s
          ierr = 0
          call star_ptr(id, s, ierr)
          if (ierr /= 0) return
          how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns


      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
          integer, intent(in) :: id, n, nz
          character (len=maxlen_profile_column_name) :: names(n)
          real(dp) :: vals(nz,n)
          integer, intent(out) :: ierr
          ! integer :: vals_nr
          type (star_info), pointer :: s
          integer :: k
          ierr = 0
          call star_ptr(id, s, ierr)
          if (ierr /= 0) return

          ! vals_nr = 1
          ! call prof_hydroEq(id, names, vals, vals_nr, nz, ierr)

      end subroutine data_for_extra_profile_columns


      integer function how_many_extra_history_header_items(id)
          integer, intent(in) :: id
          integer :: ierr
          type (star_info), pointer :: s
          ierr = 0
          call star_ptr(id, s, ierr)
          if (ierr /= 0) return
          how_many_extra_history_header_items = 0
      end function how_many_extra_history_header_items


      subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
          integer, intent(in) :: id, n
          character (len=maxlen_history_column_name) :: names(n)
          real(dp) :: vals(n)
          type(star_info), pointer :: s
          integer, intent(out) :: ierr
          ierr = 0
          call star_ptr(id,s,ierr)
          if(ierr/=0) return

          ! here is an example for adding an extra history header item
          ! also set how_many_extra_history_header_items
          ! names(1) = 'mixing_length_alpha'
          ! vals(1) = s% mixing_length_alpha
      end subroutine data_for_extra_history_header_items


      integer function how_many_extra_profile_header_items(id)
          integer, intent(in) :: id
          integer :: ierr
          type (star_info), pointer :: s
          ierr = 0
          call star_ptr(id, s, ierr)
          if (ierr /= 0) return
          how_many_extra_profile_header_items = 0
      end function how_many_extra_profile_header_items


      subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
          integer, intent(in) :: id, n
          character (len=maxlen_profile_column_name) :: names(n)
          real(dp) :: vals(n)
          type(star_info), pointer :: s
          integer, intent(out) :: ierr
          ierr = 0
          call star_ptr(id,s,ierr)
          if(ierr/=0) return

          ! here is an example for adding an extra profile header item
          ! also set how_many_extra_profile_header_items
          ! names(1) = 'mixing_length_alpha'
          ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_profile_header_items


      ! returns either keep_going or terminate.
      ! note: cannot request retry or backup; extras_check_model can do that.
      integer function extras_finish_step(id)
          integer, intent(in) :: id
          integer :: ierr
          type (star_info), pointer :: s
          ierr = 0
          call star_ptr(id, s, ierr)
          if (ierr /= 0) return
          extras_finish_step = keep_going

          if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step

      end function extras_finish_step


      subroutine extras_after_evolve(id, ierr)
          integer, intent(in) :: id
          integer, intent(out) :: ierr
          type (star_info), pointer :: s
          ierr = 0
          call star_ptr(id, s, ierr)
          if (ierr /= 0) return
      end subroutine extras_after_evolve

!!! CUSTOM

      subroutine other_adjust_mlt_gradT_fraction_Peclet(id, ierr)
          integer, intent(in) :: id
          integer, intent(out) :: ierr
          type(star_info), pointer :: s
          real(dp) :: fraction, Peclet_number, diffusivity, Hp       ! f is fraction to compose grad_T = f*grad_ad + (1-f)*grad_rad
          integer :: k
          logical, parameter :: DEBUG = .FALSE.

          ierr = 0
          call star_ptr(id, s, ierr)
          if (ierr /= 0) return

          if (s%D_mix(1) .ne. s%D_mix(1)) return  ! To ignore iterations where Dmix and gradT are NaNs

          if (s%num_conv_boundaries .lt. 1) then  ! Is zero at initialisation of the run
          if (DEBUG) then
              write(*,*) 'runstarex_gradT: skip since there are no convective boundaries'
          end if
          return
          endif

          do k= s%nz, 1, -1
              if (s%D_mix(k) <= s% min_D_mix) exit

              diffusivity = 16.0_dp * boltz_sigma * pow3(s% T(k)) / ( 3.0_dp * s% opacity(k) * pow2(s% rho(k)) * s% cp(k) )
              ! Hp = s% P(k)/(s% rho(k)*s% grav(k)) ! Pressure scale height
              Hp = s% scale_height(k) ! Pressure scale height
              Peclet_number = s% conv_vel(k) * Hp * s% mixing_length_alpha / diffusivity

              if (Peclet_number >= 100.0_dp) then
                  fraction = 1.0_dp
              else if (Peclet_number .le. 0.01_dp) then
                  fraction = 0.0_dp
              else
                  fraction = (safe_log10(Peclet_number)+2.0_dp)/4.0_dp
              end if

              s% adjust_mlt_gradT_fraction(k) = fraction
          end do

      end subroutine other_adjust_mlt_gradT_fraction_Peclet


      subroutine extended_convective_penetration(id, i, j, k_a, k_b, D, vc, ierr)
          integer, intent(in) :: id, i, j
          integer, intent(out) :: k_a, k_b
          real(dp), intent(out), dimension(:) :: D, vc
          integer, intent(out) :: ierr
          type (star_info), pointer :: s

          logical, parameter :: DEBUG = .FALSE.
          real(dp) :: f, f2, f0
          real(dp) :: D0, Delta0
          real(dp) :: w
          real(dp) :: factor
          real(dp) :: r_cb, Hp_cb
          real(dp) :: r_ob, D_ob, vc_ob
          logical  :: outward
          integer  :: dk, k, k_ob
          real(dp) :: r, dr, r_step

          ! real(dp) :: dm_core, dr_core, dr_core_div_h, r_core, m_core, rho_core_top


          ! Evaluate the overshoot diffusion coefficient D(k_a:k_b) and
          ! mixing velocity vc(k_a:k_b) at the i'th convective boundary,
          ! using the j'th set of overshoot parameters. The overshoot
          ! follows the extended convective penetration scheme description by Mathias
          ! Michielsen, "Probing the shape of the mixing profile and of the thermal
          ! structure at the convective core boundary through asteroseismology",
          ! A&A, 628, 76 (2019)

          ierr = 0
          call star_ptr(id, s, ierr)
          if (ierr /= 0) return

          ! call JermynAnders_active_penetration(s, id) !, m_core, dm_core, dr_core, dr_core_div_h, r_core, rho_core_top)
          call JermynAnders_penetration(s, id) !, m_core, dm_core, dr_core, dr_core_div_h, r_core, rho_core_top)

          ! Extract parameters
          f = dr_core_div_h        ! extend of step function (a_ov)
          f0 = s%overshoot_f0(j)
          f2 = 0.01            ! exponential decay (f_ov)

          D0 = s%overshoot_D0(j)
          Delta0 = s%overshoot_Delta0(j)

          if (f < 0.0_dp .OR. f0 <= 0.0_dp .OR. f2 < 0.0_dp) then
              write(*,*) 'ERROR: for extended convective penetration, must set f0 > 0, and f and f2 >= 0'
              write(*,*) 'see description of overshooting in star/defaults/control.defaults'
              ierr = -1
              return
          end if

          ! Evaluate convective boundary (_cb) parameters
          call star_eval_conv_bdy_r(s, i, r_cb, ierr)
          if (ierr /= 0) return

          call star_eval_conv_bdy_Hp(s, i, Hp_cb, ierr)
          if (ierr /= 0) return

          ! Evaluate overshoot boundary (_ob) parameters
          call star_eval_over_bdy_params(s, i, f0, k_ob, r_ob, D_ob, vc_ob, ierr)
          if (ierr /= 0) return

          ! Loop over cell faces, adding overshoot until D <= overshoot_D_min
          outward = s%top_conv_bdy(i)

          if (outward) then
              k_a = k_ob
              k_b = 1
              dk = -1
          else
              k_a = k_ob+1
              k_b = s%nz
              dk = 1
          endif

          if (f > 0.0_dp) then
              r_step = f*Hp_cb
          else
              r_step = 0.0_dp
          endif

          face_loop : do k = k_a, k_b, dk
              ! Evaluate the extended convective penetration factor
              r = s%r(k)
              if (outward) then
                  dr = r - r_ob
              else
                  dr = r_ob - r
              endif

              if (dr < r_step .AND. f > 0.0_dp) then  ! step factor
                  factor = 1.0_dp
              else
                  if ( f2 > 0.0_dp) then                ! exponential factor
                      factor = exp(-2.0_dp*(dr-r_step)/(f2*Hp_cb))
                  else
                      factor = 0.0_dp
                  endif
              endif

              ! Store the diffusion coefficient and velocity
              D(k) = (D0 + Delta0*D_ob)*factor
              vc(k) = (D0/D_ob + Delta0)*vc_ob*factor

              ! Check for early overshoot completion
              if (D(k) < s%overshoot_D_min) then
                  k_b = k
                  exit face_loop
              endif

          end do face_loop

          write(*,*)  'alpha_pen = ', dr_core_div_h
          !write(*,*)  'mcc = ', s% mass_conv_core

          if (DEBUG) then
              write(*,*) 'step exponential overshoot:'
              write(*,*) '  k_a, k_b   =', k_a, k_b
              write(*,*) '  r_a, r_b   =', s%r(k_a), s%r(k_b)
              write(*,*) '  r_ob, r_cb =', r_ob, r_cb
              write(*,*) '  Hp_cb      =', Hp_cb
          end if

      end subroutine extended_convective_penetration



      subroutine JermynAnders_penetration(s, id)
         use eos_def
         use star_lib
         use kap_def
         real(dp), parameter :: f = 0.86d0
         type (star_info), pointer :: s
         integer, intent(in) :: id
         real(dp), parameter :: xi = 0.6d0
         integer :: k, j, nz, ierr
         real(dp) :: Lint, delta_r, V_CZ, Favg, RHS, dr, h
         real(dp) :: Rho, T, logRho, logT, Pr
         ! real(dp), intent(out) :: dm_core, dr_core, dr_core_div_h, r_core, m_core, rho_core_top
         real(dp), dimension(num_eos_basic_results) :: res, dres_dlnRho, dres_dlnT
         real(dp) :: dres_dxa(num_eos_d_dxa_results,s% species)
         real(dp) :: kap, dlnkap_dlnRho, dlnkap_dlnT, frac_Type2
         real(dp) :: gradr(s%nz), grada(s%nz), gradL(s%nz)
         real(dp) :: kap_fracs(num_kap_fracs)

         nz = s%nz

         ! Recalculate gradR and gradA assuming the composition of the core.
         do j=1,nz
            ! Call the EOS with the composition of the convective core
            Rho = s%rho(j)
            T = s%T(j)
            logRho = log10(Rho)
            logT = log10(T)
            ierr = 0
            call star_get_eos( &
               id, 0, s%xa(:,nz), & ! k = 0 means not being called for a particular cell
               Rho, logRho, T, logT, &
               res, dres_dlnRho, dres_dlnT, &
               dres_dxa, ierr)
            grada(j) = res(i_grad_ad)

            ! Call the opacity with the composition of the convective core.
            ierr = 0
            call star_get_kap( &
               id, 0, s%zbar(nz), s%xa(:,nz), logRho, logT, &
               res(i_lnfree_e), dres_dlnRho(i_lnfree_e), dres_dlnT(i_lnfree_e), &
               res(i_eta), dres_dlnRho(i_eta), dres_dlnT(i_eta), &
               kap_fracs, kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)

            Pr = one_third*crad*T*T*T*T
            gradr(j) = s%Peos(j)*kap*s%L(j) / (16*pi*clight*s%m(j)*s%cgrav(j)*Pr)

         end do

         delta_r = 0d0
         V_CZ = 0d0
         Lint = 0d0

         ! Integrate over CZ
         do j=nz,1,-1
         if (gradr(j) < grada(j)) then
               ! Means we've hit a radiative zone
               m_core = s%m(j)
               r_core = s%r(j)
               rho_core_top = s%rho(j)
               h = s%scale_height(j)
               k = j
               exit
            end if

            dr = s%dm(j) / (4d0 * pi * pow2(s%r(j)) * s%rho(j))
            Lint = Lint + s%L_conv(j) * dr
            delta_r = delta_r + dr
            V_CZ = V_CZ + s%dm(j)/s%rho(j)

         end do

         ! Calculate target RHS
         Favg = Lint / V_CZ
         RHS = (1d0 - f) * V_CZ * Favg
         Lint = 0d0

         ! Integrate over RZ until we find the edge of the PZ
         delta_r = 0d0
         do j=min(nz,k+1),1,-1
            dr = s%dm(j) / (4d0 * pi * pow2(s%r(j)) * s%rho(j))
            delta_r = delta_r + dr
            Lint = Lint + (xi * f * 4d0 * pi * pow2(s%r(j)) * Favg + s%L(j) * (grada(j) / gradr(j) - 1d0)) * dr

            if (Lint > RHS) then
               dm_core = s%m(j) - m_core
               dr_core = delta_r
               dr_core_div_h = delta_r / h
               exit
            end if

         end do

      end subroutine JermynAnders_penetration





      end module run_star_extras
