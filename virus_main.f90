!---------------------------------------------------------------------------------
!     Simulation of viral drplets
!                                by KIYOTA OGURA(2021/1/10)
!     updated by IDA (2021/07/13)

!ToDo:      to modify module
!
!---------------------------------------------------------------------------------
      include 'csv_reader.f90'
      include 'cases_reader.f90'
      include 'virus_mod.f90'
      include 'flow_field.f90'
      include 'motion_mod.f90'
!*******************************************************************************************
PROGRAM MAIN
      use motion_mod
      use cases_reader
      implicit none

      character(7), parameter :: OS = 'Linux'!'Windows'

      integer n, vn, vnf, vfloat, nc, step_air, nc_max
      real nowtime
      character(20) :: d_start, d_stop, t_start, t_stop
!===========================================================================================
      call date_and_time(date = d_start, time = t_start)
      print*,'date = ', trim(d_start), ' time = ', trim(t_start)

      call input_condition

      call set_initial_position

      nc_max = check_cases(PATH_AIR)

      if((nc_max > 1).and.(restart >= 1)) then
            print*, 'ERROR_program:'
            print*, 'Remove '//trim(FNAME_FMT)//' or Set restart_No.= 0'
            stop
      end if

      do nc = 1, nc_max

            call reset_status

            call set_path
            
            call Set_Coefficients
      
            call initial_virus(restart)

            if(restart > 0) then
                  n_start = restart
            else
                  n_start = 0
                  call writeout(n_start)
            end if
            
            call read_nextcell
            call read_flow_field(n_start)

            print*,'*******************************************'
            print*,'             START step_loop               '
            print*,'*******************************************'

            DO n = n_start + 1, n_end

                  nowtime = real(n*dt*L_chara/U_chara)  !現在ステップ実時刻[sec]

                  call survival_check(n)

                  call set_vn_trans(vfloat)

                  !$omp parallel do private(vn)
                  DO vnf = 1, vfloat !浮遊粒子に対してのみループ
                        vn = vn_trans(vnf)
                        call evaporation(vn)
                        call VirusCalculation(vn)
                        call update_status(vn)
                  END DO
                  !$omp end parallel do 

                  if ((mod(n,interval) == 0)) then
                        print*, 'Startdate = ', trim(d_start), ' time = ', trim(t_start)
                        print*, 'Now_Step_Time=', nowtime, '[sec]'
                        print*, 'X Position=',crd_drp(1,1,1),crd_drp(1,2,1),crd_drp(1,3,1)
                        print*, 'Number of floating', count(adhesion==0)
                        print*, 'Number of calling Nearest_Cell_Serch=', num_NCS
                        num_NCS = 0
                        call writeout(n)
                  end if

                  Step_air = int(dble(n)*Rdt)          !気流計算における経過ステップ数に相当
                  if((mod(Step_air, INTERVAL_FLOW) == 0).and.(INTERVAL_FLOW > 0)) call read_flow_field(n)

            END DO

            print*,'*******************************************'
            print*,'             END step_loop                 '
            print*,'*******************************************'

            call date_and_time(date = d_stop, time = t_stop)
            print*,'date = ', d_start, ' time = ', t_start
            print*,'date = ', d_stop,  ' time = ', t_stop

            call final_result

            call deallocation_flow
            
      end do

      contains

      subroutine set_path
            character(20) :: temperature, humidity
            integer i

            if(nc_max > 1) then

                  T = get_temperature(nc)
                  RH = get_humidity(nc)
      
                  call set_case_path(PATH_AIR, nc)
                  call set_case_path2(path_out_base, nc)

            end if

            call set_dir_from_path(PATH_AIR, PATH_AIR, FNAME_FMT)

            call set_FILE_TYPE

            ! if(nc_max > 1) path_out_base = '..\' // trim(HEAD_AIR) // '_virus\'
      
            print*, 'T =', T, 'degC'
            print*, 'RH =', RH, '%'
      
            write(temperature,'(i3.3)') T
            write(humidity,'(i3.3)') RH

            i = len_trim(path_out_base)
            if(path_out_base(i:i) == '\') path_out_base(i:i) = ' '      !末尾が区切り文字であればこれを除去
            path_out = trim(path_out_base)//'_'//trim(temperature)//'_'//trim(humidity)//'\'

            select case(trim(OS))
                  case ('Linux')  !for_Linux
                        call replace_str(path_out, '\', '/' )
                        call replace_str(PATH_AIR, '\', '/' )
                        call system('mkdir -p -v '//path_out)
                        call system('cp condition_virus.txt '//path_out)

                  case ('Windows')  !for_Windows
                        call system('md '//path_out)
                        call system('copy condition_virus.txt '//path_out)

                  case default
                        print*, 'OS ERROR', OS
                        stop
                        
            end select

            print*, 'Output_Path=', path_out

            

      end subroutine set_path

      subroutine final_result
            integer n_unit

            open(newunit=n_unit, FILE= trim(path_out)//'/'//'final_result.txt',STATUS='REPLACE')
                  write(n_unit,*)'date = ', d_start, ' time = ', t_start
                  write(n_unit,*)'date = ', d_stop,  ' time = ', t_stop
                  WRITE(n_unit,*) '======================================================='
                  WRITE(n_unit,*) 'alive =', count(adhesion==0)
                  WRITE(n_unit,*) 'Step =', n !計算回数
                  write(n_unit,*) 'TIME[sec]=', nowtime
                  WRITE(n_unit,*) 'vndeath =',count(adhesion==-1) !生存率で消滅
                  WRITE(n_unit,*) 'vnadhesion =', count(adhesion >= 1) !付着したすべてのウイルス数
                  WRITE(n_unit,*) 'vnwall =',count(adhesion==1),count(adhesion==2) !かべで止まったウイルス数
                  WRITE(n_unit,*) 'vnfloor =',count(adhesion==3),count(adhesion==4) !ゆかで止まったウイルス数
                  WRITE(n_unit,*) 'vnout =',count(adhesion==5) !流出したウイルス数
                  write(n_unit,*) 'vncombination =',count(adhesion==-2) !結合回数
                  WRITE(n_unit,*) '======================================================='
            close(n_unit)
            
      end subroutine final_result


END PROGRAM MAIN