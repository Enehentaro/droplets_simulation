!---------------------------------------------------------------------------------
!     Simulation of viral drplets
!                                by KIYOTA OGURA(2021/1/10)
!     updated by IDA (2021/07/13)
!
!---------------------------------------------------------------------------------
      include 'csv_reader.f90'
      include 'cases_reader.f90'
      include 'flow_field.f90'
      include 'equation_mod.f90'
      include 'drop_motion.f90'
!*******************************************************************************************
PROGRAM MAIN
      use drop_motion_mod
      use cases_reader
      implicit none

      character(7), parameter :: OS = 'Windows'

      integer n, nc, nc_max
      double precision Step_air, now_time
      character(20) :: d_start, d_stop, t_start, t_stop
!===========================================================================================
      call date_and_time(date = d_start, time = t_start)
      print*,'date = ', trim(d_start), ' time = ', trim(t_start)

      call pre_setting    !条件TXTの読み込み

      nc_max = check_cases(PATH_AIR)      !連続実行数の取得

      if((cases_read_flag).and.(num_restart >= 1)) then  !連続実行とリスタートを同時にするとどうなるのやら
            print*, 'WARNING:Continuous execution and restart happened at the same time.'
            print*, 'Recommend to Remove '//trim(FNAME_FMT)//' or to Set restart_No.= 0'
            ! stop
      end if

      do nc = 1, nc_max

            call set_path     !パスなどの整理

            call set_coeff_drdt(T, RH)  !温湿度依存の係数の設定

            call initialization_droplet
            
            call read_nextcell      !セルの隣接関係の取得
            call read_flow_field(n_start) !流れ場の取得

            print*,'*******************************************'
            print*,'             START step_loop               '
            print*,'*******************************************'

            DO n = n_start + 1, n_end

                  now_time = real_time(n)  !現在ステップ実時刻[sec]

                  call survival_check(n)  !生存率に関する処理

                  call Calculation_Droplets

                  if ((mod(n,interval) == 0)) then
                        call standard_output
                        call output(n)  !結果出力
                  end if

                  if(INTERVAL_FLOW > 0) then
                        Step_air = dble(n)*Rdt          !気流計算における経過ステップ数に相当
                        if(mod(Step_air, dble(INTERVAL_FLOW)) == 0.d0) call read_flow_field(n)   !流れ場の更新
                  end if

            END DO

            print*,'*******************************************'
            print*,'             END step_loop                 '
            print*,'*******************************************'

            call date_and_time(date = d_stop, time = t_stop)
            print*,'date = ', d_start, ' time = ', t_start
            print*,'date = ', d_stop,  ' time = ', t_stop

            call final_result

            call deallocation_flow  !配列解放
            
      end do

      contains

      subroutine set_path
            character(20) :: temperature, humidity
            integer i

            if(cases_read_flag) then

                  T = get_temperature(nc)
                  RH = get_humidity(nc)
      
                  PATH_AIR = get_case_path(nc)
                  path_out_base = get_case_path2(nc)

            end if

            call set_dir_from_path(PATH_AIR, PATH_AIR, FNAME_FMT)

            call set_FILE_TYPE
      
            print*, 'T =', T, 'degC'
            print*, 'RH =', RH, '%'
      
            write(temperature,'(i3.3)') T
            write(humidity,'(i3.3)') RH

            i = len_trim(path_out_base)
            if(path_out_base(i:i) == '\') path_out_base(i:i) = ' '      !末尾が区切り文字であればこれを除去
            path_out = trim(path_out_base)//'_'//trim(temperature)//'_'//trim(humidity)//'\'

            select case(trim(OS))
                  case ('Linux')  !for_Linux
                        path_out =  replace_str(path_out, '\', '/' )
                        PATH_AIR = replace_str(PATH_AIR, '\', '/' )
                        call system('mkdir -p -v '//path_out)
                        call system('cp condition.txt '//path_out)

                  case ('Windows')  !for_Windows
                        call system('md '//path_out)
                        call system('copy condition.txt '//path_out)

                  case default
                        print*, 'OS ERROR', OS
                        stop
                        
            end select

            print*, 'Output_Path=', path_out

      end subroutine set_path

      subroutine standard_output
            print*, 'Startdate = ', trim(d_start), ' time = ', trim(t_start)
            print*, 'Now_Step_Time=', now_time, '[sec]'
            print*, 'Number of floating', get_drop_info('floating')
            print*, 'Number of calling Nearest_Cell_Serch=', num_NCS
            num_NCS = 0
      end subroutine standard_output

      subroutine final_result
            integer n_unit

            open(newunit=n_unit, FILE= trim(path_out)//'/'//'final_result.txt',STATUS='REPLACE')
                  write(n_unit,*)'date = ', d_start, ' time = ', t_start
                  write(n_unit,*)'date = ', d_stop,  ' time = ', t_stop
                  WRITE(n_unit,*) '======================================================='
                  WRITE(n_unit,*) 'alive =', get_drop_info('floating')
                  WRITE(n_unit,*) 'Step =', n !計算回数
                  write(n_unit,*) 'TIME[sec]=', now_time
                  WRITE(n_unit,*) 'vndeath =',get_drop_info('death') !生存率で消滅
                  WRITE(n_unit,*) 'vnadhesion =', get_drop_info('adhesion') !付着したすべてのウイルス数
                  WRITE(n_unit,*) '======================================================='
            close(n_unit)
            
      end subroutine final_result


END PROGRAM MAIN