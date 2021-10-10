!---------------------------------------------------------------------------------
!     Simulation of viral droplets
!                                by KIYOTA OGURA(2021/1/10)
!     updated by IDA
!---------------------------------------------------------------------------------
PROGRAM MAIN
      !$ use omp_lib
      use drop_motion_mod
      use cases_reader
      implicit none

      character(7), parameter :: OS = 'Windows'

      integer n, nc, nc_max
      double precision Step_air
      character(50) start_date
      !===========================================================================================
      !$OMP parallel
            !$OMP single
            !$ print *, "Num threads:", omp_get_num_threads()
            !$OMP end single
      !$OMP end parallel

      call first_setting                        !条件TXTの読み込み、飛沫初期分布の計算など

      nc_max = check_cases(PATH_AIR)            !連続実行数の取得（通常1回）

      call check_point                          !計算条件の確認のためのチェックポイント

      do nc = 1, nc_max                         !連続実行数だけループ（通常1回）

            call set_path                       !パスなどの整理

            call set_coeff_drdt(T, RH)          !温湿度依存の係数の設定

            call initialization_droplet         !初期状態を代入

            call read_flow_field(n_start)       !流れ場の取得
            call preprocess_onFlowField         !流れ場の前処理

            print*,'*******************************************'
            print*,'             START step_loop               '
            print*,'*******************************************'

            DO n = n_start + 1, n_end           !ステップ数だけループ

                  call survival_check(n)        !生存率に関する処理

                  call Calculation_Droplets     !飛沫の運動計算

                  call coalescence_check(n)     !飛沫間の合体判定

                  if ((mod(n,interval) == 0)) then
                        call standard_output
                        call output(n)          !結果出力
                  end if

                  if(INTERVAL_FLOW > 0) then
                        Step_air = dble(n)*Rdt          !気流計算における経過ステップ数に相当
                        if(mod(Step_air, dble(INTERVAL_FLOW)) == 0.d0) call read_flow_field(n)   !流れ場の更新
                  end if

            END DO

            print*,'*******************************************'
            print*,'             END step_loop                 '
            print*,'*******************************************'

            call final_result       !最終結果出力

            call deallocation_flow  !流れ場配列解放
            
      end do

      !===========================================================================================
      !以下、内部手続き

      contains

      subroutine check_point
            character(1) input
            character(10) d_start, t_start

            do
                  print*, 'Do you want to start the calculation? (y/n)'
                  read(5,*) input

                  select case(input)
                        case('y')
                              call date_and_time(date = d_start, time = t_start)
                              start_date = '[Start Date] ' &
                              //d_start(1:4)//'/'//d_start(5:6)//'/'//d_start(7:8)//' ' &
                              //t_start(1:2)//':'//t_start(3:4)//':'//t_start(5:)
                              exit

                        case('n')
                              stop

                  end select

            end do

      end subroutine

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

            call check_FILE_GRID
      
            print*, 'T =', T, 'degC'
            print*, 'RH =', RH, '%'
      
            write(temperature,'(i3.3)') T
            write(humidity,'(i3.3)') RH

            i = len_trim(path_out_base)
            if(path_out_base(i:i) == '\') path_out_base(i:i) = ' '      !末尾が区切り文字であればこれを除去
            path_out = trim(path_out_base)//'_'//trim(temperature)//'_'//trim(humidity)//'\'
            path_backup = 'backup\'

            select case(trim(OS))
                  case ('Linux')  !for_Linux
                        path_out =  replace_str(path_out, '\', '/' )
                        path_backup = replace_str(path_backup, '\', '/')
                        PATH_AIR = replace_str(PATH_AIR, '\', '/' )
                        call system('mkdir -p -v '//trim(path_out)//trim(path_backup))
                        call system('cp condition.txt '//path_out)

                  case ('Windows')  !for_Windows

                        call system('md '//trim(path_out)//trim(path_backup))
                        call system('copy condition.txt '//path_out)

                  case default
                        print*, 'OS ERROR', OS
                        stop
                        
            end select

            print*, 'Output_Path=', path_out

      end subroutine set_path

      subroutine standard_output
            print*, start_date
            print*, 'Now_Step_Time=', real_time(n), '[sec]'
            print*, 'Number of floating', get_drop_info('floating')
      end subroutine standard_output

      subroutine final_result
            integer n_unit
            character(50) end_date
            character(10) d_end, t_end
            
            call date_and_time(date = d_end, time = t_end)

            end_date = '[ END  Date] ' &
                        //d_end(1:4)//'/'//d_end(5:6)//'/'//d_end(7:8)//' ' &
                        //t_end(1:2)//':'//t_end(3:4)//':'//t_end(5:)
            print*, start_date
            print*, end_date

            open(newunit=n_unit, FILE= trim(path_out)//'final_result.txt',STATUS='REPLACE')
                  write(n_unit,*) start_date
                  write(n_unit,*) end_date
                  write(n_unit,*) '======================================================='
                  write(n_unit,'(A20, F20.8)') 'TIME[sec] =', real_time(n_end)
                  write(n_unit,'(A20, I20)') 'Step =', n_end !計算回数
                  write(n_unit,'(A20, I20)') 'alive =', get_drop_info('floating')
                  write(n_unit,'(A20, I20)') 'death =', get_drop_info('death') !生存率で消滅
                  write(n_unit,'(A20, I20)') 'coalescence =', get_drop_info('coalescence') !生存率で消滅
                  write(n_unit,'(A20, I20)') 'adhesion =', get_drop_info('adhesion') !付着したすべてのウイルス数
                  write(n_unit,*) '======================================================='
            close(n_unit)
            
      end subroutine final_result

END PROGRAM MAIN