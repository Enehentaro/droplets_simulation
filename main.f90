!---------------------------------------------------------------------------------
!     Simulation of viral droplets
!                                by KIYOTA OGURA(2021/1/10)
!     updated by IDA
!---------------------------------------------------------------------------------
PROGRAM MAIN
      !$ use omp_lib
      use drop_motion_mod
      use case_list_m
      implicit none

      character(7), parameter :: OS = 'Windows'

      integer, pointer :: n => n_time
      integer nc, nc_max
      character(50) start_date
      character(:), allocatable :: case_name
      real start_time
      !===========================================================================================
      !$OMP parallel
            !$OMP single
            !$ print *, "Num threads:", omp_get_num_threads()
            !$OMP end single
      !$OMP end parallel

      call case_check(num_case=nc_max) 

      do nc = 1, nc_max                         !実行数だけループ（通常1回）
            case_name = get_case_name(nc)
            call set_case_path(case_name)
            
            call make_directory                     !ディレクトリ作成
            
            call first_setting                        !条件TXTの読み込み、飛沫初期分布の計算など

            call initialization_droplet         !初期状態を代入

            call check_point                    !計算条件の確認および時刻計測のためのチェックポイント

            call read_flow_field                !流れ場の取得
            call preprocess_onFlowField         !流れ場の前処理

            print*,'*******************************************'
            print*,'             START step_loop               '
            print*,'*******************************************'

            DO n = n_start + 1, n_end           !ステップ数だけループ

                  call survival_check           !生存率に関する処理

                  call Calculation_Droplets     !飛沫の運動計算

                  call coalescence_check        !飛沫間の合体判定

                  if ((mod(n,interval) == 0)) call output             !出力

                  call update_flow_check        !流れ場の更新チェック

            END DO

            print*,'*******************************************'
            print*,'             END step_loop                 '
            print*,'*******************************************'

            call final_result       !最終結果出力

            call deallocation_flow  !流れ場配列解放
            call deallocation_droplet  !飛沫配列解放
            
      end do

      !プログラムここまで
      !===========================================================================================
      !以下、内部手続き

      contains

      subroutine check_point
            character(1) input
            character(10) d_start, t_start

            if(nc == 1) then
                  do
                        print*, 'Do you want to start the calculation? (y/n)'
                        read(5,*) input

                        select case(input)
                              case('y')
                                    exit

                              case('n')
                                    stop

                        end select

                  end do
            end if

            call cpu_time(start_time)
            call date_and_time(date = d_start, time = t_start)
            start_date = '[Start Date] ' &
            //d_start(1:4)//'/'//d_start(5:6)//'/'//d_start(7:8)//' ' &
            //t_start(1:2)//':'//t_start(3:4)//':'//t_start(5:)

      end subroutine

      subroutine make_directory
            use path_operator_m
            character(:), allocatable :: VTK_DIR, backup_DIR

            print*, '#', nc

            ! PATH_FlowDIR = path_list(nc)%path2FlowDIR
            ! FNAME_FMT = path_list(nc)%FlowFileName

            ! call check_FILE_GRID    !気流ファイルのタイプをチェック
      
            ! write(temperature,'(i3.3)') T
            ! write(humidity,'(i3.3)') RH

            ! i = len_trim(path_out_base)
            ! if(path_out_base(i:i) == '\') path_out_base(i:i) = ' '      !末尾が区切り文字であればこれを除去
            ! ! path_out = trim(path_out_base)//'_'//trim(temperature)//'_'//trim(humidity)//'\'

            VTK_DIR = path%VTK
            ! i = len_trim(path_out)
            ! if(path_out(i:i) /= '\') path_out(i+1:i+1) = '\'
            backup_DIR = path%backup

            select case(trim(OS))
                  case ('Linux')  !for_Linux
                        VTK_DIR =  replace_str(VTK_DIR, '\', '/' )
                        backup_DIR = replace_str(backup_DIR, '\', '/')
                        call system('mkdir -p -v '//VTK_DIR)
                        call system('mkdir -p -v '//backup_DIR)
                        ! call system('cp condition.txt '//path_out)

                  case ('Windows')  !for_Windows
                        VTK_DIR =  replace_str(VTK_DIR, '/', '\' )
                        backup_DIR = replace_str(backup_DIR, '/', '\')
                        call system('md '//VTK_DIR)
                        call system('md '//backup_DIR)
                        ! call system('copy condition.txt '//path_out)

                  case default
                        print*, 'OS ERROR', OS
                        stop
                        
            end select

      end subroutine make_directory

      subroutine output
            print*, start_date
            print*, 'Now_Step_Time=', dimensional_time(n), '[sec]'
            print*, 'Number of floating', drop_counter('floating')
            call output_droplet
      end subroutine output

      subroutine final_result
            integer n_unit
            real end_time
            character(50) end_date
            character(10) d_end, t_end
            
            call cpu_time(end_time)
            call date_and_time(date = d_end, time = t_end)

            end_date = '[ END  Date] ' &
                        //d_end(1:4)//'/'//d_end(5:6)//'/'//d_end(7:8)//' ' &
                        //t_end(1:2)//':'//t_end(3:4)//':'//t_end(5:)
            print*, start_date
            print*, end_date

            open(newunit=n_unit, FILE= case_name//'/final_result.txt',STATUS='REPLACE')
                  write(n_unit,*)'*******************************************'
                  write(n_unit,*)'*                                         *'
                  write(n_unit,*)'*             Final Results               *'
                  write(n_unit,*)'*                                         *'
                  write(n_unit,*)'*******************************************'
                  write(n_unit,'()')
                  write(n_unit,*) start_date
                  write(n_unit,*) end_date
                  write(n_unit, *) 'Erapsed Time =', end_time - start_time, '[sec]'
                  write(n_unit, *) 'Cost of Calc =', &
                        (end_time - start_time) / (dimensional_time(n_end) - dimensional_time(n_start)), '[sec/sec]'
                  write(n_unit,'(A)') '======================================================='
                  write(n_unit, '(A12, 2(F20.8,2x,A1))') 'Time[sec] =', dimensional_time(n_start), '-', dimensional_time(n_end)
                  write(n_unit, '(A12, 2(I20,2x,A1))') 'Step =', n_start, '-', n_end !計算回数
                  write(n_unit,'(A)') '======================================================='
                  write(n_unit,'(A15, I20)') 'num_droplets =', drop_counter('total')
                  write(n_unit,'(A15, I20)') 'alive =', drop_counter('floating')
                  write(n_unit,'(A15, I20)') 'death =', drop_counter('death') !生存率で消滅
                  write(n_unit,'(A15, I20)') 'coalescence =', drop_counter('coalescence') !生存率で消滅
                  write(n_unit,'(A15, I20)') 'adhesion =', drop_counter('adhesion') !付着したすべてのウイルス数
                  write(n_unit,'(A)') '======================================================='
                  write(n_unit,'(A15, F10.2)') 'Temp [degC] =', environment('Temperature')
                  write(n_unit,'(A15, F10.2)') 'RH [%] =', environment('Relative Humidity')
                  write(n_unit, *) 'Used FlowFile : ', trim(PATH_FlowDIR), trim(FNAME_FMT)

            close(n_unit)
            
      end subroutine final_result

END PROGRAM MAIN