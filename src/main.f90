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

            call set_initialDroplet         !初期状態を代入

            call check_point                    !計算条件の確認および時刻計測のためのチェックポイント

            call update_FlowField(first=.true.)                !流れ場の取得

            print*,'*******************************************'
            print*,'             START step_loop               '
            print*,'*******************************************'

            DO n = n_start + 1, n_end           !ステップ数だけループ

                  call management_droplet       !外部サブルーチンによる管理

                  call adhesion_check

                  call survival_check           !生存率に関する処理

                  call coalescence_check        !飛沫間の合体判定

                  call Calculation_Droplets     !飛沫の運動計算

                  if ((mod(n,interval) == 0)) call output             !出力

                  call update_flow_check        !流れ場の更新チェック

            END DO

            print*,'*******************************************'
            print*,'             END step_loop                 '
            print*,'*******************************************'

            call output_ResultSummary       !最終結果出力

            call deallocation_flow  !流れ場配列解放
            
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
            start_date = '[Start Date] ' // DateAndTime_string(d_start, t_start)
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
            print*, 'Now_Step_Time =', dimensional_time(n), '[sec]'
            print*, 'Number of floating :', drop_counter('floating')
            call output_droplet
      end subroutine output

      subroutine output_ResultSummary
            integer n_unit, cnt
            real end_time
            character(50) end_date, fname
            character(10) d_end, t_end
            logical existance
            
            call cpu_time(end_time)
            call date_and_time(date = d_end, time = t_end)

            end_date = '[ END  Date] ' // DateAndTime_string(d_end, t_end)
            print*, start_date
            print*, end_date

            fname = path%DIR//'ResultSummary.txt'
            inquire(file=fname, exist=existance)
            cnt = 0
            do while(existance)
                  cnt = cnt + 1
                  write(fname,'("'//path%DIR//'ResultSummary_", i0, ".txt")') cnt
                  inquire(file=fname, exist=existance)
            end do

            open(newunit=n_unit, file=fname, status='new')
                  write(n_unit,*)'*******************************************'
                  write(n_unit,*)'*                                         *'
                  write(n_unit,*)'*             Result Summary              *'
                  write(n_unit,*)'*                                         *'
                  write(n_unit,*)'*******************************************'
                  write(n_unit,'(A)') '======================================================='
                  write(n_unit,*) start_date
                  write(n_unit,*) end_date
                  write(n_unit, '(A18, F15.3, 2X, A)') 'Erapsed Time =', end_time - start_time, '[sec]'
                  write(n_unit, '(A18, F15.3, 2X, A)') 'Cost of Calc =', &
                        (end_time - start_time) / (dimensional_time(n_end) - dimensional_time(n_start)), '[sec/sec]'
                  write(n_unit,'(A)') '======================================================='
                  write(n_unit, '(A18, 2(F15.3,2x,A))') 'Time [sec] =', dimensional_time(n_start), '-', dimensional_time(n_end)
                  write(n_unit, '(A18, 2(I15,2x,A))') 'Step =', n_start, '-', n_end !計算回数
                  write(n_unit,'(A18, I15)') 'OutputInterval =', interval
                  write(n_unit,'(A)') '======================================================='
                  write(n_unit,'(A18, I15)') '#Droplets =', drop_counter('total')
                  write(n_unit,'(A18, I15)') 'floating =', drop_counter('floating')
                  write(n_unit,'(A18, I15)') 'death =', drop_counter('death') !生存率で消滅
                  write(n_unit,'(A18, I15)') 'coalescence =', drop_counter('coalescence') !生存率で消滅
                  write(n_unit,'(A18, I15)') 'adhesion =', drop_counter('adhesion') !付着したすべてのウイルス数
                  write(n_unit,'(A)') '======================================================='
                  write(n_unit,'(A18, F18.2)') 'Temp [degC] =', environment('Temperature')
                  write(n_unit,'(A18, F18.2)') 'RH [%] =', environment('Relative Humidity')
                  write(n_unit,'(A18, 2X, A)') 'Used FlowFile :', trim(PATH_FlowDIR)//trim(FNAME_FMT)

            close(n_unit)
            
      end subroutine output_ResultSummary

      function DateAndTime_string(date, time) result(string)
            character(*), intent(in) :: date, time
            character(:), allocatable :: string

            string = date(1:4)//'/'//date(5:6)//'/'//date(7:8)//' ' &
                  //time(1:2)//':'//time(3:4)//':'//time(5:6)

      end function

END PROGRAM MAIN