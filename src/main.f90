!---------------------------------------------------------------------------------
!     Simulation of viral droplets
!                                by KIYOTA OGURA(2021/1/10)
!     updated by IDA
!---------------------------------------------------------------------------------
PROGRAM MAIN
      !$ use omp_lib
      use dropletGroup_m
      use caseNameList_m
      implicit none
      integer, pointer :: n => n_time, nc => nowCase
      integer nc_max
      character(50) start_date
      real start_time
      !===========================================================================================
      !$OMP parallel
            !$OMP single
            !$ print *, "Num threads:", omp_get_num_threads()
            !$OMP end single
      !$OMP end parallel

      call case_check(num_case=nc_max) 

      do nc = 1, nc_max                         !実行数だけループ（通常1回）
            
            call create_CaseDirectory                     !ディレクトリ作成
            
            call firstSet_mainDroplet          !条件TXTの読み込み、飛沫初期分布の計算など

            call dropletManagement       !外部サブルーチンによる管理

            call output_initialDroplet

            call check_point                    !計算条件の確認および時刻計測のためのチェックポイント

            print '("*******************************************")'
            print '("            START step_loop                ")'
            print '("*******************************************")'

            DO n = n_start + 1, n_end           !ステップ数だけループ

                  call mainDroplet%adhesion_check()

                  call mainDroplet%survival_check()           !生存率に関する処理

                  call mainDroplet%coalescence_check()        !飛沫間の合体判定

                  call mainDroplet%Calculation_Droplets()     !飛沫の運動計算

                  call dropletManagement       !外部サブルーチンによる管理

                  if ((mod(n,interval) == 0)) call output             !出力

                  call update_flow_check        !流れ場の更新チェック

            END DO

            print '("*******************************************")'
            print '("             END step_loop                 ")'
            print '("*******************************************")'

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

      subroutine create_CaseDirectory
            use path_operator_m
            character(:), allocatable :: caseName

            caseName = get_caseName()

            print*, '#', nc, '[',caseName,']'

            call make_directory(caseName//'/VTK')
            call make_directory(caseName//'/backup')
            
      end subroutine create_CaseDirectory

      subroutine output
            use terminalControler_m

            print*, start_date
            print*, 'Now_Step_Time =', Time_onSimulation(n, dimension=.true.), '[sec]'
            print*, '# floating :', mainDroplet%dropletCounter('floating')
            if(refCellSearchInfo('FalseRate') >= 1) print*, '# searchFalse :', refCellSearchInfo('NumFalse')
            call output_mainDroplet(initial=.false.)
            print '("====================================================")'
            call reset_formatTC

      end subroutine output

      subroutine output_ResultSummary
            integer n_unit, cnt
            real end_time
            character(50) end_date, fname
            character(10) d_end, t_end
            logical existance
            double precision TimeStart, TimeEnd
            character(:), allocatable :: caseName
            
            call cpu_time(end_time)
            call date_and_time(date = d_end, time = t_end)

            end_date = '[ END  Date] ' // DateAndTime_string(d_end, t_end)
            print*, start_date
            print*, end_date

            caseName = get_caseName()

            fname = caseName//'/ResultSummary.txt'
            inquire(file=fname, exist=existance)
            cnt = 0
            do while(existance)
                  cnt = cnt + 1
                  write(fname,'("'//caseName//'/ResultSummary_", i0, ".txt")') cnt
                  inquire(file=fname, exist=existance)
            end do

            TimeStart = Time_onSimulation(n_start, dimension=.true.)
            TimeEnd = Time_onSimulation(n_end, dimension=.true.)

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
                        (end_time - start_time) / (TimeEnd - TimeStart), '[sec/sec]'
                  write(n_unit,'(A)') '======================================================='
                  write(n_unit, '(A18, 2(F15.3,2x,A))') 'Time [sec] =', TimeStart, '-', TimeEnd
                  write(n_unit, '(A18, 2(I15,2x,A))') 'Step =', n_start, '-', n_end !計算回数
                  write(n_unit,'(A18, I15)') 'OutputInterval =', interval
                  write(n_unit,'(A)') '======================================================='
                  write(n_unit,'(A18, I15)') '#Droplets =', mainDroplet%dropletCounter('total')
                  write(n_unit,'(A18, I15)') 'floating =', mainDroplet%dropletCounter('floating')
                  write(n_unit,'(A18, I15)') 'death =', mainDroplet%dropletCounter('death') !生存率で消滅
                  write(n_unit,'(A18, I15)') 'coalescence =', mainDroplet%dropletCounter('coalescence') !生存率で消滅
                  write(n_unit,'(A18, I15)') 'adhesion =', mainDroplet%dropletCounter('adhesion') !付着したすべてのウイルス数
                  write(n_unit,'(A)') '======================================================='
                  write(n_unit,'(A18, F18.2)') 'Temp [degC] =', environment('Temperature')
                  write(n_unit,'(A18, F18.2)') 'RH [%] =', environment('RelativeHumidity')
                  write(n_unit,'(A18, 2X, A)') 'Used FlowFile :', trim(PATH_FlowDIR)//trim(FNAME_FMT)
                  write(n_unit, '(A18, 2(I15,2x,A))') 'SearchFalseInfo :', refCellSearchInfo('NumFalse'), &
                        ' (', refCellSearchInfo('FalseRate'), '%)'

            close(n_unit)
            
      end subroutine output_ResultSummary

      function DateAndTime_string(date, time) result(string)
            character(*), intent(in) :: date, time
            character(:), allocatable :: string

            string = date(1:4)//'/'//date(5:6)//'/'//date(7:8)//' ' &
                  //time(1:2)//':'//time(3:4)//':'//time(5:6)

      end function

END PROGRAM MAIN