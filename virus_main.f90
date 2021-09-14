!---------------------------------------------------------------------------------
!     Simulation of viral drplets
!                                by KIYOTA OGURA(2021/1/10)
!     updated by IDA (2021/07/13)

!ToDo:      to modify module
!
!---------------------------------------------------------------------------------
      include 'virus_mod.f90'
      include 'flow_field.f90'
      include 'motion_mod.f90'
!*******************************************************************************************
PROGRAM MAIN
      use motion_virus
      implicit none
      integer n, vn, vnf, vfloat, PN, ios, case_num, num_programs
      integer, allocatable :: program_values(:,:)
      real nowtime
      double precision Step_air
      character(20) :: d_start, d_stop, t_start, t_stop
!===========================================================================================
      call date_and_time(date = d_start, time = t_start)
      print*,'date = ', trim(d_start), ' time = ', trim(t_start)

      call input_condition

      call read_program

      do PN = 1, num_programs

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

                  Step_air = dble(n) * Rdt          !気流計算における経過ステップ数に相当
                  if((mod(Step_air, dble(interval_flow)) == 0.0d0).and.(interval_flow > 0)) &
                        call read_flow_field(n)

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

      subroutine read_program
            character(11) FNAME
            integer n_unit

            FNAME = 'program.csv'

            open(newunit=n_unit, iostat=ios, file=FNAME, status='old')
            close(n_unit)
            if(ios /= 0) then
                  print*, 'Normal_program'
                  num_programs = 1
            else
                  call csv_reader_int(FNAME, program_values, 3)
                  num_programs = size(program_values, dim=2)
            end if

            if((num_programs > 1).and.(restart>=1)) then
                  print*, 'ERROR_program:'
                  print*, 'Remove program.csv or Set restart_No.= 0'
                  stop
            end if

      end subroutine read_program

      subroutine set_path
            character(20) :: temperature, humidity, a

            if(allocated(program_values)) then
                  case_num = program_values(1,PN)
                  T = program_values(2,PN)
                  RH = program_values(3,PN)
      
                  write(path_out_base,'("ACAP",i3.3,"_virus")') case_num
                  write(head_out,'("ac",i3.3,"_")') case_num
      
                  if (mod(case_num,10)==0) then
                        write(PATH_AIR,'("ACAP0",i2.2)') case_num/10
                        write(HEAD_AIR,'("ACAP0",i2.2)') case_num/10
                  else
                        write(PATH_AIR,'("ACAP0",i2.2,"-", i1.1)') case_num/10, case_num - case_num/10*10
                        write(HEAD_AIR,'("ACAP0",i2.2,"-", i1.1)') case_num/10, case_num - case_num/10*10
                  end if
            end if
      
            print*, 'T =', T, 'degC'
            print*, 'RH =', RH, '%'
      
            write(temperature,'(i3.3)') T
            write(humidity,'(i3.3)') RH
            write(a,'(i1.1)') PN
      
            path_out = trim(path_out_base)//'_'//trim(temperature)//'_'//trim(humidity)
            print*, 'Output_Path=', path_out
            ! if(num_programs > 1) path_out = trim(path_out)//'_'//trim(a)
            call system('mkdir -p -v '//path_out)  !サブルーチンsystem：引数文字列をコマンドとして実行する
            call system('cp condition_virus.txt '//path_out) !Windows:`copy`, Linux:`cp`

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