module motion_mod
      use virus_mod
      use flow_field
      implicit none
      integer restart, n_start, n_end
      integer, private :: LoopS, LoopF, OFFSET
      integer :: num_NCS=0    !NearestCell探索を行った回数のカウンター
      double precision Rdt    !飛沫計算と気流計算の時間間隔の比
      double precision, private :: Re     !レイノルズ数
      integer, allocatable, private :: nearcell(:)                 !飛沫近接要素ID
      integer, allocatable, private :: adhes_bound(:)              !飛沫付着境界面ID

      contains

      subroutine input_condition
            double precision DTa
            integer L, n_unit, num_virus

            OPEN(newunit=n_unit,FILE='condition_virus.txt',STATUS='OLD')
                  read(n_unit,'()')
                  read(n_unit,*) restart
                  read(n_unit,'()')
                  read(n_unit,*) n_end
                  read(n_unit,'()')
                  read(n_unit,*) dt
                  read(n_unit,'()')
                  read(n_unit,*) interval
                  read(n_unit,'()')
                  read(n_unit,'(A)') path_out_base
                  read(n_unit,*) head_out
                  read(n_unit,'()')
                  read(n_unit,*) T
                  read(n_unit,*) RH
                  read(n_unit,'()')
                  read(n_unit,*) num_virus
                  read(n_unit,'()')
                  read(n_unit,*) (center_posi(L), L=1,3)
                  read(n_unit,*) (width_posi(L), L=1,3)
                  read(n_unit,'()')
                  read(n_unit,*) (direction_g(L), L=1,3)
                  
                  read(n_unit,'()')
            
                  read(n_unit,'()')
                  read(n_unit,'(A)') PATH_AIR
                  read(n_unit,'()')
                  read(n_unit,*) DTa
                  read(n_unit,'()')
                  read(n_unit,*) OFFSET
                  read(n_unit,'()')
                  read(n_unit,*) INTERVAL_FLOW
                  read(n_unit,'()')
                  read(n_unit,*) LoopS
                  read(n_unit,*) LoopF
                  read(n_unit,'()')
                  read(n_unit,*) L_chara
                  read(n_unit,*) U_chara
                  read(n_unit,*) Roh_chara
                  read(n_unit,*) Mu_chara
  
            CLOSE(n_unit)

            Rdt = dt/DTa                       !データ読み込み時に時間軸合わせるパラメータ
            Re = (Roh_chara*U_chara*L_chara)/Mu_chara

            call allocation_virus(num_virus)
            allocate(nearcell(num_virus), source=0)
            allocate(adhes_bound(num_virus), source=0)
  
            print*, 'restart =',restart
            print*, 'n_end =',n_end
            print*, 'interval =',interval
            print*, 'interval_air =',INTERVAL_FLOW
            print*, 'loop=',loops,loopf
            print*, 'dt =',dt
            print*, 'Rdt', Rdt
            print*, 'Re =',Re

            print*, 'PATH_AIR=', PATH_AIR

  
      end subroutine input_condition

      subroutine VirusCalculation(vn) !CALCULATE VIRUS POSITION
            integer, intent(in) :: vn
            double precision  :: X(3), V(3)
            integer NCN
            double precision randble
            logical stopflag
        
            X(:) = crd_drp(:,vn,1)
            V(:) = vel_drp(:,vn,1)    
            NCN = nearcell(vn)     !前回参照セルを代入
        
            if(NCN == 0) then   !参照セルが見つかっていない（＝初期ステップ）
                  NCN = nearest_cell(X)    
                  print*, 'FirstNCN:', NCN
                  nearcell(:) = NCN !全粒子が同一セル参照と仮定して時間短縮を図る
        
            else
                  NCN = nearer_cell(X, NCN)
                  if (NCN == 0) then
                        print*, 'NCN_ERROR:', vn, NCN
                        stop
                  end if
                  if (.not.nearcell_check(X(:), NCN)) NCN = nearest_cell(X)
        
            end if
        
            nearcell(vn) = NCN    !結果を参照セル配列に記憶
 
            stopflag = .false.
            if (NoB(NCN) >= 1) stopflag = adhesion_check(vn, NCN)

            call area_check(vn, stopflag)
        
            if (stopflag) then
              
                  ! if ((X(1)>1.59d0.and.X(1)<2.01d0).and.&
                  !       (X(2)>0.25d0.and.X(2)<0.35d0).and.&
                  !       (X(3)>0.01d0.and.X(3)<0.75d0)) then
                  !       adhesion(vn) = 2       !AP_air_cleaner
                  ! else
                  !       adhesion(vn) = 1
                  ! end if
        
                  if ((X(1)>0.29d0.and.X(1)<0.31d0).and.&
                          (X(2)>2.94d0.and.X(2)<3.36d0).and.&
                          (X(3)>0.01d0.and.X(3)<0.74d0)) then
        
                        adhesion(vn) = 2       !ACAP_air_cleaner_left
        
                  else if ((X(1)>1.58d0.and.X(1)<2.02d0).and.&
                          (X(2)>5.99d0.and.X(2)<6.01d0).and.&
                          (X(3)>0.01d0.and.X(3)<0.74d0)) then
        
                        adhesion(vn) = 2       !ACAP_air_cleaner_oposit
        
                  else if ((X(1)>1.58d0.and.X(1)<2.02d0).and.&
                          (X(2)>0.29d0.and.X(2)<0.31d0).and.&
                          (X(3)>0.01d0.and.X(3)<0.74d0)) then
        
                        adhesion(vn) = 2       !ACAP_air_cleaner_under
        
        
                  else if ((X(1)>1.43d0.and.X(1)<2.18d0).and.&
                          (X(2)>0.01d0.and.X(2)<0.22d0).and.&
                          (X(3)>2.15d0.and.X(3)<2.18d0)) then
        
                        call random_number(randble)
        
                        if (randble < 0.015d0) then
                              adhesion(vn) = 3       !ACAP_air_conditioner
        
                        else
                              crd_drp(1,vn,2) = X(1)
                              crd_drp(2,vn,2) = X(2) + 0.23d0
                              crd_drp(3,vn,2) = X(3) - 0.26d0
        
                              vel_drp(:,vn,2) = 0.0d0     !速度をゼロに
        
                              return
              
                        end if
        
        
                  else
                        adhesion(vn) = 1

                  end if
        
                  vel_drp(:,vn,2) = 0.0d0     !速度をゼロに
                  crd_drp(:,vn,2) = X(:)
        
            else
        
                  vel_drp(:,vn,2) = get_velocity(V(:), VELC(:, NCN), radius(vn,2))
                  
                  crd_drp(:,vn,2) = X(:) + (V(:) + vel_drp(:,vn,2))*0.5d0*dt
              
            end if
        
            
      end subroutine VirusCalculation

      !*******************************************************************************************
      function get_velocity(vel_d, vel_a, radius_d) result(vel_d_next)
            !*******************************************************************************************
            !=====================================================================================
            double precision, intent(in) :: vel_d(3), vel_a(3), radius_d
            double precision speed_r, Re_d, Cd, Coefficient, vel_d_next(3)
            !=====================================================================================

            speed_r = norm2(vel_a(:) - vel_d(:))
            Re_d = (speed_r * 2.0d0*radius_d) * Re + 1.d-9  !ゼロ割回避のため、小さな値を足す

            Cd = (24.0d0/Re_d)*(1.0d0 + 0.15d0*(Re_d**0.687d0))

            Coefficient = (3.0d0*Cd*gumma*speed_r)/(8.0d0*radius_d)

            vel_d_next(:) = ( vel_d(:) + ( G(:) + Coefficient*vel_a(:) )*dt ) &
                              / ( 1.0d0 + Coefficient*dt )

      end function get_velocity
      !----------------------------------------------------------------------------------


    !*******************************************************************************************
      integer function nearest_cell(X) !粒子vnに最も近いセルNCNの探索
            double precision, intent(in) :: X(3)
            integer II, IIMX
            double precision, allocatable :: distance(:)
            !=====================================================================================
            num_NCS = num_NCS +1

            IIMX = size(CENC, dim=2)

            allocate(distance(IIMX))
            !↓↓↓↓　一番近いセル中心の探索
            !$omp parallel do
            DO II = 1,IIMX
                  distance(II) = norm2(CENC(:,II) - X(:))
            END DO
            !$omp end parallel do 
            !↑↑↑↑
            
            nearest_cell = minloc(distance, dim=1)   !最小値インデックス
            
      end function nearest_cell

    !----------------------------------------------------------------------------------   
    !*******************************************************************************************
      integer function nearer_cell(X, NCN)  !近セルの探索（隣接セルから）
            integer, intent(in) :: NCN
            double precision, intent(in) :: X(3)
            integer NC, IIaround, index_min
            double precision :: distancecheck(2)
            double precision, allocatable :: distance(:)
            !=====================================================================================
            nearer_cell = NCN
            allocate(distance(NCMAX))
            distancecheck(1) = norm2(CENC(:,nearer_cell)-X(:))   !注目セル重心と粒子との距離
            
            check:DO
                  distance(:) = 1.0d10     !初期値はなるべく大きくとる
            
                  DO NC = 1, NUM_NC(nearer_cell)  !全隣接セルに対してループ。
                        IIaround = NEXT_CELL(NC, nearer_cell)       !現時点で近いとされるセルの隣接セルのひとつに注目
                        IF (IIaround > 0) then
                              distance(NC) = norm2(CENC(:,IIaround)-X(:))   !注目セル重心と粒子との距離を距離配列に代入
                        END IF
            
                  END DO
            
                  distancecheck(2) = minval(distance,dim=1)     !距離配列の最小値
            
                  if(distancecheck(2) < distancecheck(1)) then !より近いセルの発見で条件満足
                        distancecheck(1) = distancecheck(2)    !最小値の更新
                        index_min = minloc(distance,dim=1)            !最小値のインデックス
                        nearer_cell = NEXT_CELL(index_min, nearer_cell)    !現時点で近いとされるセルの更新
                        if(nearer_cell==0) then
                              print*,'nearer_cell_error', nearer_cell, X(:)
                              return
                        end if

                  else  !より近いセルを発見できなかった場合
            
                        exit check     !ループ脱出
            
                  end if
            
            END DO check
            
      end function nearer_cell
    !*******************************************************************************************
                      
      logical function adhesion_check(vn, NCN)

            integer JJ, JB
            integer, intent(in) :: vn, NCN
            double precision :: r_vector(3), inner

            adhesion_check = .false.

            do JJ = 1, NoB(NCN)
                  JB = ICB(JJ, NCN)

                  r_vector(:) = crd_drp(:,vn,1) - CENF(:,JB,1)

                  inner = sum(r_vector(:)*NVECF(:,JB))

                  if (inner >= 0.0d0) then
                        adhesion_check = .true. !外向き法線ベクトルと位置ベクトルの内積が正なら付着判定
                        adhes_bound(vn) = JB     !付着した境界面番号
                  end if
            end do
    
      end function adhesion_check
 
                     
      logical function nearcell_check(X, NCN)
            double precision, intent(in) :: X(3)
            integer, intent(in) :: NCN
            double precision :: distance

            distance = norm2(X(:)-CENC(:,NCN))

            if (distance < 1.0d1*WIDC(NCN)) then
                  nearcell_check = .True.
            else
                  nearcell_check = .False.
            end if
    

      end function nearcell_check


      integer function get_num_air(n_virus)
            integer, intent(in) :: n_virus
            integer Lamda

            get_num_air = int(dble(N_virus)*RDT)    !気流計算における経過ステップ数に相当

            Lamda = LoopF - LoopS
            
            if((Lamda>0).and.(get_num_air>(LoopF-OFFSET))) get_num_air = mod(get_num_air, Lamda)

            get_num_air = get_num_air + OFFSET

      end function

      subroutine read_flow_field(n_virus)
            integer, intent(in) :: n_virus
            integer FNUM
            character(99) FNAME
            character(4) digits_fmt

            FNUM = get_num_air(n_virus)

            write(digits_fmt,'("i", i1, ".", i1)') FNAME_DIGITS, FNAME_DIGITS

            select case(FILE_TYPE)
                  case('VTK')
                        if (INTERVAL_FLOW == -1) then !定常解析
                              FNAME = trim(PATH_AIR)//trim(FNAME_FMT)
                        else
                              write(FNAME,'("'//trim(PATH_AIR)//trim(HEAD_AIR)//'",'//digits_fmt//',".vtk")') FNUM

                        end if
                        call read_VTK(FNAME)

                  case('INP')
                        if(INTERVAL_FLOW==-1) then
                              FNAME = trim(PATH_AIR)//trim(FNAME_FMT)
                        else
                              if(FNUM==0) then
                                  write(FNAME,'("'//trim(PATH_AIR)//trim(HEAD_AIR)//'",'//digits_fmt//',".inp")') 1
                              else
                                  write(FNAME,'("'//trim(PATH_AIR)//trim(HEAD_AIR)//'",'//digits_fmt//',".inp")') FNUM
                              end if
                        end if
                        call read_INP(FNAME)   !INPを読み込む(SHARP用)

                  case default
                        print*,'FILE_TYPE NG:', FILE_TYPE
                        STOP
                        
            end select
                
            MAX_CDN(1) = maxval(CDN(1,:))
            MAX_CDN(2) = maxval(CDN(2,:))
            MAX_CDN(3) = maxval(CDN(3,:))
            print*, 'MAX_coordinates=', MAX_CDN(:)
                
            MIN_CDN(1) = minval(CDN(1,:))
            MIN_CDN(2) = minval(CDN(2,:))
            MIN_CDN(3) = minval(CDN(3,:))
            print*, 'MIN_coordinates=', MIN_CDN(:)
                
            call set_gravity_center
            call boundary_set
            if(n_virus > n_start) call boundary_move

            CENF(:,:,1) = CENF(:,:,2)
                
      end subroutine read_flow_field

                                 
      subroutine boundary_move !境界面の移動に合わせて付着飛沫も移動
            integer vn, JB

            ! print*, 'CALL:boundary_move'
    
            do vn = 1, size(adhesion)
            
                  if (adhesion(vn) <= 0) cycle !付着していないならスルー
            
                  JB = adhes_bound(vn)
                  if (JB > 0) then
                        crd_drp(:,vn,1) = crd_drp(:,vn,1) + CENF(:,JB,2) - CENF(:,JB,1) !面重心の移動量と同じだけ移動
                  else
                        call area_check(vn)
            
                  end if
            
            end do

            ! print*, 'FIN:boundary_move'

      end subroutine boundary_move

      subroutine area_check(vn,check)
            integer, intent(in) :: vn
            logical,optional,intent(inout) :: check
            integer L

            do L = 1, 3
        
                  if(crd_drp(L,vn,1) < MIN_CDN(L)) then
                    crd_drp(L,vn,1) = MIN_CDN(L)
                    if(present(check)) check = .true.
                  else if(crd_drp(L,vn,1) > MAX_CDN(L)) then
                    crd_drp(L,vn,1) = MAX_CDN(L)
                    if(present(check)) check = .true.
                  end if

            end do

      end subroutine area_check

      subroutine reset_status
            nearcell(:) = 0                !飛沫近接要素ID
            adhes_bound(:) = 0              !飛沫付着境界面ID
      end subroutine reset_status


end module motion_mod