!---------------------------------------------------------------------------------
!     Simulation of viral droplets
!                                by KIYOTA OGURA(2021/1/10)
!     updated by IDA
!---------------------------------------------------------------------------------
PROGRAM MAIN
      !$ use omp_lib
      use dropletMotionSimulation
      use caseName_m
      implicit none
      integer, pointer :: nc => nowCase
      integer nc_max

      !$OMP parallel
            !$OMP single
            !$ print *, "Num threads:", omp_get_num_threads()
            !$OMP end single
      !$OMP end parallel

      call case_check(num_case=nc_max) 

      DO nc = 1, nc_max                         !実行数だけループ（通常1回）

            call create_CaseDirectory                     !ディレクトリ作成
            
            call firstSet_mainDroplet          !条件TXTの読み込み、飛沫初期分布の計算など

            call dropletManagement       !外部サブルーチンによる管理

            call checkpoint                    !計算条件の確認および時刻計測のためのチェックポイント

            call mainDropletLoop

            call output_ResultSummary       !最終結果出力
            
      END DO


      contains

      subroutine mainDropletLoop
            integer, pointer :: n => timeStep
            
            print '("*******************************************")'
            print '("            START step_loop                ")'
            print '("*******************************************")'
    
            do n = n_start + 1, n_end           !ステップ数だけループ
    
                call mainDroplet_process      !飛沫計算の一連の処理
    
                call dropletManagement       !外部サブルーチンによる管理
    
                if (mod(n, outputInterval) == 0) call periodicOutput             !出力
    
                call check_FlowFieldUpdate        !流れ場の更新チェック
    
            end do
    
            print '("*******************************************")'
            print '("             END step_loop                 ")'
            print '("*******************************************")'
    
      end subroutine
    
END PROGRAM MAIN