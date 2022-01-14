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
            
            call simulationSetUp(get_caseName())          !SetUp

            call mainDropletLoop                !mainLoop

            call output_ResultSummary       !最終結果出力
            
      END DO
    
END PROGRAM MAIN