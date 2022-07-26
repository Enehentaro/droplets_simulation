!---------------------------------------------------------------------------------
!> Simulation of Virus-Laden Droplets Behavior in AFDET
!> by KIYOTA OGURA(2021/1/10)
!> updated by YUTA IDA
!---------------------------------------------------------------------------------
PROGRAM MAIN
      !$ use omp_lib
      use dropletMotionSimulation
      use caseName_m
      implicit none
      character(50), allocatable :: caseName(:)
      integer caseID

      !$OMP parallel
            !$OMP single
            !$ print *, "Num threads:", omp_get_num_threads()
            !$OMP end single
      !$OMP end parallel

      call read_basicSettingOnSimulation

      call case_check(caseName) 

      DO caseID = 1, size(caseName)                        !実行数だけループ（通常1回）
            
            call simulationSetUp(trim(caseName(caseID)))          !SetUp

            call mainDropletLoop                !mainLoop

            call output_ResultSummary       !最終結果出力
            
      END DO
    
END PROGRAM MAIN