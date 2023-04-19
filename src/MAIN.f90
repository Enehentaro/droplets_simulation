PROGRAM MAIN
    !!author: KIYOTA OGURA, Y.Ida
    !!summary:
    !!- 流れ場ファイルを読み込み、その流れ場における飛沫の運動をシミュレーション
    !!- 並列化には対応していない

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

        call RunDropletsSimulation(trim(caseName(caseID)))
        
    END DO
    
END PROGRAM MAIN