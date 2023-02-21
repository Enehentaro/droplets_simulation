#プログラム登録用関数
function(ADDPROGRAM name)
    ### テスト用プログラムのビルド設定
    add_executable(${name} ${name}.f90)

    # モジュール参照
    target_link_libraries(${name} lib)

    # コンパイルオプション付与
    SETCOMPILEOPTIONS(${name})

endfunction(ADDPROGRAM)

ADDPROGRAM(MAIN)
ADDPROGRAM(CUBE2USG)
ADDPROGRAM(boxFlowField)
ADDPROGRAM(dropletCount)
ADDPROGRAM(getFinalDroplet)
ADDPROGRAM(dropletCount_timeSeries)
ADDPROGRAM(count_aboveDesk)
ADDPROGRAM(initial_translate)
ADDPROGRAM(array2VTK)
