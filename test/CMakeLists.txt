#テスト登録用関数
function(ADDTEST name)
    ### テスト用プログラムのビルド設定
    add_executable(${name} ${name}.f90)

    # モジュール参照
    target_link_libraries(${name} lib)

    ### テストの登録
    add_test(NAME ${name}
        COMMAND $<TARGET_FILE:${name}>
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        )

    SETCOMPILEOPTIONS(${name})

endfunction(ADDTEST)

ADDTEST(sort_test)
ADDTEST(geometry_test)
ADDTEST(cellCenterCalc_test)
ADDTEST(sort_main_vtk_ver)