function(ADDTEST name)
    ### テスト用プログラムのビルド設定
    add_executable(${name} ${name}.f90)

    # モジュールを参照するときは，ここにtarget_link_librariesやset_target_propertiesを追加する
    target_link_libraries(${name} lib)

    ### テストの登録
    add_test(NAME ${name}
        COMMAND $<TARGET_FILE:${name}> ${CMAKE_CURRENT_BINARY_DIR}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        )

endfunction(ADDTEST)

ADDTEST(sort_main_vtk_ver)