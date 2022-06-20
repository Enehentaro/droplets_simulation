add_library(src
    caseName.f90         
    dropletGenerator.f90                   
    timeKeeper.f90
    SCTfile_reader.f90    
    conditionValue.f90   
    dropletMotionSimulation.f90  
    path_operator.f90      
    unstructured_grid.f90
    adjacency_solver.f90
    filename_mod.f90             
    plot3d_operator.f90    
    vector.f90
    boxCounter.f90        
    flow_field.f90               
    simpleFile_reader.f90  
    virusDroplet.f90
    dropletEquation.f90  
    initial_translate.f90        
    terminalControler.f90  
    vtkMesh_operator.f90
)

#    main.f90         CUBE2USG.f90      boxFlowField.f90      dropletCount.f90      droplet2CSV.f90   

add_executable(main main.f90)

target_link_libraries(main src)
