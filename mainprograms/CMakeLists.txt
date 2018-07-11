
#add_executable(FemSC main.cpp ${FemSCSOURCE})

#add_executable(TestIntegrate TestIntegrate.cpp ${FemSCSOURCE})

#add_executable(TestOneDProblem TestOneDProblem.cpp ${FemSCSOURCE})

#add_executable(TestTwoDProblem TestTwoDProblem.cpp ${FemSCSOURCE})

#add_executable(TestMatrix Testmatrix.cpp ${FemSCSOURCE})

add_executable(Fem2018 main2018.cpp ${FemSOURCE2018} ${SOURCES2018})

if(USING_PZ)
 target_link_libraries(Fem2018 pz)

endif(USING_PZ)

if(USING_ARMADILLO)
 target_link_libraries(Fem2018 ${ARMADILLO_LIBRARIES})

endif(USING_ARMADILLO)

if(USING_EIGEN3)
 target_link_libraries(Fem2018 ${EIGEN_LIBRARIES})

endif(USING_EIGEN3)
