# 拉取Google Test
include(FetchContent)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG v1.14.0
)
# windows 兼容设置
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)
# 测试文件
add_executable(MyFEM_tests
        test_dense_matrix.cpp
        test_Q4.cpp
        test_Bar3D.cpp
        test_skyline_matrix.cpp
        test_solver.cpp
        test_integration_bar.cpp
        test_integration_Q4.cpp
        test_prescridedisp_bar.cpp
        test_surface_load.cpp
        test_integration_surfaceload.cpp
        test_bodyforce.cpp
        test_stress_smoothing.cpp
        test_face_mapping.cpp
        test_q4_distorted.cpp
        test_h8.cpp
        test_integration_H8.cpp
        test_expansion.cpp
        test_profile.cpp
        test_assemble_mpc.cpp
        test_mpc_beta.cpp
        test_node_dof.cpp
)
target_link_libraries(MyFEM_tests PRIVATE
        MyFEM::IO
        MyFEM::Core
        MyFEM::Model
        MyFEM::Analysis
        GTest::gtest_main
)
# 自动发现测试
include(GoogleTest)
gtest_discover_tests(MyFEM_tests)