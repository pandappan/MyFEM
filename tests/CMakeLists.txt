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
)
target_link_libraries(MyFEM_tests PRIVATE
        MyFEM::Core
        MyFEM::Model
        MyFEM::Analysis
        GTest::gtest_main
)
# 自动发现测试
include(GoogleTest)
gtest_discover_tests(MyFEM_tests)