file(GLOB tests *test.c)

foreach(src ${tests})
    if (${CERF_CPP})
        set_property(SOURCE ${src} PROPERTY LANGUAGE CXX)
    endif()
    get_filename_component(test ${src} NAME_WE)
    add_executable(${test} ${src})
    target_include_directories(${test} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../lib)
    if (${CERF_CPP})
        set_property(TARGET ${test} PROPERTY LINKER_LANGUAGE CXX)
    endif()
    target_link_libraries(${test} ${cerf_LIBRARY})
    add_test(NAME ${test} COMMAND ${test})
endforeach()
