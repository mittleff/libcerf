set(apps
    run_erfcx
    run_fwhm
    run_imwx
    run_voigt
    run_wofz
    tabulate
    bigloop
    )

foreach(app ${apps})
    if (${CERF_CPP})
        set_property(SOURCE ${app}.c PROPERTY LANGUAGE CXX)
    endif()
    add_executable(${app} ${app}.c)
    target_include_directories(${app} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../lib)
    if (${CERF_CPP})
        set_property(TARGET ${app} PROPERTY LINKER_LANGUAGE CXX)
    endif()
    target_link_libraries(${app} ${cerf_LIBRARY})
endforeach()
