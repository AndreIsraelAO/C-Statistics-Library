# Aplica warnings estritos a um target (sem vazar para quem consome a biblioteca).
function(cstats_enable_warnings target)
    if(MSVC)
        target_compile_options(${target} PRIVATE /W4)
    else()
        target_compile_options(${target} PRIVATE
            -Wall -Wextra -Wpedantic -Wshadow
            -Wconversion -Wdouble-promotion -Wstrict-prototypes
        )
    endif()
endfunction()
