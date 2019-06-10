################################################################################
include(DownloadProject)

# With CMake 3.8 and above, we can hide warnings about git being in a
# detached head by passing an extra GIT_CONFIG option
if(NOT (${CMAKE_VERSION} VERSION_LESS "3.8.0"))
    set(TRIWILD_EXTRA_OPTIONS "GIT_CONFIG advice.detachedHead=false")
else()
    set(TRIWILD_EXTRA_OPTIONS "")
endif()

# Shortcut function
function(triwild_download_project name)
    download_project(
        PROJ         ${name}
        SOURCE_DIR   ${THIRD_PARTY_DIR}/${name}
        DOWNLOAD_DIR ${THIRD_PARTY_DIR}/.cache/${name}
        QUIET
        ${TRIWILD_EXTRA_OPTIONS}
        ${ARGN}
    )
endfunction()

################################################################################

## aabbcc
function(triwild_download_aabbcc)
    triwild_download_project(aabbcc
        GIT_REPOSITORY https://github.com/Yixin-Hu/aabbcc
        GIT_TAG        91838aff841627472e78328a79e6800c7bc1678b
    )
endfunction()


## geogram
function(triwild_download_geogram)
    triwild_download_project(geogram
        GIT_REPOSITORY https://github.com/alicevision/geogram
        GIT_TAG        0ac4a5889f8eaef9372b888126deec2334128158
    )
endfunction()


## libigl
function(triwild_download_igl)
    triwild_download_project(libigl
        GIT_REPOSITORY https://github.com/libigl/libigl
        GIT_TAG        608fc010a5d65f2edede2e7b64cf09e248d76e15
    )
endfunction()


## nlopt
function(triwild_download_nlopt)
    triwild_download_project(nlopt
        GIT_REPOSITORY https://github.com/stevengj/nlopt
        GIT_TAG        37b74a8c2037eea5dc72fea7eeb9b850fa978913
    )
endfunction()





## GMP for windows
function(triwild_download_gmp_cygwin)
    triwild_download_project(gmp
        URL https://cs.nyu.edu/~exact/core/gmp/gmp-static-cygwin-4.1.tar.gz
        URL_MD5        3c1454c7bf57f71d7c90190ce133bdbe
    )
endfunction()


## GMP for windows
function(triwild_download_gmp_mingw)
    triwild_download_project(gmp
        URL https://cs.nyu.edu/~exact/core/gmp/gmp-static-mingw-4.1.tar.gz
        URL_MD5        bef87e429338a4568a6c1321fcecadcd
    )
endfunction()


## GMP for windows
function(triwild_download_gmp_vc)
    triwild_download_project(gmp
        URL https://cs.nyu.edu/~exact/core/gmp/gmp-static-vc-4.1.2.zip
        URL_MD5        06e9efda0da0259de47c5bc0f8be81e4
    )
endfunction()



