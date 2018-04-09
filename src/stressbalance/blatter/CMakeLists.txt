add_library (blatter OBJECT
  BlatterStressBalance.cc
  blatter_test.cc
  )

if (Pism_BUILD_EXTRA_EXECS)
  add_executable(blatter_test
    Blatter_implementation.c
    BlatterStressBalance.cc
    blatter_test.cc)
  target_link_libraries (blatter_test pism)
endif()
