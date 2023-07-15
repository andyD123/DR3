add_executable(lattice lattice.cpp
  americanCrankNicholsonPricer.cpp
  americanFinitDiffPricer.cpp
  americanImplicitFiniteDiff.cpp
  americanTrinomialPricer.cpp
  americanTrinomialPricerUpAndOut.cpp
  europeanBinomialPricer.cpp
  euroTrinomial.cpp
  euroTrinomialPricerWithInit.cpp
  lattice_tools.cpp
  )

target_link_libraries(lattice PUBLIC Vectorisation)

target_include_directories(lattice PUBLIC
                           "${PROJECT_BINARY_DIR}"
                           )

        