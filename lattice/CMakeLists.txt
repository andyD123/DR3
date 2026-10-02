add_executable(lattice lattice.cpp
  americanCrankNicholsonPricer.cpp
  americanFinitDiffPricer.cpp
  americanImplicitFiniteDiff.cpp
  americanTrinomialPricer.cpp
  americanTrinomialPricerUpAndOut.cpp
  europeanBinomialPricer.cpp
  euroTrinomial.cpp
  euroTrinomialPricerWithInit.cpp
  chapter3_examples.cpp
  chapter3_adi.cpp
  lattice_tools.cpp
  )

target_link_libraries(lattice PUBLIC Vectorisation)

target_include_directories(lattice PUBLIC
                           "${PROJECT_BINARY_DIR}"
                           )

        