cmake_minimum_required(VERSION 3.11)
project(CPLibrary)


#Put this in your CMakeLists.txt file (change any options from OFF to ON if you want):


find_package(Boost COMPONENTS unit_test_framework)

set(CMAKE_CXX_STANDARD 23)

set(HEADER_FILES include/linear_algebra/linear_algebra.h include/nt/modular_arithmetic.h include/nt/number_theory.h
        include/parser/StatefulParser.h
        include/polynomial/ring_extension.h include/polynomial/polynomial.h include/polynomial/fft.h include/topology/optimisation.h include/string/string.h include/data_structures/data_structures.h include/algebra/binary_operation.h
        include/topology/topology.h include/data_structures/statistic_tree.h include/nt/primality.h include/formal_verification/node.h include/formal_verification/abstract_decision_tree.h
        include/formal_verification/dynamic_decision_tree.h include/formal_verification/abstract_printer.h
        include/formal_verification/decision_tree.h include/graph/graph.h include/graph/union_find.h include/graph/2sat.h include/linear_algebra/decomposition.h include/linear_algebra/special_matrices.h include/linear_algebra/special_polynomials.h include/linear_algebra/vector.h include/linear_algebra/hilbert_geometry.h)

set(SOURCE_FILES src/parser/Grammar.cpp src/parser/LRParserBuilder.cpp src/parser/LLParser.cpp src/parser/SLRParserBuilder.cpp src/parser/LALRParserBuilder.cpp
        src/parser/StatefulParser.cpp
        src/parser/LR0ParserBuilder.cpp src/formal_verification/abstract_decision_tree.cpp src/formal_verification/abstract_printer.cpp
        src/formal_verification/dynamic_decision_tree.cpp src/formal_verification/node.cpp)
add_executable(SourceCombiner source_combiner.cpp ${HEADER_FILES})
find_package(Boost COMPONENTS iostreams program_options)
target_link_libraries(SourceCombiner Boost::program_options)
add_executable(GrammarSupport example/grammar.cpp ${SOURCE_FILES} ${HEADER_FILES})
add_executable(Parser example/parser.cpp ${SOURCE_FILES} ${HEADER_FILES})

add_executable(CustomProblem local/WinterCupExpectedIterations.cpp ${HEADER_FILES} ${SOURCE_FILES} local/WinterCupGrundy.cpp)


add_executable(CPLibrary main.cpp ${HEADER_FILES} ${SOURCE_FILES})
add_executable(MLExample example/machine_learning.cpp ${HEADER_FILES})
add_executable(PowerSeriesExample example/power_series_composition.cpp ${HEADER_FILES})


add_executable(LinearAlgebraTest tests/linear_algebra/linear_algebra_approximate_test.cpp tests/linear_algebra/linear_algebra_exact_test.cpp tests/print.h
        ${HEADER_FILES})
add_executable(NumberTheoryTest tests/number_theory/modular_arithmetic_test.cpp tests/number_theory/number_theory_test.cpp tests/print.h
        ${HEADER_FILES} include/nt/primality.h)
add_executable(PolynomialTest tests/polynomial/polynomial_test.cpp tests/print.h
        ${HEADER_FILES} tests/polynomial/polynomial_test.cpp)
add_executable(DataStructuresTest tests/data_structures/segment_tree_test.cpp tests/data_structures/fenwick_tree.cpp tests/print.h
        ${HEADER_FILES} )

add_executable(MachineLearningTest tests/ml/machine_learning_test.cpp tests/print.h
        ${HEADER_FILES} )

add_executable(TopologyTest tests/topology/topology_test.cpp tests/print.h
        ${HEADER_FILES} )
add_executable(AlgebraTest tests/algebra/abstract_algebra.cpp tests/print.h
        ${HEADER_FILES} tests/algebra/binary_operation.cpp)
include_directories(include)
target_link_libraries (LinearAlgebraTest ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} )
target_link_libraries (NumberTheoryTest ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
target_link_libraries(PolynomialTest ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
target_link_libraries(DataStructuresTest ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
target_link_libraries(MachineLearningTest ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
target_link_libraries(TopologyTest ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
target_link_libraries(AlgebraTest ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})