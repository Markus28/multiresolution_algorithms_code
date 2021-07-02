CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)

FIND_PACKAGE(deal.II 9.1.0 REQUIRED
        HINTS ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
        )
DEAL_II_INITIALIZE_CACHED_VARIABLES()

PROJECT(Biomimetics)


add_executable(main main.cpp Ellipse.cpp Ellipse.h InverseRadonTransform.cpp InverseRadonTransform.h SheppLoganAbsorption.cpp SheppLoganAbsorption.h Absorption.cpp Absorption.h BoxAbsorption.cpp BoxAbsorption.h ForwardProblem.cpp ForwardProblem.h
        PerturbedAbsorption.cpp PerturbedAbsorption.h ExponentialWave.cpp ExponentialWave.h measure.cpp measure.h InverseProblem.cpp InverseProblem.h MeshAbsorption.cpp MeshAbsorption.h MeshEvaluator.cpp MeshEvaluator.h SquareGrid.cpp SquareGrid.h FuzzyCircleAbsorption.cpp FuzzyCircleAbsorption.h)
DEAL_II_SETUP_TARGET(main)

