/*
 * Project: VeSP
 * Copyright (c) 2025.
 * Author: Zhulev Alexander Eduardovich (zhulevae)
 */

#include <iostream>
#include <vector>

#include "VeSP/VeSP.h"


LinCore::Problem_T getProblem()
{
    const std::vector<LinCore::Constraint_T> constraints{
            {{1, 0, 0}, 200, LinCore::ConstraintType::LEQ},
            {{0, 1, 0}, 200, LinCore::ConstraintType::LEQ},
            {{0, 0, 1}, 200, LinCore::ConstraintType::LEQ},
            {{1, 1, 1}, 500, LinCore::ConstraintType::LEQ},
            {{-1, 0, 0}, 0, LinCore::ConstraintType::LEQ},
            {{0, -1, 0}, 0, LinCore::ConstraintType::LEQ},
            {{0, 0, -1}, 0, LinCore::ConstraintType::LEQ}
    };
    const LinCore::Vector_T c{1, 2, 3};
    return LinCore::Problem_T::from_constraints(constraints, c);
}


int main()
{
    const LinCore::Problem_T problem = getProblem();
    const LinCore::Point_T init_point{1000, 1000, 1000};

    constexpr VeSP::Epsilons_T epsilons
    {
       .zero = 1E-11,
       .hyperplane = 1E-8,
       .projection = 1E-9,
    };

    constexpr VeSP::Parameters parameters
    {
        .objective_length = 1E+7
    };

    const LinCore::Point_T vertex = VeSP::vesp(problem, init_point, epsilons, parameters);

    std::cout << vertex << std::endl;

    return 0;
}
