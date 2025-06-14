/*
 * Project: VeSP
 * Copyright (c) 2025.
 * Author: Zhulev Alexander Eduardovich (zhulevae)
 */

#pragma once

#include "LinCore/Functions/Generative.h"
#include "LinCore/Functions/Distance.h"

#include "VeSP/Types.h"
#include "VeSP/Projections.h"


namespace VeSP
{

    inline ConstraintIndexes_T getPointBelongHyperplanesIndexes(const Problem_T& problem, const Point_T& x, const Value_T eps)
    {
        ConstraintIndexes_T indexes;
        for(size_t i = 0; i < problem.constraints_count(); ++i)
        {
            const auto constraint = problem[i];
            if (constraint.type == ConstraintType::EQL || isPointBelongToHyperplane(constraint.a, constraint.b, x, eps))
                indexes.push_back(i);
        }
        return indexes;
    }

    inline Point_T vesp(const Problem_T& problem, const Point_T& init, const Epsilons_T& eps, const Parameters& params)
    {
        const auto n = problem.coefficents_count();
        const auto obj_vector = makeLike(problem.c, params.objective_length);
        Point_T v = maxProjectionOnPolytope(problem, init, eps.projection);
        while(true)
        {
            const auto face = getPointBelongHyperplanesIndexes(problem, v, eps.hyperplane);

            if (face.size() >= n)
                break;

            Point_T u = v + obj_vector;
            const auto w = maxProjectionOnManyfold(problem, face, u, eps.projection);
            const Vector_T d = w - v;
            v = jumpingOnPolytope(problem, v, d, eps);
        }
        return v;
    }

}