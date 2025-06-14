/*
 * Project: VeSP
 * Copyright (c) 2025.
 * Author: Zhulev Alexander Eduardovich (zhulevae)
 */

#pragma once

#include <limits>
#include <stdexcept>

#include "LinCore/Functions/Arithmetic.h"
#include "LinCore/Functions/Distance.h"
#include "LinCore/Functions/Projection.h"

#include "VeSP/Types.h"


namespace VeSP
{

    inline Point_T jumpingOnPolytope(const Problem_T& problem, const Point_T& x, const Vector_T& d, const Epsilons_T& eps)
    {
        Point_T min_projection(x.size(), 0);
        Value_T min_length = std::numeric_limits<Value_T>::max();

        for (size_t i = 0; i < problem.constraints_count(); ++i)
        {
            if (problem.types[i] == ConstraintType::EQL)
                continue;

            const auto a = problem.A[i];
            const auto b = problem.b[i];

            if (isPerpendicularVector(a, d, eps.zero))
                continue;

            if(isPointBelongToHyperplane(a, b, x, eps.hyperplane))
                return x;

            if(isPointBelongToHalfspace(a, b, x, eps.hyperplane))
            {
                const auto projection = obliqueProjection(a, b, d, x);
                const auto length = euclideanNorm(projection);
                if (min_length > length)
                {
                    min_projection = projection;
                    min_length = length;
                }
                continue;
            }

            throw std::runtime_error("Point outside of polytope");
        }

        return x + min_projection;
    }

    inline Point_T maxProjectionOnPolytope(const Problem_T& problem, const Point_T& x, const Value_T& eps)
    {
        Point_T w = x;
        Point_T max_w = w;
        Value_T max_dist_prev = 0;
        Value_T delta_dist;

        do
        {
            Value_T max_dist = 0;
            for (auto i = 0; i < problem.constraints_count(); ++i)
            {
                const auto constraint = problem[i];

                if (constraint.type != ConstraintType::EQL && isPointBelongToHalfspace(constraint.a, constraint.b, w, eps))
                    continue;

                Point_T projection = orthogonalProjection(constraint.a, constraint.b, w);
                const Value_T dist = euclideanNorm(projection);
                if (dist > max_dist)
                {
                    max_w = w + projection;
                    max_dist = dist;
                }
            }

            if (max_dist == 0)
                break;

            w = max_w;
            delta_dist = fabs(max_dist - max_dist_prev);
            max_dist_prev = max_dist;

        } while (delta_dist >= eps);

        return w;
    }

    inline Point_T maxProjectionOnManyfold(const Problem_T& problem, const ConstraintIndexes_T& manifold, const Point_T& x, const Value_T& eps)
    {
        Point_T w = x;
        Point_T max_w = w;
        Value_T max_dist_prev = 0;
        Value_T delta_dist;

        do
        {
            Value_T max_dist = 0;
            for (const size_t i : manifold)
            {
                const auto constraint = problem[i];
                Point_T projection = orthogonalProjection(constraint.a, constraint.b, w);
                const Value_T dist = euclideanNorm(projection);
                if (dist > max_dist)
                {
                    max_w = w + projection;
                    max_dist = dist;
                }
            }

            if (max_dist == 0)
                break;

            w = max_w;
            delta_dist = fabs(max_dist - max_dist_prev);
            max_dist_prev = max_dist;

        } while (delta_dist >= eps);

        return w;
    }
}