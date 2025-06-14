/*
 * Project: VeSP
 * Copyright (c) 2025.
 * Author: Zhulev Alexander Eduardovich (zhulevae)
 */

#pragma once

#include <vector>

#include "LinCore/Types.h"


namespace VeSP
{
    using namespace LinCore::Types;

    struct Parameters
    {
        Value_T objective_length;
    };

    struct Epsilons_T
    {
        Value_T zero;
        Value_T hyperplane;
        Value_T projection;
    };

    using ConstraintIndexes_T = std::vector<size_t>;

}