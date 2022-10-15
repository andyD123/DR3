#pragma once
#include "pch.h"


#include "../Vectorisation/VecX/vec.h"
#include "../Vectorisation/VecX/target_name_space.h"
#include "../Vectorisation/VecX/instruction_traits.h"

//using namespace DRC::VecDb;

//using namespace DRC::VecF4F;
//using namespace DRC::VecD2D;
//using namespace DRC::VecD4D;

//using namespace DRC::VecD8D;

using namespace DRC::VecF16F;
//using namespace DRC::VecF8F;

using Numeric = InstructionTraits<VecXX::INS>::FloatType;
#include "dr3TestUtil.h"
