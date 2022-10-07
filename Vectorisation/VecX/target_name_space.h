/****************************  target_name_space.h   *******************************
* Author:        Andrew Drakeford
* Date created:  2021-04-10
* Last modified: 2021-04-10
* Version:       1.0
* Project:       DR Cubed
* Description:
*
* (c) Copyright 2019 Andrew Drakeford
* Apache License version 2.0 or later.
*****************************************************************************/
#pragma once

#include "vec.h"
#include "vec_bool.h"
#include "vec_d.h"
#include "vec_bool_d.h"
#include "vec_double.h"
#include "vec_view.h"
#include "apply_operation.h"

namespace DRC
{

	namespace VecDb
	{
		using VecxD = VecD<VecDouble>;
		using Vecx = VecD<VecDouble>;
		using VecXX = Vec<VecDouble>;
		using VecVW = VecView<VecDouble>;
		using VecBL = VecBool<VecDouble>;

	};



	namespace VecD2D
	{
		typedef VecD<Vec2d>  VecxD;
		typedef VecD<Vec2d>  Vecx;
		typedef Vec<Vec2d>  VecXX;
		//typedef VecView<Vec2d>  VecVW;
		typedef VecBool<Vec2d>	VecBL;
	};


	namespace VecD4D
	{

		typedef VecD<Vec4d>  VecxD;
		typedef VecD<Vec4d>  Vecx;
		typedef Vec<Vec4d>  VecXX;
		typedef VecView<Vec4d>  VecVW;
		typedef VecBool<Vec4d>	VecBL;

	};

	namespace VecD8D
	{
		typedef VecD<Vec8d>  VecxD;
		typedef VecD<Vec8d>  Vecx;
		typedef Vec<Vec8d>  VecXX;
		typedef VecView<Vec8d>  VecVW;
		typedef VecBool<Vec8d>	VecBL;

	};

	namespace VecF16F
	{
		typedef VecD<Vec16f>  VecxD;
		typedef VecD<Vec16f>  Vecx;
		typedef Vec<Vec16f>  VecXX;
		typedef VecView<Vec16f>  VecVW;
		typedef VecBool<Vec16f>	VecBL;
	};

	namespace VecF8F
	{
		typedef VecD<Vec8f>  VecxD;
		typedef VecD<Vec8f>  Vecx;
		typedef Vec<Vec8f>  VecXX;
		typedef VecView<Vec8f>  VecVW;
		typedef VecBool<Vec8f>	VecBL;
	};


	namespace VecF4F
	{
		typedef VecD<Vec4f>  VecxD;
		typedef VecD<Vec4f>  Vecx;
		typedef Vec<Vec4f>  VecXX;
		typedef VecView<Vec4f>  VecVW;
		typedef VecBool<Vec4f>	VecBL;
	};



}// namespace DRC