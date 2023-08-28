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
#include "span.h"

namespace DRC
{



	namespace VecLDb
	{
		using VecxD = VecD<VecLDouble>;
		using Vecx = VecD<VecLDouble>;
		using VecXX = Vec<VecLDouble>;
		using VecVW = VecView<VecLDouble>;
		using VecBL = VecBool<VecLDouble>;
		using SpanXX = Span<VecLDouble>;
		using StrdSpanXX = StridedSpan<VecLDouble>;
	};
	// experimental
	

	namespace VecDb
	{
		using VecxD = VecD<VecDouble>;
		using Vecx = VecD<VecDouble>;
		using VecXX = Vec<VecDouble>;
		using VecVW = VecView<VecDouble>;
		using VecBL = VecBool<VecDouble>;
		using SpanXX = Span<VecDouble>;
		using StrdSpanXX = StridedSpan<VecDouble>;
	};



	namespace VecD2D
	{
		using  VecxD = VecD<Vec2d>;
		using  Vecx = VecD<Vec2d>;
		using  VecXX = Vec<Vec2d>;
		using  VecVW = VecView<Vec2d>;
		using  VecBL = VecBool<Vec2d>;
		using  SpanXX = Span<Vec2d>;
		using  StrdSpanXX = StridedSpan<Vec2d>;
	};


	namespace VecD4D
	{
		using  VecxD = VecD<Vec4d>;
		using  Vecx = VecD<Vec4d>;
		using  VecXX = Vec<Vec4d>;
		using  VecVW = VecView<Vec4d>;
		using  VecBL = VecBool<Vec4d>;
		using  SpanXX = Span<Vec4d>;
		using  StrdSpanXX = StridedSpan<Vec4d>;
	};

	namespace VecD8D
	{
		using  VecxD = VecD<Vec8d>  ;
		using  Vecx = VecD<Vec8d> ;
		using  VecXX =Vec<Vec8d> ;
		using  VecVW = VecView<Vec8d> ;
		using  VecBL = VecBool<Vec8d> ;
		using  SpanXX = Span<Vec8d> ;
		using  StrdSpanXX =StridedSpan<Vec8d> ;
	};

	namespace VecF16F
	{
		using  VecxD = VecD<Vec16f>;
		using  Vecx = VecD<Vec16f>;
		using  VecXX = Vec<Vec16f>;
		using  VecVW = VecView<Vec16f>;
		using  VecBL = VecBool<Vec16f>;
		using  SpanXX = Span<Vec16f>;
		using  StrdSpanXX = StridedSpan<Vec16f>;
	};

	namespace VecF8F
	{
		using  VecxD = VecD<Vec8f>;
		using  Vecx = VecD<Vec8f>;
		using  VecXX = Vec<Vec8f>;
		using  VecVW = VecView<Vec8f>;
		using  VecBL = VecBool<Vec8f>;
		using  SpanXX = Span<Vec8f>;
		using  StrdSpanXX = StridedSpan<Vec8f>;
	};


	namespace VecF4F
	{
		using  VecxD = VecD<Vec4f>;
		using  Vecx = VecD<Vec4f>;
		using  VecXX = Vec<Vec4f>;
		using  VecVW = VecView<Vec4f>;
		using  VecBL = VecBool<Vec4f>;
		using  SpanXX = Span<Vec4f>;
		using  StrdSpanXX = StridedSpan<Vec4f>;
	};



}// namespace DRC