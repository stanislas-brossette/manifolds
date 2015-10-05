// Copyright (c) 2015 CNRS
// Authors: Pierre Gergondet, Adrien Escande

// This file is part of manifolds
// manifolds is free software: you can redistribute it
// and/or modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation, either version
// 3 of the License, or (at your option) any later version.

// manifolds is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Lesser Public License for more details.  You should have
// received a copy of the GNU Lesser General Public License along with
// manifolds. If not, see
// <http://www.gnu.org/licenses/>.

#pragma once

#ifdef WIN32
#define MANIFOLDS_DLLIMPORT __declspec(dllimport)
#define MANIFOLDS_DLLEXPORT __declspec(dllexport)
#else
#define MANIFOLDS_DLLIMPORT
#define MANIFOLDS_DLLEXPORT
#endif

#ifdef MANIFOLDS_EXPORT
#define MANIFOLDS_API MANIFOLDS_DLLEXPORT
#else
#define MANIFOLDS_API MANIFOLDS_DLLIMPORT
#endif

