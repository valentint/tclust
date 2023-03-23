/*
    SMat - Simple Matrix Classes v0.1beta
    Copyright (C) 2011 by Heinrich Fritz (heinrich_fritz@hotmail.com)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef SMAT_H
#define SMAT_H

#ifndef SMAT_FLAG_NO_INI
//#include "smat.ini.h"	//	ON ERROR / IF NOT FOUND: create an empty "smat.ini.h" file in your source directory.
						//	smat.ini.h is a user-defined header file for providing smat with additional declarations.
#endif	//	#ifndef SMAT_FLAG_NO_INI

#include "smat.def.h"
#include "smat.base.h"
#include "smat.elop.h"
#include "smat.math.h"
#include "smat.mem.h"
#include "smat.matop.h"
#include "smat.misc.h"
#include "smat.random.h"
#include "smat.sort.h"
#include "smat.stat.h"
#include "smat.meal.h"

#endif	//	#ifndef SMAT_H
