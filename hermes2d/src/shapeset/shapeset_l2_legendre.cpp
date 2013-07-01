// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#include "global.h"
#include "shapeset.h"
#include "shapeset_common.h"
#include "shapeset_l2_all.h"

namespace Hermes
{
  namespace Hermes2D
  {
    static double leg_quad_l0_l0(double x, double y)
    {
      return Legendre0(x) * Legendre0(y);
    }

    static double leg_quad_l0_l0x(double x, double y)
    {
      return Legendre0x(x) * Legendre0(y);
    }

    static double leg_quad_l0_l0y(double x, double y)
    {
      return Legendre0(x) * Legendre0x(y);
    }

    static double leg_quad_l0_l0xx(double x, double y)
    {
      return Legendre0xx(x) * Legendre0(y);
    }

    static double leg_quad_l0_l0xy(double x, double y)
    {
      return Legendre0x(x) * Legendre0x(y);
    }

    static double leg_quad_l0_l0yy(double x, double y)
    {
      return  Legendre0(x) * Legendre0xx(y);
    }

    static double leg_quad_l0_l1(double x, double y)
    {
      return Legendre0(x) * Legendre1(y);
    }

    static double leg_quad_l0_l1x(double x, double y)
    {
      return Legendre0x(x) * Legendre1(y);
    }

    static double leg_quad_l0_l1y(double x, double y)
    {
      return Legendre0(x) * Legendre1x(y);
    }

    static double leg_quad_l0_l1xx(double x, double y)
    {
      return Legendre0xx(x) * Legendre1(y);
    }

    static double leg_quad_l0_l1xy(double x, double y)
    {
      return Legendre0x(x) * Legendre1x(y);
    }

    static double leg_quad_l0_l1yy(double x, double y)
    {
      return  Legendre0(x) * Legendre1xx(y);
    }

    static double leg_quad_l0_l2(double x, double y)
    {
      return Legendre0(x) * Legendre2(y);
    }

    static double leg_quad_l0_l2x(double x, double y)
    {
      return Legendre0x(x) * Legendre2(y);
    }

    static double leg_quad_l0_l2y(double x, double y)
    {
      return Legendre0(x) * Legendre2x(y);
    }

    static double leg_quad_l0_l2xx(double x, double y)
    {
      return Legendre0xx(x) * Legendre2(y);
    }

    static double leg_quad_l0_l2xy(double x, double y)
    {
      return Legendre0x(x) * Legendre2x(y);
    }

    static double leg_quad_l0_l2yy(double x, double y)
    {
      return  Legendre0(x) * Legendre2xx(y);
    }

    static double leg_quad_l0_l3(double x, double y)
    {
      return Legendre0(x) * Legendre3(y);
    }

    static double leg_quad_l0_l3x(double x, double y)
    {
      return Legendre0x(x) * Legendre3(y);
    }

    static double leg_quad_l0_l3y(double x, double y)
    {
      return Legendre0(x) * Legendre3x(y);
    }

    static double leg_quad_l0_l3xx(double x, double y)
    {
      return Legendre0xx(x) * Legendre3(y);
    }

    static double leg_quad_l0_l3xy(double x, double y)
    {
      return Legendre0x(x) * Legendre3x(y);
    }

    static double leg_quad_l0_l3yy(double x, double y)
    {
      return  Legendre0(x) * Legendre3xx(y);
    }

    static double leg_quad_l0_l4(double x, double y)
    {
      return Legendre0(x) * Legendre4(y);
    }

    static double leg_quad_l0_l4x(double x, double y)
    {
      return Legendre0x(x) * Legendre4(y);
    }

    static double leg_quad_l0_l4y(double x, double y)
    {
      return Legendre0(x) * Legendre4x(y);
    }

    static double leg_quad_l0_l4xx(double x, double y)
    {
      return Legendre0xx(x) * Legendre4(y);
    }

    static double leg_quad_l0_l4xy(double x, double y)
    {
      return Legendre0x(x) * Legendre4x(y);
    }

    static double leg_quad_l0_l4yy(double x, double y)
    {
      return  Legendre0(x) * Legendre4xx(y);
    }

    static double leg_quad_l0_l5(double x, double y)
    {
      return Legendre0(x) * Legendre5(y);
    }

    static double leg_quad_l0_l5x(double x, double y)
    {
      return Legendre0x(x) * Legendre5(y);
    }

    static double leg_quad_l0_l5y(double x, double y)
    {
      return Legendre0(x) * Legendre5x(y);
    }

    static double leg_quad_l0_l5xx(double x, double y)
    {
      return Legendre0xx(x) * Legendre5(y);
    }

    static double leg_quad_l0_l5xy(double x, double y)
    {
      return Legendre0x(x) * Legendre5x(y);
    }

    static double leg_quad_l0_l5yy(double x, double y)
    {
      return  Legendre0(x) * Legendre5xx(y);
    }

    static double leg_quad_l0_l6(double x, double y)
    {
      return Legendre0(x) * Legendre6(y);
    }

    static double leg_quad_l0_l6x(double x, double y)
    {
      return Legendre0x(x) * Legendre6(y);
    }

    static double leg_quad_l0_l6y(double x, double y)
    {
      return Legendre0(x) * Legendre6x(y);
    }

    static double leg_quad_l0_l6xx(double x, double y)
    {
      return Legendre0xx(x) * Legendre6(y);
    }

    static double leg_quad_l0_l6xy(double x, double y)
    {
      return Legendre0x(x) * Legendre6x(y);
    }

    static double leg_quad_l0_l6yy(double x, double y)
    {
      return  Legendre0(x) * Legendre6xx(y);
    }

    static double leg_quad_l0_l7(double x, double y)
    {
      return Legendre0(x) * Legendre7(y);
    }

    static double leg_quad_l0_l7x(double x, double y)
    {
      return Legendre0x(x) * Legendre7(y);
    }

    static double leg_quad_l0_l7y(double x, double y)
    {
      return Legendre0(x) * Legendre7x(y);
    }

    static double leg_quad_l0_l7xx(double x, double y)
    {
      return Legendre0xx(x) * Legendre7(y);
    }

    static double leg_quad_l0_l7xy(double x, double y)
    {
      return Legendre0x(x) * Legendre7x(y);
    }

    static double leg_quad_l0_l7yy(double x, double y)
    {
      return  Legendre0(x) * Legendre7xx(y);
    }

    static double leg_quad_l0_l8(double x, double y)
    {
      return Legendre0(x) * Legendre8(y);
    }

    static double leg_quad_l0_l8x(double x, double y)
    {
      return Legendre0x(x) * Legendre8(y);
    }

    static double leg_quad_l0_l8y(double x, double y)
    {
      return Legendre0(x) * Legendre8x(y);
    }

    static double leg_quad_l0_l8xx(double x, double y)
    {
      return Legendre0xx(x) * Legendre8(y);
    }

    static double leg_quad_l0_l8xy(double x, double y)
    {
      return Legendre0x(x) * Legendre8x(y);
    }

    static double leg_quad_l0_l8yy(double x, double y)
    {
      return  Legendre0(x) * Legendre8xx(y);
    }

    static double leg_quad_l0_l9(double x, double y)
    {
      return Legendre0(x) * Legendre9(y);
    }

    static double leg_quad_l0_l9x(double x, double y)
    {
      return Legendre0x(x) * Legendre9(y);
    }

    static double leg_quad_l0_l9y(double x, double y)
    {
      return Legendre0(x) * Legendre9x(y);
    }

    static double leg_quad_l0_l9xx(double x, double y)
    {
      return Legendre0xx(x) * Legendre9(y);
    }

    static double leg_quad_l0_l9xy(double x, double y)
    {
      return Legendre0x(x) * Legendre9x(y);
    }

    static double leg_quad_l0_l9yy(double x, double y)
    {
      return  Legendre0(x) * Legendre9xx(y);
    }

    static double leg_quad_l0_l10(double x, double y)
    {
      return Legendre0(x) * Legendre10(y);
    }

    static double leg_quad_l0_l10x(double x, double y)
    {
      return Legendre0x(x) * Legendre10(y);
    }

    static double leg_quad_l0_l10y(double x, double y)
    {
      return Legendre0(x) * Legendre10x(y);
    }

    static double leg_quad_l0_l10xx(double x, double y)
    {
      return Legendre0xx(x) * Legendre10(y);
    }

    static double leg_quad_l0_l10xy(double x, double y)
    {
      return Legendre0x(x) * Legendre10x(y);
    }

    static double leg_quad_l0_l10yy(double x, double y)
    {
      return  Legendre0(x) * Legendre10xx(y);
    }

    static double leg_quad_l1_l0(double x, double y)
    {
      return Legendre1(x) * Legendre0(y);
    }

    static double leg_quad_l1_l0x(double x, double y)
    {
      return Legendre1x(x) * Legendre0(y);
    }

    static double leg_quad_l1_l0y(double x, double y)
    {
      return Legendre1(x) * Legendre0x(y);
    }

    static double leg_quad_l1_l0xx(double x, double y)
    {
      return Legendre1xx(x) * Legendre0(y);
    }

    static double leg_quad_l1_l0xy(double x, double y)
    {
      return Legendre1x(x) * Legendre0x(y);
    }

    static double leg_quad_l1_l0yy(double x, double y)
    {
      return  Legendre1(x) * Legendre0xx(y);
    }

    static double leg_quad_l1_l1(double x, double y)
    {
      return Legendre1(x) * Legendre1(y);
    }

    static double leg_quad_l1_l1x(double x, double y)
    {
      return Legendre1x(x) * Legendre1(y);
    }

    static double leg_quad_l1_l1y(double x, double y)
    {
      return Legendre1(x) * Legendre1x(y);
    }

    static double leg_quad_l1_l1xx(double x, double y)
    {
      return Legendre1xx(x) * Legendre1(y);
    }

    static double leg_quad_l1_l1xy(double x, double y)
    {
      return Legendre1x(x) * Legendre1x(y);
    }

    static double leg_quad_l1_l1yy(double x, double y)
    {
      return  Legendre1(x) * Legendre1xx(y);
    }

    static double leg_quad_l1_l2(double x, double y)
    {
      return Legendre1(x) * Legendre2(y);
    }

    static double leg_quad_l1_l2x(double x, double y)
    {
      return Legendre1x(x) * Legendre2(y);
    }

    static double leg_quad_l1_l2y(double x, double y)
    {
      return Legendre1(x) * Legendre2x(y);
    }

    static double leg_quad_l1_l2xx(double x, double y)
    {
      return Legendre1xx(x) * Legendre2(y);
    }

    static double leg_quad_l1_l2xy(double x, double y)
    {
      return Legendre1x(x) * Legendre2x(y);
    }

    static double leg_quad_l1_l2yy(double x, double y)
    {
      return  Legendre1(x) * Legendre2xx(y);
    }

    static double leg_quad_l1_l3(double x, double y)
    {
      return Legendre1(x) * Legendre3(y);
    }

    static double leg_quad_l1_l3x(double x, double y)
    {
      return Legendre1x(x) * Legendre3(y);
    }

    static double leg_quad_l1_l3y(double x, double y)
    {
      return Legendre1(x) * Legendre3x(y);
    }

    static double leg_quad_l1_l3xx(double x, double y)
    {
      return Legendre1xx(x) * Legendre3(y);
    }

    static double leg_quad_l1_l3xy(double x, double y)
    {
      return Legendre1x(x) * Legendre3x(y);
    }

    static double leg_quad_l1_l3yy(double x, double y)
    {
      return  Legendre1(x) * Legendre3xx(y);
    }

    static double leg_quad_l1_l4(double x, double y)
    {
      return Legendre1(x) * Legendre4(y);
    }

    static double leg_quad_l1_l4x(double x, double y)
    {
      return Legendre1x(x) * Legendre4(y);
    }

    static double leg_quad_l1_l4y(double x, double y)
    {
      return Legendre1(x) * Legendre4x(y);
    }

    static double leg_quad_l1_l4xx(double x, double y)
    {
      return Legendre1xx(x) * Legendre4(y);
    }

    static double leg_quad_l1_l4xy(double x, double y)
    {
      return Legendre1x(x) * Legendre4x(y);
    }

    static double leg_quad_l1_l4yy(double x, double y)
    {
      return  Legendre1(x) * Legendre4xx(y);
    }

    static double leg_quad_l1_l5(double x, double y)
    {
      return Legendre1(x) * Legendre5(y);
    }

    static double leg_quad_l1_l5x(double x, double y)
    {
      return Legendre1x(x) * Legendre5(y);
    }

    static double leg_quad_l1_l5y(double x, double y)
    {
      return Legendre1(x) * Legendre5x(y);
    }

    static double leg_quad_l1_l5xx(double x, double y)
    {
      return Legendre1xx(x) * Legendre5(y);
    }

    static double leg_quad_l1_l5xy(double x, double y)
    {
      return Legendre1x(x) * Legendre5x(y);
    }

    static double leg_quad_l1_l5yy(double x, double y)
    {
      return  Legendre1(x) * Legendre5xx(y);
    }

    static double leg_quad_l1_l6(double x, double y)
    {
      return Legendre1(x) * Legendre6(y);
    }

    static double leg_quad_l1_l6x(double x, double y)
    {
      return Legendre1x(x) * Legendre6(y);
    }

    static double leg_quad_l1_l6y(double x, double y)
    {
      return Legendre1(x) * Legendre6x(y);
    }

    static double leg_quad_l1_l6xx(double x, double y)
    {
      return Legendre1xx(x) * Legendre6(y);
    }

    static double leg_quad_l1_l6xy(double x, double y)
    {
      return Legendre1x(x) * Legendre6x(y);
    }

    static double leg_quad_l1_l6yy(double x, double y)
    {
      return  Legendre1(x) * Legendre6xx(y);
    }

    static double leg_quad_l1_l7(double x, double y)
    {
      return Legendre1(x) * Legendre7(y);
    }

    static double leg_quad_l1_l7x(double x, double y)
    {
      return Legendre1x(x) * Legendre7(y);
    }

    static double leg_quad_l1_l7y(double x, double y)
    {
      return Legendre1(x) * Legendre7x(y);
    }

    static double leg_quad_l1_l7xx(double x, double y)
    {
      return Legendre1xx(x) * Legendre7(y);
    }

    static double leg_quad_l1_l7xy(double x, double y)
    {
      return Legendre1x(x) * Legendre7x(y);
    }

    static double leg_quad_l1_l7yy(double x, double y)
    {
      return  Legendre1(x) * Legendre7xx(y);
    }

    static double leg_quad_l1_l8(double x, double y)
    {
      return Legendre1(x) * Legendre8(y);
    }

    static double leg_quad_l1_l8x(double x, double y)
    {
      return Legendre1x(x) * Legendre8(y);
    }

    static double leg_quad_l1_l8y(double x, double y)
    {
      return Legendre1(x) * Legendre8x(y);
    }

    static double leg_quad_l1_l8xx(double x, double y)
    {
      return Legendre1xx(x) * Legendre8(y);
    }

    static double leg_quad_l1_l8xy(double x, double y)
    {
      return Legendre1x(x) * Legendre8x(y);
    }

    static double leg_quad_l1_l8yy(double x, double y)
    {
      return  Legendre1(x) * Legendre8xx(y);
    }

    static double leg_quad_l1_l9(double x, double y)
    {
      return Legendre1(x) * Legendre9(y);
    }

    static double leg_quad_l1_l9x(double x, double y)
    {
      return Legendre1x(x) * Legendre9(y);
    }

    static double leg_quad_l1_l9y(double x, double y)
    {
      return Legendre1(x) * Legendre9x(y);
    }

    static double leg_quad_l1_l9xx(double x, double y)
    {
      return Legendre1xx(x) * Legendre9(y);
    }

    static double leg_quad_l1_l9xy(double x, double y)
    {
      return Legendre1x(x) * Legendre9x(y);
    }

    static double leg_quad_l1_l9yy(double x, double y)
    {
      return  Legendre1(x) * Legendre9xx(y);
    }

    static double leg_quad_l1_l10(double x, double y)
    {
      return Legendre1(x) * Legendre10(y);
    }

    static double leg_quad_l1_l10x(double x, double y)
    {
      return Legendre1x(x) * Legendre10(y);
    }

    static double leg_quad_l1_l10y(double x, double y)
    {
      return Legendre1(x) * Legendre10x(y);
    }

    static double leg_quad_l1_l10xx(double x, double y)
    {
      return Legendre1xx(x) * Legendre10(y);
    }

    static double leg_quad_l1_l10xy(double x, double y)
    {
      return Legendre1x(x) * Legendre10x(y);
    }

    static double leg_quad_l1_l10yy(double x, double y)
    {
      return  Legendre1(x) * Legendre10xx(y);
    }

    static double leg_quad_l2_l0(double x, double y)
    {
      return Legendre2(x) * Legendre0(y);
    }

    static double leg_quad_l2_l0x(double x, double y)
    {
      return Legendre2x(x) * Legendre0(y);
    }

    static double leg_quad_l2_l0y(double x, double y)
    {
      return Legendre2(x) * Legendre0x(y);
    }

    static double leg_quad_l2_l0xx(double x, double y)
    {
      return Legendre2xx(x) * Legendre0(y);
    }

    static double leg_quad_l2_l0xy(double x, double y)
    {
      return Legendre2x(x) * Legendre0x(y);
    }

    static double leg_quad_l2_l0yy(double x, double y)
    {
      return  Legendre2(x) * Legendre0xx(y);
    }

    static double leg_quad_l2_l1(double x, double y)
    {
      return Legendre2(x) * Legendre1(y);
    }

    static double leg_quad_l2_l1x(double x, double y)
    {
      return Legendre2x(x) * Legendre1(y);
    }

    static double leg_quad_l2_l1y(double x, double y)
    {
      return Legendre2(x) * Legendre1x(y);
    }

    static double leg_quad_l2_l1xx(double x, double y)
    {
      return Legendre2xx(x) * Legendre1(y);
    }

    static double leg_quad_l2_l1xy(double x, double y)
    {
      return Legendre2x(x) * Legendre1x(y);
    }

    static double leg_quad_l2_l1yy(double x, double y)
    {
      return  Legendre2(x) * Legendre1xx(y);
    }

    static double leg_quad_l2_l2(double x, double y)
    {
      return Legendre2(x) * Legendre2(y);
    }

    static double leg_quad_l2_l2x(double x, double y)
    {
      return Legendre2x(x) * Legendre2(y);
    }

    static double leg_quad_l2_l2y(double x, double y)
    {
      return Legendre2(x) * Legendre2x(y);
    }

    static double leg_quad_l2_l2xx(double x, double y)
    {
      return Legendre2xx(x) * Legendre2(y);
    }

    static double leg_quad_l2_l2xy(double x, double y)
    {
      return Legendre2x(x) * Legendre2x(y);
    }

    static double leg_quad_l2_l2yy(double x, double y)
    {
      return  Legendre2(x) * Legendre2xx(y);
    }

    static double leg_quad_l2_l3(double x, double y)
    {
      return Legendre2(x) * Legendre3(y);
    }

    static double leg_quad_l2_l3x(double x, double y)
    {
      return Legendre2x(x) * Legendre3(y);
    }

    static double leg_quad_l2_l3y(double x, double y)
    {
      return Legendre2(x) * Legendre3x(y);
    }

    static double leg_quad_l2_l3xx(double x, double y)
    {
      return Legendre2xx(x) * Legendre3(y);
    }

    static double leg_quad_l2_l3xy(double x, double y)
    {
      return Legendre2x(x) * Legendre3x(y);
    }

    static double leg_quad_l2_l3yy(double x, double y)
    {
      return  Legendre2(x) * Legendre3xx(y);
    }

    static double leg_quad_l2_l4(double x, double y)
    {
      return Legendre2(x) * Legendre4(y);
    }

    static double leg_quad_l2_l4x(double x, double y)
    {
      return Legendre2x(x) * Legendre4(y);
    }

    static double leg_quad_l2_l4y(double x, double y)
    {
      return Legendre2(x) * Legendre4x(y);
    }

    static double leg_quad_l2_l4xx(double x, double y)
    {
      return Legendre2xx(x) * Legendre4(y);
    }

    static double leg_quad_l2_l4xy(double x, double y)
    {
      return Legendre2x(x) * Legendre4x(y);
    }

    static double leg_quad_l2_l4yy(double x, double y)
    {
      return  Legendre2(x) * Legendre4xx(y);
    }

    static double leg_quad_l2_l5(double x, double y)
    {
      return Legendre2(x) * Legendre5(y);
    }

    static double leg_quad_l2_l5x(double x, double y)
    {
      return Legendre2x(x) * Legendre5(y);
    }

    static double leg_quad_l2_l5y(double x, double y)
    {
      return Legendre2(x) * Legendre5x(y);
    }

    static double leg_quad_l2_l5xx(double x, double y)
    {
      return Legendre2xx(x) * Legendre5(y);
    }

    static double leg_quad_l2_l5xy(double x, double y)
    {
      return Legendre2x(x) * Legendre5x(y);
    }

    static double leg_quad_l2_l5yy(double x, double y)
    {
      return  Legendre2(x) * Legendre5xx(y);
    }

    static double leg_quad_l2_l6(double x, double y)
    {
      return Legendre2(x) * Legendre6(y);
    }

    static double leg_quad_l2_l6x(double x, double y)
    {
      return Legendre2x(x) * Legendre6(y);
    }

    static double leg_quad_l2_l6y(double x, double y)
    {
      return Legendre2(x) * Legendre6x(y);
    }

    static double leg_quad_l2_l6xx(double x, double y)
    {
      return Legendre2xx(x) * Legendre6(y);
    }

    static double leg_quad_l2_l6xy(double x, double y)
    {
      return Legendre2x(x) * Legendre6x(y);
    }

    static double leg_quad_l2_l6yy(double x, double y)
    {
      return  Legendre2(x) * Legendre6xx(y);
    }

    static double leg_quad_l2_l7(double x, double y)
    {
      return Legendre2(x) * Legendre7(y);
    }

    static double leg_quad_l2_l7x(double x, double y)
    {
      return Legendre2x(x) * Legendre7(y);
    }

    static double leg_quad_l2_l7y(double x, double y)
    {
      return Legendre2(x) * Legendre7x(y);
    }

    static double leg_quad_l2_l7xx(double x, double y)
    {
      return Legendre2xx(x) * Legendre7(y);
    }

    static double leg_quad_l2_l7xy(double x, double y)
    {
      return Legendre2x(x) * Legendre7x(y);
    }

    static double leg_quad_l2_l7yy(double x, double y)
    {
      return  Legendre2(x) * Legendre7xx(y);
    }

    static double leg_quad_l2_l8(double x, double y)
    {
      return Legendre2(x) * Legendre8(y);
    }

    static double leg_quad_l2_l8x(double x, double y)
    {
      return Legendre2x(x) * Legendre8(y);
    }

    static double leg_quad_l2_l8y(double x, double y)
    {
      return Legendre2(x) * Legendre8x(y);
    }

    static double leg_quad_l2_l8xx(double x, double y)
    {
      return Legendre2xx(x) * Legendre8(y);
    }

    static double leg_quad_l2_l8xy(double x, double y)
    {
      return Legendre2x(x) * Legendre8x(y);
    }

    static double leg_quad_l2_l8yy(double x, double y)
    {
      return  Legendre2(x) * Legendre8xx(y);
    }

    static double leg_quad_l2_l9(double x, double y)
    {
      return Legendre2(x) * Legendre9(y);
    }

    static double leg_quad_l2_l9x(double x, double y)
    {
      return Legendre2x(x) * Legendre9(y);
    }

    static double leg_quad_l2_l9y(double x, double y)
    {
      return Legendre2(x) * Legendre9x(y);
    }

    static double leg_quad_l2_l9xx(double x, double y)
    {
      return Legendre2xx(x) * Legendre9(y);
    }

    static double leg_quad_l2_l9xy(double x, double y)
    {
      return Legendre2x(x) * Legendre9x(y);
    }

    static double leg_quad_l2_l9yy(double x, double y)
    {
      return  Legendre2(x) * Legendre9xx(y);
    }

    static double leg_quad_l2_l10(double x, double y)
    {
      return Legendre2(x) * Legendre10(y);
    }

    static double leg_quad_l2_l10x(double x, double y)
    {
      return Legendre2x(x) * Legendre10(y);
    }

    static double leg_quad_l2_l10y(double x, double y)
    {
      return Legendre2(x) * Legendre10x(y);
    }

    static double leg_quad_l2_l10xx(double x, double y)
    {
      return Legendre2xx(x) * Legendre10(y);
    }

    static double leg_quad_l2_l10xy(double x, double y)
    {
      return Legendre2x(x) * Legendre10x(y);
    }

    static double leg_quad_l2_l10yy(double x, double y)
    {
      return  Legendre2(x) * Legendre10xx(y);
    }

    static double leg_quad_l3_l0(double x, double y)
    {
      return Legendre3(x) * Legendre0(y);
    }

    static double leg_quad_l3_l0x(double x, double y)
    {
      return Legendre3x(x) * Legendre0(y);
    }

    static double leg_quad_l3_l0y(double x, double y)
    {
      return Legendre3(x) * Legendre0x(y);
    }

    static double leg_quad_l3_l0xx(double x, double y)
    {
      return Legendre3xx(x) * Legendre0(y);
    }

    static double leg_quad_l3_l0xy(double x, double y)
    {
      return Legendre3x(x) * Legendre0x(y);
    }

    static double leg_quad_l3_l0yy(double x, double y)
    {
      return  Legendre3(x) * Legendre0xx(y);
    }

    static double leg_quad_l3_l1(double x, double y)
    {
      return Legendre3(x) * Legendre1(y);
    }

    static double leg_quad_l3_l1x(double x, double y)
    {
      return Legendre3x(x) * Legendre1(y);
    }

    static double leg_quad_l3_l1y(double x, double y)
    {
      return Legendre3(x) * Legendre1x(y);
    }

    static double leg_quad_l3_l1xx(double x, double y)
    {
      return Legendre3xx(x) * Legendre1(y);
    }

    static double leg_quad_l3_l1xy(double x, double y)
    {
      return Legendre3x(x) * Legendre1x(y);
    }

    static double leg_quad_l3_l1yy(double x, double y)
    {
      return  Legendre3(x) * Legendre1xx(y);
    }

    static double leg_quad_l3_l2(double x, double y)
    {
      return Legendre3(x) * Legendre2(y);
    }

    static double leg_quad_l3_l2x(double x, double y)
    {
      return Legendre3x(x) * Legendre2(y);
    }

    static double leg_quad_l3_l2y(double x, double y)
    {
      return Legendre3(x) * Legendre2x(y);
    }

    static double leg_quad_l3_l2xx(double x, double y)
    {
      return Legendre3xx(x) * Legendre2(y);
    }

    static double leg_quad_l3_l2xy(double x, double y)
    {
      return Legendre3x(x) * Legendre2x(y);
    }

    static double leg_quad_l3_l2yy(double x, double y)
    {
      return  Legendre3(x) * Legendre2xx(y);
    }

    static double leg_quad_l3_l3(double x, double y)
    {
      return Legendre3(x) * Legendre3(y);
    }

    static double leg_quad_l3_l3x(double x, double y)
    {
      return Legendre3x(x) * Legendre3(y);
    }

    static double leg_quad_l3_l3y(double x, double y)
    {
      return Legendre3(x) * Legendre3x(y);
    }

    static double leg_quad_l3_l3xx(double x, double y)
    {
      return Legendre3xx(x) * Legendre3(y);
    }

    static double leg_quad_l3_l3xy(double x, double y)
    {
      return Legendre3x(x) * Legendre3x(y);
    }

    static double leg_quad_l3_l3yy(double x, double y)
    {
      return  Legendre3(x) * Legendre3xx(y);
    }

    static double leg_quad_l3_l4(double x, double y)
    {
      return Legendre3(x) * Legendre4(y);
    }

    static double leg_quad_l3_l4x(double x, double y)
    {
      return Legendre3x(x) * Legendre4(y);
    }

    static double leg_quad_l3_l4y(double x, double y)
    {
      return Legendre3(x) * Legendre4x(y);
    }

    static double leg_quad_l3_l4xx(double x, double y)
    {
      return Legendre3xx(x) * Legendre4(y);
    }

    static double leg_quad_l3_l4xy(double x, double y)
    {
      return Legendre3x(x) * Legendre4x(y);
    }

    static double leg_quad_l3_l4yy(double x, double y)
    {
      return  Legendre3(x) * Legendre4xx(y);
    }

    static double leg_quad_l3_l5(double x, double y)
    {
      return Legendre3(x) * Legendre5(y);
    }

    static double leg_quad_l3_l5x(double x, double y)
    {
      return Legendre3x(x) * Legendre5(y);
    }

    static double leg_quad_l3_l5y(double x, double y)
    {
      return Legendre3(x) * Legendre5x(y);
    }

    static double leg_quad_l3_l5xx(double x, double y)
    {
      return Legendre3xx(x) * Legendre5(y);
    }

    static double leg_quad_l3_l5xy(double x, double y)
    {
      return Legendre3x(x) * Legendre5x(y);
    }

    static double leg_quad_l3_l5yy(double x, double y)
    {
      return  Legendre3(x) * Legendre5xx(y);
    }

    static double leg_quad_l3_l6(double x, double y)
    {
      return Legendre3(x) * Legendre6(y);
    }

    static double leg_quad_l3_l6x(double x, double y)
    {
      return Legendre3x(x) * Legendre6(y);
    }

    static double leg_quad_l3_l6y(double x, double y)
    {
      return Legendre3(x) * Legendre6x(y);
    }

    static double leg_quad_l3_l6xx(double x, double y)
    {
      return Legendre3xx(x) * Legendre6(y);
    }

    static double leg_quad_l3_l6xy(double x, double y)
    {
      return Legendre3x(x) * Legendre6x(y);
    }

    static double leg_quad_l3_l6yy(double x, double y)
    {
      return  Legendre3(x) * Legendre6xx(y);
    }

    static double leg_quad_l3_l7(double x, double y)
    {
      return Legendre3(x) * Legendre7(y);
    }

    static double leg_quad_l3_l7x(double x, double y)
    {
      return Legendre3x(x) * Legendre7(y);
    }

    static double leg_quad_l3_l7y(double x, double y)
    {
      return Legendre3(x) * Legendre7x(y);
    }

    static double leg_quad_l3_l7xx(double x, double y)
    {
      return Legendre3xx(x) * Legendre7(y);
    }

    static double leg_quad_l3_l7xy(double x, double y)
    {
      return Legendre3x(x) * Legendre7x(y);
    }

    static double leg_quad_l3_l7yy(double x, double y)
    {
      return  Legendre3(x) * Legendre7xx(y);
    }

    static double leg_quad_l3_l8(double x, double y)
    {
      return Legendre3(x) * Legendre8(y);
    }

    static double leg_quad_l3_l8x(double x, double y)
    {
      return Legendre3x(x) * Legendre8(y);
    }

    static double leg_quad_l3_l8y(double x, double y)
    {
      return Legendre3(x) * Legendre8x(y);
    }

    static double leg_quad_l3_l8xx(double x, double y)
    {
      return Legendre3xx(x) * Legendre8(y);
    }

    static double leg_quad_l3_l8xy(double x, double y)
    {
      return Legendre3x(x) * Legendre8x(y);
    }

    static double leg_quad_l3_l8yy(double x, double y)
    {
      return  Legendre3(x) * Legendre8xx(y);
    }

    static double leg_quad_l3_l9(double x, double y)
    {
      return Legendre3(x) * Legendre9(y);
    }

    static double leg_quad_l3_l9x(double x, double y)
    {
      return Legendre3x(x) * Legendre9(y);
    }

    static double leg_quad_l3_l9y(double x, double y)
    {
      return Legendre3(x) * Legendre9x(y);
    }

    static double leg_quad_l3_l9xx(double x, double y)
    {
      return Legendre3xx(x) * Legendre9(y);
    }

    static double leg_quad_l3_l9xy(double x, double y)
    {
      return Legendre3x(x) * Legendre9x(y);
    }

    static double leg_quad_l3_l9yy(double x, double y)
    {
      return  Legendre3(x) * Legendre9xx(y);
    }

    static double leg_quad_l3_l10(double x, double y)
    {
      return Legendre3(x) * Legendre10(y);
    }

    static double leg_quad_l3_l10x(double x, double y)
    {
      return Legendre3x(x) * Legendre10(y);
    }

    static double leg_quad_l3_l10y(double x, double y)
    {
      return Legendre3(x) * Legendre10x(y);
    }

    static double leg_quad_l3_l10xx(double x, double y)
    {
      return Legendre3xx(x) * Legendre10(y);
    }

    static double leg_quad_l3_l10xy(double x, double y)
    {
      return Legendre3x(x) * Legendre10x(y);
    }

    static double leg_quad_l3_l10yy(double x, double y)
    {
      return  Legendre3(x) * Legendre10xx(y);
    }

    static double leg_quad_l4_l0(double x, double y)
    {
      return Legendre4(x) * Legendre0(y);
    }

    static double leg_quad_l4_l0x(double x, double y)
    {
      return Legendre4x(x) * Legendre0(y);
    }

    static double leg_quad_l4_l0y(double x, double y)
    {
      return Legendre4(x) * Legendre0x(y);
    }

    static double leg_quad_l4_l0xx(double x, double y)
    {
      return Legendre4xx(x) * Legendre0(y);
    }

    static double leg_quad_l4_l0xy(double x, double y)
    {
      return Legendre4x(x) * Legendre0x(y);
    }

    static double leg_quad_l4_l0yy(double x, double y)
    {
      return  Legendre4(x) * Legendre0xx(y);
    }

    static double leg_quad_l4_l1(double x, double y)
    {
      return Legendre4(x) * Legendre1(y);
    }

    static double leg_quad_l4_l1x(double x, double y)
    {
      return Legendre4x(x) * Legendre1(y);
    }

    static double leg_quad_l4_l1y(double x, double y)
    {
      return Legendre4(x) * Legendre1x(y);
    }

    static double leg_quad_l4_l1xx(double x, double y)
    {
      return Legendre4xx(x) * Legendre1(y);
    }

    static double leg_quad_l4_l1xy(double x, double y)
    {
      return Legendre4x(x) * Legendre1x(y);
    }

    static double leg_quad_l4_l1yy(double x, double y)
    {
      return  Legendre4(x) * Legendre1xx(y);
    }

    static double leg_quad_l4_l2(double x, double y)
    {
      return Legendre4(x) * Legendre2(y);
    }

    static double leg_quad_l4_l2x(double x, double y)
    {
      return Legendre4x(x) * Legendre2(y);
    }

    static double leg_quad_l4_l2y(double x, double y)
    {
      return Legendre4(x) * Legendre2x(y);
    }

    static double leg_quad_l4_l2xx(double x, double y)
    {
      return Legendre4xx(x) * Legendre2(y);
    }

    static double leg_quad_l4_l2xy(double x, double y)
    {
      return Legendre4x(x) * Legendre2x(y);
    }

    static double leg_quad_l4_l2yy(double x, double y)
    {
      return  Legendre4(x) * Legendre2xx(y);
    }

    static double leg_quad_l4_l3(double x, double y)
    {
      return Legendre4(x) * Legendre3(y);
    }

    static double leg_quad_l4_l3x(double x, double y)
    {
      return Legendre4x(x) * Legendre3(y);
    }

    static double leg_quad_l4_l3y(double x, double y)
    {
      return Legendre4(x) * Legendre3x(y);
    }

    static double leg_quad_l4_l3xx(double x, double y)
    {
      return Legendre4xx(x) * Legendre3(y);
    }

    static double leg_quad_l4_l3xy(double x, double y)
    {
      return Legendre4x(x) * Legendre3x(y);
    }

    static double leg_quad_l4_l3yy(double x, double y)
    {
      return  Legendre4(x) * Legendre3xx(y);
    }

    static double leg_quad_l4_l4(double x, double y)
    {
      return Legendre4(x) * Legendre4(y);
    }

    static double leg_quad_l4_l4x(double x, double y)
    {
      return Legendre4x(x) * Legendre4(y);
    }

    static double leg_quad_l4_l4y(double x, double y)
    {
      return Legendre4(x) * Legendre4x(y);
    }

    static double leg_quad_l4_l4xx(double x, double y)
    {
      return Legendre4xx(x) * Legendre4(y);
    }

    static double leg_quad_l4_l4xy(double x, double y)
    {
      return Legendre4x(x) * Legendre4x(y);
    }

    static double leg_quad_l4_l4yy(double x, double y)
    {
      return  Legendre4(x) * Legendre4xx(y);
    }

    static double leg_quad_l4_l5(double x, double y)
    {
      return Legendre4(x) * Legendre5(y);
    }

    static double leg_quad_l4_l5x(double x, double y)
    {
      return Legendre4x(x) * Legendre5(y);
    }

    static double leg_quad_l4_l5y(double x, double y)
    {
      return Legendre4(x) * Legendre5x(y);
    }

    static double leg_quad_l4_l5xx(double x, double y)
    {
      return Legendre4xx(x) * Legendre5(y);
    }

    static double leg_quad_l4_l5xy(double x, double y)
    {
      return Legendre4x(x) * Legendre5x(y);
    }

    static double leg_quad_l4_l5yy(double x, double y)
    {
      return  Legendre4(x) * Legendre5xx(y);
    }

    static double leg_quad_l4_l6(double x, double y)
    {
      return Legendre4(x) * Legendre6(y);
    }

    static double leg_quad_l4_l6x(double x, double y)
    {
      return Legendre4x(x) * Legendre6(y);
    }

    static double leg_quad_l4_l6y(double x, double y)
    {
      return Legendre4(x) * Legendre6x(y);
    }

    static double leg_quad_l4_l6xx(double x, double y)
    {
      return Legendre4xx(x) * Legendre6(y);
    }

    static double leg_quad_l4_l6xy(double x, double y)
    {
      return Legendre4x(x) * Legendre6x(y);
    }

    static double leg_quad_l4_l6yy(double x, double y)
    {
      return  Legendre4(x) * Legendre6xx(y);
    }

    static double leg_quad_l4_l7(double x, double y)
    {
      return Legendre4(x) * Legendre7(y);
    }

    static double leg_quad_l4_l7x(double x, double y)
    {
      return Legendre4x(x) * Legendre7(y);
    }

    static double leg_quad_l4_l7y(double x, double y)
    {
      return Legendre4(x) * Legendre7x(y);
    }

    static double leg_quad_l4_l7xx(double x, double y)
    {
      return Legendre4xx(x) * Legendre7(y);
    }

    static double leg_quad_l4_l7xy(double x, double y)
    {
      return Legendre4x(x) * Legendre7x(y);
    }

    static double leg_quad_l4_l7yy(double x, double y)
    {
      return  Legendre4(x) * Legendre7xx(y);
    }

    static double leg_quad_l4_l8(double x, double y)
    {
      return Legendre4(x) * Legendre8(y);
    }

    static double leg_quad_l4_l8x(double x, double y)
    {
      return Legendre4x(x) * Legendre8(y);
    }

    static double leg_quad_l4_l8y(double x, double y)
    {
      return Legendre4(x) * Legendre8x(y);
    }

    static double leg_quad_l4_l8xx(double x, double y)
    {
      return Legendre4xx(x) * Legendre8(y);
    }

    static double leg_quad_l4_l8xy(double x, double y)
    {
      return Legendre4x(x) * Legendre8x(y);
    }

    static double leg_quad_l4_l8yy(double x, double y)
    {
      return  Legendre4(x) * Legendre8xx(y);
    }

    static double leg_quad_l4_l9(double x, double y)
    {
      return Legendre4(x) * Legendre9(y);
    }

    static double leg_quad_l4_l9x(double x, double y)
    {
      return Legendre4x(x) * Legendre9(y);
    }

    static double leg_quad_l4_l9y(double x, double y)
    {
      return Legendre4(x) * Legendre9x(y);
    }

    static double leg_quad_l4_l9xx(double x, double y)
    {
      return Legendre4xx(x) * Legendre9(y);
    }

    static double leg_quad_l4_l9xy(double x, double y)
    {
      return Legendre4x(x) * Legendre9x(y);
    }

    static double leg_quad_l4_l9yy(double x, double y)
    {
      return  Legendre4(x) * Legendre9xx(y);
    }

    static double leg_quad_l4_l10(double x, double y)
    {
      return Legendre4(x) * Legendre10(y);
    }

    static double leg_quad_l4_l10x(double x, double y)
    {
      return Legendre4x(x) * Legendre10(y);
    }

    static double leg_quad_l4_l10y(double x, double y)
    {
      return Legendre4(x) * Legendre10x(y);
    }

    static double leg_quad_l4_l10xx(double x, double y)
    {
      return Legendre4xx(x) * Legendre10(y);
    }

    static double leg_quad_l4_l10xy(double x, double y)
    {
      return Legendre4x(x) * Legendre10x(y);
    }

    static double leg_quad_l4_l10yy(double x, double y)
    {
      return  Legendre4(x) * Legendre10xx(y);
    }

    static double leg_quad_l5_l0(double x, double y)
    {
      return Legendre5(x) * Legendre0(y);
    }

    static double leg_quad_l5_l0x(double x, double y)
    {
      return Legendre5x(x) * Legendre0(y);
    }

    static double leg_quad_l5_l0y(double x, double y)
    {
      return Legendre5(x) * Legendre0x(y);
    }

    static double leg_quad_l5_l0xx(double x, double y)
    {
      return Legendre5xx(x) * Legendre0(y);
    }

    static double leg_quad_l5_l0xy(double x, double y)
    {
      return Legendre5x(x) * Legendre0x(y);
    }

    static double leg_quad_l5_l0yy(double x, double y)
    {
      return  Legendre5(x) * Legendre0xx(y);
    }

    static double leg_quad_l5_l1(double x, double y)
    {
      return Legendre5(x) * Legendre1(y);
    }

    static double leg_quad_l5_l1x(double x, double y)
    {
      return Legendre5x(x) * Legendre1(y);
    }

    static double leg_quad_l5_l1y(double x, double y)
    {
      return Legendre5(x) * Legendre1x(y);
    }

    static double leg_quad_l5_l1xx(double x, double y)
    {
      return Legendre5xx(x) * Legendre1(y);
    }

    static double leg_quad_l5_l1xy(double x, double y)
    {
      return Legendre5x(x) * Legendre1x(y);
    }

    static double leg_quad_l5_l1yy(double x, double y)
    {
      return  Legendre5(x) * Legendre1xx(y);
    }

    static double leg_quad_l5_l2(double x, double y)
    {
      return Legendre5(x) * Legendre2(y);
    }

    static double leg_quad_l5_l2x(double x, double y)
    {
      return Legendre5x(x) * Legendre2(y);
    }

    static double leg_quad_l5_l2y(double x, double y)
    {
      return Legendre5(x) * Legendre2x(y);
    }

    static double leg_quad_l5_l2xx(double x, double y)
    {
      return Legendre5xx(x) * Legendre2(y);
    }

    static double leg_quad_l5_l2xy(double x, double y)
    {
      return Legendre5x(x) * Legendre2x(y);
    }

    static double leg_quad_l5_l2yy(double x, double y)
    {
      return  Legendre5(x) * Legendre2xx(y);
    }

    static double leg_quad_l5_l3(double x, double y)
    {
      return Legendre5(x) * Legendre3(y);
    }

    static double leg_quad_l5_l3x(double x, double y)
    {
      return Legendre5x(x) * Legendre3(y);
    }

    static double leg_quad_l5_l3y(double x, double y)
    {
      return Legendre5(x) * Legendre3x(y);
    }

    static double leg_quad_l5_l3xx(double x, double y)
    {
      return Legendre5xx(x) * Legendre3(y);
    }

    static double leg_quad_l5_l3xy(double x, double y)
    {
      return Legendre5x(x) * Legendre3x(y);
    }

    static double leg_quad_l5_l3yy(double x, double y)
    {
      return  Legendre5(x) * Legendre3xx(y);
    }

    static double leg_quad_l5_l4(double x, double y)
    {
      return Legendre5(x) * Legendre4(y);
    }

    static double leg_quad_l5_l4x(double x, double y)
    {
      return Legendre5x(x) * Legendre4(y);
    }

    static double leg_quad_l5_l4y(double x, double y)
    {
      return Legendre5(x) * Legendre4x(y);
    }

    static double leg_quad_l5_l4xx(double x, double y)
    {
      return Legendre5xx(x) * Legendre4(y);
    }

    static double leg_quad_l5_l4xy(double x, double y)
    {
      return Legendre5x(x) * Legendre4x(y);
    }

    static double leg_quad_l5_l4yy(double x, double y)
    {
      return  Legendre5(x) * Legendre4xx(y);
    }

    static double leg_quad_l5_l5(double x, double y)
    {
      return Legendre5(x) * Legendre5(y);
    }

    static double leg_quad_l5_l5x(double x, double y)
    {
      return Legendre5x(x) * Legendre5(y);
    }

    static double leg_quad_l5_l5y(double x, double y)
    {
      return Legendre5(x) * Legendre5x(y);
    }

    static double leg_quad_l5_l5xx(double x, double y)
    {
      return Legendre5xx(x) * Legendre5(y);
    }

    static double leg_quad_l5_l5xy(double x, double y)
    {
      return Legendre5x(x) * Legendre5x(y);
    }

    static double leg_quad_l5_l5yy(double x, double y)
    {
      return  Legendre5(x) * Legendre5xx(y);
    }

    static double leg_quad_l5_l6(double x, double y)
    {
      return Legendre5(x) * Legendre6(y);
    }

    static double leg_quad_l5_l6x(double x, double y)
    {
      return Legendre5x(x) * Legendre6(y);
    }

    static double leg_quad_l5_l6y(double x, double y)
    {
      return Legendre5(x) * Legendre6x(y);
    }

    static double leg_quad_l5_l6xx(double x, double y)
    {
      return Legendre5xx(x) * Legendre6(y);
    }

    static double leg_quad_l5_l6xy(double x, double y)
    {
      return Legendre5x(x) * Legendre6x(y);
    }

    static double leg_quad_l5_l6yy(double x, double y)
    {
      return  Legendre5(x) * Legendre6xx(y);
    }

    static double leg_quad_l5_l7(double x, double y)
    {
      return Legendre5(x) * Legendre7(y);
    }

    static double leg_quad_l5_l7x(double x, double y)
    {
      return Legendre5x(x) * Legendre7(y);
    }

    static double leg_quad_l5_l7y(double x, double y)
    {
      return Legendre5(x) * Legendre7x(y);
    }

    static double leg_quad_l5_l7xx(double x, double y)
    {
      return Legendre5xx(x) * Legendre7(y);
    }

    static double leg_quad_l5_l7xy(double x, double y)
    {
      return Legendre5x(x) * Legendre7x(y);
    }

    static double leg_quad_l5_l7yy(double x, double y)
    {
      return  Legendre5(x) * Legendre7xx(y);
    }

    static double leg_quad_l5_l8(double x, double y)
    {
      return Legendre5(x) * Legendre8(y);
    }

    static double leg_quad_l5_l8x(double x, double y)
    {
      return Legendre5x(x) * Legendre8(y);
    }

    static double leg_quad_l5_l8y(double x, double y)
    {
      return Legendre5(x) * Legendre8x(y);
    }

    static double leg_quad_l5_l8xx(double x, double y)
    {
      return Legendre5xx(x) * Legendre8(y);
    }

    static double leg_quad_l5_l8xy(double x, double y)
    {
      return Legendre5x(x) * Legendre8x(y);
    }

    static double leg_quad_l5_l8yy(double x, double y)
    {
      return  Legendre5(x) * Legendre8xx(y);
    }

    static double leg_quad_l5_l9(double x, double y)
    {
      return Legendre5(x) * Legendre9(y);
    }

    static double leg_quad_l5_l9x(double x, double y)
    {
      return Legendre5x(x) * Legendre9(y);
    }

    static double leg_quad_l5_l9y(double x, double y)
    {
      return Legendre5(x) * Legendre9x(y);
    }

    static double leg_quad_l5_l9xx(double x, double y)
    {
      return Legendre5xx(x) * Legendre9(y);
    }

    static double leg_quad_l5_l9xy(double x, double y)
    {
      return Legendre5x(x) * Legendre9x(y);
    }

    static double leg_quad_l5_l9yy(double x, double y)
    {
      return  Legendre5(x) * Legendre9xx(y);
    }

    static double leg_quad_l5_l10(double x, double y)
    {
      return Legendre5(x) * Legendre10(y);
    }

    static double leg_quad_l5_l10x(double x, double y)
    {
      return Legendre5x(x) * Legendre10(y);
    }

    static double leg_quad_l5_l10y(double x, double y)
    {
      return Legendre5(x) * Legendre10x(y);
    }

    static double leg_quad_l5_l10xx(double x, double y)
    {
      return Legendre5xx(x) * Legendre10(y);
    }

    static double leg_quad_l5_l10xy(double x, double y)
    {
      return Legendre5x(x) * Legendre10x(y);
    }

    static double leg_quad_l5_l10yy(double x, double y)
    {
      return  Legendre5(x) * Legendre10xx(y);
    }

    static double leg_quad_l6_l0(double x, double y)
    {
      return Legendre6(x) * Legendre0(y);
    }

    static double leg_quad_l6_l0x(double x, double y)
    {
      return Legendre6x(x) * Legendre0(y);
    }

    static double leg_quad_l6_l0y(double x, double y)
    {
      return Legendre6(x) * Legendre0x(y);
    }

    static double leg_quad_l6_l0xx(double x, double y)
    {
      return Legendre6xx(x) * Legendre0(y);
    }

    static double leg_quad_l6_l0xy(double x, double y)
    {
      return Legendre6x(x) * Legendre0x(y);
    }

    static double leg_quad_l6_l0yy(double x, double y)
    {
      return  Legendre6(x) * Legendre0xx(y);
    }

    static double leg_quad_l6_l1(double x, double y)
    {
      return Legendre6(x) * Legendre1(y);
    }

    static double leg_quad_l6_l1x(double x, double y)
    {
      return Legendre6x(x) * Legendre1(y);
    }

    static double leg_quad_l6_l1y(double x, double y)
    {
      return Legendre6(x) * Legendre1x(y);
    }

    static double leg_quad_l6_l1xx(double x, double y)
    {
      return Legendre6xx(x) * Legendre1(y);
    }

    static double leg_quad_l6_l1xy(double x, double y)
    {
      return Legendre6x(x) * Legendre1x(y);
    }

    static double leg_quad_l6_l1yy(double x, double y)
    {
      return  Legendre6(x) * Legendre1xx(y);
    }

    static double leg_quad_l6_l2(double x, double y)
    {
      return Legendre6(x) * Legendre2(y);
    }

    static double leg_quad_l6_l2x(double x, double y)
    {
      return Legendre6x(x) * Legendre2(y);
    }

    static double leg_quad_l6_l2y(double x, double y)
    {
      return Legendre6(x) * Legendre2x(y);
    }

    static double leg_quad_l6_l2xx(double x, double y)
    {
      return Legendre6xx(x) * Legendre2(y);
    }

    static double leg_quad_l6_l2xy(double x, double y)
    {
      return Legendre6x(x) * Legendre2x(y);
    }

    static double leg_quad_l6_l2yy(double x, double y)
    {
      return  Legendre6(x) * Legendre2xx(y);
    }

    static double leg_quad_l6_l3(double x, double y)
    {
      return Legendre6(x) * Legendre3(y);
    }

    static double leg_quad_l6_l3x(double x, double y)
    {
      return Legendre6x(x) * Legendre3(y);
    }

    static double leg_quad_l6_l3y(double x, double y)
    {
      return Legendre6(x) * Legendre3x(y);
    }

    static double leg_quad_l6_l3xx(double x, double y)
    {
      return Legendre6xx(x) * Legendre3(y);
    }

    static double leg_quad_l6_l3xy(double x, double y)
    {
      return Legendre6x(x) * Legendre3x(y);
    }

    static double leg_quad_l6_l3yy(double x, double y)
    {
      return  Legendre6(x) * Legendre3xx(y);
    }

    static double leg_quad_l6_l4(double x, double y)
    {
      return Legendre6(x) * Legendre4(y);
    }

    static double leg_quad_l6_l4x(double x, double y)
    {
      return Legendre6x(x) * Legendre4(y);
    }

    static double leg_quad_l6_l4y(double x, double y)
    {
      return Legendre6(x) * Legendre4x(y);
    }

    static double leg_quad_l6_l4xx(double x, double y)
    {
      return Legendre6xx(x) * Legendre4(y);
    }

    static double leg_quad_l6_l4xy(double x, double y)
    {
      return Legendre6x(x) * Legendre4x(y);
    }

    static double leg_quad_l6_l4yy(double x, double y)
    {
      return  Legendre6(x) * Legendre4xx(y);
    }

    static double leg_quad_l6_l5(double x, double y)
    {
      return Legendre6(x) * Legendre5(y);
    }

    static double leg_quad_l6_l5x(double x, double y)
    {
      return Legendre6x(x) * Legendre5(y);
    }

    static double leg_quad_l6_l5y(double x, double y)
    {
      return Legendre6(x) * Legendre5x(y);
    }

    static double leg_quad_l6_l5xx(double x, double y)
    {
      return Legendre6xx(x) * Legendre5(y);
    }

    static double leg_quad_l6_l5xy(double x, double y)
    {
      return Legendre6x(x) * Legendre5x(y);
    }

    static double leg_quad_l6_l5yy(double x, double y)
    {
      return  Legendre6(x) * Legendre5xx(y);
    }

    static double leg_quad_l6_l6(double x, double y)
    {
      return Legendre6(x) * Legendre6(y);
    }

    static double leg_quad_l6_l6x(double x, double y)
    {
      return Legendre6x(x) * Legendre6(y);
    }

    static double leg_quad_l6_l6y(double x, double y)
    {
      return Legendre6(x) * Legendre6x(y);
    }

    static double leg_quad_l6_l6xx(double x, double y)
    {
      return Legendre6xx(x) * Legendre6(y);
    }

    static double leg_quad_l6_l6xy(double x, double y)
    {
      return Legendre6x(x) * Legendre6x(y);
    }

    static double leg_quad_l6_l6yy(double x, double y)
    {
      return  Legendre6(x) * Legendre6xx(y);
    }

    static double leg_quad_l6_l7(double x, double y)
    {
      return Legendre6(x) * Legendre7(y);
    }

    static double leg_quad_l6_l7x(double x, double y)
    {
      return Legendre6x(x) * Legendre7(y);
    }

    static double leg_quad_l6_l7y(double x, double y)
    {
      return Legendre6(x) * Legendre7x(y);
    }

    static double leg_quad_l6_l7xx(double x, double y)
    {
      return Legendre6xx(x) * Legendre7(y);
    }

    static double leg_quad_l6_l7xy(double x, double y)
    {
      return Legendre6x(x) * Legendre7x(y);
    }

    static double leg_quad_l6_l7yy(double x, double y)
    {
      return  Legendre6(x) * Legendre7xx(y);
    }

    static double leg_quad_l6_l8(double x, double y)
    {
      return Legendre6(x) * Legendre8(y);
    }

    static double leg_quad_l6_l8x(double x, double y)
    {
      return Legendre6x(x) * Legendre8(y);
    }

    static double leg_quad_l6_l8y(double x, double y)
    {
      return Legendre6(x) * Legendre8x(y);
    }

    static double leg_quad_l6_l8xx(double x, double y)
    {
      return Legendre6xx(x) * Legendre8(y);
    }

    static double leg_quad_l6_l8xy(double x, double y)
    {
      return Legendre6x(x) * Legendre8x(y);
    }

    static double leg_quad_l6_l8yy(double x, double y)
    {
      return  Legendre6(x) * Legendre8xx(y);
    }

    static double leg_quad_l6_l9(double x, double y)
    {
      return Legendre6(x) * Legendre9(y);
    }

    static double leg_quad_l6_l9x(double x, double y)
    {
      return Legendre6x(x) * Legendre9(y);
    }

    static double leg_quad_l6_l9y(double x, double y)
    {
      return Legendre6(x) * Legendre9x(y);
    }

    static double leg_quad_l6_l9xx(double x, double y)
    {
      return Legendre6xx(x) * Legendre9(y);
    }

    static double leg_quad_l6_l9xy(double x, double y)
    {
      return Legendre6x(x) * Legendre9x(y);
    }

    static double leg_quad_l6_l9yy(double x, double y)
    {
      return  Legendre6(x) * Legendre9xx(y);
    }

    static double leg_quad_l6_l10(double x, double y)
    {
      return Legendre6(x) * Legendre10(y);
    }

    static double leg_quad_l6_l10x(double x, double y)
    {
      return Legendre6x(x) * Legendre10(y);
    }

    static double leg_quad_l6_l10y(double x, double y)
    {
      return Legendre6(x) * Legendre10x(y);
    }

    static double leg_quad_l6_l10xx(double x, double y)
    {
      return Legendre6xx(x) * Legendre10(y);
    }

    static double leg_quad_l6_l10xy(double x, double y)
    {
      return Legendre6x(x) * Legendre10x(y);
    }

    static double leg_quad_l6_l10yy(double x, double y)
    {
      return  Legendre6(x) * Legendre10xx(y);
    }

    static double leg_quad_l7_l0(double x, double y)
    {
      return Legendre7(x) * Legendre0(y);
    }

    static double leg_quad_l7_l0x(double x, double y)
    {
      return Legendre7x(x) * Legendre0(y);
    }

    static double leg_quad_l7_l0y(double x, double y)
    {
      return Legendre7(x) * Legendre0x(y);
    }

    static double leg_quad_l7_l0xx(double x, double y)
    {
      return Legendre7xx(x) * Legendre0(y);
    }

    static double leg_quad_l7_l0xy(double x, double y)
    {
      return Legendre7x(x) * Legendre0x(y);
    }

    static double leg_quad_l7_l0yy(double x, double y)
    {
      return  Legendre7(x) * Legendre0xx(y);
    }

    static double leg_quad_l7_l1(double x, double y)
    {
      return Legendre7(x) * Legendre1(y);
    }

    static double leg_quad_l7_l1x(double x, double y)
    {
      return Legendre7x(x) * Legendre1(y);
    }

    static double leg_quad_l7_l1y(double x, double y)
    {
      return Legendre7(x) * Legendre1x(y);
    }

    static double leg_quad_l7_l1xx(double x, double y)
    {
      return Legendre7xx(x) * Legendre1(y);
    }

    static double leg_quad_l7_l1xy(double x, double y)
    {
      return Legendre7x(x) * Legendre1x(y);
    }

    static double leg_quad_l7_l1yy(double x, double y)
    {
      return  Legendre7(x) * Legendre1xx(y);
    }

    static double leg_quad_l7_l2(double x, double y)
    {
      return Legendre7(x) * Legendre2(y);
    }

    static double leg_quad_l7_l2x(double x, double y)
    {
      return Legendre7x(x) * Legendre2(y);
    }

    static double leg_quad_l7_l2y(double x, double y)
    {
      return Legendre7(x) * Legendre2x(y);
    }

    static double leg_quad_l7_l2xx(double x, double y)
    {
      return Legendre7xx(x) * Legendre2(y);
    }

    static double leg_quad_l7_l2xy(double x, double y)
    {
      return Legendre7x(x) * Legendre2x(y);
    }

    static double leg_quad_l7_l2yy(double x, double y)
    {
      return  Legendre7(x) * Legendre2xx(y);
    }

    static double leg_quad_l7_l3(double x, double y)
    {
      return Legendre7(x) * Legendre3(y);
    }

    static double leg_quad_l7_l3x(double x, double y)
    {
      return Legendre7x(x) * Legendre3(y);
    }

    static double leg_quad_l7_l3y(double x, double y)
    {
      return Legendre7(x) * Legendre3x(y);
    }

    static double leg_quad_l7_l3xx(double x, double y)
    {
      return Legendre7xx(x) * Legendre3(y);
    }

    static double leg_quad_l7_l3xy(double x, double y)
    {
      return Legendre7x(x) * Legendre3x(y);
    }

    static double leg_quad_l7_l3yy(double x, double y)
    {
      return  Legendre7(x) * Legendre3xx(y);
    }

    static double leg_quad_l7_l4(double x, double y)
    {
      return Legendre7(x) * Legendre4(y);
    }

    static double leg_quad_l7_l4x(double x, double y)
    {
      return Legendre7x(x) * Legendre4(y);
    }

    static double leg_quad_l7_l4y(double x, double y)
    {
      return Legendre7(x) * Legendre4x(y);
    }

    static double leg_quad_l7_l4xx(double x, double y)
    {
      return Legendre7xx(x) * Legendre4(y);
    }

    static double leg_quad_l7_l4xy(double x, double y)
    {
      return Legendre7x(x) * Legendre4x(y);
    }

    static double leg_quad_l7_l4yy(double x, double y)
    {
      return  Legendre7(x) * Legendre4xx(y);
    }

    static double leg_quad_l7_l5(double x, double y)
    {
      return Legendre7(x) * Legendre5(y);
    }

    static double leg_quad_l7_l5x(double x, double y)
    {
      return Legendre7x(x) * Legendre5(y);
    }

    static double leg_quad_l7_l5y(double x, double y)
    {
      return Legendre7(x) * Legendre5x(y);
    }

    static double leg_quad_l7_l5xx(double x, double y)
    {
      return Legendre7xx(x) * Legendre5(y);
    }

    static double leg_quad_l7_l5xy(double x, double y)
    {
      return Legendre7x(x) * Legendre5x(y);
    }

    static double leg_quad_l7_l5yy(double x, double y)
    {
      return  Legendre7(x) * Legendre5xx(y);
    }

    static double leg_quad_l7_l6(double x, double y)
    {
      return Legendre7(x) * Legendre6(y);
    }

    static double leg_quad_l7_l6x(double x, double y)
    {
      return Legendre7x(x) * Legendre6(y);
    }

    static double leg_quad_l7_l6y(double x, double y)
    {
      return Legendre7(x) * Legendre6x(y);
    }

    static double leg_quad_l7_l6xx(double x, double y)
    {
      return Legendre7xx(x) * Legendre6(y);
    }

    static double leg_quad_l7_l6xy(double x, double y)
    {
      return Legendre7x(x) * Legendre6x(y);
    }

    static double leg_quad_l7_l6yy(double x, double y)
    {
      return  Legendre7(x) * Legendre6xx(y);
    }

    static double leg_quad_l7_l7(double x, double y)
    {
      return Legendre7(x) * Legendre7(y);
    }

    static double leg_quad_l7_l7x(double x, double y)
    {
      return Legendre7x(x) * Legendre7(y);
    }

    static double leg_quad_l7_l7y(double x, double y)
    {
      return Legendre7(x) * Legendre7x(y);
    }

    static double leg_quad_l7_l7xx(double x, double y)
    {
      return Legendre7xx(x) * Legendre7(y);
    }

    static double leg_quad_l7_l7xy(double x, double y)
    {
      return Legendre7x(x) * Legendre7x(y);
    }

    static double leg_quad_l7_l7yy(double x, double y)
    {
      return  Legendre7(x) * Legendre7xx(y);
    }

    static double leg_quad_l7_l8(double x, double y)
    {
      return Legendre7(x) * Legendre8(y);
    }

    static double leg_quad_l7_l8x(double x, double y)
    {
      return Legendre7x(x) * Legendre8(y);
    }

    static double leg_quad_l7_l8y(double x, double y)
    {
      return Legendre7(x) * Legendre8x(y);
    }

    static double leg_quad_l7_l8xx(double x, double y)
    {
      return Legendre7xx(x) * Legendre8(y);
    }

    static double leg_quad_l7_l8xy(double x, double y)
    {
      return Legendre7x(x) * Legendre8x(y);
    }

    static double leg_quad_l7_l8yy(double x, double y)
    {
      return  Legendre7(x) * Legendre8xx(y);
    }

    static double leg_quad_l7_l9(double x, double y)
    {
      return Legendre7(x) * Legendre9(y);
    }

    static double leg_quad_l7_l9x(double x, double y)
    {
      return Legendre7x(x) * Legendre9(y);
    }

    static double leg_quad_l7_l9y(double x, double y)
    {
      return Legendre7(x) * Legendre9x(y);
    }

    static double leg_quad_l7_l9xx(double x, double y)
    {
      return Legendre7xx(x) * Legendre9(y);
    }

    static double leg_quad_l7_l9xy(double x, double y)
    {
      return Legendre7x(x) * Legendre9x(y);
    }

    static double leg_quad_l7_l9yy(double x, double y)
    {
      return  Legendre7(x) * Legendre9xx(y);
    }

    static double leg_quad_l7_l10(double x, double y)
    {
      return Legendre7(x) * Legendre10(y);
    }

    static double leg_quad_l7_l10x(double x, double y)
    {
      return Legendre7x(x) * Legendre10(y);
    }

    static double leg_quad_l7_l10y(double x, double y)
    {
      return Legendre7(x) * Legendre10x(y);
    }

    static double leg_quad_l7_l10xx(double x, double y)
    {
      return Legendre7xx(x) * Legendre10(y);
    }

    static double leg_quad_l7_l10xy(double x, double y)
    {
      return Legendre7x(x) * Legendre10x(y);
    }

    static double leg_quad_l7_l10yy(double x, double y)
    {
      return  Legendre7(x) * Legendre10xx(y);
    }

    static double leg_quad_l8_l0(double x, double y)
    {
      return Legendre8(x) * Legendre0(y);
    }

    static double leg_quad_l8_l0x(double x, double y)
    {
      return Legendre8x(x) * Legendre0(y);
    }

    static double leg_quad_l8_l0y(double x, double y)
    {
      return Legendre8(x) * Legendre0x(y);
    }

    static double leg_quad_l8_l0xx(double x, double y)
    {
      return Legendre8xx(x) * Legendre0(y);
    }

    static double leg_quad_l8_l0xy(double x, double y)
    {
      return Legendre8x(x) * Legendre0x(y);
    }

    static double leg_quad_l8_l0yy(double x, double y)
    {
      return  Legendre8(x) * Legendre0xx(y);
    }

    static double leg_quad_l8_l1(double x, double y)
    {
      return Legendre8(x) * Legendre1(y);
    }

    static double leg_quad_l8_l1x(double x, double y)
    {
      return Legendre8x(x) * Legendre1(y);
    }

    static double leg_quad_l8_l1y(double x, double y)
    {
      return Legendre8(x) * Legendre1x(y);
    }

    static double leg_quad_l8_l1xx(double x, double y)
    {
      return Legendre8xx(x) * Legendre1(y);
    }

    static double leg_quad_l8_l1xy(double x, double y)
    {
      return Legendre8x(x) * Legendre1x(y);
    }

    static double leg_quad_l8_l1yy(double x, double y)
    {
      return  Legendre8(x) * Legendre1xx(y);
    }

    static double leg_quad_l8_l2(double x, double y)
    {
      return Legendre8(x) * Legendre2(y);
    }

    static double leg_quad_l8_l2x(double x, double y)
    {
      return Legendre8x(x) * Legendre2(y);
    }

    static double leg_quad_l8_l2y(double x, double y)
    {
      return Legendre8(x) * Legendre2x(y);
    }

    static double leg_quad_l8_l2xx(double x, double y)
    {
      return Legendre8xx(x) * Legendre2(y);
    }

    static double leg_quad_l8_l2xy(double x, double y)
    {
      return Legendre8x(x) * Legendre2x(y);
    }

    static double leg_quad_l8_l2yy(double x, double y)
    {
      return  Legendre8(x) * Legendre2xx(y);
    }

    static double leg_quad_l8_l3(double x, double y)
    {
      return Legendre8(x) * Legendre3(y);
    }

    static double leg_quad_l8_l3x(double x, double y)
    {
      return Legendre8x(x) * Legendre3(y);
    }

    static double leg_quad_l8_l3y(double x, double y)
    {
      return Legendre8(x) * Legendre3x(y);
    }

    static double leg_quad_l8_l3xx(double x, double y)
    {
      return Legendre8xx(x) * Legendre3(y);
    }

    static double leg_quad_l8_l3xy(double x, double y)
    {
      return Legendre8x(x) * Legendre3x(y);
    }

    static double leg_quad_l8_l3yy(double x, double y)
    {
      return  Legendre8(x) * Legendre3xx(y);
    }

    static double leg_quad_l8_l4(double x, double y)
    {
      return Legendre8(x) * Legendre4(y);
    }

    static double leg_quad_l8_l4x(double x, double y)
    {
      return Legendre8x(x) * Legendre4(y);
    }

    static double leg_quad_l8_l4y(double x, double y)
    {
      return Legendre8(x) * Legendre4x(y);
    }

    static double leg_quad_l8_l4xx(double x, double y)
    {
      return Legendre8xx(x) * Legendre4(y);
    }

    static double leg_quad_l8_l4xy(double x, double y)
    {
      return Legendre8x(x) * Legendre4x(y);
    }

    static double leg_quad_l8_l4yy(double x, double y)
    {
      return  Legendre8(x) * Legendre4xx(y);
    }

    static double leg_quad_l8_l5(double x, double y)
    {
      return Legendre8(x) * Legendre5(y);
    }

    static double leg_quad_l8_l5x(double x, double y)
    {
      return Legendre8x(x) * Legendre5(y);
    }

    static double leg_quad_l8_l5y(double x, double y)
    {
      return Legendre8(x) * Legendre5x(y);
    }

    static double leg_quad_l8_l5xx(double x, double y)
    {
      return Legendre8xx(x) * Legendre5(y);
    }

    static double leg_quad_l8_l5xy(double x, double y)
    {
      return Legendre8x(x) * Legendre5x(y);
    }

    static double leg_quad_l8_l5yy(double x, double y)
    {
      return  Legendre8(x) * Legendre5xx(y);
    }

    static double leg_quad_l8_l6(double x, double y)
    {
      return Legendre8(x) * Legendre6(y);
    }

    static double leg_quad_l8_l6x(double x, double y)
    {
      return Legendre8x(x) * Legendre6(y);
    }

    static double leg_quad_l8_l6y(double x, double y)
    {
      return Legendre8(x) * Legendre6x(y);
    }

    static double leg_quad_l8_l6xx(double x, double y)
    {
      return Legendre8xx(x) * Legendre6(y);
    }

    static double leg_quad_l8_l6xy(double x, double y)
    {
      return Legendre8x(x) * Legendre6x(y);
    }

    static double leg_quad_l8_l6yy(double x, double y)
    {
      return  Legendre8(x) * Legendre6xx(y);
    }

    static double leg_quad_l8_l7(double x, double y)
    {
      return Legendre8(x) * Legendre7(y);
    }

    static double leg_quad_l8_l7x(double x, double y)
    {
      return Legendre8x(x) * Legendre7(y);
    }

    static double leg_quad_l8_l7y(double x, double y)
    {
      return Legendre8(x) * Legendre7x(y);
    }

    static double leg_quad_l8_l7xx(double x, double y)
    {
      return Legendre8xx(x) * Legendre7(y);
    }

    static double leg_quad_l8_l7xy(double x, double y)
    {
      return Legendre8x(x) * Legendre7x(y);
    }

    static double leg_quad_l8_l7yy(double x, double y)
    {
      return  Legendre8(x) * Legendre7xx(y);
    }

    static double leg_quad_l8_l8(double x, double y)
    {
      return Legendre8(x) * Legendre8(y);
    }

    static double leg_quad_l8_l8x(double x, double y)
    {
      return Legendre8x(x) * Legendre8(y);
    }

    static double leg_quad_l8_l8y(double x, double y)
    {
      return Legendre8(x) * Legendre8x(y);
    }

    static double leg_quad_l8_l8xx(double x, double y)
    {
      return Legendre8xx(x) * Legendre8(y);
    }

    static double leg_quad_l8_l8xy(double x, double y)
    {
      return Legendre8x(x) * Legendre8x(y);
    }

    static double leg_quad_l8_l8yy(double x, double y)
    {
      return  Legendre8(x) * Legendre8xx(y);
    }

    static double leg_quad_l8_l9(double x, double y)
    {
      return Legendre8(x) * Legendre9(y);
    }

    static double leg_quad_l8_l9x(double x, double y)
    {
      return Legendre8x(x) * Legendre9(y);
    }

    static double leg_quad_l8_l9y(double x, double y)
    {
      return Legendre8(x) * Legendre9x(y);
    }

    static double leg_quad_l8_l9xx(double x, double y)
    {
      return Legendre8xx(x) * Legendre9(y);
    }

    static double leg_quad_l8_l9xy(double x, double y)
    {
      return Legendre8x(x) * Legendre9x(y);
    }

    static double leg_quad_l8_l9yy(double x, double y)
    {
      return  Legendre8(x) * Legendre9xx(y);
    }

    static double leg_quad_l8_l10(double x, double y)
    {
      return Legendre8(x) * Legendre10(y);
    }

    static double leg_quad_l8_l10x(double x, double y)
    {
      return Legendre8x(x) * Legendre10(y);
    }

    static double leg_quad_l8_l10y(double x, double y)
    {
      return Legendre8(x) * Legendre10x(y);
    }

    static double leg_quad_l8_l10xx(double x, double y)
    {
      return Legendre8xx(x) * Legendre10(y);
    }

    static double leg_quad_l8_l10xy(double x, double y)
    {
      return Legendre8x(x) * Legendre10x(y);
    }

    static double leg_quad_l8_l10yy(double x, double y)
    {
      return  Legendre8(x) * Legendre10xx(y);
    }

    static double leg_quad_l9_l0(double x, double y)
    {
      return Legendre9(x) * Legendre0(y);
    }

    static double leg_quad_l9_l0x(double x, double y)
    {
      return Legendre9x(x) * Legendre0(y);
    }

    static double leg_quad_l9_l0y(double x, double y)
    {
      return Legendre9(x) * Legendre0x(y);
    }

    static double leg_quad_l9_l0xx(double x, double y)
    {
      return Legendre9xx(x) * Legendre0(y);
    }

    static double leg_quad_l9_l0xy(double x, double y)
    {
      return Legendre9x(x) * Legendre0x(y);
    }

    static double leg_quad_l9_l0yy(double x, double y)
    {
      return  Legendre9(x) * Legendre0xx(y);
    }

    static double leg_quad_l9_l1(double x, double y)
    {
      return Legendre9(x) * Legendre1(y);
    }

    static double leg_quad_l9_l1x(double x, double y)
    {
      return Legendre9x(x) * Legendre1(y);
    }

    static double leg_quad_l9_l1y(double x, double y)
    {
      return Legendre9(x) * Legendre1x(y);
    }

    static double leg_quad_l9_l1xx(double x, double y)
    {
      return Legendre9xx(x) * Legendre1(y);
    }

    static double leg_quad_l9_l1xy(double x, double y)
    {
      return Legendre9x(x) * Legendre1x(y);
    }

    static double leg_quad_l9_l1yy(double x, double y)
    {
      return  Legendre9(x) * Legendre1xx(y);
    }

    static double leg_quad_l9_l2(double x, double y)
    {
      return Legendre9(x) * Legendre2(y);
    }

    static double leg_quad_l9_l2x(double x, double y)
    {
      return Legendre9x(x) * Legendre2(y);
    }

    static double leg_quad_l9_l2y(double x, double y)
    {
      return Legendre9(x) * Legendre2x(y);
    }

    static double leg_quad_l9_l2xx(double x, double y)
    {
      return Legendre9xx(x) * Legendre2(y);
    }

    static double leg_quad_l9_l2xy(double x, double y)
    {
      return Legendre9x(x) * Legendre2x(y);
    }

    static double leg_quad_l9_l2yy(double x, double y)
    {
      return  Legendre9(x) * Legendre2xx(y);
    }

    static double leg_quad_l9_l3(double x, double y)
    {
      return Legendre9(x) * Legendre3(y);
    }

    static double leg_quad_l9_l3x(double x, double y)
    {
      return Legendre9x(x) * Legendre3(y);
    }

    static double leg_quad_l9_l3y(double x, double y)
    {
      return Legendre9(x) * Legendre3x(y);
    }

    static double leg_quad_l9_l3xx(double x, double y)
    {
      return Legendre9xx(x) * Legendre3(y);
    }

    static double leg_quad_l9_l3xy(double x, double y)
    {
      return Legendre9x(x) * Legendre3x(y);
    }

    static double leg_quad_l9_l3yy(double x, double y)
    {
      return  Legendre9(x) * Legendre3xx(y);
    }

    static double leg_quad_l9_l4(double x, double y)
    {
      return Legendre9(x) * Legendre4(y);
    }

    static double leg_quad_l9_l4x(double x, double y)
    {
      return Legendre9x(x) * Legendre4(y);
    }

    static double leg_quad_l9_l4y(double x, double y)
    {
      return Legendre9(x) * Legendre4x(y);
    }

    static double leg_quad_l9_l4xx(double x, double y)
    {
      return Legendre9xx(x) * Legendre4(y);
    }

    static double leg_quad_l9_l4xy(double x, double y)
    {
      return Legendre9x(x) * Legendre4x(y);
    }

    static double leg_quad_l9_l4yy(double x, double y)
    {
      return  Legendre9(x) * Legendre4xx(y);
    }

    static double leg_quad_l9_l5(double x, double y)
    {
      return Legendre9(x) * Legendre5(y);
    }

    static double leg_quad_l9_l5x(double x, double y)
    {
      return Legendre9x(x) * Legendre5(y);
    }

    static double leg_quad_l9_l5y(double x, double y)
    {
      return Legendre9(x) * Legendre5x(y);
    }

    static double leg_quad_l9_l5xx(double x, double y)
    {
      return Legendre9xx(x) * Legendre5(y);
    }

    static double leg_quad_l9_l5xy(double x, double y)
    {
      return Legendre9x(x) * Legendre5x(y);
    }

    static double leg_quad_l9_l5yy(double x, double y)
    {
      return  Legendre9(x) * Legendre5xx(y);
    }

    static double leg_quad_l9_l6(double x, double y)
    {
      return Legendre9(x) * Legendre6(y);
    }

    static double leg_quad_l9_l6x(double x, double y)
    {
      return Legendre9x(x) * Legendre6(y);
    }

    static double leg_quad_l9_l6y(double x, double y)
    {
      return Legendre9(x) * Legendre6x(y);
    }

    static double leg_quad_l9_l6xx(double x, double y)
    {
      return Legendre9xx(x) * Legendre6(y);
    }

    static double leg_quad_l9_l6xy(double x, double y)
    {
      return Legendre9x(x) * Legendre6x(y);
    }

    static double leg_quad_l9_l6yy(double x, double y)
    {
      return  Legendre9(x) * Legendre6xx(y);
    }

    static double leg_quad_l9_l7(double x, double y)
    {
      return Legendre9(x) * Legendre7(y);
    }

    static double leg_quad_l9_l7x(double x, double y)
    {
      return Legendre9x(x) * Legendre7(y);
    }

    static double leg_quad_l9_l7y(double x, double y)
    {
      return Legendre9(x) * Legendre7x(y);
    }

    static double leg_quad_l9_l7xx(double x, double y)
    {
      return Legendre9xx(x) * Legendre7(y);
    }

    static double leg_quad_l9_l7xy(double x, double y)
    {
      return Legendre9x(x) * Legendre7x(y);
    }

    static double leg_quad_l9_l7yy(double x, double y)
    {
      return  Legendre9(x) * Legendre7xx(y);
    }

    static double leg_quad_l9_l8(double x, double y)
    {
      return Legendre9(x) * Legendre8(y);
    }

    static double leg_quad_l9_l8x(double x, double y)
    {
      return Legendre9x(x) * Legendre8(y);
    }

    static double leg_quad_l9_l8y(double x, double y)
    {
      return Legendre9(x) * Legendre8x(y);
    }

    static double leg_quad_l9_l8xx(double x, double y)
    {
      return Legendre9xx(x) * Legendre8(y);
    }

    static double leg_quad_l9_l8xy(double x, double y)
    {
      return Legendre9x(x) * Legendre8x(y);
    }

    static double leg_quad_l9_l8yy(double x, double y)
    {
      return  Legendre9(x) * Legendre8xx(y);
    }

    static double leg_quad_l9_l9(double x, double y)
    {
      return Legendre9(x) * Legendre9(y);
    }

    static double leg_quad_l9_l9x(double x, double y)
    {
      return Legendre9x(x) * Legendre9(y);
    }

    static double leg_quad_l9_l9y(double x, double y)
    {
      return Legendre9(x) * Legendre9x(y);
    }

    static double leg_quad_l9_l9xx(double x, double y)
    {
      return Legendre9xx(x) * Legendre9(y);
    }

    static double leg_quad_l9_l9xy(double x, double y)
    {
      return Legendre9x(x) * Legendre9x(y);
    }

    static double leg_quad_l9_l9yy(double x, double y)
    {
      return  Legendre9(x) * Legendre9xx(y);
    }

    static double leg_quad_l9_l10(double x, double y)
    {
      return Legendre9(x) * Legendre10(y);
    }

    static double leg_quad_l9_l10x(double x, double y)
    {
      return Legendre9x(x) * Legendre10(y);
    }

    static double leg_quad_l9_l10y(double x, double y)
    {
      return Legendre9(x) * Legendre10x(y);
    }

    static double leg_quad_l9_l10xx(double x, double y)
    {
      return Legendre9xx(x) * Legendre10(y);
    }

    static double leg_quad_l9_l10xy(double x, double y)
    {
      return Legendre9x(x) * Legendre10x(y);
    }

    static double leg_quad_l9_l10yy(double x, double y)
    {
      return  Legendre9(x) * Legendre10xx(y);
    }

    static double leg_quad_l10_l0(double x, double y)
    {
      return Legendre10(x) * Legendre0(y);
    }

    static double leg_quad_l10_l0x(double x, double y)
    {
      return Legendre10x(x) * Legendre0(y);
    }

    static double leg_quad_l10_l0y(double x, double y)
    {
      return Legendre10(x) * Legendre0x(y);
    }

    static double leg_quad_l10_l0xx(double x, double y)
    {
      return Legendre10xx(x) * Legendre0(y);
    }

    static double leg_quad_l10_l0xy(double x, double y)
    {
      return Legendre10x(x) * Legendre0x(y);
    }

    static double leg_quad_l10_l0yy(double x, double y)
    {
      return  Legendre10(x) * Legendre0xx(y);
    }

    static double leg_quad_l10_l1(double x, double y)
    {
      return Legendre10(x) * Legendre1(y);
    }

    static double leg_quad_l10_l1x(double x, double y)
    {
      return Legendre10x(x) * Legendre1(y);
    }

    static double leg_quad_l10_l1y(double x, double y)
    {
      return Legendre10(x) * Legendre1x(y);
    }

    static double leg_quad_l10_l1xx(double x, double y)
    {
      return Legendre10xx(x) * Legendre1(y);
    }

    static double leg_quad_l10_l1xy(double x, double y)
    {
      return Legendre10x(x) * Legendre1x(y);
    }

    static double leg_quad_l10_l1yy(double x, double y)
    {
      return  Legendre10(x) * Legendre1xx(y);
    }

    static double leg_quad_l10_l2(double x, double y)
    {
      return Legendre10(x) * Legendre2(y);
    }

    static double leg_quad_l10_l2x(double x, double y)
    {
      return Legendre10x(x) * Legendre2(y);
    }

    static double leg_quad_l10_l2y(double x, double y)
    {
      return Legendre10(x) * Legendre2x(y);
    }

    static double leg_quad_l10_l2xx(double x, double y)
    {
      return Legendre10xx(x) * Legendre2(y);
    }

    static double leg_quad_l10_l2xy(double x, double y)
    {
      return Legendre10x(x) * Legendre2x(y);
    }

    static double leg_quad_l10_l2yy(double x, double y)
    {
      return  Legendre10(x) * Legendre2xx(y);
    }

    static double leg_quad_l10_l3(double x, double y)
    {
      return Legendre10(x) * Legendre3(y);
    }

    static double leg_quad_l10_l3x(double x, double y)
    {
      return Legendre10x(x) * Legendre3(y);
    }

    static double leg_quad_l10_l3y(double x, double y)
    {
      return Legendre10(x) * Legendre3x(y);
    }

    static double leg_quad_l10_l3xx(double x, double y)
    {
      return Legendre10xx(x) * Legendre3(y);
    }

    static double leg_quad_l10_l3xy(double x, double y)
    {
      return Legendre10x(x) * Legendre3x(y);
    }

    static double leg_quad_l10_l3yy(double x, double y)
    {
      return  Legendre10(x) * Legendre3xx(y);
    }

    static double leg_quad_l10_l4(double x, double y)
    {
      return Legendre10(x) * Legendre4(y);
    }

    static double leg_quad_l10_l4x(double x, double y)
    {
      return Legendre10x(x) * Legendre4(y);
    }

    static double leg_quad_l10_l4y(double x, double y)
    {
      return Legendre10(x) * Legendre4x(y);
    }

    static double leg_quad_l10_l4xx(double x, double y)
    {
      return Legendre10xx(x) * Legendre4(y);
    }

    static double leg_quad_l10_l4xy(double x, double y)
    {
      return Legendre10x(x) * Legendre4x(y);
    }

    static double leg_quad_l10_l4yy(double x, double y)
    {
      return  Legendre10(x) * Legendre4xx(y);
    }

    static double leg_quad_l10_l5(double x, double y)
    {
      return Legendre10(x) * Legendre5(y);
    }

    static double leg_quad_l10_l5x(double x, double y)
    {
      return Legendre10x(x) * Legendre5(y);
    }

    static double leg_quad_l10_l5y(double x, double y)
    {
      return Legendre10(x) * Legendre5x(y);
    }

    static double leg_quad_l10_l5xx(double x, double y)
    {
      return Legendre10xx(x) * Legendre5(y);
    }

    static double leg_quad_l10_l5xy(double x, double y)
    {
      return Legendre10x(x) * Legendre5x(y);
    }

    static double leg_quad_l10_l5yy(double x, double y)
    {
      return  Legendre10(x) * Legendre5xx(y);
    }

    static double leg_quad_l10_l6(double x, double y)
    {
      return Legendre10(x) * Legendre6(y);
    }

    static double leg_quad_l10_l6x(double x, double y)
    {
      return Legendre10x(x) * Legendre6(y);
    }

    static double leg_quad_l10_l6y(double x, double y)
    {
      return Legendre10(x) * Legendre6x(y);
    }

    static double leg_quad_l10_l6xx(double x, double y)
    {
      return Legendre10xx(x) * Legendre6(y);
    }

    static double leg_quad_l10_l6xy(double x, double y)
    {
      return Legendre10x(x) * Legendre6x(y);
    }

    static double leg_quad_l10_l6yy(double x, double y)
    {
      return  Legendre10(x) * Legendre6xx(y);
    }

    static double leg_quad_l10_l7(double x, double y)
    {
      return Legendre10(x) * Legendre7(y);
    }

    static double leg_quad_l10_l7x(double x, double y)
    {
      return Legendre10x(x) * Legendre7(y);
    }

    static double leg_quad_l10_l7y(double x, double y)
    {
      return Legendre10(x) * Legendre7x(y);
    }

    static double leg_quad_l10_l7xx(double x, double y)
    {
      return Legendre10xx(x) * Legendre7(y);
    }

    static double leg_quad_l10_l7xy(double x, double y)
    {
      return Legendre10x(x) * Legendre7x(y);
    }

    static double leg_quad_l10_l7yy(double x, double y)
    {
      return  Legendre10(x) * Legendre7xx(y);
    }

    static double leg_quad_l10_l8(double x, double y)
    {
      return Legendre10(x) * Legendre8(y);
    }

    static double leg_quad_l10_l8x(double x, double y)
    {
      return Legendre10x(x) * Legendre8(y);
    }

    static double leg_quad_l10_l8y(double x, double y)
    {
      return Legendre10(x) * Legendre8x(y);
    }

    static double leg_quad_l10_l8xx(double x, double y)
    {
      return Legendre10xx(x) * Legendre8(y);
    }

    static double leg_quad_l10_l8xy(double x, double y)
    {
      return Legendre10x(x) * Legendre8x(y);
    }

    static double leg_quad_l10_l8yy(double x, double y)
    {
      return  Legendre10(x) * Legendre8xx(y);
    }

    static double leg_quad_l10_l9(double x, double y)
    {
      return Legendre10(x) * Legendre9(y);
    }

    static double leg_quad_l10_l9x(double x, double y)
    {
      return Legendre10x(x) * Legendre9(y);
    }

    static double leg_quad_l10_l9y(double x, double y)
    {
      return Legendre10(x) * Legendre9x(y);
    }

    static double leg_quad_l10_l9xx(double x, double y)
    {
      return Legendre10xx(x) * Legendre9(y);
    }

    static double leg_quad_l10_l9xy(double x, double y)
    {
      return Legendre10x(x) * Legendre9x(y);
    }

    static double leg_quad_l10_l9yy(double x, double y)
    {
      return  Legendre10(x) * Legendre9xx(y);
    }

    static double leg_quad_l10_l10(double x, double y)
    {
      return Legendre10(x) * Legendre10(y);
    }

    static double leg_quad_l10_l10x(double x, double y)
    {
      return Legendre10x(x) * Legendre10(y);
    }

    static double leg_quad_l10_l10y(double x, double y)
    {
      return Legendre10(x) * Legendre10x(y);
    }

    static double leg_quad_l10_l10xx(double x, double y)
    {
      return Legendre10xx(x) * Legendre10(y);
    }

    static double leg_quad_l10_l10xy(double x, double y)
    {
      return Legendre10x(x) * Legendre10x(y);
    }

    static double leg_quad_l10_l10yy(double x, double y)
    {
      return  Legendre10(x) * Legendre10xx(y);
    }

    static Shapeset::shape_fn_t leg_quad_fn[] =
    {
      leg_quad_l0_l0,   leg_quad_l0_l1,   leg_quad_l0_l2,   leg_quad_l0_l3,   leg_quad_l0_l4,
      leg_quad_l0_l5,   leg_quad_l0_l6,   leg_quad_l0_l7,   leg_quad_l0_l8,   leg_quad_l0_l9,
      leg_quad_l0_l10,   leg_quad_l1_l0,   leg_quad_l1_l1,   leg_quad_l1_l2,   leg_quad_l1_l3,
      leg_quad_l1_l4,   leg_quad_l1_l5,   leg_quad_l1_l6,   leg_quad_l1_l7,   leg_quad_l1_l8,
      leg_quad_l1_l9,   leg_quad_l1_l10,   leg_quad_l2_l0,   leg_quad_l2_l1,   leg_quad_l2_l2,
      leg_quad_l2_l3,   leg_quad_l2_l4,   leg_quad_l2_l5,   leg_quad_l2_l6,   leg_quad_l2_l7,
      leg_quad_l2_l8,   leg_quad_l2_l9,   leg_quad_l2_l10,   leg_quad_l3_l0,   leg_quad_l3_l1,
      leg_quad_l3_l2,   leg_quad_l3_l3,   leg_quad_l3_l4,   leg_quad_l3_l5,   leg_quad_l3_l6,
      leg_quad_l3_l7,   leg_quad_l3_l8,   leg_quad_l3_l9,   leg_quad_l3_l10,   leg_quad_l4_l0,
      leg_quad_l4_l1,   leg_quad_l4_l2,   leg_quad_l4_l3,   leg_quad_l4_l4,   leg_quad_l4_l5,
      leg_quad_l4_l6,   leg_quad_l4_l7,   leg_quad_l4_l8,   leg_quad_l4_l9,   leg_quad_l4_l10,
      leg_quad_l5_l0,   leg_quad_l5_l1,   leg_quad_l5_l2,   leg_quad_l5_l3,   leg_quad_l5_l4,
      leg_quad_l5_l5,   leg_quad_l5_l6,   leg_quad_l5_l7,   leg_quad_l5_l8,   leg_quad_l5_l9,
      leg_quad_l5_l10,   leg_quad_l6_l0,   leg_quad_l6_l1,   leg_quad_l6_l2,   leg_quad_l6_l3,
      leg_quad_l6_l4,   leg_quad_l6_l5,   leg_quad_l6_l6,   leg_quad_l6_l7,   leg_quad_l6_l8,
      leg_quad_l6_l9,   leg_quad_l6_l10,   leg_quad_l7_l0,   leg_quad_l7_l1,   leg_quad_l7_l2,
      leg_quad_l7_l3,   leg_quad_l7_l4,   leg_quad_l7_l5,   leg_quad_l7_l6,   leg_quad_l7_l7,
      leg_quad_l7_l8,   leg_quad_l7_l9,   leg_quad_l7_l10,   leg_quad_l8_l0,   leg_quad_l8_l1,
      leg_quad_l8_l2,   leg_quad_l8_l3,   leg_quad_l8_l4,   leg_quad_l8_l5,   leg_quad_l8_l6,
      leg_quad_l8_l7,   leg_quad_l8_l8,   leg_quad_l8_l9,   leg_quad_l8_l10,   leg_quad_l9_l0,
      leg_quad_l9_l1,   leg_quad_l9_l2,   leg_quad_l9_l3,   leg_quad_l9_l4,   leg_quad_l9_l5,
      leg_quad_l9_l6,   leg_quad_l9_l7,   leg_quad_l9_l8,   leg_quad_l9_l9,   leg_quad_l9_l10,
      leg_quad_l10_l0,   leg_quad_l10_l1,   leg_quad_l10_l2,   leg_quad_l10_l3,   leg_quad_l10_l4,
      leg_quad_l10_l5,   leg_quad_l10_l6,   leg_quad_l10_l7,   leg_quad_l10_l8,   leg_quad_l10_l9,
      leg_quad_l10_l10,
    };
    static Shapeset::shape_fn_t leg_quad_fn_dx[] =
    {
      leg_quad_l0_l0x,   leg_quad_l0_l1x,   leg_quad_l0_l2x,   leg_quad_l0_l3x,   leg_quad_l0_l4x,
      leg_quad_l0_l5x,   leg_quad_l0_l6x,   leg_quad_l0_l7x,   leg_quad_l0_l8x,   leg_quad_l0_l9x,
      leg_quad_l0_l10x,   leg_quad_l1_l0x,   leg_quad_l1_l1x,   leg_quad_l1_l2x,   leg_quad_l1_l3x,
      leg_quad_l1_l4x,   leg_quad_l1_l5x,   leg_quad_l1_l6x,   leg_quad_l1_l7x,   leg_quad_l1_l8x,
      leg_quad_l1_l9x,   leg_quad_l1_l10x,   leg_quad_l2_l0x,   leg_quad_l2_l1x,   leg_quad_l2_l2x,
      leg_quad_l2_l3x,   leg_quad_l2_l4x,   leg_quad_l2_l5x,   leg_quad_l2_l6x,   leg_quad_l2_l7x,
      leg_quad_l2_l8x,   leg_quad_l2_l9x,   leg_quad_l2_l10x,   leg_quad_l3_l0x,   leg_quad_l3_l1x,
      leg_quad_l3_l2x,   leg_quad_l3_l3x,   leg_quad_l3_l4x,   leg_quad_l3_l5x,   leg_quad_l3_l6x,
      leg_quad_l3_l7x,   leg_quad_l3_l8x,   leg_quad_l3_l9x,   leg_quad_l3_l10x,   leg_quad_l4_l0x,
      leg_quad_l4_l1x,   leg_quad_l4_l2x,   leg_quad_l4_l3x,   leg_quad_l4_l4x,   leg_quad_l4_l5x,
      leg_quad_l4_l6x,   leg_quad_l4_l7x,   leg_quad_l4_l8x,   leg_quad_l4_l9x,   leg_quad_l4_l10x,
      leg_quad_l5_l0x,   leg_quad_l5_l1x,   leg_quad_l5_l2x,   leg_quad_l5_l3x,   leg_quad_l5_l4x,
      leg_quad_l5_l5x,   leg_quad_l5_l6x,   leg_quad_l5_l7x,   leg_quad_l5_l8x,   leg_quad_l5_l9x,
      leg_quad_l5_l10x,   leg_quad_l6_l0x,   leg_quad_l6_l1x,   leg_quad_l6_l2x,   leg_quad_l6_l3x,
      leg_quad_l6_l4x,   leg_quad_l6_l5x,   leg_quad_l6_l6x,   leg_quad_l6_l7x,   leg_quad_l6_l8x,
      leg_quad_l6_l9x,   leg_quad_l6_l10x,   leg_quad_l7_l0x,   leg_quad_l7_l1x,   leg_quad_l7_l2x,
      leg_quad_l7_l3x,   leg_quad_l7_l4x,   leg_quad_l7_l5x,   leg_quad_l7_l6x,   leg_quad_l7_l7x,
      leg_quad_l7_l8x,   leg_quad_l7_l9x,   leg_quad_l7_l10x,   leg_quad_l8_l0x,   leg_quad_l8_l1x,
      leg_quad_l8_l2x,   leg_quad_l8_l3x,   leg_quad_l8_l4x,   leg_quad_l8_l5x,   leg_quad_l8_l6x,
      leg_quad_l8_l7x,   leg_quad_l8_l8x,   leg_quad_l8_l9x,   leg_quad_l8_l10x,   leg_quad_l9_l0x,
      leg_quad_l9_l1x,   leg_quad_l9_l2x,   leg_quad_l9_l3x,   leg_quad_l9_l4x,   leg_quad_l9_l5x,
      leg_quad_l9_l6x,   leg_quad_l9_l7x,   leg_quad_l9_l8x,   leg_quad_l9_l9x,   leg_quad_l9_l10x,
      leg_quad_l10_l0x,   leg_quad_l10_l1x,   leg_quad_l10_l2x,   leg_quad_l10_l3x,   leg_quad_l10_l4x,
      leg_quad_l10_l5x,   leg_quad_l10_l6x,   leg_quad_l10_l7x,   leg_quad_l10_l8x,   leg_quad_l10_l9x,
      leg_quad_l10_l10x,
    };
    static Shapeset::shape_fn_t leg_quad_fn_dy[] =
    {
      leg_quad_l0_l0y,   leg_quad_l0_l1y,   leg_quad_l0_l2y,   leg_quad_l0_l3y,   leg_quad_l0_l4y,
      leg_quad_l0_l5y,   leg_quad_l0_l6y,   leg_quad_l0_l7y,   leg_quad_l0_l8y,   leg_quad_l0_l9y,
      leg_quad_l0_l10y,   leg_quad_l1_l0y,   leg_quad_l1_l1y,   leg_quad_l1_l2y,   leg_quad_l1_l3y,
      leg_quad_l1_l4y,   leg_quad_l1_l5y,   leg_quad_l1_l6y,   leg_quad_l1_l7y,   leg_quad_l1_l8y,
      leg_quad_l1_l9y,   leg_quad_l1_l10y,   leg_quad_l2_l0y,   leg_quad_l2_l1y,   leg_quad_l2_l2y,
      leg_quad_l2_l3y,   leg_quad_l2_l4y,   leg_quad_l2_l5y,   leg_quad_l2_l6y,   leg_quad_l2_l7y,
      leg_quad_l2_l8y,   leg_quad_l2_l9y,   leg_quad_l2_l10y,   leg_quad_l3_l0y,   leg_quad_l3_l1y,
      leg_quad_l3_l2y,   leg_quad_l3_l3y,   leg_quad_l3_l4y,   leg_quad_l3_l5y,   leg_quad_l3_l6y,
      leg_quad_l3_l7y,   leg_quad_l3_l8y,   leg_quad_l3_l9y,   leg_quad_l3_l10y,   leg_quad_l4_l0y,
      leg_quad_l4_l1y,   leg_quad_l4_l2y,   leg_quad_l4_l3y,   leg_quad_l4_l4y,   leg_quad_l4_l5y,
      leg_quad_l4_l6y,   leg_quad_l4_l7y,   leg_quad_l4_l8y,   leg_quad_l4_l9y,   leg_quad_l4_l10y,
      leg_quad_l5_l0y,   leg_quad_l5_l1y,   leg_quad_l5_l2y,   leg_quad_l5_l3y,   leg_quad_l5_l4y,
      leg_quad_l5_l5y,   leg_quad_l5_l6y,   leg_quad_l5_l7y,   leg_quad_l5_l8y,   leg_quad_l5_l9y,
      leg_quad_l5_l10y,   leg_quad_l6_l0y,   leg_quad_l6_l1y,   leg_quad_l6_l2y,   leg_quad_l6_l3y,
      leg_quad_l6_l4y,   leg_quad_l6_l5y,   leg_quad_l6_l6y,   leg_quad_l6_l7y,   leg_quad_l6_l8y,
      leg_quad_l6_l9y,   leg_quad_l6_l10y,   leg_quad_l7_l0y,   leg_quad_l7_l1y,   leg_quad_l7_l2y,
      leg_quad_l7_l3y,   leg_quad_l7_l4y,   leg_quad_l7_l5y,   leg_quad_l7_l6y,   leg_quad_l7_l7y,
      leg_quad_l7_l8y,   leg_quad_l7_l9y,   leg_quad_l7_l10y,   leg_quad_l8_l0y,   leg_quad_l8_l1y,
      leg_quad_l8_l2y,   leg_quad_l8_l3y,   leg_quad_l8_l4y,   leg_quad_l8_l5y,   leg_quad_l8_l6y,
      leg_quad_l8_l7y,   leg_quad_l8_l8y,   leg_quad_l8_l9y,   leg_quad_l8_l10y,   leg_quad_l9_l0y,
      leg_quad_l9_l1y,   leg_quad_l9_l2y,   leg_quad_l9_l3y,   leg_quad_l9_l4y,   leg_quad_l9_l5y,
      leg_quad_l9_l6y,   leg_quad_l9_l7y,   leg_quad_l9_l8y,   leg_quad_l9_l9y,   leg_quad_l9_l10y,
      leg_quad_l10_l0y,   leg_quad_l10_l1y,   leg_quad_l10_l2y,   leg_quad_l10_l3y,   leg_quad_l10_l4y,
      leg_quad_l10_l5y,   leg_quad_l10_l6y,   leg_quad_l10_l7y,   leg_quad_l10_l8y,   leg_quad_l10_l9y,
      leg_quad_l10_l10y,
    };
    static Shapeset::shape_fn_t leg_quad_fn_dxx[] =
    {
      leg_quad_l0_l0xx,   leg_quad_l0_l1xx,   leg_quad_l0_l2xx,   leg_quad_l0_l3xx,   leg_quad_l0_l4xx,
      leg_quad_l0_l5xx,   leg_quad_l0_l6xx,   leg_quad_l0_l7xx,   leg_quad_l0_l8xx,   leg_quad_l0_l9xx,
      leg_quad_l0_l10xx,   leg_quad_l1_l0xx,   leg_quad_l1_l1xx,   leg_quad_l1_l2xx,   leg_quad_l1_l3xx,
      leg_quad_l1_l4xx,   leg_quad_l1_l5xx,   leg_quad_l1_l6xx,   leg_quad_l1_l7xx,   leg_quad_l1_l8xx,
      leg_quad_l1_l9xx,   leg_quad_l1_l10xx,   leg_quad_l2_l0xx,   leg_quad_l2_l1xx,   leg_quad_l2_l2xx,
      leg_quad_l2_l3xx,   leg_quad_l2_l4xx,   leg_quad_l2_l5xx,   leg_quad_l2_l6xx,   leg_quad_l2_l7xx,
      leg_quad_l2_l8xx,   leg_quad_l2_l9xx,   leg_quad_l2_l10xx,   leg_quad_l3_l0xx,   leg_quad_l3_l1xx,
      leg_quad_l3_l2xx,   leg_quad_l3_l3xx,   leg_quad_l3_l4xx,   leg_quad_l3_l5xx,   leg_quad_l3_l6xx,
      leg_quad_l3_l7xx,   leg_quad_l3_l8xx,   leg_quad_l3_l9xx,   leg_quad_l3_l10xx,   leg_quad_l4_l0xx,
      leg_quad_l4_l1xx,   leg_quad_l4_l2xx,   leg_quad_l4_l3xx,   leg_quad_l4_l4xx,   leg_quad_l4_l5xx,
      leg_quad_l4_l6xx,   leg_quad_l4_l7xx,   leg_quad_l4_l8xx,   leg_quad_l4_l9xx,   leg_quad_l4_l10xx,
      leg_quad_l5_l0xx,   leg_quad_l5_l1xx,   leg_quad_l5_l2xx,   leg_quad_l5_l3xx,   leg_quad_l5_l4xx,
      leg_quad_l5_l5xx,   leg_quad_l5_l6xx,   leg_quad_l5_l7xx,   leg_quad_l5_l8xx,   leg_quad_l5_l9xx,
      leg_quad_l5_l10xx,   leg_quad_l6_l0xx,   leg_quad_l6_l1xx,   leg_quad_l6_l2xx,   leg_quad_l6_l3xx,
      leg_quad_l6_l4xx,   leg_quad_l6_l5xx,   leg_quad_l6_l6xx,   leg_quad_l6_l7xx,   leg_quad_l6_l8xx,
      leg_quad_l6_l9xx,   leg_quad_l6_l10xx,   leg_quad_l7_l0xx,   leg_quad_l7_l1xx,   leg_quad_l7_l2xx,
      leg_quad_l7_l3xx,   leg_quad_l7_l4xx,   leg_quad_l7_l5xx,   leg_quad_l7_l6xx,   leg_quad_l7_l7xx,
      leg_quad_l7_l8xx,   leg_quad_l7_l9xx,   leg_quad_l7_l10xx,   leg_quad_l8_l0xx,   leg_quad_l8_l1xx,
      leg_quad_l8_l2xx,   leg_quad_l8_l3xx,   leg_quad_l8_l4xx,   leg_quad_l8_l5xx,   leg_quad_l8_l6xx,
      leg_quad_l8_l7xx,   leg_quad_l8_l8xx,   leg_quad_l8_l9xx,   leg_quad_l8_l10xx,   leg_quad_l9_l0xx,
      leg_quad_l9_l1xx,   leg_quad_l9_l2xx,   leg_quad_l9_l3xx,   leg_quad_l9_l4xx,   leg_quad_l9_l5xx,
      leg_quad_l9_l6xx,   leg_quad_l9_l7xx,   leg_quad_l9_l8xx,   leg_quad_l9_l9xx,   leg_quad_l9_l10xx,
      leg_quad_l10_l0xx,   leg_quad_l10_l1xx,   leg_quad_l10_l2xx,   leg_quad_l10_l3xx,   leg_quad_l10_l4xx,
      leg_quad_l10_l5xx,   leg_quad_l10_l6xx,   leg_quad_l10_l7xx,   leg_quad_l10_l8xx,   leg_quad_l10_l9xx,
      leg_quad_l10_l10xx,
    };
    static Shapeset::shape_fn_t leg_quad_fn_dxy[] =
    {
      leg_quad_l0_l0xy,   leg_quad_l0_l1xy,   leg_quad_l0_l2xy,   leg_quad_l0_l3xy,   leg_quad_l0_l4xy,
      leg_quad_l0_l5xy,   leg_quad_l0_l6xy,   leg_quad_l0_l7xy,   leg_quad_l0_l8xy,   leg_quad_l0_l9xy,
      leg_quad_l0_l10xy,   leg_quad_l1_l0xy,   leg_quad_l1_l1xy,   leg_quad_l1_l2xy,   leg_quad_l1_l3xy,
      leg_quad_l1_l4xy,   leg_quad_l1_l5xy,   leg_quad_l1_l6xy,   leg_quad_l1_l7xy,   leg_quad_l1_l8xy,
      leg_quad_l1_l9xy,   leg_quad_l1_l10xy,   leg_quad_l2_l0xy,   leg_quad_l2_l1xy,   leg_quad_l2_l2xy,
      leg_quad_l2_l3xy,   leg_quad_l2_l4xy,   leg_quad_l2_l5xy,   leg_quad_l2_l6xy,   leg_quad_l2_l7xy,
      leg_quad_l2_l8xy,   leg_quad_l2_l9xy,   leg_quad_l2_l10xy,   leg_quad_l3_l0xy,   leg_quad_l3_l1xy,
      leg_quad_l3_l2xy,   leg_quad_l3_l3xy,   leg_quad_l3_l4xy,   leg_quad_l3_l5xy,   leg_quad_l3_l6xy,
      leg_quad_l3_l7xy,   leg_quad_l3_l8xy,   leg_quad_l3_l9xy,   leg_quad_l3_l10xy,   leg_quad_l4_l0xy,
      leg_quad_l4_l1xy,   leg_quad_l4_l2xy,   leg_quad_l4_l3xy,   leg_quad_l4_l4xy,   leg_quad_l4_l5xy,
      leg_quad_l4_l6xy,   leg_quad_l4_l7xy,   leg_quad_l4_l8xy,   leg_quad_l4_l9xy,   leg_quad_l4_l10xy,
      leg_quad_l5_l0xy,   leg_quad_l5_l1xy,   leg_quad_l5_l2xy,   leg_quad_l5_l3xy,   leg_quad_l5_l4xy,
      leg_quad_l5_l5xy,   leg_quad_l5_l6xy,   leg_quad_l5_l7xy,   leg_quad_l5_l8xy,   leg_quad_l5_l9xy,
      leg_quad_l5_l10xy,   leg_quad_l6_l0xy,   leg_quad_l6_l1xy,   leg_quad_l6_l2xy,   leg_quad_l6_l3xy,
      leg_quad_l6_l4xy,   leg_quad_l6_l5xy,   leg_quad_l6_l6xy,   leg_quad_l6_l7xy,   leg_quad_l6_l8xy,
      leg_quad_l6_l9xy,   leg_quad_l6_l10xy,   leg_quad_l7_l0xy,   leg_quad_l7_l1xy,   leg_quad_l7_l2xy,
      leg_quad_l7_l3xy,   leg_quad_l7_l4xy,   leg_quad_l7_l5xy,   leg_quad_l7_l6xy,   leg_quad_l7_l7xy,
      leg_quad_l7_l8xy,   leg_quad_l7_l9xy,   leg_quad_l7_l10xy,   leg_quad_l8_l0xy,   leg_quad_l8_l1xy,
      leg_quad_l8_l2xy,   leg_quad_l8_l3xy,   leg_quad_l8_l4xy,   leg_quad_l8_l5xy,   leg_quad_l8_l6xy,
      leg_quad_l8_l7xy,   leg_quad_l8_l8xy,   leg_quad_l8_l9xy,   leg_quad_l8_l10xy,   leg_quad_l9_l0xy,
      leg_quad_l9_l1xy,   leg_quad_l9_l2xy,   leg_quad_l9_l3xy,   leg_quad_l9_l4xy,   leg_quad_l9_l5xy,
      leg_quad_l9_l6xy,   leg_quad_l9_l7xy,   leg_quad_l9_l8xy,   leg_quad_l9_l9xy,   leg_quad_l9_l10xy,
      leg_quad_l10_l0xy,   leg_quad_l10_l1xy,   leg_quad_l10_l2xy,   leg_quad_l10_l3xy,   leg_quad_l10_l4xy,
      leg_quad_l10_l5xy,   leg_quad_l10_l6xy,   leg_quad_l10_l7xy,   leg_quad_l10_l8xy,   leg_quad_l10_l9xy,
      leg_quad_l10_l10xy,
    };
    static Shapeset::shape_fn_t leg_quad_fn_dyy[] =
    {
      leg_quad_l0_l0yy,   leg_quad_l0_l1yy,   leg_quad_l0_l2yy,   leg_quad_l0_l3yy,   leg_quad_l0_l4yy,
      leg_quad_l0_l5yy,   leg_quad_l0_l6yy,   leg_quad_l0_l7yy,   leg_quad_l0_l8yy,   leg_quad_l0_l9yy,
      leg_quad_l0_l10yy,   leg_quad_l1_l0yy,   leg_quad_l1_l1yy,   leg_quad_l1_l2yy,   leg_quad_l1_l3yy,
      leg_quad_l1_l4yy,   leg_quad_l1_l5yy,   leg_quad_l1_l6yy,   leg_quad_l1_l7yy,   leg_quad_l1_l8yy,
      leg_quad_l1_l9yy,   leg_quad_l1_l10yy,   leg_quad_l2_l0yy,   leg_quad_l2_l1yy,   leg_quad_l2_l2yy,
      leg_quad_l2_l3yy,   leg_quad_l2_l4yy,   leg_quad_l2_l5yy,   leg_quad_l2_l6yy,   leg_quad_l2_l7yy,
      leg_quad_l2_l8yy,   leg_quad_l2_l9yy,   leg_quad_l2_l10yy,   leg_quad_l3_l0yy,   leg_quad_l3_l1yy,
      leg_quad_l3_l2yy,   leg_quad_l3_l3yy,   leg_quad_l3_l4yy,   leg_quad_l3_l5yy,   leg_quad_l3_l6yy,
      leg_quad_l3_l7yy,   leg_quad_l3_l8yy,   leg_quad_l3_l9yy,   leg_quad_l3_l10yy,   leg_quad_l4_l0yy,
      leg_quad_l4_l1yy,   leg_quad_l4_l2yy,   leg_quad_l4_l3yy,   leg_quad_l4_l4yy,   leg_quad_l4_l5yy,
      leg_quad_l4_l6yy,   leg_quad_l4_l7yy,   leg_quad_l4_l8yy,   leg_quad_l4_l9yy,   leg_quad_l4_l10yy,
      leg_quad_l5_l0yy,   leg_quad_l5_l1yy,   leg_quad_l5_l2yy,   leg_quad_l5_l3yy,   leg_quad_l5_l4yy,
      leg_quad_l5_l5yy,   leg_quad_l5_l6yy,   leg_quad_l5_l7yy,   leg_quad_l5_l8yy,   leg_quad_l5_l9yy,
      leg_quad_l5_l10yy,   leg_quad_l6_l0yy,   leg_quad_l6_l1yy,   leg_quad_l6_l2yy,   leg_quad_l6_l3yy,
      leg_quad_l6_l4yy,   leg_quad_l6_l5yy,   leg_quad_l6_l6yy,   leg_quad_l6_l7yy,   leg_quad_l6_l8yy,
      leg_quad_l6_l9yy,   leg_quad_l6_l10yy,   leg_quad_l7_l0yy,   leg_quad_l7_l1yy,   leg_quad_l7_l2yy,
      leg_quad_l7_l3yy,   leg_quad_l7_l4yy,   leg_quad_l7_l5yy,   leg_quad_l7_l6yy,   leg_quad_l7_l7yy,
      leg_quad_l7_l8yy,   leg_quad_l7_l9yy,   leg_quad_l7_l10yy,   leg_quad_l8_l0yy,   leg_quad_l8_l1yy,
      leg_quad_l8_l2yy,   leg_quad_l8_l3yy,   leg_quad_l8_l4yy,   leg_quad_l8_l5yy,   leg_quad_l8_l6yy,
      leg_quad_l8_l7yy,   leg_quad_l8_l8yy,   leg_quad_l8_l9yy,   leg_quad_l8_l10yy,   leg_quad_l9_l0yy,
      leg_quad_l9_l1yy,   leg_quad_l9_l2yy,   leg_quad_l9_l3yy,   leg_quad_l9_l4yy,   leg_quad_l9_l5yy,
      leg_quad_l9_l6yy,   leg_quad_l9_l7yy,   leg_quad_l9_l8yy,   leg_quad_l9_l9yy,   leg_quad_l9_l10yy,
      leg_quad_l10_l0yy,   leg_quad_l10_l1yy,   leg_quad_l10_l2yy,   leg_quad_l10_l3yy,   leg_quad_l10_l4yy,
      leg_quad_l10_l5yy,   leg_quad_l10_l6yy,   leg_quad_l10_l7yy,   leg_quad_l10_l8yy,   leg_quad_l10_l9yy,
      leg_quad_l10_l10yy,
    };
    Shapeset::shape_fn_t* leg_quad_shape_fn_table[1]     = { leg_quad_fn };
    Shapeset::shape_fn_t* leg_quad_shape_fn_table_dx[1]  = { leg_quad_fn_dx };
    Shapeset::shape_fn_t* leg_quad_shape_fn_table_dy[1]  = { leg_quad_fn_dy };
    Shapeset::shape_fn_t* leg_quad_shape_fn_table_dxx[1] = { leg_quad_fn_dxx };
    Shapeset::shape_fn_t* leg_quad_shape_fn_table_dxy[1] = { leg_quad_fn_dxy };
    Shapeset::shape_fn_t* leg_quad_shape_fn_table_dyy[1] = { leg_quad_fn_dyy };

    static int qb_0_0[] = { 0, };
    static int qb_0_1[] = { 0, 1, };
    static int qb_0_2[] = { 0, 1, 2, };
    static int qb_0_3[] = { 0, 1, 2, 3, };
    static int qb_0_4[] = { 0, 1, 2, 3, 4, };
    static int qb_0_5[] = { 0, 1, 2, 3, 4, 5, };
    static int qb_0_6[] = { 0, 1, 2, 3, 4, 5, 6, };
    static int qb_0_7[] = { 0, 1, 2, 3, 4, 5, 6, 7, };
    static int qb_0_8[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, };
    static int qb_0_9[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, };
    static int qb_0_10[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, };
    static int qb_1_0[] = { 0, 11, };
    static int qb_1_1[] = { 0, 1, 11, 12, };
    static int qb_1_2[] = { 0, 1, 2, 11, 12, 13, };
    static int qb_1_3[] = { 0, 1, 2, 3, 11, 12, 13, 14, };
    static int qb_1_4[] = { 0, 1, 2, 3, 4, 11, 12, 13, 14, 15, };
    static int qb_1_5[] = { 0, 1, 2, 3, 4, 5, 11, 12, 13, 14, 15, 16, };
    static int qb_1_6[] = { 0, 1, 2, 3, 4, 5, 6, 11, 12, 13, 14, 15, 16, 17, };
    static int qb_1_7[] = { 0, 1, 2, 3, 4, 5, 6, 7, 11, 12, 13, 14, 15, 16, 17, 18, };
    static int qb_1_8[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 12, 13, 14, 15, 16, 17, 18, 19, };
    static int qb_1_9[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, };
    static int qb_1_10[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, };
    static int qb_2_0[] = { 0, 11, 22, };
    static int qb_2_1[] = { 0, 1, 11, 12, 22, 23, };
    static int qb_2_2[] = { 0, 1, 2, 11, 12, 13, 22, 23, 24, };
    static int qb_2_3[] = { 0, 1, 2, 3, 11, 12, 13, 14, 22, 23, 24, 25, };
    static int qb_2_4[] = { 0, 1, 2, 3, 4, 11, 12, 13, 14, 15, 22, 23, 24, 25, 26, };
    static int qb_2_5[] = { 0, 1, 2, 3, 4, 5, 11, 12, 13, 14, 15, 16, 22, 23, 24, 25, 26, 27, };
    static int qb_2_6[] = { 0, 1, 2, 3, 4, 5, 6, 11, 12, 13, 14, 15, 16, 17, 22, 23, 24, 25, 26, 27, 28, };
    static int qb_2_7[] = { 0, 1, 2, 3, 4, 5, 6, 7, 11, 12, 13, 14, 15, 16, 17, 18, 22, 23, 24, 25, 26, 27, 28, 29, };
    static int qb_2_8[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 12, 13, 14, 15, 16, 17, 18, 19, 22, 23, 24, 25, 26, 27, 28, 29, 30, };
    static int qb_2_9[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, };
    static int qb_2_10[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, };
    static int qb_3_0[] = { 0, 11, 22, 33, };
    static int qb_3_1[] = { 0, 1, 11, 12, 22, 23, 33, 34, };
    static int qb_3_2[] = { 0, 1, 2, 11, 12, 13, 22, 23, 24, 33, 34, 35, };
    static int qb_3_3[] = { 0, 1, 2, 3, 11, 12, 13, 14, 22, 23, 24, 25, 33, 34, 35, 36, };
    static int qb_3_4[] = { 0, 1, 2, 3, 4, 11, 12, 13, 14, 15, 22, 23, 24, 25, 26, 33, 34, 35, 36, 37, };
    static int qb_3_5[] = { 0, 1, 2, 3, 4, 5, 11, 12, 13, 14, 15, 16, 22, 23, 24, 25, 26, 27, 33, 34, 35, 36, 37, 38, };
    static int qb_3_6[] = { 0, 1, 2, 3, 4, 5, 6, 11, 12, 13, 14, 15, 16, 17, 22, 23, 24, 25, 26, 27, 28, 33, 34, 35, 36, 37, 38, 39, };
    static int qb_3_7[] = { 0, 1, 2, 3, 4, 5, 6, 7, 11, 12, 13, 14, 15, 16, 17, 18, 22, 23, 24, 25, 26, 27, 28, 29, 33, 34, 35, 36, 37, 38, 39, 40, };
    static int qb_3_8[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 12, 13, 14, 15, 16, 17, 18, 19, 22, 23, 24, 25, 26, 27, 28, 29, 30, 33, 34, 35, 36, 37, 38, 39, 40, 41, };
    static int qb_3_9[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, };
    static int qb_3_10[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, };
    static int qb_4_0[] = { 0, 11, 22, 33, 44, };
    static int qb_4_1[] = { 0, 1, 11, 12, 22, 23, 33, 34, 44, 45, };
    static int qb_4_2[] = { 0, 1, 2, 11, 12, 13, 22, 23, 24, 33, 34, 35, 44, 45, 46, };
    static int qb_4_3[] = { 0, 1, 2, 3, 11, 12, 13, 14, 22, 23, 24, 25, 33, 34, 35, 36, 44, 45, 46, 47, };
    static int qb_4_4[] = { 0, 1, 2, 3, 4, 11, 12, 13, 14, 15, 22, 23, 24, 25, 26, 33, 34, 35, 36, 37, 44, 45, 46, 47, 48, };
    static int qb_4_5[] = { 0, 1, 2, 3, 4, 5, 11, 12, 13, 14, 15, 16, 22, 23, 24, 25, 26, 27, 33, 34, 35, 36, 37, 38, 44, 45, 46, 47, 48, 49, };
    static int qb_4_6[] = { 0, 1, 2, 3, 4, 5, 6, 11, 12, 13, 14, 15, 16, 17, 22, 23, 24, 25, 26, 27, 28, 33, 34, 35, 36, 37, 38, 39, 44, 45, 46, 47, 48, 49, 50, };
    static int qb_4_7[] = { 0, 1, 2, 3, 4, 5, 6, 7, 11, 12, 13, 14, 15, 16, 17, 18, 22, 23, 24, 25, 26, 27, 28, 29, 33, 34, 35, 36, 37, 38, 39, 40, 44, 45, 46, 47, 48, 49, 50, 51, };
    static int qb_4_8[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 12, 13, 14, 15, 16, 17, 18, 19, 22, 23, 24, 25, 26, 27, 28, 29, 30, 33, 34, 35, 36, 37, 38, 39, 40, 41, 44, 45, 46, 47, 48, 49, 50, 51, 52, };
    static int qb_4_9[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, };
    static int qb_4_10[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, };
    static int qb_5_0[] = { 0, 11, 22, 33, 44, 55, };
    static int qb_5_1[] = { 0, 1, 11, 12, 22, 23, 33, 34, 44, 45, 55, 56, };
    static int qb_5_2[] = { 0, 1, 2, 11, 12, 13, 22, 23, 24, 33, 34, 35, 44, 45, 46, 55, 56, 57, };
    static int qb_5_3[] = { 0, 1, 2, 3, 11, 12, 13, 14, 22, 23, 24, 25, 33, 34, 35, 36, 44, 45, 46, 47, 55, 56, 57, 58, };
    static int qb_5_4[] = { 0, 1, 2, 3, 4, 11, 12, 13, 14, 15, 22, 23, 24, 25, 26, 33, 34, 35, 36, 37, 44, 45, 46, 47, 48, 55, 56, 57, 58, 59, };
    static int qb_5_5[] = { 0, 1, 2, 3, 4, 5, 11, 12, 13, 14, 15, 16, 22, 23, 24, 25, 26, 27, 33, 34, 35, 36, 37, 38, 44, 45, 46, 47, 48, 49, 55, 56, 57, 58, 59, 60, };
    static int qb_5_6[] = { 0, 1, 2, 3, 4, 5, 6, 11, 12, 13, 14, 15, 16, 17, 22, 23, 24, 25, 26, 27, 28, 33, 34, 35, 36, 37, 38, 39, 44, 45, 46, 47, 48, 49, 50, 55, 56, 57, 58, 59, 60, 61, };
    static int qb_5_7[] = { 0, 1, 2, 3, 4, 5, 6, 7, 11, 12, 13, 14, 15, 16, 17, 18, 22, 23, 24, 25, 26, 27, 28, 29, 33, 34, 35, 36, 37, 38, 39, 40, 44, 45, 46, 47, 48, 49, 50, 51, 55, 56, 57, 58, 59, 60, 61, 62, };
    static int qb_5_8[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 12, 13, 14, 15, 16, 17, 18, 19, 22, 23, 24, 25, 26, 27, 28, 29, 30, 33, 34, 35, 36, 37, 38, 39, 40, 41, 44, 45, 46, 47, 48, 49, 50, 51, 52, 55, 56, 57, 58, 59, 60, 61, 62, 63, };
    static int qb_5_9[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, };
    static int qb_5_10[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, };
    static int qb_6_0[] = { 0, 11, 22, 33, 44, 55, 66, };
    static int qb_6_1[] = { 0, 1, 11, 12, 22, 23, 33, 34, 44, 45, 55, 56, 66, 67, };
    static int qb_6_2[] = { 0, 1, 2, 11, 12, 13, 22, 23, 24, 33, 34, 35, 44, 45, 46, 55, 56, 57, 66, 67, 68, };
    static int qb_6_3[] = { 0, 1, 2, 3, 11, 12, 13, 14, 22, 23, 24, 25, 33, 34, 35, 36, 44, 45, 46, 47, 55, 56, 57, 58, 66, 67, 68, 69, };
    static int qb_6_4[] = { 0, 1, 2, 3, 4, 11, 12, 13, 14, 15, 22, 23, 24, 25, 26, 33, 34, 35, 36, 37, 44, 45, 46, 47, 48, 55, 56, 57, 58, 59, 66, 67, 68, 69, 70, };
    static int qb_6_5[] = { 0, 1, 2, 3, 4, 5, 11, 12, 13, 14, 15, 16, 22, 23, 24, 25, 26, 27, 33, 34, 35, 36, 37, 38, 44, 45, 46, 47, 48, 49, 55, 56, 57, 58, 59, 60, 66, 67, 68, 69, 70, 71, };
    static int qb_6_6[] = { 0, 1, 2, 3, 4, 5, 6, 11, 12, 13, 14, 15, 16, 17, 22, 23, 24, 25, 26, 27, 28, 33, 34, 35, 36, 37, 38, 39, 44, 45, 46, 47, 48, 49, 50, 55, 56, 57, 58, 59, 60, 61, 66, 67, 68, 69, 70, 71, 72, };
    static int qb_6_7[] = { 0, 1, 2, 3, 4, 5, 6, 7, 11, 12, 13, 14, 15, 16, 17, 18, 22, 23, 24, 25, 26, 27, 28, 29, 33, 34, 35, 36, 37, 38, 39, 40, 44, 45, 46, 47, 48, 49, 50, 51, 55, 56, 57, 58, 59, 60, 61, 62, 66, 67, 68, 69, 70, 71, 72, 73, };
    static int qb_6_8[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 12, 13, 14, 15, 16, 17, 18, 19, 22, 23, 24, 25, 26, 27, 28, 29, 30, 33, 34, 35, 36, 37, 38, 39, 40, 41, 44, 45, 46, 47, 48, 49, 50, 51, 52, 55, 56, 57, 58, 59, 60, 61, 62, 63, 66, 67, 68, 69, 70, 71, 72, 73, 74, };
    static int qb_6_9[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, };
    static int qb_6_10[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, };
    static int qb_7_0[] = { 0, 11, 22, 33, 44, 55, 66, 77, };
    static int qb_7_1[] = { 0, 1, 11, 12, 22, 23, 33, 34, 44, 45, 55, 56, 66, 67, 77, 78, };
    static int qb_7_2[] = { 0, 1, 2, 11, 12, 13, 22, 23, 24, 33, 34, 35, 44, 45, 46, 55, 56, 57, 66, 67, 68, 77, 78, 79, };
    static int qb_7_3[] = { 0, 1, 2, 3, 11, 12, 13, 14, 22, 23, 24, 25, 33, 34, 35, 36, 44, 45, 46, 47, 55, 56, 57, 58, 66, 67, 68, 69, 77, 78, 79, 80, };
    static int qb_7_4[] = { 0, 1, 2, 3, 4, 11, 12, 13, 14, 15, 22, 23, 24, 25, 26, 33, 34, 35, 36, 37, 44, 45, 46, 47, 48, 55, 56, 57, 58, 59, 66, 67, 68, 69, 70, 77, 78, 79, 80, 81, };
    static int qb_7_5[] = { 0, 1, 2, 3, 4, 5, 11, 12, 13, 14, 15, 16, 22, 23, 24, 25, 26, 27, 33, 34, 35, 36, 37, 38, 44, 45, 46, 47, 48, 49, 55, 56, 57, 58, 59, 60, 66, 67, 68, 69, 70, 71, 77, 78, 79, 80, 81, 82, };
    static int qb_7_6[] = { 0, 1, 2, 3, 4, 5, 6, 11, 12, 13, 14, 15, 16, 17, 22, 23, 24, 25, 26, 27, 28, 33, 34, 35, 36, 37, 38, 39, 44, 45, 46, 47, 48, 49, 50, 55, 56, 57, 58, 59, 60, 61, 66, 67, 68, 69, 70, 71, 72, 77, 78, 79, 80, 81, 82, 83, };
    static int qb_7_7[] = { 0, 1, 2, 3, 4, 5, 6, 7, 11, 12, 13, 14, 15, 16, 17, 18, 22, 23, 24, 25, 26, 27, 28, 29, 33, 34, 35, 36, 37, 38, 39, 40, 44, 45, 46, 47, 48, 49, 50, 51, 55, 56, 57, 58, 59, 60, 61, 62, 66, 67, 68, 69, 70, 71, 72, 73, 77, 78, 79, 80, 81, 82, 83, 84, };
    static int qb_7_8[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 12, 13, 14, 15, 16, 17, 18, 19, 22, 23, 24, 25, 26, 27, 28, 29, 30, 33, 34, 35, 36, 37, 38, 39, 40, 41, 44, 45, 46, 47, 48, 49, 50, 51, 52, 55, 56, 57, 58, 59, 60, 61, 62, 63, 66, 67, 68, 69, 70, 71, 72, 73, 74, 77, 78, 79, 80, 81, 82, 83, 84, 85, };
    static int qb_7_9[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, };
    static int qb_7_10[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, };
    static int qb_8_0[] = { 0, 11, 22, 33, 44, 55, 66, 77, 88, };
    static int qb_8_1[] = { 0, 1, 11, 12, 22, 23, 33, 34, 44, 45, 55, 56, 66, 67, 77, 78, 88, 89, };
    static int qb_8_2[] = { 0, 1, 2, 11, 12, 13, 22, 23, 24, 33, 34, 35, 44, 45, 46, 55, 56, 57, 66, 67, 68, 77, 78, 79, 88, 89, 90, };
    static int qb_8_3[] = { 0, 1, 2, 3, 11, 12, 13, 14, 22, 23, 24, 25, 33, 34, 35, 36, 44, 45, 46, 47, 55, 56, 57, 58, 66, 67, 68, 69, 77, 78, 79, 80, 88, 89, 90, 91, };
    static int qb_8_4[] = { 0, 1, 2, 3, 4, 11, 12, 13, 14, 15, 22, 23, 24, 25, 26, 33, 34, 35, 36, 37, 44, 45, 46, 47, 48, 55, 56, 57, 58, 59, 66, 67, 68, 69, 70, 77, 78, 79, 80, 81, 88, 89, 90, 91, 92, };
    static int qb_8_5[] = { 0, 1, 2, 3, 4, 5, 11, 12, 13, 14, 15, 16, 22, 23, 24, 25, 26, 27, 33, 34, 35, 36, 37, 38, 44, 45, 46, 47, 48, 49, 55, 56, 57, 58, 59, 60, 66, 67, 68, 69, 70, 71, 77, 78, 79, 80, 81, 82, 88, 89, 90, 91, 92, 93, };
    static int qb_8_6[] = { 0, 1, 2, 3, 4, 5, 6, 11, 12, 13, 14, 15, 16, 17, 22, 23, 24, 25, 26, 27, 28, 33, 34, 35, 36, 37, 38, 39, 44, 45, 46, 47, 48, 49, 50, 55, 56, 57, 58, 59, 60, 61, 66, 67, 68, 69, 70, 71, 72, 77, 78, 79, 80, 81, 82, 83, 88, 89, 90, 91, 92, 93, 94, };
    static int qb_8_7[] = { 0, 1, 2, 3, 4, 5, 6, 7, 11, 12, 13, 14, 15, 16, 17, 18, 22, 23, 24, 25, 26, 27, 28, 29, 33, 34, 35, 36, 37, 38, 39, 40, 44, 45, 46, 47, 48, 49, 50, 51, 55, 56, 57, 58, 59, 60, 61, 62, 66, 67, 68, 69, 70, 71, 72, 73, 77, 78, 79, 80, 81, 82, 83, 84, 88, 89, 90, 91, 92, 93, 94, 95, };
    static int qb_8_8[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 12, 13, 14, 15, 16, 17, 18, 19, 22, 23, 24, 25, 26, 27, 28, 29, 30, 33, 34, 35, 36, 37, 38, 39, 40, 41, 44, 45, 46, 47, 48, 49, 50, 51, 52, 55, 56, 57, 58, 59, 60, 61, 62, 63, 66, 67, 68, 69, 70, 71, 72, 73, 74, 77, 78, 79, 80, 81, 82, 83, 84, 85, 88, 89, 90, 91, 92, 93, 94, 95, 96, };
    static int qb_8_9[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, };
    static int qb_8_10[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, };
    static int qb_9_0[] = { 0, 11, 22, 33, 44, 55, 66, 77, 88, 99, };
    static int qb_9_1[] = { 0, 1, 11, 12, 22, 23, 33, 34, 44, 45, 55, 56, 66, 67, 77, 78, 88, 89, 99, 100, };
    static int qb_9_2[] = { 0, 1, 2, 11, 12, 13, 22, 23, 24, 33, 34, 35, 44, 45, 46, 55, 56, 57, 66, 67, 68, 77, 78, 79, 88, 89, 90, 99, 100, 101, };
    static int qb_9_3[] = { 0, 1, 2, 3, 11, 12, 13, 14, 22, 23, 24, 25, 33, 34, 35, 36, 44, 45, 46, 47, 55, 56, 57, 58, 66, 67, 68, 69, 77, 78, 79, 80, 88, 89, 90, 91, 99, 100, 101, 102, };
    static int qb_9_4[] = { 0, 1, 2, 3, 4, 11, 12, 13, 14, 15, 22, 23, 24, 25, 26, 33, 34, 35, 36, 37, 44, 45, 46, 47, 48, 55, 56, 57, 58, 59, 66, 67, 68, 69, 70, 77, 78, 79, 80, 81, 88, 89, 90, 91, 92, 99, 100, 101, 102, 103, };
    static int qb_9_5[] = { 0, 1, 2, 3, 4, 5, 11, 12, 13, 14, 15, 16, 22, 23, 24, 25, 26, 27, 33, 34, 35, 36, 37, 38, 44, 45, 46, 47, 48, 49, 55, 56, 57, 58, 59, 60, 66, 67, 68, 69, 70, 71, 77, 78, 79, 80, 81, 82, 88, 89, 90, 91, 92, 93, 99, 100, 101, 102, 103, 104, };
    static int qb_9_6[] = { 0, 1, 2, 3, 4, 5, 6, 11, 12, 13, 14, 15, 16, 17, 22, 23, 24, 25, 26, 27, 28, 33, 34, 35, 36, 37, 38, 39, 44, 45, 46, 47, 48, 49, 50, 55, 56, 57, 58, 59, 60, 61, 66, 67, 68, 69, 70, 71, 72, 77, 78, 79, 80, 81, 82, 83, 88, 89, 90, 91, 92, 93, 94, 99, 100, 101, 102, 103, 104, 105, };
    static int qb_9_7[] = { 0, 1, 2, 3, 4, 5, 6, 7, 11, 12, 13, 14, 15, 16, 17, 18, 22, 23, 24, 25, 26, 27, 28, 29, 33, 34, 35, 36, 37, 38, 39, 40, 44, 45, 46, 47, 48, 49, 50, 51, 55, 56, 57, 58, 59, 60, 61, 62, 66, 67, 68, 69, 70, 71, 72, 73, 77, 78, 79, 80, 81, 82, 83, 84, 88, 89, 90, 91, 92, 93, 94, 95, 99, 100, 101, 102, 103, 104, 105, 106, };
    static int qb_9_8[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 12, 13, 14, 15, 16, 17, 18, 19, 22, 23, 24, 25, 26, 27, 28, 29, 30, 33, 34, 35, 36, 37, 38, 39, 40, 41, 44, 45, 46, 47, 48, 49, 50, 51, 52, 55, 56, 57, 58, 59, 60, 61, 62, 63, 66, 67, 68, 69, 70, 71, 72, 73, 74, 77, 78, 79, 80, 81, 82, 83, 84, 85, 88, 89, 90, 91, 92, 93, 94, 95, 96, 99, 100, 101, 102, 103, 104, 105, 106, 107, };
    static int qb_9_9[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, };
    static int qb_9_10[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, };
    static int qb_10_0[] = { 0, 11, 22, 33, 44, 55, 66, 77, 88, 99, 110, };
    static int qb_10_1[] = { 0, 1, 11, 12, 22, 23, 33, 34, 44, 45, 55, 56, 66, 67, 77, 78, 88, 89, 99, 100, 110, 111, };
    static int qb_10_2[] = { 0, 1, 2, 11, 12, 13, 22, 23, 24, 33, 34, 35, 44, 45, 46, 55, 56, 57, 66, 67, 68, 77, 78, 79, 88, 89, 90, 99, 100, 101, 110, 111, 112, };
    static int qb_10_3[] = { 0, 1, 2, 3, 11, 12, 13, 14, 22, 23, 24, 25, 33, 34, 35, 36, 44, 45, 46, 47, 55, 56, 57, 58, 66, 67, 68, 69, 77, 78, 79, 80, 88, 89, 90, 91, 99, 100, 101, 102, 110, 111, 112, 113, };
    static int qb_10_4[] = { 0, 1, 2, 3, 4, 11, 12, 13, 14, 15, 22, 23, 24, 25, 26, 33, 34, 35, 36, 37, 44, 45, 46, 47, 48, 55, 56, 57, 58, 59, 66, 67, 68, 69, 70, 77, 78, 79, 80, 81, 88, 89, 90, 91, 92, 99, 100, 101, 102, 103, 110, 111, 112, 113, 114, };
    static int qb_10_5[] = { 0, 1, 2, 3, 4, 5, 11, 12, 13, 14, 15, 16, 22, 23, 24, 25, 26, 27, 33, 34, 35, 36, 37, 38, 44, 45, 46, 47, 48, 49, 55, 56, 57, 58, 59, 60, 66, 67, 68, 69, 70, 71, 77, 78, 79, 80, 81, 82, 88, 89, 90, 91, 92, 93, 99, 100, 101, 102, 103, 104, 110, 111, 112, 113, 114, 115, };
    static int qb_10_6[] = { 0, 1, 2, 3, 4, 5, 6, 11, 12, 13, 14, 15, 16, 17, 22, 23, 24, 25, 26, 27, 28, 33, 34, 35, 36, 37, 38, 39, 44, 45, 46, 47, 48, 49, 50, 55, 56, 57, 58, 59, 60, 61, 66, 67, 68, 69, 70, 71, 72, 77, 78, 79, 80, 81, 82, 83, 88, 89, 90, 91, 92, 93, 94, 99, 100, 101, 102, 103, 104, 105, 110, 111, 112, 113, 114, 115, 116, };
    static int qb_10_7[] = { 0, 1, 2, 3, 4, 5, 6, 7, 11, 12, 13, 14, 15, 16, 17, 18, 22, 23, 24, 25, 26, 27, 28, 29, 33, 34, 35, 36, 37, 38, 39, 40, 44, 45, 46, 47, 48, 49, 50, 51, 55, 56, 57, 58, 59, 60, 61, 62, 66, 67, 68, 69, 70, 71, 72, 73, 77, 78, 79, 80, 81, 82, 83, 84, 88, 89, 90, 91, 92, 93, 94, 95, 99, 100, 101, 102, 103, 104, 105, 106, 110, 111, 112, 113, 114, 115, 116, 117, };
    static int qb_10_8[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 12, 13, 14, 15, 16, 17, 18, 19, 22, 23, 24, 25, 26, 27, 28, 29, 30, 33, 34, 35, 36, 37, 38, 39, 40, 41, 44, 45, 46, 47, 48, 49, 50, 51, 52, 55, 56, 57, 58, 59, 60, 61, 62, 63, 66, 67, 68, 69, 70, 71, 72, 73, 74, 77, 78, 79, 80, 81, 82, 83, 84, 85, 88, 89, 90, 91, 92, 93, 94, 95, 96, 99, 100, 101, 102, 103, 104, 105, 106, 107, 110, 111, 112, 113, 114, 115, 116, 117, 118, };
    static int qb_10_9[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, };
    static int qb_10_10[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, };

    #define NULL16 NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL

    int* leg_quad_bubble_indices[] =
    {
      qb_0_0,  qb_0_1,  qb_0_2,  qb_0_3,  qb_0_4,  qb_0_5,  qb_0_6,  qb_0_7,  qb_0_8,  qb_0_9,  qb_0_10,   NULL, NULL, NULL, NULL, NULL, NULL16,
      qb_1_0,  qb_1_1,  qb_1_2,  qb_1_3,  qb_1_4,  qb_1_5,  qb_1_6,  qb_1_7,  qb_1_8,  qb_1_9,  qb_1_10,   NULL, NULL, NULL, NULL, NULL, NULL16,
      qb_2_0,  qb_2_1,  qb_2_2,  qb_2_3,  qb_2_4,  qb_2_5,  qb_2_6,  qb_2_7,  qb_2_8,  qb_2_9,  qb_2_10,   NULL, NULL, NULL, NULL, NULL, NULL16,
      qb_3_0,  qb_3_1,  qb_3_2,  qb_3_3,  qb_3_4,  qb_3_5,  qb_3_6,  qb_3_7,  qb_3_8,  qb_3_9,  qb_3_10,   NULL, NULL, NULL, NULL, NULL, NULL16,
      qb_4_0,  qb_4_1,  qb_4_2,  qb_4_3,  qb_4_4,  qb_4_5,  qb_4_6,  qb_4_7,  qb_4_8,  qb_4_9,  qb_4_10,   NULL, NULL, NULL, NULL, NULL, NULL16,
      qb_5_0,  qb_5_1,  qb_5_2,  qb_5_3,  qb_5_4,  qb_5_5,  qb_5_6,  qb_5_7,  qb_5_8,  qb_5_9,  qb_5_10,   NULL, NULL, NULL, NULL, NULL, NULL16,
      qb_6_0,  qb_6_1,  qb_6_2,  qb_6_3,  qb_6_4,  qb_6_5,  qb_6_6,  qb_6_7,  qb_6_8,  qb_6_9,  qb_6_10,   NULL, NULL, NULL, NULL, NULL, NULL16,
      qb_7_0,  qb_7_1,  qb_7_2,  qb_7_3,  qb_7_4,  qb_7_5,  qb_7_6,  qb_7_7,  qb_7_8,  qb_7_9,  qb_7_10,   NULL, NULL, NULL, NULL, NULL, NULL16,
      qb_8_0,  qb_8_1,  qb_8_2,  qb_8_3,  qb_8_4,  qb_8_5,  qb_8_6,  qb_8_7,  qb_8_8,  qb_8_9,  qb_8_10,   NULL, NULL, NULL, NULL, NULL, NULL16,
      qb_9_0,  qb_9_1,  qb_9_2,  qb_9_3,  qb_9_4,  qb_9_5,  qb_9_6,  qb_9_7,  qb_9_8,  qb_9_9,  qb_9_10,   NULL, NULL, NULL, NULL, NULL, NULL16,
      qb_10_0,  qb_10_1,  qb_10_2,  qb_10_3,  qb_10_4,  qb_10_5,  qb_10_6,  qb_10_7,  qb_10_8,  qb_10_9,  qb_10_10,   NULL, NULL, NULL, NULL, NULL, NULL16,
    };

    #define zero16  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

    int leg_quad_bubble_count[] =
    {
      1,  2,  3,  4,  5,  6,  7,  8,  9,  10,  11,  0,  0,  0,  0,  0, zero16,
      2,  4,  6,  8,  10,  12,  14,  16,  18,  20,  22,  0,  0,  0,  0,  0, zero16,
      3,  6,  9,  12,  15,  18,  21,  24,  27,  30,  33,  0,  0,  0,  0,  0, zero16,
      4,  8,  12,  16,  20,  24,  28,  32,  36,  40,  44,  0,  0,  0,  0,  0, zero16,
      5,  10,  15,  20,  25,  30,  35,  40,  45,  50,  55,  0,  0,  0,  0,  0, zero16,
      6,  12,  18,  24,  30,  36,  42,  48,  54,  60,  66,  0,  0,  0,  0,  0, zero16,
      7,  14,  21,  28,  35,  42,  49,  56,  63,  70,  77,  0,  0,  0,  0,  0, zero16,
      8,  16,  24,  32,  40,  48,  56,  64,  72,  80,  88,  0,  0,  0,  0,  0, zero16,
      9,  18,  27,  36,  45,  54,  63,  72,  81,  90,  99,  0,  0,  0,  0,  0, zero16,
      10,  20,  30,  40,  50,  60,  70,  80,  90,  100,  110,  0,  0,  0,  0,  0, zero16,
      11,  22,  33,  44,  55,  66,  77,  88,  99,  110,  120,  0,  0,  0,  0,  0, zero16,
    };

    int leg_quad_vertex_indices[4] = { -1, -1, -1, -1 };

    static int leg_quad_edge_indices_0[22] =  { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1  };
    static int leg_quad_edge_indices_1[22] =  { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1  };
    static int leg_quad_edge_indices_2[22] =  { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1  };
    static int leg_quad_edge_indices_3[22] =  { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1  };

    int* leg_quad_edge_indices[4] =
    {
      leg_quad_edge_indices_0,
      leg_quad_edge_indices_1,
      leg_quad_edge_indices_2,
      leg_quad_edge_indices_3
    };

    #define oo H2D_MAKE_QUAD_ORDER
    #define XX(a, b) oo(a, b), oo(a, b)

    int leg_quad_index_to_order[] =
    {
      oo(0, 0),   oo(0, 1),   oo(0, 2),   oo(0, 3),   oo(0, 4),   oo(0, 5),   oo(0, 6),   oo(0, 7),   oo(0, 8),   oo(0, 9),   oo(0, 10),
      oo(1, 0),   oo(1, 1),   oo(1, 2),   oo(1, 3),   oo(1, 4),   oo(1, 5),   oo(1, 6),   oo(1, 7),   oo(1, 8),   oo(1, 9),   oo(1, 10),
      oo(2, 0),   oo(2, 1),   oo(2, 2),   oo(2, 3),   oo(2, 4),   oo(2, 5),   oo(2, 6),   oo(2, 7),   oo(2, 8),   oo(2, 9),   oo(2, 10),
      oo(3, 0),   oo(3, 1),   oo(3, 2),   oo(3, 3),   oo(3, 4),   oo(3, 5),   oo(3, 6),   oo(3, 7),   oo(3, 8),   oo(3, 9),   oo(3, 10),
      oo(4, 0),   oo(4, 1),   oo(4, 2),   oo(4, 3),   oo(4, 4),   oo(4, 5),   oo(4, 6),   oo(4, 7),   oo(4, 8),   oo(4, 9),   oo(4, 10),
      oo(5, 0),   oo(5, 1),   oo(5, 2),   oo(5, 3),   oo(5, 4),   oo(5, 5),   oo(5, 6),   oo(5, 7),   oo(5, 8),   oo(5, 9),   oo(5, 10),
      oo(6, 0),   oo(6, 1),   oo(6, 2),   oo(6, 3),   oo(6, 4),   oo(6, 5),   oo(6, 6),   oo(6, 7),   oo(6, 8),   oo(6, 9),   oo(6, 10),
      oo(7, 0),   oo(7, 1),   oo(7, 2),   oo(7, 3),   oo(7, 4),   oo(7, 5),   oo(7, 6),   oo(7, 7),   oo(7, 8),   oo(7, 9),   oo(7, 10),
      oo(8, 0),   oo(8, 1),   oo(8, 2),   oo(8, 3),   oo(8, 4),   oo(8, 5),   oo(8, 6),   oo(8, 7),   oo(8, 8),   oo(8, 9),   oo(8, 10),
      oo(9, 0),   oo(9, 1),   oo(9, 2),   oo(9, 3),   oo(9, 4),   oo(9, 5),   oo(9, 6),   oo(9, 7),   oo(9, 8),   oo(9, 9),   oo(9, 10),
      oo(10, 0),   oo(10, 1),   oo(10, 2),   oo(10, 3),   oo(10, 4),   oo(10, 5),   oo(10, 6),   oo(10, 7),   oo(10, 8),   oo(10, 9),   oo(10, 10),
    };

    //// triangle legendre shapeset /////////////////////////////////////////////////////////////////

    static double leg_tri_l0_l0(double x, double y)
    {
      return 1.0;
    }

    static double leg_tri_l0_l0x(double x, double y)
    {
      return 0.0;
    }

    static double leg_tri_l0_l0y(double x, double y)
    {
      return 0.0;
    }

    static double leg_tri_l0_l0xx(double x, double y)
    {
      return 0.0;
    }

    static double leg_tri_l0_l0xy(double x, double y)
    {
      return 0.0;
    }

    static double leg_tri_l0_l0yy(double x, double y)
    {
      return 0.0;
    }

    // number 1
    inline double leg_f1(double x, double y)
    {
      return lambda2(x, y);
    }

    inline double leg_f1_dx(double x, double y)
    {
      return lambda2x(x, y);
    }

    inline double leg_f1_dy(double x, double y)
    {
      return lambda2y(x, y);
    }

    inline double leg_f1_dxx(double x, double y)
    {
      return 0.0;
    }

    inline double leg_f1_dyy(double x, double y)
    {
      return 0.0;
    }

    inline double leg_f1_dxy(double x, double y)
    {
      return 0.0;
    }

    // number 2
    inline double leg_f2(double x, double y)
    {
      return lambda3(x, y);
    }

    inline double leg_f2_dx(double x, double y)
    {
      return lambda3x(x, y);
    }

    inline double leg_f2_dy(double x, double y)
    {
      return lambda3y(x, y);
    }

    inline double leg_f2_dxx(double x, double y)
    {
      return 0.0;
    }

    inline double leg_f2_dyy(double x, double y)
    {
      return 0.0;
    }

    inline double leg_f2_dxy(double x, double y)
    {
      return 0.0;
    }

    // number 3
    inline double leg_f3(double x, double y)
    {
      return lambda1(x, y);
    }

    inline double leg_f3_dx(double x, double y)
    {
      return lambda1x(x, y);
    }

    inline double leg_f3_dy(double x, double y)
    {
      return lambda1y(x, y);
    }

    inline double leg_f3_dxx(double x, double y)
    {
      return 0.0;
    }

    inline double leg_f3_dyy(double x, double y)
    {
      return 0.0;
    }

    inline double leg_f3_dxy(double x, double y)
    {
      return 0.0;
    }

    // ORDER 2

    // Edge functions, order 2

    // number 4
    inline double leg_f4(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      return .15309310892394863113733025466911821*a2-.15309310892394863113733025466911821*b2;
    }

    inline double leg_f4_dx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return .61237243569579452454932101867647285*a;
    }

    inline double leg_f4_dy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return .30618621784789726227466050933823642*a + .30618621784789726227466050933823642*b;
    }

    inline double leg_f4_dxx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return 1.2247448713915890490986420373529457;
    }

    inline double leg_f4_dyy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return 0.;
    }

    inline double leg_f4_dxy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return .61237243569579452454932101867647285;
    }

    // number 5
    inline double leg_f5(double x, double y)
    {
      return -.61237243569579452454932101867647285*(1. + y)*(x + 1.);
    }

    inline double leg_f5_dx(double x, double y)
    {
      return -.61237243569579452454932101867647285-.61237243569579452454932101867647285*y;
    }

    inline double leg_f5_dy(double x, double y)
    {
      return -.61237243569579452454932101867647285*x-.61237243569579452454932101867647285;
    }

    inline double leg_f5_dxx(double x, double y)
    {
      return 0.;
    }

    inline double leg_f5_dyy(double x, double y)
    {
      return 0.;
    }

    inline double leg_f5_dxy(double x, double y)
    {
      return -.61237243569579452454932101867647285;
    }

    // number 6
    inline double leg_f6(double x, double y)
    {
      return .61237243569579452454932101867647285*(1. + y)*(x + y);
    }

    inline double leg_f6_dx(double x, double y)
    {
      return .61237243569579452454932101867647285 + .61237243569579452454932101867647285*y;
    }

    inline double leg_f6_dy(double x, double y)
    {
      return .61237243569579452454932101867647285*x + 1.2247448713915890490986420373529457*y + .61237243569579452454932101867647285;
    }

    inline double leg_f6_dxx(double x, double y)
    {
      return 0.;
    }

    inline double leg_f6_dyy(double x, double y)
    {
      return 1.2247448713915890490986420373529457;
    }

    inline double leg_f6_dxy(double x, double y)
    {
      return .61237243569579452454932101867647285;
    }

    // ORDER 3

    // Edge functions, order 3

    // number 7
    inline double leg_f7_0(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      return (-.98821176880261854124965423263522453e-1*b2 + .98821176880261854124965423263522453e-1*a2)*a;
    }

    inline double leg_f7_1(double x, double y)
    {
      return -leg_f7_0(x, y);
    }

    inline double leg_f7_dx_0(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      return .59292706128157112474979253958113472*a2-.19764235376052370824993084652704491*b2;
    }

    inline double leg_f7_dx_1(double x, double y)
    {
      return -leg_f7_dx_0(x, y);
    }

    inline double leg_f7_dy_0(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      return -.98821176880261854124965423263522453e-1*b2 + (.19764235376052370824993084652704491*b + .29646353064078556237489626979056736*a)*a;
    }

    inline double leg_f7_dy_1(double x, double y)
    {
      return -leg_f7_dy_0(x, y);
    }

    inline double leg_f7_dxx_0(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return 2.3717082451262844989991701583245389*a;
    }

    inline double leg_f7_dxx_1(double x, double y)
    {
      return -leg_f7_dxx_0(x, y);
    }

    inline double leg_f7_dyy_0(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return .39528470752104741649986169305408981*a + .39528470752104741649986169305408981*b;
    }

    inline double leg_f7_dyy_1(double x, double y)
    {
      return -leg_f7_dyy_0(x, y);
    }

    inline double leg_f7_dxy_0(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return 1.1858541225631422494995850791622694*a + .39528470752104741649986169305408981*b;
    }

    inline double leg_f7_dxy_1(double x, double y)
    {
      return -leg_f7_dxy_0(x, y);
    }

    // number 8
    inline double leg_f8_0(double x, double y)
    {
      return -.79056941504209483299972338610817962*(1. + y)*y*(x + 1.);
    }

    // number 8
    inline double leg_f8_1(double x, double y)
    {
      return -leg_f8_0(x, y);
    }

    inline double leg_f8_dx_0(double x, double y)
    {
      return -.79056941504209483299972338610817962*(1. + y)*y;
    }

    inline double leg_f8_dx_1(double x, double y)
    {
      return -leg_f8_dx_0(x, y);
    }

    inline double leg_f8_dy_0(double x, double y)
    {
      return -.79056941504209483299972338610817962*(x + 1.)*(2.*y + 1.);
    }

    inline double leg_f8_dy_1(double x, double y)
    {
      return -leg_f8_dy_0(x, y);
    }

    inline double leg_f8_dxx_0(double x, double y)
    {
      return 0.;
    }

    inline double leg_f8_dxx_1(double x, double y)
    {
      return -leg_f8_dxx_0(x, y);
    }

    inline double leg_f8_dyy_0(double x, double y)
    {
      return -1.5811388300841896659994467722163592*x-1.5811388300841896659994467722163592;
    }

    inline double leg_f8_dyy_1(double x, double y)
    {
      return -leg_f8_dyy_0(x, y);
    }

    inline double leg_f8_dxy_0(double x, double y)
    {
      return -1.5811388300841896659994467722163592*y-.79056941504209483299972338610817962;
    }

    inline double leg_f8_dxy_1(double x, double y)
    {
      return -leg_f8_dxy_0(x, y);
    }

    // number 9
    inline double leg_f9_0(double x, double y)
    {
      return -.79056941504209483299972338610817962*(1. + y)*y*(x + y);
    }

    // number 9
    inline double leg_f9_1(double x, double y)
    {
      return -leg_f9_0(x, y);
    }

    inline double leg_f9_dx_0(double x, double y)
    {
      return -.79056941504209483299972338610817962*(1. + y)*y;
    }

    inline double leg_f9_dx_1(double x, double y)
    {
      return -leg_f9_dx_0(x, y);
    }

    inline double leg_f9_dy_0(double x, double y)
    {
      return -.79056941504209483299972338610817962*y*(3.*y + 2.)-.79056941504209483299972338610817962*(2.*y + 1.)*x;
    }

    inline double leg_f9_dy_1(double x, double y)
    {
      return -leg_f9_dy_0(x, y);
    }

    inline double leg_f9_dxx_0(double x, double y)
    {
      return 0.;
    }

    inline double leg_f9_dxx_1(double x, double y)
    {
      return -leg_f9_dxx_0(x, y);
    }

    inline double leg_f9_dyy_0(double x, double y)
    {
      return -1.5811388300841896659994467722163592*x-4.7434164902525689979983403166490778*y-1.5811388300841896659994467722163592;
    }

    inline double leg_f9_dyy_1(double x, double y)
    {
      return -leg_f9_dyy_0(x, y);
    }

    inline double leg_f9_dxy_0(double x, double y)
    {
      return -1.5811388300841896659994467722163592*y-.79056941504209483299972338610817962;
    }

    inline double leg_f9_dxy_1(double x, double y)
    {
      return -leg_f9_dxy_0(x, y);
    }

    // Bubble functions, order 3

    // number 10
    inline double leg_f10(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      return .48412291827592711064740817497279995*(.61237243569579452454932101867647285*a2-.61237243569579452454932101867647285*b2)*(1. + y);
    }

    inline double leg_f10_dx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return 1.1858541225631422494995850791622694*a*(1. + y);
    }

    inline double leg_f10_dy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return (.59292706128157112474979253958113472 + .59292706128157112474979253958113472*y-.29646353064078556237489626979056736*b)*b+
          (.59292706128157112474979253958113472 + .59292706128157112474979253958113472*y + .29646353064078556237489626979056736*a)*a;
    }

    inline double leg_f10_dxx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return 2.3717082451262844989991701583245389 + 2.3717082451262844989991701583245389*y;
    }

    inline double leg_f10_dyy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return 1.1858541225631422494995850791622694*a + 1.1858541225631422494995850791622694*b;
    }

    inline double leg_f10_dxy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return 1.1858541225631422494995850791622694 + 1.1858541225631422494995850791622694*y + 1.1858541225631422494995850791622694*a;
    }

    // ORDER 4

    // Edge functions, order 4

    // number 11
    inline double leg_f11(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b4 = b2*b2;
      return .14615849167085708537436518485611521e-1*b4 + (-.87695095002514251224619110913669124e-1*b2 + .73079245835428542687182592428057604e-1*a2)*a2;
    }

    inline double leg_f11_dx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      return (-.35078038001005700489847644365467650*b2 + .58463396668342834149746073942446083*a2)*a;
    }

    inline double leg_f11_dy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      return -.58463396668342834149746073942446083e-1*b3 + (-.17539019000502850244923822182733825*b2+
          (.17539019000502850244923822182733825*b + .29231698334171417074873036971223041*a)*a)*a;
    }

    inline double leg_f11_dxx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      return 3.5078038001005700489847644365467650*a2-.70156076002011400979695288730935299*b2;
    }

    inline double leg_f11_dyy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      return (.70156076002011400979695288730935299*b + .70156076002011400979695288730935299*a)*a;
    }

    inline double leg_f11_dxy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      return -.35078038001005700489847644365467650*b2 + (.70156076002011400979695288730935299*b + 1.7539019000502850244923822182733825*a)*a;
    }

    // number 12
    inline double leg_f12(double x, double y)
    {
      double y2 = y*y;
      return -.23385358667337133659898429576978433*(x + 1.)*(5.*y2-1.)*(1. + y);
    }

    inline double leg_f12_dx(double x, double y)
    {
      double y2 = y*y;
      return -.23385358667337133659898429576978433*(5.*y2-1.)*(1. + y);
    }

    inline double leg_f12_dy(double x, double y)
    {
      double y2 = y*y;
      return -.23385358667337133659898429576978433*(x + 1.)*(-1. + (10. + 15.*y)*y);
    }

    inline double leg_f12_dxx(double x, double y)
    {
      return 0.;
    }

    inline double leg_f12_dyy(double x, double y)
    {
      return -2.3385358667337133659898429576978433*(x + 1.)*(1. + 3.*y);
    }

    inline double leg_f12_dxy(double x, double y)
    {
      return .23385358667337133659898429576978433-.23385358667337133659898429576978433*(10. + 15.*y)*y;
    }

    // number 13
    inline double leg_f13(double x, double y)
    {
      double y2 = y*y;
      return .23385358667337133659898429576978433*(x + y)*(5.*y2-1.)*(1. + y);
    }

    inline double leg_f13_dx(double x, double y)
    {
      double y2 = y*y;
      return .23385358667337133659898429576978433*(5.*y2-1.)*(1. + y);
    }

    inline double leg_f13_dy(double x, double y)
    {
      double y2 = y*y;
      return -.23385358667337133659898429576978433 + .23385358667337133659898429576978433*(-2. + (15. + 20.*y)*y)*y+
          .23385358667337133659898429576978433*(-1. + (10. + 15.*y)*y)*x;
    }

    inline double leg_f13_dxx(double x, double y)
    {
      return 0.;
    }

    inline double leg_f13_dyy(double x, double y)
    {
      return -.46770717334674267319796859153956866 + .46770717334674267319796859153956866*(15. + 30.*y)*y+
          .46770717334674267319796859153956866*(5. + 15.*y)*x;
    }

    inline double leg_f13_dxy(double x, double y)
    {
      return -.23385358667337133659898429576978433 + .23385358667337133659898429576978433*(10. + 15.*y)*y;
    }

    // Bubble functions, order 4

    // number 14
    inline double leg_f14(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      return .84229416444580658121520354710938710*(.61237243569579452454932101867647285*a2-.61237243569579452454932101867647285*b2)*
          (.25000000000000000000000000000000000 + (1.5000000000000000000000000000000000 + 1.2500000000000000000000000000000000*y)*y);
    }

    inline double leg_f14_dx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return 2.0631909162161306774360682303013260*a*(.25000000000000000000000000000000000+
          (1.5000000000000000000000000000000000 + 1.2500000000000000000000000000000000*y)*y);
    }

    inline double leg_f14_dy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return (.25789886452701633467950852878766574 + 1.0315954581080653387180341151506630*(1.5000000000000000000000000000000000+
          1.2500000000000000000000000000000000*y)*y-.51579772905403266935901705757533149*(2.5000000000000000000000000000000000*y + 1.5000000000000000000000000000000000)*b)*b + (.25789886452701633467950852878766574 + 1.0315954581080653387180341151506630*(1.5000000000000000000000000000000000 + 1.2500000000000000000000000000000000*y)*y + .51579772905403266935901705757533149*(2.5000000000000000000000000000000000*y + 1.5000000000000000000000000000000000)*a)*a;
    }

    inline double leg_f14_dxx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return 1.0315954581080653387180341151506630 + 4.1263818324322613548721364606026519*(1.5000000000000000000000000000000000+
          1.2500000000000000000000000000000000*y)*y;
    }

    inline double leg_f14_dyy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return (5.1579772905403266935901705757533151*y + 3.0947863743241960161541023454519890-1.2894943226350816733975426439383287*b)*b+
          (5.1579772905403266935901705757533151*y + 3.0947863743241960161541023454519890 + 1.2894943226350816733975426439383287*a)*a;
    }

    inline double leg_f14_dxy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return .51579772905403266935901705757533151 + 2.0631909162161306774360682303013260*(1.5000000000000000000000000000000000+
          1.2500000000000000000000000000000000*y)*y + 2.0631909162161306774360682303013260*(2.5000000000000000000000000000000000*y+
          1.5000000000000000000000000000000000)*a;
    }

    // number 15
    inline double leg_f15(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      return .44370598373247120319254962186712128*(-.79056941504209483299972338610817962*b2 + .79056941504209483299972338610817962*a2)*a*(1. + y);
    }

    inline double leg_f15_dx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      return .44370598373247120319254962186712128*(4.7434164902525689979983403166490778*a2-1.5811388300841896659994467722163592*b2)*(1. + y);
    }

    inline double leg_f15_dy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      return -.35078038001005700489847644365467649*b2*(1. + y) + ((.70156076002011400979695288730935298 + .70156076002011400979695288730935298*y-
          .35078038001005700489847644365467649*b)*b + (1.0523411400301710146954293309640295 + 1.0523411400301710146954293309640295*y+
          .35078038001005700489847644365467649*a)*a)*a;
    }

    inline double leg_f15_dxx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return 8.4187291202413681175634346477122357*a*(1. + y);
    }

    inline double leg_f15_dyy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return (1.4031215200402280195939057746187060 + 1.4031215200402280195939057746187060*y-.70156076002011400979695288730935298*b)*b+
          (1.4031215200402280195939057746187060 + 1.4031215200402280195939057746187060*y + 1.4031215200402280195939057746187060*b+
          2.1046822800603420293908586619280589*a)*a;
    }

    inline double leg_f15_dxy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return (1.4031215200402280195939057746187060 + 1.4031215200402280195939057746187060*y-.70156076002011400979695288730935298*b)*b+
          (4.2093645601206840587817173238561178 + 4.2093645601206840587817173238561178*y + 2.1046822800603420293908586619280589*a)*a;
    }

    // ORDER 5

    // Edge functions, order 5

    // number 16
    inline double leg_f16_0(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b4 = b2*b2;
      return (.24859222776089561404717184605248599e-1*b4 + (-.82864075920298538015723948684161998e-1*b2 + .58004853144208976611006764078913399e-1*a2)*a2)*a;
    }

    inline double leg_f16_1(double x, double y)
    {
      return -leg_f16_0(x, y);
    }

    inline double leg_f16_dx_0(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      return .49718445552179122809434369210497199e-1*b4 + (-.49718445552179122809434369210497199*b2 + .58004853144208976611006764078913399*a2)*a2;
    }

    inline double leg_f16_dx_1(double x, double y)
    {
      return -leg_f16_dx_0(x, y);
    }

    inline double leg_f16_dy_0(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      return .24859222776089561404717184605248599e-1*b4 + (-.99436891104358245618868738420994398e-1*b3 + (-.24859222776089561404717184605248599*b2+
          (.16572815184059707603144789736832400*b + .29002426572104488305503382039456699*a)*a)*a)*a;
    }

    inline double leg_f16_dy_1(double x, double y)
    {
      return -leg_f16_dy_0(x, y);
    }

    inline double leg_f16_dxx_0(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      return (-1.9887378220871649123773747684198880*b2 + 4.6403882515367181288805411263130719*a2)*a;
    }

    inline double leg_f16_dxx_1(double x, double y)
    {
      return -leg_f16_dxx_0(x, y);
    }

    inline double leg_f16_dyy_0(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      return -.19887378220871649123773747684198880*b3 + (-.19887378220871649123773747684198880*b2+
          (.99436891104358245618868738420994398*b + .99436891104358245618868738420994398*a)*a)*a;
    }

    inline double leg_f16_dyy_1(double x, double y)
    {
      return -leg_f16_dyy_0(x, y);
    }

    inline double leg_f16_dxy_0(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      return -.19887378220871649123773747684198880*b3 + (-.99436891104358245618868738420994398*b2+
          (.99436891104358245618868738420994398*b + 2.3201941257683590644402705631565359*a)*a)*a;
    }

    inline double leg_f16_dxy_1(double x, double y)
    {
      return -leg_f16_dxy_0(x, y);
    }

    // number 17
    inline double leg_f17_0(double x, double y)
    {
      double y2 = y*y;
      return -.26516504294495532165031663578931839*(x + 1.)*(7.*y2-3.)*y*(1. + y);
    }

    // number 17
    inline double leg_f17_1(double x, double y)
    {
      return -leg_f17_0(x, y);
    }

    inline double leg_f17_dx_0(double x, double y)
    {
      double y2 = y*y;
      return -.26516504294495532165031663578931839*(7.*y2-3.)*y*(1. + y);
    }

    inline double leg_f17_dx_1(double x, double y)
    {
      return -leg_f17_dx_0(x, y);
    }

    inline double leg_f17_dy_0(double x, double y)
    {
      double y2 = y*y;
      return -.26516504294495532165031663578931839*(x + 1.)*(-3. + (-6. + (21. + 28.*y)*y)*y);
    }

    inline double leg_f17_dy_1(double x, double y)
    {
      return -leg_f17_dy_0(x, y);
    }

    inline double leg_f17_dxx_0(double x, double y)
    {
      return 0.;
    }

    inline double leg_f17_dxx_1(double x, double y)
    {
      return -leg_f17_dxx_0(x, y);
    }

    inline double leg_f17_dyy_0(double x, double y)
    {
      return -1.5909902576697319299018998147359104*(x + 1.)*(-1. + (7. + 14.*y)*y);
    }

    inline double leg_f17_dyy_1(double x, double y)
    {
      return -leg_f17_dyy_0(x, y);
    }

    inline double leg_f17_dxy_0(double x, double y)
    {
      return .79549512883486596495094990736795518-.26516504294495532165031663578931839*(-6. + (21. + 28.*y)*y)*y;
    }

    inline double leg_f17_dxy_1(double x, double y)
    {
      return -leg_f17_dxy_0(x, y);
    }

    // number 18
    inline double leg_f18_0(double x, double y)
    {
      double y2 = y*y;
      return -.26516504294495532165031663578931839*(x + y)*(7.*y2-3.)*y*(1. + y);
    }

    // number 18
    inline double leg_f18_1(double x, double y)
    {
      return -leg_f18_0(x, y);
    }

    inline double leg_f18_dx_0(double x, double y)
    {
      double y2 = y*y;
      return -.26516504294495532165031663578931839*(7.*y2-3.)*y*(1. + y);
    }

    inline double leg_f18_dx_1(double x, double y)
    {
      return -leg_f18_dx_0(x, y);
    }

    inline double leg_f18_dy_0(double x, double y)
    {
      double y2 = y*y;
      return -.26516504294495532165031663578931839*(-6. + (-9. + (28. + 35.*y)*y)*y)*y-.26516504294495532165031663578931839*(-3. + (-6. + (21. + 28.*y)*y)*y)*x;
    }

    inline double leg_f18_dy_1(double x, double y)
    {
      return -leg_f18_dy_0(x, y);
    }

    inline double leg_f18_dxx_0(double x, double y)
    {
      return 0.;
    }

    inline double leg_f18_dxx_1(double x, double y)
    {
      return -leg_f18_dxx_0(x, y);
    }

    inline double leg_f18_dyy_0(double x, double y)
    {
      return 1.5909902576697319299018998147359104-.53033008588991064330063327157863679*(-9. + (42. + 70.*y)*y)*y-
          .53033008588991064330063327157863679*(-3. + (21. + 42.*y)*y)*x;
    }

    inline double leg_f18_dyy_1(double x, double y)
    {
      return -leg_f18_dyy_0(x, y);
    }

    inline double leg_f18_dxy_0(double x, double y)
    {
      return .79549512883486596495094990736795518-.26516504294495532165031663578931839*(-6. + (21. + 28.*y)*y)*y;
    }

    inline double leg_f18_dxy_1(double x, double y)
    {
      return -leg_f18_dxy_0(x, y);
    }

    // Bubble functions, order 5

    // number 19
    inline double leg_f19(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      return 1.1023963796102460793756732306830252*(.61237243569579452454932101867647285*a2-.61237243569579452454932101867647285*b2)*
          (-.25000000000000000000000000000000000 + (.25000000000000000000000000000000000 + (2.2500000000000000000000000000000000+
          1.7500000000000000000000000000000000*y)*y)*y);
    }

    inline double leg_f19_dx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return 2.7003086243366084295691530983699986*a*(-.25000000000000000000000000000000000 + (.25000000000000000000000000000000000+
          (2.2500000000000000000000000000000000 + 1.7500000000000000000000000000000000*y)*y)*y);
    }

    inline double leg_f19_dy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return (-.33753857804207605369614413729624983 + 1.3501543121683042147845765491849993*(.25000000000000000000000000000000000+
          (2.2500000000000000000000000000000000 + 1.7500000000000000000000000000000000*y)*y)*y-.67507715608415210739228827459249967*
          (.25000000000000000000000000000000000 + (4.5000000000000000000000000000000000 + 5.2500000000000000000000000000000000*y)*y)*b)*b+
          (-.33753857804207605369614413729624983 + 1.3501543121683042147845765491849993*(.25000000000000000000000000000000000+
          (2.2500000000000000000000000000000000 + 1.7500000000000000000000000000000000*y)*y)*y + .67507715608415210739228827459249967*
          (.25000000000000000000000000000000000 + (4.5000000000000000000000000000000000 + 5.2500000000000000000000000000000000*y)*y)*a)*a;
    }

    inline double leg_f19_dxx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return -1.3501543121683042147845765491849993 + 5.4006172486732168591383061967399971*(.25000000000000000000000000000000000+
          (2.2500000000000000000000000000000000 + 1.7500000000000000000000000000000000*y)*y)*y;
    }

    inline double leg_f19_dyy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return (.67507715608415210739228827459249967 + 2.7003086243366084295691530983699986*(4.5000000000000000000000000000000000+
          5.2500000000000000000000000000000000*y)*y-.67507715608415210739228827459249967*(10.500000000000000000000000000000000*y+
          4.5000000000000000000000000000000000)*b)*b + (.67507715608415210739228827459249967 + 2.7003086243366084295691530983699986*
          (4.5000000000000000000000000000000000 + 5.2500000000000000000000000000000000*y)*y + .67507715608415210739228827459249967*
          (10.500000000000000000000000000000000*y + 4.5000000000000000000000000000000000)*a)*a;
    }

    inline double leg_f19_dxy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return -.67507715608415210739228827459249967 + 2.7003086243366084295691530983699986*(.25000000000000000000000000000000000+
          (2.2500000000000000000000000000000000 + 1.7500000000000000000000000000000000*y)*y)*y + 2.7003086243366084295691530983699986*
          (.25000000000000000000000000000000000 + (4.5000000000000000000000000000000000 + 5.2500000000000000000000000000000000*y)*y)*a;
    }

    // number 20
    inline double leg_f20(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      return .72944814386945601069179830653602616*(-.79056941504209483299972338610817962*b2 + .79056941504209483299972338610817962*a2)*a*(.75000000000000000000000000000000000 + (2.5000000000000000000000000000000000 + 1.7500000000000000000000000000000000*y)*y);
    }

    inline double leg_f20_dx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      return .72944814386945601069179830653602616*(4.7434164902525689979983403166490778*a2-1.5811388300841896659994467722163592*b2)*(.75000000000000000000000000000000000 + (2.5000000000000000000000000000000000 + 1.7500000000000000000000000000000000*y)*y);
    }

    inline double leg_f20_dy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      return -.57667939240241767253899168157713749*b2*(.75000000000000000000000000000000000 + (2.5000000000000000000000000000000000 + 1.7500000000000000000000000000000000*y)*y) + ((.86501908860362650880848752236570622 + 1.1533587848048353450779833631542750*(2.5000000000000000000000000000000000 + 1.7500000000000000000000000000000000*y)*y-.57667939240241767253899168157713749*(3.5000000000000000000000000000000000*y + 2.5000000000000000000000000000000000)*b)*b + (1.2975286329054397632127312835485594 + 1.7300381772072530176169750447314125*(2.5000000000000000000000000000000000 + 1.7500000000000000000000000000000000*y)*y + .57667939240241767253899168157713749*(3.5000000000000000000000000000000000*y + 2.5000000000000000000000000000000000)*a)*a)*a;
    }

    inline double leg_f20_dxx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return 13.840305417658024140935800357851299*a*(.75000000000000000000000000000000000 + (2.5000000000000000000000000000000000 + 1.7500000000000000000000000000000000*y)*y);
    }

    inline double leg_f20_dyy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return (1.7300381772072530176169750447314124 + 2.3067175696096706901559667263085499*(2.5000000000000000000000000000000000 + 1.7500000000000000000000000000000000*y)*y-1.1533587848048353450779833631542750*(3.5000000000000000000000000000000000*y + 2.5000000000000000000000000000000000)*b)*b + (1.7300381772072530176169750447314124 + 2.3067175696096706901559667263085499*(2.5000000000000000000000000000000000 + 1.7500000000000000000000000000000000*y)*y + (8.0735114936338474155458835420799248*y + 5.7667939240241767253899168157713747-2.0183778734084618538864708855199811*b)*b + (12.110267240450771123318825313119887*y + 8.6501908860362650880848752236570620 + 2.0183778734084618538864708855199811*a)*a)*a;
    }

    inline double leg_f20_dxy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return (1.7300381772072530176169750447314124 + 2.3067175696096706901559667263085499*(2.5000000000000000000000000000000000 + 1.7500000000000000000000000000000000*y)*y-1.1533587848048353450779833631542750*(3.5000000000000000000000000000000000*y + 2.5000000000000000000000000000000000)*b)*b + (5.1901145316217590528509251341942374 + 6.9201527088290120704679001789256498*(2.5000000000000000000000000000000000 + 1.7500000000000000000000000000000000*y)*y + 3.4600763544145060352339500894628249*(3.5000000000000000000000000000000000*y + 2.5000000000000000000000000000000000)*a)*a;
    }

    // number 21
    inline double leg_f21(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b4 = b2*b2;
      return .33145630368119415206289579473664799*(.23385358667337133659898429576978433*b4 + (-1.4031215200402280195939057746187060*b2 + 1.1692679333668566829949214788489217*a2)*a2)*(1. + y);
    }

    inline double leg_f21_dx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      return .33145630368119415206289579473664799*(-5.6124860801609120783756230984748240*b2 + 9.3541434669348534639593718307913732*a2)*a*(1. + y);
    }

    inline double leg_f21_dy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      return (-.31004898176538170982440809612960083-.31004898176538170982440809612960083*y + .77512245441345427456102024032400208e-1*b)*b3 + (-.93014694529614512947322428838880250*b2*(1. + y) + ((.93014694529614512947322428838880250 + .93014694529614512947322428838880250*y-.46507347264807256473661214419440125*b)*b + (1.5502449088269085491220404806480042 + 1.5502449088269085491220404806480042*y + .38756122720672713728051012016200104*a)*a)*a)*a;
    }

    inline double leg_f21_dxx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      return .33145630368119415206289579473664799*(56.124860801609120783756230984748240*a2-11.224972160321824156751246196949648*b2)*(1. + y);
    }

    inline double leg_f21_dyy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      return -.62009796353076341964881619225920167*b3 + ((3.7205877811845805178928971535552100 + 3.7205877811845805178928971535552100*y-1.8602938905922902589464485767776050*b)*b + (3.7205877811845805178928971535552100 + 3.7205877811845805178928971535552100*y + 1.8602938905922902589464485767776050*b + 3.1004898176538170982440809612960083*a)*a)*a;
    }

    inline double leg_f21_dxy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      return -1.8602938905922902589464485767776050*b2*(1. + y) + ((3.7205877811845805178928971535552100 + 3.7205877811845805178928971535552100*y-1.8602938905922902589464485767776050*b)*b + (9.3014694529614512947322428838880250 + 9.3014694529614512947322428838880250*y + 3.1004898176538170982440809612960083*a)*a)*a;
    }

    // ORDER 6

    // Edge functions, order 6

    // number 22
    inline double leg_f22(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b4 = b2*b2;
      double b6 = b4*b2;
      return -.22902420702262839621902490788791339e-2*b6 + (.34353631053394259432853736183187009e-1*b4 + (-.80158472457919938676658717760769688e-1*b2 + .48095083474751963205995230656461813e-1*a2)*a2)*a2;
    }

    inline double leg_f22_dx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      return (.13741452421357703773141494473274804*b4 + (-.64126777966335950941326974208615750*b2 + .57714100169702355847194276787754175*a2)*a2)*a;
    }

    inline double leg_f22_dy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      return .13741452421357703773141494473274804e-1*b5 + (.68707262106788518865707472366374018e-1*b4 + (-.13741452421357703773141494473274804*b3 + (-.32063388983167975470663487104307875*b2 + (.16031694491583987735331743552153938*b + .28857050084851177923597138393877088*a)*a)*a)*a)*a;
    }

    inline double leg_f22_dxx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      return .27482904842715407546282988946549607*b4 + (-3.8476066779801570564796184525169450*b2 + 5.7714100169702355847194276787754175*a2)*a2;
    }

    inline double leg_f22_dyy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      return (-.54965809685430815092565977893099214*b3 + (-.54965809685430815092565977893099214*b2 + (1.2825355593267190188265394841723150*b + 1.2825355593267190188265394841723150*a)*a)*a)*a;
    }

    inline double leg_f22_dxy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      return .13741452421357703773141494473274804*b4 + (-.54965809685430815092565977893099214*b3 + (-1.9238033389900785282398092262584725*b2 + (1.2825355593267190188265394841723150*b + 2.8857050084851177923597138393877088*a)*a)*a)*a;
    }

    // number 23
    inline double leg_f23(double x, double y)
    {
      double y2 = y*y;
      return -.14657549249448217358017594104826457*(x + 1.)*(1. + (-14. + 21.*y2)*y2)*(1. + y);
    }

    inline double leg_f23_dx(double x, double y)
    {
      double y2 = y*y;
      return -.14657549249448217358017594104826457*(1. + (-14. + 21.*y2)*y2)*(1. + y);
    }

    inline double leg_f23_dy(double x, double y)
    {
      double y2 = y*y;
      return -.14657549249448217358017594104826457*(x + 1.)*(1. + (-28. + (-42. + (84. + 105.*y)*y)*y)*y);
    }

    inline double leg_f23_dxx(double x, double y)
    {
      return 0.;
    }

    inline double leg_f23_dyy(double x, double y)
    {
      return -4.1041137898455008602449263493514080*(x + 1.)*(-1. + (-3. + (9. + 15.*y)*y)*y);
    }

    inline double leg_f23_dxy(double x, double y)
    {
      return -.14657549249448217358017594104826457-.14657549249448217358017594104826457*(-28. + (-42. + (84. + 105.*y)*y)*y)*y;
    }

    // number 24
    inline double leg_f24(double x, double y)
    {
      double y2 = y*y;
      return .14657549249448217358017594104826457*(x + y)*(1. + (-14. + 21.*y2)*y2)*(1. + y);
    }

    inline double leg_f24_dx(double x, double y)
    {
      double y2 = y*y;
      return .14657549249448217358017594104826457*(1. + (-14. + 21.*y2)*y2)*(1. + y);
    }

    inline double leg_f24_dy(double x, double y)
    {
      double y2 = y*y;
      return .14657549249448217358017594104826457 + .14657549249448217358017594104826457*(2. + (-42. + (-56. + (105. + 126.*y)*y)*y)*y)*y + .14657549249448217358017594104826457*(1. + (-28. + (-42. + (84. + 105.*y)*y)*y)*y)*x;
    }

    inline double leg_f24_dxx(double x, double y)
    {
      return 0.;
    }

    inline double leg_f24_dyy(double x, double y)
    {
      return .29315098498896434716035188209652914 + .29315098498896434716035188209652914*(-42. + (-84. + (210. + 315.*y)*y)*y)*y + .29315098498896434716035188209652914*(-14. + (-42. + (126. + 210.*y)*y)*y)*x;
    }

    inline double leg_f24_dxy(double x, double y)
    {
      return .14657549249448217358017594104826457 + .14657549249448217358017594104826457*(-28. + (-42. + (84. + 105.*y)*y)*y)*y;
    }

    // Bubble functions, order 6

    // number 25
    inline double leg_f25(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double y2 = y*y;
      return 1.3055824196677337863296844245952112*(.61237243569579452454932101867647285*a2-.61237243569579452454932101867647285*b2)*(-.12500000000000000000000000000000000 + (-1. + (3.5000000000000000000000000000000000 + 2.6250000000000000000000000000000000*y)*y2)*y);
    }

    inline double leg_f25_dx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double y2 = y*y;
      return 3.1980107453341565144765659865075906*a*(-.12500000000000000000000000000000000 + (-1. + (3.5000000000000000000000000000000000 + 2.6250000000000000000000000000000000*y)*y2)*y);
    }

    inline double leg_f25_dy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double y2 = y*y;
      return (-.19987567158338478215478537415672443 + 1.5990053726670782572382829932537954*(-1. + (3.5000000000000000000000000000000000 + 2.6250000000000000000000000000000000*y)*y2)*y-.79950268633353912861914149662689767*(-1. + (10.500000000000000000000000000000000 + 10.500000000000000000000000000000000*y)*y2)*b)*b + (-.19987567158338478215478537415672443 + 1.5990053726670782572382829932537954*(-1. + (3.5000000000000000000000000000000000 + 2.6250000000000000000000000000000000*y)*y2)*y + .79950268633353912861914149662689767*(-1. + (10.500000000000000000000000000000000 + 10.500000000000000000000000000000000*y)*y2)*a)*a;
    }

    inline double leg_f25_dxx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double y2 = y*y;
      return -.79950268633353912861914149662689769 + 6.3960214906683130289531319730151814*(-1. + (3.5000000000000000000000000000000000 + 2.6250000000000000000000000000000000*y)*y2)*y;
    }

    inline double leg_f25_dyy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double y2 = y*y;
      return (-3.1980107453341565144765659865075906 + 3.1980107453341565144765659865075906*(10.500000000000000000000000000000000 + 10.500000000000000000000000000000000*y)*y2-.79950268633353912861914149662689767*(21. + 31.500000000000000000000000000000000*y)*y*b)*b + (-3.1980107453341565144765659865075906 + 3.1980107453341565144765659865075906*(10.500000000000000000000000000000000 + 10.500000000000000000000000000000000*y)*y2 + .79950268633353912861914149662689767*(21. + 31.500000000000000000000000000000000*y)*y*a)*a;
    }

    inline double leg_f25_dxy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double y2 = y*y;
      return -.39975134316676956430957074831344884 + 3.1980107453341565144765659865075906*(-1. + (3.5000000000000000000000000000000000 + 2.6250000000000000000000000000000000*y)*y2)*y + 3.1980107453341565144765659865075906*(-1. + (10.500000000000000000000000000000000 + 10.500000000000000000000000000000000*y)*y2)*a;
    }

    // number 26
    inline double leg_f26(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      return .99215674164922147143810590761472265*(-.79056941504209483299972338610817962*b2 + .79056941504209483299972338610817962*a2)*a*(2. + (5. + 3.*y)*y)*y;
    }

    inline double leg_f26_dx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      return .99215674164922147143810590761472265*(4.7434164902525689979983403166490778*a2-1.5811388300841896659994467722163592*b2)*(2. + (5. + 3.*y)*y)*y;
    }

    inline double leg_f26_dy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      return -.78436877487569582622953627417361325*b2*(2. + (5. + 3.*y)*y)*y + ((1.5687375497513916524590725483472265*(2. + (5. + 3.*y)*y)*y-.78436877487569582622953627417361325*(2. + (10. + 9.*y)*y)*b)*b + (2.3531063246270874786886088225208398*(2. + (5. + 3.*y)*y)*y + .78436877487569582622953627417361325*(2. + (10. + 9.*y)*y)*a)*a)*a;
    }

    inline double leg_f26_dxx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return 18.824850597016699829508870580166718*a*(2. + (5. + 3.*y)*y)*y;
    }

    inline double leg_f26_dyy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return (3.1374750995027833049181450966944530*(2. + (5. + 3.*y)*y)*y-1.5687375497513916524590725483472265*(2. + (10. + 9.*y)*y)*b)*b + (3.1374750995027833049181450966944530*(2. + (5. + 3.*y)*y)*y + (6.2749501990055666098362901933889059 + 3.1374750995027833049181450966944530*(10. + 9.*y)*y-.78436877487569582622953627417361325*(18.*y + 10.)*b)*b + (9.4124252985083499147544352900833588 + 4.7062126492541749573772176450416795*(10. + 9.*y)*y + .78436877487569582622953627417361325*(18.*y + 10.)*a)*a)*a;
    }

    inline double leg_f26_dxy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return (3.1374750995027833049181450966944530*(2. + (5. + 3.*y)*y)*y-1.5687375497513916524590725483472265*(2. + (10. + 9.*y)*y)*b)*b + (9.4124252985083499147544352900833590*(2. + (5. + 3.*y)*y)*y + 4.7062126492541749573772176450416795*(2. + (10. + 9.*y)*y)*a)*a;
    }

    // number 27
    inline double leg_f27(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b4 = b2*b2;
      return .50111482858579572128353777664953786*(.23385358667337133659898429576978433*b4 + (-1.4031215200402280195939057746187060*b2 + 1.1692679333668566829949214788489217*a2)*a2)*(1.2500000000000000000000000000000000 + (3.5000000000000000000000000000000000 + 2.2500000000000000000000000000000000*y)*y);
    }

    inline double leg_f27_dx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      return .50111482858579572128353777664953786*(-5.6124860801609120783756230984748240*b2 + 9.3541434669348534639593718307913732*a2)*a*(1.2500000000000000000000000000000000 + (3.5000000000000000000000000000000000 + 2.2500000000000000000000000000000000*y)*y);
    }

    inline double leg_f27_dy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      return (-.58593750000000000000000000000000000-.46875000000000000000000000000000000*(3.5000000000000000000000000000000000 + 2.2500000000000000000000000000000000*y)*y + .11718750000000000000000000000000000*(4.5000000000000000000000000000000000*y + 3.5000000000000000000000000000000000)*b)*b3 + (-1.4062500000000000000000000000000000*b2*(1.2500000000000000000000000000000000 + (3.5000000000000000000000000000000000 + 2.2500000000000000000000000000000000*y)*y) + ((1.7578125000000000000000000000000000 + 1.4062500000000000000000000000000000*(3.5000000000000000000000000000000000 + 2.2500000000000000000000000000000000*y)*y-.70312500000000000000000000000000000*(4.5000000000000000000000000000000000*y + 3.5000000000000000000000000000000000)*b)*b + (2.9296875000000000000000000000000000 + 2.3437500000000000000000000000000000*(3.5000000000000000000000000000000000 + 2.2500000000000000000000000000000000*y)*y + .58593750000000000000000000000000000*(4.5000000000000000000000000000000000*y + 3.5000000000000000000000000000000000)*a)*a)*a)*a;
    }

    inline double leg_f27_dxx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      return .50111482858579572128353777664953786*(56.124860801609120783756230984748240*a2-11.224972160321824156751246196949648*b2)*(1.2500000000000000000000000000000000 + (3.5000000000000000000000000000000000 + 2.2500000000000000000000000000000000*y)*y);
    }

    inline double leg_f27_dyy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      return (-4.2187500000000000000000000000000000*y-3.2812500000000000000000000000000000 + .52734375000000000000000000000000000*b)*b3 + ((7.0312500000000000000000000000000000 + 5.6250000000000000000000000000000000*(3.5000000000000000000000000000000000 + 2.2500000000000000000000000000000000*y)*y-2.8125000000000000000000000000000000*(4.5000000000000000000000000000000000*y + 3.5000000000000000000000000000000000)*b)*b + (7.0312500000000000000000000000000000 + 5.6250000000000000000000000000000000*(3.5000000000000000000000000000000000 + 2.2500000000000000000000000000000000*y)*y + (12.656250000000000000000000000000000*y + 9.8437500000000000000000000000000000-3.1640625000000000000000000000000000*b)*b + (21.093750000000000000000000000000000*y + 16.406250000000000000000000000000000 + 2.6367187500000000000000000000000000*a)*a)*a)*a;
    }

    inline double leg_f27_dxy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      return -2.8125000000000000000000000000000000*b2*(1.2500000000000000000000000000000000 + (3.5000000000000000000000000000000000 + 2.2500000000000000000000000000000000*y)*y) + ((7.0312500000000000000000000000000000 + 5.6250000000000000000000000000000000*(3.5000000000000000000000000000000000 + 2.2500000000000000000000000000000000*y)*y-2.8125000000000000000000000000000000*(4.5000000000000000000000000000000000*y + 3.5000000000000000000000000000000000)*b)*b + (17.578125000000000000000000000000000 + 14.062500000000000000000000000000000*(3.5000000000000000000000000000000000 + 2.2500000000000000000000000000000000*y)*y + 4.6875000000000000000000000000000000*(4.5000000000000000000000000000000000*y + 3.5000000000000000000000000000000000)*a)*a)*a;
    }

    // number 28
    inline double leg_f28(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b4 = b2*b2;
      return .22642776165920997771435412745577769*(.79549512883486596495094990736795518*b4 + (-2.6516504294495532165031663578931839*b2 + 1.8561553006146872515522164505252288*a2)*a2)*a*(1. + y);
    }

    inline double leg_f28_dx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      return .22642776165920997771435412745577769*(1.5909902576697319299018998147359104*b4 + (-15.909902576697319299018998147359104*b2 + 18.561553006146872515522164505252288*a2)*a2)*(1. + y);
    }

    inline double leg_f28_dy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      return .18012218143288356533005732055603413*b4*(1. + y) + ((-.72048872573153426132022928222413650-.72048872573153426132022928222413650*y + .18012218143288356533005732055603413*b)*b3 + (-1.8012218143288356533005732055603413*b2*(1. + y) + ((1.2008145428858904355337154703735608 + 1.2008145428858904355337154703735608*y-.60040727144294521776685773518678042*b)*b + (2.1014254500503082621840020731537315 + 2.1014254500503082621840020731537315*y + .42028509001006165243680041463074629*a)*a)*a)*a)*a;
    }

    inline double leg_f28_dxx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      return .22642776165920997771435412745577769*(-63.639610306789277196075992589436414*b2 + 148.49242404917498012417731604201830*a2)*a*(1. + y);
    }

    inline double leg_f28_dyy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      return (-1.4409774514630685226404585644482730-1.4409774514630685226404585644482730*y + .36024436286576713066011464111206825*b)*b3 + ((-1.4409774514630685226404585644482730-1.4409774514630685226404585644482730*y-1.4409774514630685226404585644482730*b)*b2 + ((7.2048872573153426132022928222413650 + 7.2048872573153426132022928222413650*y-3.6024436286576713066011464111206825*b)*b + (7.2048872573153426132022928222413650 + 7.2048872573153426132022928222413650*y + 2.4016290857717808710674309407471217*b + 4.2028509001006165243680041463074629*a)*a)*a)*a;
    }

    inline double leg_f28_dxy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      return (-1.4409774514630685226404585644482730-1.4409774514630685226404585644482730*y + .36024436286576713066011464111206825*b)*b3 + (-7.2048872573153426132022928222413650*b2*(1. + y) + ((7.2048872573153426132022928222413650 + 7.2048872573153426132022928222413650*y-3.6024436286576713066011464111206825*b)*b + (16.811403600402466097472016585229852 + 16.811403600402466097472016585229852*y + 4.2028509001006165243680041463074629*a)*a)*a)*a;
    }

    // ORDER 7

    // Edge functions, order 7

    // number 29
    inline double leg_f29_0(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b4 = b2*b2;
      double b6 = b4*b2;
      return (-.62243890546786924194680470080844507e-2*b6 + (.43570723382750846936276329056591155e-1*b4 + (-.78427302088951524485297392301864079e-1*b2 + .41080967760879369968489110253357375e-1*a2)*a2)*a2)*a;
    }

    inline double leg_f29_1(double x, double y)
    {
      return -leg_f29_0(x, y);
    }

    inline double leg_f29_dx_0(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      double b6 = b5*b;
      return -.12448778109357384838936094016168901e-1*b6 + (.26142434029650508161765797433954693*b4 + (-.78427302088951524485297392301864079*b2 + .57513354865231117955884754354700324*a2)*a2)*a2;
    }

    inline double leg_f29_dx_1(double x, double y)
    {
      return -leg_f29_dx_0(x, y);
    }

    inline double leg_f29_dy_0(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      double b6 = b5*b;
      return -.62243890546786924194680470080844507e-2*b6 + (.37346334328072154516808282048506704e-1*b5 + (.13071217014825254080882898716977346*b4 + (-.17428289353100338774510531622636462*b3 + (-.39213651044475762242648696150932039*b2 + (.15685460417790304897059478460372816*b + .28756677432615558977942377177350162*a)*a)*a)*a)*a)*a;
    }

    inline double leg_f29_dy_1(double x, double y)
    {
      return -leg_f29_dy_0(x, y);
    }

    inline double leg_f29_dxx_0(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      return (1.0456973611860203264706318973581877*b4 + (-6.2741841671161219588237913841491263*b2 + 6.9016025838277341547061705225640389*a2)*a2)*a;
    }

    inline double leg_f29_dxx_1(double x, double y)
    {
      return -leg_f29_dxx_0(x, y);
    }

    inline double leg_f29_dyy_0(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      return .74692668656144309033616564097013408e-1*b5 + (.74692668656144309033616564097013408e-1*b4 + (-1.0456973611860203264706318973581877*b3 + (-1.0456973611860203264706318973581877*b2 + (1.5685460417790304897059478460372816*b + 1.5685460417790304897059478460372816*a)*a)*a)*a)*a;
    }

    inline double leg_f29_dyy_1(double x, double y)
    {
      return -leg_f29_dyy_0(x, y);
    }

    inline double leg_f29_dxy_0(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      return .74692668656144309033616564097013408e-1*b5 + (.52284868059301016323531594867909386*b4 + (-1.0456973611860203264706318973581877*b3 + (-3.1370920835580609794118956920745631*b2 + (1.5685460417790304897059478460372816*b + 3.4508012919138670773530852612820195*a)*a)*a)*a)*a;
    }

    inline double leg_f29_dxy_1(double x, double y)
    {
      return -leg_f29_dxy_0(x, y);
    }

    // number 30
    inline double leg_f30_0(double x, double y)
    {
      double y2 = y*y;
      return -.15934435979977452593838200340696194*(x + 1.)*(5. + (-30. + 33.*y2)*y2)*y*(1. + y);
    }

    // number 30
    inline double leg_f30_1(double x, double y)
    {
      return -leg_f30_0(x, y);
    }

    inline double leg_f30_dx_0(double x, double y)
    {
      double y2 = y*y;
      return -.15934435979977452593838200340696194*(5. + (-30. + 33.*y2)*y2)*y*(1. + y);
    }

    inline double leg_f30_dx_1(double x, double y)
    {
      return -leg_f30_dx_0(x, y);
    }

    inline double leg_f30_dy_0(double x, double y)
    {
      double y2 = y*y;
      return -.15934435979977452593838200340696194*(x + 1.)*(5. + (10. + (-90. + (-120. + (165. + 198.*y)*y)*y)*y)*y);
    }

    inline double leg_f30_dy_1(double x, double y)
    {
      return -leg_f30_dy_0(x, y);
    }

    inline double leg_f30_dxx_0(double x, double y)
    {
      return 0.;
    }

    inline double leg_f30_dxx_1(double x, double y)
    {
      return -leg_f30_dxx_0(x, y);
    }

    inline double leg_f30_dyy_0(double x, double y)
    {
      return -1.5934435979977452593838200340696194*(x + 1.)*(1. + (-18. + (-36. + (66. + 99.*y)*y)*y)*y);
    }

    inline double leg_f30_dyy_1(double x, double y)
    {
      return -leg_f30_dyy_0(x, y);
    }

    inline double leg_f30_dxy_0(double x, double y)
    {
      return -.79672179899887262969191001703480969-.15934435979977452593838200340696194*(10. + (-90. + (-120. + (165. + 198.*y)*y)*y)*y)*y;
    }

    inline double leg_f30_dxy_1(double x, double y)
    {
      return -leg_f30_dxy_0(x, y);
    }

    // number 31
    inline double leg_f31_0(double x, double y)
    {
      double y2 = y*y;
      return -.15934435979977452593838200340696194*(x + y)*(5. + (-30. + 33.*y2)*y2)*y*(1. + y);
    }

    // number 31
    inline double leg_f31_1(double x, double y)
    {
      return -leg_f31_0(x, y);
    }

    inline double leg_f31_dx_0(double x, double y)
    {
      double y2 = y*y;
      return -.15934435979977452593838200340696194*(5. + (-30. + 33.*y2)*y2)*y*(1. + y);
    }

    inline double leg_f31_dx_1(double x, double y)
    {
      return -leg_f31_dx_0(x, y);
    }

    inline double leg_f31_dy_0(double x, double y)
    {
      double y2 = y*y;
      return -.15934435979977452593838200340696194*(10. + (15. + (-120. + (-150. + (198. + 231.*y)*y)*y)*y)*y)*y-.15934435979977452593838200340696194*(5. + (10. + (-90. + (-120. + (165. + 198.*y)*y)*y)*y)*y)*x;
    }

    inline double leg_f31_dy_1(double x, double y)
    {
      return -leg_f31_dy_0(x, y);
    }

    inline double leg_f31_dxx_0(double x, double y)
    {
      return 0.;
    }

    inline double leg_f31_dxx_1(double x, double y)
    {
      return -leg_f31_dxx_0(x, y);
    }

    inline double leg_f31_dyy_0(double x, double y)
    {
      return -1.5934435979977452593838200340696194-.31868871959954905187676400681392388*(15. + (-180. + (-300. + (495. + 693.*y)*y)*y)*y)*y-.31868871959954905187676400681392388*(5. + (-90. + (-180. + (330. + 495.*y)*y)*y)*y)*x;
    }

    inline double leg_f31_dyy_1(double x, double y)
    {
      return -leg_f31_dyy_0(x, y);
    }

    inline double leg_f31_dxy_0(double x, double y)
    {
      return -.79672179899887262969191001703480969-.15934435979977452593838200340696194*(10. + (-90. + (-120. + (165. + 198.*y)*y)*y)*y)*y;
    }

    inline double leg_f31_dxy_1(double x, double y)
    {
      return -leg_f31_dxy_0(x, y);
    }

    // Bubble functions, order 7

    // number 32
    inline double leg_f32(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      return 1.4752421108802058983765559207693054*(.61237243569579452454932101867647285*a2-.61237243569579452454932101867647285*b2)*(.12500000000000000000000000000000000 + (-.37500000000000000000000000000000000 + (-2.7500000000000000000000000000000000 + (-.75000000000000000000000000000000000 + (5.6250000000000000000000000000000000 + 4.1250000000000000000000000000000000*y)*y)*y)*y)*y);
    }

    inline double leg_f32_dx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return 3.6135904187228682497009747173648780*a*(.12500000000000000000000000000000000 + (-.37500000000000000000000000000000000 + (-2.7500000000000000000000000000000000 + (-.75000000000000000000000000000000000 + (5.6250000000000000000000000000000000 + 4.1250000000000000000000000000000000*y)*y)*y)*y)*y);
    }

    inline double leg_f32_dy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return (.22584940117017926560631091983530488 + 1.8067952093614341248504873586824390*(-.37500000000000000000000000000000000 + (-2.7500000000000000000000000000000000 + (-.75000000000000000000000000000000000 + (5.6250000000000000000000000000000000 + 4.1250000000000000000000000000000000*y)*y)*y)*y)*y-.90339760468071706242524367934121949*(-.37500000000000000000000000000000000 + (-5.5000000000000000000000000000000000 + (-2.2500000000000000000000000000000000 + (22.500000000000000000000000000000000 + 20.625000000000000000000000000000000*y)*y)*y)*y)*b)*b + (.22584940117017926560631091983530488 + 1.8067952093614341248504873586824390*(-.37500000000000000000000000000000000 + (-2.7500000000000000000000000000000000 + (-.75000000000000000000000000000000000 + (5.6250000000000000000000000000000000 + 4.1250000000000000000000000000000000*y)*y)*y)*y)*y + .90339760468071706242524367934121949*(-.37500000000000000000000000000000000 + (-5.5000000000000000000000000000000000 + (-2.2500000000000000000000000000000000 + (22.500000000000000000000000000000000 + 20.625000000000000000000000000000000*y)*y)*y)*y)*a)*a;
    }

    inline double leg_f32_dxx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return .90339760468071706242524367934121949 + 7.2271808374457364994019494347297560*(-.37500000000000000000000000000000000 + (-2.7500000000000000000000000000000000 + (-.75000000000000000000000000000000000 + (5.6250000000000000000000000000000000 + 4.1250000000000000000000000000000000*y)*y)*y)*y)*y;
    }

    inline double leg_f32_dyy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return (-1.3550964070210755936378655190118293 + 3.6135904187228682497009747173648780*(-5.5000000000000000000000000000000000 + (-2.2500000000000000000000000000000000 + (22.500000000000000000000000000000000 + 20.625000000000000000000000000000000*y)*y)*y)*y-.90339760468071706242524367934121949*(-5.5000000000000000000000000000000000 + (-4.5000000000000000000000000000000000 + (67.500000000000000000000000000000000 + 82.500000000000000000000000000000000*y)*y)*y)*b)*b + (-1.3550964070210755936378655190118293 + 3.6135904187228682497009747173648780*(-5.5000000000000000000000000000000000 + (-2.2500000000000000000000000000000000 + (22.500000000000000000000000000000000 + 20.625000000000000000000000000000000*y)*y)*y)*y + .90339760468071706242524367934121949*(-5.5000000000000000000000000000000000 + (-4.5000000000000000000000000000000000 + (67.500000000000000000000000000000000 + 82.500000000000000000000000000000000*y)*y)*y)*a)*a;
    }

    inline double leg_f32_dxy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return .45169880234035853121262183967060975 + 3.6135904187228682497009747173648780*(-.37500000000000000000000000000000000 + (-2.7500000000000000000000000000000000 + (-.75000000000000000000000000000000000 + (5.6250000000000000000000000000000000 + 4.1250000000000000000000000000000000*y)*y)*y)*y)*y + 3.6135904187228682497009747173648780*(-.37500000000000000000000000000000000 + (-5.5000000000000000000000000000000000 + (-2.2500000000000000000000000000000000 + (22.500000000000000000000000000000000 + 20.625000000000000000000000000000000*y)*y)*y)*y)*a;
    }

    // number 33
    inline double leg_f33(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      return 1.2208540490103741136884064212095202*(-.79056941504209483299972338610817962*b2 + .79056941504209483299972338610817962*a2)*a*(-.34375000000000000000000000000000000 + (-.62500000000000000000000000000000000 + (3.9375000000000000000000000000000000 + (9.3750000000000000000000000000000000 + 5.1562500000000000000000000000000000*y)*y)*y)*y);
    }

    inline double leg_f33_dx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      return 1.2208540490103741136884064212095202*(4.7434164902525689979983403166490778*a2-1.5811388300841896659994467722163592*b2)*(-.34375000000000000000000000000000000 + (-.62500000000000000000000000000000000 + (3.9375000000000000000000000000000000 + (9.3750000000000000000000000000000000 + 5.1562500000000000000000000000000000*y)*y)*y)*y);
    }

    inline double leg_f33_dy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      return -.96516987137790443929991476509124903*b2*(-.34375000000000000000000000000000000 + (-.62500000000000000000000000000000000 + (3.9375000000000000000000000000000000 + (9.3750000000000000000000000000000000 + 5.1562500000000000000000000000000000*y)*y)*y)*y) + ((-.66355428657230930201869140100023374 + 1.9303397427558088785998295301824982*(-.62500000000000000000000000000000000 + (3.9375000000000000000000000000000000 + (9.3750000000000000000000000000000000 + 5.1562500000000000000000000000000000*y)*y)*y)*y-.96516987137790443929991476509124903*(-.62500000000000000000000000000000000 + (7.8750000000000000000000000000000000 + (28.125000000000000000000000000000000 + 20.625000000000000000000000000000000*y)*y)*y)*b)*b + (-.99533142985846395302803710150035059 + 2.8955096141337133178997442952737472*(-.62500000000000000000000000000000000 + (3.9375000000000000000000000000000000 + (9.3750000000000000000000000000000000 + 5.1562500000000000000000000000000000*y)*y)*y)*y + .96516987137790443929991476509124903*(-.62500000000000000000000000000000000 + (7.8750000000000000000000000000000000 + (28.125000000000000000000000000000000 + 20.625000000000000000000000000000000*y)*y)*y)*a)*a)*a;
    }

    inline double leg_f33_dxx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return 23.164076913069706543197954362189978*a*(-.34375000000000000000000000000000000 + (-.62500000000000000000000000000000000 + (3.9375000000000000000000000000000000 + (9.3750000000000000000000000000000000 + 5.1562500000000000000000000000000000*y)*y)*y)*y);
    }

    inline double leg_f33_dyy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return (-1.3271085731446186040373828020004674 + 3.8606794855116177571996590603649962*(-.62500000000000000000000000000000000 + (3.9375000000000000000000000000000000 + (9.3750000000000000000000000000000000 + 5.1562500000000000000000000000000000*y)*y)*y)*y-1.9303397427558088785998295301824982*(-.62500000000000000000000000000000000 + (7.8750000000000000000000000000000000 + (28.125000000000000000000000000000000 + 20.625000000000000000000000000000000*y)*y)*y)*b)*b + (-1.3271085731446186040373828020004674 + 3.8606794855116177571996590603649962*(-.62500000000000000000000000000000000 + (3.9375000000000000000000000000000000 + (9.3750000000000000000000000000000000 + 5.1562500000000000000000000000000000*y)*y)*y)*y + (-2.4129246784447610982497869127281226 + 3.8606794855116177571996590603649962*(7.8750000000000000000000000000000000 + (28.125000000000000000000000000000000 + 20.625000000000000000000000000000000*y)*y)*y-.96516987137790443929991476509124903*(7.8750000000000000000000000000000000 + (56.250000000000000000000000000000000 + 61.875000000000000000000000000000000*y)*y)*b)*b + (-3.6193870176671416473746803690921840 + 5.7910192282674266357994885905474944*(7.8750000000000000000000000000000000 + (28.125000000000000000000000000000000 + 20.625000000000000000000000000000000*y)*y)*y + .96516987137790443929991476509124903*(7.8750000000000000000000000000000000 + (56.250000000000000000000000000000000 + 61.875000000000000000000000000000000*y)*y)*a)*a)*a;
    }

    inline double leg_f33_dxy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return (-1.3271085731446186040373828020004674 + 3.8606794855116177571996590603649962*(-.62500000000000000000000000000000000 + (3.9375000000000000000000000000000000 + (9.3750000000000000000000000000000000 + 5.1562500000000000000000000000000000*y)*y)*y)*y-1.9303397427558088785998295301824982*(-.62500000000000000000000000000000000 + (7.8750000000000000000000000000000000 + (28.125000000000000000000000000000000 + 20.625000000000000000000000000000000*y)*y)*y)*b)*b + (-3.9813257194338558121121484060014024 + 11.582038456534853271598977181094989*(-.62500000000000000000000000000000000 + (3.9375000000000000000000000000000000 + (9.3750000000000000000000000000000000 + 5.1562500000000000000000000000000000*y)*y)*y)*y + 5.7910192282674266357994885905474944*(-.62500000000000000000000000000000000 + (7.8750000000000000000000000000000000 + (28.125000000000000000000000000000000 + 20.625000000000000000000000000000000*y)*y)*y)*a)*a;
    }

    // number 34
    inline double leg_f34(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b4 = b2*b2;
      return .67031968203506862686507519221230759*(.23385358667337133659898429576978433*b4 + (-1.4031215200402280195939057746187060*b2 + 1.1692679333668566829949214788489217*a2)*a2)*(.58333333333333333333333333333333333 + (4.7500000000000000000000000000000000 + (8.7500000000000000000000000000000000 + 4.5833333333333333333333333333333333*y)*y)*y);
    }

    inline double leg_f34_dx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      return .67031968203506862686507519221230759*(-5.6124860801609120783756230984748240*b2 + 9.3541434669348534639593718307913732*a2)*a*(.58333333333333333333333333333333333 + (4.7500000000000000000000000000000000 + (8.7500000000000000000000000000000000 + 4.5833333333333333333333333333333333*y)*y)*y);
    }

    inline double leg_f34_dy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      return (-.36576554434386081824112622088898979-.62702664744661854555621637866683965*(4.7500000000000000000000000000000000 + (8.7500000000000000000000000000000000 + 4.5833333333333333333333333333333333*y)*y)*y + .15675666186165463638905409466670991*(4.7500000000000000000000000000000000 + (17.500000000000000000000000000000000 + 13.750000000000000000000000000000000*y)*y)*b)*b3 + (-1.8810799423398556366686491360005189*b2*(.58333333333333333333333333333333333 + (4.7500000000000000000000000000000000 + (8.7500000000000000000000000000000000 + 4.5833333333333333333333333333333333*y)*y)*y) + ((1.0972966330315824547233786626669694 + 1.8810799423398556366686491360005189*(4.7500000000000000000000000000000000 + (8.7500000000000000000000000000000000 + 4.5833333333333333333333333333333333*y)*y)*y-.94053997116992781833432456800025947*(4.7500000000000000000000000000000000 + (17.500000000000000000000000000000000 + 13.750000000000000000000000000000000*y)*y)*b)*b + (1.8288277217193040912056311044449490 + 3.1351332372330927277810818933341982*(4.7500000000000000000000000000000000 + (8.7500000000000000000000000000000000 + 4.5833333333333333333333333333333333*y)*y)*y + .78378330930827318194527047333354953*(4.7500000000000000000000000000000000 + (17.500000000000000000000000000000000 + 13.750000000000000000000000000000000*y)*y)*a)*a)*a)*a;
    }

    inline double leg_f34_dxx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      return .67031968203506862686507519221230759*(56.124860801609120783756230984748240*a2-11.224972160321824156751246196949648*b2)*(.58333333333333333333333333333333333 + (4.7500000000000000000000000000000000 + (8.7500000000000000000000000000000000 + 4.5833333333333333333333333333333333*y)*y)*y);
    }

    inline double leg_f34_dyy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      return (-5.9567531507428761827840555973349767-1.2540532948932370911124327573336793*(17.500000000000000000000000000000000 + 13.750000000000000000000000000000000*y)*y + .15675666186165463638905409466670991*(27.500000000000000000000000000000000*y + 17.500000000000000000000000000000000)*b)*b3 + ((4.3891865321263298188935146506678776 + 7.5243197693594225466745965440020758*(4.7500000000000000000000000000000000 + (8.7500000000000000000000000000000000 + 4.5833333333333333333333333333333333*y)*y)*y-3.7621598846797112733372982720010379*(4.7500000000000000000000000000000000 + (17.500000000000000000000000000000000 + 13.750000000000000000000000000000000*y)*y)*b)*b + (4.3891865321263298188935146506678776 + 7.5243197693594225466745965440020758*(4.7500000000000000000000000000000000 + (8.7500000000000000000000000000000000 + 4.5833333333333333333333333333333333*y)*y)*y + (17.870259452228628548352166792004930 + 3.7621598846797112733372982720010379*(17.500000000000000000000000000000000 + 13.750000000000000000000000000000000*y)*y-.94053997116992781833432456800025947*(27.500000000000000000000000000000000*y + 17.500000000000000000000000000000000)*b)*b + (29.783765753714380913920277986674883 + 6.2702664744661854555621637866683965*(17.500000000000000000000000000000000 + 13.750000000000000000000000000000000*y)*y + .78378330930827318194527047333354953*(27.500000000000000000000000000000000*y + 17.500000000000000000000000000000000)*a)*a)*a)*a;
    }

    inline double leg_f34_dxy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      return -3.7621598846797112733372982720010379*b2*(.58333333333333333333333333333333333 + (4.7500000000000000000000000000000000 + (8.7500000000000000000000000000000000 + 4.5833333333333333333333333333333333*y)*y)*y) + ((4.3891865321263298188935146506678776 + 7.5243197693594225466745965440020758*(4.7500000000000000000000000000000000 + (8.7500000000000000000000000000000000 + 4.5833333333333333333333333333333333*y)*y)*y-3.7621598846797112733372982720010379*(4.7500000000000000000000000000000000 + (17.500000000000000000000000000000000 + 13.750000000000000000000000000000000*y)*y)*b)*b + (10.972966330315824547233786626669694 + 18.810799423398556366686491360005189*(4.7500000000000000000000000000000000 + (8.7500000000000000000000000000000000 + 4.5833333333333333333333333333333333*y)*y)*y + 6.2702664744661854555621637866683965*(4.7500000000000000000000000000000000 + (17.500000000000000000000000000000000 + 13.750000000000000000000000000000000*y)*y)*a)*a)*a;
    }

    // number 35
    inline double leg_f35(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b4 = b2*b2;
      return .32110057093916629088182291342744104*(.79549512883486596495094990736795518*b4 + (-2.6516504294495532165031663578931839*b2 + 1.8561553006146872515522164505252288*a2)*a2)*a*(1.7500000000000000000000000000000000 + (4.5000000000000000000000000000000000 + 2.7500000000000000000000000000000000*y)*y);
    }

    inline double leg_f35_dx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      return .32110057093916629088182291342744104*(1.5909902576697319299018998147359104*b4 + (-15.909902576697319299018998147359104*b2 + 18.561553006146872515522164505252288*a2)*a2)*(1.7500000000000000000000000000000000 + (4.5000000000000000000000000000000000 + 2.7500000000000000000000000000000000*y)*y);
    }

    inline double leg_f35_dy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      return .25543394004820110678174896037405033*b4*(1.7500000000000000000000000000000000 + (4.5000000000000000000000000000000000 + 2.7500000000000000000000000000000000*y)*y) + ((-1.7880375803374077474722427226183523-1.0217357601928044271269958414962013*(4.5000000000000000000000000000000000 + 2.7500000000000000000000000000000000*y)*y + .25543394004820110678174896037405033*(5.5000000000000000000000000000000000*y + 4.5000000000000000000000000000000000)*b)*b3 + (-2.5543394004820110678174896037405033*b2*(1.7500000000000000000000000000000000 + (4.5000000000000000000000000000000000 + 2.7500000000000000000000000000000000*y)*y) + ((2.9800626338956795791204045376972539 + 1.7028929336546740452116597358270022*(4.5000000000000000000000000000000000 + 2.7500000000000000000000000000000000*y)*y-.85144646682733702260582986791350110*(5.5000000000000000000000000000000000*y + 4.5000000000000000000000000000000000)*b)*b + (5.2151096093174392634607079409701943 + 2.9800626338956795791204045376972539*(4.5000000000000000000000000000000000 + 2.7500000000000000000000000000000000*y)*y + .59601252677913591582408090753945079*(5.5000000000000000000000000000000000*y + 4.5000000000000000000000000000000000)*a)*a)*a)*a)*a;
    }

    inline double leg_f35_dxx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      return .32110057093916629088182291342744104*(-63.639610306789277196075992589436414*b2 + 148.49242404917498012417731604201830*a2)*a*(1.7500000000000000000000000000000000 + (4.5000000000000000000000000000000000 + 2.7500000000000000000000000000000000*y)*y);
    }

    inline double leg_f35_dyy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      return (-3.5760751606748154949444854452367047-2.0434715203856088542539916829924027*(4.5000000000000000000000000000000000 + 2.7500000000000000000000000000000000*y)*y + .51086788009640221356349792074810066*(5.5000000000000000000000000000000000*y + 4.5000000000000000000000000000000000)*b)*b3 + ((-3.5760751606748154949444854452367047-2.0434715203856088542539916829924027*(4.5000000000000000000000000000000000 + 2.7500000000000000000000000000000000*y)*y + (-11.239093362120848698396954256458215*y-9.1956218417352398441429625734658122 + 1.4048866702651060872996192820572768*b)*b)*b2 + ((17.880375803374077474722427226183523 + 10.217357601928044271269958414962013*(4.5000000000000000000000000000000000 + 2.7500000000000000000000000000000000*y)*y-5.1086788009640221356349792074810066*(5.5000000000000000000000000000000000*y + 4.5000000000000000000000000000000000)*b)*b + (17.880375803374077474722427226183523 + 10.217357601928044271269958414962013*(4.5000000000000000000000000000000000 + 2.7500000000000000000000000000000000*y)*y + (18.731822270201414497328257094097024*y + 15.326036402892066406904937622443020-4.6829555675503536243320642735242559*b)*b + (32.780688972852475370324449914669793*y + 26.820563705061116212083640839275286 + 3.2780688972852475370324449914669793*a)*a)*a)*a)*a;
    }

    inline double leg_f35_dxy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      return (-3.5760751606748154949444854452367047-2.0434715203856088542539916829924027*(4.5000000000000000000000000000000000 + 2.7500000000000000000000000000000000*y)*y + .51086788009640221356349792074810066*(5.5000000000000000000000000000000000*y + 4.5000000000000000000000000000000000)*b)*b3 + (-10.217357601928044271269958414962013*b2*(1.7500000000000000000000000000000000 + (4.5000000000000000000000000000000000 + 2.7500000000000000000000000000000000*y)*y) + ((17.880375803374077474722427226183523 + 10.217357601928044271269958414962013*(4.5000000000000000000000000000000000 + 2.7500000000000000000000000000000000*y)*y-5.1086788009640221356349792074810066*(5.5000000000000000000000000000000000*y + 4.5000000000000000000000000000000000)*b)*b + (41.720876874539514107685663527761555 + 23.840501071165436632963236301578031*(4.5000000000000000000000000000000000 + 2.7500000000000000000000000000000000*y)*y + 5.9601252677913591582408090753945079*(5.5000000000000000000000000000000000*y + 4.5000000000000000000000000000000000)*a)*a)*a)*a;
    }

    // number 36
    inline double leg_f36(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b4 = b2*b2;
      double b6 = b4*b2;
      return .14636714058746245794953655752202643*(-.14657549249448217358017594104826457*b6 + (2.1986323874172326037026391157239686*b4 + (-5.1301422373068760753061579366892600*b2 + 3.0780853423841256451836947620135560*a2)*a2)*a2)*(1. + y);
    }

    inline double leg_f36_dx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      return .14636714058746245794953655752202643*(8.7945295496689304148105564628958743*b4 + (-41.041137898455008602449263493514080*b2 + 36.937024108609507742204337144162672*a2)*a2)*a*(1. + y);
    }

    inline double leg_f36_dy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      return (.12872301429969852374331755029475975 + .12872301429969852374331755029475975*y-.21453835716616420623886258382459958e-1*b)*b5 + (.64361507149849261871658775147379875*b4*(1. + y) + ((-1.2872301429969852374331755029475975-1.2872301429969852374331755029475975*y + .32180753574924630935829387573689937*b)*b3 + (-3.0035370003262988873440761735443942*b2*(1. + y) + ((1.5017685001631494436720380867721971 + 1.5017685001631494436720380867721971*y-.75088425008157472183601904338609854*b)*b + (2.7031833002936689986096685561899547 + 2.7031833002936689986096685561899547*y + .45053055004894483310161142603165912*a)*a)*a)*a)*a)*a;
    }

    inline double leg_f36_dxx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      return .14636714058746245794953655752202643*(17.589059099337860829621112925791749*b4 + (-246.24682739073005161469558096108448*b2 + 369.37024108609507742204337144162672*a2)*a2)*(1. + y);
    }

    inline double leg_f36_dyy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      return .25744602859939704748663510058951950*b5 + ((-5.1489205719879409497327020117903900-5.1489205719879409497327020117903900*y + 1.2872301429969852374331755029475975*b)*b3 + ((-5.1489205719879409497327020117903900-5.1489205719879409497327020117903900*y-2.5744602859939704748663510058951950*b)*b2 + ((12.014148001305195549376304694177577 + 12.014148001305195549376304694177577*y-6.0070740006525977746881523470887883*b)*b + (12.014148001305195549376304694177577 + 12.014148001305195549376304694177577*y + 3.0035370003262988873440761735443942*b + 5.4063666005873379972193371123799095*a)*a)*a)*a)*a;
    }

    inline double leg_f36_dxy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      return 1.2872301429969852374331755029475975*b4*(1. + y) + ((-5.1489205719879409497327020117903900-5.1489205719879409497327020117903900*y + 1.2872301429969852374331755029475975*b)*b3 + (-18.021222001957793324064457041266365*b2*(1. + y) + ((12.014148001305195549376304694177577 + 12.014148001305195549376304694177577*y-6.0070740006525977746881523470887883*b)*b + (27.031833002936689986096685561899547 + 27.031833002936689986096685561899547*y + 5.4063666005873379972193371123799095*a)*a)*a)*a)*a;
    }

    // ORDER 8

    // Edge functions, order 8

    // number 37
    inline double leg_f37(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b4 = b2*b2;
      double b6 = b4*b2;
      double b8 = b6*b2;
      return .41787914848721779896314222930969401e-3*b8 + (-.11700616157642098370967982420671432e-1*b6 + (.52652772709389442669355920893021445e-1*b4 + (-.77224066640437849248388683976431453e-1*b2 + .35854030940203287151037603274771746e-1*a2)*a2)*a2)*a2;
    }

    inline double leg_f37_dx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      double b6 = b5*b;
      double b7 = b6*b;
      return (-.46802464630568393483871929682685729e-1*b6 + (.42122218167511554135484736714417156*b4 + (-.92668879968525419098066420771717743*b2 + .57366449504325259441660165239634793*a2)*a2)*a2)*a;
    }

    inline double leg_f37_dy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      double b6 = b5*b;
      double b7 = b6*b;
      return -.33430331878977423917051378344775521e-2*b7 + (-.23401232315284196741935964841342864e-1*b6 + (.70203696945852590225807894524028593e-1*b5 + (.21061109083755777067742368357208578*b4 + (-.21061109083755777067742368357208578*b3 + (-.46334439984262709549033210385858872*b2 + (.15444813328087569849677736795286291*b + .28683224752162629720830082619817397*a)*a)*a)*a)*a)*a)*a;
    }

    inline double leg_f37_dxx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      double b6 = b5*b;
      return -.93604929261136786967743859365371458e-1*b6 + (2.5273330900506932481290842028650294*b4 + (-9.2668879968525419098066420771717743*b2 + 8.0313029306055363218324231335488711*a2)*a2)*a2;
    }

    inline double leg_f37_dyy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      double b6 = b5*b;
      return (.28081478778341036090323157809611437*b5 + (.28081478778341036090323157809611437*b4 + (-1.6848887267004621654193894685766862*b3 + (-1.6848887267004621654193894685766862*b2 + (1.8533775993705083819613284154343549*b + 1.8533775993705083819613284154343549*a)*a)*a)*a)*a)*a;
    }

    inline double leg_f37_dxy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      double b6 = b5*b;
      return -.46802464630568393483871929682685729e-1*b6 + (.28081478778341036090323157809611437*b5 + (1.2636665450253466240645421014325147*b4 + (-1.6848887267004621654193894685766862*b3 + (-4.6334439984262709549033210385858872*b2 + (1.8533775993705083819613284154343549*b + 4.0156514653027681609162115667744355*a)*a)*a)*a)*a)*a;
    }

    // number 38
    inline double leg_f38(double x, double y)
    {
      double y2 = y*y;
      return -.21395412402545551306912882140656333e-1*(x + 1.)*(-5. + (135. + (-495. + 429.*y2)*y2)*y2)*(1. + y);
    }

    inline double leg_f38_dx(double x, double y)
    {
      double y2 = y*y;
      return -.21395412402545551306912882140656333e-1*(-5. + (135. + (-495. + 429.*y2)*y2)*y2)*(1. + y);
    }

    inline double leg_f38_dy(double x, double y)
    {
      double y2 = y*y;
      return -.21395412402545551306912882140656333e-1*(x + 1.)*(-5. + (270. + (405. + (-1980. + (-2475. + (2574. + 3003.*y)*y)*y)*y)*y)*y);
    }

    inline double leg_f38_dxx(double x, double y)
    {
      return 0.;
    }

    inline double leg_f38_dyy(double x, double y)
    {
      return -.38511742324581992352443187853181400*(x + 1.)*(15. + (45. + (-330. + (-550. + (715. + 1001.*y)*y)*y)*y)*y);
    }

    inline double leg_f38_dxy(double x, double y)
    {
      return .10697706201272775653456441070328166-.21395412402545551306912882140656333e-1*(270. + (405. + (-1980. + (-2475. + (2574. + 3003.*y)*y)*y)*y)*y)*y;
    }

    // number 39
    inline double leg_f39(double x, double y)
    {
      double y2 = y*y;
      return .21395412402545551306912882140656333e-1*(x + y)*(-5. + (135. + (-495. + 429.*y2)*y2)*y2)*(1. + y);
    }

    inline double leg_f39_dx(double x, double y)
    {
      double y2 = y*y;
      return .21395412402545551306912882140656333e-1*(-5. + (135. + (-495. + 429.*y2)*y2)*y2)*(1. + y);
    }

    inline double leg_f39_dy(double x, double y)
    {
      double y2 = y*y;
      return -.10697706201272775653456441070328166 + .21395412402545551306912882140656333e-1*(-10. + (405. + (540. + (-2475. + (-2970. + (3003. + 3432.*y)*y)*y)*y)*y)*y)*y + .21395412402545551306912882140656333e-1*(-5. + (270. + (405. + (-1980. + (-2475. + (2574. + 3003.*y)*y)*y)*y)*y)*y)*x;
    }

    inline double leg_f39_dxx(double x, double y)
    {
      return 0.;
    }

    inline double leg_f39_dyy(double x, double y)
    {
      return -.21395412402545551306912882140656333 + .42790824805091102613825764281312666e-1*(405. + (810. + (-4950. + (-7425. + (9009. + 12012.*y)*y)*y)*y)*y)*y + .42790824805091102613825764281312666e-1*(135. + (405. + (-2970. + (-4950. + (6435. + 9009.*y)*y)*y)*y)*y)*x;
    }

    inline double leg_f39_dxy(double x, double y)
    {
      return -.10697706201272775653456441070328166 + .21395412402545551306912882140656333e-1*(270. + (405. + (-1980. + (-2475. + (2574. + 3003.*y)*y)*y)*y)*y)*y;
    }

    // Bubble functions, order 8

    // number 40
    inline double leg_f40(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      return 1.6233099319400270136250045367308206*(.61237243569579452454932101867647285*a2-.61237243569579452454932101867647285*b2)*(.78125000000000000000000000000000000e-1 + (.78125000000000000000000000000000000 + (-.70312500000000000000000000000000000 + (-6.5625000000000000000000000000000000 + (-2.5781250000000000000000000000000000 + (9.2812500000000000000000000000000000 + 6.7031250000000000000000000000000000*y)*y)*y)*y)*y)*y);
    }

    inline double leg_f40_dx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return 3.9762810276451551143588958519087301*a*(.78125000000000000000000000000000000e-1 + (.78125000000000000000000000000000000 + (-.70312500000000000000000000000000000 + (-6.5625000000000000000000000000000000 + (-2.5781250000000000000000000000000000 + (9.2812500000000000000000000000000000 + 6.7031250000000000000000000000000000*y)*y)*y)*y)*y)*y);
    }

    inline double leg_f40_dy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return (.15532347764238887165464436921518476 + 1.9881405138225775571794479259543650*(.78125000000000000000000000000000000 + (-.70312500000000000000000000000000000 + (-6.5625000000000000000000000000000000 + (-2.5781250000000000000000000000000000 + (9.2812500000000000000000000000000000 + 6.7031250000000000000000000000000000*y)*y)*y)*y)*y)*y-.99407025691128877858972396297718251*(.78125000000000000000000000000000000 + (-1.4062500000000000000000000000000000 + (-19.687500000000000000000000000000000 + (-10.312500000000000000000000000000000 + (46.406250000000000000000000000000000 + 40.218750000000000000000000000000000*y)*y)*y)*y)*y)*b)*b + (.15532347764238887165464436921518476 + 1.9881405138225775571794479259543650*(.78125000000000000000000000000000000 + (-.70312500000000000000000000000000000 + (-6.5625000000000000000000000000000000 + (-2.5781250000000000000000000000000000 + (9.2812500000000000000000000000000000 + 6.7031250000000000000000000000000000*y)*y)*y)*y)*y)*y + .99407025691128877858972396297718251*(.78125000000000000000000000000000000 + (-1.4062500000000000000000000000000000 + (-19.687500000000000000000000000000000 + (-10.312500000000000000000000000000000 + (46.406250000000000000000000000000000 + 40.218750000000000000000000000000000*y)*y)*y)*y)*y)*a)*a;
    }

    inline double leg_f40_dxx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return .62129391056955548661857747686073905 + 7.9525620552903102287177917038174600*(.78125000000000000000000000000000000 + (-.70312500000000000000000000000000000 + (-6.5625000000000000000000000000000000 + (-2.5781250000000000000000000000000000 + (9.2812500000000000000000000000000000 + 6.7031250000000000000000000000000000*y)*y)*y)*y)*y)*y;
    }

    inline double leg_f40_dyy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return (3.1064695528477774330928873843036954 + 3.9762810276451551143588958519087301*(-1.4062500000000000000000000000000000 + (-19.687500000000000000000000000000000 + (-10.312500000000000000000000000000000 + (46.406250000000000000000000000000000 + 40.218750000000000000000000000000000*y)*y)*y)*y)*y-.99407025691128877858972396297718251*(-1.4062500000000000000000000000000000 + (-39.375000000000000000000000000000000 + (-30.937500000000000000000000000000000 + (185.62500000000000000000000000000000 + 201.09375000000000000000000000000000*y)*y)*y)*y)*b)*b + (3.1064695528477774330928873843036954 + 3.9762810276451551143588958519087301*(-1.4062500000000000000000000000000000 + (-19.687500000000000000000000000000000 + (-10.312500000000000000000000000000000 + (46.406250000000000000000000000000000 + 40.218750000000000000000000000000000*y)*y)*y)*y)*y + .99407025691128877858972396297718251*(-1.4062500000000000000000000000000000 + (-39.375000000000000000000000000000000 + (-30.937500000000000000000000000000000 + (185.62500000000000000000000000000000 + 201.09375000000000000000000000000000*y)*y)*y)*y)*a)*a;
    }

    inline double leg_f40_dxy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return .31064695528477774330928873843036954 + 3.9762810276451551143588958519087301*(.78125000000000000000000000000000000 + (-.70312500000000000000000000000000000 + (-6.5625000000000000000000000000000000 + (-2.5781250000000000000000000000000000 + (9.2812500000000000000000000000000000 + 6.7031250000000000000000000000000000*y)*y)*y)*y)*y)*y + 3.9762810276451551143588958519087301*(.78125000000000000000000000000000000 + (-1.4062500000000000000000000000000000 + (-19.687500000000000000000000000000000 + (-10.312500000000000000000000000000000 + (46.406250000000000000000000000000000 + 40.218750000000000000000000000000000*y)*y)*y)*y)*y)*a;
    }

    // number 41
    inline double leg_f41(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      return 1.4187993420613458492725943245537961*(-.79056941504209483299972338610817962*b2 + .79056941504209483299972338610817962*a2)*a*(-.62500000000000000000000000000000000e-1 + (-1.8125000000000000000000000000000000 + (-3.1250000000000000000000000000000000 + (6.8750000000000000000000000000000000 + (17.187500000000000000000000000000000 + 8.9375000000000000000000000000000000*y)*y)*y)*y)*y);
    }

    inline double leg_f41_dx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      return 1.4187993420613458492725943245537961*(4.7434164902525689979983403166490778*a2-1.5811388300841896659994467722163592*b2)*(-.62500000000000000000000000000000000e-1 + (-1.8125000000000000000000000000000000 + (-3.1250000000000000000000000000000000 + (6.8750000000000000000000000000000000 + (17.187500000000000000000000000000000 + 8.9375000000000000000000000000000000*y)*y)*y)*y)*y);
    }

    inline double leg_f41_dy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      return -1.1216593659155472035363028403674680*b2*(-.62500000000000000000000000000000000e-1 + (-1.8125000000000000000000000000000000 + (-3.1250000000000000000000000000000000 + (6.8750000000000000000000000000000000 + (17.187500000000000000000000000000000 + 8.9375000000000000000000000000000000*y)*y)*y)*y)*y) + ((-.14020742073944340044203785504593349 + 2.2433187318310944070726056807349359*(-1.8125000000000000000000000000000000 + (-3.1250000000000000000000000000000000 + (6.8750000000000000000000000000000000 + (17.187500000000000000000000000000000 + 8.9375000000000000000000000000000000*y)*y)*y)*y)*y-1.1216593659155472035363028403674680*(-1.8125000000000000000000000000000000 + (-6.2500000000000000000000000000000000 + (20.625000000000000000000000000000000 + (68.750000000000000000000000000000000 + 44.687500000000000000000000000000000*y)*y)*y)*y)*b)*b + (-.21031113110916510066305678256890025 + 3.3649780977466416106089085211024040*(-1.8125000000000000000000000000000000 + (-3.1250000000000000000000000000000000 + (6.8750000000000000000000000000000000 + (17.187500000000000000000000000000000 + 8.9375000000000000000000000000000000*y)*y)*y)*y)*y + 1.1216593659155472035363028403674680*(-1.8125000000000000000000000000000000 + (-6.2500000000000000000000000000000000 + (20.625000000000000000000000000000000 + (68.750000000000000000000000000000000 + 44.687500000000000000000000000000000*y)*y)*y)*y)*a)*a)*a;
    }

    inline double leg_f41_dxx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return 26.919824781973132884871268168819232*a*(-.62500000000000000000000000000000000e-1 + (-1.8125000000000000000000000000000000 + (-3.1250000000000000000000000000000000 + (6.8750000000000000000000000000000000 + (17.187500000000000000000000000000000 + 8.9375000000000000000000000000000000*y)*y)*y)*y)*y);
    }

    inline double leg_f41_dyy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return (-.28041484147888680088407571009186700 + 4.4866374636621888141452113614698719*(-1.8125000000000000000000000000000000 + (-3.1250000000000000000000000000000000 + (6.8750000000000000000000000000000000 + (17.187500000000000000000000000000000 + 8.9375000000000000000000000000000000*y)*y)*y)*y)*y-2.2433187318310944070726056807349359*(-1.8125000000000000000000000000000000 + (-6.2500000000000000000000000000000000 + (20.625000000000000000000000000000000 + (68.750000000000000000000000000000000 + 44.687500000000000000000000000000000*y)*y)*y)*y)*b)*b + (-.28041484147888680088407571009186700 + 4.4866374636621888141452113614698719*(-1.8125000000000000000000000000000000 + (-3.1250000000000000000000000000000000 + (6.8750000000000000000000000000000000 + (17.187500000000000000000000000000000 + 8.9375000000000000000000000000000000*y)*y)*y)*y)*y + (-8.1320304028877172256381955926641427 + 4.4866374636621888141452113614698719*(-6.2500000000000000000000000000000000 + (20.625000000000000000000000000000000 + (68.750000000000000000000000000000000 + 44.687500000000000000000000000000000*y)*y)*y)*y-1.1216593659155472035363028403674680*(-6.2500000000000000000000000000000000 + (41.250000000000000000000000000000000 + (206.25000000000000000000000000000000 + 178.75000000000000000000000000000000*y)*y)*y)*b)*b + (-12.198045604331575838457293388996214 + 6.7299561954932832212178170422048078*(-6.2500000000000000000000000000000000 + (20.625000000000000000000000000000000 + (68.750000000000000000000000000000000 + 44.687500000000000000000000000000000*y)*y)*y)*y + 1.1216593659155472035363028403674680*(-6.2500000000000000000000000000000000 + (41.250000000000000000000000000000000 + (206.25000000000000000000000000000000 + 178.75000000000000000000000000000000*y)*y)*y)*a)*a)*a;
    }

    inline double leg_f41_dxy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return (-.28041484147888680088407571009186700 + 4.4866374636621888141452113614698719*(-1.8125000000000000000000000000000000 + (-3.1250000000000000000000000000000000 + (6.8750000000000000000000000000000000 + (17.187500000000000000000000000000000 + 8.9375000000000000000000000000000000*y)*y)*y)*y)*y-2.2433187318310944070726056807349359*(-1.8125000000000000000000000000000000 + (-6.2500000000000000000000000000000000 + (20.625000000000000000000000000000000 + (68.750000000000000000000000000000000 + 44.687500000000000000000000000000000*y)*y)*y)*y)*b)*b + (-.84124452443666040265222713027560098 + 13.459912390986566442435634084409616*(-1.8125000000000000000000000000000000 + (-3.1250000000000000000000000000000000 + (6.8750000000000000000000000000000000 + (17.187500000000000000000000000000000 + 8.9375000000000000000000000000000000*y)*y)*y)*y)*y + 6.7299561954932832212178170422048078*(-1.8125000000000000000000000000000000 + (-6.2500000000000000000000000000000000 + (20.625000000000000000000000000000000 + (68.750000000000000000000000000000000 + 44.687500000000000000000000000000000*y)*y)*y)*y)*a)*a;
    }

    // number 42
    inline double leg_f42(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b4 = b2*b2;
      return .82946008478893335842946748358641783*(.23385358667337133659898429576978433*b4 + (-1.4031215200402280195939057746187060*b2 + 1.1692679333668566829949214788489217*a2)*a2)*(-.31250000000000000000000000000000000 + (1.7500000000000000000000000000000000 + (12.375000000000000000000000000000000 + (19.250000000000000000000000000000000 + 8.9375000000000000000000000000000000*y)*y)*y)*y);
    }

    inline double leg_f42_dx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      return .82946008478893335842946748358641783*(-5.6124860801609120783756230984748240*b2 + 9.3541434669348534639593718307913732*a2)*a*(-.31250000000000000000000000000000000 + (1.7500000000000000000000000000000000 + (12.375000000000000000000000000000000 + (19.250000000000000000000000000000000 + 8.9375000000000000000000000000000000*y)*y)*y)*y);
    }

    inline double leg_f42_dy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      return (.24246526978786345616698149680517773-.77588886332116305973434078977656874*(1.7500000000000000000000000000000000 + (12.375000000000000000000000000000000 + (19.250000000000000000000000000000000 + 8.9375000000000000000000000000000000*y)*y)*y)*y + .19397221583029076493358519744414218*(1.7500000000000000000000000000000000 + (24.750000000000000000000000000000000 + (57.750000000000000000000000000000000 + 35.750000000000000000000000000000000*y)*y)*y)*b)*b3 + (-2.3276665899634891792030223693297061*b2*(-.31250000000000000000000000000000000 + (1.7500000000000000000000000000000000 + (12.375000000000000000000000000000000 + (19.250000000000000000000000000000000 + 8.9375000000000000000000000000000000*y)*y)*y)*y) + ((-.72739580936359036850094449041553317 + 2.3276665899634891792030223693297061*(1.7500000000000000000000000000000000 + (12.375000000000000000000000000000000 + (19.250000000000000000000000000000000 + 8.9375000000000000000000000000000000*y)*y)*y)*y-1.1638332949817445896015111846648531*(1.7500000000000000000000000000000000 + (24.750000000000000000000000000000000 + (57.750000000000000000000000000000000 + 35.750000000000000000000000000000000*y)*y)*y)*b)*b + (-1.2123263489393172808349074840258886 + 3.8794443166058152986717039488828435*(1.7500000000000000000000000000000000 + (12.375000000000000000000000000000000 + (19.250000000000000000000000000000000 + 8.9375000000000000000000000000000000*y)*y)*y)*y + .96986107915145382466792598722071090*(1.7500000000000000000000000000000000 + (24.750000000000000000000000000000000 + (57.750000000000000000000000000000000 + 35.750000000000000000000000000000000*y)*y)*y)*a)*a)*a)*a;
    }

    inline double leg_f42_dxx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      return .82946008478893335842946748358641783*(56.124860801609120783756230984748240*a2-11.224972160321824156751246196949648*b2)*(-.31250000000000000000000000000000000 + (1.7500000000000000000000000000000000 + (12.375000000000000000000000000000000 + (19.250000000000000000000000000000000 + 8.9375000000000000000000000000000000*y)*y)*y)*y);
    }

    inline double leg_f42_dyy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      return (-2.7156110216240707090701927642179905-1.5517777266423261194686815795531374*(24.750000000000000000000000000000000 + (57.750000000000000000000000000000000 + 35.750000000000000000000000000000000*y)*y)*y + .19397221583029076493358519744414218*(24.750000000000000000000000000000000 + (115.50000000000000000000000000000000 + 107.25000000000000000000000000000000*y)*y)*b)*b3 + ((-2.9095832374543614740037779616621327 + 9.3106663598539567168120894773188247*(1.7500000000000000000000000000000000 + (12.375000000000000000000000000000000 + (19.250000000000000000000000000000000 + 8.9375000000000000000000000000000000*y)*y)*y)*y-4.6553331799269783584060447386594123*(1.7500000000000000000000000000000000 + (24.750000000000000000000000000000000 + (57.750000000000000000000000000000000 + 35.750000000000000000000000000000000*y)*y)*y)*b)*b + (-2.9095832374543614740037779616621327 + 9.3106663598539567168120894773188247*(1.7500000000000000000000000000000000 + (12.375000000000000000000000000000000 + (19.250000000000000000000000000000000 + 8.9375000000000000000000000000000000*y)*y)*y)*y + (8.1468330648722121272105782926539716 + 4.6553331799269783584060447386594123*(24.750000000000000000000000000000000 + (57.750000000000000000000000000000000 + 35.750000000000000000000000000000000*y)*y)*y-1.1638332949817445896015111846648531*(24.750000000000000000000000000000000 + (115.50000000000000000000000000000000 + 107.25000000000000000000000000000000*y)*y)*b)*b + (13.578055108120353545350963821089953 + 7.7588886332116305973434078977656874*(24.750000000000000000000000000000000 + (57.750000000000000000000000000000000 + 35.750000000000000000000000000000000*y)*y)*y + .96986107915145382466792598722071090*(24.750000000000000000000000000000000 + (115.50000000000000000000000000000000 + 107.25000000000000000000000000000000*y)*y)*a)*a)*a)*a;
    }

    inline double leg_f42_dxy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      return -4.6553331799269783584060447386594123*b2*(-.31250000000000000000000000000000000 + (1.7500000000000000000000000000000000 + (12.375000000000000000000000000000000 + (19.250000000000000000000000000000000 + 8.9375000000000000000000000000000000*y)*y)*y)*y) + ((-2.9095832374543614740037779616621327 + 9.3106663598539567168120894773188247*(1.7500000000000000000000000000000000 + (12.375000000000000000000000000000000 + (19.250000000000000000000000000000000 + 8.9375000000000000000000000000000000*y)*y)*y)*y-4.6553331799269783584060447386594123*(1.7500000000000000000000000000000000 + (24.750000000000000000000000000000000 + (57.750000000000000000000000000000000 + 35.750000000000000000000000000000000*y)*y)*y)*b)*b + (-7.2739580936359036850094449041553317 + 23.276665899634891792030223693297061*(1.7500000000000000000000000000000000 + (12.375000000000000000000000000000000 + (19.250000000000000000000000000000000 + 8.9375000000000000000000000000000000*y)*y)*y)*y + 7.7588886332116305973434078977656874*(1.7500000000000000000000000000000000 + (24.750000000000000000000000000000000 + (57.750000000000000000000000000000000 + 35.750000000000000000000000000000000*y)*y)*y)*a)*a)*a;
    }

    // number 43
    inline double leg_f43(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b4 = b2*b2;
      return .41851680570743669240797557420461342*(.79549512883486596495094990736795518*b4 + (-2.6516504294495532165031663578931839*b2 + 1.8561553006146872515522164505252288*a2)*a2)*a*(1.5000000000000000000000000000000000 + (8.5000000000000000000000000000000000 + (13.500000000000000000000000000000000 + 6.5000000000000000000000000000000000*y)*y)*y);
    }

    inline double leg_f43_dx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      return .41851680570743669240797557420461342*(1.5909902576697319299018998147359104*b4 + (-15.909902576697319299018998147359104*b2 + 18.561553006146872515522164505252288*a2)*a2)*(1.5000000000000000000000000000000000 + (8.5000000000000000000000000000000000 + (13.500000000000000000000000000000000 + 6.5000000000000000000000000000000000*y)*y)*y);
    }

    inline double leg_f43_dy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      return .33292808027579391902410751141834787*b4*(1.5000000000000000000000000000000000 + (8.5000000000000000000000000000000000 + (13.500000000000000000000000000000000 + 6.5000000000000000000000000000000000*y)*y)*y) + ((-1.9975684816547635141446450685100873-1.3317123211031756760964300456733915*(8.5000000000000000000000000000000000 + (13.500000000000000000000000000000000 + 6.5000000000000000000000000000000000*y)*y)*y + .33292808027579391902410751141834787*(8.5000000000000000000000000000000000 + (27. + 19.500000000000000000000000000000000*y)*y)*b)*b3 + (-3.3292808027579391902410751141834787*b2*(1.5000000000000000000000000000000000 + (8.5000000000000000000000000000000000 + (13.500000000000000000000000000000000 + 6.5000000000000000000000000000000000*y)*y)*y) + ((3.3292808027579391902410751141834788 + 2.2195205351719594601607167427889859*(8.5000000000000000000000000000000000 + (13.500000000000000000000000000000000 + 6.5000000000000000000000000000000000*y)*y)*y-1.1097602675859797300803583713944929*(8.5000000000000000000000000000000000 + (27. + 19.500000000000000000000000000000000*y)*y)*b)*b + (5.8262414048263935829218814498210879 + 3.8841609365509290552812542998807252*(8.5000000000000000000000000000000000 + (13.500000000000000000000000000000000 + 6.5000000000000000000000000000000000*y)*y)*y + .77683218731018581105625085997614505*(8.5000000000000000000000000000000000 + (27. + 19.500000000000000000000000000000000*y)*y)*a)*a)*a)*a)*a;
    }

    inline double leg_f43_dxx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      return .41851680570743669240797557420461342*(-63.639610306789277196075992589436414*b2 + 148.49242404917498012417731604201830*a2)*a*(1.5000000000000000000000000000000000 + (8.5000000000000000000000000000000000 + (13.500000000000000000000000000000000 + 6.5000000000000000000000000000000000*y)*y)*y);
    }

    inline double leg_f43_dyy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      return (-3.9951369633095270282892901370201747-2.6634246422063513521928600913467831*(8.5000000000000000000000000000000000 + (13.500000000000000000000000000000000 + 6.5000000000000000000000000000000000*y)*y)*y + .66585616055158783804821502283669576*(8.5000000000000000000000000000000000 + (27. + 19.500000000000000000000000000000000*y)*y)*b)*b3 + ((-3.9951369633095270282892901370201747-2.6634246422063513521928600913467831*(8.5000000000000000000000000000000000 + (13.500000000000000000000000000000000 + 6.5000000000000000000000000000000000*y)*y)*y + (-22.639109458753986493639310776447656-2.6634246422063513521928600913467831*(27. + 19.500000000000000000000000000000000*y)*y + .33292808027579391902410751141834787*(39.*y + 27.)*b)*b)*b2 + ((19.975684816547635141446450685100873 + 13.317123211031756760964300456733915*(8.5000000000000000000000000000000000 + (13.500000000000000000000000000000000 + 6.5000000000000000000000000000000000*y)*y)*y-6.6585616055158783804821502283669576*(8.5000000000000000000000000000000000 + (27. + 19.500000000000000000000000000000000*y)*y)*b)*b + (19.975684816547635141446450685100873 + 13.317123211031756760964300456733915*(8.5000000000000000000000000000000000 + (13.500000000000000000000000000000000 + 6.5000000000000000000000000000000000*y)*y)*y + (37.731849097923310822732184627412759 + 4.4390410703439189203214334855779718*(27. + 19.500000000000000000000000000000000*y)*y-1.1097602675859797300803583713944929*(39.*y + 27.)*b)*b + (66.030735921365793939781323097972327 + 7.7683218731018581105625085997614505*(27. + 19.500000000000000000000000000000000*y)*y + .77683218731018581105625085997614505*(39.*y + 27.)*a)*a)*a)*a)*a;
    }

    inline double leg_f43_dxy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      return (-3.9951369633095270282892901370201747-2.6634246422063513521928600913467831*(8.5000000000000000000000000000000000 + (13.500000000000000000000000000000000 + 6.5000000000000000000000000000000000*y)*y)*y + .66585616055158783804821502283669576*(8.5000000000000000000000000000000000 + (27. + 19.500000000000000000000000000000000*y)*y)*b)*b3 + (-13.317123211031756760964300456733915*b2*(1.5000000000000000000000000000000000 + (8.5000000000000000000000000000000000 + (13.500000000000000000000000000000000 + 6.5000000000000000000000000000000000*y)*y)*y) + ((19.975684816547635141446450685100873 + 13.317123211031756760964300456733915*(8.5000000000000000000000000000000000 + (13.500000000000000000000000000000000 + 6.5000000000000000000000000000000000*y)*y)*y-6.6585616055158783804821502283669576*(8.5000000000000000000000000000000000 + (27. + 19.500000000000000000000000000000000*y)*y)*b)*b + (46.609931238611148663375051598568703 + 31.073287492407432442250034399045802*(8.5000000000000000000000000000000000 + (13.500000000000000000000000000000000 + 6.5000000000000000000000000000000000*y)*y)*y + 7.7683218731018581105625085997614505*(8.5000000000000000000000000000000000 + (27. + 19.500000000000000000000000000000000*y)*y)*a)*a)*a)*a;
    }

    // number 44
    inline double leg_f44(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b4 = b2*b2;
      double b6 = b4*b2;
      return .19784376250916950649135326318598913*(-.14657549249448217358017594104826457*b6 + (2.1986323874172326037026391157239686*b4 + (-5.1301422373068760753061579366892600*b2 + 3.0780853423841256451836947620135560*a2)*a2)*a2)*(2.2500000000000000000000000000000000 + (5.5000000000000000000000000000000000 + 3.2500000000000000000000000000000000*y)*y);
    }

    inline double leg_f44_dx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      return .19784376250916950649135326318598913*(8.7945295496689304148105564628958743*b4 + (-41.041137898455008602449263493514080*b2 + 36.937024108609507742204337144162672*a2)*a2)*a*(2.2500000000000000000000000000000000 + (5.5000000000000000000000000000000000 + 3.2500000000000000000000000000000000*y)*y);
    }

    inline double leg_f44_dy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      return (.39148713351102899664245105172126804 + .17399428156045733184108935632056358*(5.5000000000000000000000000000000000 + 3.2500000000000000000000000000000000*y)*y-.28999046926742888640181559386760596e-1*(6.5000000000000000000000000000000000*y + 5.5000000000000000000000000000000000)*b)*b5 + (.86997140780228665920544678160281788*b4*(2.2500000000000000000000000000000000 + (5.5000000000000000000000000000000000 + 3.2500000000000000000000000000000000*y)*y) + ((-3.9148713351102899664245105172126804-1.7399428156045733184108935632056358*(5.5000000000000000000000000000000000 + 3.2500000000000000000000000000000000*y)*y + .43498570390114332960272339080140894*(6.5000000000000000000000000000000000*y + 5.5000000000000000000000000000000000)*b)*b3 + (-4.0598665697440044096254183141464834*b2*(2.2500000000000000000000000000000000 + (5.5000000000000000000000000000000000 + 3.2500000000000000000000000000000000*y)*y) + ((4.5673498909620049608285956034147938 + 2.0299332848720022048127091570732417*(5.5000000000000000000000000000000000 + 3.2500000000000000000000000000000000*y)*y-1.0149666424360011024063545785366209*(6.5000000000000000000000000000000000*y + 5.5000000000000000000000000000000000)*b)*b + (8.2212298037316089294914720861466287 + 3.6538799127696039686628764827318351*(5.5000000000000000000000000000000000 + 3.2500000000000000000000000000000000*y)*y + .60897998546160066144381274712197252*(6.5000000000000000000000000000000000*y + 5.5000000000000000000000000000000000)*a)*a)*a)*a)*a)*a;
    }

    inline double leg_f44_dxx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      return .19784376250916950649135326318598913*(17.589059099337860829621112925791749*b4 + (-246.24682739073005161469558096108448*b2 + 369.37024108609507742204337144162672*a2)*a2)*(2.2500000000000000000000000000000000 + (5.5000000000000000000000000000000000 + 3.2500000000000000000000000000000000*y)*y);
    }

    inline double leg_f44_dyy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      return (2.2619256602859453139341616321673265*y + 1.9139370971650306502519829195261993-.18849380502382877616118013601394388*b)*b5 + ((-15.659485340441159865698042068850722-6.9597712624182932736435742528225433*(5.5000000000000000000000000000000000 + 3.2500000000000000000000000000000000*y)*y + 1.7399428156045733184108935632056358*(6.5000000000000000000000000000000000*y + 5.5000000000000000000000000000000000)*b)*b3 + ((-15.659485340441159865698042068850722-6.9597712624182932736435742528225433*(5.5000000000000000000000000000000000 + 3.2500000000000000000000000000000000*y)*y + (-22.619256602859453139341616321673265*y-19.139370971650306502519829195261993 + 2.8274070753574316424177020402091581*b)*b)*b2 + ((36.538799127696039686628764827318350 + 16.239466278976017638501673256585934*(5.5000000000000000000000000000000000 + 3.2500000000000000000000000000000000*y)*y-8.1197331394880088192508366282929670*(6.5000000000000000000000000000000000*y + 5.5000000000000000000000000000000000)*b)*b + (36.538799127696039686628764827318350 + 16.239466278976017638501673256585934*(5.5000000000000000000000000000000000 + 3.2500000000000000000000000000000000*y)*y + (26.389132703336028662565219041952143*y + 22.329266133592024252939800727805659-6.5972831758340071656413047604880354*b)*b + (47.500438866004851592617394275513857*y + 40.192679040465643655291641310050186 + 3.9583699055004042993847828562928213*a)*a)*a)*a)*a)*a;
    }

    inline double leg_f44_dxy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      return 1.7399428156045733184108935632056358*b4*(2.2500000000000000000000000000000000 + (5.5000000000000000000000000000000000 + 3.2500000000000000000000000000000000*y)*y) + ((-15.659485340441159865698042068850722-6.9597712624182932736435742528225433*(5.5000000000000000000000000000000000 + 3.2500000000000000000000000000000000*y)*y + 1.7399428156045733184108935632056358*(6.5000000000000000000000000000000000*y + 5.5000000000000000000000000000000000)*b)*b3 + (-24.359199418464026457752509884878900*b2*(2.2500000000000000000000000000000000 + (5.5000000000000000000000000000000000 + 3.2500000000000000000000000000000000*y)*y) + ((36.538799127696039686628764827318350 + 16.239466278976017638501673256585934*(5.5000000000000000000000000000000000 + 3.2500000000000000000000000000000000*y)*y-8.1197331394880088192508366282929670*(6.5000000000000000000000000000000000*y + 5.5000000000000000000000000000000000)*b)*b + (82.212298037316089294914720861466287 + 36.538799127696039686628764827318351*(5.5000000000000000000000000000000000 + 3.2500000000000000000000000000000000*y)*y + 7.3077598255392079373257529654636703*(6.5000000000000000000000000000000000*y + 5.5000000000000000000000000000000000)*a)*a)*a)*a)*a;
    }

    // number 45
    inline double leg_f45(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b4 = b2*b2;
      double b6 = b4*b2;
      return .91069207987015857989057523397342891e-1*(-.79672179899887262969191001703480969*b6 + (5.5770525929921084078433701192436678*b4 + (-10.038694667385795134118066214638602*b2 + 5.2583638733925593559666061124297439*a2)*a2)*a2)*a*(1. + y);
    }

    inline double leg_f45_dx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      double b6 = b5*b;
      return .91069207987015857989057523397342891e-1*(-1.5934435979977452593838200340696194*b6 + (33.462315557952650447060220715462007*b4 + (-100.38694667385795134118066214638602*b2 + 73.617094227495830983532485574016415*a2)*a2)*a2)*(1. + y);
    }

    inline double leg_f45_dy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      double b6 = b5*b;
      return -.72556823220817774297504599550421409e-1*b6*(1. + y) + ((.43534093932490664578502759730252845 + .43534093932490664578502759730252845*y-.72556823220817774297504599550421409e-1*b)*b5 + (1.5236932876371732602475965905588496*b4*(1. + y) + ((-2.0315910501828976803301287874117995-2.0315910501828976803301287874117995*y + .50789776254572442008253219685294985*b)*b3 + (-4.5710798629115197807427897716765487*b2*(1. + y) + ((1.8284319451646079122971159086706195 + 1.8284319451646079122971159086706195*y-.91421597258230395614855795433530975*b)*b + (3.3521252328017811725447124992294690 + 3.3521252328017811725447124992294690*y + .47887503325739731036353035703278130*a)*a)*a)*a)*a)*a)*a;
    }

    inline double leg_f45_dxx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      return .91069207987015857989057523397342891e-1*(133.84926223181060178824088286184803*b4 + (-803.09557339086361072944529717108816*b2 + 883.40513072994997180238982688819698*a2)*a2)*a*(1. + y);
    }

    inline double leg_f45_dyy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      return (.87068187864981329157005519460505690 + .87068187864981329157005519460505690*y-.14511364644163554859500919910084282*b)*b5 + ((.87068187864981329157005519460505690 + .87068187864981329157005519460505690*y + .87068187864981329157005519460505690*b)*b4 + ((-12.189546301097386081980772724470796-12.189546301097386081980772724470796*y + 3.0473865752743465204951931811176992*b)*b3 + ((-12.189546301097386081980772724470796-12.189546301097386081980772724470796*y-4.0631821003657953606602575748235988*b)*b2 + ((18.284319451646079122971159086706195 + 18.284319451646079122971159086706195*y-9.1421597258230395614855795433530975*b)*b + (18.284319451646079122971159086706195 + 18.284319451646079122971159086706195*y + 3.6568638903292158245942318173412390*b + 6.7042504656035623450894249984589381*a)*a)*a)*a)*a)*a;
    }

    inline double leg_f45_dxy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      return (.87068187864981329157005519460505690 + .87068187864981329157005519460505690*y-.14511364644163554859500919910084282*b)*b5 + (6.0947731505486930409903863622353983*b4*(1. + y) + ((-12.189546301097386081980772724470796-12.189546301097386081980772724470796*y + 3.0473865752743465204951931811176992*b)*b3 + (-36.568638903292158245942318173412390*b2*(1. + y) + ((18.284319451646079122971159086706195 + 18.284319451646079122971159086706195*y-9.1421597258230395614855795433530975*b)*b + (40.225502793621374070536549990753629 + 40.225502793621374070536549990753629*y + 6.7042504656035623450894249984589381*a)*a)*a)*a)*a)*a;
    }

    // ORDER 9

    // Edge functions, order 9

    // number 46
    inline double leg_f46_0(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b4 = b2*b2;
      double b6 = b4*b2;
      double b8 = b6*b2;
      return (.15570321374480096167037609154822953e-2*b8 + (-.18684385649376115400445130985787544e-1*b6 + (.61658472642941180821468932253098894e-1*b4 + (-.76339061367450985778961535170503392e-1*b2 + .31807942236437910741233972987709747e-1*a2)*a2)*a2)*a2)*a;
    }

    inline double leg_f46_1(double x, double y)
    {
      return -leg_f46_0(x, y);
    }

    inline double leg_f46_dx_0(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      double b6 = b5*b;
      double b7 = b6*b;
      double b8 = b7*b;
      return .31140642748960192334075218309645906e-2*b8 + (-.11210631389625669240267078591472526*b6 + (.61658472642941180821468932253098894*b4 + (-1.0687468591443138009054614923870475*b2 + .57254296025588239334221151377877544*a2)*a2)*a2)*a2;
    }

    inline double leg_f46_dx_1(double x, double y)
    {
      return -leg_f46_dx_0(x, y);
    }

    inline double leg_f46_dy_0(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      double b6 = b5*b;
      double b7 = b6*b;
      double b8 = b7*b;
      return .15570321374480096167037609154822953e-2*b8 + (-.12456257099584076933630087323858362e-1*b7 + (-.56053156948128346201335392957362631e-1*b6 + (.11210631389625669240267078591472526*b5 + (.30829236321470590410734466126549447*b4 + (-.24663389057176472328587572901239557*b3 + (-.53437342957215690045273074619352374*b2 + (.15267812273490197155792307034100678*b + .28627148012794119667110575688938772*a)*a)*a)*a)*a)*a)*a)*a;
    }

    inline double leg_f46_dy_1(double x, double y)
    {
      return -leg_f46_dy_0(x, y);
    }

    inline double leg_f46_dxx_0(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      double b6 = b5*b;
      double b7 = b6*b;
      return (-.44842525558502676961068314365890104*b6 + (4.9326778114352944657175145802479115*b4 + (-12.824962309731765610865537908644570*b2 + 9.1606873640941182934753842204604070*a2)*a2)*a2)*a;
    }

    inline double leg_f46_dxx_1(double x, double y)
    {
      return -leg_f46_dxx_0(x, y);
    }

    inline double leg_f46_dyy_0(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      double b6 = b5*b;
      double b7 = b6*b;
      return -.24912514199168153867260174647716725e-1*b7 + (-.24912514199168153867260174647716725e-1*b6 + (.67263788337754015441602471548835157*b5 + (.67263788337754015441602471548835157*b4 + (-2.4663389057176472328587572901239557*b3 + (-2.4663389057176472328587572901239557*b2 + (2.1374937182886276018109229847740950*b + 2.1374937182886276018109229847740950*a)*a)*a)*a)*a)*a)*a;
    }

    inline double leg_f46_dyy_1(double x, double y)
    {
      return -leg_f46_dyy_0(x, y);
    }

    inline double leg_f46_dxy_0(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      double b6 = b5*b;
      double b7 = b6*b;
      return -.24912514199168153867260174647716725e-1*b7 + (-.22421262779251338480534157182945052*b6 + (.67263788337754015441602471548835157*b5 + (2.4663389057176472328587572901239557*b4 + (-2.4663389057176472328587572901239557*b3 + (-6.4124811548658828054327689543222849*b2 + (2.1374937182886276018109229847740950*b + 4.5803436820470591467376921102302035*a)*a)*a)*a)*a)*a)*a;
    }

    inline double leg_f46_dxy_1(double x, double y)
    {
      return -leg_f46_dxy_0(x, y);
    }

    // number 47
    inline double leg_f47_0(double x, double y)
    {
      double y2 = y*y;
      return -.22777155839239454964352159677912434e-1*(x + 1.)*(-35. + (385. + (-1001. + 715.*y2)*y2)*y2)*y*(1. + y);
    }

    // number 47
    inline double leg_f47_1(double x, double y)
    {
      return -leg_f47_0(x, y);
    }

    inline double leg_f47_dx_0(double x, double y)
    {
      double y2 = y*y;
      return -.22777155839239454964352159677912434e-1*(-35. + (385. + (-1001. + 715.*y2)*y2)*y2)*y*(1. + y);
    }

    inline double leg_f47_dx_1(double x, double y)
    {
      return -leg_f47_dx_0(x, y);
    }

    inline double leg_f47_dy_0(double x, double y)
    {
      double y2 = y*y;
      return -.22777155839239454964352159677912434e-1*(x + 1.)*(-35. + (-70. + (1155. + (1540. + (-5005. + (-6006. + (5005. + 5720.*y)*y)*y)*y)*y)*y)*y);
    }

    inline double leg_f47_dy_1(double x, double y)
    {
      return -leg_f47_dy_0(x, y);
    }

    inline double leg_f47_dxx_0(double x, double y)
    {
      return 0.;
    }

    inline double leg_f47_dxx_1(double x, double y)
    {
      return -leg_f47_dxx_0(x, y);
    }

    inline double leg_f47_dyy_0(double x, double y)
    {
      return -1.5944009087467618475046511774538704*(x + 1.)*(-1. + (33. + (66. + (-286. + (-429. + (429. + 572.*y)*y)*y)*y)*y)*y);
    }

    inline double leg_f47_dyy_1(double x, double y)
    {
      return -leg_f47_dyy_0(x, y);
    }

    inline double leg_f47_dxy_0(double x, double y)
    {
      return .79720045437338092375232558872693520-.22777155839239454964352159677912434e-1*(-70. + (1155. + (1540. + (-5005. + (-6006. + (5005. + 5720.*y)*y)*y)*y)*y)*y)*y;
    }

    inline double leg_f47_dxy_1(double x, double y)
    {
      return -leg_f47_dxy_0(x, y);
    }

    // number 48
    inline double leg_f48_0(double x, double y)
    {
      double y2 = y*y;
      return -.22777155839239454964352159677912434e-1*(x + y)*(-35. + (385. + (-1001. + 715.*y2)*y2)*y2)*y*(1. + y);
    }

    // number 48
    inline double leg_f48_1(double x, double y)
    {
      return -leg_f48_0(x, y);
    }

    inline double leg_f48_dx_0(double x, double y)
    {
      double y2 = y*y;
      return -.22777155839239454964352159677912434e-1*(-35. + (385. + (-1001. + 715.*y2)*y2)*y2)*y*(1. + y);
    }

    inline double leg_f48_dx_1(double x, double y)
    {
      return -leg_f48_dx_0(x, y);
    }

    inline double leg_f48_dy_0(double x, double y)
    {
      double y2 = y*y;
      return -.22777155839239454964352159677912434e-1*(-70. + (-105. + (1540. + (1925. + (-6006. + (-7007. + (5720. + 6435.*y)*y)*y)*y)*y)*y)*y)*y-.22777155839239454964352159677912434e-1*(-35. + (-70. + (1155. + (1540. + (-5005. + (-6006. + (5005. + 5720.*y)*y)*y)*y)*y)*y)*y)*x;
    }

    inline double leg_f48_dy_1(double x, double y)
    {
      return -leg_f48_dy_0(x, y);
    }

    inline double leg_f48_dxx_0(double x, double y)
    {
      return 0.;
    }

    inline double leg_f48_dxx_1(double x, double y)
    {
      return -leg_f48_dxx_0(x, y);
    }

    inline double leg_f48_dyy_0(double x, double y)
    {
      return 1.5944009087467618475046511774538704-.45554311678478909928704319355824868e-1*(-105. + (2310. + (3850. + (-15015. + (-21021. + (20020. + 25740.*y)*y)*y)*y)*y)*y)*y-.45554311678478909928704319355824868e-1*(-35. + (1155. + (2310. + (-10010. + (-15015. + (15015. + 20020.*y)*y)*y)*y)*y)*y)*x;
    }

    inline double leg_f48_dyy_1(double x, double y)
    {
      return -leg_f48_dyy_0(x, y);
    }

    inline double leg_f48_dxy_0(double x, double y)
    {
      return .79720045437338092375232558872693520-.22777155839239454964352159677912434e-1*(-70. + (1155. + (1540. + (-5005. + (-6006. + (5005. + 5720.*y)*y)*y)*y)*y)*y)*y;
    }

    inline double leg_f48_dxy_1(double x, double y)
    {
      return -leg_f48_dxy_0(x, y);
    }

    // Bubble functions, order 9

    // number 49
    inline double leg_f49(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      return 1.7563249754401486360699529803715482*(.61237243569579452454932101867647285*a2-.61237243569579452454932101867647285*b2)*(-.78125000000000000000000000000000000e-1 + (.39062500000000000000000000000000000 + (3.0468750000000000000000000000000000 + (-.85937500000000000000000000000000000 + (-14.609375000000000000000000000000000 + (-6.7031250000000000000000000000000000 + (15.640625000000000000000000000000000 + 11.171875000000000000000000000000000*y)*y)*y)*y)*y)*y)*y);
    }

    inline double leg_f49_dx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return 4.3021000123345612730902482891395519*a*(-.78125000000000000000000000000000000e-1 + (.39062500000000000000000000000000000 + (3.0468750000000000000000000000000000 + (-.85937500000000000000000000000000000 + (-14.609375000000000000000000000000000 + (-6.7031250000000000000000000000000000 + (15.640625000000000000000000000000000 + 11.171875000000000000000000000000000*y)*y)*y)*y)*y)*y)*y);
    }

    inline double leg_f49_dy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return (-.16805078173181879973008782379451375 + 2.1510500061672806365451241445697760*(.39062500000000000000000000000000000 + (3.0468750000000000000000000000000000 + (-.85937500000000000000000000000000000 + (-14.609375000000000000000000000000000 + (-6.7031250000000000000000000000000000 + (15.640625000000000000000000000000000 + 11.171875000000000000000000000000000*y)*y)*y)*y)*y)*y)*y-1.0755250030836403182725620722848880*(.39062500000000000000000000000000000 + (6.0937500000000000000000000000000000 + (-2.5781250000000000000000000000000000 + (-58.437500000000000000000000000000000 + (-33.515625000000000000000000000000000 + (93.843750000000000000000000000000000 + 78.203125000000000000000000000000000*y)*y)*y)*y)*y)*y)*b)*b + (-.16805078173181879973008782379451375 + 2.1510500061672806365451241445697760*(.39062500000000000000000000000000000 + (3.0468750000000000000000000000000000 + (-.85937500000000000000000000000000000 + (-14.609375000000000000000000000000000 + (-6.7031250000000000000000000000000000 + (15.640625000000000000000000000000000 + 11.171875000000000000000000000000000*y)*y)*y)*y)*y)*y)*y + 1.0755250030836403182725620722848880*(.39062500000000000000000000000000000 + (6.0937500000000000000000000000000000 + (-2.5781250000000000000000000000000000 + (-58.437500000000000000000000000000000 + (-33.515625000000000000000000000000000 + (93.843750000000000000000000000000000 + 78.203125000000000000000000000000000*y)*y)*y)*y)*y)*y)*a)*a;
    }

    inline double leg_f49_dxx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return -.67220312692727519892035129517805500 + 8.6042000246691225461804965782791037*(.39062500000000000000000000000000000 + (3.0468750000000000000000000000000000 + (-.85937500000000000000000000000000000 + (-14.609375000000000000000000000000000 + (-6.7031250000000000000000000000000000 + (15.640625000000000000000000000000000 + 11.171875000000000000000000000000000*y)*y)*y)*y)*y)*y)*y;
    }

    inline double leg_f49_dyy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return (1.6805078173181879973008782379451375 + 4.3021000123345612730902482891395519*(6.0937500000000000000000000000000000 + (-2.5781250000000000000000000000000000 + (-58.437500000000000000000000000000000 + (-33.515625000000000000000000000000000 + (93.843750000000000000000000000000000 + 78.203125000000000000000000000000000*y)*y)*y)*y)*y)*y-1.0755250030836403182725620722848880*(6.0937500000000000000000000000000000 + (-5.1562500000000000000000000000000000 + (-175.31250000000000000000000000000000 + (-134.06250000000000000000000000000000 + (469.21875000000000000000000000000000 + 469.21875000000000000000000000000000*y)*y)*y)*y)*y)*b)*b + (1.6805078173181879973008782379451375 + 4.3021000123345612730902482891395519*(6.0937500000000000000000000000000000 + (-2.5781250000000000000000000000000000 + (-58.437500000000000000000000000000000 + (-33.515625000000000000000000000000000 + (93.843750000000000000000000000000000 + 78.203125000000000000000000000000000*y)*y)*y)*y)*y)*y + 1.0755250030836403182725620722848880*(6.0937500000000000000000000000000000 + (-5.1562500000000000000000000000000000 + (-175.31250000000000000000000000000000 + (-134.06250000000000000000000000000000 + (469.21875000000000000000000000000000 + 469.21875000000000000000000000000000*y)*y)*y)*y)*y)*a)*a;
    }

    inline double leg_f49_dxy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return -.33610156346363759946017564758902751 + 4.3021000123345612730902482891395519*(.39062500000000000000000000000000000 + (3.0468750000000000000000000000000000 + (-.85937500000000000000000000000000000 + (-14.609375000000000000000000000000000 + (-6.7031250000000000000000000000000000 + (15.640625000000000000000000000000000 + 11.171875000000000000000000000000000*y)*y)*y)*y)*y)*y)*y + 4.3021000123345612730902482891395519*(.39062500000000000000000000000000000 + (6.0937500000000000000000000000000000 + (-2.5781250000000000000000000000000000 + (-58.437500000000000000000000000000000 + (-33.515625000000000000000000000000000 + (93.843750000000000000000000000000000 + 78.203125000000000000000000000000000*y)*y)*y)*y)*y)*y)*a;
    }

    // number 50
    inline double leg_f50(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      return 1.5921853230058714151559241260993711*(-.79056941504209483299972338610817962*b2 + .79056941504209483299972338610817962*a2)*a*(.20312500000000000000000000000000000 + (.31250000000000000000000000000000000e-1 + (-6.0156250000000000000000000000000000 + (-10.312500000000000000000000000000000 + (11.171875000000000000000000000000000 + (31.281250000000000000000000000000000 + 15.640625000000000000000000000000000*y)*y)*y)*y)*y)*y);
    }

    inline double leg_f50_dx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      return 1.5921853230058714151559241260993711*(4.7434164902525689979983403166490778*a2-1.5811388300841896659994467722163592*b2)*(.20312500000000000000000000000000000 + (.31250000000000000000000000000000000e-1 + (-6.0156250000000000000000000000000000 + (-10.312500000000000000000000000000000 + (11.171875000000000000000000000000000 + (31.281250000000000000000000000000000 + 15.640625000000000000000000000000000*y)*y)*y)*y)*y)*y);
    }

    inline double leg_f50_dy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      return -1.2587330194473605815215838648165324*b2*(.20312500000000000000000000000000000 + (.31250000000000000000000000000000000e-1 + (-6.0156250000000000000000000000000000 + (-10.312500000000000000000000000000000 + (11.171875000000000000000000000000000 + (31.281250000000000000000000000000000 + 15.640625000000000000000000000000000*y)*y)*y)*y)*y)*y) + ((.51136028915049023624314344508171630 + 2.5174660388947211630431677296330648*(.31250000000000000000000000000000000e-1 + (-6.0156250000000000000000000000000000 + (-10.312500000000000000000000000000000 + (11.171875000000000000000000000000000 + (31.281250000000000000000000000000000 + 15.640625000000000000000000000000000*y)*y)*y)*y)*y)*y-1.2587330194473605815215838648165324*(.31250000000000000000000000000000000e-1 + (-12.031250000000000000000000000000000 + (-30.937500000000000000000000000000000 + (44.687500000000000000000000000000000 + (156.40625000000000000000000000000000 + 93.843750000000000000000000000000000*y)*y)*y)*y)*y)*b)*b + (.76704043372573535436471516762257442 + 3.7761990583420817445647515944495971*(.31250000000000000000000000000000000e-1 + (-6.0156250000000000000000000000000000 + (-10.312500000000000000000000000000000 + (11.171875000000000000000000000000000 + (31.281250000000000000000000000000000 + 15.640625000000000000000000000000000*y)*y)*y)*y)*y)*y + 1.2587330194473605815215838648165324*(.31250000000000000000000000000000000e-1 + (-12.031250000000000000000000000000000 + (-30.937500000000000000000000000000000 + (44.687500000000000000000000000000000 + (156.40625000000000000000000000000000 + 93.843750000000000000000000000000000*y)*y)*y)*y)*y)*a)*a)*a;
    }

    inline double leg_f50_dxx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return 30.209592466736653956518012755596777*a*(.20312500000000000000000000000000000 + (.31250000000000000000000000000000000e-1 + (-6.0156250000000000000000000000000000 + (-10.312500000000000000000000000000000 + (11.171875000000000000000000000000000 + (31.281250000000000000000000000000000 + 15.640625000000000000000000000000000*y)*y)*y)*y)*y)*y);
    }

    inline double leg_f50_dyy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return (1.0227205783009804724862868901634326 + 5.0349320777894423260863354592661295*(.31250000000000000000000000000000000e-1 + (-6.0156250000000000000000000000000000 + (-10.312500000000000000000000000000000 + (11.171875000000000000000000000000000 + (31.281250000000000000000000000000000 + 15.640625000000000000000000000000000*y)*y)*y)*y)*y)*y-2.5174660388947211630431677296330648*(.31250000000000000000000000000000000e-1 + (-12.031250000000000000000000000000000 + (-30.937500000000000000000000000000000 + (44.687500000000000000000000000000000 + (156.40625000000000000000000000000000 + 93.843750000000000000000000000000000*y)*y)*y)*y)*y)*b)*b + (1.0227205783009804724862868901634326 + 5.0349320777894423260863354592661295*(.31250000000000000000000000000000000e-1 + (-6.0156250000000000000000000000000000 + (-10.312500000000000000000000000000000 + (11.171875000000000000000000000000000 + (31.281250000000000000000000000000000 + 15.640625000000000000000000000000000*y)*y)*y)*y)*y)*y + (.15734162743092007269019798310206655 + 5.0349320777894423260863354592661295*(-12.031250000000000000000000000000000 + (-30.937500000000000000000000000000000 + (44.687500000000000000000000000000000 + (156.40625000000000000000000000000000 + 93.843750000000000000000000000000000*y)*y)*y)*y)*y-1.2587330194473605815215838648165324*(-12.031250000000000000000000000000000 + (-61.875000000000000000000000000000000 + (134.06250000000000000000000000000000 + (625.62500000000000000000000000000000 + 469.21875000000000000000000000000000*y)*y)*y)*y)*b)*b + (.23601244114638010903529697465309982 + 7.5523981166841634891295031888991943*(-12.031250000000000000000000000000000 + (-30.937500000000000000000000000000000 + (44.687500000000000000000000000000000 + (156.40625000000000000000000000000000 + 93.843750000000000000000000000000000*y)*y)*y)*y)*y + 1.2587330194473605815215838648165324*(-12.031250000000000000000000000000000 + (-61.875000000000000000000000000000000 + (134.06250000000000000000000000000000 + (625.62500000000000000000000000000000 + 469.21875000000000000000000000000000*y)*y)*y)*y)*a)*a)*a;
    }

    inline double leg_f50_dxy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return (1.0227205783009804724862868901634326 + 5.0349320777894423260863354592661295*(.31250000000000000000000000000000000e-1 + (-6.0156250000000000000000000000000000 + (-10.312500000000000000000000000000000 + (11.171875000000000000000000000000000 + (31.281250000000000000000000000000000 + 15.640625000000000000000000000000000*y)*y)*y)*y)*y)*y-2.5174660388947211630431677296330648*(.31250000000000000000000000000000000e-1 + (-12.031250000000000000000000000000000 + (-30.937500000000000000000000000000000 + (44.687500000000000000000000000000000 + (156.40625000000000000000000000000000 + 93.843750000000000000000000000000000*y)*y)*y)*y)*y)*b)*b + (3.0681617349029414174588606704902977 + 15.104796233368326978259006377798388*(.31250000000000000000000000000000000e-1 + (-6.0156250000000000000000000000000000 + (-10.312500000000000000000000000000000 + (11.171875000000000000000000000000000 + (31.281250000000000000000000000000000 + 15.640625000000000000000000000000000*y)*y)*y)*y)*y)*y + 7.5523981166841634891295031888991943*(.31250000000000000000000000000000000e-1 + (-12.031250000000000000000000000000000 + (-30.937500000000000000000000000000000 + (44.687500000000000000000000000000000 + (156.40625000000000000000000000000000 + 93.843750000000000000000000000000000*y)*y)*y)*y)*y)*a)*a;
    }

    // number 51
    inline double leg_f51(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b4 = b2*b2;
      return .97509833420206208638416899890309576*(.23385358667337133659898429576978433*b4 + (-1.4031215200402280195939057746187060*b2 + 1.1692679333668566829949214788489217*a2)*a2)*(-.43750000000000000000000000000000000 + (-2.6875000000000000000000000000000000 + (2.6250000000000000000000000000000000 + (27.625000000000000000000000000000000 + (39.812500000000000000000000000000000 + 17.062500000000000000000000000000000*y)*y)*y)*y)*y);
    }

    inline double leg_f51_dx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      return .97509833420206208638416899890309576*(-5.6124860801609120783756230984748240*b2 + 9.3541434669348534639593718307913732*a2)*a*(-.43750000000000000000000000000000000 + (-2.6875000000000000000000000000000000 + (2.6250000000000000000000000000000000 + (27.625000000000000000000000000000000 + (39.812500000000000000000000000000000 + 17.062500000000000000000000000000000*y)*y)*y)*y)*y);
    }

    inline double leg_f51_dy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      return (.39905292492166838819716431718425096-.91212097124952774445066129642114507*(-2.6875000000000000000000000000000000 + (2.6250000000000000000000000000000000 + (27.625000000000000000000000000000000 + (39.812500000000000000000000000000000 + 17.062500000000000000000000000000000*y)*y)*y)*y)*y + .22803024281238193611266532410528627*(-2.6875000000000000000000000000000000 + (5.2500000000000000000000000000000000 + (82.875000000000000000000000000000000 + (159.25000000000000000000000000000000 + 85.312500000000000000000000000000000*y)*y)*y)*y)*b)*b3 + (-2.7363629137485832333519838892634352*b2*(-.43750000000000000000000000000000000 + (-2.6875000000000000000000000000000000 + (2.6250000000000000000000000000000000 + (27.625000000000000000000000000000000 + (39.812500000000000000000000000000000 + 17.062500000000000000000000000000000*y)*y)*y)*y)*y) + ((-1.1971587747650051645914929515527529 + 2.7363629137485832333519838892634352*(-2.6875000000000000000000000000000000 + (2.6250000000000000000000000000000000 + (27.625000000000000000000000000000000 + (39.812500000000000000000000000000000 + 17.062500000000000000000000000000000*y)*y)*y)*y)*y-1.3681814568742916166759919446317176*(-2.6875000000000000000000000000000000 + (5.2500000000000000000000000000000000 + (82.875000000000000000000000000000000 + (159.25000000000000000000000000000000 + 85.312500000000000000000000000000000*y)*y)*y)*y)*b)*b + (-1.9952646246083419409858215859212548 + 4.5606048562476387222533064821057254*(-2.6875000000000000000000000000000000 + (2.6250000000000000000000000000000000 + (27.625000000000000000000000000000000 + (39.812500000000000000000000000000000 + 17.062500000000000000000000000000000*y)*y)*y)*y)*y + 1.1401512140619096805633266205264313*(-2.6875000000000000000000000000000000 + (5.2500000000000000000000000000000000 + (82.875000000000000000000000000000000 + (159.25000000000000000000000000000000 + 85.312500000000000000000000000000000*y)*y)*y)*y)*a)*a)*a)*a;
    }

    inline double leg_f51_dxx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      return .97509833420206208638416899890309576*(56.124860801609120783756230984748240*a2-11.224972160321824156751246196949648*b2)*(-.43750000000000000000000000000000000 + (-2.6875000000000000000000000000000000 + (2.6250000000000000000000000000000000 + (27.625000000000000000000000000000000 + (39.812500000000000000000000000000000 + 17.062500000000000000000000000000000*y)*y)*y)*y)*y);
    }

    inline double leg_f51_dyy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      return (4.9026502204662116264223044682636546-1.8242419424990554889013225928422901*(5.2500000000000000000000000000000000 + (82.875000000000000000000000000000000 + (159.25000000000000000000000000000000 + 85.312500000000000000000000000000000*y)*y)*y)*y + .22803024281238193611266532410528627*(5.2500000000000000000000000000000000 + (165.75000000000000000000000000000000 + (477.75000000000000000000000000000000 + 341.25000000000000000000000000000000*y)*y)*y)*b)*b3 + ((-4.7886350990600206583659718062110115 + 10.945451654994332933407935557053741*(-2.6875000000000000000000000000000000 + (2.6250000000000000000000000000000000 + (27.625000000000000000000000000000000 + (39.812500000000000000000000000000000 + 17.062500000000000000000000000000000*y)*y)*y)*y)*y-5.4727258274971664667039677785268703*(-2.6875000000000000000000000000000000 + (5.2500000000000000000000000000000000 + (82.875000000000000000000000000000000 + (159.25000000000000000000000000000000 + 85.312500000000000000000000000000000*y)*y)*y)*y)*b)*b + (-4.7886350990600206583659718062110115 + 10.945451654994332933407935557053741*(-2.6875000000000000000000000000000000 + (2.6250000000000000000000000000000000 + (27.625000000000000000000000000000000 + (39.812500000000000000000000000000000 + 17.062500000000000000000000000000000*y)*y)*y)*y)*y + (-14.707950661398634879266913404790964 + 5.4727258274971664667039677785268703*(5.2500000000000000000000000000000000 + (82.875000000000000000000000000000000 + (159.25000000000000000000000000000000 + 85.312500000000000000000000000000000*y)*y)*y)*y-1.3681814568742916166759919446317176*(5.2500000000000000000000000000000000 + (165.75000000000000000000000000000000 + (477.75000000000000000000000000000000 + 341.25000000000000000000000000000000*y)*y)*y)*b)*b + (-24.513251102331058132111522341318273 + 9.1212097124952774445066129642114507*(5.2500000000000000000000000000000000 + (82.875000000000000000000000000000000 + (159.25000000000000000000000000000000 + 85.312500000000000000000000000000000*y)*y)*y)*y + 1.1401512140619096805633266205264313*(5.2500000000000000000000000000000000 + (165.75000000000000000000000000000000 + (477.75000000000000000000000000000000 + 341.25000000000000000000000000000000*y)*y)*y)*a)*a)*a)*a;
    }

    inline double leg_f51_dxy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      return -5.4727258274971664667039677785268703*b2*(-.43750000000000000000000000000000000 + (-2.6875000000000000000000000000000000 + (2.6250000000000000000000000000000000 + (27.625000000000000000000000000000000 + (39.812500000000000000000000000000000 + 17.062500000000000000000000000000000*y)*y)*y)*y)*y) + ((-4.7886350990600206583659718062110115 + 10.945451654994332933407935557053741*(-2.6875000000000000000000000000000000 + (2.6250000000000000000000000000000000 + (27.625000000000000000000000000000000 + (39.812500000000000000000000000000000 + 17.062500000000000000000000000000000*y)*y)*y)*y)*y-5.4727258274971664667039677785268703*(-2.6875000000000000000000000000000000 + (5.2500000000000000000000000000000000 + (82.875000000000000000000000000000000 + (159.25000000000000000000000000000000 + 85.312500000000000000000000000000000*y)*y)*y)*y)*b)*b + (-11.971587747650051645914929515527529 + 27.363629137485832333519838892634352*(-2.6875000000000000000000000000000000 + (2.6250000000000000000000000000000000 + (27.625000000000000000000000000000000 + (39.812500000000000000000000000000000 + 17.062500000000000000000000000000000*y)*y)*y)*y)*y + 9.1212097124952774445066129642114507*(-2.6875000000000000000000000000000000 + (5.2500000000000000000000000000000000 + (82.875000000000000000000000000000000 + (159.25000000000000000000000000000000 + 85.312500000000000000000000000000000*y)*y)*y)*y)*a)*a)*a;
    }

    // number 52
    inline double leg_f52(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b4 = b2*b2;
      return .51404418387524420707255933386925244*(.79549512883486596495094990736795518*b4 + (-2.6516504294495532165031663578931839*b2 + 1.8561553006146872515522164505252288*a2)*a2)*a*(.21875000000000000000000000000000000 + (7.1250000000000000000000000000000000 + (26.812500000000000000000000000000000 + (34.125000000000000000000000000000000 + 14.218750000000000000000000000000000*y)*y)*y)*y);
    }

    inline double leg_f52_dx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      return .51404418387524420707255933386925244*(1.5909902576697319299018998147359104*b4 + (-15.909902576697319299018998147359104*b2 + 18.561553006146872515522164505252288*a2)*a2)*(.21875000000000000000000000000000000 + (7.1250000000000000000000000000000000 + (26.812500000000000000000000000000000 + (34.125000000000000000000000000000000 + 14.218750000000000000000000000000000*y)*y)*y)*y);
    }

    inline double leg_f52_dy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      return .40891964427865092013487337229616315*b4*(.21875000000000000000000000000000000 + (7.1250000000000000000000000000000000 + (26.812500000000000000000000000000000 + (34.125000000000000000000000000000000 + 14.218750000000000000000000000000000*y)*y)*y)*y) + ((-.35780468874381955511801420075914278-1.6356785771146036805394934891846526*(7.1250000000000000000000000000000000 + (26.812500000000000000000000000000000 + (34.125000000000000000000000000000000 + 14.218750000000000000000000000000000*y)*y)*y)*y + .40891964427865092013487337229616315*(7.1250000000000000000000000000000000 + (53.625000000000000000000000000000000 + (102.37500000000000000000000000000000 + 56.875000000000000000000000000000000*y)*y)*y)*b)*b3 + (-4.0891964427865092013487337229616315*b2*(.21875000000000000000000000000000000 + (7.1250000000000000000000000000000000 + (26.812500000000000000000000000000000 + (34.125000000000000000000000000000000 + 14.218750000000000000000000000000000*y)*y)*y)*y) + ((.59634114790636592519669033459857128 + 2.7261309618576728008991558153077544*(7.1250000000000000000000000000000000 + (26.812500000000000000000000000000000 + (34.125000000000000000000000000000000 + 14.218750000000000000000000000000000*y)*y)*y)*y-1.3630654809288364004495779076538772*(7.1250000000000000000000000000000000 + (53.625000000000000000000000000000000 + (102.37500000000000000000000000000000 + 56.875000000000000000000000000000000*y)*y)*y)*b)*b + (1.0435970088361403690942080855474997 + 4.7707291832509274015735226767885701*(7.1250000000000000000000000000000000 + (26.812500000000000000000000000000000 + (34.125000000000000000000000000000000 + 14.218750000000000000000000000000000*y)*y)*y)*y + .95414583665018548031470453535771402*(7.1250000000000000000000000000000000 + (53.625000000000000000000000000000000 + (102.37500000000000000000000000000000 + 56.875000000000000000000000000000000*y)*y)*y)*a)*a)*a)*a)*a;
    }

    inline double leg_f52_dxx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      return .51404418387524420707255933386925244*(-63.639610306789277196075992589436414*b2 + 148.49242404917498012417731604201830*a2)*a*(.21875000000000000000000000000000000 + (7.1250000000000000000000000000000000 + (26.812500000000000000000000000000000 + (34.125000000000000000000000000000000 + 14.218750000000000000000000000000000*y)*y)*y)*y);
    }

    inline double leg_f52_dyy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      return (-.71560937748763911023602840151828554-3.2713571542292073610789869783693052*(7.1250000000000000000000000000000000 + (26.812500000000000000000000000000000 + (34.125000000000000000000000000000000 + 14.218750000000000000000000000000000*y)*y)*y)*y + .81783928855730184026974674459232632*(7.1250000000000000000000000000000000 + (53.625000000000000000000000000000000 + (102.37500000000000000000000000000000 + 56.875000000000000000000000000000000*y)*y)*y)*b)*b3 + ((-.71560937748763911023602840151828554-3.2713571542292073610789869783693052*(7.1250000000000000000000000000000000 + (26.812500000000000000000000000000000 + (34.125000000000000000000000000000000 + 14.218750000000000000000000000000000*y)*y)*y)*y + (-23.308419723883102447687782220881300-3.2713571542292073610789869783693052*(53.625000000000000000000000000000000 + (102.37500000000000000000000000000000 + 56.875000000000000000000000000000000*y)*y)*y + .40891964427865092013487337229616315*(53.625000000000000000000000000000000 + (204.75000000000000000000000000000000 + 170.62500000000000000000000000000000*y)*y)*b)*b)*b2 + ((3.5780468874381955511801420075914278 + 16.356785771146036805394934891846526*(7.1250000000000000000000000000000000 + (26.812500000000000000000000000000000 + (34.125000000000000000000000000000000 + 14.218750000000000000000000000000000*y)*y)*y)*y-8.1783928855730184026974674459232632*(7.1250000000000000000000000000000000 + (53.625000000000000000000000000000000 + (102.37500000000000000000000000000000 + 56.875000000000000000000000000000000*y)*y)*y)*b)*b + (3.5780468874381955511801420075914278 + 16.356785771146036805394934891846526*(7.1250000000000000000000000000000000 + (26.812500000000000000000000000000000 + (34.125000000000000000000000000000000 + 14.218750000000000000000000000000000*y)*y)*y)*y + (38.847366206471837412812970368135500 + 5.4522619237153456017983116306155087*(53.625000000000000000000000000000000 + (102.37500000000000000000000000000000 + 56.875000000000000000000000000000000*y)*y)*y-1.3630654809288364004495779076538772*(53.625000000000000000000000000000000 + (204.75000000000000000000000000000000 + 170.62500000000000000000000000000000*y)*y)*b)*b + (67.982890861325715472422698144237124 + 9.5414583665018548031470453535771402*(53.625000000000000000000000000000000 + (102.37500000000000000000000000000000 + 56.875000000000000000000000000000000*y)*y)*y + .95414583665018548031470453535771402*(53.625000000000000000000000000000000 + (204.75000000000000000000000000000000 + 170.62500000000000000000000000000000*y)*y)*a)*a)*a)*a)*a;
    }

    inline double leg_f52_dxy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      return (-.71560937748763911023602840151828554-3.2713571542292073610789869783693052*(7.1250000000000000000000000000000000 + (26.812500000000000000000000000000000 + (34.125000000000000000000000000000000 + 14.218750000000000000000000000000000*y)*y)*y)*y + .81783928855730184026974674459232632*(7.1250000000000000000000000000000000 + (53.625000000000000000000000000000000 + (102.37500000000000000000000000000000 + 56.875000000000000000000000000000000*y)*y)*y)*b)*b3 + (-16.356785771146036805394934891846526*b2*(.21875000000000000000000000000000000 + (7.1250000000000000000000000000000000 + (26.812500000000000000000000000000000 + (34.125000000000000000000000000000000 + 14.218750000000000000000000000000000*y)*y)*y)*y) + ((3.5780468874381955511801420075914278 + 16.356785771146036805394934891846526*(7.1250000000000000000000000000000000 + (26.812500000000000000000000000000000 + (34.125000000000000000000000000000000 + 14.218750000000000000000000000000000*y)*y)*y)*y-8.1783928855730184026974674459232632*(7.1250000000000000000000000000000000 + (53.625000000000000000000000000000000 + (102.37500000000000000000000000000000 + 56.875000000000000000000000000000000*y)*y)*y)*b)*b + (8.3487760706891229527536646843799979 + 38.165833466007419212588181414308561*(7.1250000000000000000000000000000000 + (26.812500000000000000000000000000000 + (34.125000000000000000000000000000000 + 14.218750000000000000000000000000000*y)*y)*y)*y + 9.5414583665018548031470453535771402*(7.1250000000000000000000000000000000 + (53.625000000000000000000000000000000 + (102.37500000000000000000000000000000 + 56.875000000000000000000000000000000*y)*y)*y)*a)*a)*a)*a;
    }

    // number 53
    inline double leg_f53(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b4 = b2*b2;
      double b6 = b4*b2;
      return .25148937907519400465337169927700054*(-.14657549249448217358017594104826457*b6 + (2.1986323874172326037026391157239686*b4 + (-5.1301422373068760753061579366892600*b2 + 3.0780853423841256451836947620135560*a2)*a2)*a2)*(2.7500000000000000000000000000000000 + (13.250000000000000000000000000000000 + (19.250000000000000000000000000000000 + 8.7500000000000000000000000000000000*y)*y)*y);
    }

    inline double leg_f53_dx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      return .25148937907519400465337169927700054*(8.7945295496689304148105564628958743*b4 + (-41.041137898455008602449263493514080*b2 + 36.937024108609507742204337144162672*a2)*a2)*a*(2.7500000000000000000000000000000000 + (13.250000000000000000000000000000000 + (19.250000000000000000000000000000000 + 8.7500000000000000000000000000000000*y)*y)*y);
    }

    inline double leg_f53_dy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      return (.60822596331878833691165842189837477 + .22117307757046848614969397159940900*(13.250000000000000000000000000000000 + (19.250000000000000000000000000000000 + 8.7500000000000000000000000000000000*y)*y)*y-.36862179595078081024948995266568168e-1*(13.250000000000000000000000000000000 + (38.500000000000000000000000000000000 + 26.250000000000000000000000000000000*y)*y)*b)*b5 + (1.1058653878523424307484698579970450*b4*(2.7500000000000000000000000000000000 + (13.250000000000000000000000000000000 + (19.250000000000000000000000000000000 + 8.7500000000000000000000000000000000*y)*y)*y) + ((-6.0822596331878833691165842189837477-2.2117307757046848614969397159940900*(13.250000000000000000000000000000000 + (19.250000000000000000000000000000000 + 8.7500000000000000000000000000000000*y)*y)*y + .55293269392617121537423492899852251*(13.250000000000000000000000000000000 + (38.500000000000000000000000000000000 + 26.250000000000000000000000000000000*y)*y)*b)*b3 + (-5.1607051433109313434928593373195434*b2*(2.7500000000000000000000000000000000 + (13.250000000000000000000000000000000 + (19.250000000000000000000000000000000 + 8.7500000000000000000000000000000000*y)*y)*y) + ((7.0959695720525305973026815888143723 + 2.5803525716554656717464296686597717*(13.250000000000000000000000000000000 + (19.250000000000000000000000000000000 + 8.7500000000000000000000000000000000*y)*y)*y-1.2901762858277328358732148343298859*(13.250000000000000000000000000000000 + (38.500000000000000000000000000000000 + 26.250000000000000000000000000000000*y)*y)*b)*b + (12.772745229694555075144826859865870 + 4.6446346289798382091435734035875891*(13.250000000000000000000000000000000 + (19.250000000000000000000000000000000 + 8.7500000000000000000000000000000000*y)*y)*y + .77410577149663970152392890059793149*(13.250000000000000000000000000000000 + (38.500000000000000000000000000000000 + 26.250000000000000000000000000000000*y)*y)*a)*a)*a)*a)*a)*a;
    }

    inline double leg_f53_dxx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      return .25148937907519400465337169927700054*(17.589059099337860829621112925791749*b4 + (-246.24682739073005161469558096108448*b2 + 369.37024108609507742204337144162672*a2)*a2)*(2.7500000000000000000000000000000000 + (13.250000000000000000000000000000000 + (19.250000000000000000000000000000000 + 8.7500000000000000000000000000000000*y)*y)*y);
    }

    inline double leg_f53_dyy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      return (5.8610865556174148829668902473843387 + .44234615514093697229938794319881801*(38.500000000000000000000000000000000 + 26.250000000000000000000000000000000*y)*y-.36862179595078081024948995266568168e-1*(52.500000000000000000000000000000000*y + 38.500000000000000000000000000000000)*b)*b5 + ((-24.329038532751533476466336875934991-8.8469231028187394459877588639763602*(13.250000000000000000000000000000000 + (19.250000000000000000000000000000000 + 8.7500000000000000000000000000000000*y)*y)*y + 2.2117307757046848614969397159940900*(13.250000000000000000000000000000000 + (38.500000000000000000000000000000000 + 26.250000000000000000000000000000000*y)*y)*b)*b3 + ((-24.329038532751533476466336875934991-8.8469231028187394459877588639763602*(13.250000000000000000000000000000000 + (19.250000000000000000000000000000000 + 8.7500000000000000000000000000000000*y)*y)*y + (-58.610865556174148829668902473843387-4.4234615514093697229938794319881801*(38.500000000000000000000000000000000 + 26.250000000000000000000000000000000*y)*y + .55293269392617121537423492899852251*(52.500000000000000000000000000000000*y + 38.500000000000000000000000000000000)*b)*b)*b2 + ((56.767756576420244778421452710514978 + 20.642820573243725373971437349278174*(13.250000000000000000000000000000000 + (19.250000000000000000000000000000000 + 8.7500000000000000000000000000000000*y)*y)*y-10.321410286621862686985718674639087*(13.250000000000000000000000000000000 + (38.500000000000000000000000000000000 + 26.250000000000000000000000000000000*y)*y)*b)*b + (56.767756576420244778421452710514978 + 20.642820573243725373971437349278174*(13.250000000000000000000000000000000 + (19.250000000000000000000000000000000 + 8.7500000000000000000000000000000000*y)*y)*y + (68.379343148869840301280386219483951 + 5.1607051433109313434928593373195434*(38.500000000000000000000000000000000 + 26.250000000000000000000000000000000*y)*y-1.2901762858277328358732148343298859*(52.500000000000000000000000000000000*y + 38.500000000000000000000000000000000)*b)*b + (123.08281766796571254230469519507111 + 9.2892692579596764182871468071751781*(38.500000000000000000000000000000000 + 26.250000000000000000000000000000000*y)*y + .77410577149663970152392890059793149*(52.500000000000000000000000000000000*y + 38.500000000000000000000000000000000)*a)*a)*a)*a)*a)*a;
    }

    inline double leg_f53_dxy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      return 2.2117307757046848614969397159940900*b4*(2.7500000000000000000000000000000000 + (13.250000000000000000000000000000000 + (19.250000000000000000000000000000000 + 8.7500000000000000000000000000000000*y)*y)*y) + ((-24.329038532751533476466336875934991-8.8469231028187394459877588639763602*(13.250000000000000000000000000000000 + (19.250000000000000000000000000000000 + 8.7500000000000000000000000000000000*y)*y)*y + 2.2117307757046848614969397159940900*(13.250000000000000000000000000000000 + (38.500000000000000000000000000000000 + 26.250000000000000000000000000000000*y)*y)*b)*b3 + (-30.964230859865588060957156023917261*b2*(2.7500000000000000000000000000000000 + (13.250000000000000000000000000000000 + (19.250000000000000000000000000000000 + 8.7500000000000000000000000000000000*y)*y)*y) + ((56.767756576420244778421452710514978 + 20.642820573243725373971437349278174*(13.250000000000000000000000000000000 + (19.250000000000000000000000000000000 + 8.7500000000000000000000000000000000*y)*y)*y-10.321410286621862686985718674639087*(13.250000000000000000000000000000000 + (38.500000000000000000000000000000000 + 26.250000000000000000000000000000000*y)*y)*b)*b + (127.72745229694555075144826859865870 + 46.446346289798382091435734035875891*(13.250000000000000000000000000000000 + (19.250000000000000000000000000000000 + 8.7500000000000000000000000000000000*y)*y)*y + 9.2892692579596764182871468071751781*(13.250000000000000000000000000000000 + (38.500000000000000000000000000000000 + 26.250000000000000000000000000000000*y)*y)*a)*a)*a)*a)*a;
    }

    // number 54
    inline double leg_f54(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b4 = b2*b2;
      double b6 = b4*b2;
      return .11865250751057162342943393039327360*(-.79672179899887262969191001703480969*b6 + (5.5770525929921084078433701192436678*b4 + (-10.038694667385795134118066214638602*b2 + 5.2583638733925593559666061124297439*a2)*a2)*a2)*a*(2.7500000000000000000000000000000000 + (6.5000000000000000000000000000000000 + 3.7500000000000000000000000000000000*y)*y);
    }

    inline double leg_f54_dx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      double b6 = b5*b;
      return .11865250751057162342943393039327360*(-1.5934435979977452593838200340696194*b6 + (33.462315557952650447060220715462007*b4 + (-100.38694667385795134118066214638602*b2 + 73.617094227495830983532485574016415*a2)*a2)*a2)*(2.7500000000000000000000000000000000 + (6.5000000000000000000000000000000000 + 3.7500000000000000000000000000000000*y)*y);
    }

    inline double leg_f54_dy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      double b6 = b5*b;
      return -.94533039239549870023101308425754951e-1*b6*(2.7500000000000000000000000000000000 + (6.5000000000000000000000000000000000 + 3.7500000000000000000000000000000000*y)*y) + ((1.5597951474525728553811715890249567 + .56719823543729922013860785055452969*(6.5000000000000000000000000000000000 + 3.7500000000000000000000000000000000*y)*y-.94533039239549870023101308425754951e-1*(7.5000000000000000000000000000000000*y + 6.5000000000000000000000000000000000)*b)*b5 + (1.9851938240305472704851274769408539*b4*(2.7500000000000000000000000000000000 + (6.5000000000000000000000000000000000 + 3.7500000000000000000000000000000000*y)*y) + ((-7.2790440214453399917788007487831311-2.6469250987073963606468366359211386*(6.5000000000000000000000000000000000 + 3.7500000000000000000000000000000000*y)*y + .66173127467684909016170915898028464*(7.5000000000000000000000000000000000*y + 6.5000000000000000000000000000000000)*b)*b3 + (-5.9555814720916418114553824308225617*b2*(2.7500000000000000000000000000000000 + (6.5000000000000000000000000000000000 + 3.7500000000000000000000000000000000*y)*y) + ((6.5511396193008059926009206739048179 + 2.3822325888366567245821529723290247*(6.5000000000000000000000000000000000 + 3.7500000000000000000000000000000000*y)*y-1.1911162944183283622910764861645123*(7.5000000000000000000000000000000000*y + 6.5000000000000000000000000000000000)*b)*b + (12.010422635384810986435021235492166 + 4.3674264128672039950672804492698786*(6.5000000000000000000000000000000000 + 3.7500000000000000000000000000000000*y)*y + .62391805898102914215246863560998266*(7.5000000000000000000000000000000000*y + 6.5000000000000000000000000000000000)*a)*a)*a)*a)*a)*a)*a;
    }

    inline double leg_f54_dxx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      return .11865250751057162342943393039327360*(133.84926223181060178824088286184803*b4 + (-803.09557339086361072944529717108816*b2 + 883.40513072994997180238982688819698*a2)*a2)*a*(2.7500000000000000000000000000000000 + (6.5000000000000000000000000000000000 + 3.7500000000000000000000000000000000*y)*y);
    }

    inline double leg_f54_dyy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      return (3.1195902949051457107623431780499134 + 1.1343964708745984402772157011090594*(6.5000000000000000000000000000000000 + 3.7500000000000000000000000000000000*y)*y-.18906607847909974004620261685150989*(7.5000000000000000000000000000000000*y + 6.5000000000000000000000000000000000)*b)*b5 + ((3.1195902949051457107623431780499134 + 1.1343964708745984402772157011090594*(6.5000000000000000000000000000000000 + 3.7500000000000000000000000000000000*y)*y + (8.5079735315594883020791177583179455*y + 7.3735770606848898618019020572088860-.70899779429662402517325981319316211*b)*b)*b4 + ((-43.674264128672039950672804492698787-15.881550592244378163881019815526832*(6.5000000000000000000000000000000000 + 3.7500000000000000000000000000000000*y)*y + 3.9703876480610945409702549538817078*(7.5000000000000000000000000000000000*y + 6.5000000000000000000000000000000000)*b)*b3 + ((-43.674264128672039950672804492698787-15.881550592244378163881019815526832*(6.5000000000000000000000000000000000 + 3.7500000000000000000000000000000000*y)*y + (-39.703876480610945409702549538817079*y-34.410026283196152688408876266974801 + 4.9629845600763681762128186923521348*b)*b)*b2 + ((65.511396193008059926009206739048179 + 23.822325888366567245821529723290247*(6.5000000000000000000000000000000000 + 3.7500000000000000000000000000000000*y)*y-11.911162944183283622910764861645123*(7.5000000000000000000000000000000000*y + 6.5000000000000000000000000000000000)*b)*b + (65.511396193008059926009206739048179 + 23.822325888366567245821529723290247*(6.5000000000000000000000000000000000 + 3.7500000000000000000000000000000000*y)*y + (35.733488832549850868732294584935370*y + 30.969023654876537419567988640277321-8.9333722081374627171830736462338426*b)*b + (65.511396193008059926009206739048179*y + 56.776543367273651935874645840508422 + 4.6793854423577185661435147670748699*a)*a)*a)*a)*a)*a)*a;
    }

    inline double leg_f54_dxy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      return (3.1195902949051457107623431780499134 + 1.1343964708745984402772157011090594*(6.5000000000000000000000000000000000 + 3.7500000000000000000000000000000000*y)*y-.18906607847909974004620261685150989*(7.5000000000000000000000000000000000*y + 6.5000000000000000000000000000000000)*b)*b5 + (7.9407752961221890819405099077634157*b4*(2.7500000000000000000000000000000000 + (6.5000000000000000000000000000000000 + 3.7500000000000000000000000000000000*y)*y) + ((-43.674264128672039950672804492698787-15.881550592244378163881019815526832*(6.5000000000000000000000000000000000 + 3.7500000000000000000000000000000000*y)*y + 3.9703876480610945409702549538817078*(7.5000000000000000000000000000000000*y + 6.5000000000000000000000000000000000)*b)*b3 + (-47.644651776733134491643059446580494*b2*(2.7500000000000000000000000000000000 + (6.5000000000000000000000000000000000 + 3.7500000000000000000000000000000000*y)*y) + ((65.511396193008059926009206739048179 + 23.822325888366567245821529723290247*(6.5000000000000000000000000000000000 + 3.7500000000000000000000000000000000*y)*y-11.911162944183283622910764861645123*(7.5000000000000000000000000000000000*y + 6.5000000000000000000000000000000000)*b)*b + (144.12507162461773183722025482590599 + 52.409116954406447940807365391238543*(6.5000000000000000000000000000000000 + 3.7500000000000000000000000000000000*y)*y + 8.7348528257344079901345608985397572*(7.5000000000000000000000000000000000*y + 6.5000000000000000000000000000000000)*a)*a)*a)*a)*a)*a;
    }

    // number 55
    inline double leg_f55(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b4 = b2*b2;
      double b6 = b4*b2;
      double b8 = b6*b2;
      return .55090590344733555051405113120224709e-1*(.10697706201272775653456441070328167*b8 + (-2.9953577363563771829678034996918866*b6 + (13.479109813603697323355115748613490*b4 + (-19.769361059952089407587503097966452*b2 + 9.1786319206920415106656264383415669*a2)*a2)*a2)*a2)*(1. + y);
    }

    inline double leg_f55_dx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      double b6 = b5*b;
      double b7 = b6*b;
      return .55090590344733555051405113120224709e-1*(-11.981430945425508731871213998767547*b6 + (107.83287850882957858684092598890792*b4 + (-237.23233271942507289105003717559742*b2 + 146.85811073107266417065002301346507*a2)*a2)*a2)*a*(1. + y);
    }

    inline double leg_f55_dy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      double b6 = b5*b;
      double b7 = b6*b;
      return (-.47147435997010740107742172080116283e-1-.47147435997010740107742172080116283e-1*y + .58934294996263425134677715100145354e-2*b)*b7 + (-.33003205197907518075419520456081398*b6*(1. + y) + ((.99009615593722554226258561368244194 + .99009615593722554226258561368244194*y-.16501602598953759037709760228040699*b)*b5 + (2.9702884678116766267877568410473258*b4*(1. + y) + ((-2.9702884678116766267877568410473258-2.9702884678116766267877568410473258*y + .74257211695291915669693921026183146*b)*b3 + (-6.5346346291856885789330650503041168*b2*(1. + y) + ((2.1782115430618961929776883501013723 + 2.1782115430618961929776883501013723*y-1.0891057715309480964888441750506861*b)*b + (4.0452500085435215012442783644739771 + 4.0452500085435215012442783644739771*y + .50565625106794018765553479555924713*a)*a)*a)*a)*a)*a)*a)*a;
    }

    inline double leg_f55_dxx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      double b6 = b5*b;
      double b7 = b6*b;
      return .55090590344733555051405113120224709e-1*(-23.962861890851017463742427997535093*b6 + (646.99727105297747152104555593344752*b4 + (-2372.3233271942507289105003717559742*b2 + 2056.0135502350172983891003221885110*a2)*a2)*a2)*(1. + y);
    }

    inline double leg_f55_dyy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      double b6 = b5*b;
      double b7 = b6*b;
      return -.94294871994021480215484344160232566e-1*b7 + ((3.9603846237489021690503424547297678 + 3.9603846237489021690503424547297678*y-.66006410395815036150839040912162796*b)*b5 + ((3.9603846237489021690503424547297678 + 3.9603846237489021690503424547297678*y + 1.9801923118744510845251712273648839*b)*b4 + ((-23.762307742493413014302054728378607-23.762307742493413014302054728378607*y + 5.9405769356233532535755136820946517*b)*b3 + ((-23.762307742493413014302054728378607-23.762307742493413014302054728378607*y-5.9405769356233532535755136820946517*b)*b2 + ((26.138538516742754315732260201216467 + 26.138538516742754315732260201216467*y-13.069269258371377157866130100608234*b)*b + (26.138538516742754315732260201216467 + 26.138538516742754315732260201216467*y + 4.3564230861237923859553767002027445*b + 8.0905000170870430024885567289479542*a)*a)*a)*a)*a)*a)*a;
    }

    inline double leg_f55_dxy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      double b6 = b5*b;
      double b7 = b6*b;
      return -.66006410395815036150839040912162796*b6*(1. + y) + ((3.9603846237489021690503424547297678 + 3.9603846237489021690503424547297678*y-.66006410395815036150839040912162796*b)*b5 + (17.821730806870059760726541046283955*b4*(1. + y) + ((-23.762307742493413014302054728378607-23.762307742493413014302054728378607*y + 5.9405769356233532535755136820946517*b)*b3 + (-65.346346291856885789330650503041168*b2*(1. + y) + ((26.138538516742754315732260201216467 + 26.138538516742754315732260201216467*y-13.069269258371377157866130100608234*b)*b + (56.633500119609301017419897102635679 + 56.633500119609301017419897102635679*y + 8.0905000170870430024885567289479542*a)*a)*a)*a)*a)*a)*a;
    }

    // ORDER 10

    // Edge functions, order 10

    // number 56
    inline double leg_f56(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b4 = b2*b2;
      double b6 = b4*b2;
      double b8 = b6*b2;
      double b10 = b8*b2;
      return -.82303806344571752837660497036322993e-4*b10 + (.37036712855057288776947223666345347e-2*b8 + (-.27160256093708678436427964021986588e-1*b6 + (.70616665843642563934712706457165128e-1*b4 + (-.75660713403902747072906471204105494e-1*b2 + .28582936174807704449764666899328742e-1*a2)*a2)*a2)*a2)*a2;
    }

    inline double leg_f56_dx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      double b6 = b5*b;
      double b7 = b6*b;
      double b8 = b7*b;
      double b9 = b8*b;
      return (.14814685142022915510778889466538139e-1*b8 + (-.21728204874966942749142371217589270*b6 + (.84739999012371076721655247748598154*b4 + (-1.2105714144624439531665035392656879*b2 + .57165872349615408899529333798657485*a2)*a2)*a2)*a2)*a;
    }

    inline double leg_f56_dy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      double b6 = b5*b;
      double b7 = b6*b;
      double b8 = b7*b;
      double b9 = b8*b;
      return .82303806344571752837660497036322993e-3*b9 + (.74073425710114577553894447332690694e-2*b8 + (-.29629370284045831021557778933076278e-1*b7 + (-.10864102437483471374571185608794635*b6 + (.16296153656225207061856778413191953*b5 + (.42369999506185538360827623874299077*b4 + (-.28246666337457025573885082582866051*b3 + (-.60528570723122197658325176963284395*b2 + (.15132142680780549414581294240821099*b + .28582936174807704449764666899328742*a)*a)*a)*a)*a)*a)*a)*a)*a;
    }

    inline double leg_f56_dxx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      double b6 = b5*b;
      double b7 = b6*b;
      double b8 = b7*b;
      return .29629370284045831021557778933076278e-1*b8 + (-1.3036922924980165649485422730553562*b6 + (8.4739999012371076721655247748598154*b4 + (-16.947999802474215344331049549719631*b2 + 10.289857022930773601915280083758347*a2)*a2)*a2)*a2;
    }

    inline double leg_f56_dyy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      double b6 = b5*b;
      double b7 = b6*b;
      double b8 = b7*b;
      return (-.11851748113618332408623111573230511*b7 + (-.11851748113618332408623111573230511*b6 + (1.3036922924980165649485422730553562*b5 + (1.3036922924980165649485422730553562*b4 + (-3.3895999604948430688662099099439261*b3 + (-3.3895999604948430688662099099439261*b2 + (2.4211428289248879063330070785313758*b + 2.4211428289248879063330070785313758*a)*a)*a)*a)*a)*a)*a)*a;
    }

    inline double leg_f56_dxy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      double b6 = b5*b;
      double b7 = b6*b;
      double b8 = b7*b;
      return .14814685142022915510778889466538139e-1*b8 + (-.11851748113618332408623111573230511*b7 + (-.65184614624900828247427113652767811*b6 + (1.3036922924980165649485422730553562*b5 + (4.2369999506185538360827623874299077*b4 + (-3.3895999604948430688662099099439261*b3 + (-8.4739999012371076721655247748598154*b2 + (2.4211428289248879063330070785313758*b + 5.1449285114653868009576400418791736*a)*a)*a)*a)*a)*a)*a)*a;
    }

    // number 57
    inline double leg_f57(double x, double y)
    {
      double y2 = y*y;
      return -.12039871099548782129394906995027821e-1*(x + 1.)*(7. + (-308. + (2002. + (-4004. + 2431.*y2)*y2)*y2)*y2)*(1. + y);
    }

    inline double leg_f57_dx(double x, double y)
    {
      double y2 = y*y;
      return -.12039871099548782129394906995027821e-1*(7. + (-308. + (2002. + (-4004. + 2431.*y2)*y2)*y2)*y2)*(1. + y);
    }

    inline double leg_f57_dy(double x, double y)
    {
      double y2 = y*y;
      return -.12039871099548782129394906995027821e-1*(x + 1.)*(7. + (-616. + (-924. + (8008. + (10010. + (-24024. + (-28028. + (19448. + 21879.*y)*y)*y)*y)*y)*y)*y)*y);
    }

    inline double leg_f57_dxx(double x, double y)
    {
      return 0.;
    }

    inline double leg_f57_dyy(double x, double y)
    {
      return -1.0595086567602928273867518155624482*(x + 1.)*(-7. + (-21. + (273. + (455. + (-1365. + (-1911. + (1547. + 1989.*y)*y)*y)*y)*y)*y)*y);
    }

    inline double leg_f57_dxy(double x, double y)
    {
      return -.84279097696841474905764348965194744e-1-.12039871099548782129394906995027821e-1*(-616. + (-924. + (8008. + (10010. + (-24024. + (-28028. + (19448. + 21879.*y)*y)*y)*y)*y)*y)*y)*y;
    }

    // number 58
    inline double leg_f58(double x, double y)
    {
      double y2 = y*y;
      return .12039871099548782129394906995027821e-1*(x + y)*(7. + (-308. + (2002. + (-4004. + 2431.*y2)*y2)*y2)*y2)*(1. + y);
    }

    inline double leg_f58_dx(double x, double y)
    {
      double y2 = y*y;
      return .12039871099548782129394906995027821e-1*(7. + (-308. + (2002. + (-4004. + 2431.*y2)*y2)*y2)*y2)*(1. + y);
    }

    inline double leg_f58_dy(double x, double y)
    {
      double y2 = y*y;
      return .84279097696841474905764348965194744e-1 + .12039871099548782129394906995027821e-1*(14. + (-924. + (-1232. + (10010. + (12012. + (-28028. + (-32032. + (21879. + 24310.*y)*y)*y)*y)*y)*y)*y)*y)*y + .12039871099548782129394906995027821e-1*(7. + (-616. + (-924. + (8008. + (10010. + (-24024. + (-28028. + (19448. + 21879.*y)*y)*y)*y)*y)*y)*y)*y)*x;
    }

    inline double leg_f58_dxx(double x, double y)
    {
      return 0.;
    }

    inline double leg_f58_dyy(double x, double y)
    {
      return .16855819539368294981152869793038949 + .24079742199097564258789813990055641e-1*(-924. + (-1848. + (20020. + (30030. + (-84084. + (-112112. + (87516. + 109395.*y)*y)*y)*y)*y)*y)*y)*y + .24079742199097564258789813990055641e-1*(-308. + (-924. + (12012. + (20020. + (-60060. + (-84084. + (68068. + 87516.*y)*y)*y)*y)*y)*y)*y)*x;
    }

    inline double leg_f58_dxy(double x, double y)
    {
      return .84279097696841474905764348965194744e-1 + .12039871099548782129394906995027821e-1*(-616. + (-924. + (8008. + (10010. + (-24024. + (-28028. + (19448. + 21879.*y)*y)*y)*y)*y)*y)*y)*y;
    }

    // Bubble functions, order 10

    // number 59
    inline double leg_f59(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double y2 = y*y;
      return 1.8781826691835577799972514353411927*(.61237243569579452454932101867647285*a2-.61237243569579452454932101867647285*b2)*(-.54687500000000000000000000000000000e-1 + (-.65625000000000000000000000000000000 + (1.2031250000000000000000000000000000 + (9.6250000000000000000000000000000000 + (-31.281250000000000000000000000000000 + (-15.640625000000000000000000000000000 + (26.812500000000000000000000000000000 + 18.992187500000000000000000000000000*y)*y)*y)*y2)*y)*y)*y);
    }

    inline double leg_f59_dx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double y2 = y*y;
      return 4.6005891832382558280858377040928971*a*(-.54687500000000000000000000000000000e-1 + (-.65625000000000000000000000000000000 + (1.2031250000000000000000000000000000 + (9.6250000000000000000000000000000000 + (-31.281250000000000000000000000000000 + (-15.640625000000000000000000000000000 + (26.812500000000000000000000000000000 + 18.992187500000000000000000000000000*y)*y)*y)*y2)*y)*y)*y);
    }

    inline double leg_f59_dy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double y2 = y*y;
      return (-.12579736047917105779922212472129016 + 2.3002945916191279140429188520464486*(-.65625000000000000000000000000000000 + (1.2031250000000000000000000000000000 + (9.6250000000000000000000000000000000 + (-31.281250000000000000000000000000000 + (-15.640625000000000000000000000000000 + (26.812500000000000000000000000000000 + 18.992187500000000000000000000000000*y)*y)*y)*y2)*y)*y)*y-1.1501472958095639570214594260232243*(-.65625000000000000000000000000000000 + (2.4062500000000000000000000000000000 + (28.875000000000000000000000000000000 + (-156.40625000000000000000000000000000 + (-93.843750000000000000000000000000000 + (187.68750000000000000000000000000000 + 151.93750000000000000000000000000000*y)*y)*y)*y2)*y)*y)*b)*b + (-.12579736047917105779922212472129016 + 2.3002945916191279140429188520464486*(-.65625000000000000000000000000000000 + (1.2031250000000000000000000000000000 + (9.6250000000000000000000000000000000 + (-31.281250000000000000000000000000000 + (-15.640625000000000000000000000000000 + (26.812500000000000000000000000000000 + 18.992187500000000000000000000000000*y)*y)*y)*y2)*y)*y)*y + 1.1501472958095639570214594260232243*(-.65625000000000000000000000000000000 + (2.4062500000000000000000000000000000 + (28.875000000000000000000000000000000 + (-156.40625000000000000000000000000000 + (-93.843750000000000000000000000000000 + (187.68750000000000000000000000000000 + 151.93750000000000000000000000000000*y)*y)*y)*y2)*y)*y)*a)*a;
    }

    inline double leg_f59_dxx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double y2 = y*y;
      return -.50318944191668423119688849888516062 + 9.2011783664765116561716754081857942*(-.65625000000000000000000000000000000 + (1.2031250000000000000000000000000000 + (9.6250000000000000000000000000000000 + (-31.281250000000000000000000000000000 + (-15.640625000000000000000000000000000 + (26.812500000000000000000000000000000 + 18.992187500000000000000000000000000*y)*y)*y)*y2)*y)*y)*y;
    }

    inline double leg_f59_dyy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double y2 = y*y;
      return (-3.0191366515001053871813309933109638 + 4.6005891832382558280858377040928971*(2.4062500000000000000000000000000000 + (28.875000000000000000000000000000000 + (-156.40625000000000000000000000000000 + (-93.843750000000000000000000000000000 + (187.68750000000000000000000000000000 + 151.93750000000000000000000000000000*y)*y)*y)*y2)*y)*y-1.1501472958095639570214594260232243*(2.4062500000000000000000000000000000 + (57.750000000000000000000000000000000 + (-625.62500000000000000000000000000000 + (-469.21875000000000000000000000000000 + (1126.1250000000000000000000000000000 + 1063.5625000000000000000000000000000*y)*y)*y)*y2)*y)*b)*b + (-3.0191366515001053871813309933109638 + 4.6005891832382558280858377040928971*(2.4062500000000000000000000000000000 + (28.875000000000000000000000000000000 + (-156.40625000000000000000000000000000 + (-93.843750000000000000000000000000000 + (187.68750000000000000000000000000000 + 151.93750000000000000000000000000000*y)*y)*y)*y2)*y)*y + 1.1501472958095639570214594260232243*(2.4062500000000000000000000000000000 + (57.750000000000000000000000000000000 + (-625.62500000000000000000000000000000 + (-469.21875000000000000000000000000000 + (1126.1250000000000000000000000000000 + 1063.5625000000000000000000000000000*y)*y)*y)*y2)*y)*a)*a;
    }

    inline double leg_f59_dxy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double y2 = y*y;
      return -.25159472095834211559844424944258031 + 4.6005891832382558280858377040928971*(-.65625000000000000000000000000000000 + (1.2031250000000000000000000000000000 + (9.6250000000000000000000000000000000 + (-31.281250000000000000000000000000000 + (-15.640625000000000000000000000000000 + (26.812500000000000000000000000000000 + 18.992187500000000000000000000000000*y)*y)*y)*y2)*y)*y)*y + 4.6005891832382558280858377040928971*(-.65625000000000000000000000000000000 + (2.4062500000000000000000000000000000 + (28.875000000000000000000000000000000 + (-156.40625000000000000000000000000000 + (-93.843750000000000000000000000000000 + (187.68750000000000000000000000000000 + 151.93750000000000000000000000000000*y)*y)*y)*y2)*y)*y)*a;
    }

    // number 60
    inline double leg_f60(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      return 1.7465387901340515412774319981019878*(-.79056941504209483299972338610817962*b2 + .79056941504209483299972338610817962*a2)*a*(.62500000000000000000000000000000000e-1 + (1.5625000000000000000000000000000000 + (1.5000000000000000000000000000000000 + (-16.250000000000000000000000000000000 + (-28.437500000000000000000000000000000 + (17.062500000000000000000000000000000 + (56.875000000000000000000000000000000 + 27.625000000000000000000000000000000*y)*y)*y)*y)*y)*y)*y);
    }

    inline double leg_f60_dx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      return 1.7465387901340515412774319981019878*(4.7434164902525689979983403166490778*a2-1.5811388300841896659994467722163592*b2)*(.62500000000000000000000000000000000e-1 + (1.5625000000000000000000000000000000 + (1.5000000000000000000000000000000000 + (-16.250000000000000000000000000000000 + (-28.437500000000000000000000000000000 + (17.062500000000000000000000000000000 + (56.875000000000000000000000000000000 + 27.625000000000000000000000000000000*y)*y)*y)*y)*y)*y)*y);
    }

    inline double leg_f60_dy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      return -1.3807601496646051572657795977932815*b2*(.62500000000000000000000000000000000e-1 + (1.5625000000000000000000000000000000 + (1.5000000000000000000000000000000000 + (-16.250000000000000000000000000000000 + (-28.437500000000000000000000000000000 + (17.062500000000000000000000000000000 + (56.875000000000000000000000000000000 + 27.625000000000000000000000000000000*y)*y)*y)*y)*y)*y)*y) + ((.17259501870807564465822244972416018 + 2.7615202993292103145315591955865629*(1.5625000000000000000000000000000000 + (1.5000000000000000000000000000000000 + (-16.250000000000000000000000000000000 + (-28.437500000000000000000000000000000 + (17.062500000000000000000000000000000 + (56.875000000000000000000000000000000 + 27.625000000000000000000000000000000*y)*y)*y)*y)*y)*y)*y-1.3807601496646051572657795977932815*(1.5625000000000000000000000000000000 + (3. + (-48.750000000000000000000000000000000 + (-113.75000000000000000000000000000000 + (85.312500000000000000000000000000000 + (341.25000000000000000000000000000000 + 193.37500000000000000000000000000000*y)*y)*y)*y)*y)*y)*b)*b + (.25889252806211346698733367458624028 + 4.1422804489938154717973387933798444*(1.5625000000000000000000000000000000 + (1.5000000000000000000000000000000000 + (-16.250000000000000000000000000000000 + (-28.437500000000000000000000000000000 + (17.062500000000000000000000000000000 + (56.875000000000000000000000000000000 + 27.625000000000000000000000000000000*y)*y)*y)*y)*y)*y)*y + 1.3807601496646051572657795977932815*(1.5625000000000000000000000000000000 + (3. + (-48.750000000000000000000000000000000 + (-113.75000000000000000000000000000000 + (85.312500000000000000000000000000000 + (341.25000000000000000000000000000000 + 193.37500000000000000000000000000000*y)*y)*y)*y)*y)*y)*a)*a)*a;
    }

    inline double leg_f60_dxx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return 33.138243591950523774378710347038755*a*(.62500000000000000000000000000000000e-1 + (1.5625000000000000000000000000000000 + (1.5000000000000000000000000000000000 + (-16.250000000000000000000000000000000 + (-28.437500000000000000000000000000000 + (17.062500000000000000000000000000000 + (56.875000000000000000000000000000000 + 27.625000000000000000000000000000000*y)*y)*y)*y)*y)*y)*y);
    }

    inline double leg_f60_dyy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return (.34519003741615128931644489944832037 + 5.5230405986584206290631183911731259*(1.5625000000000000000000000000000000 + (1.5000000000000000000000000000000000 + (-16.250000000000000000000000000000000 + (-28.437500000000000000000000000000000 + (17.062500000000000000000000000000000 + (56.875000000000000000000000000000000 + 27.625000000000000000000000000000000*y)*y)*y)*y)*y)*y)*y-2.7615202993292103145315591955865629*(1.5625000000000000000000000000000000 + (3. + (-48.750000000000000000000000000000000 + (-113.75000000000000000000000000000000 + (85.312500000000000000000000000000000 + (341.25000000000000000000000000000000 + 193.37500000000000000000000000000000*y)*y)*y)*y)*y)*y)*b)*b + (.34519003741615128931644489944832037 + 5.5230405986584206290631183911731259*(1.5625000000000000000000000000000000 + (1.5000000000000000000000000000000000 + (-16.250000000000000000000000000000000 + (-28.437500000000000000000000000000000 + (17.062500000000000000000000000000000 + (56.875000000000000000000000000000000 + 27.625000000000000000000000000000000*y)*y)*y)*y)*y)*y)*y + (8.6297509354037822329111224862080092 + 5.5230405986584206290631183911731259*(3. + (-48.750000000000000000000000000000000 + (-113.75000000000000000000000000000000 + (85.312500000000000000000000000000000 + (341.25000000000000000000000000000000 + 193.37500000000000000000000000000000*y)*y)*y)*y)*y)*y-1.3807601496646051572657795977932815*(3. + (-97.500000000000000000000000000000000 + (-341.25000000000000000000000000000000 + (341.25000000000000000000000000000000 + (1706.2500000000000000000000000000000 + 1160.2500000000000000000000000000000*y)*y)*y)*y)*y)*b)*b + (12.944626403105673349366683729312014 + 8.2845608979876309435946775867596891*(3. + (-48.750000000000000000000000000000000 + (-113.75000000000000000000000000000000 + (85.312500000000000000000000000000000 + (341.25000000000000000000000000000000 + 193.37500000000000000000000000000000*y)*y)*y)*y)*y)*y + 1.3807601496646051572657795977932815*(3. + (-97.500000000000000000000000000000000 + (-341.25000000000000000000000000000000 + (341.25000000000000000000000000000000 + (1706.2500000000000000000000000000000 + 1160.2500000000000000000000000000000*y)*y)*y)*y)*y)*a)*a)*a;
    }

    inline double leg_f60_dxy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      return (.34519003741615128931644489944832037 + 5.5230405986584206290631183911731259*(1.5625000000000000000000000000000000 + (1.5000000000000000000000000000000000 + (-16.250000000000000000000000000000000 + (-28.437500000000000000000000000000000 + (17.062500000000000000000000000000000 + (56.875000000000000000000000000000000 + 27.625000000000000000000000000000000*y)*y)*y)*y)*y)*y)*y-2.7615202993292103145315591955865629*(1.5625000000000000000000000000000000 + (3. + (-48.750000000000000000000000000000000 + (-113.75000000000000000000000000000000 + (85.312500000000000000000000000000000 + (341.25000000000000000000000000000000 + 193.37500000000000000000000000000000*y)*y)*y)*y)*y)*y)*b)*b + (1.0355701122484538679493346983449611 + 16.569121795975261887189355173519378*(1.5625000000000000000000000000000000 + (1.5000000000000000000000000000000000 + (-16.250000000000000000000000000000000 + (-28.437500000000000000000000000000000 + (17.062500000000000000000000000000000 + (56.875000000000000000000000000000000 + 27.625000000000000000000000000000000*y)*y)*y)*y)*y)*y)*y + 8.2845608979876309435946775867596891*(1.5625000000000000000000000000000000 + (3. + (-48.750000000000000000000000000000000 + (-113.75000000000000000000000000000000 + (85.312500000000000000000000000000000 + (341.25000000000000000000000000000000 + 193.37500000000000000000000000000000*y)*y)*y)*y)*y)*y)*a)*a;
    }

    // number 61
    inline double leg_f61(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b4 = b2*b2;
      double y2 = y*y;
      return 1.1071816700297210425642255794583961*(.23385358667337133659898429576978433*b4 + (-1.4031215200402280195939057746187060*b2 + 1.1692679333668566829949214788489217*a2)*a2)*(.83333333333333333333333333333333333e-1 + (-2.6250000000000000000000000000000000 + (-12.187500000000000000000000000000000 + (56.875000000000000000000000000000000 + (79.625000000000000000000000000000000 + 32.229166666666666666666666666666667*y)*y)*y2)*y)*y);
    }

    inline double leg_f61_dx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double y2 = y*y;
      return 1.1071816700297210425642255794583961*(-5.6124860801609120783756230984748240*b2 + 9.3541434669348534639593718307913732*a2)*a*(.83333333333333333333333333333333333e-1 + (-2.6250000000000000000000000000000000 + (-12.187500000000000000000000000000000 + (56.875000000000000000000000000000000 + (79.625000000000000000000000000000000 + 32.229166666666666666666666666666667*y)*y)*y2)*y)*y);
    }

    inline double leg_f61_dy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double y2 = y*y;
      return (-.86306134878487797796370546040008498e-1-1.0356736185418535735564465524801020*(-2.6250000000000000000000000000000000 + (-12.187500000000000000000000000000000 + (56.875000000000000000000000000000000 + (79.625000000000000000000000000000000 + 32.229166666666666666666666666666667*y)*y)*y2)*y)*y + .25891840463546339338911163812002550*(-2.6250000000000000000000000000000000 + (-24.375000000000000000000000000000000 + (227.50000000000000000000000000000000 + (398.12500000000000000000000000000000 + 193.37500000000000000000000000000000*y)*y)*y2)*y)*b)*b3 + (-3.1070208556255607206693396574403060*b2*(.83333333333333333333333333333333333e-1 + (-2.6250000000000000000000000000000000 + (-12.187500000000000000000000000000000 + (56.875000000000000000000000000000000 + (79.625000000000000000000000000000000 + 32.229166666666666666666666666666667*y)*y)*y2)*y)*y) + ((.25891840463546339338911163812002550 + 3.1070208556255607206693396574403060*(-2.6250000000000000000000000000000000 + (-12.187500000000000000000000000000000 + (56.875000000000000000000000000000000 + (79.625000000000000000000000000000000 + 32.229166666666666666666666666666667*y)*y)*y2)*y)*y-1.5535104278127803603346698287201530*(-2.6250000000000000000000000000000000 + (-24.375000000000000000000000000000000 + (227.50000000000000000000000000000000 + (398.12500000000000000000000000000000 + 193.37500000000000000000000000000000*y)*y)*y2)*y)*b)*b + (.43153067439243898898185273020004250 + 5.1783680927092678677822327624005100*(-2.6250000000000000000000000000000000 + (-12.187500000000000000000000000000000 + (56.875000000000000000000000000000000 + (79.625000000000000000000000000000000 + 32.229166666666666666666666666666667*y)*y)*y2)*y)*y + 1.2945920231773169669455581906001275*(-2.6250000000000000000000000000000000 + (-24.375000000000000000000000000000000 + (227.50000000000000000000000000000000 + (398.12500000000000000000000000000000 + 193.37500000000000000000000000000000*y)*y)*y2)*y)*a)*a)*a)*a;
    }

    inline double leg_f61_dxx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double y2 = y*y;
      return 1.1071816700297210425642255794583961*(56.124860801609120783756230984748240*a2-11.224972160321824156751246196949648*b2)*(.83333333333333333333333333333333333e-1 + (-2.6250000000000000000000000000000000 + (-12.187500000000000000000000000000000 + (56.875000000000000000000000000000000 + (79.625000000000000000000000000000000 + 32.229166666666666666666666666666667*y)*y)*y2)*y)*y);
    }

    inline double leg_f61_dyy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double y2 = y*y;
      return (5.4372864973447312611713444005205354-2.0713472370837071471128931049602040*(-24.375000000000000000000000000000000 + (227.50000000000000000000000000000000 + (398.12500000000000000000000000000000 + 193.37500000000000000000000000000000*y)*y)*y2)*y + .25891840463546339338911163812002550*(-24.375000000000000000000000000000000 + (682.50000000000000000000000000000000 + (1592.5000000000000000000000000000000 + 966.87500000000000000000000000000000*y)*y)*y2)*b)*b3 + ((1.0356736185418535735564465524801020 + 12.428083422502242882677358629761224*(-2.6250000000000000000000000000000000 + (-12.187500000000000000000000000000000 + (56.875000000000000000000000000000000 + (79.625000000000000000000000000000000 + 32.229166666666666666666666666666667*y)*y)*y2)*y)*y-6.2140417112511214413386793148806119*(-2.6250000000000000000000000000000000 + (-24.375000000000000000000000000000000 + (227.50000000000000000000000000000000 + (398.12500000000000000000000000000000 + 193.37500000000000000000000000000000*y)*y)*y2)*y)*b)*b + (1.0356736185418535735564465524801020 + 12.428083422502242882677358629761224*(-2.6250000000000000000000000000000000 + (-12.187500000000000000000000000000000 + (56.875000000000000000000000000000000 + (79.625000000000000000000000000000000 + 32.229166666666666666666666666666667*y)*y)*y2)*y)*y + (-16.311859492034193783514033201561606 + 6.2140417112511214413386793148806119*(-24.375000000000000000000000000000000 + (227.50000000000000000000000000000000 + (398.12500000000000000000000000000000 + 193.37500000000000000000000000000000*y)*y)*y2)*y-1.5535104278127803603346698287201530*(-24.375000000000000000000000000000000 + (682.50000000000000000000000000000000 + (1592.5000000000000000000000000000000 + 966.87500000000000000000000000000000*y)*y)*y2)*b)*b + (-27.186432486723656305856722002602676 + 10.356736185418535735564465524801020*(-24.375000000000000000000000000000000 + (227.50000000000000000000000000000000 + (398.12500000000000000000000000000000 + 193.37500000000000000000000000000000*y)*y)*y2)*y + 1.2945920231773169669455581906001275*(-24.375000000000000000000000000000000 + (682.50000000000000000000000000000000 + (1592.5000000000000000000000000000000 + 966.87500000000000000000000000000000*y)*y)*y2)*a)*a)*a)*a;
    }

    inline double leg_f61_dxy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double y2 = y*y;
      return -6.2140417112511214413386793148806119*b2*(.83333333333333333333333333333333333e-1 + (-2.6250000000000000000000000000000000 + (-12.187500000000000000000000000000000 + (56.875000000000000000000000000000000 + (79.625000000000000000000000000000000 + 32.229166666666666666666666666666667*y)*y)*y2)*y)*y) + ((1.0356736185418535735564465524801020 + 12.428083422502242882677358629761224*(-2.6250000000000000000000000000000000 + (-12.187500000000000000000000000000000 + (56.875000000000000000000000000000000 + (79.625000000000000000000000000000000 + 32.229166666666666666666666666666667*y)*y)*y2)*y)*y-6.2140417112511214413386793148806119*(-2.6250000000000000000000000000000000 + (-24.375000000000000000000000000000000 + (227.50000000000000000000000000000000 + (398.12500000000000000000000000000000 + 193.37500000000000000000000000000000*y)*y)*y2)*y)*b)*b + (2.5891840463546339338911163812002550 + 31.070208556255607206693396574403060*(-2.6250000000000000000000000000000000 + (-12.187500000000000000000000000000000 + (56.875000000000000000000000000000000 + (79.625000000000000000000000000000000 + 32.229166666666666666666666666666667*y)*y)*y2)*y)*y + 10.356736185418535735564465524801020*(-2.6250000000000000000000000000000000 + (-24.375000000000000000000000000000000 + (227.50000000000000000000000000000000 + (398.12500000000000000000000000000000 + 193.37500000000000000000000000000000*y)*y)*y2)*y)*a)*a)*a;
    }

    // number 62
    inline double leg_f62(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b4 = b2*b2;
      return .60500365815888498128641266942176687*(.79549512883486596495094990736795518*b4 + (-2.6516504294495532165031663578931839*b2 + 1.8561553006146872515522164505252288*a2)*a2)*a*(-.75000000000000000000000000000000000 + (-.75000000000000000000000000000000000 + (21. + (70. + (78.750000000000000000000000000000000 + 29.750000000000000000000000000000000*y)*y)*y)*y)*y);
    }

    inline double leg_f62_dx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      return .60500365815888498128641266942176687*(1.5909902576697319299018998147359104*b4 + (-15.909902576697319299018998147359104*b2 + 18.561553006146872515522164505252288*a2)*a2)*(-.75000000000000000000000000000000000 + (-.75000000000000000000000000000000000 + (21. + (70. + (78.750000000000000000000000000000000 + 29.750000000000000000000000000000000*y)*y)*y)*y)*y);
    }

    inline double leg_f62_dy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      return .48127746299266741539338459887275870*b4*(-.75000000000000000000000000000000000 + (-.75000000000000000000000000000000000 + (21. + (70. + (78.750000000000000000000000000000000 + 29.750000000000000000000000000000000*y)*y)*y)*y)*y) + ((1.4438323889780022461801537966182761-1.9251098519706696615735383954910348*(-.75000000000000000000000000000000000 + (21. + (70. + (78.750000000000000000000000000000000 + 29.750000000000000000000000000000000*y)*y)*y)*y)*y + .48127746299266741539338459887275870*(-.75000000000000000000000000000000000 + (42. + (210. + (315. + 148.75000000000000000000000000000000*y)*y)*y)*y)*b)*b3 + (-4.8127746299266741539338459887275870*b2*(-.75000000000000000000000000000000000 + (-.75000000000000000000000000000000000 + (21. + (70. + (78.750000000000000000000000000000000 + 29.750000000000000000000000000000000*y)*y)*y)*y)*y) + ((-2.4063873149633370769669229943637935 + 3.2085164199511161026225639924850580*(-.75000000000000000000000000000000000 + (21. + (70. + (78.750000000000000000000000000000000 + 29.750000000000000000000000000000000*y)*y)*y)*y)*y-1.6042582099755580513112819962425290*(-.75000000000000000000000000000000000 + (42. + (210. + (315. + 148.75000000000000000000000000000000*y)*y)*y)*y)*b)*b + (-4.2111778011858398846921152401366384 + 5.6149037349144531795894869868488513*(-.75000000000000000000000000000000000 + (21. + (70. + (78.750000000000000000000000000000000 + 29.750000000000000000000000000000000*y)*y)*y)*y)*y + 1.1229807469828906359178973973697703*(-.75000000000000000000000000000000000 + (42. + (210. + (315. + 148.75000000000000000000000000000000*y)*y)*y)*y)*a)*a)*a)*a)*a;
    }

    inline double leg_f62_dxx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      return .60500365815888498128641266942176687*(-63.639610306789277196075992589436414*b2 + 148.49242404917498012417731604201830*a2)*a*(-.75000000000000000000000000000000000 + (-.75000000000000000000000000000000000 + (21. + (70. + (78.750000000000000000000000000000000 + 29.750000000000000000000000000000000*y)*y)*y)*y)*y);
    }

    inline double leg_f62_dyy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      return (2.8876647779560044923603075932365521-3.8502197039413393231470767909820695*(-.75000000000000000000000000000000000 + (21. + (70. + (78.750000000000000000000000000000000 + 29.750000000000000000000000000000000*y)*y)*y)*y)*y + .96255492598533483078676919774551738*(-.75000000000000000000000000000000000 + (42. + (210. + (315. + 148.75000000000000000000000000000000*y)*y)*y)*y)*b)*b3 + ((2.8876647779560044923603075932365521-3.8502197039413393231470767909820695*(-.75000000000000000000000000000000000 + (21. + (70. + (78.750000000000000000000000000000000 + 29.750000000000000000000000000000000*y)*y)*y)*y)*y + (2.8876647779560044923603075932365521-3.8502197039413393231470767909820695*(42. + (210. + (315. + 148.75000000000000000000000000000000*y)*y)*y)*y + .48127746299266741539338459887275870*(42. + (420. + (945. + 595.*y)*y)*y)*b)*b)*b2 + ((-14.438323889780022461801537966182761 + 19.251098519706696615735383954910348*(-.75000000000000000000000000000000000 + (21. + (70. + (78.750000000000000000000000000000000 + 29.750000000000000000000000000000000*y)*y)*y)*y)*y-9.6255492598533483078676919774551738*(-.75000000000000000000000000000000000 + (42. + (210. + (315. + 148.75000000000000000000000000000000*y)*y)*y)*y)*b)*b + (-14.438323889780022461801537966182761 + 19.251098519706696615735383954910348*(-.75000000000000000000000000000000000 + (21. + (70. + (78.750000000000000000000000000000000 + 29.750000000000000000000000000000000*y)*y)*y)*y)*y + (-4.8127746299266741539338459887275868 + 6.4170328399022322052451279849701158*(42. + (210. + (315. + 148.75000000000000000000000000000000*y)*y)*y)*y-1.6042582099755580513112819962425290*(42. + (420. + (945. + 595.*y)*y)*y)*b)*b + (-8.4223556023716797693842304802732771 + 11.229807469828906359178973973697703*(42. + (210. + (315. + 148.75000000000000000000000000000000*y)*y)*y)*y + 1.1229807469828906359178973973697703*(42. + (420. + (945. + 595.*y)*y)*y)*a)*a)*a)*a)*a;
    }

    inline double leg_f62_dxy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      return (2.8876647779560044923603075932365521-3.8502197039413393231470767909820695*(-.75000000000000000000000000000000000 + (21. + (70. + (78.750000000000000000000000000000000 + 29.750000000000000000000000000000000*y)*y)*y)*y)*y + .96255492598533483078676919774551738*(-.75000000000000000000000000000000000 + (42. + (210. + (315. + 148.75000000000000000000000000000000*y)*y)*y)*y)*b)*b3 + (-19.251098519706696615735383954910348*b2*(-.75000000000000000000000000000000000 + (-.75000000000000000000000000000000000 + (21. + (70. + (78.750000000000000000000000000000000 + 29.750000000000000000000000000000000*y)*y)*y)*y)*y) + ((-14.438323889780022461801537966182761 + 19.251098519706696615735383954910348*(-.75000000000000000000000000000000000 + (21. + (70. + (78.750000000000000000000000000000000 + 29.750000000000000000000000000000000*y)*y)*y)*y)*y-9.6255492598533483078676919774551738*(-.75000000000000000000000000000000000 + (42. + (210. + (315. + 148.75000000000000000000000000000000*y)*y)*y)*y)*b)*b + (-33.689422409486719077536921921093108 + 44.919229879315625436715895894790810*(-.75000000000000000000000000000000000 + (21. + (70. + (78.750000000000000000000000000000000 + 29.750000000000000000000000000000000*y)*y)*y)*y)*y + 11.229807469828906359178973973697703*(-.75000000000000000000000000000000000 + (42. + (210. + (315. + 148.75000000000000000000000000000000*y)*y)*y)*y)*a)*a)*a)*a;
    }

    // number 63
    inline double leg_f63(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b4 = b2*b2;
      double b6 = b4*b2;
      return .30530019864318345929552039857796281*(-.14657549249448217358017594104826457*b6 + (2.1986323874172326037026391157239686*b4 + (-5.1301422373068760753061579366892600*b2 + 3.0780853423841256451836947620135560*a2)*a2)*a2)*(1.5000000000000000000000000000000000 + (16.500000000000000000000000000000000 + (48.750000000000000000000000000000000 + (55. + 21.250000000000000000000000000000000*y)*y)*y)*y);
    }

    inline double leg_f63_dx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      return .30530019864318345929552039857796281*(8.7945295496689304148105564628958743*b4 + (-41.041137898455008602449263493514080*b2 + 36.937024108609507742204337144162672*a2)*a2)*a*(1.5000000000000000000000000000000000 + (16.500000000000000000000000000000000 + (48.750000000000000000000000000000000 + (55. + 21.250000000000000000000000000000000*y)*y)*y)*y);
    }

    inline double leg_f63_dy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      return (.40274574277309068430839343011964242 + .26849716184872712287226228674642828*(16.500000000000000000000000000000000 + (48.750000000000000000000000000000000 + (55. + 21.250000000000000000000000000000000*y)*y)*y)*y-.44749526974787853812043714457738048e-1*(16.500000000000000000000000000000000 + (97.500000000000000000000000000000000 + (165. + 85.*y)*y)*y)*b)*b5 + (1.3424858092436356143613114337321414*b4*(1.5000000000000000000000000000000000 + (16.500000000000000000000000000000000 + (48.750000000000000000000000000000000 + (55. + 21.250000000000000000000000000000000*y)*y)*y)*y) + ((-4.0274574277309068430839343011964242-2.6849716184872712287226228674642828*(16.500000000000000000000000000000000 + (48.750000000000000000000000000000000 + (55. + 21.250000000000000000000000000000000*y)*y)*y)*y + .67124290462181780718065571686607068*(16.500000000000000000000000000000000 + (97.500000000000000000000000000000000 + (165. + 85.*y)*y)*y)*b)*b3 + (-6.2649337764702995336861200240833266*b2*(1.5000000000000000000000000000000000 + (16.500000000000000000000000000000000 + (48.750000000000000000000000000000000 + (55. + 21.250000000000000000000000000000000*y)*y)*y)*y) + ((4.6987003323527246502645900180624949 + 3.1324668882351497668430600120416633*(16.500000000000000000000000000000000 + (48.750000000000000000000000000000000 + (55. + 21.250000000000000000000000000000000*y)*y)*y)*y-1.5662334441175748834215300060208316*(16.500000000000000000000000000000000 + (97.500000000000000000000000000000000 + (165. + 85.*y)*y)*y)*b)*b + (8.4576605982349043704762620325124910 + 5.6384403988232695803175080216749940*(16.500000000000000000000000000000000 + (48.750000000000000000000000000000000 + (55. + 21.250000000000000000000000000000000*y)*y)*y)*y + .93974006647054493005291800361249898*(16.500000000000000000000000000000000 + (97.500000000000000000000000000000000 + (165. + 85.*y)*y)*y)*a)*a)*a)*a)*a)*a;
    }

    inline double leg_f63_dxx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      return .30530019864318345929552039857796281*(17.589059099337860829621112925791749*b4 + (-246.24682739073005161469558096108448*b2 + 369.37024108609507742204337144162672*a2)*a2)*(1.5000000000000000000000000000000000 + (16.500000000000000000000000000000000 + (48.750000000000000000000000000000000 + (55. + 21.250000000000000000000000000000000*y)*y)*y)*y);
    }

    inline double leg_f63_dyy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      return (8.8604063410079950547846554626321332 + .53699432369745424574452457349285655*(97.500000000000000000000000000000000 + (165. + 85.*y)*y)*y-.44749526974787853812043714457738048e-1*(97.500000000000000000000000000000000 + (330. + 255.*y)*y)*b)*b5 + ((-16.109829710923627372335737204785697-10.739886473949084914890491469857131*(16.500000000000000000000000000000000 + (48.750000000000000000000000000000000 + (55. + 21.250000000000000000000000000000000*y)*y)*y)*y + 2.6849716184872712287226228674642828*(16.500000000000000000000000000000000 + (97.500000000000000000000000000000000 + (165. + 85.*y)*y)*y)*b)*b3 + ((-16.109829710923627372335737204785697-10.739886473949084914890491469857131*(16.500000000000000000000000000000000 + (48.750000000000000000000000000000000 + (55. + 21.250000000000000000000000000000000*y)*y)*y)*y + (-88.604063410079950547846554626321332-5.3699432369745424574452457349285655*(97.500000000000000000000000000000000 + (165. + 85.*y)*y)*y + .67124290462181780718065571686607068*(97.500000000000000000000000000000000 + (330. + 255.*y)*y)*b)*b)*b2 + ((37.589602658821797202116720144499959 + 25.059735105881198134744480096333306*(16.500000000000000000000000000000000 + (48.750000000000000000000000000000000 + (55. + 21.250000000000000000000000000000000*y)*y)*y)*y-12.529867552940599067372240048166653*(16.500000000000000000000000000000000 + (97.500000000000000000000000000000000 + (165. + 85.*y)*y)*y)*b)*b + (37.589602658821797202116720144499959 + 25.059735105881198134744480096333306*(16.500000000000000000000000000000000 + (48.750000000000000000000000000000000 + (55. + 21.250000000000000000000000000000000*y)*y)*y)*y + (103.37140731175994230582098039737489 + 6.2649337764702995336861200240833266*(97.500000000000000000000000000000000 + (165. + 85.*y)*y)*y-1.5662334441175748834215300060208316*(97.500000000000000000000000000000000 + (330. + 255.*y)*y)*b)*b + (186.06853316116789615047776471527480 + 11.276880797646539160635016043349988*(97.500000000000000000000000000000000 + (165. + 85.*y)*y)*y + .93974006647054493005291800361249898*(97.500000000000000000000000000000000 + (330. + 255.*y)*y)*a)*a)*a)*a)*a)*a;
    }

    inline double leg_f63_dxy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      return 2.6849716184872712287226228674642828*b4*(1.5000000000000000000000000000000000 + (16.500000000000000000000000000000000 + (48.750000000000000000000000000000000 + (55. + 21.250000000000000000000000000000000*y)*y)*y)*y) + ((-16.109829710923627372335737204785697-10.739886473949084914890491469857131*(16.500000000000000000000000000000000 + (48.750000000000000000000000000000000 + (55. + 21.250000000000000000000000000000000*y)*y)*y)*y + 2.6849716184872712287226228674642828*(16.500000000000000000000000000000000 + (97.500000000000000000000000000000000 + (165. + 85.*y)*y)*y)*b)*b3 + (-37.589602658821797202116720144499958*b2*(1.5000000000000000000000000000000000 + (16.500000000000000000000000000000000 + (48.750000000000000000000000000000000 + (55. + 21.250000000000000000000000000000000*y)*y)*y)*y) + ((37.589602658821797202116720144499959 + 25.059735105881198134744480096333306*(16.500000000000000000000000000000000 + (48.750000000000000000000000000000000 + (55. + 21.250000000000000000000000000000000*y)*y)*y)*y-12.529867552940599067372240048166653*(16.500000000000000000000000000000000 + (97.500000000000000000000000000000000 + (165. + 85.*y)*y)*y)*b)*b + (84.576605982349043704762620325124910 + 56.384403988232695803175080216749940*(16.500000000000000000000000000000000 + (48.750000000000000000000000000000000 + (55. + 21.250000000000000000000000000000000*y)*y)*y)*y + 11.276880797646539160635016043349988*(16.500000000000000000000000000000000 + (97.500000000000000000000000000000000 + (165. + 85.*y)*y)*y)*a)*a)*a)*a)*a;
    }

    // number 64
    inline double leg_f64(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b4 = b2*b2;
      double b6 = b4*b2;
      return .14753991249720048155874466028543705*(-.79672179899887262969191001703480969*b6 + (5.5770525929921084078433701192436678*b4 + (-10.038694667385795134118066214638602*b2 + 5.2583638733925593559666061124297439*a2)*a2)*a2)*a*(4.3333333333333333333333333333333333 + (19. + (26. + 11.333333333333333333333333333333333*y)*y)*y);
    }

    inline double leg_f64_dx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      double b6 = b5*b;
      return .14753991249720048155874466028543705*(-1.5934435979977452593838200340696194*b6 + (33.462315557952650447060220715462007*b4 + (-100.38694667385795134118066214638602*b2 + 73.617094227495830983532485574016415*a2)*a2)*a2)*(4.3333333333333333333333333333333333 + (19. + (26. + 11.333333333333333333333333333333333*y)*y)*y);
    }

    inline double leg_f64_dy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      double b6 = b5*b;
      return -.11754826450890581801454183188941421*b6*(4.3333333333333333333333333333333333 + (19. + (26. + 11.333333333333333333333333333333333*y)*y)*y) + ((3.0562548772315512683780876291247695 + .70528958705343490808725099133648527*(19. + (26. + 11.333333333333333333333333333333333*y)*y)*y-.11754826450890581801454183188941421*(19. + (52. + 34.*y)*y)*b)*b5 + (2.4685135546870221783053784696776984*b4*(4.3333333333333333333333333333333333 + (19. + (26. + 11.333333333333333333333333333333333*y)*y)*y) + ((-14.262522760413905919097742269248924-3.2913514062493629044071712929035978*(19. + (26. + 11.333333333333333333333333333333333*y)*y)*y + .82283785156234072610179282322589947*(19. + (52. + 34.*y)*y)*b)*b3 + (-7.4055406640610665349161354090330953*b2*(4.3333333333333333333333333333333333 + (19. + (26. + 11.333333333333333333333333333333333*y)*y)*y) + ((12.836270484372515327187968042324031 + 2.9622162656244266139664541636132381*(19. + (26. + 11.333333333333333333333333333333333*y)*y)*y-1.4811081328122133069832270818066191*(19. + (52. + 34.*y)*y)*b)*b + (23.533162554682944766511274744260725 + 5.4307298203114487922718326332909365*(19. + (26. + 11.333333333333333333333333333333333*y)*y)*y + .77581854575877839889597609047013379*(19. + (52. + 34.*y)*y)*a)*a)*a)*a)*a)*a)*a;
    }

    inline double leg_f64_dxx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      return .14753991249720048155874466028543705*(133.84926223181060178824088286184803*b4 + (-803.09557339086361072944529717108816*b2 + 883.40513072994997180238982688819698*a2)*a2)*a*(4.3333333333333333333333333333333333 + (19. + (26. + 11.333333333333333333333333333333333*y)*y)*y);
    }

    inline double leg_f64_dyy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      return (6.1125097544631025367561752582495390 + 1.4105791741068698161745019826729705*(19. + (26. + 11.333333333333333333333333333333333*y)*y)*y-.23509652901781163602908366377882842*(19. + (52. + 34.*y)*y)*b)*b5 + ((6.1125097544631025367561752582495390 + 1.4105791741068698161745019826729705*(19. + (26. + 11.333333333333333333333333333333333*y)*y)*y + (26.801004308030526507315537670786440 + 1.4105791741068698161745019826729705*(52. + 34.*y)*y-.11754826450890581801454183188941421*(68.*y + 52.)*b)*b)*b4 + ((-85.575136562483435514586453615493545-19.748108437496177426443027757421587*(19. + (26. + 11.333333333333333333333333333333333*y)*y)*y + 4.9370271093740443566107569393553969*(19. + (52. + 34.*y)*y)*b)*b3 + ((-85.575136562483435514586453615493545-19.748108437496177426443027757421587*(19. + (26. + 11.333333333333333333333333333333333*y)*y)*y + (-125.07135343747579036747250913033672-6.5827028124987258088143425858071959*(52. + 34.*y)*y + .82283785156234072610179282322589947*(68.*y + 52.)*b)*b)*b2 + ((128.36270484372515327187968042324031 + 29.622162656244266139664541636132381*(19. + (26. + 11.333333333333333333333333333333333*y)*y)*y-14.811081328122133069832270818066191*(19. + (52. + 34.*y)*y)*b)*b + (128.36270484372515327187968042324031 + 29.622162656244266139664541636132381*(19. + (26. + 11.333333333333333333333333333333333*y)*y)*y + (112.56421809372821133072525821730305 + 5.9244325312488532279329083272264761*(52. + 34.*y)*y-1.4811081328122133069832270818066191*(68.*y + 52.)*b)*b + (206.36773317183505410632964006505559 + 10.861459640622897584543665266581873*(52. + 34.*y)*y + .77581854575877839889597609047013379*(68.*y + 52.)*a)*a)*a)*a)*a)*a)*a;
    }

    inline double leg_f64_dxy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      return (6.1125097544631025367561752582495390 + 1.4105791741068698161745019826729705*(19. + (26. + 11.333333333333333333333333333333333*y)*y)*y-.23509652901781163602908366377882842*(19. + (52. + 34.*y)*y)*b)*b5 + (9.8740542187480887132215138787107937*b4*(4.3333333333333333333333333333333333 + (19. + (26. + 11.333333333333333333333333333333333*y)*y)*y) + ((-85.575136562483435514586453615493545-19.748108437496177426443027757421587*(19. + (26. + 11.333333333333333333333333333333333*y)*y)*y + 4.9370271093740443566107569393553969*(19. + (52. + 34.*y)*y)*b)*b3 + (-59.244325312488532279329083272264761*b2*(4.3333333333333333333333333333333333 + (19. + (26. + 11.333333333333333333333333333333333*y)*y)*y) + ((128.36270484372515327187968042324031 + 29.622162656244266139664541636132381*(19. + (26. + 11.333333333333333333333333333333333*y)*y)*y-14.811081328122133069832270818066191*(19. + (52. + 34.*y)*y)*b)*b + (282.39795065619533719813529693112870 + 65.168757843737385507261991599491238*(19. + (26. + 11.333333333333333333333333333333333*y)*y)*y + 10.861459640622897584543665266581873*(19. + (52. + 34.*y)*y)*a)*a)*a)*a)*a)*a;
    }

    // number 65
    inline double leg_f65(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b4 = b2*b2;
      double b6 = b4*b2;
      double b8 = b6*b2;
      return .69728750379731214109841784629056465e-1*(.10697706201272775653456441070328167*b8 + (-2.9953577363563771829678034996918866*b6 + (13.479109813603697323355115748613490*b4 + (-19.769361059952089407587503097966452*b2 + 9.1786319206920415106656264383415669*a2)*a2)*a2)*a2)*(3.2500000000000000000000000000000000 + (7.5000000000000000000000000000000000 + 4.2500000000000000000000000000000000*y)*y);
    }

    inline double leg_f65_dx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      double b6 = b5*b;
      double b7 = b6*b;
      return .69728750379731214109841784629056465e-1*(-11.981430945425508731871213998767547*b6 + (107.83287850882957858684092598890792*b4 + (-237.23233271942507289105003717559742*b2 + 146.85811073107266417065002301346507*a2)*a2)*a2)*a*(3.2500000000000000000000000000000000 + (7.5000000000000000000000000000000000 + 4.2500000000000000000000000000000000*y)*y);
    }

    inline double leg_f65_dy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      double b6 = b5*b;
      double b7 = b6*b;
      return (-.19394379818950552502986780397515149-.59675014827540161547651631992354305e-1*(7.5000000000000000000000000000000000 + 4.2500000000000000000000000000000000*y)*y + .74593768534425201934564539990442882e-2*(8.5000000000000000000000000000000000*y + 7.5000000000000000000000000000000000)*b)*b7 + (-.41772510379278113083356142394648014*b6*(3.2500000000000000000000000000000000 + (7.5000000000000000000000000000000000 + 4.2500000000000000000000000000000000*y)*y) + ((4.0728197619796160256272238834781813 + 1.2531753113783433925006842718394404*(7.5000000000000000000000000000000000 + 4.2500000000000000000000000000000000*y)*y-.20886255189639056541678071197324007*(8.5000000000000000000000000000000000*y + 7.5000000000000000000000000000000000)*b)*b5 + (3.7595259341350301775020528155183212*b4*(3.2500000000000000000000000000000000 + (7.5000000000000000000000000000000000 + 4.2500000000000000000000000000000000*y)*y) + ((-12.218459285938848076881671650434544-3.7595259341350301775020528155183212*(7.5000000000000000000000000000000000 + 4.2500000000000000000000000000000000*y)*y + .93988148353375754437551320387958031*(8.5000000000000000000000000000000000*y + 7.5000000000000000000000000000000000)*b)*b3 + (-8.2709570550970663905045161941403067*b2*(3.2500000000000000000000000000000000 + (7.5000000000000000000000000000000000 + 4.2500000000000000000000000000000000*y)*y) + ((8.9602034763551552563798925436519989 + 2.7569856850323554635015053980467689*(7.5000000000000000000000000000000000 + 4.2500000000000000000000000000000000*y)*y-1.3784928425161777317507526990233844*(8.5000000000000000000000000000000000*y + 7.5000000000000000000000000000000000)*b)*b + (16.640377884659574047562657581067998 + 5.1201162722029458607885100249439995*(7.5000000000000000000000000000000000 + 4.2500000000000000000000000000000000*y)*y + .64001453402536823259856375311799993*(8.5000000000000000000000000000000000*y + 7.5000000000000000000000000000000000)*a)*a)*a)*a)*a)*a)*a)*a;
    }

    inline double leg_f65_dxx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      double b6 = b5*b;
      double b7 = b6*b;
      return .69728750379731214109841784629056465e-1*(-23.962861890851017463742427997535093*b6 + (646.99727105297747152104555593344752*b4 + (-2372.3233271942507289105003717559742*b2 + 2056.0135502350172983891003221885110*a2)*a2)*a2)*(3.2500000000000000000000000000000000 + (7.5000000000000000000000000000000000 + 4.2500000000000000000000000000000000*y)*y);
    }

    inline double leg_f65_dyy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      double b6 = b5*b;
      double b7 = b6*b;
      return (-1.0144752520681827463100777438700232*y-.89512522241310242321477447988531460 + .63404703254261421644379858991876449e-1*b)*b7 + ((16.291279047918464102508895533912725 + 5.0127012455133735700027370873577615*(7.5000000000000000000000000000000000 + 4.2500000000000000000000000000000000*y)*y-.83545020758556226166712284789296028*(8.5000000000000000000000000000000000*y + 7.5000000000000000000000000000000000)*b)*b5 + ((16.291279047918464102508895533912725 + 5.0127012455133735700027370873577615*(7.5000000000000000000000000000000000 + 4.2500000000000000000000000000000000*y)*y + (21.303980293431837672511632621270488*y + 18.797629670675150887510264077591607-1.7753316911193198060426360517725406*b)*b)*b4 + ((-97.747674287510784615053373203476350-30.076207473080241420016422524146569*(7.5000000000000000000000000000000000 + 4.2500000000000000000000000000000000*y)*y + 7.5190518682700603550041056310366424*(8.5000000000000000000000000000000000*y + 7.5000000000000000000000000000000000)*b)*b3 + ((-97.747674287510784615053373203476350-30.076207473080241420016422524146569*(7.5000000000000000000000000000000000 + 4.2500000000000000000000000000000000*y)*y + (-63.911940880295513017534897863811461*y-56.392889012025452662530792232774819 + 7.9889926100369391271918622329764326*b)*b)*b2 + ((107.52244171626186307655871052382399 + 33.083828220388265562018064776561226*(7.5000000000000000000000000000000000 + 4.2500000000000000000000000000000000*y)*y-16.541914110194132781009032388280613*(8.5000000000000000000000000000000000*y + 7.5000000000000000000000000000000000)*b)*b + (107.52244171626186307655871052382399 + 33.083828220388265562018064776561226*(7.5000000000000000000000000000000000 + 4.2500000000000000000000000000000000*y)*y + (46.868756645550042879525591766795072*y + 41.354785275485331952522580970701534-11.717189161387510719881397941698768*b)*b + (87.041976627450079633404670424047990*y + 76.801744083044187911827650374159992 + 5.4401235392156299770877919015029993*a)*a)*a)*a)*a)*a)*a)*a;
    }

    inline double leg_f65_dxy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      double b6 = b5*b;
      double b7 = b6*b;
      return -.83545020758556226166712284789296028*b6*(3.2500000000000000000000000000000000 + (7.5000000000000000000000000000000000 + 4.2500000000000000000000000000000000*y)*y) + ((16.291279047918464102508895533912725 + 5.0127012455133735700027370873577615*(7.5000000000000000000000000000000000 + 4.2500000000000000000000000000000000*y)*y-.83545020758556226166712284789296028*(8.5000000000000000000000000000000000*y + 7.5000000000000000000000000000000000)*b)*b5 + (22.557155604810181065012316893109926*b4*(3.2500000000000000000000000000000000 + (7.5000000000000000000000000000000000 + 4.2500000000000000000000000000000000*y)*y) + ((-97.747674287510784615053373203476350-30.076207473080241420016422524146569*(7.5000000000000000000000000000000000 + 4.2500000000000000000000000000000000*y)*y + 7.5190518682700603550041056310366424*(8.5000000000000000000000000000000000*y + 7.5000000000000000000000000000000000)*b)*b3 + (-82.709570550970663905045161941403067*b2*(3.2500000000000000000000000000000000 + (7.5000000000000000000000000000000000 + 4.2500000000000000000000000000000000*y)*y) + ((107.52244171626186307655871052382399 + 33.083828220388265562018064776561226*(7.5000000000000000000000000000000000 + 4.2500000000000000000000000000000000*y)*y-16.541914110194132781009032388280613*(8.5000000000000000000000000000000000*y + 7.5000000000000000000000000000000000)*b)*b + (232.96529038523403666587720613495197 + 71.681627810841242051039140349215992*(7.5000000000000000000000000000000000 + 4.2500000000000000000000000000000000*y)*y + 10.240232544405891721577020049887999*(8.5000000000000000000000000000000000*y + 7.5000000000000000000000000000000000)*a)*a)*a)*a)*a)*a)*a;
    }

    // number 66
    inline double leg_f66(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b4 = b2*b2;
      double b6 = b4*b2;
      double b8 = b6*b2;
      return .32612178450642258898018265094570013e-1*(.79720045437338092375232558872693519*b8 + (-9.5664054524805710850279070647232223*b6 + (31.569137993185884580592093313586634*b4 + (-39.085599420134904718828306007297737*b2 + 16.285666425056210299511794169707390*a2)*a2)*a2)*a2)*a*(1. + y);
    }

    inline double leg_f66_dx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      double b6 = b5*b;
      double b7 = b6*b;
      double b8 = b7*b;
      return .32612178450642258898018265094570013e-1*(1.5944009087467618475046511774538704*b8 + (-57.398432714883426510167442388339334*b6 + (315.69137993185884580592093313586634*b4 + (-547.19839188188866606359628410216831*b2 + 293.14199565101178539121229505473303*a2)*a2)*a2)*a2)*(1. + y);
    }

    inline double leg_f66_dy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      double b6 = b5*b;
      double b7 = b6*b;
      double b8 = b7*b;
      return .25998443478957790700562193556377897e-1*b8*(1. + y) + ((-.20798754783166232560449754845102318-.20798754783166232560449754845102318*y + .25998443478957790700562193556377897e-1*b)*b7 + (-.93594396524248046522023896802960429*b6*(1. + y) + ((1.8718879304849609304404779360592086 + 1.8718879304849609304404779360592086*y-.31198132174749348840674632267653476*b)*b5 + (5.1476918088336425587113143241628238*b4*(1. + y) + ((-4.1181534470669140469690514593302588-4.1181534470669140469690514593302588*y + 1.0295383617667285117422628648325647*b)*b3 + (-8.9226658019783137684329448285488942*b2*(1. + y) + ((2.5493330862795182195522699510139698 + 2.5493330862795182195522699510139698*y-1.2746665431397591097761349755069849*b)*b + (4.7799995367740966616605061581511934 + 4.7799995367740966616605061581511934*y + .53111105964156629574005623979457703*a)*a)*a)*a)*a)*a)*a)*a)*a;
    }

    inline double leg_f66_dxx(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      double b6 = b5*b;
      double b7 = b6*b;
      return .32612178450642258898018265094570013e-1*(-229.59373085953370604066976955335733*b6 + (2525.5310394548707664473674650869307*b4 + (-6566.3807025826639927631554092260198*b2 + 4690.2719304161885662593967208757284*a2)*a2)*a2)*a*(1. + y);
    }

    inline double leg_f66_dyy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      double b6 = b5*b;
      double b7 = b6*b;
      return (-.41597509566332465120899509690204634-.41597509566332465120899509690204634*y + .51996886957915581401124387112755795e-1*b)*b7 + ((-.41597509566332465120899509690204634-.41597509566332465120899509690204634*y-.41597509566332465120899509690204634*b)*b6 + ((11.231327582909765582642867616355252 + 11.231327582909765582642867616355252*y-1.8718879304849609304404779360592086*b)*b5 + ((11.231327582909765582642867616355252 + 11.231327582909765582642867616355252*y + 3.7437758609699218608809558721184173*b)*b4 + ((-41.181534470669140469690514593302588-41.181534470669140469690514593302588*y + 10.295383617667285117422628648325647*b)*b3 + ((-41.181534470669140469690514593302588-41.181534470669140469690514593302588*y-8.2363068941338280939381029186605177*b)*b2 + ((35.690663207913255073731779314195579 + 35.690663207913255073731779314195579*y-17.845331603956627536865889657097789*b)*b + (35.690663207913255073731779314195579 + 35.690663207913255073731779314195579*y + 5.0986661725590364391045399020279395*b + 9.5599990735481933233210123163023868*a)*a)*a)*a)*a)*a)*a)*a;
    }

    inline double leg_f66_dxy(double x, double y)
    {
      double a = 1.0 + 2.0*x + y;
      double b = 1.0-y;
      double a2 = a*a;
      double b2 = b*b;
      double b3 = b2*b;
      double b4 = b3*b;
      double b5 = b4*b;
      double b6 = b5*b;
      double b7 = b6*b;
      return (-.41597509566332465120899509690204634-.41597509566332465120899509690204634*y + .51996886957915581401124387112755795e-1*b)*b7+
          (-3.7437758609699218608809558721184173*b6*(1. + y) + ((11.231327582909765582642867616355252 + 11.231327582909765582642867616355252*y-
          1.8718879304849609304404779360592086*b)*b5 + (41.181534470669140469690514593302588*b4*(1. + y) + ((-41.181534470669140469690514593302588-
          41.181534470669140469690514593302588*y + 10.295383617667285117422628648325647*b)*b3 + (-107.07198962373976522119533794258673*b2*(1. + y)+
          ((35.690663207913255073731779314195579 + 35.690663207913255073731779314195579*y-17.845331603956627536865889657097789*b)*b+
          (76.479992588385546586568098530419093 + 76.479992588385546586568098530419093*y+
          9.5599990735481933233210123163023868*a)*a)*a)*a)*a)*a)*a)*a;
    }

    static Shapeset::shape_fn_t leg_tri_fn[] =
    {
      leg_tri_l0_l0,   leg_f1,    leg_f2,    leg_f3,    leg_f4,    leg_f5,    leg_f6,    leg_f7_0,
        leg_f8_0,    leg_f9_0,    leg_f10,   leg_f11,
      leg_f12,   leg_f13,   leg_f14,   leg_f15,   leg_f16_0,  leg_f17_0,
       leg_f18_0,  leg_f19,   leg_f20,   leg_f21,   leg_f22,
      leg_f23,   leg_f24,   leg_f25,   leg_f26,   leg_f27,   leg_f28,   leg_f29_0,
       leg_f30_0,  leg_f31_0,  leg_f32,   leg_f33,
      leg_f34,   leg_f35,   leg_f36,   leg_f37,   leg_f38,   leg_f39,   leg_f40,
      leg_f41,   leg_f42,   leg_f43,   leg_f44,   leg_f45,   leg_f46_0, 
      leg_f47_0,  leg_f48_0,  leg_f49,   leg_f50,   leg_f51,
      leg_f52,   leg_f53,   leg_f54,   leg_f55,   leg_f56,   leg_f57,   leg_f58,
      leg_f59,   leg_f60,   leg_f61,   leg_f62,   leg_f63,   leg_f64,   leg_f65,
      leg_f66
    };
    static Shapeset::shape_fn_t leg_tri_fn_dx[] =
    {
      leg_tri_l0_l0x,   leg_f1_dx,   leg_f2_dx,   leg_f3_dx,   leg_f4_dx,   leg_f5_dx,   leg_f6_dx,   leg_f7_dx_0,
       leg_f8_dx_0,  leg_f9_dx_0,  leg_f10_dx,  leg_f11_dx,
      leg_f12_dx,  leg_f13_dx,  leg_f14_dx,  leg_f15_dx,  leg_f16_dx_0,  leg_f17_dx_0,
       leg_f18_dx_0,  leg_f19_dx,  leg_f20_dx,  leg_f21_dx,  leg_f22_dx,
      leg_f23_dx,  leg_f24_dx,  leg_f25_dx,  leg_f26_dx,  leg_f27_dx,  leg_f28_dx,  leg_f29_dx_0,
       leg_f30_dx_0,  leg_f31_dx_0,  leg_f32_dx,  leg_f33_dx,
      leg_f34_dx,  leg_f35_dx,  leg_f36_dx,  leg_f37_dx,  leg_f38_dx,  leg_f39_dx,  leg_f40_dx,
      leg_f41_dx,  leg_f42_dx,  leg_f43_dx,  leg_f44_dx,  leg_f45_dx,  leg_f46_dx_0, 
      leg_f47_dx_0,  leg_f48_dx_0,  leg_f49_dx,  leg_f50_dx,  leg_f51_dx,
      leg_f52_dx,  leg_f53_dx,  leg_f54_dx,  leg_f55_dx,  leg_f56_dx,  leg_f57_dx,  leg_f58_dx,
      leg_f59_dx,  leg_f60_dx,  leg_f61_dx,  leg_f62_dx,  leg_f63_dx,  leg_f64_dx,  leg_f65_dx,
      leg_f66_dx
    };
    static Shapeset::shape_fn_t leg_tri_fn_dy[] =
    {
      leg_tri_l0_l0y,   leg_f1_dy,   leg_f2_dy,   leg_f3_dy,   leg_f4_dy,   leg_f5_dy,   leg_f6_dy,   leg_f7_dy_0,
       leg_f8_dy_0,  leg_f9_dy_0,  leg_f10_dy,  leg_f11_dy,
      leg_f12_dy,  leg_f13_dy,  leg_f14_dy,  leg_f15_dy,  leg_f16_dy_0,  leg_f17_dy_0,
       leg_f18_dy_0,  leg_f19_dy,  leg_f20_dy,  leg_f21_dy,  leg_f22_dy,
      leg_f23_dy,  leg_f24_dy,  leg_f25_dy,  leg_f26_dy,  leg_f27_dy,  leg_f28_dy,  leg_f29_dy_0,
       leg_f30_dy_0,  leg_f31_dy_0,  leg_f32_dy,  leg_f33_dy,
      leg_f34_dy,  leg_f35_dy,  leg_f36_dy,  leg_f37_dy,  leg_f38_dy,  leg_f39_dy,  leg_f40_dy,
      leg_f41_dy,  leg_f42_dy,  leg_f43_dy,  leg_f44_dy,  leg_f45_dy,  leg_f46_dy_0, 
      leg_f47_dy_0,  leg_f48_dy_0,  leg_f49_dy,  leg_f50_dy,  leg_f51_dy,
      leg_f52_dy,  leg_f53_dy,  leg_f54_dy,  leg_f55_dy,  leg_f56_dy,  leg_f57_dy,  leg_f58_dy,
      leg_f59_dy,  leg_f60_dy,  leg_f61_dy,  leg_f62_dy,  leg_f63_dy,  leg_f64_dy,  leg_f65_dy,
      leg_f66_dy
    };
    static Shapeset::shape_fn_t leg_tri_fn_dxx[] =
    {
      leg_tri_l0_l0xx,   leg_f1_dxx,   leg_f2_dxx,   leg_f3_dxx,   leg_f4_dxx,   leg_f5_dxx,   leg_f6_dxx,   leg_f7_dxx_0,
      leg_f7_dxx_1, leg_f8_dxx_0, leg_f8_dxx_1, leg_f9_dxx_0, leg_f9_dxx_1, leg_f10_dxx,  leg_f11_dxx,
      leg_f12_dxx,  leg_f13_dxx,  leg_f14_dxx,  leg_f15_dxx,  leg_f16_dxx_0, leg_f16_dxx_1, leg_f17_dxx_0,
      leg_f17_dxx_1, leg_f18_dxx_0, leg_f18_dxx_1, leg_f19_dxx,  leg_f20_dxx,  leg_f21_dxx,  leg_f22_dxx,
      leg_f23_dxx,  leg_f24_dxx,  leg_f25_dxx,  leg_f26_dxx,  leg_f27_dxx,  leg_f28_dxx,  leg_f29_dxx_0,
      leg_f29_dxx_1, leg_f30_dxx_0, leg_f30_dxx_1, leg_f31_dxx_0, leg_f31_dxx_1, leg_f32_dxx,  leg_f33_dxx,
      leg_f34_dxx,  leg_f35_dxx,  leg_f36_dxx,  leg_f37_dxx,  leg_f38_dxx,  leg_f39_dxx,  leg_f40_dxx,
      leg_f41_dxx,  leg_f42_dxx,  leg_f43_dxx,  leg_f44_dxx,  leg_f45_dxx,  leg_f46_dxx_0, leg_f46_dxx_1,
      leg_f47_dxx_0, leg_f47_dxx_1, leg_f48_dxx_0, leg_f48_dxx_1, leg_f49_dxx,  leg_f50_dxx,  leg_f51_dxx,
      leg_f52_dxx,  leg_f53_dxx,  leg_f54_dxx,  leg_f55_dxx,  leg_f56_dxx,  leg_f57_dxx,  leg_f58_dxx,
      leg_f59_dxx,  leg_f60_dxx,  leg_f61_dxx,  leg_f62_dxx,  leg_f63_dxx,  leg_f64_dxx,  leg_f65_dxx,
      leg_f66_dxx
    };
    static Shapeset::shape_fn_t leg_tri_fn_dxy[] =
    {
      leg_tri_l0_l0xy,   leg_f1_dxy,   leg_f2_dxy,   leg_f3_dxy,   leg_f4_dxy,   leg_f5_dxy,   leg_f6_dxy,   leg_f7_dxy_0,
      leg_f7_dxy_1, leg_f8_dxy_0, leg_f8_dxy_1, leg_f9_dxy_0, leg_f9_dxy_1, leg_f10_dxy,  leg_f11_dxy,
      leg_f12_dxy,  leg_f13_dxy,  leg_f14_dxy,  leg_f15_dxy,  leg_f16_dxy_0, leg_f16_dxy_1, leg_f17_dxy_0,
      leg_f17_dxy_1, leg_f18_dxy_0, leg_f18_dxy_1, leg_f19_dxy,  leg_f20_dxy,  leg_f21_dxy,  leg_f22_dxy,
      leg_f23_dxy,  leg_f24_dxy,  leg_f25_dxy,  leg_f26_dxy,  leg_f27_dxy,  leg_f28_dxy,  leg_f29_dxy_0,
      leg_f29_dxy_1, leg_f30_dxy_0, leg_f30_dxy_1, leg_f31_dxy_0, leg_f31_dxy_1, leg_f32_dxy,  leg_f33_dxy,
      leg_f34_dxy,  leg_f35_dxy,  leg_f36_dxy,  leg_f37_dxy,  leg_f38_dxy,  leg_f39_dxy,  leg_f40_dxy,
      leg_f41_dxy,  leg_f42_dxy,  leg_f43_dxy,  leg_f44_dxy,  leg_f45_dxy,  leg_f46_dxy_0, leg_f46_dxy_1,
      leg_f47_dxy_0, leg_f47_dxy_1, leg_f48_dxy_0, leg_f48_dxy_1, leg_f49_dxy,  leg_f50_dxy,  leg_f51_dxy,
      leg_f52_dxy,  leg_f53_dxy,  leg_f54_dxy,  leg_f55_dxy,  leg_f56_dxy,  leg_f57_dxy,  leg_f58_dxy,
      leg_f59_dxy,  leg_f60_dxy,  leg_f61_dxy,  leg_f62_dxy,  leg_f63_dxy,  leg_f64_dxy,  leg_f65_dxy,
      leg_f66_dxy
    };
    static Shapeset::shape_fn_t leg_tri_fn_dyy[] =
    {
      leg_tri_l0_l0yy,   leg_f1_dyy,   leg_f2_dyy,   leg_f3_dyy,   leg_f4_dyy,   leg_f5_dyy,   leg_f6_dyy,   leg_f7_dyy_0,
      leg_f7_dyy_1, leg_f8_dyy_0, leg_f8_dyy_1, leg_f9_dyy_0, leg_f9_dyy_1, leg_f10_dyy,  leg_f11_dyy,
      leg_f12_dyy,  leg_f13_dyy,  leg_f14_dyy,  leg_f15_dyy,  leg_f16_dyy_0, leg_f16_dyy_1, leg_f17_dyy_0,
      leg_f17_dyy_1, leg_f18_dyy_0, leg_f18_dyy_1, leg_f19_dyy,  leg_f20_dyy,  leg_f21_dyy,  leg_f22_dyy,
      leg_f23_dyy,  leg_f24_dyy,  leg_f25_dyy,  leg_f26_dyy,  leg_f27_dyy,  leg_f28_dyy,  leg_f29_dyy_0,
      leg_f29_dyy_1, leg_f30_dyy_0, leg_f30_dyy_1, leg_f31_dyy_0, leg_f31_dyy_1, leg_f32_dyy,  leg_f33_dyy,
      leg_f34_dyy,  leg_f35_dyy,  leg_f36_dyy,  leg_f37_dyy,  leg_f38_dyy,  leg_f39_dyy,  leg_f40_dyy,
      leg_f41_dyy,  leg_f42_dyy,  leg_f43_dyy,  leg_f44_dyy,  leg_f45_dyy,  leg_f46_dyy_0, leg_f46_dyy_1,
      leg_f47_dyy_0, leg_f47_dyy_1, leg_f48_dyy_0, leg_f48_dyy_1, leg_f49_dyy,  leg_f50_dyy,  leg_f51_dyy,
      leg_f52_dyy,  leg_f53_dyy,  leg_f54_dyy,  leg_f55_dyy,  leg_f56_dyy,  leg_f57_dyy,  leg_f58_dyy,
      leg_f59_dyy,  leg_f60_dyy,  leg_f61_dyy,  leg_f62_dyy,  leg_f63_dyy,  leg_f64_dyy,  leg_f65_dyy,
      leg_f66_dyy
    };
    Shapeset::shape_fn_t* leg_tri_shape_fn_table[1]     = { leg_tri_fn };
    Shapeset::shape_fn_t* leg_tri_shape_fn_table_dx[1]  = { leg_tri_fn_dx };
    Shapeset::shape_fn_t* leg_tri_shape_fn_table_dy[1]  = { leg_tri_fn_dy };
    Shapeset::shape_fn_t* leg_tri_shape_fn_table_dxx[1]  = { leg_tri_fn_dxx };
    Shapeset::shape_fn_t* leg_tri_shape_fn_table_dxy[1]  = { leg_tri_fn_dxy };
    Shapeset::shape_fn_t* leg_tri_shape_fn_table_dyy[1]  = { leg_tri_fn_dyy };

    static int qb_0[] = { 0, };
    static int qb_1[] = { 1, 2, 3, };
    static int qb_2[] = { 1, 2, 3, 4, 5, 6, };
    static int qb_3[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
    static int qb_4[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, };
    static int qb_5[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, };
    static int qb_6[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, };
    static int qb_7[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36 };
    static int qb_8[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, };
    static int qb_9[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55 };
    static int qb_10[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66 };

    int* leg_tri_bubble_indices[11] = {  qb_0,   qb_1,   qb_2,   qb_3,   qb_4,   qb_5,   qb_6,   qb_7,   qb_8,   qb_9,   qb_10 };

    int leg_tri_bubble_count[11] = { 1,  3,  6,  10,  15,  21,  28,  36,  45,  55,  66 };

    int leg_tri_vertex_indices[4] = { -1, -1, -1, -1 };

    static int leg_tri_edge_indices_0[22] =  { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1  };
    static int leg_tri_edge_indices_1[22] =  { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1  };
    static int leg_tri_edge_indices_2[22] =  { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1  };
    static int leg_tri_edge_indices_3[22] =  { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1  };

    int* leg_tri_edge_indices[4] =
    {
      leg_tri_edge_indices_0,
      leg_tri_edge_indices_1,
      leg_tri_edge_indices_2,
      leg_tri_edge_indices_3
    };

    int leg_tri_index_to_order[] =
    {
      0, 
      1, 1, 1, 
      2, 2, 2, 
      3, 3, 3, 3, 
      4, 4, 4, 4, 4,
      5, 5, 5, 5, 5, 5, 
      6, 6, 6, 6, 6, 6, 6, 
      7, 7, 7, 7, 7, 7, 7, 7, 
      8, 8, 8, 8, 8, 8, 8, 8, 8, 
      9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 
      10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    static Shapeset::shape_fn_t** leg_shape_fn_table[2] =
    {
      leg_tri_shape_fn_table,
      leg_quad_shape_fn_table
    };

    static Shapeset::shape_fn_t** leg_shape_fn_table_dx[2] =
    {
      leg_tri_shape_fn_table_dx,
      leg_quad_shape_fn_table_dx
    };

    static Shapeset::shape_fn_t** leg_shape_fn_table_dy[2] =
    {
      leg_tri_shape_fn_table_dy,
      leg_quad_shape_fn_table_dy
    };

    static Shapeset::shape_fn_t** leg_shape_fn_table_dxx[2] =
    {
      NULL,
      leg_quad_shape_fn_table_dxx
    };

    static Shapeset::shape_fn_t** leg_shape_fn_table_dyy[2] =
    {
      NULL,
      leg_quad_shape_fn_table_dyy
    };

    static Shapeset::shape_fn_t** leg_shape_fn_table_dxy[2] =
    {
      NULL,
      leg_quad_shape_fn_table_dxy
    };

    static int* leg_vertex_indices[2] =
    {
      leg_tri_vertex_indices,
      leg_quad_vertex_indices
    };

    static int** leg_edge_indices[2] =
    {
      leg_tri_edge_indices,
      leg_quad_edge_indices
    };

    static int** leg_bubble_indices[2] =
    {
      leg_tri_bubble_indices,
      leg_quad_bubble_indices
    };

    static int* leg_bubble_count[2] =
    {
      leg_tri_bubble_count,
      leg_quad_bubble_count
    };

    static int* leg_index_to_order[2] =
    {
      leg_tri_index_to_order,
      leg_quad_index_to_order
    };

    int L2ShapesetLegendre::get_max_index(ElementMode2D mode) { return max_index[mode]; }

    L2ShapesetLegendre::L2ShapesetLegendre()
    {
      shape_table[0] = leg_shape_fn_table;
      shape_table[1] = leg_shape_fn_table_dx;
      shape_table[2] = leg_shape_fn_table_dy;
      shape_table[3] = leg_shape_fn_table_dxx;
      shape_table[4] = leg_shape_fn_table_dyy;
      shape_table[5] = leg_shape_fn_table_dxy;

      vertex_indices = leg_vertex_indices;
      edge_indices = leg_edge_indices;
      bubble_indices = leg_bubble_indices;
      bubble_count = leg_bubble_count;
      index_to_order = leg_index_to_order;

      ref_vert[0][0][0] = -1.0;
      ref_vert[0][0][1] = -1.0;
      ref_vert[0][1][0] =  1.0;
      ref_vert[0][1][1] = -1.0;
      ref_vert[0][2][0] = -1.0;
      ref_vert[0][2][1] =  1.0;

      ref_vert[1][0][0] = -1.0;
      ref_vert[1][0][1] = -1.0;
      ref_vert[1][1][0] =  1.0;
      ref_vert[1][1][1] = -1.0;
      ref_vert[1][2][0] =  1.0;
      ref_vert[1][2][1] =  1.0;
      ref_vert[1][3][0] = -1.0;
      ref_vert[1][3][1] =  1.0;

      max_order = 10;
      min_order = 0;
      num_components = 1;

      ebias = 0;

      comb_table = NULL;
    }
    
    const int L2ShapesetLegendre::max_index[2] = { 66, 120 };
  }
}