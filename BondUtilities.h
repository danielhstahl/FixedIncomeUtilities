#ifndef __BONDUTILITIES_H_INCLUDED__
#define __BONDUTILITIES_H_INCLUDED__
#include <cmath>

auto convertLiborToContinuous(
  const auto&, //continuously compounded yield at point t.
  const auto& //point t
);
auto convertLiborToBond(
  const auto&, //libor rate over (0, t)
  const auto& //t
);
auto convertBondToLibor(
  const auto&, //bond from 0 to t
  const auto& //t
);
#include "BondUtilities.hpp"

#endif
