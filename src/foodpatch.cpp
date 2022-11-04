//
// Created by aleferna on 31-10-2022.
//

#include <algorithm>
#include "foodpatch.h"
#include "dish.h"

FoodPatch::FoodPatch(
  Dish *owner,
  int centerx,
  int centery,
  int side,
  int food_start
) : centerx(centerx), centery(centery), side(side), food_start(food_start), owner(owner) {
  // TODO: Add some random noise so patches look more interesting
  food_spots_left = side * side;
  sigma = new int[side * side];
  fill_n(sigma, side * side, food_start);
}

FoodPatch::~FoodPatch() {
  delete[] sigma;
}
