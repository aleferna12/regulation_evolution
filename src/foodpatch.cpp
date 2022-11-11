//
// Created by aleferna on 31-10-2022.
//

#include <algorithm>
#include "foodpatch.h"
#include "dish.h"

FoodPatch::FoodPatch(
        Dish *owner,
        int id,
        int x,
        int y,
        int length,
        int food_per_spot
) : id(id), x(x), y(y), length(length), food_per_spot(food_per_spot), owner(owner) {
  initSigma();
}

FoodPatch::FoodPatch(
        const FoodPatch &fp
) : id(fp.id),
    x(fp.x),
    y(fp.y),
    length(fp.length),
    food_per_spot(fp.food_per_spot),
    food_left(fp.food_left),
    owner(fp.owner),
    empty(fp.empty),
    removed(fp.removed) {
}

// Initialize the FoodPatch on both the FoodPlane of dish and its internal plane
void FoodPatch::initSigma() {
  for (int i = max(2, x); i < min(owner->SizeX() - 2, x + length); ++i) {
    // Even though the actual border is only 0 and size - 1, for some reason
    // AmoebaMove also prevents cells from reaching pos 1 and size - 2
    // which is why this loop is smaller
    for (int j = max(2, y); j < min(owner->SizeY() - 2, y + length); ++j) {
      if (owner->FoodPlane->Sigma(i, j) == -1) {
        owner->FoodPlane->setSigma(i, j, id);
        food_left += food_per_spot;
      }
    }
  }
  // Just to be sure
  if (food_left == 0) {
    cerr << "tried to initialize FoodPatch at an invalid position: " << x << ", " << y;
    exit(1);
  }
}

// Valid check to guarantee that there has been no miscount of food_left
void FoodPatch::checkEmpty() {
  empty = true;
  for (int i = max(2, x); i < min(owner->SizeX() - 2, x + length); ++i) {
    // Even though the actual border is only 0 and size - 1, for some reason
    // AmoebaMove also prevents cells from reaching pos 1 and size - 2
    // which is why this loop is smaller
    for (int j = max(2, y); j < min(owner->SizeY() - 2, y + length); ++j) {
      if (owner->FoodPlane->Sigma(i, j) == id) {
        cerr << "FoodPatch has no food left but is not empty";
        exit(1);
      }
    }
  }
}

int FoodPatch::consumeFood(int i, int j) {
  if (owner->FoodPlane->Sigma(i, j) == id) {
    cout << i << j << endl;
    owner->FoodPlane->setSigma(i, j, -1);
    food_left -= food_per_spot;
    if (food_left == 0)
      checkEmpty();
    return food_per_spot;
  }
  return 0;
}
