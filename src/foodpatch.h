//
// Created by aleferna on 31-10-2022.
//

#ifndef REGULATION_EVOLUTION_FOODPATCH_H
#define REGULATION_EVOLUTION_FOODPATCH_H

#include <vector>
#include <stdexcept>
#include "boundingbox.h"

class Dish;

class FoodPatch {
private:
    int centerx;
    int centery;
    int side;
    int food_start;
    int food_spots_left;
    int *sigma;
    Dish *owner;
    BoundingBox grad_box{};

public:
    FoodPatch(Dish *owner, int centerx, int centery, int side, int food_start);

    virtual ~FoodPatch();

    int getCenterX() const {
      return centerx;
    }

    void setCenterX(int val) {
      centerx = val;
    }

    int getCenterY() const {
      return centery;
    }

    void setCenterY(int val) {
      centery = val;
    }

    int getSide() const {
      return side;
    }

    void setSide(int val) {
      side = val;
    }

    int getFoodStart() const {
      return food_start;
    }

    int getFoodSpotsLeft() const {
      return food_spots_left;
    }

    void decrementFoodSpotsLeft() {
      food_spots_left -= 1;
    }

    int getSigma(int x, int y) const {
      return sigma[x * side + y];
    }

    void setSigma(int x, int y, int val) {
      sigma[x * side + y] = val;
    }

    // TODO
    BoundingBox &getGradBox() {
      throw logic_error("grad boxes are yet to be implemented");
      return grad_box;
    }
};


#endif //REGULATION_EVOLUTION_FOODPATCH_H
