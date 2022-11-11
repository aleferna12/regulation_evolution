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
    int id;
    int x;
    int y;
    int length;
    int food_per_spot;
    int food_left = 0;
    Dish *owner;
    BoundingBox grad_box{};

public:
    FoodPatch(Dish *owner, int id, int x, int y, int length, int food_per_spot);

    FoodPatch(const FoodPatch &fp);

    int getX() const {
      return x;
    }

    int getY() const {
      return y;
    }

    double getCenterX() const {
      return x + length / (double) 2;
    }

    double getCenterY() const {
      return y + length / (double) 2;
    }

    int getId() const {
      return id;
    }

    int getLength() const {
      return length;
    }

    int getFoodPerSpot() const {
      return food_per_spot;
    }

    unsigned int getFoodLeft() const {
      return food_left;
    }

    // TODO
    BoundingBox &getGradBox() {
      throw logic_error("grad boxes are yet to be implemented");
      return grad_box;
    }

    void initSigma();

    int consumeFood(int gi, int gj);

    // Whether there is any food spots for this patch
    // If true we can remove FoodPatch from the ChemPlane
    // If false, removed is also false
    bool empty = false;

    // Whether FoodPatch has been removed from the ChemPlane
    // If true, we can recycle the position on fpatches
    // If true, empty is also true
    bool removed = false;

    void checkEmpty();
};


#endif //REGULATION_EVOLUTION_FOODPATCH_H