//
// Created by aleferna on 31-10-2022.
//

#include <algorithm>
#include "foodpatch.h"
#include "dish.h"
#include "misc.h"

FoodPatch::FoodPatch(
    Dish *owner,
    int id,
    int x,
    int y,
    int length,
    int food_per_spot,
    int *sigmas
) : id(id), x(x), y(y), length(length), food_per_spot(food_per_spot), owner(owner) {
    sigma = new int[length * length]{};
    initSigmas(sigmas);
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
    sigma = new int[length * length];
    copy(fp.sigma, fp.sigma + fp.length * fp.length, sigma);
}

FoodPatch &FoodPatch::operator=(const FoodPatch &fp) {
    if (this == &fp)
        return *this;
    if (length != fp.length)
        throw runtime_error("tried assigning FoodPatches of different length");

    id = fp.id;
    x = fp.x;
    y = fp.y;
    length = fp.length;
    food_per_spot = fp.food_per_spot;
    food_left = fp.food_left;
    owner = fp.owner;
    empty = fp.empty;
    removed = fp.removed;

    copy(fp.sigma, fp.sigma + fp.length * fp.length, sigma);

    return *this;
}

FoodPatch::~FoodPatch() {
    delete[] sigma;
}

// Initialize the FoodPatch on both the FoodPlane of dish and its internal plane
void FoodPatch::initSigmas(const int *sigmas) {
    int minx = max(2, x);
    int miny = max(2, y);
    int maxx = min(owner->SizeX() - 2, x + length);
    int maxy = min(owner->SizeY() - 2, y + length);
    int center = length / 2;
    for (int gi = minx; gi < maxx; ++gi) {
        int i = getLocalX(gi);
        // Even though the actual border is only 0 and size - 1, for some reason
        // AmoebaMove also prevents cells from reaching pos 1 and size - 2
        // which is why this loop is smaller
        for (int gj = miny; gj < maxy; ++gj) {
            int j = getLocalY(gj);
            int value;
            if (sigmas == nullptr)
                value = food_per_spot;
            else
                value = sigmas[i * length + j];
            double dist_to_center = dist(center, center, i + .5, j + .5);
            if (owner->FoodPlane->Sigma(gi, gj) == -1 and value != 0 and dist_to_center <= length / 2.) {
                owner->FoodPlane->setSigma(gi, gj, id);
                setSigma(i, j, value);
                food_left += value;
            }
        }
    }
    // Just to be sure
    if (food_left == 0) {
        empty = true;
        cerr << "Warning: tried to initialize FoodPatch at an invalid position: " << x << ", " << y;
    }
}

// Valid check to guarantee that there has been no miscount of food_left
void FoodPatch::checkEmpty() {
    int minx = max(2, x);
    int miny = max(2, y);
    int maxx = min(owner->SizeX() - 2, x + length);
    int maxy = min(owner->SizeY() - 2, y + length);
    for (int gi = minx; gi < maxx; ++gi) {
        for (int gj = miny; gj < maxy; ++gj) {
            if (owner->FoodPlane->Sigma(gi, gj) == id) {
                throw runtime_error("FoodPatch has no food left but is not empty");
            }
        }
    }
    empty = true;
}

int FoodPatch::consumeFood(int gi, int gj) {
    if (owner->FoodPlane->Sigma(gi, gj) != id)
        throw runtime_error("consumeFood called with incorrect FoodPatch id");

    int i = getLocalX(gi);
    int j = getLocalY(gj);
    int sigma_at = getSigma(i, j);
    if (sigma_at > 0) {
        setSigma(i, j, sigma_at - 1);
        food_left -= 1;
        if (sigma_at == 1) {
            owner->FoodPlane->setSigma(gi, gj, -1);
            if (food_left == 0 and not empty)
                checkEmpty();
        }
        return 1;
    }
    return 0;
}
