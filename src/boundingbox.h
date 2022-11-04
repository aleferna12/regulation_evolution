//
// Created by aleferna on 31-10-2022.
//

#ifndef REGULATION_EVOLUTION_BOUNDINGBOX_H
#define REGULATION_EVOLUTION_BOUNDINGBOX_H

#include <utility>

using namespace std;

class BoundingBox {
private:
    int minx;
    int miny;
    int maxx;
    int maxy;

public:
    BoundingBox() = default;

    BoundingBox(int minx, int miny, int maxx, int maxy) : minx(minx), miny(miny), maxx(maxx), maxy(maxy) {}

    int getMinX() const {
      return minx;
    }

    void setMinX(int val) {
      minx = val;
    }

    int getMinY() const {
      return miny;
    }

    void setMinY(int val) {
      miny = val;
    }

    int getMaxX() const {
      return maxx;
    }

    void setMaxX(int val) {
      maxx = val;
    }

    int getMaxY() const {
      return maxy;
    }

    void setMaxY(int val) {
      maxy = val;
    }

    int getSideX() const {
      return maxx - minx;
    }

    int getSideY() const {
      return maxy - miny;
    }

    int getArea() const {
      return getSideX() * getSideY();
    }

    pair<int, int> getOverlapLengths(const BoundingBox &box) const;

    bool overlaps(const BoundingBox &box) const;

    int getOverlapArea(const BoundingBox &box) const;
};

#endif //REGULATION_EVOLUTION_BOUNDINGBOX_H
