//
// Created by aleferna on 31-10-2022.
//

#include <algorithm>
#include "boundingbox.h"

pair<int, int> BoundingBox::getOverlapLengths(const BoundingBox &box) const {
    int lenx = max(minx, box.minx) - min(maxx, box.maxx);
    int leny = max(miny, box.miny) - min(maxy, box.maxy);
    return {lenx, leny};
}

int BoundingBox::getOverlapArea(const BoundingBox &box) const {
    auto lens = getOverlapLengths(box);
    return max(0, lens.first * lens.second);
}

bool BoundingBox::overlaps(const BoundingBox &box) const {
    auto lens = getOverlapLengths(box);
    return lens.first > 0 and lens.second > 0;
}