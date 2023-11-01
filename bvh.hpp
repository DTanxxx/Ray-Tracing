#ifndef BVH_H
#define BVH_H

#include "constants.hpp"
#include "hittable.hpp"
#include "hittable_list.hpp"
#include <algorithm>

// generic bounding box comparator that returns true if the first argument is less than the second
// i.e. the first hittable's bounding box has a smaller coordinate along a given axis than the second
// hittable's bounding box's along that axis
inline bool box_compare(const shared_ptr<hittable> a, const shared_ptr<hittable> b, int axis) {
    aabb box_a;
    aabb box_b;

    if (!a->bounding_box(0, 0, box_a) || !b->bounding_box(0, 0, box_b)) {
        std::cerr << "No bounding box in bvh_node constructor.\n";
    }

    return box_a.min().e[axis] < box_b.min().e[axis];
}

// x axis specific comparator
bool box_x_compare(const shared_ptr<hittable> a, const shared_ptr<hittable> b) {
    return box_compare(a, b, 0);
}

// y axis specific comparator
bool box_y_compare(const shared_ptr<hittable> a, const shared_ptr<hittable> b) {
    return box_compare(a, b, 1);
}

// z axis specific comparator
bool box_z_compare(const shared_ptr<hittable> a, const shared_ptr<hittable> b) {
    return box_compare(a, b, 2);
}

class bvh_node : public hittable {
    public:
        // constructors
        bvh_node();
        bvh_node(const hittable_list& list, double time0, double time1)
            : bvh_node(list.objects, 0, list.objects.size(), time0, time1)
        {}
        // the central constructor: divide the list of hittable objects into two sublists,
        // one for each bvh child node (and repeat recursively), so that the two bvh child nodes
        // have smaller bounding boxes than their parent's bounding box
        bvh_node(const std::vector<shared_ptr<hittable>>& src_objects, size_t start, 
                 size_t end, double time0, double time1) {
            // create a modifiable array of the source scene objects
            auto objects = src_objects;

            // randomly choose an axis to split the list of primitives along
            int axis = random_int(0, 2);
            auto comparator = (axis == 0) ? box_x_compare :
                              (axis == 1) ? box_y_compare :
                                            box_z_compare;

            size_t object_span = end - start;

            if (object_span == 1) {
                // when the list of objects contains just one primitive, duplicate it so we have
                // one per subtree - this way we do not have to check for null pointers
                left = right = objects[start];
            }
            else if (object_span == 2) {
                // when the list of objects contains two primitives, we put one in each subtree and end the recursion
                if (comparator(objects[start], objects[start + 1])) {
                    left = objects[start];
                    right = objects[start + 1];
                }
                else {
                    left = objects[start + 1];
                    right = objects[start];
                }
            }
            else {
                // sort the primitives using std::sort, depending on which axis we chose to split along
                std::sort(objects.begin() + start, objects.begin() + end, comparator);

                // put half of primitives into each bvh subtree
                auto mid = start + object_span / 2;
                left = make_shared<bvh_node>(objects, start, mid, time0, time1);
                right = make_shared<bvh_node>(objects, mid, end, time0, time1);
            }

            aabb box_left, box_right;

            // check if there is a bounding box for left and right children bvh nodes, in case 
            // an infinite plane is passed into the tree which does not have a bounding box
            if (!left->bounding_box(time0, time1, box_left) || !right->bounding_box(time0, time1, box_right)) {
                std::cerr << "No bounding box in bvh_node constructor.\n";
            }

            box = surrounding_box(box_left, box_right);
        }

        virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override;
        virtual bool bounding_box(double time0, double time1, aabb& output_box) const override;

    public:
        shared_ptr<hittable> left;  // pointer to left generic hittable child
        shared_ptr<hittable> right;  // pointer to right generic hittable child
        aabb box;  // bounding box for this bvh node
};

// check whether the box for the bvh node is hit, and if so, check the children and sort out any details
bool bvh_node::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    if (!box.hit(r, t_min, t_max)) {
        // box for the bvh node is not hit, return false
        return false;
    }

    // check if the ray hits the left or right child
    bool hit_left = left->hit(r, t_min, t_max, rec);
    bool hit_right = right->hit(r, t_min, hit_left ? rec.t : t_max, rec);  // if the ray hits left child, check if the right child can be hit with a shorter ray (by setting t_max to rec.t)

    return hit_left || hit_right;
}

bool bvh_node::bounding_box(double time0, double time1, aabb& output_box) const {
    output_box = box;
    return true;
}

#endif
