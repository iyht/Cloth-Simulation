#include <collision_detection_cloth_sphere.h>
#include <iostream>
void collision_detection_cloth_sphere(std::vector<unsigned int> &cloth_index, std::vector<Eigen::Vector3d> &normals, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::Vector3d> center, double radius) {


    // for each vertex x, check if || x - c ||^2 < r^2, where c is center, r is radius
    cloth_index.clear();
    normals.clear();

    for(int i = 0; i < q.rows(); i += 3)
    {
        Eigen::Vector3d q_i = q.segment(i, 3);
        // check
        if((q_i - center).squaredNorm() <= (radius*radius))
        {
            cloth_index.push_back(i);
            Eigen::Vector3d n = (q_i - center) /(q_i - center).norm();
            normals.push_back(n);
        }
    }



}