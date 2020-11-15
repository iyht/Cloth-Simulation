#include <velocity_filter_cloth_sphere.h>
#include <iostream>


void velocity_filter_cloth_sphere(Eigen::VectorXd &qdot, const std::vector<unsigned int> &indices, 
                                  const std::vector<Eigen::Vector3d> &normals) {

    for(int i = 0; i < indices.size(); i++)
    {
        Eigen::Vector3d v_i = qdot.segment(indices[i], 3);
        Eigen::Vector3d n_i = normals[i];
        //std::cout << "n\n" << n_i << std::endl;

        double nv = n_i.transpose()*v_i;
        double alpha = 0.0 <= nv ? 0.0 : nv;
        qdot.segment(indices[i], 3) -= alpha*n_i;
    }
    

}