#include <dV_cloth_gravity_dq.h>
#include <iostream>

void dV_cloth_gravity_dq(Eigen::VectorXd &fg, Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::Vector3d> g) {
    Eigen::VectorXd g_all(M.rows());
    for(int i = 0; i < M.rows(); i += 3)
    {
        g_all.segment(i, 3) = g;
    }

    fg = -M*g_all;
    std::cout << "M\n" << M.block(0,0,20,20) << std::endl;
    std::cout << "fg\n" << fg.block(0,0,20,1) << std::endl;
    //exit(0);

}
