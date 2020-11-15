#include <dV_cloth_gravity_dq.h>
#include <iostream>

void dV_cloth_gravity_dq(Eigen::VectorXd &fg, Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::Vector3d> g) {
    Eigen::VectorXd g_all(M.rows());
    for(int i = 0; i < M.rows(); i += 3)
    {
        g_all.segment(i, 3) = g;
    }

    fg = -M*g_all;

}
