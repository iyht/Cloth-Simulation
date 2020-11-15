#include <dphi_cloth_triangle_dX.h>

//compute 3x3 deformation gradient 
void dphi_cloth_triangle_dX(Eigen::Matrix3d &dphi, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {
    Eigen::Matrix32d T;
    T.setZero();
    T.col(0) = (V.row(element(1)) - V.row(element(0))).transpose();
    T.col(1) = (V.row(element(2)) - V.row(element(0))).transpose();


    dphi.setZero();
    Eigen::MatrixXd bottom(2, 3);
    bottom << (T.transpose()*T).inverse()*T.transpose();
    dphi.row(0) = -bottom.colwise().sum();
    dphi.block(1, 0, 2, 3) = bottom;

    
}