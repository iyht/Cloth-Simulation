#include <V_membrane_corotational.h>

//Allowed to use libigl SVD or Eigen SVD for this part
void V_membrane_corotational(double &energy, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::Matrix3d> dX, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double area, 
                          double mu, double lambda) {

    // get the deformation gradient

    // q element
    Eigen::Vector9d X, x;
    X << V.row(element(0)).transpose(), V.row(element(1)).transpose(), V.row(element(2)).transpose();
    x << q.segment(3*element(0),3), q.segment(3*element(1),3), q.segment(3*element(2),3);

    Eigen::Vector3d delta_X2= X.segment(6, 3) - X.segment(0, 3);
    Eigen::Vector3d delta_X1= X.segment(3, 3) - X.segment(0, 3);
    Eigen::Vector3d N = (delta_X1).cross(delta_X2);
    N = N/N.norm();

    Eigen::Vector3d delta_x2= x.segment(6, 3) - x.segment(0, 3);
    Eigen::Vector3d delta_x1= x.segment(3, 3) - x.segment(0, 3);
    Eigen::Vector3d n = (delta_x1).cross(delta_x2);
    n = n/n.norm();
    //std::cout << "n\n" << n << std::endl;

    Eigen::Matrix34d x_mat;
    x_mat.setZero();

    x_mat.block(0, 0, 3, 1) = q.segment(3*element(0),3);
    x_mat.block(0, 1, 3, 1) = q.segment(3*element(1),3);
    x_mat.block(0, 2, 3, 1) = q.segment(3*element(2),3);
    x_mat.block(0, 3, 3, 1) = n;
    //std::cout << "x_mat\n" << x_mat << std::endl;

    Eigen::Matrix43d dX_and_N;
    dX_and_N.setZero();
    dX_and_N.block(0, 0, 3, 3) = dX;
    dX_and_N.block(3, 0, 1, 3) = N.transpose();
    //std::cout << "dX_and_N\n" << dX_and_N << std::endl;

    Eigen::Matrix3d F = x_mat * dX_and_N;


    // SVD of the deformation gradient to get the pricipal stretch
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(F, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::Matrix3d svd_U, svd_V;
    Eigen::Vector3d svd_S;
    svd_U << svd.matrixU();
    svd_V << svd.matrixV();
    svd_S << svd.singularValues();

    // build the Co-Rotational Linear Elasticity
    energy = area * (mu * ((svd_S[0] - 1.0)*(svd_S[0] - 1.0) + (svd_S[1] - 1.0)*(svd_S[1] - 1.0) + (svd_S[2] - 1.0)*(svd_S[2] - 1.0)) +
             0.5*lambda*(svd_S[0] + svd_S[1] + svd_S[2] - 3.0)*(svd_S[0] + svd_S[1] + svd_S[2] - 3.0));



}
