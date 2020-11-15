#include <dV_membrane_corotational_dq.h>
#include <dphi_cloth_triangle_dX.h>
#include <iostream>

void dV_membrane_corotational_dq(Eigen::Vector9d &dV, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::Matrix3d> dX, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double area, 
                          double mu, double lambda) {

    //Deformation Gradient
    Eigen::Matrix3d F;
    Eigen::Matrix3d U;
    Eigen::Vector3d S; 
    Eigen::Matrix3d W; 

    //TODO: SVD Here

    // q element
    Eigen::Vector9d X, x;
    X << V.row(element(0)).transpose(), V.row(element(1)).transpose(), V.row(element(2)).transpose();
    x << q.segment(3*element(0),3), q.segment(3*element(1),3), q.segment(3*element(2),3);

    Eigen::Vector3d delta_X2= X.segment(6, 3) - X.segment(0, 3);
    Eigen::Vector3d delta_X1= X.segment(3, 3) - X.segment(0, 3);
    Eigen::Vector3d N = (delta_X1).cross(delta_X2);
    N = (1.0/N.norm()) * N;

    Eigen::Vector3d delta_x2= x.segment(6, 3) - x.segment(0, 3);
    Eigen::Vector3d delta_x1= x.segment(3, 3) - x.segment(0, 3);
    Eigen::Vector3d n = (delta_x1).cross(delta_x2);
    double n_norm =  n.norm();
    n = (1.0/n_norm) * n;

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

    F = x_mat * dX_and_N;

    //std::cout << "F\n" << F << std::endl;
    //exit(0);
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(F, Eigen::ComputeThinU | Eigen::ComputeThinV);
    U << svd.matrixU();
    W << svd.matrixV();
    S << svd.singularValues();




    //Eigen::JacobiSVD<Eigen::MatrixXd> svd(dX, Eigen::ComputeThinU | Eigen::ComputeThinV);
    //U << svd.matrixU();
    //W << svd.matrixV();
    //S << svd.singularValues();

    //Fix for inverted elements (thanks to Danny Kaufman)
    double det = S[0]*S[1];
    
     if(det <= -1e-10)
    {
        if(S[0] < 0) S[0] *= -1;
        if(S[1] < 0) S[1] *= -1;
        if(S[2] < 0) S[2] *= -1;
    }
    
    if(U.determinant() <= 0)
    {
        U(0, 2) *= -1;
        U(1, 2) *= -1;
        U(2, 2) *= -1;
    }
    
    if(W.determinant() <= 0)
    {
        W(0, 2) *= -1;
        W(1, 2) *= -1;
        W(2, 2) *= -1;
    }
    
    //TODO: energy model gradient
    Eigen::Vector3d dphi_dsigma;
    dphi_dsigma.setZero();
    dphi_dsigma(0) = 2.0*mu*(-1.0+S[0])+lambda*(-3.0+S[0]+S[1]+S[2]);
    dphi_dsigma(1) = 2.0*mu*(-1.0+S[1])+lambda*(-3.0+S[0]+S[1]+S[2]);
    dphi_dsigma(2) = 2.0*mu*(-1.0+S[2])+lambda*(-3.0+S[0]+S[1]+S[2]);
    //dphi_dsigma << lambda*(S.sum() - 3.0) + (mu * 2.0 * (S.array()-1.0).sum());
    Eigen::Matrix<double, 3, 3, Eigen::RowMajor> dphi_dF = U * dphi_dsigma.asDiagonal() * W.transpose();

    // q element
    //Eigen::Vector9d x, X;
    //X << V.row(element(0)).transpose(), V.row(element(1)).transpose(), V.row(element(2)).transpose();
    //x << q.segment(3*element(0),3), q.segment(3*element(1),3), q.segment(3*element(2),3);

    //Eigen::Vector3d delta_X2= X.segment(6, 3) - X.segment(0, 3);
    //Eigen::Vector3d delta_X1= X.segment(3, 3) - X.segment(0, 3);
    //Eigen::Vector3d N = (delta_X2).cross(delta_X1);
    //N = N/N.norm();

    //Eigen::Vector3d delta_x2= x.segment(6, 3) - x.segment(0, 3);
    //Eigen::Vector3d delta_x1= x.segment(3, 3) - x.segment(0, 3);
    //Eigen::Vector3d n = (delta_x2).cross(delta_x1);
    //n = n/n.norm();

    Eigen::Matrix39d n_gradient;
    n_gradient.setZero();
    n_gradient(0,0) = -n(0)*n(2)*(delta_x1(1)-delta_x2(1))+n(0)*n(1)*(delta_x1(2)-delta_x2(2));
    n_gradient(0,1) = -(n(0)*n(0)-1.0)*(delta_x1(2)-delta_x2(2))+n(0)*n(2)*(delta_x1(0)-delta_x2(0));
    n_gradient(0,2) = (n(0)*n(0)-1.0)*(delta_x1(1)-delta_x2(1))-n(0)*n(1)*(delta_x1(0)-delta_x2(0));
    n_gradient(0,3) = n(0)*n(1)*delta_x2(2)-n(0)*n(2)*delta_x2(1);
    n_gradient(0,4) = -delta_x2(2)*(n(0)*n(0)-1.0)+n(0)*n(2)*delta_x2(0);
    n_gradient(0,5) = delta_x2(1)*(n(0)*n(0)-1.0)-n(0)*n(1)*delta_x2(0);
    n_gradient(0,6) = -n(0)*n(1)*delta_x1(2)+n(0)*n(2)*delta_x1(1);
    n_gradient(0,7) = delta_x1(2)*(n(0)*n(0)-1.0)-n(0)*n(2)*delta_x1(0);
    n_gradient(0,8) = -delta_x1(1)*(n(0)*n(0)-1.0)+n(0)*n(1)*delta_x1(0);
    n_gradient(1,0) = (n(1)*n(1)-1.0)*(delta_x1(2)-delta_x2(2))-n(1)*n(2)*(delta_x1(1)-delta_x2(1));
    n_gradient(1,1) = n(1)*n(2)*(delta_x1(0)-delta_x2(0))-n(0)*n(1)*(delta_x1(2)-delta_x2(2));
    n_gradient(1,2) = -(n(1)*n(1)-1.0)*(delta_x1(0)-delta_x2(0))+n(0)*n(1)*(delta_x1(1)-delta_x2(1));
    n_gradient(1,3) = delta_x2(2)*(n(1)*n(1)-1.0)-n(1)*n(2)*delta_x2(1);
    n_gradient(1,4) = -n(0)*n(1)*delta_x2(2)+n(1)*n(2)*delta_x2(0);
    n_gradient(1,5) = -delta_x2(0)*(n(1)*n(1)-1.0)+n(0)*n(1)*delta_x2(1);
    n_gradient(1,6) = -delta_x1(2)*(n(1)*n(1)-1.0)+n(1)*n(2)*delta_x1(1);
    n_gradient(1,7) = n(0)*n(1)*delta_x1(2)-n(1)*n(2)*delta_x1(0);
    n_gradient(1,8) = delta_x1(0)*(n(1)*n(1)-1.0)-n(0)*n(1)*delta_x1(1);
    n_gradient(2,0) = -(n(2)*n(2)-1.0)*(delta_x1(1)-delta_x2(1))+n(1)*n(2)*(delta_x1(2)-delta_x2(2));
    n_gradient(2,1) = (n(2)*n(2)-1.0)*(delta_x1(0)-delta_x2(0))-n(0)*n(2)*(delta_x1(2)-delta_x2(2));
    n_gradient(2,2) = -n(1)*n(2)*(delta_x1(0)-delta_x2(0))+n(0)*n(2)*(delta_x1(1)-delta_x2(1));
    n_gradient(2,3) = -delta_x2(1)*(n(2)*n(2)-1.0)+n(1)*n(2)*delta_x2(2);
    n_gradient(2,4) = delta_x2(0)*(n(2)*n(2)-1.0)-n(0)*n(2)*delta_x2(2);
    n_gradient(2,5) = n(0)*n(2)*delta_x2(1)-n(1)*n(2)*delta_x2(0);
    n_gradient(2,6) = delta_x1(1)*(n(2)*n(2)-1.0)-n(1)*n(2)*delta_x1(2);
    n_gradient(2,7) = -delta_x1(0)*(n(2)*n(2)-1.0)+n(0)*n(2)*delta_x1(2);
    n_gradient(2,8) = -n(0)*n(2)*delta_x1(1)+n(1)*n(2)*delta_x1(0);
    n_gradient = (1.0/n_norm)*n_gradient;

    Eigen::Matrix3d dphi;
    dphi.setZero();
    dphi_cloth_triangle_dX(dphi, V, element, n); // n is not really use here
    Eigen::MatrixXd B(9, 9);
    B.setZero();
    B.block(0, 0, 3, 1) = dphi.block(0, 0, 1, 3).transpose();
    B.block(3, 1, 3, 1) = dphi.block(0, 0, 1, 3).transpose();
    B.block(6, 2, 3, 1) = dphi.block(0, 0, 1, 3).transpose();

    B.block(0, 3, 3, 1) = dphi.block(1, 0, 1, 3).transpose();
    B.block(3, 4, 3, 1) = dphi.block(1, 0, 1, 3).transpose();
    B.block(6, 5, 3, 1) = dphi.block(1, 0, 1, 3).transpose();

    B.block(0, 6, 3, 1) = dphi.block(2, 0, 1, 3).transpose();
    B.block(3, 7, 3, 1) = dphi.block(2, 0, 1, 3).transpose();
    B.block(6, 8, 3, 1) = dphi.block(2, 0, 1, 3).transpose();

    Eigen::Matrix93d N_mat;
    N_mat.setZero();
    N_mat.block(0, 0, 3, 1) = N;
    N_mat.block(3, 1, 3, 1) = N;
    N_mat.block(6, 2, 3, 1) = N;

    Eigen::Matrix99d dF_dq = B + N_mat*n_gradient;
    //Eigen::Vector9d dphi_dF_flatten;
    Eigen::Map<Eigen::Vector9d> dphi_dF_flatten(dphi_dF.data(), dphi_dF.size());

    dV = area*dF_dq.transpose()*dphi_dF_flatten;



}
