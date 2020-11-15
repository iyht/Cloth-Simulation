#include <d2V_membrane_corotational_dq2.h>
#include <dV_membrane_corotational_dq.h>
#include <dphi_cloth_triangle_dX.h>
#include <iostream>

void d2V_membrane_corotational_dq2(Eigen::Matrix99d &H, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::Matrix3d> dX,
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double area, 
                          double mu, double lambda) {
    

    //SVD = USW^T
    Eigen::Matrix3d U;
    Eigen::Vector3d S; 
    Eigen::Matrix3d W; 
    Eigen::Matrix3d F; //deformation gradient
    
    double tol = 1e-5;
    
    //Compute SVD of F here
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
    double n_norm = n.norm();
    n = (1.0/n_norm) * n;

    Eigen::Matrix34d x_mat;
    x_mat.setZero();

    x_mat.block(0, 0, 3, 1) = q.segment(3*element(0),3);
    x_mat.block(0, 1, 3, 1) = q.segment(3*element(1),3);
    x_mat.block(0, 2, 3, 1) = q.segment(3*element(2),3);
    x_mat.block(0, 3, 3, 1) = n;

    Eigen::Matrix43d dX_and_N;
    dX_and_N.setZero();
    dX_and_N.block(0, 0, 3, 3) = dX;
    dX_and_N.block(3, 0, 1, 3) = N.transpose();

    F = x_mat * dX_and_N;

    //exit(0);
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(F, Eigen::ComputeThinU | Eigen::ComputeThinV);
    U << svd.matrixU();
    W << svd.matrixV();
    S << svd.singularValues();

    //deal with singularity in the svd gradient
    if(std::fabs(S[0] - S[1]) < tol || std::fabs(S[1] - S[2]) < tol || std::fabs(S[0] - S[2]) < tol) {
        F += Eigen::Matrix3d::Random()*tol;
        Eigen::JacobiSVD<Eigen::Matrix3d> svd2(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
        U = svd2.matrixU();
        W = svd2.matrixV();
        S = svd2.singularValues();
    }
    
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

    //TODO: compute H, the hessian of the corotational energy
    Eigen::Tensor3333d dU, dV;
    Eigen::Tensor333d dS;
    dsvd(dU, dS, dV, F);




    //dphi_dsigma
    Eigen::Vector3d dphi_dsigma;
    dphi_dsigma(0) = 2.0*mu*(-1.0+S[0])+lambda*(-3.0+S[0]+S[1]+S[2]);
    dphi_dsigma(1) = 2.0*mu*(-1.0+S[1])+lambda*(-3.0+S[0]+S[1]+S[2]);
    dphi_dsigma(2) = 2.0*mu*(-1.0+S[2])+lambda*(-3.0+S[0]+S[1]+S[2]);

    //d2phi_dsigma2
    Eigen::Matrix3d d2phi_dsigma2;
    d2phi_dsigma2.setIdentity();
    d2phi_dsigma2 *= 2.0*mu;
    d2phi_dsigma2 += Eigen::Matrix3d::Constant(lambda);
    //std::cout << "d2phi_dsigma2\n" << d2phi_dsigma2 << std::endl;

    Eigen::Vector3d tmp_Vec;
    Eigen::Matrix3d row_Mat;
    Eigen::Matrix99d d2phi_dF2;

    for(int r = 0; r <3; ++r) {
        for(int s = 0; s<3; ++s) {
            tmp_Vec  = d2phi_dsigma2*dS[r][s];
            row_Mat = (dU[r][s]*dphi_dsigma.asDiagonal()*W.transpose()  + U*tmp_Vec.asDiagonal()*W.transpose()) + U*dphi_dsigma.asDiagonal()*dV[r][s].transpose();
            row_Mat.transposeInPlace();
            d2phi_dF2.row(3*r + s) = Eigen::Map<Eigen::Matrix<double, 1,9> >(row_Mat.data(), 9);
        }
    }

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

    H = area*dF_dq.transpose()* d2phi_dF2 * dF_dq;


    //fix errant eigenvalues
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix99d> es(H);
    
    Eigen::MatrixXd DiagEval = es.eigenvalues().real().asDiagonal();
    Eigen::MatrixXd Evec = es.eigenvectors().real();
    
    for (int i = 0; i < 9; ++i) {
        if (es.eigenvalues()[i]<1e-6) {
            DiagEval(i,i) = 1e-3;
        }
    }
    H = Evec * DiagEval * Evec.transpose();


}
