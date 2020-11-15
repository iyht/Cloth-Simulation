#include <assemble_forces.h>
#include <iostream>

void assemble_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::MatrixXd> qdot, Eigen::Ref<const Eigen::MatrixXd> dX,
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F, Eigen::Ref<const Eigen::VectorXd> a0,
                     double mu, double lambda) {

    int r = q.rows();
    f.resize(r);
    f.setZero();

    for(int i = 0; i < F.rows(); i++)
    {
        Eigen::Vector3i index_q0, index_q1, index_q2, index_q3; // store the index to the vector q to reach [q0_x, q0_y, q0_z]
        index_q0 << F(i, 0)*3, F(i, 0)*3+1, F(i, 0)*3+2;
        index_q1 << F(i, 1)*3, F(i, 1)*3+1, F(i, 1)*3+2;
        index_q2 << F(i, 2)*3, F(i, 2)*3+1, F(i, 2)*3+2;
        //index_q3 << F(i, 3)*3, F(i, 3)*3+1, F(i, 3)*3+2;


        Eigen::Vector9d dV_i;
        Eigen::Matrix<double, 1,9> tmp_row;
        tmp_row = dX.row(i); //ei is the triangle index.
        Eigen::Map<const Eigen::Matrix3d> tmp_dX(tmp_row.data());
        //std::cout << "a0(i)\n" << a0(i) << std::endl;
        //std::cout << "F.row(i)\n" << F.row(i) << std::endl;
        dV_membrane_corotational_dq(dV_i, q, tmp_dX, V, F.row(i), a0(i), mu, lambda);

        f(index_q0(0)) -= dV_i(0);
        f(index_q0(1)) -= dV_i(1);
        f(index_q0(2)) -= dV_i(2);

        f(index_q1(0)) -= dV_i(3);
        f(index_q1(1)) -= dV_i(4);
        f(index_q1(2)) -= dV_i(5);

        f(index_q2(0)) -= dV_i(6);
        f(index_q2(1)) -= dV_i(7);
        f(index_q2(2)) -= dV_i(8);

        //f(index_q3(0)) -= dV_i(9);
        //f(index_q3(1)) -= dV_i(10);
        //f(index_q3(2)) -= dV_i(11);
    }

    //f.setZero();
        
       
}
