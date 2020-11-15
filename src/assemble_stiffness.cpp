#include <assemble_stiffness.h>
#include <iostream>
void assemble_stiffness(Eigen::SparseMatrixd &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::MatrixXd> dX,
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F, Eigen::Ref<const Eigen::VectorXd> a0, 
                     double mu, double lambda) {

    int r = qdot.rows(); // r = 3n
    K.setZero();
    K.resize(r, r);

    typedef Eigen::Triplet<double> Trip;
    std::vector<Trip> tripleList;
    tripleList.reserve(r*r);

    for(int i = 0; i < F.rows(); i++)
    {

        Eigen::Matrix99d H_i;

        Eigen::Matrix<double, 1,9> tmp_row;
        tmp_row = dX.row(i); //ei is the triangle index.
        Eigen::Map<const Eigen::Matrix3d> dX_row(tmp_row.data());


        d2V_membrane_corotational_dq2(H_i, q, dX_row, V, F.row(i), a0(i), mu, lambda);


        Eigen::Matrix3d H_i00 = H_i.block(0, 0, 3, 3);
        Eigen::Matrix3d H_i01 = H_i.block(0, 3, 3, 3);
        Eigen::Matrix3d H_i02 = H_i.block(0, 6, 3, 3);

        Eigen::Matrix3d H_i10 = H_i.block(3, 0, 3, 3);
        Eigen::Matrix3d H_i11 = H_i.block(3, 3, 3, 3);
        Eigen::Matrix3d H_i12 = H_i.block(3, 6, 3, 3);

        Eigen::Matrix3d H_i20 = H_i.block(6, 0, 3, 3);
        Eigen::Matrix3d H_i21 = H_i.block(6, 3, 3, 3);
        Eigen::Matrix3d H_i22 = H_i.block(6, 6, 3, 3);


        for(int ii = 0; ii < 3; ii++)
        {
            for (int jj = 0; jj < 3; jj++)
            {

                tripleList.push_back(Trip(3*F(i, 0)+ii, 3*F(i, 0)+jj, -H_i00.coeff(ii, jj)));
                tripleList.push_back(Trip(3*F(i, 0)+ii, 3*F(i, 1)+jj, -H_i01.coeff(ii, jj)));
                tripleList.push_back(Trip(3*F(i, 0)+ii, 3*F(i, 2)+jj, -H_i02.coeff(ii, jj)));


                tripleList.push_back(Trip(3*F(i, 1)+ii, 3*F(i, 0)+jj, -H_i10.coeff(ii, jj)));
                tripleList.push_back(Trip(3*F(i, 1)+ii, 3*F(i, 1)+jj, -H_i11.coeff(ii, jj)));
                tripleList.push_back(Trip(3*F(i, 1)+ii, 3*F(i, 2)+jj, -H_i12.coeff(ii, jj)));

                tripleList.push_back(Trip(3*F(i, 2)+ii, 3*F(i, 0)+jj, -H_i20.coeff(ii, jj)));
                tripleList.push_back(Trip(3*F(i, 2)+ii, 3*F(i, 1)+jj, -H_i21.coeff(ii, jj)));
                tripleList.push_back(Trip(3*F(i, 2)+ii, 3*F(i, 2)+jj, -H_i22.coeff(ii, jj)));


            }
        }

    }
    K.setFromTriplets(tripleList.begin(), tripleList.end());
       
        
    };
