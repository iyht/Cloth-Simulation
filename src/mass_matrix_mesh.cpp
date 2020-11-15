#include <mass_matrix_mesh.h>
#include <iostream>

void mass_matrix_mesh(Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::VectorXd> q, 
                         Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F,
                         double density, Eigen::Ref<const Eigen::VectorXd> areas) {
    int r = q.rows(); // r = 3n
    M.setZero();
    M.resize(r, r);

    typedef Eigen::Triplet<double> Trip;
    std::vector<Trip> tripleList;
    tripleList.reserve(r*r);

    for(int i = 0; i < F.rows(); i++)
    {


        // single triangle
        Eigen::Matrix99d M_i;
        double mass = density*areas(i);

        double c0 = (1.0/2.0)*mass;
        double c1 = (1.0/4.0)*mass;

        //setup the mass matrix
        //m_massMatrix
        //Point Indices
        unsigned int p0 = 0;
        unsigned int p1 = 3;
        unsigned int p2 = 6;

        //Assemble these bad boys //really big 4x4 block matrix
        M_i.block(p0,p0, 3,3) = c0*Eigen::Matrix<double,3,3>::Identity();
        M_i.block(p0,p1, 3,3) = c1*Eigen::Matrix<double,3,3>::Identity();
        M_i.block(p0,p2, 3,3) = c1*Eigen::Matrix<double,3,3>::Identity();

        M_i.block(p1,p0, 3,3) = c1*Eigen::Matrix<double,3,3>::Identity();
        M_i.block(p1,p1, 3,3) = c0*Eigen::Matrix<double,3,3>::Identity();
        M_i.block(p1,p2, 3,3) = c1*Eigen::Matrix<double,3,3>::Identity();

        M_i.block(p2,p0, 3,3) = c1*Eigen::Matrix<double,3,3>::Identity();
        M_i.block(p2,p1, 3,3) = c1*Eigen::Matrix<double,3,3>::Identity();
        M_i.block(p2,p2, 3,3) = c0*Eigen::Matrix<double,3,3>::Identity();

        // single triangle
        Eigen::Matrix3d M_i00 = M_i.block(0, 0, 3, 3);
        Eigen::Matrix3d M_i01 = M_i.block(0, 3, 3, 3);
        Eigen::Matrix3d M_i02 = M_i.block(0, 6, 3, 3);

        Eigen::Matrix3d M_i10 = M_i.block(3, 0, 3, 3);
        Eigen::Matrix3d M_i11 = M_i.block(3, 3, 3, 3);
        Eigen::Matrix3d M_i12 = M_i.block(3, 6, 3, 3);

        Eigen::Matrix3d M_i20 = M_i.block(6, 0, 3, 3);
        Eigen::Matrix3d M_i21 = M_i.block(6, 3, 3, 3);
        Eigen::Matrix3d M_i22 = M_i.block(6, 6, 3, 3);


        for(int ii = 0; ii < 3; ii++)
        {
            for (int jj = 0; jj < 3; jj++)
            {

                tripleList.push_back(Trip(3*F(i, 0)+ii, 3*F(i, 0)+jj, M_i00.coeff(ii, jj)));
                tripleList.push_back(Trip(3*F(i, 0)+ii, 3*F(i, 1)+jj, M_i01.coeff(ii, jj)));
                tripleList.push_back(Trip(3*F(i, 0)+ii, 3*F(i, 2)+jj, M_i02.coeff(ii, jj)));


                tripleList.push_back(Trip(3*F(i, 1)+ii, 3*F(i, 0)+jj, M_i10.coeff(ii, jj)));
                tripleList.push_back(Trip(3*F(i, 1)+ii, 3*F(i, 1)+jj, M_i11.coeff(ii, jj)));
                tripleList.push_back(Trip(3*F(i, 1)+ii, 3*F(i, 2)+jj, M_i12.coeff(ii, jj)));

                tripleList.push_back(Trip(3*F(i, 2)+ii, 3*F(i, 0)+jj, M_i20.coeff(ii, jj)));
                tripleList.push_back(Trip(3*F(i, 2)+ii, 3*F(i, 1)+jj, M_i21.coeff(ii, jj)));
                tripleList.push_back(Trip(3*F(i, 2)+ii, 3*F(i, 2)+jj, M_i22.coeff(ii, jj)));

            }
        }

    }
    M.setFromTriplets(tripleList.begin(), tripleList.end());

}
 
