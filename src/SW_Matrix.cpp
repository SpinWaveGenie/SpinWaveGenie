//
//  CreateMatrix.cpp
//  Spin Wave Fit
//
//  Created by Hahn, Steven E. on 1/7/13.
//  Copyright (c) 2013 Oak Ridge National Laboratory. All rights reserved.
//

#include "SW_Matrix.h"
#include <iomanip>

using namespace Eigen;
using namespace std;

/*int SW_Matrix::mod(int K)
{
    //need to ensure that modulus is greater than or equal to zero
    //http://www.worldgeekz.com/2012/05/cc-positive-modulo-function.html
    return (K%P + P)%P;
}*/

void SW_Matrix::set_parr(vector<SW_sublattice> SLin,vector<double>& Xin)
{
    X = Xin;
    SL = SLin;
    
    M= (int)SL.size();
    N=2*M;

    LN.resize(N,N); LN.setZero();
    
    SS.resize(N);
    SS.setZero();
    for (int j=0;j<M;j++)
    {
        SS(j) = 1.0;
        SS(j+M) = -1.0;
    }
}

vector<double> SW_Matrix::get_parr()
{
    return X;
}

void SW_Matrix::CreateMatrix_exchange( double KXP, double KYP, double KZP, MatrixXi &interactions)
{
    double Sr, z_rs;
    Matrix3d F;
    Vector3d K;
    vector<Vector3d> neighbors;
    
    K[0] = KXP*2.0*M_PI;
    K[1] = KYP*2.0*M_PI;
    K[2] = KZP*2.0*M_PI;
    
    for (int r=0;r<M;r++)
    {
        Sr = SL[r].get_sublattice()[0];
        for (int s=0;s<M;s++)
        {
            //cout << r << '\t' << s << endl;
            if (interactions(r,s) >= 0)
            {
                F = SL[r].get_rot_matrix()*SL[s].get_inv_matrix();
                //cout << q << "\t" << s << endl << F << endl;
                //cout << endl;
                neighbors = SL[r].get_neighbors(s);
                z_rs = (double) neighbors.size(); //avoid double-counting bonds
            
                complex<double> gamma_rs (0.0,0.0);
                for (int i=0;i<neighbors.size();i++)
                {
                    double dot_prod = K.dot(neighbors[i]);
                    gamma_rs += complex<double> (cos(dot_prod),-1.0*sin(dot_prod));
                }
                gamma_rs /= z_rs; //force gamma_rs(k=0) = 1.0
                //cout << "gamma_rs(" << r << "," << s << ")= " << gamma_rs << endl;
                
                complex<double> G1 (F(0,0) + F(1,1),F(1,0)-F(0,1));
                G1 *= -0.5;
                
                complex<double> G2 (F(0,0) - F(1,1),-F(1,0)-F(0,1));
                G2 *= -0.5;
                
                cout << "G1= " << G1 << endl;
		        cout << "G2= " << G2 << endl;
                {
                    LN(r,r) += 0.25*z_rs*X[interactions(r,s)]*Sr*F(2,2);
                    LN(s,s) += 0.25*z_rs*X[interactions(r,s)]*Sr*F(2,2);
                    LN(r+M,r+M) += 0.25*z_rs*X[interactions(r,s)]*Sr*F(2,2);
                    LN(s+M,s+M) += 0.25*z_rs*X[interactions(r,s)]*Sr*F(2,2);
                    
                    LN(r,s) += 0.25*z_rs*X[interactions(r,s)]*Sr*gamma_rs*G1;
                    LN(r+M,s+M) += 0.25*z_rs*X[interactions(r,s)]*Sr*conj(gamma_rs)*conj(G1);
                    LN(s,r) += 0.25*z_rs*X[interactions(r,s)]*Sr*gamma_rs*conj(G1);
                    LN(s+M,r+M) += 0.25*z_rs*X[interactions(r,s)]*Sr*conj(gamma_rs)*G1;
                    
                    LN(r,s+M) += 0.25*z_rs*X[interactions(r,s)]*Sr*conj(gamma_rs)*conj(G2);
                    LN(s,r+M) += 0.25*z_rs*X[interactions(r,s)]*Sr*gamma_rs*conj(G2);
                    
                    LN(r+M,s) += 0.25*z_rs*X[interactions(r,s)]*Sr*conj(gamma_rs)*G2;
                    LN(s+M,r) += 0.25*z_rs*X[interactions(r,s)]*Sr*gamma_rs*G2;
                }
                
            }
        }
    }
}


void SW_Matrix::CreateMatrix_DMy( double KXP, double KYP, double KZP, MatrixXi &interactions)
{
    /* diagonal terms in exchange*/
    complex<double> XI (0.0,1.0);
    double tmp;
    complex<double> tmp1,tmp2,tmp3;
    Vector3d K;
    double Sr, theta_r, phi_r,Ss, theta_s, phi_s, z_rs;
    vector<Vector3d> neighbors;

    
    K[0] = KXP*2.0*M_PI;
    K[1] = KYP*2.0*M_PI;
    K[2] = KZP*2.0*M_PI;
     
    for (int r=0;r<M;r++)
    {
        Sr = SL[r].get_sublattice()[0];
        theta_r = SL[r].get_sublattice()[1];
        phi_r = SL[r].get_sublattice()[2];
        for (int s=0;s<M;s++)
        {
            //cout << r << '\t' << s << endl;
            if (interactions(r,s) >= 0)
            {
                Ss = SL[s].get_sublattice()[0];
                theta_s = SL[s].get_sublattice()[1];
                phi_s = SL[s].get_sublattice()[2];
                
                //vector<Vector3d> neighbors = SL[r].get_neighbors(s);
                
                //double z_rs = (double) neighbors.size(); //avoid double-counting bonds
                
                neighbors = SL[r].get_neighbors(s);
                z_rs = (double) neighbors.size(); //avoid double-counting bonds
                
                complex<double> gamma_rs (0.0,0.0);
                for (int i=0;i<neighbors.size();i++)
                {
                    double dot_prod = K.dot(neighbors[i]);
                    gamma_rs += complex<double> (cos(dot_prod),-1.0*sin(dot_prod));
                }
                gamma_rs /= z_rs; //force gamma_rs(k=0) = 1.0
                //cout << "gamma_rs(" << r << "," << s << ")= " << gamma_rs << endl;

                tmp = 0.5*X[interactions(r,s)]*Sr*z_rs*(sin(theta_r)*cos(theta_s)*cos(phi_r) - cos(theta_r)*sin(theta_s)*cos(phi_s));
                LN(r,r) -= tmp;
                LN(r+M,r+M) -= tmp;
                LN(s,s) -= tmp;
                LN(s+M,s+M) -= tmp;
                
                tmp1 = cos(theta_r)*sin(theta_s)*cos(phi_r)-sin(theta_r)*cos(theta_s)*cos(phi_s);
                tmp2 = sin(theta_r)*sin(phi_s);
                tmp3 = sin(theta_s)*sin(phi_r);
                
                //cout << "r= " << r << ", s= " << s << endl;
                //cout << tmp1 << '\t' << tmp2 << '\t' << tmp3 << endl;
                
                LN(r,s) -= 0.25*X[interactions(r,s)]*Sr*z_rs*conj(gamma_rs)*(tmp1 - XI*tmp2 - XI*tmp3);
                LN(s,r) -= 0.25*X[interactions(r,s)]*Sr*z_rs*gamma_rs*(tmp1 + XI*tmp2 + XI*tmp3);
                LN(r+M,s+M) -= 0.25*X[interactions(r,s)]*z_rs*conj(gamma_rs)*Sr*(tmp1 + XI*tmp2 + XI*tmp3);
                LN(s+M,r+M) -= 0.25*X[interactions(r,s)]*z_rs*gamma_rs*Sr*(tmp1 - XI*tmp2 - XI*tmp3);
                
                LN(r+M,s) -= 0.25*X[interactions(r,s)]*z_rs*conj(gamma_rs)*Sr*(tmp1 - XI*tmp2 + XI*tmp3);
                LN(s+M,r) -= 0.25*X[interactions(r,s)]*z_rs*gamma_rs*Sr*(tmp1 - XI*tmp2 + XI*tmp3);
                
                LN(r,s+M) -= 0.25*X[interactions(r,s)]*z_rs*conj(gamma_rs)*Sr*(tmp1 + XI*tmp2 - XI*tmp3);
                LN(s,r+M) -= 0.25*X[interactions(r,s)]*z_rs*gamma_rs*Sr*(tmp1 + XI*tmp2 - XI*tmp3);
            }
        }
    }
}

void SW_Matrix::CreateMatrix_DMz( double KXP, double KYP, double KZP, MatrixXi &interactions)
{
    /* diagonal terms in exchange*/
    complex<double> XI (0.0,1.0);
    double tmp;
    complex<double> tmp1,tmp2,tmp3,tmp4;
    double Sr, theta_r, phi_r,Ss, theta_s, phi_s, z_rs;
    vector<Vector3d> neighbors;
    Vector3d K;
    
    K[0] = KXP*2.0*M_PI;
    K[1] = KYP*2.0*M_PI;
    K[2] = KZP*2.0*M_PI;
    
    for (int r=0;r<M;r++)
    {
        Sr = SL[r].get_sublattice()[0];
        theta_r = SL[r].get_sublattice()[1];
        phi_r = SL[r].get_sublattice()[2];
        for (int s=0;s<M;s++)
        {
            //cout << r << '\t' << s << endl;
            if (interactions(r,s) >= 0)
            {
                Ss = SL[s].get_sublattice()[0];
                theta_s = SL[s].get_sublattice()[1];
                phi_s = SL[s].get_sublattice()[2];
                
                neighbors = SL[r].get_neighbors(s);
                z_rs = (double) neighbors.size();
                
                complex<double> gamma_rs (0.0,0.0);
                for (int i=0;i<neighbors.size();i++)
                {
                    double dot_prod = K.dot(neighbors[i]);
                    gamma_rs += complex<double> (cos(dot_prod),-1.0*sin(dot_prod));
                }
                gamma_rs /= z_rs; //force gamma_rs(k=0) = 1.0
                
                tmp = 0.5*X[interactions(r,s)]*Sr*z_rs*sin(theta_r)*sin(theta_s)*sin(phi_r-phi_s);
                //cout << tmp << endl;
                LN(r,r) -= tmp;
                LN(r+M,r+M) -= tmp;
                LN(s,s) -= tmp;
                LN(s+M,s+M) -= tmp;
                
                tmp1 = cos(theta_r)*cos(theta_s)*sin(phi_r-phi_s);
                tmp2 = cos(theta_r)*cos(phi_r-phi_s);
                tmp3 = cos(theta_s)*cos(phi_r-phi_s);
                tmp4 = sin(phi_r-phi_s);
                
                //cout << "r= " << r << ", s= " << s << endl;
                //cout << tmp1 << '\t' << tmp2 << '\t' << tmp3 << '\t' << tmp4 <<  endl;
                //cout << phi_r-phi_s << endl;
                
                
                LN(r,s) -= 0.25*X[interactions(r,s)]*Sr*z_rs*conj(gamma_rs)*(-tmp1 - XI*tmp2 - XI*tmp3 - tmp4);
                LN(s,r) -= 0.25*X[interactions(r,s)]*Sr*z_rs*gamma_rs*(-tmp1 + XI*tmp2 + XI*tmp3 - tmp4);
                
                LN(r+M,s+M) -= 0.25*X[interactions(r,s)]*Sr*z_rs*conj(gamma_rs)*(-tmp1 + XI*tmp2 + XI*tmp3 - tmp4);
                LN(s+M,r+M) -= 0.25*X[interactions(r,s)]*Sr*z_rs*gamma_rs*(-tmp1 - XI*tmp2 - XI*tmp3 - tmp4);
                
                LN(r+M,s) -= 0.25*X[interactions(r,s)]*Sr*z_rs*conj(gamma_rs)*(-tmp1 - XI*tmp2 + XI*tmp3 + tmp4);
                LN(s+M,r) -= 0.25*X[interactions(r,s)]*Sr*z_rs*gamma_rs*(-tmp1 - XI*tmp2 + XI*tmp3 + tmp4);
                
                LN(r,s+M) -= 0.25*X[interactions(r,s)]*Sr*z_rs*conj(gamma_rs)*(-tmp1 + XI*tmp2 - XI*tmp3 + tmp4);
                LN(s,r+M) -= 0.25*X[interactions(r,s)]*Sr*z_rs*gamma_rs*(-tmp1 + XI*tmp2 - XI*tmp3 + tmp4);
            }
        }
    }
}

void SW_Matrix::CreateMatrix_anis_z()
{
    for (int q=0;q<M;q++)
    {
        double Sq = SL[q].get_sublattice()[0];
        double theta_q = SL[q].get_sublattice()[1];
        LN(q,q) -= 0.5*X[2]*Sq*(1.0 - 3.0*pow(cos(theta_q),2));
        LN(q+M,q+M) -= 0.5*X[2]*Sq*(1.0 - 3.0*pow(cos(theta_q),2));
        LN(q,q+M) -= 0.5*X[2]*Sq*pow(sin(theta_q),2);
        LN(q+M,q) -= 0.5*X[2]*Sq*pow(sin(theta_q),2);
    }
}

void SW_Matrix::CreateMatrix_anis_x()
{
    complex<double> XI (0.0,1.0);
    
    for (int r=0;r<M;r++)
    {
        double S = SL[r].get_sublattice()[0];
        double theta = SL[r].get_sublattice()[1];
        double phi = SL[r].get_sublattice()[2];
        //cout << r << endl << UI << endl;
        LN(r,r)     -= 0.5*X[3]*S*(pow(cos(theta),2)*pow(cos(phi),2)+pow(sin(phi),2)-2.0*pow(sin(theta),2)*pow(cos(phi),2));
        LN(r+M,r+M) -= 0.5*X[3]*S*(pow(cos(theta),2)*pow(cos(phi),2)+pow(sin(phi),2)-2.0*pow(sin(theta),2)*pow(cos(phi),2));
        LN(r+M,r) -= 0.5*X[3]*S*pow(cos(theta)*cos(phi)+XI*sin(phi),2);
        LN(r,r+M) -= 0.5*X[3]*S*pow(cos(theta)*cos(phi)-XI*sin(phi),2);
    }
    
}

void SW_Matrix::CreateMatrix_bfield()
{
    for (int q=0;q<M;q++)
    {
        double Sq = SL[q].get_sublattice()[0];
        double theta_q = SL[q].get_sublattice()[1];
        LN(q,q) -= 0.5*X[3]*Sq*cos(theta_q);
        LN(q+M,q+M) -= 0.5*X[3]*Sq*cos(theta_q);
    }
}

void SW_Matrix::CreateMatrix_YFeO3( double KXP, double KYP, double KZP)
{
    double J,D,Kx,Kz,z,S,theta;
    double eta_q,gamma_q;
    double A_q,B_q,C_q,D_q;
    double KX = KXP*2.0*M_PI;
    double KY = KYP*2.0*M_PI;
    double KZ = KZP*2.0*M_PI;

    J = -X[0]/2.0;
    D = X[4];
    Kx = X[3];
    Kz = X[2];
    z=6.0;
    S=SL[0].get_sublattice()[0];
    double delta = abs(D/(2.0*J));
    //cout << delta << endl;

    theta = SL[0].get_sublattice()[1]; //M_PI/2.0 - 0.01098;
    
    cout << M_PI/2.0 - theta << endl;

    gamma_q = (cos(KZ/2.0) + cos((KX+KY)/2.0) + cos((KX-KY)/2.0))/3.0;

    A_q = -z*J*S*cos(2.0*theta) - 0.5*z*D*S*sin(2.0*theta) - 0.5*Kx*S*(pow(cos(theta),2) - 2.0*pow(sin(theta),2)) -0.5*Kz*S*(1.0 - 3.0*pow(cos(theta),2));

    B_q = -0.5*Kx*S*pow(cos(theta),2) -0.5*Kz*S*pow(sin(theta),2);

    C_q = -z*J*S*pow(cos(theta),2)*gamma_q - 0.25*z*D*S*sin(2.0*theta)*gamma_q;

    D_q = z*J*S*pow(sin(theta),2)*gamma_q - 0.25*z*D*S*sin(2.0*theta)*gamma_q;

    MatrixXcd ML;
    ML.resize(4,4);
    ML.setZero();

    ML(0,0) = A_q;
    ML(0,1) = C_q;
    ML(0,2) = B_q;
    ML(0,3) = D_q;
    ML(1,0) = C_q;
    ML(1,1) = A_q;
    ML(1,2) = D_q;
    ML(1,3) = B_q;

    ML.block(M,M,M,M) = ML.block(0,0,M,M);
    ML.block(M,0,M,M) = ML.block(0,M,M,M);

    cout << ML << endl;
    cout << LN << endl;
    
    //LN = ML;
    
    //cout << sqrt(24.0*J*S*S*(2.0*(Kx-Kz))) << '\t';
    //cout << sqrt(24.0*J*S*S*(6.0*abs(D)*tan(abs(delta)) + 2.0*Kx)) << endl;
    
    cout << sqrt(24.0*J*S*2.0*S*(Kx-Kz)) << endl;
    cout << sqrt(24.0*J*S*2.0*Kx*S) << endl;
}

void SW_Matrix::Calc_Eigenvalues()
{
    int i;
    //MatrixXcd test;
    // testing if LN is a normal matrix
    // http://en.wikipedia.org/wiki/Normal_matrix
    //test = LN.adjoint()*LN - LN*LN.adjoint();
    //cout << test.squaredNorm() << "\n";
    //    for (int L1=0;L1<N;L1++)
    //    {
    //        for (int L2=0;L2<N;L2++)
    //        {
    //        if (abs(test(L1,L2)) > 1.0E-6)
    //            cout << "AR and AI matrices: " << L1 << "   " << L2 << "    " <<test(L1,L2) << "\n";
    //        }
    //    }
    
    
    //LN.block(M,0,M,M) = LN.block(0,M,M,M);
    //LN.block(M,M,M,M) = LN.block(0,0,M,M);

    //cout << LN << endl;
    
    //cout << "M Matrix:" << endl;
    //cout << LN;
    //cout << endl;
    //cout << LN - LN.adjoint() << endl;
    //cout << endl;

    LN.block(0,M,M,M) *= -1.0;
    LN.block(M,M,M,M) *= -1.0;
    
    LN = LN*2.0;
    
    ces.compute(LN);
    if (ces.info() != Success)
        cout << ces.info() << endl;
    
    cout << ces.eigenvalues() << endl;
    
    //
    //     Test eigenvalue condition
    //
    
    for (i=0;i<N;i++)
    {
        complex<double> lambda = ces.eigenvalues()[i];
        VectorXcd v = ces.eigenvectors().col(i);
        VectorXcd tmp = lambda * v - LN*v;
        complex<double> zero_test = tmp.dot(tmp);
        if (abs(zero_test) > 1.0e-6)
        {
            cout << "Consider the first eigenvalue, lambda = " << lambda << endl;
            cout << "i = " << i << ", eigenvalue condition = " << tmp.dot(tmp) << endl;
        }
    }
    
    //cout << "XX(0) = " << endl;
    //for (int i=0;i<N;i++)
    //{
    //    //cout << abs(ces.eigenvectors().col(39)[i]) << endl;
    //
    //}
    //cout << " " << endl;
    //
    // Test orthogonality condition
    //
    
    /*MatrixXcd ortho_test = ces.eigenvectors().adjoint()*SS.asDiagonal()*ces.eigenvectors();
    for (int L1=0;L1<N;L1++)
    {
        for (int L2=0;L2<N;L2++)
        {
            if (L1 != L2 && abs(ortho_test(L1,L2)) > 1.0E-3)
                cout << "AR and AI matrices: " << L1 << " " << L2 << " " << ces.eigenvalues()[L1] << " " << ces.eigenvalues()[L2] << " " << ortho_test(L1,L2) << "\n";
        }
    }*/
}

bool sortweights(pair<double, pair<complex<double>,VectorXcd> > a, pair<double, pair<complex<double>,VectorXcd> > b)
{
    // sortweights from largest to smallest
    return b.first < a.first;
}

void SW_Matrix::Calc_Weights()
{
    MatrixXcd XX;
    //VectorXcd tmp;
    VectorXd AL(N);
    XY.resize(N,N);
    WW.resize(N);
    XX = ces.eigenvectors().adjoint();
    
    //cout << "eigenvectors" << endl;
    //cout << ces.eigenvectors() << endl;
    

    //MatrixXcd ortho_test = XX*SS.asDiagonal()*XX.adjoint();
    MatrixXi IPR;
    int IR;
    int IFL;
    IPR.resize(N,N);
    //IPR.setZero();
    
    for(int ito=0;ito<50;ito++)
    {
        cout << "iteration # " << ito << endl;
        
        MatrixXcd ortho_test = XX*SS.asDiagonal()*XX.adjoint();
        IPR.setZero();
        IR = 1;
        for (int L1=0;L1<N;L1++)
        {
            for (int L2=0;L2<N;L2++)
            {
                if (L1 != L2 && abs(ortho_test(L1,L2)) > 1.0E-5)
                {
                    IPR(L1,L2) = 1;
                    IR = -1;
                    cout << "AR and AI matrices: " << L1 << " " << L2 << " " << ces.eigenvalues()[L1] << " " << ces.eigenvalues()[L2] << " " << ortho_test(L1,L2) << "\n";
                }
            }
        }
        if (IR != -1)
            break;
        
        for (int L1=0;L1<N;L1++)
        {
            IFL = 0;
            for (int L2=0;L2<N;L2++)
            {
                if (L2 > L1 && IPR(L1,L2) != 0 && IFL == 0)
                {
                    for(int J=0;J<N;J++)
                        {
                            XX(L1,J) = XX(L1,J) - XX(L2,J)*ortho_test(L1,L2)/ortho_test(L2,L2);
                        }
                    IFL = 1;
                }
            }
        }
    }

    
    
    
    
    for (int L1=0;L1<N;L1++)
    {
        //VectorXcd tmp1 = ces.eigenvectors().col(L1).array()*SS.array()*ces.eigenvectors().col(L1).conjugate().array();
        VectorXcd tmp1 = XX.row(L1).array()*SS.transpose().array()*XX.row(L1).conjugate().array();
        AL(L1) = tmp1.sum().real();
        //tmp = ces.eigenvectors().col(L1).array()*SS.array();
        //AL(L1) = tmp.dot(ces.eigenvectors().col(L1)).real();
        //cout << ces.eigenvalues()[L1] << "\t" << AL(L1) << endl;
        XX.row(L1) /= sqrt(abs(AL(L1)));
        //cout << "XX= " << XX.row(L1) << endl;
    }


    //
    // Reorder the XX's by the weights
    //
    
    for(int L1=0;L1<N;L1++)
    {
        pair<double, pair<complex<double>,VectorXcd> > tmp;
        tmp.first = AL(L1);
        tmp.second.first = ces.eigenvalues()[L1];
        tmp.second.second = XX.row(L1);
        eigen.push_back(tmp);
    }
    
    sort(eigen.begin(),eigen.end(),sortweights);
    
    
    for(int L1=0;L1<N;L1++)
    {
        WW(L1) = eigen[L1].second.first.real(); //eigenvalue
        //cout << "Eigenvalues= " << WW(L1) << endl;
        XY.row(L1) = eigen[L1].second.second;   //eigenvector
    }
    //cout << XY << endl;
    //cout << LN << endl << endl << XX << endl << endl;
    
    eigen.erase(eigen.begin(),eigen.end());

    //;eigen.erase (0,N);

    //cout << WW << endl;
    //cout << endl;

    /*C     Now reorder the XXs
    C
    L1=1
    L2=M+1
    DO 750 L=1,N
      if(AL(L).gt.0.) then
        DO 751 J=1,N
        XY(L1,J)=XX(L,J)
751     CONTINUE
        WW(L1)=WR(L)
C       write(6,851) L,L1,AL(L),WW(L1)
851     format(' Positive AL: ',2I5,2E15.7)
        L1=L1+1
      else
        if(L2.le.N) then
        DO 752 J=1,N
        XY(L2,J)=XX(L,J)
752     CONTINUE
        WW(L2)=WR(L)
C       write(6,852) L,L2,AL(L),WW(L2)
852     format(' Negative AL: ',2I5,2E15.7)
        endif
        L2=L2+1
      endif
750   CONTINUE*/

    
    //
    //Evalue inverse of XY or XIN
    //
    
    XIN = XY.adjoint();
    XIN.block(0,M,M,M) = -1.0*XIN.block(0,M,M,M);
    XIN.block(M,0,M,M) = -1.0*XIN.block(M,0,M,M);
    
    //cout << "XY= " << endl << XY << endl;
    //cout << "XIN= " << endl << XIN << endl ;
    
    //cout << endl;
    //cout << (XY.inverse()*XY).block(0,0,5,5)  << endl ;
    //cout << XIN*XY << endl;
    /*C
     C     Evaluate inverse of XY or XIN
     C
      DO 801 L1=1,N
      DO 802 L2=1,N
        XIN(L1,L2)=CONJG(XY(L2,L1))
        if(L1.le.M.and.L2.gt.M) then
          XIN(L1,L2)=-XIN(L1,L2)
          endif
        if(L1.gt.M.and.L2.le.M) then
          XIN(L1,L2)=-XIN(L1,L2)
        endif
802   CONTINUE
801   CONTINUE*/
    
    
}

void SW_Matrix::Calc_Intensities()
{
    complex<double> XI (0.0,1.0);
    double S = SL[0].get_sublattice()[0];
    vector<Matrix3d> V_array;
    Matrix3d V_r,V_s;
    MatrixXcd Intensities(M,3); Intensities.setZero();
    
    MatrixXcd C,SS;
    
    C.resize(2*M,2*M);
    SS.resize(3,3);
    SS.setZero();
    
    for (int i=0;i<M;i++)
    {
        V_array.push_back(SL[i].get_inv_matrix());
    }
    
    
    for(int L=0;L<2*M;L++)
    {
        for(int alpha=0;alpha<3;alpha++) //r
        {
            for( int beta=0;beta<3;beta++) // s
            {
                for(int r=0;r<M;r++) //r
                {
                    V_r = V_array[r];
                    for( int s=0;s<M;s++) // s
                    {
                        V_s = V_array[s];
                    
                        C(s,r) = 0.25*(V_r(alpha,0)+XI*V_r(alpha,1))*(V_s(beta,0)-XI*V_s(beta,1));
                        C(s,r+M) = 0.25*(V_r(alpha,0)-XI*V_r(alpha,1))*(V_s(beta,0)-XI*V_s(beta,1));
                        C(s+M,r+M) = 0.25*(V_r(alpha,0)-XI*V_r(alpha,1))*(V_s(beta,0)+XI*V_s(beta,1));
                        C(s+M,r) = 0.25*(V_r(alpha,0)+XI*V_r(alpha,1))*(V_s(beta,0)+XI*V_s(beta,1));
                    }
                }
                //cout << C << endl << endl;
                for(int r=0;r<M;r++) //r
                {
                    for( int s=0;s<M;s++) // s
                    {
                        SS(alpha,beta) -= conj(XIN(r,L))*XIN(s,L)*C(s,r);
                    }
                }
            }
        }
        std::cout << fixed;
        cout << setprecision(5) << SS*S/M*4.0 << endl << endl;
        SS.setZero();
    }




    for(int L=0;L<M;L++) //n
    {
        for(int L1=0;L1<3;L1++) //alpha
        {
            for( int L2=0;L2<M;L2++) // r
            {
                complex<double> Intensities_r = (V_array[L2](L1,0) - XI*V_array[L2](L1,1)) * XIN(L2,L+M)
                              +  (V_array[L2](L1,0) + XI*V_array[L2](L1,1)) * XIN(L2+M,L+M);
                Intensities(L,L1) += Intensities_r;
            }
        }
    }
    
    
    
    //cout << XIN << endl;
    
    /*for(int L1=0;L1<3;L1++)
        {
        for( int L2=0;L2<M;L2++)
        {
            Intensities.block(0,L1,M,1) += (V_array[L2](L1,0) - XI*V_array[L2](L1,1)) * XIN.block(L2, M, 1, M).transpose()
            +(V_array[L2](L1,0) + XI*V_array[L2](L1,1)) * XIN.block(L2+M,M,1,M).transpose();
        }
    }*/
       
    Intensities = Intensities.array().conjugate()*Intensities.array();
    Intensities *= S/(4.0*M);
    
    SXX = Intensities.col(0).real();
    SYY = Intensities.col(1).real();
    SZZ = Intensities.col(2).real();
    WP = WW.segment(M,M).array().abs();
    //cout << "WP= " << WP << endl;
}

bool evalues_equal(double a, double b)
{
    // remove eigenvalues that are equal
    double EPS = 1.0e-5;
    return fabs(a-b) < EPS;
}

void SW_Matrix::Unique_Solutions()
{
    double EPS = 1.0e-5;
    int VP_pos;
    vector<double>::iterator last;
    vector<double>::iterator uniq;
    vector<double> WP_unique;
    // Find unique eigenvalues
    for (int i=0;i<M;i++)
    {
        WP_unique.push_back(WP(i));
    }
    //reverse order
    //sort(WP_unique.begin(),WP_unique.end());
    
    /*for (int i=0;i<M;i++)
    {
        cout << WP_unique[i] << "\t" << WP[i] << endl;
    }*/
    
    WP_unique.erase(unique(WP_unique.begin(),WP_unique.end(),evalues_equal),WP_unique.end());
    // map WP_unique to VP
    NU = (int)WP_unique.size();
    VP.resize(NU);
    //Map<VectorXd> VP(&WP_unique[0],NU);
    for (int i=0;i<NU;i++)
        VP(i) = WP_unique[i];

    //cout << VP << endl;
    //TPM.resize(NU); TPM.setZero();
    //TMP.resize(NU); TMP.setZero();
    TXX.resize(NU); TXX.setZero();
    TYY.resize(NU); TYY.setZero();
    TZZ.resize(NU); TZZ.setZero();
    
    for (int i=0;i<M;i++)
    {
        VP_pos = NU; //set position to a nonsense value
        for (int j=0;j<NU;j++)
        {
            if (abs(WP(i) - VP(j)) < EPS)
            {
                VP_pos = j;
                TXX(j) += SXX(i);
                TYY(j) += SYY(i);
                TZZ(j) += SZZ(i);
                break;
            }
        }
        if (VP_pos== NU)
            cout << "error finding unique value" << endl;
    }
    //for (int i=0;i<NU;i++)
    //    cout << VP(i) << "\t" << TPM(i) << "\t" << TMP(i) << "\t" << TZZ(i) << endl;
    //TXX.array().abs();
    //TYY.array().abs();
    //TZZ.array().abs();
}

void SW_Matrix::Signif_Solutions(double KXP,double KYP,double KZP)
{
    double CXX,CYY,CZZ;
    double KX = KXP*2.0*M_PI;
    double KY = KYP*2.0*M_PI;
    double KZ = KZP*2.0*M_PI;
    double ETS = 1.0e-25;
    VectorXd SJN(NU),ZJN(NU);
    IM = 0; MI = 0;
    VI.resize(0);
    SVI.resize(0);
    
    for (int k=0;k<NU;k++)
    {
        //if (TXX(k) > ETS || TYY(k) > ETS || TZZ(k) > ETS )
        {
            VI.push_back(VP(k));
            CXX = TXX(k);
            CYY = TYY(k);
            CZZ = TZZ(k);
            cout << "CXX= " << CXX << "\t CYY= " << CYY << "\t CZZ= " << CZZ << endl;
            SVI.push_back(CXX + CYY + CZZ - (pow(KX,2)*CXX + pow(KY,2)*CYY + pow(KZ/sqrt(2.0),2)*CZZ)/(pow(KX,2)+pow(KY,2)+pow(KZ/sqrt(2.0),2)));
            MI++;
        }
    }

    VI.resize(MI);
    SVI.resize(MI);
    //cout << "Numerical Result" << endl;
    //for (int i=0;i<MI;i++)
    //{
    //    cout << KZP << "\t" << VI[i] << "\t" << SVI[i] << endl;
    //}
}

vector<double> SW_Matrix::Get_Frequencies()
{
    return VI;
}

vector<double> SW_Matrix::Get_Intensities()
{
    return SVI;
}

