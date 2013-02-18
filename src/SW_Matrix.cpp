//
//  CreateMatrix.cpp
//  Spin Wave Fit
//
//  Created by Hahn, Steven E. on 1/7/13.
//  Copyright (c) 2013 Oak Ridge National Laboratory. All rights reserved.
//

#include "SW_Matrix.h"

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

void SW_Matrix::CreateMatrix_exchange( double KXP, double KYP, double KZP)
{
    /* diagonal terms in exchange*/
    
    MatrixXi interactions;
    
    interactions << 0,1,
                    1,0;
    
    
    
    
    for (int q=0;q<M;q++)
    {
        for (int s=0;s<M;s++)
        {
            if (interactions(q,s) == 1)
            {
                
            }
                
            
        }
        for (int r=0;r<M;r++)
        {
            if (interactions(r,q) == 1)
            {
                
            }
            
        }
    }
    
    
}

void SW_Matrix::CreateMatrix_YFeO3( double KXP, double KYP, double KZP)
{
    double J,D,Kx,Kz,z,S,theta;
    double eta_q,gamma_q;
    double A_q,B_q,C_q,D_q;
    double KX = KXP*2.0*PI;
    double KY = KYP*2.0*PI;
    double KZ = KZP*2.0*PI;
    
    J = X[0];
    D = X[1];
    Kx = X[2];
    Kz = X[3];
    z=6.0;
    S=SL[0].get_sublattice()[0];
    cout << S << endl;
    
    theta = SL[0].get_sublattice()[1]; //PI/2.0 - 0.01098;

    cout << theta << endl;
 
    gamma_q = (cos(KZ/2.0) + cos((KX+KY)/2.0) + cos((KX-KY)/2.0))/3.0;
    eta_q =   (cos(KZ/2.0) + cos((KX+KY)/2.0) + cos((KX-KY)/2.0))/3.0;
    
    A_q = -2.0*z*J*S*cos(2.0*theta)+z*D*S*sin(2.0*theta) - Kx*S*(-2.0+3.0*pow(cos(theta),2))
        -Kz*S*(-2.0+3.0*pow(sin(theta),2));
    
    B_q = -0.5*Kx*S*pow(cos(theta),2) -0.5*Kz*S*pow(sin(theta),2);
    
    C_q = 0.5*z*D*S*sin(2.0*theta)*eta_q - z*J*S*(cos(2.0*theta)+1.0)*gamma_q;
    
    D_q = 0.5*z*D*S*sin(2.0*theta)*eta_q - z*J*S*(cos(2.0*theta)-1.0)*gamma_q;
    
    LN(0,0) = A_q;
    LN(0,1) = C_q;
    LN(0,2) = 2.0*B_q;
    LN(0,3) = D_q;
    LN(1,0) = C_q;
    LN(1,1) = A_q;
    LN(1,2) = D_q;
    LN(1,3) = 2.0*B_q;
    
    LN.block(M,0,M,M) = -1.0*LN.block(0,M,M,M);
    LN.block(M,M,M,M) = -1.0*LN.block(0,0,M,M);
    
    //cout << LN << endl;
    
    //cout << sqrt(24.0*J*S*2.0*S*(Kx-Kz)) << endl;
    //cout << sqrt(24.0*J*S*2.0*Kx*S) << endl;
}

void SW_Matrix::CreateMatrix_AFM(double KXP, double KYP, double KZP)
{
    double J,D,z,S;
    double gamma_q;
    double A_q,C_q;
    double KX = KXP*2.0*PI;
    double KY = KYP*2.0*PI;
    double KZ = KZP*2.0*PI;
    LN.resize(N,N); LN.setZero();
    SS.resize(N);  SS.setZero();
    
    J = X[0];
    D = X[1];
    z=6.0;
    S=5.0/2.0;
    
    for (int j=0;j<M;j++)
    {
        SS(j) = 1.0;
        SS(j+M) = -1.0;
    }
    
    gamma_q = (cos(KX) + cos(KY) + cos(KZ))/3.0;
    
    A_q = -2.0*z*J*S + D*S;
    
    C_q = 2.0*z*J*S*gamma_q;
    
    LN(0,0) = A_q;
    LN(0,3) = C_q;
    LN(1,1) = A_q;
    LN(1,2) = C_q;
    
    LN.block(M,0,M,M) = -1.0*LN.block(0,M,M,M);
    LN.block(M,M,M,M) = -1.0*LN.block(0,0,M,M);
    
    LN *= 1.0;
    
    cout << sqrt(pow(2.0*J*z*S,2)*(1.0-pow(gamma_q,2))+4*abs(J)*z*S*S*D+S*S*D*D) << endl;
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
    
    ces.compute(LN);
    if (ces.info() != Success)
        cout << ces.info() << endl;

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
    //    cout << abs(ces.eigenvectors().col(39)[i]) << endl;
    //}
    //cout << " " << endl;
    //
    // Test orthogonality condition
    //
    
    MatrixXcd ortho_test = ces.eigenvectors().adjoint()*SS.asDiagonal()*ces.eigenvectors();
    for (int L1=0;L1<N;L1++)
    {
        for (int L2=0;L2<N;L2++)
        {
            if (L1 != L2 && abs(ortho_test(L1,L2)) > 1.0E-3)
                cout << "AR and AI matrices: " << L1 << " " << L2 << " " << ces.eigenvalues()[L1] << " " << ces.eigenvalues()[L2] << " " << ortho_test(L1,L2) << "\n";
        }
    }

    //cout << ces.eigenvalues() << endl;
    
    //
    //     Normalize the XXs
    //
}

bool sortweights(pair<double, pair<complex<double>,VectorXcd> >& a, pair<double, pair<complex<double>,VectorXcd> >& b)
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

    for (int L1=0;L1<N;L1++)
    {
        //VectorXcd tmp1 = ces.eigenvectors().col(L1).array()*SS.array()*ces.eigenvectors().col(L1).conjugate().array();
        VectorXcd tmp1 = XX.row(L1).array()*SS.transpose().array()*XX.row(L1).conjugate().array();
        AL(L1) = tmp1.sum().real();
        //tmp = ces.eigenvectors().col(L1).array()*SS.array();
        //AL(L1) = tmp.dot(ces.eigenvectors().col(L1)).real();
        //cout << ces.eigenvalues()[L1] << "\t" << AL(L1) << endl;
        XX.row(L1) /= sqrt(abs(AL(L1)));
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

    //cout << (XIN*XY).diagonal() << endl;
    //cout << (XIN*XY).block(0,0,5,5) << endl ;
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

void SW_Matrix::Rotation_Matrix()
{
    complex<double> XI (0.0,1.0);
    double S = SL[0].get_sublattice()[0];
    vector<Matrix3d> V_array;
    Matrix3d V1,V2;
    MatrixXcd Intensities(M,3); Intensities.setZero();
    
    V_array.push_back(SL[0].get_inv_matrix());
    V_array.push_back(SL[1].get_inv_matrix());
   
    /*for(int L=0;L<M;L++)
    {
        for(int L1=0;L1<3;L1++)
        {
            for( int L2=0;L2<M;L2++)
            {
                Intensities(L,L1) += (V_array[L2](L1,0) - XI*V_array[L2](L1,1)) * XIN(L2,L+M)
                                    +(V_array[L2](L1,0) + XI*V_array[L2](L1,1)) * XIN(L2+M,L+M);
            }
        }
    }*/
    
    for(int L1=0;L1<3;L1++)
        {
        for( int L2=0;L2<M;L2++)
        {
            Intensities.block(0,L1,M,1) += (V_array[L2](L1,0) - XI*V_array[L2](L1,1)) * XIN.block(L2, M, 1, M).transpose()
            +(V_array[L2](L1,0) + XI*V_array[L2](L1,1)) * XIN.block(L2+M,M,1,M).transpose();
        }
    }
       
    Intensities = Intensities.array().conjugate()*Intensities.array();
    Intensities *= S/8.0;
    
    SXX = Intensities.col(0).real();
    SYY = Intensities.col(1).real();
    SZZ = Intensities.col(2).real();
    WP = WW.segment(0,M).array().abs();
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
    double KX = KXP*2.0*PI;
    double KY = KYP*2.0*PI;
    double KZ = KZP*2.0*PI;
    double ETS = 1.0e-2;
    VectorXd SJN(NU),ZJN(NU);
    IM = 0; MI = 0;
    VI.resize(NU); VI.setZero();
    SVI.resize(NU); SVI.setZero();
    /*for (int k=0;k<NU;k++)
    {
        SJN(k) = (TPM(k) + TMP(k))/2.0;
        ZJN(k) = TZZ(k);
        if (SJN(k) > ETS || ZJN(k) > ETS )
        {
            VI(IM) = VP(k);
            CXX = ZJN(k);
            CYY = SJN(k) - ZJN(k);
            CZZ = ZJN(k);
            SVI(IM) = CXX + CYY + CZZ - (pow(KZ,2)*CXX + pow(KX,2)*CYY+pow(KY,2)*CZZ)/(pow(KX,2)+pow(KY,2)+pow(KZ,2));
            if (abs(SVI(IM)) < 1.0e-6)
            {
                SVI(IM) = 0.0;
            }
            IM++;
            MI++;
        }
    }*/
    for (int k=0; k<NU;k++)
    {
        if (TXX(k) > ETS || TYY(k) > ETS || TZZ(k) > ETS )
        {
            VI(k) = VP(k);
            CXX = TXX(k);
            CYY = TYY(k);
            CZZ = TZZ(k);
            SVI(IM) = CXX + CYY + CZZ - (pow(KZ,2)*CXX + pow(KX,2)*CYY+pow(KY,2)*CZZ)/(pow(KX,2)+pow(KY,2)+pow(KZ,2));
            IM++;
            MI++;
        }
    }
    
    
    //VI.resize(MI);
    //SVI.resize(MI);
    for (int i=0;i<MI;i++)
    {
        cout << KY/2.0/PI << "\t" << VI(i) << "\t" << SVI(i) << endl;
    }
}