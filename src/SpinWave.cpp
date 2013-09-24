#include "SpinWave.h"

using namespace Eigen;
using namespace std;

SpinWave::SpinWave()
{
    
}

SpinWave::SpinWave(boost::shared_ptr<Cell>& cell_in)
{

    cell = cell_in;
    M = (int) cell->size();
    N = 2*M;
    
    LN.resize(N,N); LN.setZero();
    
    SS.resize(N);
    SS.setZero();
    
    for (int j=0;j<M;j++)
    {
        SS(j) = 1.0;
    }
    
    for (int j=M;j<2*M;j++)
    {
        SS(j) = -1.0;
    }
}

void SpinWave::Clear_Matrix()
{
    LN.setZero();
}

void SpinWave::Set_Kpoint(double KX, double KY, double KZ)
{
    KXP = KX;
    KYP = KY;
    KZP = KZ;
}

void SpinWave::Calc_Eigenvalues()
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
    
    
    //cout << LN << endl;
    
    ces.compute(LN);
    if (ces.info() != Success)
        cout << ces.info() << endl;
    
    //cout << ces.eigenvalues() << endl;
    
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

void SpinWave::Calc_Weights()
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
        //cout << "iteration # " << ito << endl;
        
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
    
    
    VectorXcd tmp1(N);
    for (int L1=0;L1<N;L1++)
    {
        //VectorXcd tmp1 = ces.eigenvectors().col(L1).array()*SS.array()*ces.eigenvectors().col(L1).conjugate().array();
        //cout << XX.row(L1).array() << endl;
        //cout << XX.row(L1).array()*SS.transpose().array()*XX.row(L1).conjugate().array() << endl;

        tmp1 = XX.row(L1).array()*SS.transpose().array()*XX.row(L1).conjugate().array();
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
        results tmp;
        tmp.weight = AL(L1);
        tmp.eigenvalue = ces.eigenvalues()[L1];
        tmp.eigenvector = XX.row(L1);
        eigenresults.push_back(tmp);
    }
    
    /*for(int L1=0;L1<N;L1++)
     {
     pair<double, pair<complex<double>,VectorXcd> > tmp;
     tmp.first = AL(L1);
     tmp.second.first = ces.eigenvalues()[L1];
     tmp.second.second = XX.row(L1);
     eigen.push_back(tmp);
     }*/
    
    sort(eigenresults.begin(),eigenresults.end());
    
    for(int L1=0;L1<N;L1++)
    {
        WW(L1) = eigenresults[L1].eigenvalue.real(); //eigenvalue
        //cout << "Eigenvalues= " << WW(L1) << endl;
        XY.row(L1) = eigenresults[L1].eigenvector;   //eigenvector
    }
    
    //cout << XY << endl;
    //cout << LN << endl << endl << XX << endl << endl;
    
    eigenresults.erase(eigenresults.begin(),eigenresults.end());
    
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

void SpinWave::Calc_Intensities()
{
    complex<double> XI (0.0,1.0);
    double KX = KXP;//*2.0*M_PI;
    double KY = KYP;//*2.0*M_PI;
    double KZ = KZP;//*2.0*M_PI;
    vector<Matrix3d> V_array;
    Matrix3d V_r,V_s;
    MatrixXcd Intensities(M,3); Intensities.setZero();
    VectorXd SXX,SYY,SZZ;
    MatrixXcd C,SS;
    
    C.resize(2*M,2*M);
    SS.resize(3,3);
    SS.setZero();
    

    
    //for(sl.First();!sl.IsDone();sl.Next())
    //{
    //    V_array.push_back(sl.CurrentItem()->get_inv_matrix());
    //}
        
    //typedef std::map<> >::iterator it_type;
    //for(it_type iterator = m.begin(); iterator != m.end(); iterator++) {
    
    /*for(int L=0;L<2*M;L++)
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
     for(int r=0;r<2*M;r++) //r
     {
     for( int s=0;s<2*M;s++) // s
     {
     SS(alpha,beta) -= conj(XIN(r,L))*XIN(s,L)*C(s,r);
     }
     }
     }
     }
     std::cout << fixed;
     cout << setprecision(5) << -1.0*SS*S/M << endl << endl;
     SS.setZero();
     }*/
    
    //CellIter sl(cell);
    //double S = sl.CurrentItem()->getMoment()[0];
    
    SublatticeIterator sl = cell->begin();
    double S = (*(*sl)->getMoment())[0];


    for(int L=0;L<M;L++) //n
    {
        for(int L1=0;L1<3;L1++) //alpha
        {
            //sl.First();
            sl = cell->begin();
            for( int L2=0;L2<M;L2++) // r
            {
                V_r = (*sl)->getInverseMatrix();
                complex<double> Intensities_r = (V_r(L1,0) - XI*V_r(L1,1)) * XIN(L2,L+M)
                +  (V_r(L1,0) + XI*V_r(L1,1)) * XIN(L2+M,L+M);
                Intensities(L,L1) += Intensities_r;
                ++sl;
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
    
    VI.clear();
    SVI.clear();
    
    for (int i=0;i<M;i++)
    {
        VI.push_back(abs(WW[i+M]));
        //cout << "SXX= " << SXX[i] << "\t SYY= " << SYY[i] << "\t SZZ= " << SZZ[i] << endl;
        SVI.push_back(SXX(i) + SYY(i) + SZZ(i) - (pow(KX,2)*SXX(i) + pow(KY,2)*SYY(i) + pow(KZ,2)*SZZ(i))/(pow(KX,2)+pow(KY,2)+pow(KZ,2)));
    }
    
    //WP = WW.segment(M,M).array().abs();
    //cout << "WP= " << WP << endl;
}

bool evalues_equal(double a, double b)
{
    // remove eigenvalues that are equal
    double EPS = 1.0e-5;
    return fabs(a-b) < EPS;
}

void SpinWave::Unique_Solutions()
{
    double EPS = 1.0e-5;
    int VP_pos;
    vector<double>::iterator last;
    vector<double>::iterator uniq;
    vector<double> VI_unique,SVI_unique;
    // Find unique eigenvalues
    VI_unique = VI;
    
    VI_unique.erase(unique(VI_unique.begin(),VI_unique.end(),evalues_equal),VI_unique.end());
    
    NU = (int)VI_unique.size();
    VI_unique.resize(NU);
    
    for (int i=0;i<NU;i++)
    {
        SVI_unique.push_back(0.0);
    }
    
    for (int i=0;i<M;i++)
    {
        VP_pos = NU; //set position to a nonsense value
        for (int j=0;j<NU;j++)
        {
            if (abs(VI[i] - VI_unique[j]) < EPS)
            {
                VP_pos = j;
                SVI_unique[j] += SVI[i];
                break;
            }
        }
        if (VP_pos== NU)
            cout << "error finding unique value" << endl;
    }
    
    VI = VI_unique;
    SVI = SVI_unique;
    
    //for (int i=0;i<NU;i++)
    //    cout << VP(i) << "\t" << TPM(i) << "\t" << TMP(i) << "\t" << TZZ(i) << endl;
    //TXX.array().abs();
    //TYY.array().abs();
    //TZZ.array().abs();
}

void SpinWave::Signif_Solutions()
{
    double ETS = 1.0e-5;
    vector<double> VI_signif,SVI_signif;
    IM = 0; MI = 0;
    VI.resize(0);
    SVI.resize(0);
    
    for (int k=0;k<NU;k++)
    {
        if (SVI[k] > ETS )
        {
            VI_signif.push_back(VI[k]);
            SVI_signif.push_back(SVI[k]);
        }
    }
    
    VI = VI_signif;
    SVI = SVI_signif;
    
    //SVI.resize(MI);
    //cout << "Numerical Result" << endl;
    //for (int i=0;i<MI;i++)
    //{
    //    cout << KZP << "\t" << VI[i] << "\t" << SVI[i] << endl;
    //}
}

void SpinWave::Calc()
{
    Calc_Eigenvalues();
    Calc_Weights();
    Calc_Intensities();
    //Unique_Solutions();
    //Signif_Solutions();
}

vector<double> SpinWave::Get_Frequencies()
{
    return VI;
}

vector<double> SpinWave::Get_Intensities()
{
    return SVI;
}

/*VectorXd SpinWave::Get_Evector(double E_min, double E_max, double E_points)
{
    
}*/
