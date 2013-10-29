#include "SpinWave.h"
#include <Eigen/Cholesky>
#include <iomanip>
#include <unordered_map>


using namespace Eigen;
using namespace std;


double form_factor(string type, double x,double y,double z)
{
    Vector4d big_F,little_f,eval_points,weights;
    std::map<string,std::vector<double> > coeff;

    std::vector<double> tmp = {0.408600, 28.810900,  0.607700,  8.543700, -0.029500,  0.276800,  0.012300};
    coeff.emplace("V0", tmp);
    /*tmp = {0.444400, 32.647900,  0.568300,  9.097100, -0.228500,  0.021800,  0.215000};
    coeff.emplace("V1",tmp);
    tmp = {0.408500, 23.852600,  0.609100,  8.245600, -0.167600,  0.041500,  0.149600};
    coeff.emplace("V2",tmp);
    tmp = {0.359800, 19.336399,  0.663200,  7.617200, -0.306400,  0.029600,  0.283500};
    coeff.emplace("V3",tmp);
    tmp = {0.310600, 16.816000,  0.719800,  7.048700, -0.052100,  0.302000,  0.022100};
    coeff.emplace("V4",tmp);
    tmp = {0.243800, 24.962900,  0.147200, 15.672800,  0.618900,  6.540300, -0.010500};
    coeff.emplace("MN0",tmp);
    tmp = {-0.013800,  0.421300,  0.423100, 24.667999,  0.590500,  6.654500, -0.001000};
    coeff.emplace("MN1",tmp);
    tmp = {0.422000, 17.684000,  0.594800,  6.005000,  0.004300, -0.609000, -0.021900};
    coeff.emplace("MN2",tmp);
    tmp = {0.419800, 14.282900,  0.605400,  5.468900,  0.924100, -0.008800, -0.949800};
    coeff.emplace("MN3",tmp);
    tmp = {0.376000, 12.566100,  0.660200,  5.132900, -0.037200,  0.563000,  0.001100};
    coeff.emplace("MN4",tmp);*/
    //coeff[FE0]   0.070600, 35.008499,  0.358900, 15.358300,  0.581900,  5.560600, -0.011400
    //coeff[FE1]   0.125100, 34.963299,  0.362900, 15.514400,  0.522300,  5.591400, -0.010500
    //coeff[FE2]   0.026300, 34.959702,  0.366800, 15.943500,  0.618800,  5.593500, -0.011900
    tmp = {0.397200, 13.244200,  0.629500,  4.903400, -0.031400,  0.349600,  0.004400};
    coeff.emplace("FE3",tmp);
    //coeff[FE4]   0.378200, 11.380000,  0.655600,  4.592000, -0.034600,  0.483300,  0.000500
    
    //big_F << 0.3972,0.6295,-0.0314,0.0044;
    //little_f << 13.2442,4.9034,0.3496,0.0;
    
    vector<double> F = coeff[type];
    double s = sqrt(pow(x,2) + pow(y,2) + pow(z,2))/(4.0*M_PI);
    double f_Q = F[6];
    for(int k=0;k<3;k++)
    {
        f_Q += F[2*k]*exp(-1.0*F[2*k+1]*pow(s,2));
    }
    return f_Q;
}

SpinWave::SpinWave()
{
    
}

SpinWave::SpinWave(Cell& cell_in)
{

    cell = cell_in;
    M = cell.size();
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
    VI.clear();
    SVI.clear();
}

void SpinWave::Set_Kpoint(double KX, double KY, double KZ)
{
    KXP = KX;
    KYP = KY;
    KZP = KZ;
}

void SpinWave::Calc_Eigenvalues()
{
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
    
    /*cout << KXP << KYP << KZP << endl;
    LLT<MatrixXcd> lltOfA(LN);
    MatrixXcd H = lltOfA.matrixL();
    MatrixXcd Ht = H.transpose();
    
    MatrixXcd Ltest = H*LN*Ht;
    
    ces.compute(Ltest);
    if (ces.info() != Success)
        cout << ces.info() << endl;
    
    cout << ces.eigenvalues().transpose() << endl;
    */

    //cout << "LN= " << endl;
    //cout << fixed << setprecision(4) << LN << endl;
    
    ces.compute(LN);
    if (ces.info() != Success)
        cout << ces.info() << endl;
    
    //cout << ces.eigenvalues().transpose() << endl << endl;
    
    //
    //     Test eigenvalue condition
    //
    
    for (int i=0;i<N;i++)
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
    MatrixXcdRowMajor XX;
    //VectorXcd tmp;
    //XY.resize(N,N);
    WW.resize(N);
    XX = ces.eigenvectors().adjoint();
    
    //cout << "eigenvectors" << endl;
    //cout << ces.eigenvectors() << endl;
    
    
    //MatrixXcd ortho_test = XX*SS.asDiagonal()*XX.adjoint();
    MatrixXi IPR(N,N);
    VectorXcd TEST(N);
    int IR;
    int IFL;
    //IPR.setZero();
    int maxIterations = 50;
    for(int ito=0;ito<maxIterations;ito++)
    {
        //cout << "ito= " << ito << endl;
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
                    //cout << "AR and AI matrices: " << L1 << " " << L2 << " " << ces.eigenvalues()[L1] << " " << ces.eigenvalues()[L2] << " " << ortho_test(L1,L2) << "\n";
                }
            }
        }
        //cout << ortho_test << endl;
        //cout << "IR= " << IR << endl;
        if (IR != -1)
        {
            TEST = ortho_test.diagonal();
            break;
        }
        for (int L1=0;L1<N;L1++)
        {
            IFL = 0;
            for (int L2=0;L2<N;L2++)
            {
                if (L2 > L1 && IPR(L1,L2) != 0 && IFL == 0)
                {
                    for(int J=0;J<N;J++)
                    {
                        XX(L1,J) -= XX(L2,J)*ortho_test(L1,L2)/ortho_test(L2,L2);
                    }
                    IFL = 1;
                }
            }
        }
        if (ito==maxIterations-1)
        {
            TEST = ortho_test.diagonal();
            cout << "Error calculating frequencies" << endl;
        }
    }
    
    vector<results> AL(N);
    for (int L1=0;L1<N;L1++)
    {
        AL[L1].weight = TEST[L1].real();
        AL[L1].index = L1;
        XX.row(L1) /= sqrt(abs(AL[L1].weight));
    }
    
    //
    // Reorder the XX's by the weights
    //
    
    sort(AL.begin(),AL.end());

    //cout << "AL.index= " << endl;
    for(int L1=0;L1<N;L1++)
    {
        WW(L1) = ces.eigenvalues()[AL[L1].index].real(); //eigenvalue
        //cout << AL[L1].index << " ";

    }
    //cout << endl;
    
    
    //Swap rows to reflect ordering of eigenvalues.
    //The swap moves row L1 to a new position and the index must be
    //updated to reflect this.
    int old_index;
    for(int L1=0;L1<N;L1++)
    {
        for(int L2=L1;L2<N;L2++)
        {
            if( L1 == AL[L2].index)
            {
                old_index = L2;
                break;
            }
        }
        XX.row(L1).swap(XX.row(AL[L1].index));   //eigenvector
        AL[old_index].index = AL[L1].index;
        //AL[L1].index = L1;
    }
    

    /*cout <<"Eigenvalues: " << endl ;
    cout << ces.eigenvalues().transpose() << endl;
    cout << "Weights: " << endl;
    cout << TEST.transpose() << endl;
    cout <<"Sorted Eigenvalues: " << endl;
    cout << WW.transpose() << endl;
    cout << "Eigenvectors: " << endl;
    cout << ces.eigenvectors() << endl;
    cout << "Normalized eigenvectors:" << endl;
    cout << XX << endl;
    */
    //eigenresults.erase(eigenresults.begin(),eigenresults.end());
    
    //
    //Evalue inverse of XY or XIN
    //
    
    XIN = XX.adjoint();
    XIN.block(0,M,M,M) *= -1.0;//*XIN.block(0,M,M,M);
    XIN.block(M,0,M,M) *= -1.0;//*XIN.block(M,0,M,M);
    
    //cout << "XY= " << endl << XY << endl;
    //cout << "XIN= " << endl << XIN << endl ;
    
    //cout << endl;
    //cout << (XY.inverse()*XY).block(0,0,5,5)  << endl ;
    //cout << XIN*XY << endl;
}

void SpinWave::Calc_Intensities()
{
    complex<double> XI (0.0,1.0);
    double KX = KXP;
    double KY = KYP;
    double KZ = KZP;
    //vector<Matrix3d> V_array;
    Matrix3d V_r;//,V_s;
    double S_r,formFactor;
    ArrayXXcd Intensities(M,3); Intensities.setZero();
    VectorXd SXX,SYY,SZZ;
    MatrixXcd C,SS;
    
    C.resize(2*M,2*M);
    //SS.resize(3,3);
    //SS.setZero();

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
    
    // sl = cell->begin();
    //double S = ;

    long L2 = 0;
    for (SublatticeIterator sl = cell.begin(); sl!=cell.end();++sl) //r
    {
        V_r = (*sl->getInverseMatrix());
        S_r = sl->getMoment();
        string type = sl->getType();
        formFactor = form_factor(type,KX,KY,KZ);
        for(int L=0;L<M;L++) //n
        {
            for(int L1=0;L1<3;L1++) //alpha
            {
                complex<double> Intensities_r = (V_r(L1,0) - XI*V_r(L1,1)) * XIN(L2,L+M)
                +  (V_r(L1,0) + XI*V_r(L1,1)) * XIN(L2+M,L+M);
                Intensities(L,L1) += sqrt(S_r)*formFactor*Intensities_r;
            }
        }
        L2++;
    }
    
    Intensities *= Intensities.conjugate();
    Intensities *= 1.0/(4.0*M);
    
    SXX = Intensities.col(0).real();
    SYY = Intensities.col(1).real();
    SZZ = Intensities.col(2).real();
    
    VI.reserve(M);
    SVI.reserve(M);
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
    //VI.resize(0);
    //SVI.resize(0);
    
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
    this->Calc_Eigenvalues();
    this->Calc_Weights();
    this->Calc_Intensities();
    //this->Unique_Solutions();
    //this->Signif_Solutions();
}

vector<double> SpinWave::Get_Frequencies()
{
    return VI;
}

vector<double> SpinWave::Get_Intensities()
{
    return SVI;
}
