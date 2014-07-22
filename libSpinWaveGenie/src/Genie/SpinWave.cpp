
#include <Eigen/Cholesky>
#include <iomanip>
#include <cmath>
#include <random>

#include "SpinWaveGenie/Genie/SpinWave.h"
#include "SpinWaveGenie/Genie/Neighbors.h"
#include "SpinWaveGenie/Containers/Matrices.h"

using namespace Eigen;
using namespace std;

namespace SpinWaveGenie
{


SpinWave::SpinWave()
{
    
}

SpinWave::SpinWave(Cell& cell_in, boost::ptr_vector<Interaction> interactions_in)
{
    cell = cell_in;
    interactions = interactions_in;
    M = cell.size();
    N = 2*M;
    
    LN.setZero(N,N);
    
    SS.setZero(N);
    
    for (int j=0;j<M;j++)
    {
        SS(j) = 1.0;
    }
    
    for (int j=M;j<2*M;j++)
    {
        SS(j) = -1.0;
    }
}



void SpinWave::createMatrix(double KX,double KY,double KZ)
{
    Vector3 K;
    Matrix3 recip;
    recip = cell.getReciprocalVectors();
    K << KX,KY,KZ;
    K = K.transpose()*recip;
    this->KXP = KX;
    this->KYP = KY;
    this->KZP = KZ;
    clearMatrix();
    boost::ptr_vector<Interaction>::iterator iter;
    for (iter = interactions.begin(); iter != interactions.end(); iter++)
    {
        iter->updateMatrix(K,LN);
    }
}

void SpinWave::clearMatrix()
{
    LN.setZero();
    VI.clear();
}

const Cell& SpinWave::getCell() const
{
    return cell;
}

void SpinWave::calculateEigenvalues()
{
    //cout << LN << endl << endl;
    
    LN.block(0,M,M,M) *= -1.0;
    LN.block(M,M,M,M) *= -1.0;
    
    LN = LN*2.0;
    
    for(int i=0;i<2*M;i++)
    {
        for( int j=0;j<2*M;j++)
        {
            if (std::abs(LN(i,j)) < 1.0e-10)
            {
                LN(i,j) = complex<double>(0.0,0.0);
            }
        }
    }

    ces.compute(LN);
    if (ces.info() != Success)
        cout << ces.info() << endl;
    
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
}

bool IsPositive (results i) { return i.weight > 0; }

void SpinWave::calculateWeights()
{
    MatrixXcdRowMajor XX;
    WW.resize(N);
    XX = ces.eigenvectors().adjoint();
    
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
    
    
    
    partition(AL.begin(), AL.end(),IsPositive);

    for(int L1=0;L1<N;L1++)
    {
        WW(L1) = ces.eigenvalues()[AL[L1].index].real(); //eigenvalue
    }
    
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
    }
   
    //
    //Evalue inverse of XY or XIN
    //
    
    XIN = XX.adjoint();
    XIN.block(0,M,M,M) *= -1.0;
    XIN.block(M,0,M,M) *= -1.0;
    
}

void SpinWave::calculateIntensities()
{
    complex<double> XI (0.0,1.0);
    double KX = KXP;
    double KY = KYP;
    double KZ = KZP;
    Matrix3 V_r;//,V_s;
    double S_r,ff;
    ArrayXXcd Intensities(M,3); Intensities.setZero();
    VectorXd SXX,SYY,SZZ;
    long L2 = 0;
    for (Cell::Iterator sl = cell.begin(); sl!=cell.end();++sl) //r
    {
        V_r = sl->getInverseMatrix();
        S_r = sl->getMoment();
        formFactor.setType(sl->getType());
        ff = formFactor.getFormFactor(KX,KY,KZ);
        /*if (sl->getType() == "MN2")
        {
            formFactor.setType("CO2");
            ff = 0.6*ff+ 0.4*formFactor.getFormFactor(KX,KY,KZ);
        }*/
        /*if (sl->getType() == "CO2")
        {
            formFactor.setType("MN2");
            ff = 0.8*ff+ 0.2*formFactor.getFormFactor(KX,KY,KZ);
        }*/
        for(int L=0;L<M;L++) //n
        {
            for(int L1=0;L1<3;L1++) //alpha
            {
                complex<double> Intensities_r = (V_r(L1,0) - XI*V_r(L1,1)) * XIN(L2,L+M)
                +  (V_r(L1,0) + XI*V_r(L1,1)) * XIN(L2+M,L+M);
                Intensities(L,L1) += sqrt(S_r)*ff*Intensities_r;
            }
        }
        L2++;
    }
    
    Intensities *= Intensities.conjugate();
    Intensities *= 1.0/(4.0*M);
    
    SXX = Intensities.col(0).real();
    SYY = Intensities.col(1).real();
    SZZ = Intensities.col(2).real();
    
    for (int i=0;i<M;i++)
    {
        Point pt;
        pt.frequency = abs(WW[i+M]);
        //cout << SXX(i) << " " << SYY(i) << " " << SZZ(i) << endl;
        pt.intensity = SXX(i) + SYY(i) + SZZ(i) - (pow(KX,2)*SXX(i) + pow(KY,2)*SYY(i) + pow(KZ,2)*SZZ(i))/(pow(KX,2)+pow(KY,2)+pow(KZ,2));
        VI.insert(pt);
    }
}

void SpinWave::calculate()
{
    this->calculateEigenvalues();
    this->calculateWeights();
    this->calculateIntensities();
}

Results SpinWave::getPoints()
{
    return VI;
}

}
