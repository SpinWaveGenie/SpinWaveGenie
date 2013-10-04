#include "SW_Builder.h"
#include <iostream>

using namespace std;

SW_Builder::SW_Builder(boost::shared_ptr<Cell>& cell_in)
{
    cell = cell_in;
    SpinWave temp(cell);
    SW = temp;
}

void SW_Builder::Add_Interaction(Interaction * in)
{
    in->calcConstantValues(this->cell);
    interactions.push_back(in);
    interactions.sort();

}

Eigen::VectorXcd SW_Builder::checkFirstOrderTerms()
{
    int M=0;
    for (SublatticeIterator sl=cell->begin(); sl!=cell->end(); ++sl)
    {
        M++;
    }
    Eigen::VectorXcd firstOrder;
    firstOrder.setZero(2*M);
    boost::ptr_vector<Interaction>::iterator iter;
    for (iter = interactions.begin(); iter != interactions.end(); iter++)
    {
        iter->checkFirstOrderTerms(this->cell,firstOrder);
    }
    return firstOrder;
}


SpinWave SW_Builder::Create_Element(double KX,double KY,double KZ)
{
    Eigen::Vector3d K;
    Eigen::Matrix3d recip;
    recip = cell->getReciprocalVectors();
    K << KX,KY,KZ;
    //cout << "K before " << K.transpose() << endl;
    K = recip*K;
    //cout << "K after " << K.transpose() << endl;
    SW.Set_Kpoint(K[0],K[1],K[2]);
    SW.Clear_Matrix();
    boost::ptr_vector<Interaction>::iterator iter;
    int quad = 0;
    for (iter = interactions.begin(); iter != interactions.end(); iter++)
    {
        iter->calcChangingValues(SW.cell,K);
        iter->Update_Matrix(K,SW.cell,SW.LN,quad);
        iter->Update_Matrix(K,SW.cell,SW.LN,quad+1);
        iter->Update_Matrix(K,SW.cell,SW.LN,quad+2);
        iter->Update_Matrix(K,SW.cell,SW.LN,quad+3);
    }
    
    //cout << "LN" << endl;
    //cout << SW.LN << endl;
    return SW;
}