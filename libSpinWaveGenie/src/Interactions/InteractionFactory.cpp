//
//  InteractionFactory.cpp
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 3/27/14.
//
//

#include "SpinWaveGenie/Interactions/InteractionFactory.h"
#include "SpinWaveGenie/Interactions/ExchangeInteraction.h"
#include "SpinWaveGenie/Interactions/ExchangeInteractionSameSublattice.h"
#include "SpinWaveGenie/Interactions/DM_Y_Interaction.h"
#include "SpinWaveGenie/Interactions/DM_Z_Interaction.h"
#include "SpinWaveGenie/Interactions/AnisotropyInteraction.h"
#include "SpinWaveGenie/Interactions/MagneticFieldInteraction.h"

namespace SpinWaveGenie
{

std::unique_ptr<Interaction> InteractionFactory::getExchange(std::string name, double value, std::string sl_r, std::string sl_s, double min, double max)
{
    if (sl_r.compare(sl_s) == 0)
        return std::move(std::unique_ptr<Interaction>(new ExchangeInteractionSameSublattice(name,value,sl_r,min,max)));
    else
        return std::move(std::unique_ptr<Interaction>(new ExchangeInteraction(name,value,sl_r,sl_s,min,max)));
}

std::unique_ptr<Interaction> InteractionFactory::getDzyaloshinskiiMoriya(std::string name, double value, Vector3 unitVector, std::string sl_r, std::string sl_s, double min, double max)
{
    double x,y,z;
    x = unitVector[0];
    y = unitVector[1];
    z = unitVector[2];
    
    if (x < 0.01 && y > 0.99 && z < 0.01)
    {
        return std::move(std::unique_ptr<Interaction>(new DM_Y_Interaction(name,value,sl_r,sl_s,min,max)));
    }
    else if (x < 0.01 && y < 0.01 && z > 0.99)
    {
        return std::move(std::unique_ptr<Interaction>( new DM_Z_Interaction(name,value,sl_r,sl_s,min,max)));
    }
    else
    {
        // a general DM interaction has not yet been implemented
        return std::move(std::unique_ptr<Interaction>( new DM_Z_Interaction(name,value,sl_r,sl_s,min,max)));
    }
}

std::unique_ptr<Interaction> InteractionFactory::getAnisotropy(std::string name, double value, Vector3 unitVector, std::string sl_r)
{
    return std::move(std::unique_ptr<Interaction>(new AnisotropyInteraction(name,value,unitVector,sl_r)));
}

std::unique_ptr<Interaction> InteractionFactory::getMagneticField(std::string name, double value, Vector3 unitVector, std::string sl_r)
{
    return std::move(std::unique_ptr<Interaction>(new MagneticFieldInteraction(name,value,unitVector,sl_r)));
}
    
}



