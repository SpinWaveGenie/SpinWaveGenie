#include "CommonFunctions.h"

using namespace std;
using namespace SpinWaveGenie;

void alongKAxis()
{
    SpinWave SW = createModel();
    
    PointsAlongLine Line;
    Line.setFirstPoint(2.0,-1.5,-3.0);
    Line.setFinalPoint(2.0, 1.5,-3.0);
    Line.setNumberPoints(201);
    ThreeVectors<double> kPoints = Line.getPoints();
    
    Energies energies(0.0, 80.0, 201);
    
    TwoDimGaussian resinfo;
    resinfo.a = 1109.0;
    resinfo.b = 0.0;
    resinfo.c = 0.48;
    resinfo.tol = 5.0e-3;
    resinfo.direction = Vector3(0.0,1.0,0.0);
    
    unique_ptr<SpinWavePlot> res(new TwoDimensionResolutionFunction(resinfo, SW, energies));
    
    HKLDirections direction;
    direction.addDirection(0.0,0.0,1.0,0.2);
    direction.addDirection(1.0,0.0,0.0,0.2);
    //direction.addDirection(0.0,1.0,0.0,0.05);
    
    unique_ptr<SpinWavePlot> asdf(new IntegrateAxes(std::move(res),direction,1.0e-2));

    TwoDimensionalCut twodimcut;
    twodimcut.setFilename("YFeO3_2_K_m3");
    twodimcut.setPlotObject(move(asdf));
    twodimcut.setPoints(kPoints);
    twodimcut.save();

    Line.setFirstPoint(2.0,-1.5,-2.0);
    Line.setFinalPoint(2.0, 1.5,-2.0);
    Line.setNumberPoints(201);
    kPoints = Line.getPoints();

    twodimcut.setFilename("YFeO3_2_K_m2");
    twodimcut.setPoints(kPoints);
    twodimcut.save();
}

void alongLAxis()
{
    SpinWave SW = createModel();

    PointsAlongLine Line;
    Line.setFirstPoint(3.0,0.0,-4.0);
    Line.setFinalPoint(3.0,0.0, 4.0);
    Line.setNumberPoints(101);
    ThreeVectors<double> kPoints = Line.getPoints();

    Energies energies(0.0, 80.0, 201);

    TwoDimGaussian resinfo;
    resinfo.a = 579.7;
    resinfo.b = -20.0;
    resinfo.c = 1.3;
    resinfo.tol = 7.0e-3;
    resinfo.direction = Vector3(0.0,0.0,1.0);

    unique_ptr<SpinWavePlot> res(new TwoDimensionResolutionFunction(resinfo, SW, energies));

    HKLDirections direction;
    direction.addDirection(1.0,0.0,0.0,0.2);
    direction.addDirection(0.0,1.0,0.0,0.2);
    //direction.addDirection(0.0,0.0,1.0,0.05);
    
    unique_ptr<SpinWavePlot> asdf(new IntegrateAxes(std::move(res),direction,7.0e-3));
   
    TwoDimensionalCut twodimcut;
    twodimcut.setFilename("YFeO3_3_0_L");
    twodimcut.setPlotObject(move(asdf));
    twodimcut.setPoints(kPoints);
    twodimcut.save();
   
    Line.setFirstPoint(3.0,1.0,-4.0);
    Line.setFinalPoint(3.0,1.0, 4.0);
    Line.setNumberPoints(201);
    kPoints = Line.getPoints();
    twodimcut.setFilename("YFeO3_3_1_L");
    twodimcut.setPoints(kPoints);
    //twodimcut.save();
}

int main()
{
    //alongKAxis();
    alongLAxis();
    return 0;
}
