#include <functional>
#include "Initializer.h"
#include "Interactions/InteractionFactory.h"
#include "Interactions/ExchangeInteraction.h"
#include "Interactions/AnisotropyInteraction.h"
#include "Interactions/DM_Y_Interaction.h"
#include "Interactions/DM_Z_Interaction.h"
#include "SpinWaveDispersion.h"
#include "TwoDimensionCut.h"
#include "PointsAlongLine.h"
#include "Containers/ThreeVectors.h"
#include "OneDimensionalFactory.h"
#include "OneDimensionalShapes.h"
#include "EnergyResolutionFunction.h"

using namespace std;
using namespace boost;

double stringToDouble(string value)
{
    algorithm::trim(value); // get rid of surrounding whitespace
    return lexical_cast<double>(value);
}

Init::Init(string filename)
{
    read_input(filename);
}

void Init::read_input(string filename)
{    //Sublattice Fe1,Fe2;
    
    m_parser_map["cell"] = bind(&Init::parseCrystalNode, this, _1);
    m_parser_map["interactions"] = bind(&Init::parseInteractionNode, this, _1);
    m_parser_map["dispersion"] = bind(&Init::parseDispersion, this, _1);
    m_parser_map["TwoDimensionCut"] = bind(&Init::parseTwoDimensionCut, this, _1);
    
    pugi::xml_document doc;
    
    pugi::xml_parse_result result = doc.load_file(filename.c_str());
    std::cout << "Load result: " << result.description() << std::endl;
    assert(result);
    
    
    pugi::xml_node tools = doc.child("spinwavegenie");
    for (pugi::xml_node_iterator it = tools.begin(); it != tools.end(); ++it)
    {    
        string name = it->name();
        cout << "name" << name << endl;
        pugi::xml_node node = (*it);
        std::map< std::string, std::function< void (const pugi::xml_node&) > >::iterator parser;
        parser = m_parser_map.find(name);
        if (parser != m_parser_map.end())
            parser->second(node);
        else
            cout << "Notice: Parser for node <" << name << "> not defined, ignoring" << endl;
    }
}

void Init::save_input(string filename)
{
    cout << "Not yet implemented" << endl;
}

void Init::parseCrystalNode(const pugi::xml_node &node)
{
    pugi::xml_node lattice_parameters = node.child("latticevectors");
    double a,b,c,alpha,beta,gamma;
    
    a = stringToDouble(lattice_parameters.child_value("a"));
    b = stringToDouble(lattice_parameters.child_value("b"));
    c = stringToDouble(lattice_parameters.child_value("c"));
    alpha = stringToDouble(lattice_parameters.child_value("alpha"));
    beta = stringToDouble(lattice_parameters.child_value("beta"));
    gamma = stringToDouble(lattice_parameters.child_value("gamma"));
    
    cout << a << " " << b << " " << c << " " << alpha << " " << beta << " " << gamma << endl;

    unit_cell.setBasisVectors(a,b,c,alpha,beta,gamma);

    pugi::xml_node moments = node.child("moments");

    string temp;
    
    for (pugi::xml_node tool = moments.child("sublattice"); tool; tool = tool.next_sibling("sublattice"))
    {
        Sublattice new_sl;
            
        string name(tool.child_value("name"));
        algorithm::trim(name);
        cout << name << endl;
        new_sl.setName(name);
        
        string type(tool.child_value("type"));
        algorithm::trim(type);
        cout << type << endl;
        new_sl.setType(type);
        
        double S,theta,phi;
        
        S = stringToDouble(tool.child("moment").child_value("magnitude"));
        theta = stringToDouble(tool.child("moment").child_value("theta"))*M_PI/180.0;
        phi = stringToDouble(tool.child("moment").child_value("phi"))*M_PI/180.0;
        
        cout << S << '\t' << theta << '\t' << phi << endl;

        new_sl.setMoment(S,theta,phi);

        unit_cell.addSublattice(new_sl);
        
        pugi::xml_node atomicpositions = tool.child("atomicpositions");
        for (pugi::xml_node position = atomicpositions.child("position"); position; position = position.next_sibling("position"))
        {
            double x,y,z;
            
            x = stringToDouble(position.child_value("x"));
            y = stringToDouble(position.child_value("y"));
            z = stringToDouble(position.child_value("z"));
            
            cout << x << '\t' << y << '\t' << z << endl;
           
            unit_cell.addAtom(name,x,y,z);
        }
    }
}

void Init::parseInteractionNode(const pugi::xml_node &node)
{
     SW_Builder buildertemp(unit_cell);
     builder = buildertemp;
     InteractionFactory factory;
    
     pugi::xml_node exchange = node.child("Exchange");
     for (pugi::xml_node tool = exchange.child("group"); tool; tool = tool.next_sibling("group"))
     {
         cout << "Exchange: " << endl;
         
         string name = tool.child_value("name");
         algorithm::trim(name); // get rid of surrounding whitespace

         double value,min,max;
         
         value = stringToDouble(tool.child_value("value"));
         min = stringToDouble(tool.child_value("mindist"));
         max = stringToDouble(tool.child_value("maxdist"));
         
         pugi::xml_node pairs = tool.child("pairs");
         for (pugi::xml_node pair = pairs.child("pair"); pair; pair = pair.next_sibling("pair"))
         {
             string atom1 = pair.child_value("name1");
             string atom2 = pair.child_value("name2");
             
             cout << name << " " << value <<" " << atom1 << " " <<  atom2 << " " << min << " " << " " << max << endl;
             builder.addInteraction(factory.getExchange(name,value,atom1,atom2,min,max));

         }
     }
    
    pugi::xml_node DzyaloshinskyMoriya = node.child("DzyaloshinskyMoriya");
    for (pugi::xml_node tool = DzyaloshinskyMoriya.child("group"); tool; tool = tool.next_sibling("group"))
    {
        cout << "Dzyaloshinsky-Moriya: " << endl;
        
        string name = tool.child_value("name");
        algorithm::trim(name); // get rid of surrounding whitespace
        
        double value,min,max;
        value = stringToDouble(tool.child_value("value"));
        min = stringToDouble(tool.child_value("min"));
        max = stringToDouble(tool.child_value("max"));
        
        double x,y,z;
        
        x = stringToDouble(tool.child("direction").child_value("x"));
        y = stringToDouble(tool.child("direction").child_value("y"));
        z = stringToDouble(tool.child("direction").child_value("z"));
        
        cout << x << '\t' << y << '\t' << z << endl;
        
        Vector3 direction(x,y,z);
        
        pugi::xml_node pairs = tool.child("pairs");
        for (pugi::xml_node pair = pairs.child("pair"); pair; pair = pair.next_sibling("pair"))
        {
            string atom1 = pair.child_value("name1");
            string atom2 = pair.child_value("name2");
            
            cout << name << " " << value <<" " << atom1 << " " <<  atom2 << " " << min << " " << " " << max << endl;
            builder.addInteraction( factory.getDzyaloshinskiiMoriya(name,value,direction,atom1,atom2,min,max) );
        }
    }
    
    pugi::xml_node anisotropy = node.child("Anisotropy");
    
    for (pugi::xml_node tool = anisotropy.child("group"); tool; tool = tool.next_sibling("group"))
    {
        cout << "Anisotropy: " << endl;
        
        string identifier = tool.child_value("name");
        algorithm::trim(identifier); // get rid of surrounding whitespace
        
        double value,x,y,z;
        
        value = stringToDouble(tool.child_value("value"));
        x = stringToDouble(tool.child("direction").child_value("x"));
        y = stringToDouble(tool.child("direction").child_value("y"));
        z = stringToDouble(tool.child("direction").child_value("z"));
    
        cout << x << '\t' << y << '\t' << z << endl;
            
        pugi::xml_node sublattices = tool.child("sublattices");
        for (pugi::xml_node name = sublattices.child("name"); name; name = name.next_sibling("name"))
        {
            string atom1 = name.child_value();
            Vector3 direction(x,y,z);
            cout << "direction= " << direction.transpose() << endl;
            builder.addInteraction(factory.getAnisotropy(identifier,value,direction,atom1));
        }
    }
}

void Init::parseDispersion(const pugi::xml_node &node)
{
    
    cout << "Dispersion:" << endl;
    string filename = node.child_value("filename");
    string filetype = node.child_value("filetype");
    
    SpinWaveDispersion Dispersion;
    Dispersion.setFilename(filename);
    
    SpinWaveDispersion::Options PrintOptions;
    string temp;
    bool value;
    
    temp = node.child("filetype").child_value("printposition");
    value = boost::lexical_cast<bool, std::string>(temp);
    PrintOptions = SpinWaveDispersion::Options::PrintPosition;
    cout << "print position: "<< value << endl;
    Dispersion.setOptions(PrintOptions,value);
    
    temp = node.child("filetype").child_value("printfrequency");
    value = boost::lexical_cast<bool, std::string>(temp);
    PrintOptions = SpinWaveDispersion::Options::PrintFrequency;
    cout << "print frequency: "<< value << endl;
    Dispersion.setOptions(PrintOptions,value);
    
    temp = node.child("filetype").child_value("printintensity");
    value = boost::lexical_cast<bool, std::string>(temp);
    PrintOptions = SpinWaveDispersion::Options::PrintIntensity;
    cout << "print intensity: "<< value << endl;
    Dispersion.setOptions(PrintOptions,value);
    
    pugi::xml_node lines = node.child("lines");
    for (pugi::xml_node group = lines.child("group");group; group = group.next_sibling("group"))
    {
        double x0,y0,z0,x1,y1,z1,NumberPoints;
        
        NumberPoints = stringToDouble(group.child_value("numberpoints"));
        x0 = stringToDouble(group.child("firstpoint").child_value("x"));
        y0 = stringToDouble(group.child("firstpoint").child_value("y"));
        z0 = stringToDouble(group.child("firstpoint").child_value("z"));
        x1 = stringToDouble(group.child("lastpoint").child_value("x"));
        y1 = stringToDouble(group.child("lastpoint").child_value("y"));
        z1 = stringToDouble(group.child("lastpoint").child_value("z"));
        
        cout << x0 << " " << y0 << " " << z0 << " " << x1 << " " << y1 << "  " << z1 << endl;
        
        PointsAlongLine Line;
        Line.setFirstPoint(x0,y0,z0);
        Line.setFinalPoint(x1,y1,z1);
        Line.setNumberPoints(NumberPoints);
        Dispersion.setPoints(Line.getPoints());
    }
    
    Dispersion.setGenie(builder.Create_Element());
    Dispersion.save();
}

void Init::parseTwoDimensionCut(const pugi::xml_node &node)
{
    
    cout << "Two Dimension Cut:" << endl;
    string filename = node.child_value("filename");
    string filetype = node.child_value("filetype");
    
    TwoDimensionCut Cut;
    Cut.setFilename(filename);
    
    cout << filename << endl;
    
    {
        pugi::xml_node group = node.child("setkpoints");
        double x0,y0,z0,x1,y1,z1,NumberPoints;
        string temp;
        
        NumberPoints = stringToDouble(group.child_value("numberpoints"));
        x0 = stringToDouble(group.child("firstpoint").child_value("x"));
        y0 = stringToDouble(group.child("firstpoint").child_value("y"));
        z0 = stringToDouble(group.child("firstpoint").child_value("z"));
        x1 = stringToDouble(group.child("lastpoint").child_value("x"));
        y1 = stringToDouble(group.child("lastpoint").child_value("y"));
        z1 = stringToDouble(group.child("lastpoint").child_value("z"));
        
        cout << x0 << " " << y0 << " " << z0 << " " << x1 << " " << y1 << "  " << z1 << endl;
        
        PointsAlongLine Line;
        Line.setFirstPoint(x0,y0,z0);
        Line.setFinalPoint(x1,y1,z1);
        Line.setNumberPoints(NumberPoints);
        Cut.setPoints(Line.getPoints());
        

    }
        
    {
        pugi::xml_node group = node.child("setenergypoints");
        double MinEnergy,MaxEnergy,NumberPoints;
        string temp;
        
        NumberPoints = stringToDouble(group.child_value("numberpoints"));
        MinEnergy = stringToDouble(group.child_value("firstpoint"));
        MaxEnergy = stringToDouble(group.child_value("lastpoint"));
        
        cout << MinEnergy << " " << MaxEnergy << " " << NumberPoints << endl;
        
        PointsAlongLine Line;

        pugi::xml_node type = node.child("type");
        
        pugi::xml_node Gaussian = type.child("OneDimensionGaussian");
        pugi::xml_node Lorentzian = type.child("OneDimensionLorentzian");
        pugi::xml_node PseudoVoigt = type.child("OneDimensionPseudoVoigt");
        
        OneDimensionalFactory factory;
        
        unique_ptr<OneDimensionalShapes> resinfo;
        
        if (Gaussian)
        {
            double fwhm,tolerance;
            fwhm = stringToDouble(Gaussian.child_value("fwhm"));
            tolerance = stringToDouble(Gaussian.child_value("tol"));
            
            resinfo = factory.getGaussian(fwhm,tolerance);
            cout << "Gaussian resolution function set" << endl;
        }
        else if(Lorentzian)
        {
            double fwhm,tolerance;
            fwhm = stringToDouble(Lorentzian.child_value("fwhm"));
            tolerance = stringToDouble(Lorentzian.child_value("tol"));

            resinfo = factory.getLorentzian(fwhm,tolerance);
            
            cout << "Lorentzian resolution function set" << endl;
        }
        else if(PseudoVoigt)
        {
            double eta,fwhm,tolerance;
            eta = stringToDouble(PseudoVoigt.child_value("eta"));
            fwhm = stringToDouble(PseudoVoigt.child_value("fwhm"));
            tolerance = stringToDouble(PseudoVoigt.child_value("tol"));
            
            resinfo = factory.getPseudoVoigt(eta,fwhm,tolerance);
            
            cout << "Pseudo-Voigt resolution function set" << endl;
        }
        else
        {
            cout << "RESOLUTION FUNCTION NOT SET!!!" << endl;
        }
        
        SpinWave SW = builder.Create_Element();
        unique_ptr<SpinWavePlot> res(new EnergyResolutionFunction(move(resinfo), SW, MinEnergy, MaxEnergy, NumberPoints));
        Cut.setEnergyPoints(MinEnergy,MaxEnergy,NumberPoints);
        Cut.setPlotObject(move(res));
    }
    Cut.save();


}

Cell Init::get_cell()
{
    return unit_cell;
}

SW_Builder Init::get_builder()
{
    return builder;
}


