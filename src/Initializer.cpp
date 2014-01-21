#include <functional>
#include "Initializer.h"
#include "Interactions/Exch_Interaction.h"
#include "Interactions/AnisotropyInteraction.h"
#include "Interactions/DM_Y_Interaction.h"
#include "Interactions/DM_Z_Interaction.h"
#include "SpinWaveDispersion.h"
#include "TwoDimensionCut.h"
#include "PointsAlongLine.h"
#include "Containers/Positions.h"


using namespace std;
using namespace boost;

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
        cout << name << endl;
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
    
    string temp = lattice_parameters.child_value("a");
    algorithm::trim(temp); // get rid of surrounding whitespace
    a = lexical_cast<double>(temp);
    
    temp = lattice_parameters.child_value("b");
    algorithm::trim(temp); // get rid of surrounding whitespace
    b = lexical_cast<double>(temp);
    
    temp = lattice_parameters.child_value("c");
    algorithm::trim(temp); // get rid of surrounding whitespace
    c = lexical_cast<double>(temp);
    
    temp = lattice_parameters.child_value("alpha");
    algorithm::trim(temp); // get rid of surrounding whitespace
    alpha = lexical_cast<double>(temp);
    
    temp = lattice_parameters.child_value("beta");
    algorithm::trim(temp); // get rid of surrounding whitespace
    beta = lexical_cast<double>(temp);
    
    temp = lattice_parameters.child_value("gamma");
    algorithm::trim(temp); // get rid of surrounding whitespace
    gamma = lexical_cast<double>(temp);
    
    cout << a << " " << b << " " << c << " " << alpha << " " << beta << " " << gamma << endl;

    unit_cell.setBasisVectors(a,b,c,alpha,beta,gamma);

    pugi::xml_node moments = node.child("moments");

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
        
        temp = tool.child("moment").child_value("magnitude");
        algorithm::trim(temp); // get rid of surrounding whitespace
        S = lexical_cast<double>(temp);
        
        temp = tool.child("moment").child_value("theta");
        algorithm::trim(temp); // get rid of surrounding whitespace
        theta = lexical_cast<double>(temp)*M_PI/180.0;
        
        temp = tool.child("moment").child_value("phi");
        algorithm::trim(temp); // get rid of surrounding whitespace
        phi = lexical_cast<double>(temp)*M_PI/180.0;
        
        cout << S << '\t' << theta << '\t' << phi << endl;

        new_sl.setMoment(S,theta,phi);

        unit_cell.addSublattice(new_sl);
        
        pugi::xml_node atomicpositions = tool.child("atomicpositions");
        for (pugi::xml_node position = atomicpositions.child("position"); position; position = position.next_sibling("position"))
        {
            double x,y,z;
            temp = position.child_value("x");
            algorithm::trim(temp); // get rid of surrounding whitespace
            x = lexical_cast<double>(temp);
            
            temp = position.child_value("y");
            algorithm::trim(temp); // get rid of surrounding whitespace
            y = lexical_cast<double>(temp);
            
            temp = position.child_value("z");
            algorithm::trim(temp); // get rid of surrounding whitespace
            z = lexical_cast<double>(temp);
            
            cout << x << '\t' << y << '\t' << z << endl;
           
            unit_cell.addAtom(name,x,y,z);
            
        }
    }
}

void Init::parseInteractionNode(const pugi::xml_node &node)
{
     SW_Builder buildertemp(unit_cell);
     builder = buildertemp;
    
     pugi::xml_node exchange = node.child("Exchange");
     for (pugi::xml_node tool = exchange.child("group"); tool; tool = tool.next_sibling("group"))
     {
         cout << "Exchange: " << endl;
         
         string name = tool.child_value("name");
         algorithm::trim(name); // get rid of surrounding whitespace
         
         string temp = tool.child_value("value");
         algorithm::trim(temp); // get rid of surrounding whitespace
         double value = lexical_cast<double>(temp);
         
         temp = tool.child_value("mindist");
         algorithm::trim(temp); // get rid of surrounding whitespace
         double min = lexical_cast<double>(temp);
         
         temp = tool.child_value("maxdist");
         algorithm::trim(temp); // get rid of surrounding whitespace
         double max = lexical_cast<double>(temp);

         pugi::xml_node pairs = tool.child("pairs");
         for (pugi::xml_node pair = pairs.child("pair"); pair; pair = pair.next_sibling("pair"))
         {
             string atom1 = pair.child_value("name1");
             string atom2 = pair.child_value("name2");
             
             cout << name << " " << value <<" " << atom1 << " " <<  atom2 << " " << min << " " << " " << max << endl;
             builder.Add_Interaction(new Exch_Interaction(name,value,atom1,atom2,min,max));
         }
     }
    
    pugi::xml_node DzyaloshinskyMoriya = node.child("DzyaloshinskyMoriya");
    for (pugi::xml_node tool = DzyaloshinskyMoriya.child("group"); tool; tool = tool.next_sibling("group"))
    {
        cout << "Dzyaloshinsky-Moriya: " << endl;
        
        string name = tool.child_value("name");
        algorithm::trim(name); // get rid of surrounding whitespace
        
        string temp = tool.child_value("value");
        algorithm::trim(temp); // get rid of surrounding whitespace
        double value = lexical_cast<double>(temp);
        
        temp = tool.child_value("min");
        algorithm::trim(temp); // get rid of surrounding whitespace
        double min = lexical_cast<double>(temp);
        
        temp = tool.child_value("max");
        algorithm::trim(temp); // get rid of surrounding whitespace
        double max = lexical_cast<double>(temp);
        
        double x,y,z;
        
        temp = tool.child("direction").child_value("x");
        algorithm::trim(temp); // get rid of surrounding whitespace
        x = lexical_cast<double>(temp);
        
        temp = tool.child("direction").child_value("y");
        algorithm::trim(temp); // get rid of surrounding whitespace
        y = lexical_cast<double>(temp);
        
        temp = tool.child("direction").child_value("z");
        algorithm::trim(temp); // get rid of surrounding whitespace
        z = lexical_cast<double>(temp);
        
        cout << x << '\t' << y << '\t' << z << endl;
        
        pugi::xml_node pairs = tool.child("pairs");
        for (pugi::xml_node pair = pairs.child("pair"); pair; pair = pair.next_sibling("pair"))
        {
            string atom1 = pair.child_value("name1");
            string atom2 = pair.child_value("name2");
            
            cout << name << " " << value <<" " << atom1 << " " <<  atom2 << " " << min << " " << " " << max << endl;
            
            if ( x > 0.99 && y < 0.01 && z < 0.01)
            {
                cout << "Error: Dzyaloshinsky-Moriya along X not yet implemented!" << endl;
            }
            else if (x < 0.01 && y > 0.99 && z < 0.01)
            {
                builder.Add_Interaction(new DM_Y_Interaction(name,value,atom1,atom2,min,max));
            }
            else if (x < 0.01 && y < 0.01 && z > 0.99)
            {
                builder.Add_Interaction(new DM_Z_Interaction(name,value,atom1,atom2,min,max));
            }
        }
    }
    
    pugi::xml_node anisotropy = node.child("Anisotropy");
    
    for (pugi::xml_node tool = anisotropy.child("group"); tool; tool = tool.next_sibling("group"))
    {
        cout << "Anisotropy: " << endl;
        
        string identifier = tool.child_value("name");
        algorithm::trim(identifier); // get rid of surrounding whitespace
        
        string temp = tool.child_value("value");
        algorithm::trim(temp); // get rid of surrounding whitespace
        double value = lexical_cast<double>(temp);
        
        double x,y,z;
        
        temp = tool.child("direction").child_value("x");
        algorithm::trim(temp); // get rid of surrounding whitespace
        x = lexical_cast<double>(temp);
        
        temp = tool.child("direction").child_value("y");
        algorithm::trim(temp); // get rid of surrounding whitespace
        y = lexical_cast<double>(temp);
        
        temp = tool.child("direction").child_value("z");
        algorithm::trim(temp); // get rid of surrounding whitespace
        z = lexical_cast<double>(temp);
    
        cout << x << '\t' << y << '\t' << z << endl;
            
        pugi::xml_node sublattices = tool.child("sublattices");
        for (pugi::xml_node name = sublattices.child("name"); name; name = name.next_sibling("name"))
        {
            string atom1 = name.child_value();
            Vector3 direction(x,y,z);
            //cout "direction= " << direction.transpose() << endl;
            builder.Add_Interaction(new AnisotropyInteraction(identifier,value,direction,atom1));
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
        string temp;
        
        temp = group.child_value("numberpoints");
        algorithm::trim(temp); // get rid of surrounding whitespace
        NumberPoints = lexical_cast<double>(temp);
        
        temp = group.child("firstpoint").child_value("x");
        algorithm::trim(temp); // get rid of surrounding whitespace
        x0 = lexical_cast<double>(temp);
        
        temp = group.child("firstpoint").child_value("y");
        algorithm::trim(temp); // get rid of surrounding whitespace
        y0 = lexical_cast<double>(temp);
        
        temp = group.child("firstpoint").child_value("z");
        algorithm::trim(temp); // get rid of surrounding whitespace
        z0 = lexical_cast<double>(temp);
        
        temp = group.child("lastpoint").child_value("x");
        algorithm::trim(temp); // get rid of surrounding whitespace
        x1 = lexical_cast<double>(temp);
        
        temp = group.child("lastpoint").child_value("y");
        algorithm::trim(temp); // get rid of surrounding whitespace
        y1 = lexical_cast<double>(temp);
        
        temp = group.child("lastpoint").child_value("z");
        algorithm::trim(temp); // get rid of surrounding whitespace
        z1 = lexical_cast<double>(temp);
        
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
        
        temp = group.child_value("numberpoints");
        algorithm::trim(temp); // get rid of surrounding whitespace
        NumberPoints = lexical_cast<double>(temp);
        
        temp = group.child("firstpoint").child_value("x");
        algorithm::trim(temp); // get rid of surrounding whitespace
        x0 = lexical_cast<double>(temp);
        
        temp = group.child("firstpoint").child_value("y");
        algorithm::trim(temp); // get rid of surrounding whitespace
        y0 = lexical_cast<double>(temp);
        
        temp = group.child("firstpoint").child_value("z");
        algorithm::trim(temp); // get rid of surrounding whitespace
        z0 = lexical_cast<double>(temp);
        
        temp = group.child("lastpoint").child_value("x");
        algorithm::trim(temp); // get rid of surrounding whitespace
        x1 = lexical_cast<double>(temp);
        
        temp = group.child("lastpoint").child_value("y");
        algorithm::trim(temp); // get rid of surrounding whitespace
        y1 = lexical_cast<double>(temp);
        
        temp = group.child("lastpoint").child_value("z");
        algorithm::trim(temp); // get rid of surrounding whitespace
        z1 = lexical_cast<double>(temp);
        
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
        
        temp = group.child_value("numberpoints");
        algorithm::trim(temp); // get rid of surrounding whitespace
        NumberPoints = lexical_cast<double>(temp);
        
        temp = group.child_value("firstpoint");
        algorithm::trim(temp); // get rid of surrounding whitespace
        MinEnergy = lexical_cast<double>(temp);
        
        temp = group.child_value("lastpoint");
        algorithm::trim(temp); // get rid of surrounding whitespace
        MaxEnergy = lexical_cast<double>(temp);
        
        cout << MinEnergy << " " << MaxEnergy << " " << NumberPoints << endl;
        
        PointsAlongLine Line;
        Cut.setEnergyPoints(MinEnergy,MaxEnergy,NumberPoints);
    }

    OneDimGaussian resinfo;
    resinfo.fwhm = 1.0;
    resinfo.tol = 1.0e-5;
    resinfo.SW = builder.Create_Element();
    Cut.setConvolutionObject(resinfo);
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


