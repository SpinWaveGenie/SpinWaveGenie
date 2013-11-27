#include <functional>
#include "Initializer.h"
#include "Exch_Interaction.h"
#include "Anis_Z_Interaction.h"
#include "Anis_Y_Interaction.h"
#include "Anis_X_Interaction.h"
#include "AnisotropyInteraction.h"
#include "DM_Y_Interaction.h"
#include "DM_Z_Interaction.h"

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
    m_parser_map["kpoints"] = bind(&Init::parseKpointsNode, this, _1);
    
    pugi::xml_document doc;
    
    pugi::xml_parse_result result = doc.load_file(filename.c_str());
    std::cout << "Load result: " << result.description() << std::endl;
    assert(result);
    pugi::xml_node tools = doc.child("spinwavegenie");
    
    for (pugi::xml_node_iterator it = tools.begin(); it != tools.end(); ++it)
    {    
        string name = it->name();
        cout << name << endl;

        pugi::xml_node node = tools.child(name.c_str());
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
             if (atom1 != atom2)
             {
                 cout << atom2 << " " << atom1 << endl;
                 builder.Add_Interaction(new Exch_Interaction(name,value,atom2,atom1,min,max));
             }
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

void Init::parseKpointsNode(const pugi::xml_node &node)
{
    double x,y,z,x0,y0,z0,x1,y1,z1;
    if(strncmp(node.attribute("type").value(),"line",4) == 0)
    {
        double npoints = node.attribute("points").as_double();
        istringstream parser;
        parser.str(node.child_value());
        while(parser.good())
        {
            parser >> x0 >> y0 >> z0;
            parser >> x1 >> y1 >> z1;
            if(parser.good())
            {
                for(int n=0;n!=npoints;n++)
                {
                    x = x0 + (x1-x0)*n/(npoints-1);
                    y = y0 + (y1-y0)*n/(npoints-1);
                    z = z0 + (z1-z0)*n/(npoints-1);
                    //SpinWave* test = builder.Create_Element(x,y,z);
                    //cout << x << " " << y << " " << z << endl;
                    //test->Calc();
                    //vector<double> frequencies = test->Get_Frequencies();
                    /*cout << endl;
                    for(int i=0;i!=frequencies.size();i++)
                    {
                        cout << frequencies[i] << endl;
                    }*/
                }
            }
        }
    }
}

Cell Init::get_cell()
{
    return unit_cell;
}

SW_Builder Init::get_builder()
{
    return builder;
}


