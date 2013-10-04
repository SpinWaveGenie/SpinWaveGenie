#include "Initializer.h"
#include "Exch_Interaction.h"
#include "Anis_Z_Interaction.h"
#include "Anis_Y_Interaction.h"
#include "Anis_X_Interaction.h"
#include "DM_Y_Interaction.h"

using namespace std;
using namespace boost;

void Init::read_input(string filename)
{    //Sublattice Fe1,Fe2;
    
    m_parser_map["crystal"] = bind(&Init::parseCrystalNode, this, _1);
    m_parser_map["interaction"] = bind(&Init::parseInteractionNode, this, _1);
    m_parser_map["kpoints"] = bind(&Init::parseKpointsNode, this, _1);
    
    pugi::xml_document doc;
    
    pugi::xml_parse_result result = doc.load_file(filename.c_str());
    std::cout << "Load result: " << result.description() << std::endl;
    assert(result);
    pugi::xml_node tools = doc.child("spin_wave_genie");
    
    for (pugi::xml_node_iterator it = tools.begin(); it != tools.end(); ++it)
    {    
        string name = it->name();
        pugi::xml_node node = tools.child(name.c_str());
        std::map< std::string, boost::function< void (const pugi::xml_node&) > >::iterator parser;
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
    unit_cell = boost::shared_ptr<Cell>(new Cell);
    istringstream parser;
    double scale = node.child("basevect").attribute("scale").as_double();
    cout << scale << endl;
    
    /*parser.str(node.child_value("basevect"));
    Eigen::Matrix3d basis;
    int i = 0;
    while(parser.good())
    {
        double x,y,z;
        parser >> x >> y >> z;
        if (parser.good())
        {
            //cout << x << '\t' << y << '\t' << z << endl;
            basis(i,0) = x;
            basis(i,1) = y;
            basis(i,2) = z;
            i++;
        }
    }
    cout << basis << endl;
            
    unit_cell->setBasisVectors(scale,basis);*/
    pugi::xml_node lattice_parameters = node.child("basevect");
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

    unit_cell->setBasisVectors(a,b,c,alpha,beta,gamma);

    for (pugi::xml_node tool = node.child("sublattice"); tool; tool = tool.next_sibling("sublattice"))
    {
        Sublattice* new_sl(new Sublattice());
            
        string name(tool.child_value("name"));
        name.erase(std::remove_if(name.begin(),name.end(), ::isspace), name.end());
        cout << name << endl;
        new_sl->setName(name);
                
        parser.clear();
        parser.str(tool.child_value("moment"));
        double S,theta,phi;
        parser >> S >> theta >> phi;
        cout << S << '\t' << theta << '\t' << phi << endl;
        new_sl->setMoment(S,theta*M_PI/180.0,phi*M_PI/180.0);

        unit_cell->addSublattice(name,std::move(new_sl));
        parser.clear();
        parser.str(tool.child_value("position"));
        while(parser.good())
        {
            double x,y,z;
            parser >> x >> y >> z;
            if (parser.good())
            {
                unit_cell->addAtom(name,x,y,z);
                cout << x << '\t' << y << '\t' << z << endl;
            }
        }
    }
}

void Init::parseInteractionNode(const pugi::xml_node &node)
{
    istringstream parser;
    double min,max;
    string atom1,atom2;
    builder = boost::shared_ptr<SW_Builder>(new SW_Builder(unit_cell));
    for (pugi::xml_node tool = node.child("exchange"); tool; tool = tool.next_sibling("exchange"))
    {
        double value = tool.attribute("value").as_double();
        parser.clear();
        parser.str(tool.child_value());
        while(parser.good())
        {
            parser >> atom1 >> atom2 >> min >> max;
            if(parser.good())
            {
                cout << atom1 << atom2 << endl;
                builder->Add_Interaction(new Exch_Interaction(value,atom1,atom2,min,max));
                if (atom1 != atom2)
                {
                    cout << atom2 << atom1 << endl;
                    builder->Add_Interaction(new Exch_Interaction(value,atom2,atom1,min,max));
                }
            }
        }
    }
    
    for (pugi::xml_node tool = node.child("DM_y"); tool; tool = tool.next_sibling("DM_y"))
    {
        double value = tool.attribute("value").as_double();
        parser.clear();
        parser.str(tool.child_value());
        while(parser.good())
        {
            parser >> atom1 >> atom2 >> min >> max;
            if(parser.good())
            {
                cout << atom1 << atom2 << endl;
                builder->Add_Interaction(new DM_Y_Interaction(value,atom2,atom1,min,max));
            }
        }
    }
    
    for (pugi::xml_node tool = node.child("anis_z"); tool; tool = tool.next_sibling("anis_z"))
    {
        double value = tool.attribute("value").as_double();
        parser.clear();
        parser.str(tool.child_value());
        while(parser.good())
        {
            parser >> atom1;
            if(parser.good())
            {
                cout << atom1 << endl;
                builder->Add_Interaction(new Anis_Z_Interaction(value,atom1));
            }
        }
    }
    
    for (pugi::xml_node tool = node.child("anis_y"); tool; tool = tool.next_sibling("anis_y"))
    {
        double value = tool.attribute("value").as_double();
        parser.clear();
        parser.str(tool.child_value());
        while(parser.good())
        {
            parser >> atom1;
            if(parser.good())
            {
                cout << atom1 << endl;
                builder->Add_Interaction(new Anis_Y_Interaction(value,atom1));
            }
        }
        
    }
    
    for (pugi::xml_node tool = node.child("anis_x"); tool; tool = tool.next_sibling("anis_x"))
    {
        double value = tool.attribute("value").as_double();
        parser.clear();
        parser.str(tool.child_value());
        while(parser.good())
        {
            parser >> atom1;
            if(parser.good())
            {
                cout << atom1 << endl;
                builder->Add_Interaction(new Anis_X_Interaction(value,atom1));
            }
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
        //ProgressBar pbar(npoints);
        //pbar.start();
        while(parser.good())
        {
            parser >> x0 >> y0 >> z0;
            parser >> x1 >> y1 >> z1;
            if(parser.good())
            {
                for(int n=0;n!=npoints;n++)
                {
                    //pbar.update(n);
                    x = x0 + (x1-x0)*n/(npoints-1);
                    y = y0 + (y1-y0)*n/(npoints-1);
                    z = z0 + (z1-z0)*n/(npoints-1);
                    SpinWave test = builder->Create_Element(x,y,z);
                    //cout << x << " " << y << " " << z << endl;
                    test.Calc();
                    vector<double> frequencies = test.Get_Frequencies();
                    /*cout << endl;
                    for(int i=0;i!=frequencies.size();i++)
                    {
                        cout << frequencies[i] << endl;
                    }*/
                }
            }
        }
        //pbar.finish();
    }
}

boost::shared_ptr<Cell> Init::get_cell()
{
    return unit_cell;
}

boost::shared_ptr<SW_Builder> Init::get_builder()
{
    return builder;
}


