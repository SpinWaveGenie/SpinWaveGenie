#ifndef __Initializer__
#define __Initializer__
#include <map>

#include "Cell/Cell.h"
#include "Genie/SpinWaveBuilder.h"
#include "External/pugixml.hpp"

class Init
{
public:
    Init(std::string filename);
    //! reads input from xml file "filename"
    void read_input(std::string filename);
    void save_input(std::string filename);
    // saves Cell into xml file "filename"
    // not yet implemented
    SpinWaveGenie::Cell get_cell();
    SpinWaveGenie::SpinWaveBuilder get_builder();
    
private:
    std::map< std::string, std::function< void (const pugi::xml_node&) > > m_parser_map; //!< Map for dispatching parsers based on node type
    SpinWaveGenie::Cell unit_cell;
    SpinWaveGenie::SpinWaveBuilder builder;
    void parseCrystalNode(const pugi::xml_node &node);
    void parseInteractionNode(const pugi::xml_node &node);
    void parseDispersion(const pugi::xml_node &node);
    void parseTwoDimensionCut( const pugi::xml_node &node);
};

#endif /* defined(__Initializer__) */
