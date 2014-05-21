#ifndef __Initializer__
#define __Initializer__

#include <iostream>
#include <sstream>
#include <new>
#include <functional>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/trim.hpp>
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
    Cell get_cell();
    SpinWaveBuilder get_builder();
    
private:
    std::map< std::string, std::function< void (const pugi::xml_node&) > > m_parser_map; //!< Map for dispatching parsers based on node type
    Cell unit_cell;
    SpinWaveBuilder builder;
    void parseCrystalNode(const pugi::xml_node &node);
    void parseInteractionNode(const pugi::xml_node &node);
    void parseDispersion(const pugi::xml_node &node);
    void parseTwoDimensionCut( const pugi::xml_node &node);
};

#endif /* defined(__Initializer__) */
