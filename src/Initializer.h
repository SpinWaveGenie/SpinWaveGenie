#ifndef __Initializer__
#define __Initializer__

#include <iostream>
#include <sstream>
#include <new>
#include <functional>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/trim.hpp>
#include "extern/pugixml.hpp"
#include "Cell.h"
#include "SW_Builder.h"

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
    SW_Builder get_builder();
    
private:
    std::map< std::string, std::function< void (const pugi::xml_node&) > > m_parser_map; //!< Map for dispatching parsers based on node type
    Cell unit_cell;
    SW_Builder builder;
    void parseCrystalNode(const pugi::xml_node &node);
    void parseInteractionNode(const pugi::xml_node &node);
    void parseKpointsNode(const pugi::xml_node &node);
};

#endif /* defined(__Initializer__) */
