#pragma once

typedef robin_hood::unordered_map< uint32_t, std::string > TMap;

inline TMap load_map( std::string map_file )
{
    TMap          map;
    std::string   line;
    std::ifstream infile;
    infile.open( map_file );
    while ( std::getline( infile, line, '\n' ) )
    {
        std::istringstream         stream_line( line );
        std::vector< std::string > fields;
        std::string                field;
        while ( std::getline( stream_line, field, '\t' ) )
            fields.push_back( field );
        // target <tab> binid
        map[std::stoul( fields[1] )] = fields[0];
    }
    infile.close();
    return map;
}
