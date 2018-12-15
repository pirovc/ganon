#include <ganon-build/Config.hpp>
#include <ganon-build/GanonBuild.hpp>

#include <catch2/catch.hpp>

#include <fstream>
#include <iterator>
#include <streambuf>
#include <string>

namespace detail
{

bool filesAreEqual( const std::string& a, const std::string& b )
{
    std::ifstream stream{ a };
    std::string   file1{ std::istreambuf_iterator< char >( stream ), std::istreambuf_iterator< char >() };

    stream = std::ifstream{ b };
    std::string file2{ std::istreambuf_iterator< char >( stream ), std::istreambuf_iterator< char >() };

    return file1 == file2;
}

} // namespace detail

SCENARIO( "ganon-build test" )
{
    Config config;
    config.seqid_bin_file     = "bacteria_seqid_bin.txt";
    config.output_filter_file = "test_output.filter";
    config.filter_size        = 15797760;
    config.kmer_size          = 19;
    config.hash_functions     = 3;
    config.reference_files    = {
        "bacteria_NC_010333.1.fasta.gz",
        "bacteria_NC_017163.1.fasta.gz",
        "bacteria_NC_017164.1.fasta.gz",
        "bacteria_NC_017543.1.fasta.gz"
    };
    config.threads       = 1;
    config.build_threads = 1;

    REQUIRE( GanonBuild::run( config ) );
    REQUIRE( detail::filesAreEqual( config.output_filter_file, "build_output.filter" ) );
}
