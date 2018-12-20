#pragma once

#include <fstream>
#include <iterator>
#include <streambuf>
#include <string>

namespace aux
{

inline bool filesAreEqual( const std::string& file1, const std::string& file2 )
{
    std::ifstream     stream1{ file1 };
    const std::string data1{ std::istreambuf_iterator< char >{ stream1 }, std::istreambuf_iterator< char >{} };

    std::ifstream     stream2{ file2 };
    const std::string data2{ std::istreambuf_iterator< char >{ stream2 }, std::istreambuf_iterator< char >{} };

    return data1 == data2;
}

inline bool fileIsEmpty( const std::string& file )
{
    std::ifstream stream{ file };
    return stream.peek() == std::ifstream::traits_type::eof();
}

} // namespace aux
