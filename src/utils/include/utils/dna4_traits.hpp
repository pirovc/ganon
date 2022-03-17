#pragma once
#include <seqan3/io/sequence_file/input.hpp>
struct dna4_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};
