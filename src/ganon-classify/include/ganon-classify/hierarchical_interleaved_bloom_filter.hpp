// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides raptor::hierarchical_interleaved_bloom_filter.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <ranges>

#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>

//#ifndef RAPTOR_HIBF_HAS_COUNT
//#define RAPTOR_HIBF_HAS_COUNT 0
//#endif
#define RAPTOR_HIBF_HAS_COUNT 1

namespace raptor
{

/*!\brief The HIBF binning directory. A data structure that efficiently answers set-membership queries for multiple
 *        bins.
 * \tparam data_layout_mode_ Indicates whether the underlying data type is compressed. See
 *                           [seqan3::data_layout](https://docs.seqan.de/seqan/3.0.3/group__submodule__dream__index.html#gae9cb143481c46a1774b3cdf5d9fdb518).
 * \see [seqan3::interleaved_bloom_filter][1]
 * \details
 *
 * This class improves the [seqan3::interleaved_bloom_filter][1] by adding additional bookkeeping that allows
 * to establish a hierarchical structure. This structure can then be used to split or merge user bins and distribute
 * them over a variable number of technical bins. In the [seqan3::interleaved_bloom_filter][1], the number of user bins
 * and technical bins is always the same. This causes performance degradation when there are many user bins or the user
 * bins are unevenly distributed.
 *
 * # Terminology
 *
 * ## User Bin
 * The user may impose a structure on his sequence data in the form of logical groups (e.g. species). When querying the
 * (H)IBF, the user is interested in an answer that differentiates between these groups.
 *
 * ## Technical Bin
 * A Technical Bin represents an actual bin in the binning directory. In the IBF, it stores its kmers in a single Bloom
 * Filter (which is interleaved with all the other BFs). In the HIBF each of these bins could be merged or splitted,
 * thus it differs to the original user bins.
 *
 * ## Layout of the HIBF
 * The relationship between user bins and technical bins in the HIBF, i.e. which are split or merged, how substructures
 * of merged bins (lower-level IBFs) look like is called layout.
 *
 * # Hierarchical Interleaved Bloom Filter (HIBF)
 *
 * In constrast to the [seqan3::interleaved_bloom_filter][1], the user bins may be split across multiple technical bins,
 * or multiple user bins may be merged into one technical bin. When merging multiple user bins, the HIBF stores another
 * IBF that is built over the user bins constituting the merged bin. This lower-level IBF can then be used to further
 * distinguish between merged bins.
 *
 * In this example layout, user bin 1 was split into two technical bins. Bins 3, 4, and 5 were merged into a single
 * technical bin, and another IBF was added for the merged bin.
 * \image html hibf.svg width=40%
 *
 * The individual IBFs may have a different number of technical bins and differ in their sizes, allowing an efficient
 * distribution of the user bins.
 *
 * ## Querying
 * To query the Hierarchical Interleaved Bloom Filter for values, call
 * raptor::hierarchical_interleaved_bloom_filter::membership_agent() and use the returned
 * raptor::hierarchical_interleaved_bloom_filter::membership_agent.
 * In contrast to the [seqan3::interleaved_bloom_filter][1], the result will consist of indices of user bins.
 *
 * To count the occurrences in each user bin of a range of values in the Hierarchical Interleaved Bloom Filter, call
 * raptor::hierarchical_interleaved_bloom_filter::counting_agent() and use
 * the returned raptor::hierarchical_interleaved_bloom_filter::counting_agent_type.
 *
 * ## Thread safety
 *
 * The Interleaved Bloom Filter promises the basic thread-safety by the STL that all
 * calls to `const` member functions are safe from multiple threads (as long as no thread calls
 * a non-`const` member function at the same time).
 *
 * [1]: https://docs.seqan.de/seqan/3.0.3/classseqan3_1_1interleaved__bloom__filter.html
 */
template < seqan3::data_layout data_layout_mode_ = seqan3::data_layout::uncompressed >
class hierarchical_interleaved_bloom_filter
{
public:
    // Forward declaration
    class user_bins;

    // Forward declaration
    class membership_agent;

#if RAPTOR_HIBF_HAS_COUNT
    // Forward declaration
    template < std::integral value_t >
    class counting_agent_type;
#endif

    //!\brief Indicates whether the Interleaved Bloom Filter is compressed.
    static constexpr seqan3::data_layout data_layout_mode = data_layout_mode_;

    //!\brief The type of an individual Bloom filter.
    using ibf_t = seqan3::interleaved_bloom_filter< data_layout_mode_ >;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    hierarchical_interleaved_bloom_filter()                                               = default; //!< Defaulted.
    hierarchical_interleaved_bloom_filter( hierarchical_interleaved_bloom_filter const& ) = default; //!< Defaulted.
    hierarchical_interleaved_bloom_filter& operator=( hierarchical_interleaved_bloom_filter const& ) =
        default;                                                                                //!< Defaulted.
    hierarchical_interleaved_bloom_filter( hierarchical_interleaved_bloom_filter&& ) = default; //!< Defaulted.
    hierarchical_interleaved_bloom_filter& operator=( hierarchical_interleaved_bloom_filter&& ) =
        default;                                        //!< Defaulted.
    ~hierarchical_interleaved_bloom_filter() = default; //!< Defaulted.

    //!\}

    //!\brief The individual interleaved Bloom filters.
    std::vector< ibf_t > ibf_vector;

    /*!\brief Stores for each bin in each IBF of the HIBF the ID of the next IBF.
     * \details
     * Assume we look up a bin `b` in IBF `i`, i.e. `next_ibf_id[i][b]`.
     * If `i` is returned, there is no lower level IBF, bin `b` is hence not a merged bin.
     * If `j != i` is returned, there is a lower level IBF, bin `b` is a merged bin, and `j` is the ID of the lower
     * level IBF in ibf_vector.
     */
    std::vector< std::vector< int64_t > > next_ibf_id;

    //!\brief The underlying user bins.
    user_bins user_bins;

    //!\brief Returns a membership_agent to be used for counting.
    membership_agent membership_agent() const
    {
        return typename hierarchical_interleaved_bloom_filter< data_layout_mode >::membership_agent{ *this };
    }

#if RAPTOR_HIBF_HAS_COUNT
    /*!\brief Returns a counting_agent_type to be used for counting.
     * \tparam value_t The type to use for the counters; must model std::integral.
     */
    template < std::integral value_t = uint16_t >
    counting_agent_type< value_t > counting_agent() const
    {
        return counting_agent_type< value_t >{ *this };
    }
#endif

    /*!\cond DEV
     * \brief Serialisation support function.
     * \tparam archive_t Type of `archive`; must satisfy seqan3::cereal_archive.
     * \param[in] archive The archive being serialised from/to.
     *
     * \attention These functions are never called directly.
     * \sa https://docs.seqan.de/seqan/3.2.0/group__io.html#serialisation
     */
    template < seqan3::cereal_archive archive_t >
    void CEREAL_SERIALIZE_FUNCTION_NAME( archive_t& archive )
    {
        archive( ibf_vector );
        archive( next_ibf_id );
        archive( user_bins );
    }
    //!\endcond
};

/*!\brief Bookkeeping for user and technical bins.
 */
template < seqan3::data_layout data_layout_mode >
class hierarchical_interleaved_bloom_filter< data_layout_mode >::user_bins
{
private:
    //!\brief Contains filenames of all user bins.
    std::vector< std::string > user_bin_filenames;

    /*!\brief Stores for each bin in each IBF of the HIBF the ID of the filename.
     * \details
     * Assume we look up a bin `b` in IBF `i`, i.e. `ibf_bin_to_filename_position[i][b]`.
     * If `-1` is returned, bin `b` is a merged bin, and there is no filename, we need to look into the lower level IBF.
     * Otherwise, the returned value `j` can be used to access the corresponding filename `user_bin_filenames[j]`.
     */
    std::vector< std::vector< int64_t > > ibf_bin_to_filename_position{};

public:
    //!\brief Returns the number of managed user bins.
    size_t num_user_bins() const noexcept
    {
        return user_bin_filenames.size();
    }

    //!\brief Changes the number of managed IBFs.
    void set_ibf_count( size_t const size )
    {
        ibf_bin_to_filename_position.resize( size );
    }

    //!\brief Changes the number of managed user bins.
    void set_user_bin_count( size_t const size )
    {
        user_bin_filenames.resize( size );
    }

    /*!\brief Returns a vector containing user bin indices for each bin in the `idx`th IBF.
     * \param idx The id of the x-th IBF.
     *
     * \details
     *
     * ### Example
     *
     * \include test/snippet/hibf/bin_indices_of_ibf.cpp
     */
    std::vector< int64_t >& bin_indices_of_ibf( size_t const idx )
    {
        return ibf_bin_to_filename_position[idx];
    }

    /*!\brief Returns the filename of the `idx`th user bin.
     * \param idx The id of the x-th user bin.
     *
     * \details
     *
     * ### Example
     *
     * \include test/snippet/hibf/filename_of_user_bin.cpp
     */
    std::string& filename_of_user_bin( size_t const idx )
    {
        return user_bin_filenames[idx];
    }

    //!\brief For a pair `(a,b)`, returns a const reference to the filename of the user bin at IBF `a`, bin `b`.
    std::string const& operator[]( std::pair< size_t, size_t > const& index_pair ) const
    {
        return user_bin_filenames[ibf_bin_to_filename_position[index_pair.first][index_pair.second]];
    }

    /*!\brief Returns a view over the user bin filenames for the `ibf_idx`th IBF.
     *        An empty string is returned for merged bins.
     */
    auto operator[]( size_t const ibf_idx ) const
    {
        return ibf_bin_to_filename_position[ibf_idx] | std::views::transform( [this]( int64_t i ) {
                   if ( i == -1 )
                       return std::string{};
                   else
                       return user_bin_filenames[i];
               } );
    }

    //!\brief Returns the filename index of the `ibf_idx`th IBF for bin `bin_idx`.
    int64_t filename_index( size_t const ibf_idx, size_t const bin_idx ) const
    {
        return ibf_bin_to_filename_position[ibf_idx][bin_idx];
    }

    /*!\brief Writes all filenames to a stream. Index and filename are tab-separated.
     * \details
     * 0	\<path_to_user_bin_0\>
     * 1	\<path_to_user_bin_1\>
     */
    template < typename stream_t >
    void write_filenames( stream_t& out_stream ) const
    {
        size_t      position{};
        std::string line{};
        for ( auto const& filename : user_bin_filenames )
        {
            line.clear();
            line = '#';
            line += std::to_string( position );
            line += '\t';
            line += filename;
            line += '\n';
            out_stream << line;
            ++position;
        }
    }

    /*!\cond DEV
     * \brief Serialisation support function.
     * \tparam archive_t Type of `archive`; must satisfy seqan3::cereal_archive.
     * \param[in] archive The archive being serialised from/to.
     *
     * \attention These functions are never called directly.
     * \sa https://docs.seqan.de/seqan/3.2.0/group__io.html#serialisation
     */
    template < typename archive_t >
    void serialize( archive_t& archive )
    {
        archive( user_bin_filenames );
        archive( ibf_bin_to_filename_position );
    }
    //!\endcond
};

/*!\brief Manages membership queries for the raptor::hierarchical_interleaved_bloom_filter.
 * \see raptor::hierarchical_interleaved_bloom_filter::user_bins::filename_of_user_bin
 * \details
 * In contrast to the [seqan3::interleaved_bloom_filter][1], the result will consist of indices of user bins.
 */
template < seqan3::data_layout data_layout_mode > // TODO: value_t as template?
class hierarchical_interleaved_bloom_filter< data_layout_mode >::membership_agent
{
private:
    //!\brief The type of the augmented hierarchical_interleaved_bloom_filter.
    using hibf_t = hierarchical_interleaved_bloom_filter< data_layout_mode >;

    //!\brief A pointer to the augmented hierarchical_interleaved_bloom_filter.
    hibf_t const* const hibf_ptr{ nullptr };

    //!\brief Helper for recursive membership querying.
    template < std::ranges::forward_range value_range_t >
    void bulk_contains_impl( value_range_t&& values, int64_t const ibf_idx, size_t const threshold )
    {
        auto  agent  = hibf_ptr->ibf_vector[ibf_idx].template counting_agent< uint16_t >();
        auto& result = agent.bulk_count( values );

        uint16_t sum{};

        for ( size_t bin{}; bin < result.size(); ++bin )
        {
            sum += result[bin];

            auto const current_filename_index = hibf_ptr->user_bins.filename_index( ibf_idx, bin );

            if ( current_filename_index < 0 ) // merged bin
            {
                if ( sum >= threshold )
                    bulk_contains_impl( values, hibf_ptr->next_ibf_id[ibf_idx][bin], threshold );
                sum = 0u;
            }
            else if ( bin + 1u == result.size() || // last bin
                      current_filename_index
                          != hibf_ptr->user_bins.filename_index( ibf_idx, bin + 1 ) ) // end of split bin
            {
                if ( sum >= threshold )
                    result_buffer.emplace_back( current_filename_index );
                sum = 0u;
            }
        }
    }

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    membership_agent()                          = default;            //!< Defaulted.
    membership_agent( membership_agent const& ) = default;            //!< Defaulted.
    membership_agent& operator=( membership_agent const& ) = default; //!< Defaulted.
    membership_agent( membership_agent&& )                 = default; //!< Defaulted.
    membership_agent& operator=( membership_agent&& ) = default;      //!< Defaulted.
    ~membership_agent()                               = default;      //!< Defaulted.

    /*!\brief Construct a membership_agent for an existing hierarchical_interleaved_bloom_filter.
     * \private
     * \param hibf The hierarchical_interleaved_bloom_filter.
     */
    explicit membership_agent( hibf_t const& hibf )
    : hibf_ptr( std::addressof( hibf ) )
    {
    }
    //!\}

    //!\brief Stores the result of bulk_contains().
    std::vector< int64_t > result_buffer;

    /*!\name Lookup
     * \{
     */
    /*!\brief Determines set membership of given values, and returns the user bin indices of occurrences.
     * \param[in] values The values to process; must model std::ranges::forward_range.
     * \param[in] threshold Report a user bin if there are at least this many hits.
     *
     * \attention The result of this function must always be bound via reference, e.g. `auto &`, to prevent copying.
     * \attention Sequential calls to this function invalidate the previously returned reference.
     *
     * \details
     *
     * ### Thread safety
     *
     * Concurrent invocations of this function are not thread safe, please create a
     * raptor::hierarchical_interleaved_bloom_filter::membership_agent for each thread.
     */
    template < std::ranges::forward_range value_range_t >
    [[nodiscard]] std::vector< int64_t > const& bulk_contains( value_range_t&& values,
                                                               size_t const    threshold ) & noexcept
    {
        assert( hibf_ptr != nullptr );

        static_assert( std::ranges::forward_range< value_range_t >, "The values must model forward_range." );
        static_assert( std::unsigned_integral< std::ranges::range_value_t< value_range_t > >,
                       "An individual value must be an unsigned integral." );

        result_buffer.clear();

        bulk_contains_impl( values, 0, threshold );

        std::ranges::sort( result_buffer ); // TODO: necessary?

        return result_buffer;
    }

    // `bulk_contains` cannot be called on a temporary, since the object the returned reference points to
    // is immediately destroyed.
    template < std::ranges::range value_range_t >
    [[nodiscard]] std::vector< int64_t > const& bulk_contains( value_range_t&& values,
                                                               size_t const    threshold ) && noexcept = delete;
    //!\}
};

#if RAPTOR_HIBF_HAS_COUNT
/*!\brief Manages counting ranges of values for the raptor::hierarchical_interleaved_bloom_filter.
 */
template < seqan3::data_layout data_layout_mode >
template < std::integral value_t >
class hierarchical_interleaved_bloom_filter< data_layout_mode >::counting_agent_type
{
private:
    //!\brief The type of the augmented hierarchical_interleaved_bloom_filter.
    using hibf_t = hierarchical_interleaved_bloom_filter< data_layout_mode >;

    //!\brief A pointer to the augmented hierarchical_interleaved_bloom_filter.
    hibf_t const* const hibf_ptr{ nullptr };

    //!\brief Helper for recursive bulk counting.
    template < std::ranges::forward_range value_range_t >
    void bulk_count_impl( value_range_t&& values, int64_t const ibf_idx, size_t const threshold )
    {
        auto  agent  = hibf_ptr->ibf_vector[ibf_idx].template counting_agent< value_t >();
        auto& result = agent.bulk_count( values );

        value_t sum{};

        for ( size_t bin{}; bin < result.size(); ++bin )
        {
            sum += result[bin];
            auto const current_filename_index = hibf_ptr->user_bins.filename_index( ibf_idx, bin );

            if ( current_filename_index < 0 ) // merged bin
            {
                if ( sum >= threshold )
                    bulk_count_impl( values, hibf_ptr->next_ibf_id[ibf_idx][bin], threshold );
                sum = 0u;
            }
            else if ( bin + 1u == result.size() || // last bin
                      current_filename_index
                          != hibf_ptr->user_bins.filename_index( ibf_idx, bin + 1 ) ) // end of split bin
            {
                if ( sum >= threshold )
                    result_buffer[current_filename_index] = sum;
                sum = 0u;
            }
        }
    }

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    counting_agent_type()                             = default;            //!< Defaulted.
    counting_agent_type( counting_agent_type const& ) = default;            //!< Defaulted.
    counting_agent_type& operator=( counting_agent_type const& ) = default; //!< Defaulted.
    counting_agent_type( counting_agent_type&& )                 = default; //!< Defaulted.
    counting_agent_type& operator=( counting_agent_type&& ) = default;      //!< Defaulted.
    ~counting_agent_type()                                  = default;      //!< Defaulted.

    /*!\brief Construct a counting_agent_type for an existing hierarchical_interleaved_bloom_filter.
     * \private
     * \param hibf The hierarchical_interleaved_bloom_filter.
     */
    explicit counting_agent_type( hibf_t const& hibf )
    : hibf_ptr( std::addressof( hibf ) )
    , result_buffer( hibf_ptr->user_bins.num_user_bins() )
    {
    }
    //!\}

    //!\brief Stores the result of bulk_count().
    seqan3::counting_vector< value_t > result_buffer;

    /*!\name Counting
     * \{
     */
    /*!\brief Counts the occurrences in each bin for all values in a range.
     * \tparam value_range_t The type of the range of values. Must model std::ranges::forward_range. The reference type
     *                       must model std::unsigned_integral.
     * \param[in] values The range of values to process.
     * \param[in] threshold Do not recurse into merged bins with less than this many hits. Default: 1.
     *
     * \attention The result of this function must always be bound via reference, e.g. `auto &`, to prevent copying.
     * \attention Sequential calls to this function invalidate the previously returned reference.
     *
     * \details
     *
     * ### Thread safety
     *
     * Concurrent invocations of this function are not thread safe, please create a
     * raptor::hierarchical_interleaved_bloom_filter::counting_agent_type for each thread.
     */
    template < std::ranges::forward_range value_range_t >
    [[nodiscard]] seqan3::counting_vector< value_t > const& bulk_count( value_range_t&& values,
                                                                        size_t const    threshold = 1u ) & noexcept
    {
        assert( hibf_ptr != nullptr );
        assert( threshold > 0u );
        assert( result_buffer.size() == hibf_ptr->user_bins.num_user_bins() );

        static_assert( std::ranges::forward_range< value_range_t >, "The values must model forward_range." );
        static_assert( std::unsigned_integral< std::ranges::range_value_t< value_range_t > >,
                       "An individual value must be an unsigned integral." );

        std::ranges::fill( result_buffer, static_cast< value_t >( 0u ) );

        bulk_count_impl( values, 0, threshold );

        return result_buffer;
    }

    // `bulk_count` cannot be called on a temporary, since the object the returned reference points to
    // is immediately destroyed.
    template < std::ranges::range value_range_t >
    [[nodiscard]] seqan3::counting_vector< value_t > const& bulk_count(
        value_range_t&& values, size_t const threshold = 1u ) && noexcept = delete;
    //!\}
};
#endif // RAPTOR_HIBF_HAS_COUNT

} // namespace raptor
