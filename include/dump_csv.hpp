#pragma once
// GPL v3 License
// Copyright 2023--present Flowy developers
#include <exception>
#include <istream>
#include <iterator>
#include <sstream>
#include <string>
#include <utility>

#include <xtensor/xtensor.hpp>
#include <xtensor/xtensor_config.hpp>

namespace Flowy::Utility
{
/**
 * @brief Dump tensor to CSV. This is just a copy of the xtensor::dump_csv,
 * but with the ability to specify a custom delimiter
 *
 * @param stream the output stream to write the CSV encoded values
 * @param e the tensor expression to serialize
 * @param delimiter the delimiter, defaults to ','
 */
template<class E>
void dump_csv( std::ostream & stream, const xt::xexpression<E> & e, const char & delimiter = ',' )
{
    using size_type = typename E::size_type;
    const E & ex    = e.derived_cast();
    if( ex.dimension() != 2 )
    {
        XTENSOR_THROW( std::runtime_error, "Only 2-D expressions can be serialized to CSV" );
    }
    size_type nbrows = ex.shape()[0], nbcols = ex.shape()[1];
    auto st = ex.stepper_begin( ex.shape() );
    for( size_type r = 0; r != nbrows; ++r )
    {
        for( size_type c = 0; c != nbcols; ++c )
        {
            stream << *st;
            if( c != nbcols - 1 )
            {
                st.step( 1 );
                stream << delimiter;
            }
            else
            {
                st.reset( 1 );
                st.step( 0 );
                stream << std::endl;
            }
        }
    }
}
} // namespace Flowy::Utility
