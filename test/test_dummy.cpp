#include "flowtastic.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>

TEST_CASE( "DummyT3st", "[dummy]" )
{
    REQUIRE( Flowtastic::i_am_true() );
}