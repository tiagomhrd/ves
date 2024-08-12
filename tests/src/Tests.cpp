#include <iostream>
#include "ves.h"
#include <catch_amalgamated.hpp>

TEST_CASE("Truth") {
    REQUIRE(true);
}

int main(int argc, char* argv[]) {
	int result = Catch::Session().run(argc, argv);
	return result;
}