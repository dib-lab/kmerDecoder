#include "doctest.h"
#include <string>

std::string test_func(){
    return "pending..";
}

TEST_CASE("testing the factorial function") {
CHECK(test_func() == "pending..");
}