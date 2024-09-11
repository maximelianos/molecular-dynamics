// Include a library file to make sure proper includes are set
#include <gtest/gtest.h>
#include "Eigen/Core"
#include "hello.h"

// Demonstrate some basic assertions.
TEST(HelloTest, BasicAssertions) {
    // Expect two strings not to be equal.
    EXPECT_STRNE("hello", "world");
    // Expect equality.
    EXPECT_EQ(7 * 6, 42);
    // Testing if we can call a function from our MD library
    hello_eigen();
}

TEST(HelloTest, IntToString) {
    int x = 10;
    std::string s;
    s = std::to_string(x);
    EXPECT_EQ(s, "10");
}

TEST(HelloTest, Norm3D) {
    Eigen::Vector3d a;
    a << 1, 0, 0;
    Eigen::Vector3d b;
    b << 0, 1, 0;
    // norm = sqrt(1^2 + 1^2 + 0) = sqrt(2)
    ASSERT_FLOAT_EQ((a - b).norm(), std::sqrt(2));
}


