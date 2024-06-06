#include "hello.h"
#include <sys/param.h>
#include <Eigen/Dense>
#include <iostream>

using Eigen::MatrixXd;

std::string get_working_path()
{
    char temp[MAXPATHLEN];
    return ( getcwd(temp, sizeof(temp)) ? std::string( temp ) : std::string("") );
}

void hello_eigen() {
    std::cout << get_working_path() << std::endl;

    MatrixXd m(2, 2);
    m(0, 0) = 3;
    m(1, 0) = 2.5;
    m(0, 1) = -1;
    m(1, 1) = m(1, 0) + m(0, 1);

    // C++11 style loop over elements of an Eigen matrix
    for (auto&& value : m.reshaped())
        value += 1.;

    // Output matrix values
    std::cout << m << std::endl;
}
