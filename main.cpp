#include <iostream>
#include "joint_tap.h"

int main() {
    auto [graphA, graphT, assistance] = read_file_joint("../tests/graphs.txt");
    JointTAP solver(graphA, graphT, assistance);
    cout << solver.solve() << endl;
    return 0;
}