#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <cmath>
#include "utils.h"
using namespace std;

tuple<Graph, Graph, PairsSet> read_file_joint(string path_file) {
    Graph graphA;
    Graph graphT;
    PairsSet ps;

    ifstream file(path_file);
    if (!file.is_open()) {
        cout << "File not found: " << path_file << endl;
        exit(1);
    }
    string line;

    while (std::getline(file, line) && line != "-") {
        std::istringstream iss(line);
        int u, len;
        iss >> u >> len;
        vector<AdjEntry> adjList;
        for (int i = 0; i < len; i++) {
            int v;
            double l;
            iss >> v >> l;
            adjList.emplace_back(v, l);
        }
        graphA.push_back(adjList);
    }

    while (std::getline(file, line) && line != "-") {
        std::istringstream iss(line);
        int u, len;
        iss >> u >> len;
        vector<AdjEntry> adjList;
        for (int i = 0; i < len; i++) {
            int v;
            double l;
            iss >> v >> l;
            adjList.emplace_back(v, l);
        }
        graphT.push_back(adjList);
    }

    while (std::getline(file, line) && line != "-") {
        std::istringstream iss(line);
        int u, v;
        iss >> u >> v;
        ps.insert({u,v});
    }

    return {graphA, graphT, ps};
}