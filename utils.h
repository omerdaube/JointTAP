#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <climits>
using namespace std;

const int NO_COMPUTE = INT_MAX;
const double EPS = 1e-6;

class TimeRewardPair {
public:
    double time;
    double reward;

    TimeRewardPair(double time, double reward) {
        this->time = time;
        this->reward = reward;
    }

    TimeRewardPair(const pair<double, double>& p) {
        this->time = p.first;
        this->reward = p.second;
    }

    bool operator==(const TimeRewardPair& other) const {
        return time == other.time && reward == other.reward;
    }
    bool operator<(const TimeRewardPair& other) const {
        if (time < other.time)
            return true;
        if (time > other.time)
            return false;
        return reward < other.reward;
    }
};

typedef vector<TimeRewardPair> ArrivalsList;

class Interval {
public:
    double start;
    double end;
    ArrivalsList arrivals;
    Interval(double start, double end, ArrivalsList arrivals)
            : start(start), end(end), arrivals(std::move(arrivals)) {}
    Interval(double start, double end)
            : start(start), end(end), arrivals() {}

    Interval(const Interval&) = default;
    Interval& operator=(const Interval&) = default;
};

class Edge {
public:
    int u;
    int v;
    double l;

    bool operator==(const Edge& other) const {
        return u == other.u && v == other.v && l == other.l;
    }

    Edge(int u, int v, double l)
            : u(u), v(v), l(l) {}
};

typedef vector<Edge> Path;
typedef vector<Interval> VertexIntervalList;
typedef map<int, VertexIntervalList> Intervals;

typedef vector<TimeRewardPair> PairList;
typedef set<double> TimeSet;




class AdjEntry {
public:
    int v;
    double l;

    AdjEntry(int v, double l) : v(v), l(l) {}
};
typedef vector<vector<AdjEntry>> Graph;
struct pair_hash {
    size_t operator()(const pair<int, int>& p) const {
        return hash<int>()(p.first) ^ (hash<int>()(p.second) << 1);
    }
};
using PairsSet = unordered_set<std::pair<int, int>, pair_hash>;
const double INF = 1e15;


pair<Intervals, Path> read_file_otp(string path_file, bool sim);
tuple<Graph, Graph, PairsSet> read_file_joint(string path_file);
tuple<Path, Graph, PairsSet> read_file_task(string path_file);
tuple<Graph, Intervals> read_file_assistance(string path_file);

class ValueHistory {
    vector<vector<double>> value_history;

public:
    explicit ValueHistory(size_t graph_size) {
        for (int i = 0; i < graph_size; i++) {
            value_history.push_back(vector<double>({-2 * EPS}));
        }
    }

    double last(int u) {
        return value_history[u].back();
    }

    void add(int u, double value) {
        value_history[u].push_back(value);
    }

    void removeLast(int u) {
        value_history[u].pop_back();
    }

};

class PathHash {
public:
    size_t operator()(const Path& path) const {
        size_t hash = 0;
        std::hash<int> intHasher;
        std::hash<double> doubleHasher;
        for (const auto& edge : path) {
            hash ^= intHasher(edge.u) ^ (intHasher(edge.v) << 1) ^ (doubleHasher(edge.l) << 2);
        }
        return hash;
    }
};

class PathEqual {
public:
    bool operator()(const Path& path1, const Path& path2) const {
        return path1 == path2;
    }
};

using TimePair = pair<double,double>; // (entry, exit)
using TimingList = vector<TimePair>;
using Config = std::vector<std::vector<double>>;
std::tuple<Config, Config, double> loadConfigs(const std::string &filepath);
void write_configs_to_json(const Config& assistance,
                           const Config& task,
                           const std::string& filepath);
#endif // UTILS_H
