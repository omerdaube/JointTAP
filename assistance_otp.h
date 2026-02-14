#ifndef DTAP_CPP_ASSISTANCE_OTP_H
#define DTAP_CPP_ASSISTANCE_OTP_H


#include "utils.h"
#include "incremental_otp.h"
#include "lp_bound.h"
#include <iostream>
#include <map>
#include <algorithm>
#include <chrono>

typedef vector<vector<double>> Matrix; // Represents a 2D matrix of doubles.

/**
 * Function to find the shortest path from a starting node to all other nodes in the graph using Dijkstra's algorithm.
 * @param graph The graph to perform the algorithm on.
 * @param start The starting node.
 * @param dist The vector to store the shortest distances.
 */
void dijkstra(Graph& graph, int start, vector<double>& dist);

/**
 * Represents an interval associated with a vertex.
 */
class VInterval : public Interval {
public:
    int vertex; // The vertex associated with the interval.

    /**
     * Constructor initializing the vertex and the interval with a start and end time.
     * @param vertex The vertex associated with the interval.
     * @param start The start time of the interval.
     * @param end The end time of the interval.
     */
    VInterval(int vertex, double start, double end) : Interval(start, end) {
        this->vertex = vertex;
    }

    /**
     * Constructor initializing the vertex and the interval with an existing interval.
     * @param vertex The vertex associated with the interval.
     * @param interval The existing interval.
     */
    VInterval(int vertex, Interval interval) : Interval(interval.start, interval.end) {
        this->vertex = vertex;
    }

    bool operator<(const VInterval& other) const {
        // is smaller if i end first
        if (end<other.end)
            return true;
        if (end > other.end)
            return false;
        // if we end at the same time, im smaller if i start first
        if (start < other.start)
            return true;
        if (start > other.start)
            return false;
        // if we start and end at the same time, im smaller if i am at a smaller vertex
        return vertex < other.vertex;
    }
};

/**
 * A class containing a method to bound the reward obtainable starting from a given node at a given time.
 */
class GraphBound {
public:
    Graph graph; // The graph associated with this bound.
    Intervals intervals; // The intervals associated with this bound.

    Matrix min_dist; // The minimum distances between all pairs of nodes excluding half of the first and last edge.
    Matrix min_dist_no_last; // The minimum distances between all pairs of nodes excluding half of the last edge.

    map<VInterval, double> UBi; // a mapping from an interval to its corresponding UB_I value.

    /**
     * Constructor initializing the graph and intervals and computing all other variables.
     * @param graph The graph to initialize with.
     * @param intervals The intervals to initialize with.
     */
    GraphBound(Graph& graph, Intervals& intervals);

    /**
     * Function to calculate the upper bound for a given node at a given time [UB(u,t)].
     * @param u The node to calculate the upper bound for.
     * @param time The time at which to calculate the upper bound.
     * @return The calculated upper bound.
     */
    double UB(int u, double time);

private:
    /**
     * Function to compute the minimum distances between all pairs of nodes.
     */
    void computeDists();

    /**
     * Function to compute UB_I.
     */
    void computeUBI();
};

typedef set<TimeRewardPair> PairSet; // Represents a set of time-reward pairs.


/**
 * Wraps the GraphBound class to implement UB on a path.
 * Contains methods for caching rewards and determining whether a node should be pruned.
 */
class GraphBoundWrapper {
    double max_inner_UB;
    GraphBound bound; // The GraphBound object to wrap.
    map<int, PairSet> cache; // A cache of UB(u,t) values

    /**
     * Function to add a reward for a given node at a given time to the cache.
     * @param u The node to add the reward for.
     * @param time The time at which to add the reward.
     * @param reward The reward to add.
     */
    void addToCache(int u, double time, double reward);

    /**
     * Function to retrieve a reward for a given node at a given time from the cache.
     * @param u The node to retrieve the reward for.
     * @param time The time at which to retrieve the reward.
     * @return A pair containing the retrieved reward and a boolean indicating whether the reward was found in the cache.
     */
    pair<double, bool> getFromCache(int u, double time);

public:
    double get_max_inner_UB();

    /**
     * Constructor initializing the graph and intervals.
     * @param graph The graph to initialize with.
     * @param intervals The intervals to initialize with.
     */
    GraphBoundWrapper(Graph& graph, Intervals& intervals);

    /**
     * Function to determine whether a node should be pruned based on UB(pi).
     * @param u The node to check for pruning.
     * @param bound_pairs The list of bound pairs computed by the updatedOTP algorithm.
     * @param max_reward The maximum reward obtained until now.
     * @return True if the node should be pruned, false otherwise.
     */
    bool shouldPrune(int u, PairList& bound_pairs, double max_reward);

};

class AbstractAssistanceOTP {
protected:
    Graph& graph; // The graph to solve the AbstractAssistanceOTP for.
    Intervals intervals; // The intervals for each vertex in the graph.

    double max_value; // The maximal value found by the AbstractAssistanceOTP solver.
    Path max_path; // The path resulting in the maximal value found by the AbstractAssistanceOTP solver.

    int path_count; // The number of paths evaluated by the AbstractAssistanceOTP solver.

    std::chrono::high_resolution_clock::time_point start_time; // The start time for the AbstractAssistanceOTP solver.
    const double time_limit = 1800; // The time limit for the AbstractAssistanceOTP solver.


public:
    /**
     * Constructs a new AbstractAssistanceOTP object.
     *
     * @param graph The graph to solve the AbstractAssistanceOTP for.
     * @param intervals The intervals for each vertex in the graph.
     */
    AbstractAssistanceOTP(Graph& graph, Intervals& intervals) : graph(graph), intervals(intervals), max_value(0), path_count(0) {};

    /**
     * Runs the AbstractAssistanceOTP solver.
     *
     * @return The maximum reward obtained by the AbstractAssistanceOTP solver.
     */
    virtual double solve() = 0;

    /**
     * Prints the best path found by the AbstractAssistanceOTP solver.
     */
    void printBestPath() {
        if (max_path.empty()) {
            cout << "No path found" << endl;
            cout << "path: " << path_count << endl;
            return;
        }
        cout << max_path.front().u;
        for (auto& edge: max_path) {
            cout << " -> " << edge.v;
        }
        cout << endl;
        cout << "path: " << path_count << endl;
    }

    /**
     * Returns the number of paths evaluated by the AbstractAssistanceOTP solver.
     *
     * @return The number of paths evaluated by the AbstractAssistanceOTP solver.
     */
    int getPathCount() {
        return path_count;
    }

    /**
         * Returns the best path found by the AbstractAssistanceOTP solver.
         *
         * @return The best path found by the AbstractAssistanceOTP solver.
         */
    Path get_best_path() {
        return max_path;
    }
};

/**
 * A Branch and Bound solver for the AbstractAssistanceOTP problem
 */
class AssistanceOTP : public AbstractAssistanceOTP {
    vector<double> min_time_vertex;
    PathSolver path_solver;
    double eps; // The epsilon value for the AssistanceOTP solver.
    double max_interval_size; // The maximum interval size for the AssistanceOTP solver.
    double inner_UB;

    GraphBoundWrapper bound; // The upper bound for the AssistanceOTP solver.
    ValueHistory value_history; // The value history for the AssistanceOTP solver.

    Intervals added_intervals;

    int dummy_vertex; // A dummy vertex.


    /**
     * Preprocesses the intervals.
     */
    void preprocess_intervals();

    /**
     * Preprocesses the graph.
     */
    void preprocess_graph();

    /**
     * Runs the branch and bound algorithm.
     *
     * @param path The path to start the branch and bound algorithm from.
     * @param path_length the length of the path.
     */
    void branchAndBound(Path& path, IncrementalOTP inc_solver, double path_length);

    /**
     * Runs the OTP algorithm on a path.
     *
     * @param path The path to solve the OTP algorithm on.
     * @param w_min The weight of the last edge to add to the path.
     * @return A pair containing the reward of an optimal timing profile, and the PairList needed to compute the upper bound of a path.
     */
    pair<double, PairList *> runOTP(Path& path, IncrementalOTP& inc_solver, double w_min = 0);

public:
    /**
     * Constructs a new BranchAndBound object.
     *
     * @param graph The graph to solve the AssistanceOTP for.
     * @param intervals The intervals for each vertex in the graph.
     * @param previous_max The maximum value obtained on the previous graph.
     * @param eps The epsilon approximation value.
     * @param max_interval_size The maximum interval size.
     */
    AssistanceOTP(Graph& graph, Intervals& intervals, double previous_max = 0, bool use_path_solver = false, double eps = 0,
                   double max_interval_size = 1);

    /**
     * Runs the Branch and Bound AssistanceOTP solver.
     *
     * @return The maximum reward obtained by the Branch and Bound solver.
     */
    double solve() override;
    double get_best_UB();


    void add_intervals(Intervals new_intervals);
    void update_max_value(double new_max_value);

};



#endif //DTAP_CPP_ASSISTANCE_OTP_H
