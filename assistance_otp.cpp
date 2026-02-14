#include "assistance_otp.h"


void dijkstra(Graph& graph, int start, vector<double>& dist) {
    set<pair<double, int>> s;
    s.insert({0, start});
    dist[start] = 0;
    while (!s.empty()) {
        auto it = s.begin();
        int u = it->second;
        s.erase(it);
        for (auto& adj: graph[u]) {
            int v = adj.v;
            double w = adj.l;
            if (dist[v] > dist[u] + w) {
                s.erase({dist[v], v});
                dist[v] = dist[u] + w;
                s.insert({dist[v], v});
            }
        }
    }
}


GraphBound::GraphBound(Graph& graph, Intervals& intervals) : graph(graph), intervals(intervals) {
    min_dist.resize(graph.size(), vector<double>(graph.size(), INF));
    min_dist_no_last.resize(graph.size(), vector<double>(graph.size(), INF));
    computeDists();
    computeUBI();
}

void GraphBound::computeDists() {
    // real distances
    Matrix distances = Matrix(graph.size(), vector<double>(graph.size(), INF));
    for (int u = 0; u < graph.size(); u++) {
        dijkstra(graph, u, distances[u]);
    }

    // distances without half of last edge
    min_dist_no_last = distances;
    for (int s = 0; s < graph.size(); s++) {
        for (int u = 0; u < graph.size(); u++) {
            for (auto& edge: graph[u]) {
                int v = edge.v;
                double w = edge.l;
                min_dist_no_last[s][v] = min(min_dist_no_last[s][v], distances[s][u] + w / 2);
            }
        }
    }

    // distances without half of first and last edge
    min_dist = min_dist_no_last;
    for (int u = 0; u < graph.size(); u++) {
        for (auto& edge: graph[u]) {
            for (int t = 0; t < graph.size(); t++) {
                int v = edge.v;
                double w = edge.l;
                min_dist[u][t] = min(min_dist[u][t], min_dist_no_last[v][t] + w / 2);
            }
        }
    }

    // set distances to 0 for all edges
    for (int u = 0; u < graph.size(); u++) {
        for (auto edge: graph[u])
            min_dist[u][edge.v] = 0;
    }
}

void GraphBound::computeUBI() {
    // create list of (interval, vertex) pairs
    vector<VInterval> v_intervals;
    for (int u = 0; u < intervals.size(); u++) {
        for (auto& interval: intervals[u]) {
            v_intervals.emplace_back(u, interval);
        }
    }
    // sort by interval.end from largest to smallest
    sort(v_intervals.begin(), v_intervals.end(),
         [](const VInterval& a, const VInterval& b) { return a.end > b.end; }
    );

    // initialize UBi
    for (auto v_interval: v_intervals) {
        UBi[v_interval] = v_interval.end - v_interval.start;
    }

    // iterate over all intervals
    bool changed = true; // flag to check if we need to iterate again
    while (changed) {
        changed = false;

        for (auto v_interval: v_intervals) {
            for (auto next_interval: v_intervals) {
                // consider only intervals that are after v_interval
                if (v_interval.start >= next_interval.end - EPS) {
                    continue;
                }

                auto dist = min_dist[v_interval.vertex][next_interval.vertex];
                auto exit_time = min(max(v_interval.start, next_interval.start - dist), v_interval.end);
                auto arrival_time = exit_time + dist;

                auto reward_in_interval = exit_time - v_interval.start;
                auto future_reward = UBi[next_interval] - max(0.0, arrival_time - next_interval.start);
                auto reward = reward_in_interval + future_reward;
                if (reward > UBi[v_interval] + EPS) {
                    UBi[v_interval] = reward;
                    changed = true;
                }
            }
        }
    }
}

double GraphBound::UB(int u, double time) {
    double bound = 0;
    // iterate over v_intervals
    for (auto pair: UBi) {
        auto v_interval = pair.first;
        auto reward = pair.second;

        auto distance = min_dist_no_last[u][v_interval.vertex];
        if (time + distance < v_interval.end + EPS) {
            auto overlap = max(time + distance - v_interval.start, 0.0);
            bound = max(bound, reward - overlap);
        }
    }

    return bound;
}


GraphBoundWrapper::GraphBoundWrapper(Graph& graph, Intervals& intervals) : bound(graph, intervals), max_inner_UB(0.0) {
    // initialize cache
    for (int u = 0; u < graph.size(); u++) {
        cache[u] = PairSet();
    }
}

void GraphBoundWrapper::addToCache(int u, double time, double reward) {
    cache[u].insert({time, reward});
}

pair<double, bool> GraphBoundWrapper::getFromCache(int u, double time) {
    // if cache is empty, return false
    if (cache[u].empty())
        return {1, false};

    // find the closest time that is less than or equal to time
    auto it = cache[u].lower_bound({time, 0});
    if (it == cache[u].begin())
        return {1, false};
    it--;
    return {it->reward, it->time == time};
}


bool GraphBoundWrapper::shouldPrune(int u, PairList& bound_pairs, double max_reward) {
    max_inner_UB = 0;
    for (auto& pair: bound_pairs)
        max_inner_UB = max(max_inner_UB, pair.reward + bound.UB(u, pair.time));

    for (auto& pair: bound_pairs) {
        auto exit_time = pair.time;
        auto reward = pair.reward;

        auto entry = getFromCache(u, exit_time);
        if (reward + entry.first <= max_reward)
            continue;
        // if entry is accurate, we cant prune
        if (entry.second) {
            return false;
        }

        // compute bound
        auto real_bound = bound.UB(u, exit_time);
        addToCache(u, exit_time, real_bound);
        if (reward + real_bound > max_reward + EPS) {
            return false;
        }
    }

    return true;
}

double GraphBoundWrapper::get_max_inner_UB() {
    return max_inner_UB;
}


AssistanceOTP::AssistanceOTP(Graph& graph, Intervals& intervals, double previous_max, bool use_path_solver, double eps,
                               double max_interval_size) : AbstractAssistanceOTP(graph, intervals),
                                                           value_history(graph.size()), eps(eps),
                                                           bound(graph, intervals),
                                                           max_interval_size(max_interval_size),
                                                           dummy_vertex((int) graph.size()), path_solver() {
    max_value = previous_max;
    preprocess_intervals();
    preprocess_graph();
    intervals[dummy_vertex] = {};
}


void AssistanceOTP::preprocess_intervals() {
    // if interval length is more than max_interval_size, split it
    for (auto& [v, vec]: intervals) {
        vector<Interval> temp_list;
        for (auto& interval: vec) {
            auto start = interval.start, end = interval.end;
            while (end - start > max_interval_size) {
                temp_list.push_back({Interval(start, start + max_interval_size)});
                start += max_interval_size;
            }
            if (end - start > EPS)
                temp_list.push_back({Interval(start, end)});
        }
        vec = temp_list;
    }
    if (min_time_vertex.empty())
        return;

    for (auto& [v, vec] : intervals) {
        double min_time = min_time_vertex[v];
        std::vector<Interval> new_vec;

        for (const auto& interval : vec) {
            if (interval.end <= min_time) {
                continue;
            }
            if (interval.start < min_time) {
                new_vec.emplace_back(min_time, interval.end);
            } else {
                new_vec.push_back(interval);
            }
        }
        vec = std::move(new_vec);
    }

}

void AssistanceOTP::preprocess_graph() {
    // sort adj list from longest edge to shortest
    for (auto& vec: graph) {
        sort(vec.begin(), vec.end(), [this](const AdjEntry& a, const AdjEntry& b) {
            if (intervals[a.v].size() == intervals[b.v].size())
                return a.l > b.l;
            return intervals[a.v].size() > intervals[b.v].size();
        });
    }
}

pair<double, PairList *> AssistanceOTP::runOTP(Path& path, IncrementalOTP& inc_solver, double w_min) {
    path_count++;

    int last_vertex = path.empty() ? 0 : path.back().v;
    double last_weight = path.empty() ? 0.0 : path.back().l;

    if (!path.empty())
        inc_solver.add_vertex(last_vertex, last_weight);
    for (const auto& interval : intervals[last_vertex])
        inc_solver.add_interval(last_vertex, interval.start, interval.end);
    auto rwd = inc_solver.solve();

    path.emplace_back(last_vertex, dummy_vertex, w_min);
    auto pl = inc_solver.get_pairs();
    path.pop_back();

    return {rwd, pl};
}


double AssistanceOTP::solve() {
    path_count = 0;
    start_time = std::chrono::high_resolution_clock::now();
    inner_UB = 0;
    auto path = Path();
    IncrementalOTP inc_solver({}, {});
    branchAndBound(path, inc_solver, 0);
    added_intervals = {};
    return max_value;
}

void AssistanceOTP::branchAndBound(Path& path, IncrementalOTP inc_solver, double path_length) {
    // check if time limit is reached
    auto current_time = std::chrono::high_resolution_clock::now();
    auto duration =
            (double) std::chrono::duration_cast<std::chrono::milliseconds>(current_time - start_time).count() / 1000.;
    if (duration > time_limit) {
        cout << "AssistanceOTP: Time limit reached" << endl;
        return;
    }

    int last_vertex = path.empty() ? -1 : path.back().u;
    int current_vertex = path.empty() ? 0 : path.back().v;

    // min edge length from current vertex (this can be added to the path when running OTP)
    double w_min = graph[current_vertex].empty() ? 0 : graph[current_vertex].back().l;

    // solve OTP
    auto pair = runOTP(path, inc_solver, w_min);
    double value = pair.first;
    PairList *bound_pairs = pair.second;

    // prune unnecessary loops using ValueHistory to keep track of
    if (value_history.last(current_vertex) > value - EPS) {
        return;
    }

    // update real value
    value_history.add(current_vertex, value);

    // update max reward
    if (value > max_value) {
        max_value = value;
        max_path = path;
    }

    // iterate over all adjacent edges
    for (auto& adj: graph[current_vertex]) {
        // if vertex is not reachable
        if (path_length + adj.l > 1)
            continue;

        // prune unnecessary loops
        if (adj.v == last_vertex && intervals[adj.v].empty())
            continue;

        // prune if upper bound is less than max value
        bool prune = bound.shouldPrune(adj.v, *bound_pairs, eps + max_value); // NOTICE: was (1 + eps) * max_value
        inner_UB = max(inner_UB, bound.get_max_inner_UB());
        if (prune)
            continue;

        // add edge and solve recursive call
        path.emplace_back(current_vertex, adj.v, adj.l);
        branchAndBound(path, inc_solver, path_length + adj.l);
        path.pop_back();
    }

    // update value history
    value_history.removeLast(current_vertex);

    delete bound_pairs;
}

void AssistanceOTP::add_intervals(Intervals new_intervals) {
    if (!min_time_vertex.empty()) {
        for (auto &[v, vec]: new_intervals) {
            double min_time = min_time_vertex[v];
            std::vector<Interval> new_vec;

            for (const auto &interval: vec) {
                if (interval.end <= min_time) {
                    continue; // Discard interval
                }
                if (interval.start < min_time) {
                    new_vec.push_back({min_time, interval.end}); // Clamp start
                } else {
                    new_vec.push_back(interval); // Keep as-is
                }
            }

            vec = std::move(new_vec); // Replace with filtered/updated intervals
        }
    }

    added_intervals = new_intervals;

    for (auto& [v, new_list] : new_intervals) {
        auto& current_list = intervals[v];

        current_list.insert(current_list.end(), new_list.begin(), new_list.end());

        VertexIntervalList merged;
        for (const Interval& interval : current_list) {
            if (!merged.empty() && std::abs(merged.back().end - interval.start) <= EPS) {
                merged.back().end = interval.end;
            } else {
                merged.push_back(interval);
            }
        }
        current_list = std::move(merged);
    }


    bound = GraphBoundWrapper(graph, intervals);
    preprocess_graph();
    preprocess_intervals();
}

void AssistanceOTP::update_max_value(double new_max_value) {
    max_value = new_max_value;
}
