#include "incremental_otp.h"
#include <algorithm>

// Constructor: builds the internal representation of the given path and intervals,
// initializes timing structures and mappings, and prepares the incremental OTP state.
IncrementalOTP::IncrementalOTP(Path path_given, const Intervals& intervals) {
    path = path_given;
    last_reward = 0;
    min_dirty_vertex = NO_COMPUTE;
    if (path_given.empty()){
        internal_to_external.push_back(0);
        external_to_internal.resize(1);
        external_to_internal[0].push_back(0);
        ell_minus.push_back(0.0);
        ell_plus.push_back(0.0);
        add_interval(0, 0.0, 0.0);
        return;
    }

    int max_vertex_id = 0;
    for (auto &e: path)
        max_vertex_id = max({max_vertex_id, e.u, e.v});

    external_to_internal.resize(max_vertex_id + 1);

    int next_id = 0;
    internal_to_external.push_back(path.front().u);
    external_to_internal[path.front().u].push_back(next_id);
    for (auto &e: path){
        internal_to_external.push_back(e.v);
        external_to_internal[e.v].push_back(++next_id);
    }

    int N = internal_to_external.size();
    ell_minus.resize(N);
    ell_plus.resize(N);

    double cumulative_length = 0;
    for (int i = 0; i < N; ++i) {
        ell_minus[i] = cumulative_length - 0.5 * (i > 0 ? path[i - 1].l : 0);
        cumulative_length += path[i].l;
        ell_plus[i] = cumulative_length - 0.5 * path[i].l;
    }

    add_interval(0, 0.0, 0.0);
    for (auto & [idx, bucket] : intervals){
        for (const auto& interval : bucket)
            add_interval(idx, interval.start, interval.end);
    }
    for (auto & [idx, bucket] : this->intervals)
        for (auto & interval : bucket)
            interval.arrivals[0].reward = intermediate_reward(0, 0, idx, interval.arrivals[0].time);

    int index = 0;
    path.front().u = 0;
    for (int i = 0; i < path.size() - 1; i++) {
        path[i].v = ++index;
        path[i+1].u = index;
    }
    path.back().v = ++index;
}

// Computes the optimal reward incrementally by updating only affected intervals
// and returning the cached maximum reward.
double IncrementalOTP::solve() {
    if (min_dirty_vertex == 0)
        for (auto& [x,y] : intervals)
            for (auto& v : y){
                double adjusted_alpha = std::max(v.start, ell_dist_rel(x, 0));
                double reward = intermediate_reward(0, 0, x, adjusted_alpha);

                v.arrivals.clear();
                v.arrivals.emplace_back(adjusted_alpha, reward);
            }

    if (min_dirty_vertex != NO_COMPUTE){
        update_intervals();
        min_dirty_vertex = NO_COMPUTE;
        last_reward = max_reward();
    }
    return last_reward;
}

// Generates and returns all feasible (time, reward) pairs that represent optimal
// completion tradeoffs for the current intervals and path.
PairList* IncrementalOTP::get_pairs() {

    auto *pairs = new PairList();

    if (path.empty()) return pairs;

    double first_edge = path.front().l;
    intervals[path.front().u].insert(intervals[path.front().u].begin(), {first_edge / 2, first_edge / 2});
    intervals[path.back().v].emplace_back(1, 1);


    const int LAST = static_cast<int>(internal_to_external.size()) - 1;
    const double dist_to_end = ell_minus[LAST];

    for (int v = 0; v < LAST; ++v) {
        const auto vit = intervals.find(v);
        if (vit == intervals.end()) continue;

        const double dist_left = dist_to_end - ell_plus[v];

        const double exit_direct = ell_plus[v];
        double best_direct = 0.0;

        for (const auto &I : vit->second) {
            for (const auto &arr : I.arrivals) {
                if (arr.time > exit_direct + EPS) break;
                const double local = std::max(0.0, std::min(exit_direct, I.end) - std::max(I.start, arr.time));
                best_direct = std::max(best_direct, arr.reward + local);
            }
        }
        const double t_end_direct = exit_direct + dist_left;
        if (t_end_direct <= 1.0 + EPS)
            pairs->emplace_back(t_end_direct, best_direct);

        for (const auto &I : vit->second) {
            const double t_end = I.start + dist_left;
            if (t_end > 1.0 + EPS) continue;

            double best = 0.0;
            for (const auto &arr : I.arrivals) {
                if (arr.time > I.end + EPS) break;
                const double local = I.end - std::max(I.start, arr.time);
                best = std::max(best, arr.reward + local);
            }
            pairs->emplace_back(t_end, best);
        }
    }

    std::sort(pairs->begin(), pairs->end(),
              [](const TimeRewardPair& a, const TimeRewardPair& b){ return a.time < b.time - EPS; });

    intervals[path.front().u].erase(intervals[path.front().u].begin());
    intervals[path.back().v].pop_back();

    return pairs;

}

// Computes the accumulated reward obtained when traveling from a source vertex
// to a destination vertex within a given time window.
double IncrementalOTP::intermediate_reward(int src, double leave_time, int dest, double arrive_time) {
    double cumulative_reward = 0;
    double curr_time = arrive_time;
    double prev_time;

    auto it_low  = intervals.lower_bound(src);
    auto it_high = intervals.lower_bound(dest);

    using MapIt     = Intervals::const_iterator;
    auto rit_begin  = std::reverse_iterator<MapIt>(it_high);
    auto rit_end    = std::reverse_iterator<MapIt>(it_low);

    int last_idx = dest;
    for (auto rit = rit_begin; rit != rit_end; ++rit) {
        prev_time = curr_time - ell_dist_abs(last_idx, rit->first);
        curr_time = curr_time - ell_dist_abs(last_idx, rit->first + 1);
        last_idx = rit->first;
        if (last_idx <= src)
            break;
        if (prev_time < leave_time) {
            prev_time = leave_time;
            cumulative_reward += reward(last_idx, prev_time, curr_time);
            break;
        }
        cumulative_reward += reward(last_idx, prev_time, curr_time);
        curr_time = prev_time;
    }
    return cumulative_reward;
}

// Computes the reward collected at a single vertex between entry and exit times
// by summing overlaps with its intervals.
double IncrementalOTP::reward(int i, double t_entry, double t_exit) const {
    auto curr_I = intervals.find(i);
    if (curr_I == intervals.end()) return 0;

    const auto& vil = curr_I->second;
    double r = 0;
    for (const auto& iv : vil)
        r += std::max(0.0, std::min(t_exit, iv.end) - std::max(t_entry, iv.start));

    return r;
}

// Removes a specific interval from all internal copies of a given vertex and
// marks affected vertices for recomputation.
void IncrementalOTP::remove_interval(int vertex_idx, double alpha, double beta) {
    for (int internal : external_to_internal[vertex_idx]) {
        auto &vec = intervals[internal];
        vec.erase(
                std::remove_if(vec.begin(), vec.end(),
                               [&](const Interval &interval) {
                                   return std::abs(interval.start - alpha) < EPS &&
                                          std::abs(interval.end - beta) < EPS;
                               }),
                vec.end()
        );

        if (vec.empty()) {
            intervals.erase(internal);
        }

        min_dirty_vertex = std::min(min_dirty_vertex, internal);
    }
}

// Adds a new interval to all internal copies of a vertex, initializes its arrival
// states, and marks affected vertices as dirty.
void IncrementalOTP::add_interval(int vertex_idx, double alpha, double beta) {
    for (int internal : external_to_internal[vertex_idx]) {
        if (ell_dist_rel(internal, 0) > 1.0) return;

        double adjusted_alpha = std::max(alpha, ell_dist_rel(internal, 0));
        double reward = intermediate_reward(0, 0, internal, adjusted_alpha);

        ArrivalsList arrivals;
        arrivals.reserve(20);
        arrivals.emplace_back(adjusted_alpha, reward);

        auto &vec = intervals[internal];

        vec.erase(
                std::remove_if(vec.begin(), vec.end(),
                               [&](const Interval &interval) {
                                   return std::abs(interval.start - alpha) < EPS &&
                                          std::abs(interval.end   - beta)  < EPS;
                               }),
                vec.end()
        );

        vec.emplace_back(alpha, beta, arrivals);

        min_dirty_vertex = std::min(min_dirty_vertex, internal);
    }
}
// Propagates reward from a source interval to a destination interval by computing
// feasible arrival times and updating dominating arrival states.
void IncrementalOTP::add_reward(double start_time, double arrive_time, double reward_before,
                                const Interval &src_interval, Interval &dest_interval, int src, int dest) {
    double travel           = ell_dist(src, dest, src_interval);
    double dest_wait_arrive = (src == dest) ? 0.0 : 0.5 * path[dest - 1].l;
    double dest_wait_leave  = (src == dest) ? 0.0 : 0.5 * path[dest].l;

    double leave_time   = arrive_time - travel;
    double end_time     = arrive_time + dest_wait_arrive + dest_wait_leave;

    if (arrive_time - EPS < ell_dist_rel(dest, 0) ||
        (leave_time - EPS < start_time && src != dest))
        return;

    double reward_src   = max(0.0, min(leave_time, src_interval.end) - max(src_interval.start, start_time));
    double reward_trans = intermediate_reward(src, leave_time, dest, arrive_time);
    double reward_dest  = max(0.0, min(end_time, dest_interval.end) - max(dest_interval.start, arrive_time));

    add_dominating(dest_interval.arrivals, end_time, reward_before + reward_src + reward_trans + reward_dest);
}

// Inserts a new arrival state into a list while removing dominated states and
// maintaining time sorted order.
void IncrementalOTP::add_dominating(ArrivalsList &arrivals, double t, double r) {
    if (arrivals.empty()) {
        arrivals.emplace_back(t, r);
        return;
    }

    auto it = arrivals.begin();
    for (; it != arrivals.end(); ++it) {
        if (!(it->time < t)) break;
    }

    for (size_t i = 1; i < arrivals.size(); ++i) {
        if (arrivals[i].time < arrivals[i-1].time) {
            std::sort(arrivals.begin(), arrivals.end(),
                      [](auto &a, auto &b){ return a.time < b.time; });
            break;
        }
    }

    bool dominated = false;
    for (auto &p : arrivals) {
        if (p.reward >= r && p.time <= t) {
            dominated = true;
        }
    }
    for (auto &p : arrivals) {
        if (p.reward >= r && p.time <= t) {
            dominated = true;
        }
    }
    if (dominated) return;

    auto jt = it;
    while (jt != arrivals.end()) {
        if (jt->reward <= r && jt->time >= t - EPS) {
            ++jt;
        } else {
            ++jt;
        }
    }

    std::sort(arrivals.begin(), arrivals.end(),
              [](auto &a, auto &b){ return a.time < b.time; });

    arrivals.emplace(it, t, r);

    for (size_t i = 1; i < arrivals.size(); ++i) {
        if (arrivals[i].time < arrivals[i-1].time) {
            std::stable_sort(arrivals.begin(), arrivals.end(),
                             [](auto &a, auto &b){ return a.time < b.time; });
            break;
        }
    }
}

// Returns the maximum achievable reward over all intervals and arrival states.
double IncrementalOTP::max_reward() {
    double maxReward = 0.0;

    for (auto& [idx, vil] : intervals)
        for (const auto& interval : vil)
            for (const auto& arrival : interval.arrivals)
                maxReward = std::max(maxReward, arrival.reward + interval.end - std::clamp(arrival.time, interval.start, interval.end));

    return maxReward;
}

// Recomputes interval arrival states starting from the earliest dirty vertex
// using dynamic programming style propagation.
void IncrementalOTP::update_intervals() {
    if (min_dirty_vertex == NO_COMPUTE) return;

    Intervals prev_intervals = intervals, curr_intervals;
    auto limit_dirty = prev_intervals.lower_bound(min_dirty_vertex);
    for (auto it = prev_intervals.begin(); it != limit_dirty; ++it)
        curr_intervals[it->first] = it->second;

    for (auto dest_it = limit_dirty; dest_it != prev_intervals.end(); ++dest_it) {
        int dest = dest_it->first;
        auto limit = curr_intervals.upper_bound(dest);
        for (auto& dest_interval : dest_it->second) {
            auto& dest_list = curr_intervals[dest];
            for (auto src_it = curr_intervals.begin(); src_it != limit; ++src_it){
                int src = src_it->first;
                for (const auto& src_interval : src_it->second) {
                    double travel = ell_dist(src, dest, src_interval);
                    for (const auto& arrival : src_interval.arrivals) {
                        double t0 = arrival.time, r0 = arrival.reward;
                        add_reward(t0, dest_interval.start, r0, src_interval, dest_interval, src, dest);
                        add_reward(t0, src_interval.end + travel, r0, src_interval, dest_interval, src, dest);
                    }
                }
            }
            dest_list.push_back(dest_interval);
        }
    }

    intervals = curr_intervals;
}

// Appends a new vertex and edge to the path and updates internal timing structures.
void IncrementalOTP::add_vertex(int v, double l) {
    int last_internal = (path.empty()) ? 0 : path.back().v;
    int new_internal  = static_cast<int>(internal_to_external.size());

    internal_to_external.push_back(v);
    if (v >= static_cast<int>(external_to_internal.size()))
        external_to_internal.resize(v + 1);
    external_to_internal[v].push_back(new_internal);

    path.emplace_back(last_internal, new_internal, l);

    ell_plus[ell_plus.size()-1] = ell_plus.back() + 0.5 * l;
    double cumulative_length = ell_plus.back();
    ell_minus.push_back(cumulative_length);
    ell_plus.push_back(cumulative_length + 0.5 * l);
}