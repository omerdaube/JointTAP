#include "joint_tap.h"

// Constructor: initializes the JointTAP instance with the assistant and task graphs,
// assistance pairs, and upper bound structure, then preprocesses the task graph.
JointTAP::JointTAP(Graph graphA, Graph graphT, PairsSet assistance)
        : graphA(graphA), graphT(graphT), assistance(assistance), max_reward(0),
          value_history(graphT.size()), current_intervals(),
          UBo(this->graphT, this->graphA, this->assistance) {
    preprocess_graph();
}

// Checks whether assistant vertex vA is allowed to assist task vertex vT
// according to the predefined assistance pairs.
bool JointTAP::can_assist(int vT, int vA) {
    return assistance.count({vA, vT});
}

// Solves the JointTAP problem using branch and bound by initializing the search
// and returning the maximum reward found.
double JointTAP::solve() {
    max_reward = 0;
    current_intervals.clear();
    auto pathT = Path();
    Intervals empty_intervals;
    AssistanceOTP ll_solver(graphA, empty_intervals);  // persistent instance
    BnB(pathT, 0, ll_solver);
    return max_reward;
}

// Performs a recursive branch and bound search over task paths, uses the low level
// solver to evaluate rewards, prunes suboptimal branches, and updates the best solution.
void JointTAP::BnB(Path& pathT, double path_length, AssistanceOTP& ll_solver) {
    int current_vertex = pathT.empty() ? 0 : pathT.back().v;

    ll_solver.update_max_value(0);
    auto value = ll_solver.solve();
    if (value + UBo.getUB(current_vertex, 1 - path_length) <= max_reward) {
        return;
    }

    if (value_history.last(current_vertex) > value - EPS)
        return;
    value_history.add(current_vertex, value);

    if (value > max_reward) {
        max_reward = value;
        best_pathA = ll_solver.get_best_path();
        best_pathT = pathT;
    }

    for (auto& adj: graphT[current_vertex]) {
        if (path_length + adj.l > 1.0)
            continue;

        pathT.emplace_back(current_vertex, adj.v, adj.l);

        AssistanceOTP next_solver = ll_solver;
        Intervals new_intervals = add_current_interval(pathT, path_length + adj.l);
        next_solver.add_intervals(new_intervals);
        BnB(pathT, path_length + adj.l, next_solver);
        pathT.pop_back();
    }

    value_history.removeLast(current_vertex);
}

// Computes all assistance time intervals induced by a given task path by splitting
// each edge into two halves and assigning intervals to compatible assistant vertices.
Intervals JointTAP::calculate_intervals(Path pathT) {
    Intervals result;
    if (pathT.empty()) return result;

    double acc_len = 0.0;

    for (const auto& e : pathT) {
        const int  u   = e.u;
        const int  v   = e.v;
        const double L = e.l;

        const double mid = acc_len + 0.5 * L;

        double s_u = acc_len;
        double e_u = mid;
        if (s_u < 1.0) {
            if (e_u > 1.0) e_u = 1.0;
            for (int a = 0; a < static_cast<int>(graphA.size()); ++a)
                if (can_assist(u, a))
                    result[a].emplace_back(s_u, e_u);
        }

        double s_v = mid;
        double e_v = acc_len + L;
        if (s_v < 1.0) {
            if (e_v > 1.0) e_v = 1.0;
            for (int a = 0; a < static_cast<int>(graphA.size()); ++a)
                if (can_assist(v, a))
                    result[a].emplace_back(s_v, e_v);
        }

        acc_len += L;
        if (acc_len >= 1.0) break;
    }

    return result;
}

// Computes and returns only the new assistance intervals introduced by the most
// recently added edge in the current task path.
Intervals JointTAP::add_current_interval(Path &pathT, double path_length) {
    Intervals new_intervals;

    int n = pathT.size();
    if (n == 0) return {};

    const auto& last_edge = pathT.back();
    int u = last_edge.u;
    int v = last_edge.v;
    double len = last_edge.l;

    double mid = path_length - len + len * 0.5;

    double start_u = path_length - len;
    double end_u   = mid;
    if (start_u < 1.0) {
        if (end_u > 1.0) end_u = 1.0;
        for (int j = 0; j < graphA.size(); ++j)
            if (can_assist(u, j))
                if (end_u - start_u > EPS)
                    new_intervals[j].emplace_back(start_u, end_u);
    }

    double start_v = mid;
    double end_v   = path_length;
    if (start_v < 1.0) {
        if (end_v > 1.0) end_v = 1.0;
        for (int j = 0; j < graphA.size(); ++j)
            if (can_assist(v, j))
                if (end_v - start_v > EPS)
                    new_intervals[j].emplace_back(start_v, end_v);
    }

    return new_intervals;
}

// Solves the JointTAP problem using a full depth first search without pruning
// and returns the best reward found.
double JointTAP::solve_dfs() {
    max_reward_dfs = 0;
    auto pathT = Path();
    DFS(pathT, 0);
    return max_reward_dfs;
}

// Recursively explores all feasible task paths with DFS, evaluates each path
// with the low level solver, and tracks the best reward and path.
void JointTAP::DFS(Path& pathT, double path_length) {
    int current_vertex = pathT.empty() ? 0 : pathT.back().v;

    auto intervals = calculate_intervals(pathT);
    AssistanceOTP ll_solver(graphA, intervals);  // persistent instance
    auto value = ll_solver.solve();
    if (value > max_reward_dfs){
        max_reward_dfs = value;
        max_path_dfs = pathT;
    }

    for (auto& adj: graphT[current_vertex]) {
        if (path_length + adj.l > 1.0)
            continue;

        pathT.emplace_back(current_vertex, adj.v, adj.l);
        DFS(pathT, path_length + adj.l);
        pathT.pop_back();
    }
}

// Preprocesses the task graph by sorting outgoing edges of each vertex according
// to how frequently their destinations appear in the assistance pairs.
void JointTAP::preprocess_graph() {
    std::unordered_map<int, int> freq;
    for (const auto& p : assistance)
        ++freq[p.second];

    for (auto& vec: graphT) {
        sort(vec.begin(), vec.end(), [&freq](const AdjEntry& a, const AdjEntry& b) {
            int fa = freq.count(a.v) ? freq.at(a.v) : 0;
            int fb = freq.count(b.v) ? freq.at(b.v) : 0;
            if (fa == fb) return a.v < b.v;
            return fa > fb;
        });
    }
}

// Returns the best assistant path and task path found during the search.
tuple<Path, Path> JointTAP::get_paths() {
    return {best_pathA, best_pathT};
}