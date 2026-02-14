#ifndef DTAP_CPP_INCREMENTAL_OTP_H
#define DTAP_CPP_INCREMENTAL_OTP_H


#include "utils.h"

class IncrementalOTP {
private:
    Path path;
    Intervals intervals;
    vector<double> ell_minus;
    vector<double> ell_plus;
    std::vector<int> internal_to_external;
    std::vector<std::vector<int>> external_to_internal;
    double last_reward;

    inline double ell_dist_abs(int i, int j) {
        return ell_minus[i] - ell_minus[j];
    }

    inline double ell_dist_rel(int i, int j) {
        return (i == j) ? 0 : ell_minus[i] - ell_plus[j];
    }

    inline double ell_dist(int src, int dest, const Interval& src_interval) {
        return ell_dist_rel(dest, src);
    }

    void add_reward(double start_time, double arrive_time, double reward_before,
                    const Interval& src_interval, Interval& dest_interval,
                    int src, int dest);

    double reward(int i, double t_entry, double t_exit) const;

    double intermediate_reward(int src, double leave_time,
                               int dest, double arrive_time);

    void add_dominating(ArrivalsList& arrivals, double t, double r);

    double max_reward();

    void update_intervals();

public:
    int min_dirty_vertex;

    IncrementalOTP(Path path, const Intervals& intervals);
    void add_interval(int vertex_idx, double alpha, double beta);
    void remove_interval(int vertex_idx, double alpha, double beta);
    double solve();
    PairList* get_pairs();
    void add_vertex(int v, double l);
};

typedef std::unordered_map<Path, IncrementalOTP, PathHash, PathEqual> PathSolver;

#endif //DTAP_CPP_INCREMENTAL_OTP_H
