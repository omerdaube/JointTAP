#ifndef DTAP_CPP_JOINT_TAP_H
#define DTAP_CPP_JOINT_TAP_H

#include "utils.h"
#include "lp_bound.h"
#include "assistance_otp.h"
#include <queue>

class JointTAP {
private:
    Graph graphA;
    Graph graphT;
    PairsSet assistance;
    double max_reward;
    double max_reward_dfs;
    Path max_path_dfs;
    ValueHistory value_history;
    Intervals current_intervals;
    Path best_pathA;
    Path best_pathT;
    UBOracle UBo;

    bool can_assist(int vT, int vA);
    void BnB(Path& pathT, double path_length, AssistanceOTP& ll_solver);
    void DFS(Path& pathT, double path_length);
    Intervals add_current_interval(Path &pathT, double path_length);
    void preprocess_graph();

public:
    JointTAP(Graph graphA, Graph graphT, PairsSet assistance);
    double solve();
    double solve_dfs();
    tuple<Path,Path> get_paths();
    Intervals calculate_intervals(Path pathT);
};

#endif //DTAP_CPP_JOINT_TAP_H
