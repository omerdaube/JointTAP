#ifndef DTAP_CPP_LP_BOUND_H
#define DTAP_CPP_LP_BOUND_H

#include "utils.h"
#include <iostream>
#include <random>
#include <unordered_map>
#include <queue>
#include <vector>
#include <utility>
#include <cmath>
#include "utils.h"

const int sID = -1;
const int tID = -2;

inline constexpr long long makeKey(int u, int v)
{
    return ( (static_cast<long long>(u) << 32) ^
             (static_cast<long long>(v) & 0xFFFFFFFFLL) );
}

// Computes all pairs shortest path distances on graphA using Floyd Warshall
// and returns the distance matrix as nested hash maps.
std::unordered_map<int,std::unordered_map<int,double>> buildAndRunFloydWarshall(const Graph& graphA);

struct EdgeInfo {
    int    to;
    double lenT;
    double lenA;
    double reward;
};

// Builds the raw product graph of the task and assistant graphs and returns
// adjacency information together with mappings between product nodes and pairs.
std::pair< std::unordered_map<int,std::vector<int>>,
        std::pair< std::unordered_map<int,std::pair<int,int>>,
                std::unordered_map<long long,int>> >
                build_product_graph(const Graph& graphT,const Graph& graphA, bool newIdea = false);

// Constructs an annotated product graph by attaching traversal lengths and
// rewards to edges based on assistance feasibility and distances.
unordered_map<int, vector<EdgeInfo>>
buildAnnotatedProd(const Graph& graphT,
                   const Graph& graphA,
                   const unordered_map<int, vector<int>>& prodRaw,
                   const unordered_map<int, pair<int,int>>& idToPair,
                   const PairsSet& assistance,
                   const unordered_map<int, unordered_map<int,double>>& deltaEllA,
                   bool newIdea);

// Reverses all edges in an annotated product graph.
unordered_map<int, vector<EdgeInfo>> reverseAnnotatedProd(const unordered_map<int, vector<EdgeInfo>>& adj);

// Computes a topological ordering of a directed acyclic annotated graph.
vector<int> topoSort(const unordered_map<int, vector<EdgeInfo>>& adj);

// class to calculate the UB
class UBOracle {
public:
    // UBOracle constructor: preprocesses graphs and builds internal data structures
    // required for efficient upper bound computation.
    UBOracle(const Graph& graphT,
             const Graph& graphA,
             const PairsSet& assistance)
            : T_(graphT)
            , A_(graphA)
    {
        distA0_ = dijkstra(A_, /*src=*/0);
        okA_.reserve(distA0_.size());
        for (double d : distA0_) okA_.push_back(d <= 1.0 + 1e-9);

        auto deltaEllA = buildAndRunFloydWarshall(A_);
        auto [prodRaw, maps] = build_product_graph(T_, A_);
        idToPair_ = maps.first;
        pairToID_ = maps.second;

        fullRev_   = reverseAnnotatedProd(
                buildAnnotatedProd(T_, A_, prodRaw,
                                   idToPair_, assistance, deltaEllA, false));
        revTopo_   = topoSort(fullRev_);

        distT_.resize(T_.size());
    }

    // Computes an admissible upper bound on achievable reward from a task vertex
    // within the given remaining time budget.
    double getUB(int vTsrc, double time)
    {
        const std::vector<double>& dT = ensureDistT(vTsrc);

        std::unordered_set<int> keep = {sID, tID};
        keep.reserve(idToPair_.size() / 4);

        for (const auto& [id, uv] : idToPair_)
            if (okA_[uv.second] && dT[uv.first] <= time + 1e-9)
                keep.insert(id);

        std::unordered_map<int,double> front;
        front[tID] = 0.0;
        double best = -1.0;

        for (int u : revTopo_) {
            if (!keep.count(u)) continue;
            auto it = front.find(u);
            if (it == front.end()) continue;

            double acc = it->second;

            auto itPair = idToPair_.find(u);
            if (itPair != idToPair_.end() && itPair->second.first == vTsrc)
                best = std::max(best, acc);

            for (const EdgeInfo& e : fullRev_[u])
                if (keep.count(e.to))
                    front[e.to] = std::max(front[e.to], acc + e.reward);
        }
        return min(best, time);
    }

private:
    // Runs Dijkstra shortest path from a source vertex on a graph and returns distances.
    static std::vector<double> dijkstra(const Graph& G, int src)
    {
        const double INF = 1e15;
        using Q = std::pair<double,int>;
        std::vector<double> dist(G.size(), INF);
        std::priority_queue<Q, std::vector<Q>, std::greater<Q>> pq;
        dist[src]=0.0; pq.emplace(0.0,src);
        while(!pq.empty()){
            auto[du,u]=pq.top(); pq.pop();
            if(du!=dist[u]) continue;
            for(const AdjEntry&e:G[u])
                if(dist[e.v]>du+e.l){
                    dist[e.v]=du+e.l; pq.emplace(dist[e.v],e.v);
                }
        }
        return dist;
    }

    // Ensures shortest path distances from a task vertex are computed and cached.
    const std::vector<double>& ensureDistT(int vT)
    {
        if (distT_[vT].empty())
            distT_[vT] = dijkstra(T_, vT);
        return distT_[vT];
    }

    const Graph& T_;
    const Graph& A_;

    std::vector<double> distA0_;
    std::vector<char>   okA_;

    std::unordered_map<int,std::pair<int,int>> idToPair_;
    std::unordered_map<long long,int>          pairToID_;
    std::unordered_map<int,std::vector<EdgeInfo>> fullRev_;
    std::vector<int>    revTopo_;

    std::vector<std::vector<double>> distT_;
};


#endif //DTAP_CPP_LP_BOUND_H
