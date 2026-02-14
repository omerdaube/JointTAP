#include "lp_bound.h"

using namespace std;

std::unordered_map<int,std::unordered_map<int,double>> buildAndRunFloydWarshall(const Graph& graphA)
{
    const double INF = 1e15;
    int n = static_cast<int>(graphA.size());

    std::vector<std::pair<int,int>> edgeList;
    for (int u=0; u<n; ++u)
        for (const AdjEntry& e : graphA[u])
            edgeList.emplace_back(u,e.v);

    int m = static_cast<int>(edgeList.size());
    int totalNodes = n + m;

    std::vector<std::vector<double>>
            dist(totalNodes, std::vector<double>(totalNodes, INF));

    for (int i=0;i<totalNodes;++i) dist[i][i]=0.0;

    for (int idx=0; idx<m; ++idx)
    {
        int u=edgeList[idx].first, v=edgeList[idx].second;
        double w = 0.0;
        for (const AdjEntry& e : graphA[u])
            if (e.v==v){ w=e.l; break; }

        int eNode = n+idx;
        dist[u][eNode] = dist[eNode][v] = w*0.5;
    }

    for (int k=0;k<totalNodes;++k)
        for (int i=0;i<totalNodes;++i) if (dist[i][k]<INF)
                for (int j=0;j<totalNodes;++j) if (dist[k][j]<INF)
                        if (dist[i][j] > dist[i][k]+dist[k][j])
                            dist[i][j] = dist[i][k]+dist[k][j];

    std::unordered_map<int,std::unordered_map<int,double>> delta;
    for (int u=0;u<n;++u)
        for (int v=0;v<n;++v)
        {
            double best=INF;
            for (int idx=0; idx<m; ++idx)
                if (edgeList[idx].second==v)
                    best = std::min(best, dist[u][n+idx]);
            if (best<INF/2) delta[u][v]=best;
        }
    return delta;
}

std::vector<int> bfs_reachable(const Graph& graphA,int start)
{
    int n = static_cast<int>(graphA.size());
    std::vector<int> res;
    std::vector<char> vis(n,0);
    std::queue<int> Q;  Q.push(start);  vis[start]=1;

    while (!Q.empty()){
        int u=Q.front(); Q.pop(); res.push_back(u);
        for (const AdjEntry& e : graphA[u])
            if (!vis[e.v]){ vis[e.v]=1; Q.push(e.v);}
    }
    return res;
}

std::pair< std::unordered_map<int,std::vector<int>>,
        std::pair< std::unordered_map<int,std::pair<int,int>>,
                std::unordered_map<long long,int>> >
                build_product_graph(const Graph& graphT,const Graph& graphA, bool newIdea)
{
    int nT = static_cast<int>(graphT.size());
    int nA = static_cast<int>(graphA.size());

    std::vector<std::vector<int>> reachA(nA);
    for (int v=0;v<nA;++v) reachA[v]=bfs_reachable(graphA,v);

    std::vector<std::vector<int>> reachT(nT);
    for (int v=0;v<nT;++v) reachT[v]=bfs_reachable(graphT,v);

    std::unordered_map<long long,int> pairToID;
    std::unordered_map<int,std::pair<int,int>> idToPair;
    std::unordered_map<int,std::vector<int>>  adjProd;
    adjProd[sID]={tID}; adjProd[tID]={};

    int nextID=0;
    for (int u=0;u<nT;++u)
        for (int v=0;v<nA;++v){
            long long key = makeKey(u,v);
            pairToID[key]=nextID;
            idToPair[nextID]={u,v};
            adjProd[nextID]={};
            ++nextID;
        }

    for (auto& kv : idToPair)
        adjProd[kv.first].push_back(tID);

    if (!newIdea){
        for (int u=0;u<nT;++u)
            for (const AdjEntry& eT : graphT[u])
                for (int v=0;v<nA;++v){
                    int fromID = pairToID[ makeKey(u ,v) ];
                    for (int w : reachA[v]){
                        int toID = pairToID[ makeKey(eT.v , w) ];
                        adjProd[fromID].push_back(toID);
                    }
                }
    } else {
        for (int u = 0; u < nT; ++u) {
            for (int v = 0; v < nA; ++v) {
                int fromID = pairToID[makeKey(u, v)];
                for (int u2 : reachT[u]) {
                    for (int v2 : reachA[v]) {
                        int toID = pairToID[makeKey(u2, v2)];
                        adjProd[fromID].push_back(toID);
                    }
                }
            }
        }
    }

    return {adjProd,{idToPair,pairToID}};
}

unordered_map<int, vector<EdgeInfo>>
buildAnnotatedProd(const Graph&                       graphT,
                   const Graph&                       graphA,
                   const unordered_map<int, vector<int>>& prodRaw,
                   const unordered_map<int, pair<int,int>>& idToPair,
                   const PairsSet&                    assistance,
                   const unordered_map<int, unordered_map<int,double>>& deltaEllA,
                   bool newIdea) {
    auto visible = [&](int id) -> bool {
        if (id == sID || id == tID) return false;
        auto [u, v] = idToPair.at(id);
        return assistance.count({v, u}) > 0;
    };

    auto ellT = [&](int u, int uP) -> double {
        for (const AdjEntry& e : graphT[u])
            if (e.v == uP) return e.l;
        return numeric_limits<double>::infinity();
    };

    unordered_map<int, vector<EdgeInfo>> adj;
    adj.reserve(prodRaw.size());

    for (auto const& [from, rawOut] : prodRaw) {
        auto& out = adj[from];

        for (int to : rawOut) {
            double lenT = 0.0, lenA = 0.0, rew = 0.0;

            if (from != sID && from != tID && to != sID && to != tID) {
                auto [u,  v ] = idToPair.at(from);
                auto [uP, vP] = idToPair.at(to);

                lenT = ellT(u, uP);
                if (!std::isfinite(lenT)) continue;

                auto itV = deltaEllA.find(v);
                if (itV == deltaEllA.end()) continue;
                auto itCost = itV->second.find(vP);
                if (itCost == itV->second.end()) continue;

                lenA = itCost->second;

                double A_from = visible(from) ? 1.0 : 0.0;
                double A_to   = visible(to)   ? 1.0 : 0.0;
                rew = 0.5 * (A_from + A_to) * lenT;

                if (newIdea)
                    rew = max(0.0, rew - lenA);

                if (newIdea && (A_to != 1 || A_from != 1))
                    continue;
            }

            out.push_back({to, lenT, lenA, rew});
        }
    }
    return adj;
}

unordered_map<int, vector<EdgeInfo>> reverseAnnotatedProd(const unordered_map<int, vector<EdgeInfo>>& adj) {
    unordered_map<int, vector<EdgeInfo>> rev;
    for (const auto& [u, edges] : adj)
        for (const auto& e : edges) {
            rev[e.to].push_back({u, e.lenT, e.lenA, e.reward});
        }
    return rev;
}

vector<int> topoSort(const unordered_map<int, vector<EdgeInfo>>& adj) {
    unordered_map<int, int> indeg;
    for (auto const& kv : adj) indeg[kv.first];
    for (auto const& kv : adj)
        for (auto const& e : kv.second) indeg[e.to] += 1;

    queue<int> Q;
    for (auto const& [v, d] : indeg) if (d == 0) Q.push(v);

    vector<int> order; order.reserve(adj.size());
    while (!Q.empty()) {
        int u = Q.front(); Q.pop(); order.push_back(u);
        auto it = adj.find(u);
        if (it == adj.end()) continue;
        for (auto const& e : it->second)
            if (--indeg[e.to] == 0) Q.push(e.to);
    }
    return order;
}