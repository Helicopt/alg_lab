#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <alg_lab/defines.hpp>
#include <alg_lab/structures.hpp>

namespace toka
{
    template <typename Nd, typename Ed, typename T>
    std::vector<T> spfa(LinkGraph<Nd, Ed> &G, int src, T zero, T inf, T (*dist)(edge_t<Ed>) = nullptr)
    {
        if (dist == nullptr)
            dist = [](edge_t<Ed> _e) -> T { return (T)_e.w; };
        std::vector<T> D(G.size(), inf);
        D[src] = zero;
        std::vector<bool> inq(G.size(), false);
        std::queue<int> q;
        q.push(src);
        inq[src] = true;
        while (!q.empty())
        {
            int x = q.front();
            q.pop();
            inq[x] = false;
            auto &node = G[x];
            for (auto &e : G.get_edges(node))
            {
                T _d = D[x] + dist(e);
                if (_d < D[e.adj])
                {
                    D[e.adj] = _d;
                    if (!inq[e.adj])
                    {
                        q.push(e.adj);
                        inq[e.adj] = true;
                    }
                }
            }
        }
        return D;
    }
} // namespace toka

#endif
