#ifndef DEF_HPP
#define DEF_HPP
#ifndef MODP7
#define MODP7 1000000007
#else
static_assert(MODP7 == 1000000007);
#endif

#ifndef ONLINE_JUDGE
#define DEBUG(tag, A, ...) std::cerr << #tag << ":\n" \
                                     << A << std::endl
#else
#define DEBUG(...)
#endif

#include <bits/stdc++.h>
namespace toka
{
    typedef long long LL;
    typedef std::pair<int, int> PA;
    typedef std::pair<LL, LL> PL;
    typedef std::pair<double, double> PF;
    typedef std::vector<int> VI;
    typedef std::vector<long long> VL;

    void redirect_in(std::string infile = "in.txt")
    {
#ifndef ONLINE_JUDGE
        freopen(infile.c_str(), "r", stdin);
#endif
    }

    void redirect_out(std::string outfile = "out.txt")
    {
#ifndef ONLINE_JUDGE
        freopen(outfile.c_str(), "w", stdout);
#endif
    }

} // namespace toka

#endif