#ifndef MATHS_HPP
#define MATHS_HPP

#include <alg_lab/defines.hpp>

namespace toka
{
    template <typename T>
    T gcd(T x, T y)
    {
        while (y)
        {
            T r = x % y;
            x = y;
            y = r;
        }
        return x;
    }

    template <typename T>
    T lcm(T x, T y)
    {
        T g = gcd(x, y);
        assert(g != 0);
        return x * y / g;
    }

    template <typename T, typename M = unsigned long>
    T qpow(T a, LL b, T one, M m = 0)
    {
        T ret = one;
        T k = a;
        while (b)
        {
            if (b & 1)
                ret = ret * k;
            k = k * k;
            b >>= 1;
            if (m > 0)
            {
                k = k % m;
                ret = ret % m;
            }
        }
        return ret;
    }

    void get_primes(const unsigned long n, VI &ret)
    {
        ret.clear();
        std::vector<bool> flags(n + 1, true);
        for (unsigned long i = 2; i <= n; ++i)
        {
            if (flags[i])
            {
                ret.emplace_back(i);
            }
            for (auto k : ret)
            {
                if ((LL)i * k > n)
                    break;
                flags[i * k] = false;
                if (i % k == 0)
                    break;
            }
        }
    }

    template <typename T, typename M>
    M divide_large(const std::vector<T> &n, M m, std::vector<T> &q, unsigned long long base = 10)
    {
        signed long L = n.size() - 1;
        unsigned long long cum = 0;
        q.clear();
        q.resize(n.size(), 0);
        while (L >= 0)
        {
            cum = cum * base + n[L];
            if (cum >= m)
            {
                q[L] = cum / m;
                cum = cum % m;
            }
            L -= 1;
        }
        L = n.size() - 1;
        while (L >= 0 && q[L] == 0)
        {
            q.pop_back();
            L -= 1;
        }
        return (M)cum;
    }

} // namespace toka

#endif