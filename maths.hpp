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

    template <typename T>
    T qpow(T a, LL b, T one, int m = 0)
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

    void get_primes(const int n, VI &ret)
    {
        ret.clear();
        std::vector<bool> flags(n + 1, true);
        for (int i = 2; i <= n; ++i)
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
} // namespace toka

#endif