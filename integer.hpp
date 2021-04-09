#ifndef INTEGER_HPP
#define INTEGER_HPP

#include <alg_lab/defines.hpp>
#include <alg_lab/maths.hpp>

#define CHECK_INT_TYPE(T_, X_) typename std::enable_if< \
    std::is_same<char, T_>::value ||                    \
        std::is_same<unsigned char, T_>::value ||       \
        std::is_same<short, T_>::value ||               \
        std::is_same<unsigned short, T_>::value ||      \
        std::is_same<int, T_>::value ||                 \
        std::is_same<unsigned int, T_>::value ||        \
        std::is_same<long, T_>::value ||                \
        std::is_same<unsigned long, T_>::value ||       \
        std::is_same<long long, T_>::value ||           \
        std::is_same<unsigned long long, T_>::value,    \
    X_>::type

namespace toka
{

    class Integer
    {
    public:
        Integer(const long long &n)
        {
            if (n >= 0)
                pos = 1;
            else
                pos = -1;
            auto nn = (unsigned long long)abs(n);
            data.clear();
            data.push_back(nn & 0xFFFFFFFFU);
            data.push_back(nn >> 32);
            clear_zeros();
        }

        Integer(const unsigned long long &n)
        {
            pos = 1;
            auto nn = (unsigned long long)n;
            data.clear();
            data.push_back(nn & 0xFFFFFFFFU);
            data.push_back(nn >> 32);
            clear_zeros();
        }

        Integer(const long &n)
        {
            if (n >= 0)
                pos = 1;
            else
                pos = -1;
            auto nn = (unsigned long long)abs(n);
            data.clear();
            data.push_back(nn & 0xFFFFFFFFU);
            data.push_back(nn >> 32);
            clear_zeros();
        }

        Integer(const unsigned long &n)
        {
            pos = 1;
            auto nn = (unsigned long long)n;
            data.clear();
            data.push_back(nn & 0xFFFFFFFFU);
            data.push_back(nn >> 32);
            clear_zeros();
        }

        Integer(const int &n)
        {
            if (n >= 0)
                pos = 1;
            else
                pos = -1;
            auto nn = (unsigned long long)abs(n);
            data.clear();
            data.push_back(nn & 0xFFFFFFFFU);
            data.push_back(nn >> 32);
            clear_zeros();
        }

        Integer(const unsigned int &n)
        {
            pos = 1;
            auto nn = (unsigned long long)n;
            data.clear();
            data.push_back(nn & 0xFFFFFFFFU);
            data.push_back(nn >> 32);
            clear_zeros();
        }

        Integer(const std::string &n)
        {
            int start = 0;
            if (n[0] == '-')
            {
                pos = -1;
                start = 1;
            }
            else
                pos = 1;
            VI k;
            auto l = n.size();
            for (size_t i = 0; i < l - start; ++i)
            {
                k.push_back(n[l - 1 - i] - '0');
            }
            data.clear();
            VI q;
            while (k.size())
            {
                auto r = divide_large(k, 1llu << 32, q);
                data.push_back(r);
                k.swap(q);
            }
            clear_zeros();
        }

        Integer &_sub(const Integer &other)
        {
            *this += other.neg();
            return *this;
        }

        Integer sub(const Integer &other) const
        {
            return other.neg() += *this;
        }

        Integer &_add(const Integer &other)
        {
            if (pos == other.pos)
            {
                size_t n = std::max(data.size(), other.data.size());
                size_t m = std::min(data.size(), other.data.size());
                data.resize(n + 1, 0);
                auto cum = 0ull;
                for (size_t i = 0; i < m; ++i)
                {
                    cum += data[i];
                    cum += other.data[i];
                    data[i] = cum & 0xFFFFFFFFU;
                    cum >>= 32;
                }
                if (m >= other.data.size())
                    for (size_t i = m; i < n; i++)
                    {
                        cum += data[i];
                        data[i] = cum & 0xFFFFFFFFU;
                        cum >>= 32;
                    }
                else
                    for (size_t i = m; i < n; i++)
                    {
                        cum += other.data[i];
                        data[i] = cum & 0xFFFFFFFFU;
                        cum >>= 32;
                    }
                if (cum > 0)
                {
                    data[n] = cum;
                }
            }
            else
            {
                size_t n = std::max(data.size(), other.data.size());
                data.resize(n + 1, 0);
                const Integer *a, *b;
                if (pos < 0)
                {
                    a = &other;
                    b = this;
                }
                else
                {
                    b = &other;
                    a = this;
                }
                auto cum = 0ull;
                for (size_t i = 0; i < n; ++i)
                {
                    if (((unsigned long long)(a->data[i])) >= cum + b->data[i])
                    {
                        data[i] = a->data[i] - b->data[i] - cum;
                        cum = 0;
                    }
                    else
                    {
                        data[i] = (1ull << 32) + a->data[i] - b->data[i] - cum;
                        cum = 1;
                    }
                }
                if (cum > 0)
                {
                    pos = -1;
                    cum = 1;
                    for (size_t i = 0; i < n; ++i)
                    {
                        cum += ~data[i];
                        data[i] = cum & 0xFFFFFFFF;
                        cum >>= 32;
                    }
                    if (cum > 0)
                        data[n] = 1;
                }
                else
                {
                    pos = 1;
                }
            }
            clear_zeros();
            return *this;
        }

        Integer add(const Integer &other) const
        {
            Integer ret = *this;
            return ret += other;
        }

        Integer mul(const Integer &other) const
        {
            auto ret = *this;
            return ret._mul(other);
        }

        Integer &_mul(const Integer &other)
        {
            pos *= other.pos;
            Integer ret("0");
            ret.data.resize(data.size() * other.data.size() + 3, 0);
            for (size_t i = 0; i < data.size(); ++i)
            {
                for (size_t j = 0; j < other.data.size(); ++j)
                {
                    auto cum = 0llu;
                    auto f = 1llu * data[i] * other.data[j];
                    cum += ret.data[i + j];
                    cum += f & 0xFFFFFFFFU;
                    ret.data[i + j] = cum & 0xFFFFFFFFU;
                    cum >>= 32;
                    cum += ret.data[i + j + 1];
                    cum += f >> 32;
                    ret.data[i + j + 1] = cum & 0xFFFFFFFFU;
                    cum >>= 32;
                    size_t p = i + j + 2;
                    while (cum > 0)
                    {
                        cum += ret.data[p];
                        ret.data[p] = cum & 0xFFFFFFFFU;
                        cum >>= 32;
                        p += 1;
                    }
                }
            }
            data.swap(ret.data);
            clear_zeros();
            return *this;
        }

        Integer div(const Integer &other) const
        {
            auto ret = *this;
            return ret._div(other);
        }

        void divmod(const Integer &other, Integer &q, Integer &r) const
        {
            if (other.data.size() == 0)
            {
                throw std::runtime_error("div zero error");
            }
            if (other.data.size() == 1)
            {
                auto _r = divide_large(data, other.data[0], q.data, 1ull << 32);
                q.pos = pos * other.pos;
                r.data.resize(1, 0);
                r.data[0] = _r;
                r.pos = pos;
                r.clear_zeros();
            }
            else
            {
                signed long L = data.size() - 1;
                Integer cum(0);
                q.data.clear();
                q.data.resize(data.size(), 0);
                while (L >= 0)
                {
                    cum *= 1ull << 32;
                    cum += data[L];
                    if (cum >= other)
                    {
                        unsigned long long lb = 0;
                        unsigned long long rb = 1ull << 32;
                        while (lb + 1 < rb)
                        {
                            auto m = lb + rb >> 1;
                            if (m * other <= cum)
                                lb = m;
                            else
                                rb = m;
                        }
                        q.data[L] = lb;
                        cum -= lb * other;
                    }
                    L -= 1;
                }
                L = data.size() - 1;
                while (L >= 0 && q.data[L] == 0)
                {
                    q.data.pop_back();
                    L -= 1;
                }
                q.pos = pos * other.pos;
                r.data.swap(cum.data);
                r.pos = pos;
            }
        }

        Integer &_div(const Integer &other)
        {
            Integer q("0"), r("0");
            divmod(other, q, r);
            pos = q.pos;
            data.swap(q.data);
            return *this;
        }

        Integer mod(const Integer &other) const
        {
            auto ret = *this;
            return ret._mod(other);
        }

        Integer &_mod(const Integer &other)
        {
            Integer q("0"), r("0");
            divmod(other, q, r);
            pos = r.pos;
            data.swap(r.data);
            return *this;
        }

        bool equal(const Integer &other) const
        {
            if (pos != other.pos || data.size() != other.data.size())
                return false;
            for (size_t i = 0; i < data.size(); ++i)
                if (data[i] != other.data[i])
                    return false;
            return true;
        }

        bool lessequal(const Integer &other) const
        {
            return less(other) || equal(other);
        }

        bool less(const Integer &other) const
        {
            if (pos < 0 && other.pos > 0)
                return true;
            if (pos > 0 && other.pos < 0)
                return false;
            bool inv = false;
            if (pos < 0)
                inv = true;
            if (data.size() < other.data.size())
                return inv ^ true;
            if (data.size() > other.data.size())
                return inv ^ false;
            auto l = data.size();
            for (size_t i = 1; i <= l; ++i)
            {
                if (data[l - i] < other.data[l - i])
                {
                    return inv ^ true;
                }
                if (data[l - i] > other.data[l - i])
                {
                    return inv ^ false;
                }
            }
            return false;
        }

        Integer neg() const
        {
            auto ret = *this;
            ret.pos = -ret.pos;
            return ret;
        }

        Integer operator-() const
        {
            return this->neg();
        }

        Integer _neg()
        {
            pos = -pos;
            return *this;
        }

        template <typename T>
        Integer operator+(const T &other) const
        {
            return add(other);
        }

        template <typename T>
        Integer &operator+=(const T &other)
        {
            return _add(other);
        }

        template <typename T>
        Integer operator-(const T &other) const
        {
            return sub(other);
        }

        template <typename T>
        Integer &operator-=(const T &other)
        {
            return _sub(other);
        }

        template <typename T>
        Integer operator*(const T &other) const
        {
            return mul(other);
        }

        template <typename T>
        Integer &operator*=(const T &other)
        {
            return _mul(other);
        }

        template <typename T>
        Integer operator/(const T &other) const
        {
            return div(other);
        }

        template <typename T>
        Integer &operator/=(const T &other)
        {
            return _div(other);
        }

        template <typename T>
        Integer operator%(const T &other) const
        {
            return mod(other);
        }

        template <typename T>
        Integer &operator%=(const T &other)
        {
            return _mod(other);
        }

        template <typename T>
        bool operator==(const T &other) const
        {
            return equal(other);
        }

        template <typename T>
        bool operator!=(const T &other) const
        {
            return !equal(other);
        }

        template <typename T>
        bool operator<(const T &other) const
        {
            return less(other);
        }

        template <typename T>
        bool operator<=(const T &other) const
        {
            return lessequal(other);
        }

        template <typename T>
        bool operator>(const T &other) const
        {
            return !lessequal(other);
        }

        template <typename T>
        bool operator>=(const T &other) const
        {
            return !less(other);
        }

        template <typename T>
        friend CHECK_INT_TYPE(T, Integer)
        operator+(const T x, const Integer &y)
        {
            return Integer(x)._add(y);
        }

        template <typename T>
        friend CHECK_INT_TYPE(T, Integer)
        operator-(const T x, const Integer &y)
        {
            return Integer(x)._sub(y);
        }

        template <typename T>
        friend CHECK_INT_TYPE(T, Integer)
        operator*(const T x, const Integer &y)
        {
            return Integer(x)._mul(y);
        }

        template <typename T>
        friend CHECK_INT_TYPE(T, Integer)
        operator/(const T x, const Integer &y)
        {
            return Integer(x)._div(y);
        }

        template <typename T>
        friend CHECK_INT_TYPE(T, Integer)
        operator%(const T x, const Integer &y)
        {
            return Integer(x)._mod(y);
        }

        template <typename T>
        friend CHECK_INT_TYPE(T, bool)
        operator==(const T x, const Integer &y)
        {
            return Integer(x).equal(y);
        }

        template <typename T>
        friend CHECK_INT_TYPE(T, bool)
        operator!=(const T x, const Integer &y)
        {
            return !(Integer(x).equal(y));
        }

        template <typename T>
        friend CHECK_INT_TYPE(T, bool)
        operator<(const T x, const Integer &y)
        {
            return Integer(x).less(y);
        }

        template <typename T>
        friend CHECK_INT_TYPE(T, bool)
        operator<=(const T x, const Integer &y)
        {
            return Integer(x).lessequal(y);
        }

        template <typename T>
        friend CHECK_INT_TYPE(T, bool)
        operator>(const T x, const Integer &y)
        {
            return !(Integer(x).lessequal(y));
        }

        template <typename T>
        friend CHECK_INT_TYPE(T, bool)
        operator>=(const T x, const Integer &y)
        {
            return !(Integer(x).less(y));
        }

        std::string to_string(unsigned long base = 10)
        {
            clear_zeros();
            std::string s = "";
            auto tmp = data;
            std::vector<unsigned long> q;
            while (tmp.size())
            {
                auto r = divide_large(tmp, base, q, 1ull << 32);
                auto c = get_char(r);
                s = c + s;
                tmp.swap(q);
            }
            if (s.size() == 0)
                s = "0";
            else if (pos < 0)
                s = '-' + s;
            return s;
        }

        friend std::ostream &operator<<(std::ostream &stream, Integer A)
        {
            stream << A.to_string();
            return stream;
        }

    private:
        std::vector<unsigned long> data;
        int pos;
        char get_char(unsigned int x)
        {
            if (x < 10)
            {
                return '0' + x;
            }
            else
            {
                return 'a' - 10 + x;
            }
        }
        void clear_zeros()
        {
            while (data.size() > 0 && data.back() == 0)
                data.pop_back();
            if (data.size() == 0 && pos < 0)
                pos = 1;
        }
    };
} // namespace toka

#endif
