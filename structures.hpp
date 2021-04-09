#ifndef STRUCT_HPP
#define STRUCT_HPP

#include <alg_lab/defines.hpp>

namespace toka
{
    template <typename T>
    class Matrix
    {
    private:
        unsigned int _h, _w;
        std::shared_ptr<T[]> _d;
        void _allocate(unsigned int h, unsigned int w)
        {
            _h = h;
            _w = w;
            _d = static_cast<std::shared_ptr<T[]>>(new T[_h * _w]);
        }

    public:
        Matrix(unsigned int h, unsigned int w)
        {
            _allocate(h, w);
        }
        Matrix(unsigned int h, unsigned int w, const T &fill_)
        {
            _allocate(h, w);
            for (size_t i = 0; i < _h * _w; ++i)
                _d[i] = fill_;
        }

        Matrix copy()
        {
            Matrix ret = Matrix(_h, _w, 0);
            for (size_t i = 0; i < _h * _w; ++i)
            {
                ret._d[i] = _d[i];
            }
            // memcpy(ret._d.get(), _d.get(), _h * _w * sizeof(T));
            return ret;
        }

        T &operator[](PA index)
        {
            return _d[index.first * _w + index.second];
        }

        friend std::ostream &operator<<(std::ostream &stream, Matrix A)
        {
            VI width;
            for (unsigned int i = 0; i < A._w; ++i)
            {
                size_t tmp = 0;
                for (unsigned int j = 0; j < A._h; ++j)
                {
                    T k = A[PA(j, i)];
                    std::stringstream ts;
                    ts << k;
                    std::string s = ts.str();
                    tmp = std::max(tmp, s.length());
                }
                width.emplace_back(tmp);
            }
            for (unsigned int i = 0; i < A._h; ++i)
            {
                stream << "[";
                for (unsigned int j = 0; j < A._w; ++j)
                {
                    T k = A[PA(i, j)];
                    std::stringstream ts;
                    ts << k;
                    std::string s = ts.str();
                    while (s.length() < width[j] + 2)
                        s = " " + s;
                    stream << s;
                }
                stream << " ]" << std::endl;
            }
            stream << "(Matrix of " << A._h << "x" << A._w << ")";
            return stream;
        }

        Matrix ONE()
        {
            Matrix ret = Matrix(_h, _w, 0);
            for (unsigned int i = 0; i < std::min(_h, _w); ++i)
            {
                ret[PA(i, i)] = 1;
            }
            return ret;
        }

        static Matrix ONE(unsigned int n)
        {
            Matrix ret = Matrix(n, n, 0);
            for (unsigned int i = 0; i < n; ++i)
            {
                ret[PA(i, i)] = 1;
            }
            return ret;
        }

        Matrix _add(Matrix other)
        {
            assert(_h == other._h && _w == other._w);
            for (unsigned int i = 0; i < _h; ++i)
            {
                for (unsigned int j = 0; j < _w; ++j)
                {
                    operator[](PA(i, j)) += other[PA(i, j)];
                }
            }
            return *this;
        }

        Matrix operator+(Matrix other)
        {
            assert(_h == other._h && _w == other._w);
            auto ret = copy();
            ret._add(other);
            return ret;
        }

        Matrix _mul(Matrix other, Matrix dst)
        {
            assert(_h == other._h && _w == other._w && _w == other._w);
            assert(dst._h = _h && dst._w == other._w);
            for (unsigned int i = 0; i < _h; ++i)
            {
                for (unsigned int j = 0; j < other._w; ++j)
                {
                    dst[PA(i, j)] = 0;
                    for (unsigned int k = 0; k < _w; ++k)
                    {
                        dst[PA(i, j)] += operator[](PA(i, k)) * other[PA(k, j)];
                    }
                }
            }
            return dst;
        }

        Matrix operator*(Matrix other)
        {
            assert(_w == other._h);
            auto ret = Matrix(_h, other._w, 0);
            this->_mul(other, ret);
            return ret;
        }

        Matrix _mod(T mod)
        {
            for (unsigned int i = 0; i < _h; ++i)
            {
                for (unsigned int j = 0; j < _w; ++j)
                {
                    operator[](PA(i, j)) %= mod;
                }
            }
            return *this;
        }

        Matrix operator%(T mod)
        {
            auto ret = copy();
            ret._mod(mod);
            return ret;
        }
    };

    template <typename T>
    struct edge_t
    {
        int adj, nxt;
        T w;
    };

    template <typename T>
    struct node_t
    {
        int i;
        T v;
    };

    template <typename Nd, typename Ed>
    class LinkGraph
    {
    private:
        int _n, _m;
        int _tail;
        std::shared_ptr<node_t<Nd>[]> _g;
        std::shared_ptr<edge_t<Ed>[]> _e;

    public:
        class EdgeSet
        {
        private:
            std::shared_ptr<edge_t<Ed>[]> _eptr;
            int _pos;

        public:
            class Iterator
            {
            private:
                int _pos;
                std::shared_ptr<edge_t<Ed>[]> _eptr;

            public:
                Iterator(std::shared_ptr<edge_t<Ed>[]> e_ptr, int pos)
                {
                    _eptr = e_ptr;
                    _pos = pos;
                }
                Iterator &operator++()
                {
                    _pos = _eptr[_pos].nxt;
                    return *this;
                }
                bool operator!=(const Iterator &other)
                {
                    return _pos != other._pos || _eptr.get() != other._eptr.get();
                }
                edge_t<Ed> &operator*()
                {
                    return _eptr[_pos];
                }
            };
            EdgeSet(std::shared_ptr<edge_t<Ed>[]> e_ptr, int pos)
            {
                _eptr = e_ptr;
                _pos = pos;
            }
            Iterator begin()
            {
                return Iterator(_eptr, _pos);
            }
            Iterator end()
            {
                return Iterator(_eptr, 0);
            }
        };
        class Iterator
        {
        private:
            int _pos = 0;
            LinkGraph *_ptr;

        public:
            Iterator(LinkGraph *G, int pos)
            {
                _ptr = G;
                _pos = pos;
            }

            Iterator &operator++()
            {
                _pos++;
                return *this;
            }
            bool operator!=(const Iterator &other)
            {
                return _ptr != other._ptr || _pos != other._pos;
            }
            node_t<Nd> &operator*()
            {
                return _ptr->_g[_pos];
            }
        };

        LinkGraph(const int &n, const int &m)
        {
            _n = n;
            _m = m;
            _g = static_cast<std::shared_ptr<node_t<Nd>[]>>(new node_t<Nd>[n]);
            _tail = 0;
            _e = static_cast<std::shared_ptr<edge_t<Ed>[]>>(new edge_t<Ed>[m + 1]);
            for (int i = 0; i < n; ++i)
                _g[i].i = 0;
        }

        edge_t<Ed> &add_edge(int x, int y)
        {
            ++_tail;
            _e[_tail].adj = y;
            _e[_tail].nxt = _g[x].i;
            _g[x].i = _tail;
            return _e[_tail];
        }

        edge_t<Ed> &add_edge(int x, int y, Ed w)
        {
            edge_t<Ed> &e = add_edge(x, y);
            e.w = w;
            return e;
        }

        Iterator begin()
        {
            return Iterator(this, 0);
        }

        Iterator end()
        {
            return Iterator(this, _n);
        }

        friend std::ostream &operator<<(std::ostream &stream, LinkGraph A)
        {
            stream << "[*LinkGraph* of " << A._n << " nodes and " << A._m << " edges (capacity)]";
            return stream;
        }

        node_t<Nd> &operator[](int pos)
        {
            return _g[pos];
        }

        EdgeSet get_edges(const node_t<Nd> &node)
        {
            return EdgeSet(_e, node.i);
        }

        size_t size()
        {
            return _n;
        }
    };

    template <typename T>
    class SegTree
    {
    private:
        int _l, _r;
        std::shared_ptr<T[]> _d;

        virtual void _update(int l, int r, T &d, const T &c) = 0;
        virtual void _push(int l, int r, T &d, const T &c) = 0;
        virtual void _clear(int l, int r, T &d) = 0;
        virtual void _pull(int l, int r, T &d, const T &lv, const T &rv) = 0;
        virtual T _merge(int l, int r, const T &lv, const T &rv) = 0;

        void _insert(int u, int l, int r, int x, int y, const T &c)
        {
            if (l > r || x > y || x > r || y < l)
                return;
            x = std::max(x, l);
            y = std::min(y, r);
            if (l == x && r == y)
            {
                _update(l, r, _d[u], c);
                return;
            }
            int ll = u << 1, rr = ll + 1;
            int m = l + r >> 1;
            _push(l, m, _d[ll], _d[u]);
            _push(m + 1, r, _d[rr], _d[u]);
            _clear(l, r, _d[u]);
            if (y <= m)
                _insert(ll, l, m, x, y, c);
            else if (x > m)
                _insert(rr, m + 1, r, x, y, c);
            else
            {
                _insert(ll, l, m, x, m, c);
                _insert(rr, m + 1, r, m + 1, y, c);
            }
            _pull(l, r, _d[u], _d[ll], _d[rr]);
        }

        T _find(int u, int l, int r, int x, int y, const T &c)
        {
            if (l > r || x > y || x > r || y < l)
                return c;
            x = std::max(x, l);
            y = std::min(y, r);
            if (l == x && r == y)
            {
                return _d[u];
            }
            int ll = u << 1, rr = ll + 1;
            int m = l + r >> 1;
            _push(l, m, _d[ll], _d[u]);
            _push(m + 1, r, _d[rr], _d[u]);
            _clear(l, r, _d[u]);
            T ret;
            if (y <= m)
                ret = _find(ll, l, m, x, y, c);
            else if (x > m)
                ret = _find(rr, m + 1, r, x, y, c);
            else
            {
                auto lv = _find(ll, l, m, x, m, c);
                auto rv = _find(rr, m + 1, r, m + 1, y, c);
                ret = _merge(l, r, lv, rv);
            }
            _pull(l, r, _d[u], _d[ll], _d[rr]);
            return ret;
        }

    public:
        SegTree(int l, int r)
        {
            assert(l <= r);
            _l = l;
            _r = r;
            _d = static_cast<std::shared_ptr<T[]>>(new T[(_r - _l + 1) << 2]);
        }

        size_t size()
        {
            return _r - _l + 1;
        }

        void init(const T &c)
        {
            for (int i = 0; i < (size() << 2); ++i)
            {
                _d[i] = c;
            }
        }

        void insert(int l, int r, const T &c)
        {
            _insert(1, _l, _r, l, r, c);
        }

        void insert(int l, const T &c)
        {
            insert(l, l, c);
        }

        T find(int l, int r, const T &c)
        {
            return _find(1, _l, _r, l, r, c);
        }

        T find(int l, const T &c)
        {
            return find(l, l, c);
        }
    };

    struct SumSegTreeNode
    {
        int tag;
        int sum;
    };

    template <typename T>
    class _SumSegTree : public toka::SegTree<T>
    {
    private:
        void _update(int l, int r, T &d, const T &c)
        {
            d.tag = c.tag;
            d.sum = d.tag * (r - l + 1);
        }

        void _push(int l, int r, T &d, const T &c)
        {
            if (c.tag != -1)
            {
                d.tag = c.tag;
                d.sum = d.tag * (r - l + 1);
            }
        }
        void _clear(int l, int r, T &d)
        {
            d.tag = -1;
        }
        void _pull(int l, int r, T &d, const T &lv, const T &rv)
        {
            d.sum = lv.sum + rv.sum;
        }
        T _merge(int l, int r, const T &lv, const T &rv)
        {
            T ret;
            ret.sum = lv.sum + rv.sum;
            ret.tag = -1;
            return ret;
        }

    public:
        using toka::SegTree<T>::SegTree;
    };

    class SumSegTree : public _SumSegTree<SumSegTreeNode>
    {
    public:
        using _SumSegTree::_SumSegTree;
    };

    template <typename T, size_t maxsize>
    struct TrieNode
    {
        size_t MAXSIZE = maxsize;
        T data;
        size_t cnt, ends;
        TrieNode<T, maxsize> *ch[maxsize];
        TrieNode()
        {
            memset(ch, 0, sizeof(ch));
            reset();
        }
        void reset()
        {
            cnt = 0;
            ends = 0;
            free();
        }
        void free()
        {
            for (size_t i = 0; i < maxsize; ++i)
            {
                if (ch[i] != nullptr)
                {
                    ch[i]->free();
                }
            }
            memset(ch, 0, sizeof(ch));
        }
        ~TrieNode()
        {
            free();
        }
    };

    template <typename T, size_t maxsize>
    class TriePath
    {
    public:
        class Iterator
        {
        private:
            TrieNode<T, maxsize> *_pos;
            const VI *_seq;
            int _ind;

        public:
            Iterator(const VI &seq, TrieNode<T, maxsize> *root, int index = -1)
            {
                _ind = index;
                _pos = root;
                _seq = &seq;
            }
            Iterator &operator++()
            {
                _ind += 1;
                auto ind = (_ind < _seq->size()) ? (*_seq)[_ind] : -1;
                _pos = (ind >= 0) ? _pos->ch[ind] : nullptr;
                return *this;
            }
            bool operator!=(const Iterator &other)
            {
                return _pos != other._pos || _ind != other._ind;
            }
            T &operator*()
            {
                return _pos->data;
            }
        };
        TriePath(const VI &path, TrieNode<T, maxsize> *root)
        {
            seq = path;
            _root = root;
        }
        Iterator begin()
        {
            return Iterator(seq, _root);
        }
        Iterator end()
        {
            return Iterator(seq, nullptr, seq.size());
        }

    private:
        VI seq;
        TrieNode<T, maxsize> *_root;
    };

    template <typename T, size_t maxsize>
    class Trie
    {
        TrieNode<T, maxsize> root;

    public:
        void insert(const VI &seq, const T &data)
        {
            TrieNode<T, maxsize> *ptr = &root;
            int nxt = 0;
            while (nxt < seq.size())
            {
                auto ind = seq[nxt];
                if (ind < 0 || ind >= maxsize)
                {
                    throw std::runtime_error("index exceed error");
                }
                if (ptr->ch[ind] == nullptr)
                {
                    ptr->ch[ind] = new TrieNode<T, maxsize>();
                }
                ptr->ch[ind]->cnt += 1;
                if (nxt == seq.size() - 1)
                {
                    ptr->ch[ind]->ends += 1;
                }
                insert_update(nxt, ind, nxt == seq.size() - 1, ptr->ch[ind]->data, data);
                ptr = ptr->ch[ind];
                nxt += 1;
            }
        }

        void remove(const VI &seq, const T &data)
        {
            TrieNode<T, maxsize> *ptr = &root;
            int nxt = 0;
            while (nxt < seq.size())
            {
                auto ind = seq[nxt];
                if (ind < 0 || ind >= maxsize)
                {
                    throw std::runtime_error("index exceed error");
                }
                if (ptr->ch[ind] == nullptr)
                {
                    return;
                }
                ptr->ch[ind]->cnt -= 1;
                if (nxt == seq.size() - 1)
                {
                    ptr->ch[ind]->ends -= 1;
                }
                remove_update(nxt, ind, nxt == seq.size() - 1, ptr->ch[ind]->data, data);
                ptr = ptr->ch[ind];
                nxt += 1;
            }
        }

        TriePath<T, maxsize> search(const VI &path)
        {
            return TriePath<T, maxsize>(path, &root);
        }

    private:
        void insert_update(int ind, int val, bool is_end, T &node, const T &data)
        {
            node += data;
        }
        void remove_update(int ind, int val, bool is_end, T &node, const T &data)
        {
            node -= data;
        }
    };

    template <typename T, size_t maxsize>
    class ACAutomata : Trie<T, maxsize>
    {
    public:
        ACAutomata()
        {
        }
        void build(const std::vector<VI> &words)
        {
        }
        void build(const std::vector<std::string> &words)
        {
        }

    private:
    };

} // namespace toka

#endif