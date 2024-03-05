#include <bits/stdc++.h>
using namespace std;
#include <atcoder/all>
using namespace atcoder;

// macro
#define rep(i, n) for (LL i = 0; i < (LL)(n); i++)
#define rep1(i, n) for (LL i = 1; i < (LL)(n+1); i++)
#define PRINT_YN(b) cout << (b ? "Yes" : "No") << endl;

// type
using LL = long long;
using PII = pair<int, int>;
using PFF = pair<double, double>;
using VI = vector<int>;
using VLL = vector<LL>;
using VS = vector<string>;
using VB = vector<bool>;
using VVI = vector<vector<int>>;
using VVLL = vector<vector<LL>>;
using VVB = vector<vector<bool>>;
using VVVI = vector<vector<vector<int>>>;

struct HashPair {

    //注意 constがいる
    template<class T1, class T2>
    size_t operator()(const pair<T1, T2> &p) const {

        //first分をハッシュ化する
        auto hash1 = hash<T1>{}(p.first);

        //second分をハッシュ化する
        auto hash2 = hash<T2>{}(p.second);

        //重複しないようにハッシュ処理
        size_t seed = 0;
        seed ^= hash1 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= hash2 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};

// func
// template<class T, class K>
// bool contains(T collection, K elem) {
//   return collection.find(elem) != collection.end();
// }

void loadKernel(int i) { (void)i; }

template<class Head, class... Tail>
void loadKernel(int i, Head&& head, Tail&&... tail) {
  cin >> head[i];
  loadKernel(i, std::forward<Tail>(tail)...);
}

template<class... Dst>
void loadVec(const int n, Dst&&... dst) {
  rep(i, n) loadKernel(i, std::forward<Dst>(dst)...);
}

void loadVar() {}

template<class Head, class... Tail>
void loadVar(Head&& head, Tail&&... tail) {
  cin >> head;
  loadVar(std::forward<Tail>(tail)...);
}

string notationKernel() {
  return "";
}

template<class Head, class... Tail>
string notationKernel(Head&& head, Tail&&... tail) {
  string s = "," + to_string(head);
  return s + notation(std::forward<Tail>(tail)...);
}

template<class Head, class... Tail>
string notation(Head&& head, Tail&&... tail) {
  string s = to_string(head);
  return s + notation(std::forward<Tail>(tail)...);
}

template<class T>
void printVec(const T& v) {
  bool first = true;
  for(auto itr = v.begin(); itr != v.end(); itr++) {
    if(!first) cout << " ";
    first = false;
    cout << *itr;
  }
  cout <<endl;
}

template <typename T>
T Ceil(T a,T b) {
  return (a + b - 1) / b;
}

// 素因数分解 (1~n)
/* usage */
// PrimeFact pf(n);
// pf.get(x)
template <typename T>
struct PrimeFact {  
    vector<T> spf;
    PrimeFact(T N) { init(N); }
    void init(T N) { // 前処理。spf を求める
        spf.assign(N + 1, 0);
        iota(spf.begin(), spf.end(), 0);
        for (T i = 2; i * i <= N; i++) {
            if (spf[i] != i) continue;
            for (T j = i * i; j <= N; j += i) {
                if (spf[j] == j) spf[j] = i;
            }
        }
    }
    map<T, T> get(T n) { // nの素因数分解を求める
        map<T, T> m;
        while (n != 1) {
            m[spf[n]]++;
            n /= spf[n];
        }
        return m;
    }
};

/*  prime_factor(n)
    入力：整数 n
    出力：nの素因数分解
    計算量：O(√n)前後
*/
template <typename T>
map<T, T> prime_factor(T n) {
    map<T, T> ret;
    for (T i = 2; i * i <= n; i++) {
        if (n % i != 0) continue;
        T tmp = 0;
        while (n % i == 0) {
            tmp++;
            n /= i;
        }
        ret[i] = tmp;
    }
    if (n != 1) ret[n] = 1;
    return ret;
}


// 素数の集合を得る
template <typename T>
struct Prime {
  T N;
  list<T> primes;

  Prime(T n) { init(n); }

  void init(T n) {
    N = n;
    primes.clear();

    vector<bool> sieve(N+1, true);
    for(T i = 2; i <= N; i++) {
      if(!sieve[i]) continue;
      primes.push_back(i);
      for(T j = i * 2; j <= N; j+=i) sieve[j] = false;
    }
  }

  list<T>& get() { return primes; }
};

template <typename T>
class IDTable {
public:
  map<T, T> id_table;
  map<T, T> id_table_;
  IDTable() {}
  IDTable(const vector<T>& v) {
    rep(i, v.size()) {
      id_table[i] = v[i];    
      id_table_[v[i]] = i;
    }
  }
  T forward(const T n) { return id_table[n]; }
  T backward(const T n) { return id_table_[n]; }
};
////////////////////////////////////////////////////

// int Calc(VI X, VI Y, VI Z) {
//   int v = 1;
//   v *= max(0, *max_element(X.begin(), X.end()) - *min_element(X.begin(), X.end()));
//   v *= max(0, *max_element(Y.begin(), Y.end()) - *min_element(Y.begin(), Y.end()));
//   v *= max(0, *max_element(Z.begin(), Z.end()) - *min_element(Z.begin(), Z.end()));
// 
//   return v;
// }

int f(int a1, int b1, int c1, int a2, int b2, int c2) {
    int res = 1;
    res *= max(0, min(a1, a2) + 7 - max(a1, a2));
    res *= max(0, min(b1, b2) + 7 - max(b1, b2));
    res *= max(0, min(c1, c2) + 7 - max(c1, c2));
    return res;
}
int f(int a1, int b1, int c1, int a2, int b2, int c2, int a3, int b3, int c3) {
    int res = 1;
    res *= max(0, min({a1, a2, a3}) + 7 - max({a1, a2, a3}));
    res *= max(0, min({b1, b2, b3}) + 7 - max({b1, b2, b3}));
    res *= max(0, min({c1, c2, c3}) + 7 - max({c1, c2, c3}));
    return res;
}

void solve() {
  int V1, V2, V3;
  loadVar(V1, V2, V3);

  int x1, y1, z1;
  x1 = y1 = z1 = 0;
  for(int x2 = -7; x2 <= 7; x2++) {
    for(int y2 = -7; y2 <= 7; y2++) {
      for(int z2 = -7; z2 <= 7; z2++) {
        for(int x3 = -7; x3 <= 7; x3++) {
          for(int y3 = -7; y3 <= 7; y3++) {
            for(int z3 = -7; z3 <= 7; z3++) {
              // int v3 = Calc({x1,x2,x3}, {y1,y2,y3}, {z1, z2, z3});
              // int v2 = Calc({x1,x2}, {y1,y2}, {z1, z2}) + 
              //          Calc({x1,x3}, {y1,y3}, {z1, z3}) +
              //          Calc({x2,x3}, {y2,y3}, {z2, z3}) - v3 * 3;
              // int v1 = 7 * 7 * 7 * 3 - V2 * 2 - V3 * 3;
              int v3 = f(0, 0, 0, x2, y2, z2, x3, y3, z3);
              int v2 = f(0, 0, 0, x2, y2, z2) + f(0, 0, 0, x3, y3, z3) + f(x2, y2, z2, x3, y3, z3) - v3 * 3;
              int v1 = 3 * 7 * 7 * 7 - v2 * 2 - v3 * 3;
              if(v1 == V1 && v2 == V2 && v3 == V3) {
                PRINT_YN(true);
                cout << x1 << " " << y1 << " " << z1 << " "
                     << x2 << " " << y2 << " " << z2 << " "
                     << x3 << " " << y3 << " " << z3 << endl;
                return;
              }
            }
          }
        }
      }
    }
  }
  PRINT_YN(false);
}

int main(int argc, char *argv[]) {
#ifdef ONLINE_JUDGE 
  const int n_testcase = 1;  // Don't change here!!
#else
  int tmp = 1;
  cin >> tmp;
  const int n_testcase = tmp;  // Don't change here!!
#endif

  rep(i, n_testcase) {
    solve();
  }

  return 0;
}