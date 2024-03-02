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

// 素因数分解
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

////////////////////////////////////////////////////
int N, M;
VI U(M), V(M);
vector<list<int>> D;
VI colors;

int B = 0, R = 0;

bool checkBiominal(const int i, const int color) {
  colors[i] = color;
  switch(color) {
    case  1: B++; break;
    case -1: R++; break;
  }

  for(auto e: D[i]) {
    if(colors[e] == color) return false;
    if(colors[e] != 0) continue;
    if(!checkBiominal(e, -color)) return false;
  }

  return true;
}

void solve() {
  loadVar(N, M);
  U.resize(M, 0); V.resize(M, 0); loadVec(M, U, V);
  D.resize(N);
  colors.resize(N, 0);

  dsu S(N);
  rep(i, M) {
    int u = U[i] - 1;
    int v = V[i] - 1;
    D[u].push_back(v);
    D[v].push_back(u);
    S.merge(u, v);
  }

  const auto& G = S.groups();
  const int gs = S.groups().size();
  LL ans = 0;

  for(int i = 0; i < gs -1; i++) {
    for(int j = i+1; j < gs; j++) {
      ans += G[i].size() * G[j].size();
    }
  }

  rep(i, gs) {
    B = 0; R = 0;
    if(!checkBiominal(G[i][0], 1)) {
      cout << 0 << endl;
      return;
    }
    ans += B * R;
    
    LL s_edge = 0;
    for(auto e: G[i]) {
      s_edge += D[e].size();
    }
    ans -= s_edge / 2;
  }

  cout << ans << endl;
}

int main() {
#ifdef ONLINE_JUDGE 
  const int n_testcase = 1;  // Don't change here!!
#else
  const int n_testcase = 3;  // Change here for test cases.
#endif

  rep(i, n_testcase) {
    solve();
#ifndef ONLINE_JUDGE 
    cout << endl;
#endif
  }

  return 0;
}