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
////////////////////////////////////////////////////

int h, w;
VVI A;

int flip(bool flip, int elem) {
  if(!flip) return elem;
  return 1 - elem;
}

bool check_elem(const bool flip1, const bool flip2, const bool flip3, const int i, const int j) {
  const int e = flip(flip1, A[i][j]);
  bool flag = false;
  if(i > 0) flag |= (e == flip(flip3, A[i-1][j]));
  if(i < h-1) flag |= (e == flip(flip2, A[i+1][j]));
  if(j > 0) flag |= (e == flip(flip1, A[i][j-1]));
  if(j < w-1) flag |= (e == flip(flip1, A[i][j+1]));
  return flag;
}

bool check_row(const bool flip1, const bool flip2, const bool flip3, const int i) {
  bool flag = true;
  rep(e, w) flag &= check_elem(flip1, flip2, flip3, i, e);
  return flag;
}

void solve() {
  cin >> h >> w;

  A.resize(h);
  rep(i, h) A[i].resize(w);
  rep(i, h) loadVec(w, A[i]);
  
  VVVI dp(h, VVI(2, VI(2, INT_MAX)));
  dp[0][0][0] = 0;
  dp[0][0][1] = 1;
  
  rep1(i, h-1) {
    rep(j, 2) { // i-1 flip or not
      rep(k, 2) { // i flip or not
        if(check_row(j, k, 0, i-1)) {
          if(i != h-1 || check_row(k, 0, j, i)) dp[i][j][k]= min(dp[i][j][k], dp[i-1][0][j]) + k;
        }
        if(check_row(j, k, 1, i-1)) {
          if(i != h-1 || check_row(k, 0, j, i)) dp[i][j][k]= min(dp[i][j][k], dp[i-1][1][j]) + k;
        }
      }
    }
  }
  int ans = INT_MAX;
  if(check_row(0, 0, 0, h-1)) ans = min(ans, dp[h-1][0][0]);
  if(check_row(0, 0, 1, h-1)) ans = min(ans, dp[h-1][1][0]);
  if(check_row(1, 0, 0, h-1)) ans = min(ans, dp[h-1][0][1]);
  if(check_row(1, 0, 1, h-1)) ans = min(ans, dp[h-1][1][1]);
  if(ans == INT_MAX) ans = -1;
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