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
  if (v.size() == 0)  return;
  cout << v[0];
  rep1(i, v.size()-1) cout << " " << v[i];
  cout <<endl;
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

void solve() {
  int N; loadVar(N);
  VI A(N), B(N); loadVec(N, A, B);

  set<int> nums;
  rep(i, N) {
    int a_ = A[i] - 1, b_ = B[i] - 1;
    nums.insert(a_); nums.insert(b_);
  }

  map<int, int> id_table;
  map<int, int> id_table_;
  id_table[0] = 0;
  id_table_[0] = 0;
  int i = 1;
  for(auto e:nums) {
    if(e == 0) continue;
    id_table[i] = e; 
    id_table_[e] = i;
    i++;
  }

  int tmp = nums.size() + (*nums.begin() != 0 ? 1 : 0);
  dsu U(tmp);

  rep(i, N) {
    int a_ = A[i] - 1, b_ = B[i] - 1;
    U.merge(id_table_[a_], id_table_[b_]);
  }

  int ans = 1;
  for(int i = tmp-1; i> 0; i--) {
    if(U.same(0, i)) {
      ans = id_table[i] + 1;
      break;
    }
  }
  // map<int,list<int>> branchs;
  // rep(i, N) {
  //   int a_ = A[i] - 1, b_ = B[i] - 1;
  //   branchs[a_].push_back(b_);
  //   branchs[b_].push_back(a_);
  // }

  // unordered_set<int> used;
  // int ans = 0;
  // queue<int> q;
  // q.push(0);
  // while(!q.empty()) {
  //   int tmp = q.front(); q.pop();
  //   used.insert(tmp);
  //   ans = max(tmp, ans);
  //   for(auto&& e: branchs[tmp]) {
  //     if(!used.count(e)) q.push(e);
  //   }
  // }

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