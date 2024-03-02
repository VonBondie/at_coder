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

bool check(const int H, const int W,
           const int h, const int w,
           const char c, const int cost) {
  if(h >= H) return false;
  if(h < 0) return false;
  if(w >= W) return false;
  if(w < 0) return false;
  if(c == '#') return false;
  if(cost != INT_MAX) return false;

  return true;
}

VI dykstra(const pair<int, int> s, const VS& A, const int n_node) {
  VI ans(n_node, INT_MAX);
  const int H = A.size();
  const int W = A[0].size();
  
  VVI field(H, VI(W, INT_MAX));
  queue<pair<int,int>> Q;

  field[s.first][s.second] = 0;
  Q.push(s);
  while(!Q.empty()) {
    pair<int, int> current = Q.front();
    Q.pop();
    int h = current.first;
    int w = current.second;

    queue<pair<int, int>> tmp_q;
    tmp_q.push(pair<int,int>(h+1, w));
    tmp_q.push(pair<int,int>(h-1, w));
    tmp_q.push(pair<int,int>(h, w+1));
    tmp_q.push(pair<int,int>(h, w-1));

    while(!tmp_q.empty()){
      auto next = tmp_q.front();
      tmp_q.pop();
      int h_ = next.first;
      int w_ = next.second;
      if(check(H, W, h_, w_, A[h_][w_], field[h_][w_])) {
        field[h_][w_] = field[h][w] + 1;
        if(A[h_][w_] == 'S') ans[0] = field[h_][w_];
        if(A[h_][w_] == 'G') ans[n_node -1] = field[h_][w_];
        if(A[h_][w_] < 0) ans[-A[h_][w_]] = field[h_][w_];
        Q.push(next);
      }
    }
  }
  return ans;
}

VVI DP(const VVI& graph, const int capacity) {
  VVI dp(pow(2,graph.size()), VI(graph.size() + 1, INT_MAX));

  dp[1][1] = 0;
  int i = 1;
  queue<pair<int, int>> Q;
  Q.push(pair<int, int>(1, 1));

  while(!Q.empty()) {
    pair<int, int> current = Q.front();
    Q.pop();
    int status = current.first;
    int pos = current.second;

    rep1(i, graph.size() -1) {
      if((status >> i) & 1) continue;
      if(graph[pos][i] == INT_MAX) continue;
      dp[status | (1 << i)][i] = graph[pos][i] + dp[status][pos];
      Q.push(pair<int, int>(status | (1 << i), i));
    }
  }

  return dp;
}

void solve() {
  int H, W, T; loadVar(H, W, T);
  VS A(H); loadVec(H, A);

  vector<pair<int, int>> foods;
  pair<int, int> start;
  pair<int, int> goal;

  int food_id = 1;
  rep(i, H) {
    rep(j, W) {
      if(A[i][j] == 'o') {
        foods.emplace_back(i, j);
        A[i][j] = -food_id;
        food_id++;
      }
      if(A[i][j] == 'S') start = pair<int, int>(i, j);
      if(A[i][j] == 'G') goal = pair<int, int>(i, j);
    }
  }

  const int n_node = 2+foods.size();
  VVI graph(n_node, VI(n_node));

  graph[0] = dykstra(start, A, n_node);
  rep(i, foods.size()) {
    graph[i+1] = dykstra(foods[i], A, n_node);
  }

  auto dp = DP(graph, T);
  int ans = -1;
  rep(i, pow(2, n_node)) {
    if(dp[i][n_node -1] != INT_MAX) max(ans, __builtin_popcount(i) - 2);
  }

  cout << ans << endl;
}

int main(int argc, char *argv[]) {
#ifdef ONLINE_JUDGE 
  const int n_testcase = 1;  // Don't change here!!
  const bool flag_endl = false;
#else
  int tmp = 1;
  bool flag_endl = false;
  if(argc >= 2) {
    ifstream ifs(argv[1]);
    ifs >> tmp;
    flag_endl = true;
  }
  const int n_testcase = tmp;  // Don't change here!!
#endif

  rep(i, n_testcase) {
    solve();
    if(flag_endl)  cout << endl;
  }

  return 0;
}