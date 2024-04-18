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

template <class T>
vector<T> setToVector(const set<T>& s) {
  vector<T> vec(s.begin(), s.end());
  return vec;
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
int H, W;
VVI A;
vector<array<int, 3>> Data;

list<int> getReachable(int i) {
  VVI T(H, VI(W));

  list<int> reachables;
  queue<array<int, 3>> q;
  q.push(Data[i]);
  T[Data[i][0]][Data[i][1]] = Data[i][2];
  if(Data[i][2] == 0) return reachables;

  int dy[] = {-1, 0, 1, 0};
  int dx[] = {0, -1, 0, 1};
  while(!q.empty()) {
    array<int, 3> current = q.front();
    q.pop();
    rep(i, 4) {
      int ny = current[0] + dy[i];
      int nx = current[1] + dx[i];
      int ne = current[2] - 1;
      
      if(ny >= 0 && nx >= 0 && ny < H && nx < W && A[ny][nx] != -1) {
        if(T[ny][nx] != 0) continue;
        T[ny][nx] = ne;
        if(A[ny][nx] != -2) {
          reachables.push_back(A[ny][nx]);
        }
        if(ne > 0) {
          q.push(array<int, 3>{ny, nx, ne});
        }
      }
    }
  }

  return reachables;
}

void solve() {
  loadVar(H, W);
  A.resize(H);
  int sx, sy, tx, ty;
  rep(i, H) {
    A[i].resize(W);
    rep(j, W) {
      char c;
      cin >> c;
      switch(c) {
        case 'S':
          sy = i; sx = j;
          break;
        case 'T':
          ty = i; tx = j;
          break;
        case '#':
          A[i][j] = -1;
          break;
        case '.':
          A[i][j] = -2;
      }
    }
  }
  int N;
  loadVar(N);
  Data.reserve(N+2);
  Data.push_back(array<int, 3>{ty, tx, 0});
  Data.push_back(array<int, 3>{sy, sx, 0});
  rep(i, N) {
    int r, c, e;
    cin >> r >> c >> e;
    r--;c--;
    if(r == sy && c == sx) {
      Data[1][2] = e;
      continue;
    }
    Data.push_back(array<int, 3>{r, c, e});
  }
  rep(i, Data.size()) {
    int y = Data[i][0];
    int x = Data[i][1];
    A[y][x] = i;
  }

  queue<int> q;
  q.push(1);
  dsu T(Data.size());
  while(!q.empty()) {
    int c = q.front();
    q.pop();
    for(const auto& e: getReachable(c)){
      if(T.same(c, e)) continue;
      T.merge(c, e);
      q.push(e);
    }
  }
  PRINT_YN(T.same(0, 1));
}

int main(int argc, char *argv[]) {
#ifdef ONLINE_JUDGE 
  const int n_testcase = 1;  // Don't change here!!
#else
  ifstream ifs(argv[1]);
  int tmp = 1;
  ifs >> tmp;
  const int n_testcase = tmp;  // Don't change here!!
#endif

  rep(i, n_testcase) {
    solve();
    Data.clear();
  }

  return 0;
}
