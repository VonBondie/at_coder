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

void solve() {
  int N;
  loadVar(N);
  VS S(N);
  loadVec(N,S);

  int y1 = -1, x1 = -1;
  int y2 = -1, x2 = -1;
  
  rep(i, N) {
    rep(j, N) {
      if(S[i][j] != 'P') continue;;
      if(x1 == -1) y1 = i, x1 = j;
      else y2 = i, x2 = j;
    }
  }

  // vector<vector<vector<vector<int>>>> T(N);
  vector T(N, vector(N, vector(N, vector(N, INTMAX_MAX))));

  T[y1][x1][y2][x2] = 0;
  queue<array<int, 4>> que;
  que.push({y1, x1, y2, x2});

  int dy[] = {-1,  0, 1, 0};
  int dx[] = { 0, -1, 0, 1};

  while(!que.empty()) {
    array<int, 4> current = que.front();
    que.pop();
    rep(i, 4) {
      array<int, 4> next = current;
      next[0] += dy[i]; next[1] += dx[i];
      if(next[0] < 0 || next[0] >= N || next[1] < 0 || next[1] >= N || S[next[0]][next[1]] == '#') {
        next[0] = current[0]; next[1] = current[1];
      }

      next[2] += dy[i]; next[3] += dx[i];
      if(next[2] < 0 || next[2] >= N || next[3] < 0 || next[3] >= N || S[next[2]][next[3]] == '#') {
        next[2] = current[2]; next[3] = current[3];
      }

      if(T[next[0]][next[1]][next[2]][next[3]] == INTMAX_MAX) {
        T[next[0]][next[1]][next[2]][next[3]] = T[current[0]][current[1]][current[2]][current[3]] + 1;
        que.push(next);
      }
    }
  }

  long ans = INTMAX_MAX;
  rep(i, N) {
    rep(j, N) {
     ans = min(ans, T[i][j][i][j]);
    }
  }
  cout << (ans == INTMAX_MAX ? -1 : ans) << endl;

}

int main(int argc, char *argv[]) {
#ifdef ONLINE_JUDGE 
  const int n_testcase = 1;  // Don't change here!!
#else
  ifstream ifs("case_num.txt");
  int tmp = 1;
  ifs >> tmp;
  const int n_testcase = tmp;  // Don't change here!!
#endif

  rep(i, n_testcase) {
    solve();
  }

  return 0;
}
