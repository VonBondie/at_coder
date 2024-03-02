#include <bits/stdc++.h>
using namespace std;
#include <atcoder/all>
using namespace atcoder;

// macro
#define rep(i, n) for (int i = 0; i < (int)(n); i++)
#define rep1(i, n) for (int i = 1; i < (int)(n+1); i++)
#define PRINT_YN(b) cout << (b ? "Yes" : "No") << endl;

// type
using LL = long long;
using PII = pair<int, int>;
using PFF = pair<double, double>;
using VI = vector<int>;

// func
template<class T, class K>
bool contains(T collection, K elem) {
  return collection.find(elem) != collection.end();
}

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

void solve() {
  int n, m, x;
  cin >> n;
  VI a(n);
  loadVec(n, a);

  cin >> m;
  unordered_set<int> b;
  rep(i, m) {
    int tmp;
    cin >> tmp;
    b.insert(tmp);
  }
  cin >> x;

  sort(a.begin(), a.end());
  vector<bool> steps(x+1, false);
  steps[0] = true;
  rep(i, x - 1) {
    if(!steps[i]) continue;
    for(auto&& ax: a) {
      if (i+ax <= x && !contains(b, i+ax)) steps[i+ax] = true; 
    }
  }

  PRINT_YN(steps[x]);
}

int main() {
#ifdef ONLINE_JUDGE 
  const int n_testcase = 1;  // Don't change here!!
#else
  const int n_testcase = 3;  // Change here for test cases.
#endif

  rep(i, n_testcase) solve();
  return 0;
}