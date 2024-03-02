#include <bits/stdc++.h>
using namespace std;
#include <atcoder/all>
using namespace atcoder;

// macro
#define rep(i, n) for (int i = 0; i < (int)(n); i++)
#define rep1(i, n) for (int i = 1; i < (int)(n+1); i++)
#define PRINT_YN(b) cout << b ? "Yes" : "No" << endl;

// type
using LL = long long;
using PII = pair<int, int>;
using PFF = pair<double, double>;

void solve() {
  int n, m;
  cin >> n >> m;

  vector<int> a(n), b(m);
  rep(i, n) cin >> a[i];
  rep(i, m) cin >> b[i];

  LL total = 0;
  for(auto&& idx: b) {
    total += a[idx - 1];
  }

  LL ans = total;

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
  }

  return 0;
}