#include <bits/stdc++.h>
using namespace std;
#include <atcoder/all>
using namespace atcoder;
#define rep(i, n) for (int i = 0; i < (int)(n); i++)
#define rep1(i, n) for (int i = 1; i < (int)(n+1); i++)
using ll = long long;
using p = pair<int, int>;
using pf = pair<double, double>;

void solve() {
  int n;
  cin >> n;
  vector<int> x(5 * n);
  for(auto&& a: x) cin >> a;

  sort(x.begin(), x.end());
  rep(i, n) x.erase(x.begin());
  rep(i, n) x.erase(x.end() - 1);

  ll score = 0;
  for(auto&& a: x) score += a;

  double ans = (double)score / x.size();
  cout << ans << endl;
}

int main() {
#ifdef ONLINE_JUDGE 
  const int n_testcase = 1;  // Don't change here!!
#else
  const int n_testcase = 2;  // Change here for test cases.
#endif

  rep(i, n_testcase) {
    solve();
  }

  return 0;
}