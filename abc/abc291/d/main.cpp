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
  const ll div = 998244353;
  int n;
  cin >> n;

  vector<pair<ll, ll>> cards(n);
  rep(i, n) cin >> cards[i].first >> cards[i].second;

  vector<ll> na(n, 0), nb(n, 0);

  na[0] = 1;
  nb[0] = 1;

  rep1(i, n-1) {
    auto card = cards[i];
    auto before = cards[i-1];
    
    if(cards[i].first != before.first) na[i] += na[i-1] % div;
    if(cards[i].first != before.second) na[i] += nb[i-1] % div;
    if(cards[i].second != before.first) nb[i] += na[i-1] % div;
    if(cards[i].second != before.second) nb[i] += nb[i-1] % div;
  }

  ll total = na[n-1] + nb[n-1];

  ll ans = total % div;

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