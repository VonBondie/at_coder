#include <bits/stdc++.h>
using namespace std;
#include <atcoder/all>
using namespace atcoder;
#define rep(i, n) for (int i = 0; i < (int)(n); i++)
#define rep1(i, n) for (int i = 1; i < (int)(n+1); i++)
using ll = long long;
using p = pair<int, int>;
using pf = pair<double, double>;

int main() {
  std::string s;
  cin >> s;

  ll num = 0;
  for(auto&& c : s) {
    num *= 26;
    num += c - 'A' + 1;
  }

  ll ans = num;

  cout << ans << endl;
  return 0;
}