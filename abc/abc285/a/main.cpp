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
  const int max = 16;
  int a, b;
  cin >> a >> b;

  bool flag = false;
  if(b == a * 2 || b == a * 2 + 1) flag = true;

  std::string ans = flag ? "Yes" : "No";
  cout << ans << endl;
  return 0;
}