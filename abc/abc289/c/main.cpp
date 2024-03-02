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
using VI = vector<int>;

int n, m;
vector<VI> s;

int check(const unordered_set<int> before_set, int i) {
  unordered_set<int> new_set = before_set;
  for(auto&& elem: s[i]) {
    new_set.insert(elem);  
  }

  if(i == m - 1) {
    int total = 0;
    total += new_set.size() == n ? 1 : 0;
    total += before_set.size() == n ? 1 : 0;
    return total;
  }

  int total = 0;
  total += check(new_set, i+1);
  total += check(before_set, i+1);
  return total;
}

void solve() {
  cin >> n >> m;

  s.clear();
  s.resize(m);

  rep(i, m) {
    int ci;
    cin >> ci;
    s[i].resize(ci);
    rep(j, ci) cin >> s[i][j];
  }

  unordered_set<int> current_set;
  int ans = check(current_set, 0);

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