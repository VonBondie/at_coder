#include <bits/stdc++.h>
using namespace std;
#include <atcoder/all>
using namespace atcoder;
#define rep(i, n) for (int i = 0; i < (int)(n); i++)
#define rep1(i, n) for (int i = 1; i < (int)(n+1); i++)
using LL = long long;
using PII = pair<int, int>;
using PFF = pair<double, double>;

unordered_map<string, pair<string, bool>> users;

bool check(pair<string, bool>& value, const string target) {
  if(value.second) return false;
  value.second = true;

  if(users.find(value.first) != users.end()) {
    auto next_value = users[value.first];
    if(target == next_value.first) return true;
    if(check(next_value, target)) return true;
  }

  return false;
}

void solve() {
  int n;
  cin >> n;
  users.clear();
  rep(i, n) {
    string before, after;
    cin >> before >> after;
    users[before] = pair(after, false);
  }

  bool flag = false;

  for(auto&& user: users) {
    flag |= check(user.second, user.first);
    user.second.second = true;
    if(flag) break;
  }

  string ans = !flag ? "Yes" : "No";
  cout << ans << endl;
}

int main() {
#ifdef ONLINE_JUDGE 
  int n_testcase = 1;
#else
  int n_testcase = 3;
#endif

  while(n_testcase--) solve();

  return 0;
}