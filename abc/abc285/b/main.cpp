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
  int n;
  string s;
  cin >> n >> s;

  std::vector<int> ans;
  rep1(i, n-1) {
    int j = 0;
    bool flag = false;
    for(j = 0; j + i < n; j++) {
      if(s[j] == s[j+i]) {
        flag = false; 
        break;
      }
    }

    ans.push_back(j + 1 + ((flag) ? 0 : -1));
  }

  for(auto&& a : ans) {
    cout << a << endl;
  }
  return 0;
}