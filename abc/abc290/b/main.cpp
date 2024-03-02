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
  int n, k;
  cin >> n >> k;

  string s;
  cin >> s;

  string t;
  for(auto&& c: s) {
    if(k > 0 && c == 'o') {
      t += 'o';
      k--;
    } else {
      t += 'x';
    }
  }

  cout << t << endl;
}

int main() {
#ifdef ONLINE_JUDGE 
  const int n_testcase = 1;  // Don't change here!!
#else
  const int n_testcase = 1;  // Change here for test cases.
#endif

  rep(i, n_testcase) {
    solve();
  }

  return 0;
}