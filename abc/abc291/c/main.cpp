#include <bits/stdc++.h>
using namespace std;
#include <atcoder/all>
using namespace atcoder;
#define rep(i, n) for (int i = 0; i < (int)(n); i++)
#define rep1(i, n) for (int i = 1; i < (int)(n+1); i++)
using ll = long long;
using p = pair<int, int>;
using pf = pair<double, double>;

string toString(p data) {
  string tmp;
  tmp += to_string(data.first);
  tmp += ",";
  tmp += to_string(data.second);
  return tmp;
}

void solve() {
  int n;
  cin >> n;
  string move;
  cin >> move;

  p cod = p(0, 0);
  unordered_map<string, int> cordinates;
  cordinates[toString(cod)] = 0;

  rep1(i, n) {
    switch(move[i-1]) {
      case 'R' : cod.first += 1; break;
      case 'L' : cod.first -= 1; break;
      case 'U' : cod.second += 1; break;
      case 'D' : cod.second -= 1; break;
    }
    if(cordinates.find(toString(cod)) != cordinates.end()) {
      cout << "Yes" << endl;
      return;
    }

    cordinates[toString(cod)] = i;
  }
  
  cout << "No" << endl;
  return;
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