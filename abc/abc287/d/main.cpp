#include <bits/stdc++.h>
using namespace std;
#include <atcoder/all>
using namespace atcoder;

// macro
#define rep(i, n) for (int i = 0; i < (int)(n); i++)
#define rep1(i, n) for (int i = 1; i < (int)(n+1); i++)
#define PRINT_YN(b) cout << (b ? "Yes" : "No") << endl;

// type
using LL = long long;
using PII = pair<int, int>;
using PFF = pair<double, double>;
using VI = vector<int>;
using VS = vector<string>;

struct HashPair {

    //注意 constがいる
    template<class T1, class T2>
    size_t operator()(const pair<T1, T2> &p) const {

        //first分をハッシュ化する
        auto hash1 = hash<T1>{}(p.first);

        //second分をハッシュ化する
        auto hash2 = hash<T2>{}(p.second);

        //重複しないようにハッシュ処理
        size_t seed = 0;
        seed ^= hash1 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= hash2 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};

// func
template<class T, class K>
bool contains(T collection, K elem) {
  return collection.find(elem) != collection.end();
}

void loadKernel(int i) { (void)i; }

template<class Head, class... Tail>
void loadKernel(int i, Head&& head, Tail&&... tail) {
  cin >> head[i];
  loadKernel(i, std::forward<Tail>(tail)...);
}

template<class... Dst>
void loadVec(const int n, Dst&&... dst) {
  rep(i, n) loadKernel(i, std::forward<Dst>(dst)...);
}

string notationKernel() {
  return "";
}

template<class Head, class... Tail>
string notationKernel(Head&& head, Tail&&... tail) {
  string s = "," + to_string(head);
  return s + notation(std::forward<Tail>(tail)...);
}

template<class Head, class... Tail>
string notation(Head&& head, Tail&&... tail) {
  string s = to_string(head);
  return s + notation(std::forward<Tail>(tail)...);
}
////////////////////////////////////////////////////

bool isAllTrue(const vector<bool>& vec) {
  return all_of(vec.begin(), vec.end(), [](bool x){ return x; });
}

void solve() {
  string s,t;
  cin >> s >> t;

  const int lt = t.size();
  vector<bool> check(lt);
  string sd = s.substr(s.size() - lt, lt);

  rep(i, lt) check[i] = (sd[i] == t[i] || sd[i] == '?' || t[i] == '?');
  PRINT_YN(isAllTrue(check));

  int total_false = 0;
  rep(i, lt) if(!check[i]) total_false++;

  rep1(i, lt) {
    const int size_back = lt - i;
    sd[i-1] = s[i-1];
    const bool current_change = (sd[i-1] == t[i-1] || sd[i-1] == '?' || t[i-1] == '?');
    if(!check[i - 1] && current_change) total_false--;
    if(check[i - 1] && !current_change) total_false++;

    PRINT_YN(total_false == 0);
  }
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