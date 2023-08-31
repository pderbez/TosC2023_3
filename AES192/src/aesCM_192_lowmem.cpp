#include <iostream>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <map>

#include "SysOfEqs.hpp"

using namespace std;

vector<vector<unsigned>> initMC16() {
  vector<vector<unsigned>> MC16 (5);
  MC16[0] = vector<unsigned> ({0});
  MC16[1] = vector<unsigned> ({4});
  MC16[2] = vector<unsigned> ({3, 4});
  MC16[3] = vector<unsigned> ({2, 3, 4});
  MC16[4] = vector<unsigned> ({1, 2, 3, 4});
  return MC16;
}

vector<uint8_t> initPop5() {
  vector<uint8_t> count;
  for (unsigned s = 0; s < 5*5*5*5; ++s) {
    unsigned pow5 = 1;
    unsigned a = 0;
    for (unsigned c = 0; c < 4; ++c) {
      a += (s/pow5)%5;
      pow5 *= 5;
    }
    count.emplace_back(a);
  }
  return count;
}

void updateSR(vector<uint8_t> & T, uint8_t const global_bound) {
  static unsigned const n_states = 5*5*5*5;
  static unsigned const n_keys = 5*5*5*5*5*5;
  static unsigned const mypow[8] = {1, 5, 5*5, 5*5*5, 5*5*5*5, 5*5*5*5*5, 5*5*5*5*5*5, 5*5*5*5*5*5*5};
  static vector<uint8_t> const count = initPop5();
  static vector<unsigned> const shiftRows ({0, 1, 2, 3, 7, 4, 5, 6, 10, 11, 8, 9, 13, 14, 15, 12});

  vector<uint8_t> TT (n_states*n_keys, global_bound);

  for (unsigned state = 0; state < n_states; ++state) {

    unsigned mystate[4];
    for (unsigned c = 0; c < 4; ++c) mystate[c] = (state/mypow[c])%5;
    sort(&mystate[0], &mystate[0] + 4);
    unsigned smallest_state = 0;
    for (unsigned c = 0; c < 4; ++c) {smallest_state *= 5; smallest_state += mystate[c];}

    if (state != smallest_state) {
      #pragma omp parallel for
      for (unsigned key = 0; key < n_keys; ++key) {
        TT[key*n_states + state] = TT[key*n_states + smallest_state];
      }
    }
    else {
      set<unsigned> all_states;
      all_states.emplace(0);
      for (unsigned c = 0; c < 4; ++c) {
        auto const & u = mystate[c];
        set<unsigned> tmp;
        for (unsigned x = 0; x < 16; ++x) {
          if (__builtin_popcount(x) != u) continue;
          for (auto y : all_states) {
            y += ((x >> 0) & 1)*mypow[c];
            y += ((x >> 1) & 1)*mypow[(c+1)%4];
            y += ((x >> 2) & 1)*mypow[(c+2)%4];
            y += ((x >> 3) & 1)*mypow[(c+3)%4];
            tmp.emplace(y);
          }
        }
        swap(tmp, all_states);
      }

      #pragma omp parallel for
      for (unsigned key = 0; key < n_keys; ++key) {
        auto & dst = TT[key*n_states + state];
        for (auto const & s : all_states) {
          auto const & src = T[key*n_states + s];
          dst = min(dst, src);
        }
      }
    }
  }
  swap(T, TT);
}

vector<unsigned> inv_updateSR(vector<uint8_t> const & T, uint8_t const bound, unsigned state, unsigned key, vector<uint8_t> const & valX) {
  static unsigned const n_states = 5*5*5*5;
  static unsigned const mypow[4] = {1, 5, 5*5, 5*5*5};

  vector<unsigned> res;

  {
    set<unsigned> all_states;
    all_states.emplace(0);
    for (unsigned c = 0; c < 4; ++c) {
      auto u = (state/mypow[c])%5;
      set<unsigned> tmp;
      for (unsigned x = 0; x < 16; ++x) {
        if (__builtin_popcount(x) != u) continue;
        if (valX[4*0 + ((c+0)%4)] != 2 && valX[4*0 + ((c+0)%4)] != ((x >> 0) & 1)) continue;
        if (valX[4*1 + ((c+1)%4)] != 2 && valX[4*1 + ((c+1)%4)] != ((x >> 1) & 1)) continue;
        if (valX[4*2 + ((c+2)%4)] != 2 && valX[4*2 + ((c+2)%4)] != ((x >> 2) & 1)) continue;
        if (valX[4*3 + ((c+3)%4)] != 2 && valX[4*3 + ((c+3)%4)] != ((x >> 3) & 1)) continue;
        for (auto y : all_states) {
          y += ((x >> 0) & 1)*mypow[c];
          y += ((x >> 1) & 1)*mypow[(c+1)%4];
          y += ((x >> 2) & 1)*mypow[(c+2)%4];
          y += ((x >> 3) & 1)*mypow[(c+3)%4];
          tmp.emplace(y);
        }
      }
      swap(tmp, all_states);
    }

    for (auto s : all_states) {
      auto const & src = T[key*n_states + s];
      if (src <= bound) res.emplace_back(key*n_states + s);
    }
  }
  return res;
}

void updateMC_ARK(vector<uint8_t> & T, unsigned const col, unsigned const colk, uint8_t const global_bound) {
  static unsigned const n_states = 5*5*5*5;
  static unsigned const n_keys = 5*5*5*5*5*5;
  static unsigned const mypow[8] = {1, 5, 5*5, 5*5*5, 5*5*5*5, 5*5*5*5*5, 5*5*5*5*5*5, 5*5*5*5*5*5*5};
  static const auto MC16 = initMC16();

  vector<uint8_t> TT (n_states*n_keys, global_bound);

  #pragma omp parallel for
  for (unsigned key_tmp = 0; key_tmp < n_keys/5; ++key_tmp) {
    unsigned const & kk = (key_tmp/mypow[colk])%5;
    for (int k = 0; k <= 4; ++k) {
      unsigned const & key = key_tmp + k*mypow[colk] + kk*(mypow[5] - mypow[colk]);
      int x = (k == 0) ? 1 : 0;
      for (auto z : MC16[x+k]) {
        for (unsigned state_tmp = 0; state_tmp < n_states/5; ++state_tmp) {
          unsigned const & ss = (state_tmp/mypow[col])%5;
          unsigned const & state_src = state_tmp + z*mypow[col] + ss*(mypow[3] - mypow[col]);
          unsigned const & state_dst = state_tmp + x*mypow[col] + ss*(mypow[3] - mypow[col]);
          uint8_t const & src = T[key*n_states + state_src] + x;
          if (src >= global_bound) continue;
          auto & dst = TT[key*n_states + state_dst];
          dst = min(dst, src);
        }
      }
      int z = 4 - (x+k);
      for (++x; x <= 4; ++x) {
        for (unsigned state_tmp = 0; state_tmp < n_states/5; ++state_tmp) {
          unsigned const & ss = (state_tmp/mypow[col])%5;
          unsigned const & state_dst = state_tmp + x*mypow[col] + ss*(mypow[3] - mypow[col]);
          unsigned const & state_src1 = state_tmp + (x-1)*mypow[col] + ss*(mypow[3] - mypow[col]);
          if (z == 0) TT[key*n_states + state_dst] = TT[key*n_states + state_src1] + 1;
          else {
            unsigned const & state_src2 = state_tmp + z*mypow[col] + ss*(mypow[3] - mypow[col]);
            TT[key*n_states + state_dst] = min(TT[key*n_states + state_src1] + 1, T[key*n_states + state_src2] + x);
          }
        }
        if (z > 0) --z;
      }
      for (unsigned state_tmp = 0; state_tmp < n_states/5; ++state_tmp) {
        unsigned const & ss = (state_tmp/mypow[col])%5;
        unsigned const & state_src = state_tmp + ss*(mypow[3] - mypow[col]);
        unsigned const & state_dst = state_tmp + k*mypow[col] + ss*(mypow[3] - mypow[col]);
        uint8_t const & src = T[key*n_states + state_src] + k;
        if (src >= global_bound) continue;
        auto & dst = TT[key*n_states + state_dst];
        dst = min(dst, src);
      }
    }
  }
  swap(T, TT);
}

vector<unsigned> inv_updateMC_ARK(vector<uint8_t> const & T, unsigned const col_start, uint8_t const bound, unsigned const state, unsigned const key) {
  static unsigned const n_states = 5*5*5*5;
  static unsigned const mypow[8] = {1, 5, 5*5, 5*5*5, 5*5*5*5, 5*5*5*5*5, 5*5*5*5*5*5, 5*5*5*5*5*5*5};
  static const auto MC16 = initMC16();

  vector<unsigned> res(1, state);

  for (unsigned col = 0; col < 4; ++col) {
    unsigned const & colk = (col_start + col)%6;
    unsigned const & k = (key/mypow[colk])%5;
    vector<unsigned> tmp;
    for (auto s : res) {
      unsigned const & x = (s/mypow[col])%5;

      auto const & pos = s - x*mypow[col];

      {
        for (int z : MC16[min(x+k, 4u)]) {
          auto const & full = pos + z*mypow[col];
          tmp.emplace_back(full);
        }

        if (x == k && x != 0) tmp.emplace_back(pos);
      }
    }
    swap(res, tmp);
  }

  unsigned i = 0, n = res.size();
  while (i < n) {
    res[i] = key*n_states + res[i];
    if (T[res[i]] <= bound) ++i;
    else res[i] = res[--n];
  }

  res.resize(n);


  return res;
}

void updateKey192Column(unsigned const col, vector<uint8_t> & T, uint8_t const global_bound) {
  static unsigned const n_states = 5*5*5*5;
  static unsigned const n_keys = 5*5*5*5*5*5;
  static unsigned const mypow[8] = {1, 5, 5*5, 5*5*5, 5*5*5*5, 5*5*5*5*5, 5*5*5*5*5*5, 5*5*5*5*5*5*5};
  //static vector<uint8_t> const count = initPop5();
  //static const auto MC16 = initMC16();

  vector<uint8_t> TT (n_states*n_keys, global_bound);

  unsigned col2 = (col == 0) ? 5 : col-1;

  vector<vector<unsigned>> possibleCol (5*5);
  for (int k0 = 0; k0 <= 4; ++k0) {
    for (int k1 = 0; k1 <= 4; ++k1) {
      for (int k = abs(k0 - k1); k <= min(4, k0+k1); ++k) {
        possibleCol[k0 + 5*k1].emplace_back(k);
      }
    }
  }

  #pragma omp parallel for
  for (unsigned key_tmp = 0; key_tmp < n_keys/5; ++key_tmp) {
    unsigned const & n7 = (key_tmp/mypow[col])%5;
    vector<unsigned> all_keys;
    all_keys.reserve(5);
    for (unsigned n0 = 0; n0 <= 4; ++n0) {
      unsigned key = key_tmp + n0*mypow[col] + n7*(mypow[5] - mypow[col]);
      unsigned const & n1 = (key/mypow[col2])%5;
      unsigned const & key2 = key - n0*mypow[col];
      for (auto x : possibleCol[n0 + 5*n1]) all_keys.emplace_back((key2 + x*mypow[col])*n_states);
      //if (all_keys.empty()) continue;
      key *= n_states;

      for (unsigned state = 0; state < n_states; ++state) {
        auto const & src = T[key + state];
        if (src >= global_bound) continue;
        //cout << all_keys.size() << " - " << n0 << " - " << n1 << endl;
        for (auto k : all_keys) {
          auto & dst = TT[k + state];
          dst = min(dst, src);
        }
      }
      all_keys.clear();
    }
  }

  swap(T,TT);
}

vector<unsigned> inv_updateKey192Column(unsigned const col, unsigned key, vector<unsigned> & res) {
  static unsigned const mypow[8] = {1, 5, 5*5, 5*5*5, 5*5*5*5, 5*5*5*5*5, 5*5*5*5*5*5, 5*5*5*5*5*5*5};
  //static vector<uint8_t> const count = initPop5();
  //static const auto MC16 = initMC16();

  {
    unsigned const & col2 = (col == 0) ? 5 : col-1;
    unsigned const & n0 = (key/mypow[col])%5;
    unsigned const & n1 = (key/mypow[col2])%5;
    key -= n0*mypow[col];

    for (unsigned k = (n0 >= n1) ? n0 - n1 : n1 - n0; k <= min(n0+n1, 4u); ++k) {
      unsigned const & full = key + k*mypow[col];
      res.emplace_back(full);
    }
  }
  return res;
}

vector<unsigned> inv_updateKey192(int col_start, unsigned n_c, vector<uint8_t> const & T, uint8_t const bound, unsigned state, unsigned key) {
  static unsigned const n_states = 5*5*5*5;

  vector<unsigned> res (1, key);

  for (unsigned c = 0; c < n_c; ++c) {
    vector<unsigned> tmp;
    for (auto k : res) inv_updateKey192Column((col_start-c + 6)%6, k, tmp);
    swap(res, tmp);
  }

  unsigned i = 0, n = res.size();
  while (i < n) {
    res[i] = res[i]*n_states + state;
    if (T[res[i]] <= bound) ++i;
    else res[i] = res[--n];
  }

  res.resize(n);

  return res;
}




unsigned searchOnLine(unsigned l, uint8_t x, vector<vector<uint8_t>> const & valX, vector<vector<uint8_t>> const & valK, Matrix const & mat) {
  for (unsigned c = 0; c < mat.nbcols; ++c) {
    if (mat(l, c) != 0) {
      int uu = abs(mat.getColumns(c));
      int ii = uu%16;
      int vv = (uu/16)%4;
      int rr = (uu/16)/4;
      if (vv == 2 && valX[rr][ii] == x) return c;
      else if (vv == 1 && valK[rr][ii] == x) return c;
    }
  }
  return mat.nbcols;
}

void set0Mat(int uval, vector<vector<uint8_t>> & valX, vector<vector<uint8_t>> & valK, unsigned & line1, unsigned & line2, Matrix & mat) {
  if (mat.setAsPivot(uval, line1, line2)) {
    unsigned tmp = searchOnLine(line1, 2, valX, valK, mat);
    if (tmp == mat.nbcols) mat.swapLines(line1, --line2);
    else {
      mat.swapLineColumn(line1,tmp);
      mat.eraseColumn(tmp);
    }
  }
  else {
    for (unsigned ll = 0; ll < mat.nbcols; ++ll) {
      if (mat.getColumns(ll) == uval) {
        mat.eraseColumn(ll);
        break;
      }
    }
  }
}

pair<bool, bool> updateColumns_X(int r, set<int> & set_x, set<int> & set_sr, vector<vector<uint8_t>> & valX, vector<vector<uint8_t>> & valK, unsigned & line1, unsigned & line2, Matrix & mat, vector<vector<uint8_t>> const & valColX, vector<vector<uint8_t>> const & valColSR, vector<vector<uint8_t>> const & valColK) {
  pair<bool, bool> res = make_pair(false, true);
  bool done_x[4] = {false, false, false, false};
  bool done_sr[4] = {false, false, false, false};

  while (!set_x.empty() || !set_sr.empty()) {
    if (!set_x.empty()) {
      auto it = set_x.begin();
      int cc = *it;
      set_x.erase(it);
      if (valColX[r][cc] == 5) continue;
      unsigned cpt[3] = {0,0,0};
      for (unsigned ll = 0; ll < 4; ++ll) cpt[valX[r][4*ll + cc]] += 1;
      if (cpt[0] > 4-valColX[r][cc] || cpt[1] > valColX[r][cc]) {res.second = false; return res;}
      if (cpt[2] == 0) {done_x[cc] = true; continue;}
      if (cpt[0] == 4-valColX[r][cc]) {
        done_x[cc] = true;
        for (unsigned ll = 0; ll < 4; ++ll) {
          if (valX[r][4*ll + cc] != 2) continue;
          int uuval2 = 16*(4*r + 2) + 4*ll + cc;
          if (mat.setAsPivot(uuval2, line1, line2)) ++line1;
          if (mat.setAsPivot(-uuval2, line1, line2)) ++line1;
          valX[r][4*ll + cc] = 1;
          //cout << r << ": " << 4*ll + cc << " (X,1)" << endl;
          if (!done_sr[(cc - ll + 4)%4]) set_sr.emplace((cc - ll + 4)%4);
        }
      }
      else if (cpt[1] == valColX[r][cc]) {
        done_x[cc] = true;
        for (unsigned ll = 0; ll < 4; ++ll) {
          if (valX[r][4*ll + cc] != 2) continue;
          int uuval2 = 16*(4*r + 2) + 4*ll + cc;
          set0Mat(uuval2, valX, valK, line1, line2, mat);
          set0Mat(-uuval2, valX, valK, line1, line2, mat);
          valX[r][4*ll + cc] = 0;
          //cout << r << ": " << 4*ll + cc << " (X,0)" << endl;
          res.first = true;
          if (!done_sr[(cc - ll + 4)%4]) set_sr.emplace((cc - ll + 4)%4);
        }
      }
    }
    else {
      auto it = set_sr.begin();
      int cc = *it;
      set_sr.erase(it);
      if (valColSR[r][cc] == 5) continue;
      unsigned cpt[3] = {0,0,0};
      for (unsigned ll = 0; ll < 4; ++ll) cpt[valX[r][4*ll + ((cc+ll)%4)]] += 1;
      if (cpt[0] > 4-valColSR[r][cc] || cpt[1] > valColSR[r][cc]) {res.second = false; return res;}
      if (cpt[2] == 0) {done_sr[cc] = true; continue;}
      if (cpt[0] == 4-valColSR[r][cc]) {
        done_sr[cc] = true;
        for (unsigned ll = 0; ll < 4; ++ll) {
          if (valX[r][4*ll + ((cc+ll)%4)] != 2) continue;
          int uuval2 = 16*(4*r + 2) + 4*ll + ((cc+ll)%4);
          if (mat.setAsPivot(uuval2, line1, line2)) ++line1;
          if (mat.setAsPivot(-uuval2, line1, line2)) ++line1;
          valX[r][4*ll + ((cc+ll)%4)] = 1;
          //cout << r << ": " << 4*ll + ((cc+ll)%4) << " (SR,1)" << endl;
          if (!done_x[(cc+ll)%4]) set_x.emplace(((cc+ll)%4));
        }
      }
      else if (cpt[1] == valColSR[r][cc]) {
        done_sr[cc] = true;
        for (unsigned ll = 0; ll < 4; ++ll) {
          if (valX[r][4*ll + ((cc+ll)%4)] != 2) continue;
          int uuval2 = 16*(4*r + 2) + 4*ll + ((cc+ll)%4);
          set0Mat(uuval2, valX, valK, line1, line2, mat);
          set0Mat(-uuval2, valX, valK, line1, line2, mat);
          valX[r][4*ll + ((cc+ll)%4)] = 0;
          //cout << r << ": " << 4*ll + ((cc+ll)%4) << " (SR,0)" << endl;
          res.first = true;
          if (!done_x[(cc+ll)%4]) set_x.emplace(((cc+ll)%4));
        }
      }
    }
  }
  return res;
}

pair<bool, bool> updateColumns(int round_type, int r, int pos, vector<vector<uint8_t>> & valX, vector<vector<uint8_t>> & valK, unsigned & line1, unsigned & line2, Matrix & mat, vector<vector<uint8_t>> const & valColX, vector<vector<uint8_t>> const & valColSR, vector<vector<uint8_t>> const & valColK) {
  pair<bool, bool> res = make_pair(false, true);
  if (round_type == 2) {
    set<int> set_x;
    set_x.emplace(pos%4);
    set<int> set_sr;
    set_sr.emplace(((pos%4) - (pos/4) + 4)%4);
    return updateColumns_X(r, set_x, set_sr, valX, valK, line1, line2, mat, valColX, valColSR, valColK);
  }
  else {
    int cc = pos%4;
    if (valColK[r][cc] == 5) return res;
    unsigned cpt[3] = {0,0,0};
    for (unsigned ll = 0; ll < 4; ++ll) cpt[valK[r][4*ll + cc]] += 1;
    if (cpt[0] > 4-valColK[r][cc] || cpt[1] > valColK[r][cc]) {res.second = false; return res;}
    if (cpt[2] == 0) return res;
    if (cpt[0] == 4-valColK[r][cc]) {
      for (unsigned ll = 0; ll < 4; ++ll) {
        if (valK[r][4*ll + cc] != 2) continue;
        int uuval2 = 16*(4*r + 1) + 4*ll + cc;
        if (mat.setAsPivot(uuval2, line1, line2)) ++line1;
        if (mat.setAsPivot(-uuval2, line1, line2)) ++line1;
        valK[r][4*ll + cc] = 1;
        //cout << r << ": " << 4*ll + cc << " (K,1)" << endl;
      }
    }
    else if (cpt[1] == valColK[r][cc]) {
      for (unsigned ll = 0; ll < 4; ++ll) {
        if (valK[r][4*ll + cc] != 2) continue;
        int uuval2 = 16*(4*r + 1) + 4*ll + cc;
        set0Mat(uuval2, valX, valK, line1, line2, mat);
        set0Mat(-uuval2, valX, valK, line1, line2, mat);
        valK[r][4*ll + cc] = 0;
        //cout << r << ": " << 4*ll + cc << " (K,0)" << endl;
        res.first = true;
      }
    }
  }
  return res;
}

bool propagateZERO(vector<vector<uint8_t>> & valX, vector<vector<uint8_t>> & valK, unsigned & line1, unsigned & line2, Matrix & mat, vector<vector<uint8_t>> const & valColX, vector<vector<uint8_t>> const & valColSR, vector<vector<uint8_t>> const & valColK) {
  bool res = true;
  {
    //cout << "start" << flush;
    unsigned l = line1;
    while (l < line2) {
      unsigned cc = searchOnLine(l,2, valX, valK, mat);
      if (cc == mat.nbcols) {
        int uuval = mat.getFront(l);
        int uu = abs(uuval);
        int ii = uu%16;
        int vv = (uu/16)%4;
        int rr = (uu/16)/4;
        if (vv == 2) valX[rr][ii] = 0;
        else valK[rr][ii] = 0;
        mat.swapLines(l, --line2);
        unsigned ll = line1;
        while (ll < line2) {
          if (mat.getFront(ll) == -uuval) {
            cc = searchOnLine(ll, 2, valX, valK, mat);
            if (cc == mat.nbcols) mat.swapLines(ll, --line2);
            else {
              mat.swapLineColumn(ll,cc);
              mat.eraseColumn(cc);
            }
            break;
          }
          else ++ll;
        }
        if (ll == line2) {
          for (ll = 0; ll < mat.nbcols; ++ll) {
            if (mat.getColumns(ll) == -uuval) {
              mat.eraseColumn(ll);
              break;
            }
          }
        }
        if (!updateColumns(vv, rr, ii, valX, valK, line1, line2, mat, valColX, valColSR, valColK).second) return false;
        l = line1;
      }
      else ++l;
    }
  }
  return res;
}

bool propagateONE(vector<vector<uint8_t>> & valX, vector<vector<uint8_t>> & valK, unsigned line1_start, unsigned & line1, unsigned & line2, Matrix & mat, vector<vector<uint8_t>> const & valColX, vector<vector<uint8_t>> const & valColSR, vector<vector<uint8_t>> const & valColK) {
  if (line1 == line1_start) return true;
  unsigned l = line1_start;
  while (l < line1) {
    unsigned cc = 0;
    unsigned tmp = 0;
    unsigned x = 0;
    while (cc < mat.nbcols) {
      if (mat(l, cc) != 0) {
        int uu = abs(mat.getColumns(cc));
        int ii = uu%16;
        int vv = (uu/16)%4;
        int rr = (uu/16)/4;
        if (vv == 2) {
          if (valX[rr][ii] != 0) {tmp += 3 - valX[rr][ii]; x = cc;}
        }
        else {
          if (valK[rr][ii] != 0) {tmp += 3 - valK[rr][ii]; x = cc;}
        }
        if (tmp > 1) break;
      }
      ++cc;
    }
    if (tmp == 1) {
      auto uuval = mat.getColumns(x);
      int uu = abs(uuval);
      int ii = uu%16;
      int vv = (uu/16)%4;
      int rr = (uu/16)/4;
      if (vv == 2) valX[rr][ii] = 1;
      else valK[rr][ii] = 1;
      if (mat.setColumnAsPivot(x, line1, line2)) ++line1;
      else ++l;
      if (mat.setAsPivot(-uuval, line1, line2)) ++line1;
      auto pbool = updateColumns(vv, rr, ii, valX, valK, line1, line2, mat, valColX, valColSR, valColK);
      if (!pbool.second) return false;
      if (pbool.first) {
        if (!propagateZERO(valX, valK, line1, line2, mat, valColX, valColSR, valColK)) return false;
        return propagateONE(valX, valK, 0, line1, line2, mat, valColX, valColSR, valColK);
      }
    }
    else {
      if (tmp == 0) {return false;}
      else ++l;
    }
  }
  if (l == line2) return true;
  set<vector<pair<GFElement, int>>> mymap;
  for (l = line1_start; l < line1; ++l) {
    vector<pair<GFElement, int>> myvec;
    myvec.reserve(mat.nbcols);
    for (unsigned c = 0; c < mat.nbcols; ++c) {
      if (mat(l, c) != 0) {
        myvec.emplace_back(mat(l, c), mat.getColumns(c));
      }
    }
    if (myvec.empty()) return false;
    auto coef = myvec[0].first.getInverse();
    for (auto & p : myvec) p.first *= coef;
    mymap.emplace(move(myvec));
  }
  vector<int> toprocess;
  vector<pair<GFElement, int>> myvec;
  myvec.reserve(mat.nbcols);
  for (l = line1; l < line2; ++l) {
    myvec.clear();
    for (unsigned c = 0; c < mat.nbcols; ++c) {
      if (mat(l, c) != 0) {
        myvec.emplace_back(mat(l, c), mat.getColumns(c));
      }
    }
    if (myvec.empty()) continue;
    auto coef = myvec[0].first.getInverse();
    for (auto & p : myvec) p.first *= coef;
    if (mymap.count(myvec) != 0) toprocess.emplace_back(mat.getFront(l));
  }
  if (toprocess.empty()) return true;

  for (auto x : toprocess) {
    if (mat.setAsPivot(x, line1, line2)) ++line1;
    int uu = abs(x);
    int ii = uu%16;
    int vv = (uu/16)%4;
    int rr = (uu/16)/4;
    if (vv == 2) valX[rr][ii] = 1;
    else valK[rr][ii] = 1;
    auto pbool = updateColumns(vv, rr, ii, valX, valK, line1, line2, mat, valColX, valColSR, valColK);
    if (!pbool.second) return false;
    if (pbool.first) {
      for (auto xx : toprocess) {
        if (mat.setAsPivot(-xx, line1, line2)) ++line1;
        if (xx == x) break;
      }
      if (!propagateZERO(valX, valK, line1, line2, mat, valColX, valColSR, valColK)) return false;
      return propagateONE(valX, valK, 0, line1, line2, mat, valColX, valColSR, valColK);
    }
  }
  line1_start = line1;
  for (auto x : toprocess) {
    if (mat.setAsPivot(-x, line1, line2)) ++line1;
  }
  return propagateONE(valX, valK, line1_start, line1, line2, mat, valColX, valColSR, valColK);
}

bool propagateONE(vector<vector<uint8_t>> & valX, vector<vector<uint8_t>> & valK, unsigned & line1, unsigned & line2, Matrix & mat, vector<vector<uint8_t>> const & valColX, vector<vector<uint8_t>> const & valColSR, vector<vector<uint8_t>> const & valColK) {
  return propagateONE(valX, valK, 0, line1, line2, mat, valColX, valColSR, valColK);
}


bool updateStateVar(uint8_t x, unsigned r, unsigned l, unsigned c, Matrix & mat, unsigned & line1, unsigned & line2, vector<vector<uint8_t>> & valX, vector<vector<uint8_t>> & valK, vector<vector<uint8_t>> const & valColX, vector<vector<uint8_t>> const & valColSR, vector<vector<uint8_t>> const & valColK) {
  if (valX[r][4*l + c] != 2) return valX[r][4*l + c] == x;
  int uval = 16*(4*r + 2) + 4*l + c;
  if (x == 0) {
    set0Mat(uval, valX, valK, line1, line2, mat);
    set0Mat(-uval, valX, valK, line1, line2, mat);
    valX[r][4*l + c] = 0;
    //cout << r << ": " << 4*l + c << " (0)" << endl;
    //getchar();
    auto pbool = updateColumns(2, r, 4*l + c, valX, valK, line1, line2, mat, valColX, valColSR, valColK);
    if (!pbool.second) return false;
    if (!propagateZERO(valX, valK, line1, line2, mat, valColX, valColSR, valColK)) return false;
  }
  else {
    if (mat.setAsPivot(uval, line1, line2)) ++line1;
    if (mat.setAsPivot(-uval, line1, line2)) ++line1;
    valX[r][4*l + c] = 1;
    //cout << r << ": " << 4*l + c << " (1)" << endl;
    auto pbool = updateColumns(2, r, 4*l + c, valX, valK, line1, line2, mat, valColX, valColSR, valColK);
    if (!pbool.second) return false;
    if (pbool.first && !propagateZERO(valX, valK, line1, line2, mat, valColX, valColSR, valColK)) return false;
  }
  return propagateONE(valX, valK, line1, line2, mat, valColX, valColSR, valColK);
}

bool updateKeyVar(uint8_t x, unsigned r, unsigned l, unsigned c, Matrix & mat, unsigned & line1, unsigned & line2, vector<vector<uint8_t>> & valX, vector<vector<uint8_t>> & valK, vector<vector<uint8_t>> const & valColX, vector<vector<uint8_t>> const & valColSR, vector<vector<uint8_t>> const & valColK) {
  if (valK[r][4*l + c] != 2) return valK[r][4*l + c] == x;
  int uval = 16*(4*r + 1) + 4*l + c;
  if (x == 0) {
    set0Mat(uval, valX, valK, line1, line2, mat);
    set0Mat(-uval, valX, valK, line1, line2, mat);
    valK[r][4*l + c] = 0;
    auto pbool = updateColumns(1, r, 4*l + c, valX, valK, line1, line2, mat, valColX, valColSR, valColK);
    if (!pbool.second) return false;
    if (!propagateZERO(valX, valK, line1, line2, mat, valColX, valColSR, valColK)) return false;
  }
  else {
    if (mat.setAsPivot(uval, line1, line2)) ++line1;
    if (mat.setAsPivot(-uval, line1, line2)) ++line1;
    valK[r][4*l + c] = 1;
    auto pbool = updateColumns(1, r, 4*l + c, valX, valK, line1, line2, mat, valColX, valColSR, valColK);
    if (!pbool.second) return false;
    if (pbool.first && !propagateZERO(valX, valK, line1, line2, mat, valColX, valColSR, valColK)) return false;
  }
  return propagateONE(valX, valK, line1, line2, mat, valColX, valColSR, valColK);
}


bool updateColX(unsigned r, unsigned c, vector<vector<uint8_t>> & valX, vector<vector<uint8_t>> & valK, unsigned & line1, unsigned & line2, Matrix & mat, vector<vector<uint8_t>> & valColX, vector<vector<uint8_t>> & valColSR, vector<vector<uint8_t>> & valColK) {
  unsigned cpt[3] = {0,0,0};
  for (unsigned l = 0; l < 4; ++l) cpt[valX[r][4*l+c]] += 1;
  if (cpt[0] > 4-valColX[r][c] || cpt[1] > valColX[r][c]) return false;
  if (cpt[2] == 0) return true;

  if (cpt[1] == valColX[r][c]) {
    for (unsigned l = 0; l < 4; ++l) {
      if (valX[r][4*l + c] != 2) continue;
      if (!updateStateVar(0, r, l, c, mat, line1, line2, valX, valK, valColX, valColSR, valColK)) return false;
    }
  }
  else if (cpt[0] == 4 - valColX[r][c]) {
    for (unsigned l = 0; l < 4; ++l) {
      if (valX[r][4*l + c] != 2) continue;
      if (!updateStateVar(1, r, l, c, mat, line1, line2, valX, valK, valColX, valColSR, valColK)) return false;
    }
  }
  return true;
}

bool updateColSR(unsigned r, unsigned c, vector<vector<uint8_t>> & valX, vector<vector<uint8_t>> & valK, unsigned & line1, unsigned & line2, Matrix & mat, vector<vector<uint8_t>> & valColX, vector<vector<uint8_t>> & valColSR, vector<vector<uint8_t>> & valColK) {
  unsigned cpt[3] = {0,0,0};
  for (unsigned l = 0; l < 4; ++l) cpt[valX[r][4*l+((c+l)%4)]] += 1;
  if (cpt[0] > 4-valColSR[r][c] || cpt[1] > valColSR[r][c]) return false;
  if (cpt[2] == 0) return true;
  if (cpt[1] == valColSR[r][c]) {
    for (unsigned l = 0; l < 4; ++l) {
      if (valX[r][4*l + ((c+l)%4)] != 2) continue;
      if (!updateStateVar(0, r, l, (c+l)%4, mat, line1, line2, valX, valK, valColX, valColSR, valColK)) return false;
    }
  }
  else if (cpt[0] == 4 - valColSR[r][c]) {
    for (unsigned l = 0; l < 4; ++l) {
      if (valX[r][4*l + ((c+l)%4)] != 2) continue;
      if (!updateStateVar(1, r, l, (c+l)%4, mat, line1, line2, valX, valK, valColX, valColSR, valColK)) return false;
    }
  }
  return true;
}

bool updateColK(unsigned r, unsigned c, vector<vector<uint8_t>> & valX, vector<vector<uint8_t>> & valK, unsigned & line1, unsigned & line2, Matrix & mat, vector<vector<uint8_t>> & valColX, vector<vector<uint8_t>> & valColSR, vector<vector<uint8_t>> & valColK) {
  unsigned cpt[3] = {0,0,0};
  for (unsigned l = 0; l < 4; ++l) cpt[valK[r][4*l+c]] += 1;
  if (cpt[0] > 4-valColK[r][c] || cpt[1] > valColK[r][c]) return false;
  if (cpt[2] == 0) return true;
  if (cpt[1] == valColK[r][c]) {
    for (unsigned l = 0; l < 4; ++l) {
      if (valK[r][4*l + c] != 2) continue;
      if (!updateKeyVar(0, r, l, c, mat, line1, line2, valX, valK, valColX, valColSR, valColK)) return false;
    }
  }
  else if (cpt[0] == 4 - valColK[r][c]) {
    for (unsigned l = 0; l < 4; ++l) {
      if (valK[r][4*l + c] != 2) continue;
      if (!updateKeyVar(1, r, l, c, mat, line1, line2, valX, valK, valColX, valColSR, valColK)) return false;
    }
  }
  return true;
}



pair<bool,bool> constraintMC(unsigned rk, unsigned ck, unsigned deck, unsigned r1, unsigned c1, unsigned r2, unsigned c2, Matrix & mat, unsigned & line1, unsigned & line2, vector<vector<uint8_t>> & valX, vector<vector<uint8_t>> & valK, vector<vector<uint8_t>> & valColX, vector<vector<uint8_t>> & valColK, vector<vector<uint8_t>> & valColSR) {
  pair<bool, bool> res = make_pair(true, false);
  //return res;
  if (valColK[r1][c1] != 0 && valColK[r2][c2] != 0) {
    int x0 = 4;
    for (unsigned l = 0; l < 4; ++l) if (valX[r1-1][4*l + ((c1 + l)%4)] == 0 && valX[r2-1][4*l + ((c2 + l)%4)] == 0) --x0;
    int x1 = valColSR[r1-1][c1] + valColSR[r2-1][c2];
    for (unsigned l = 0; l < 4; ++l) if (valX[r1-1][4*l + ((c1 + l)%4)] == 1 && valX[r2-1][4*l + ((c2 + l)%4)] == 1) --x1;
    unsigned x = min(x0, x1);
    int y0 = 4;
    for (unsigned l = 0; l < 4; ++l) if (valX[r1][4*l + c1] == 0 && valX[r2][4*l + c2] == 0 && valK[rk][4*((l+deck)%4) + ck] == 0) --y0;
    int y1 = valColX[r1][c1] + valColX[r2][c2] + valColK[rk][ck];
    for (unsigned l = 0; l < 4; ++l) {
      int tmp = 0;
      if (valX[r1][4*l + c1] == 1) tmp += 1;
      if (valX[r2][4*l + c2] == 1) tmp += 1;
      if (valK[rk][4*((l+deck)%4) + ck] == 1) tmp += 1;
      if (tmp != 0) y1 -= (tmp-1);
    }
    unsigned y = min(y0, y1);
    //cout << x0 << " " << x1 << " - " << y0 << " " << y1 << endl;
    if (x + y <= 4) {
      if (valColSR[r1-1][c1] != valColSR[r2-1][c2]) {res.first = false; return res;}
      if (valColX[r1][c1] + valColX[r2][c2] < valColK[rk][ck] || abs(valColX[r1][c1] - valColX[r2][c2]) > valColK[rk][ck]) {res.first = false; return res;}
      for (unsigned l = 0; l < 4; ++l) if (valX[r1-1][4*l + ((c1 + l)%4)] + valX[r2-1][4*l + ((c2+l)%4)] == 1) {res.first = false; return res;}
      for (unsigned l = 0; l < 4; ++l) {
        if (valX[r1-1][4*l + ((c1 + l)%4)] != valX[r2-1][4*l + ((c2+l)%4)]) {
          if (valX[r1-1][4*l + ((c1 + l)%4)] == 2) {
            if (!updateStateVar(valX[r2-1][4*l + ((c2+l)%4)], r1-1, l, ((c1 + l)%4), mat, line1, line2, valX, valK, valColX, valColSR, valColK)) {res.first = false; return res;}
            res.second = true;
          }
          else {
            if (!updateStateVar(valX[r1-1][4*l + ((c1+l)%4)], r2-1, l, ((c2 + l)%4), mat, line1, line2, valX, valK, valColX, valColSR, valColK)) {res.first = false; return res;}
            res.second = true;
          }
        }
      }
      {
        for (unsigned l = 0; l < 4; ++l) if (valX[r1][4*l + c1] + valX[r2][4*l + c2] + valK[rk][4*((l+deck)%4) + ck] == 1) {res.first = false; return res;}
        for (unsigned l = 0; l < 4; ++l) {
          if (valK[rk][4*((l+deck)%4) + ck] == 0 && valX[r1][4*l + c1] != valX[r2][4*l + c2]) {
            if (valX[r1][4*l + c1] == 2) {
              if (!updateStateVar(valX[r2][4*l + c2], r1, l, c1, mat, line1, line2, valX, valK, valColX, valColSR, valColK)) {res.first = false; return res;}
              res.second = true;
            }
            else {
              if (!updateStateVar(valX[r1][4*l + c1], r2, l, c2, mat, line1, line2, valX, valK, valColX, valColSR, valColK)) {res.first = false; return res;}
              res.second = true;
            }
          }
          else if (valX[r2][4*l + c2] == 0 && valX[r1][4*l + c1] != valK[rk][4*((l+deck)%4) + ck]) {
            if (valX[r1][4*l + c1] == 2) {
              if (!updateStateVar(valK[rk][4*((l+deck)%4) + ck], r1, l, c1, mat, line1, line2, valX, valK, valColX, valColSR, valColK)) {res.first = false; return res;}
              res.second = true;
            }
            else {
              if (!updateKeyVar(valX[r1][4*l + c1], rk, ((l+deck)%4), ck, mat, line1, line2, valX, valK, valColX, valColSR, valColK)) {res.first = false; return res;}
              res.second = true;
            }
          }
          else if (valX[r1][4*l + c1] == 0 && valX[r2][4*l + c2] != valK[rk][4*((l+deck)%4) + ck]) {
            if (valX[r2][4*l + c2] == 2) {
              if (!updateStateVar(valK[rk][4*((l+deck)%4) + ck], r2, l, c2, mat, line1, line2, valX, valK, valColX, valColSR, valColK)) {res.first = false; return res;}
              res.second = true;
            }
            else {
              if (!updateKeyVar(valX[r2][4*l + c2], rk, ((l+deck)%4), ck, mat, line1, line2, valX, valK, valColX, valColSR, valColK)) {res.first = false; return res;}
              res.second = true;
            }
          }
        }
      }
    }
  }
  return res;
}

bool findBestTrail1(unsigned state_key, vector<vector<uint8_t>> const & T, uint8_t & global_bound, uint8_t current_bound, int step, Matrix & mat, unsigned line1, unsigned line2, vector<vector<uint8_t>> & valX, vector<vector<uint8_t>> & valK, vector<vector<uint8_t>> & valColX, vector<vector<uint8_t>> & valColK, vector<vector<uint8_t>> & valColSR) {
  static unsigned const n_states = 5*5*5*5;
  static unsigned const n_keys = 5*5*5*5*5*5;
  static unsigned const mypow[8] = {1, 5, 5*5, 5*5*5, 5*5*5*5, 5*5*5*5*5, 5*5*5*5*5*5, 5*5*5*5*5*5*5};

  static unsigned cpt = 0;

  int mod_step = step%3;


  vector<unsigned> next_state_key;

  if (step < -1) {
    for (unsigned r = 0; r < valColK.size(); ++r) {
      for (unsigned c = 0; c < 4; ++c) {
        if (valColK[r][c] > 0 && valColK[r][c] < 4) {
          for (unsigned l = 0; l < 4; ++l) {
            if (valK[r][4*l + c] == 2) {
              {
                auto m = mat;
                auto vX = valX;
                auto vK = valK;
                auto l1 = line1;
                auto l2 = line2;
                if (updateKeyVar(0, r, l, c, m, l1, l2, vX, vK, valColX, valColSR, valColK)) {
                  if (findBestTrail1(state_key, T, global_bound, current_bound, step, m, l1, l2, vX, vK, valColX, valColK, valColSR)) return true;
                }
              }
              if (updateKeyVar(1, r, l, c, mat, line1, line2, valX, valK, valColX, valColSR, valColK)) {
                return findBestTrail1(state_key, T, global_bound, current_bound, step, mat, line1, line2, valX, valK, valColX, valColK, valColSR);
              }
              else return false;
            }
          }
        }
      }
    }
    for (unsigned r = 0; r < valColX.size(); ++r) {
      for (unsigned l = 0; l < 4; ++l) {
        for (unsigned c = 0; c < 4; ++c) {
          if (valX[r][4*l + c] == 2) {
            {
              auto m = mat;
              auto vX = valX;
              auto vK = valK;
              auto l1 = line1;
              auto l2 = line2;
              if (updateStateVar(0, r, l, c, m, l1, l2, vX, vK, valColX, valColSR, valColK)) {
                if (findBestTrail1(state_key, T, global_bound, current_bound, step, m, l1, l2, vX, vK, valColX, valColK, valColSR)) return true;
              }
            }
            if (updateStateVar(1, r, l, c, mat, line1, line2, valX, valK, valColX, valColSR, valColK)) {
              return findBestTrail1(state_key, T, global_bound, current_bound, step, mat, line1, line2, valX, valK, valColX, valColK, valColSR);
            }
            else return false;
          }
        }
      }
    }
  }


  if (step == -1) {
    for (unsigned c = 0; c < 4; ++c) current_bound += valColSR[0][c];
    //current_bound += valColK[0][3];
    if (current_bound != global_bound) return false;
    --step;

    for (unsigned i = 0; i < 16; ++i) {
      if (valK[0][i] != 0 || valX[0][i] != 0 || valK[1][i] != 0) {
        auto m = mat;
        auto vX = valX;
        auto vK = valK;
        auto l1 = line1;
        auto l2 = line2;
        return findBestTrail1(state_key, T, global_bound, current_bound, step, m, l1, l2, vX, vK, valColX, valColK, valColSR);
      }
    }

    return false;
  }


  if (step < 0) {
    return true;
  }

  unsigned r = (step+1)/3;

  unsigned colk = (4*((r+1)%3))%6;

  if (mod_step == 2 && r+1 < valColK.size()) {
    bool restart = true;
    while (restart) {
      restart = false;
      for (unsigned rr = r; rr <= r+1; ++rr) {
        for (unsigned c1 = 0; c1 < 4; ++c1) {
          if (rr+1 >= valColK.size() || valColK[rr+1][c1] == 0) continue;
          unsigned pos1 = (4*(rr+1) + c1)%6;
          if (pos1 == 0) {
            unsigned rk = rr+1;
            unsigned ck = c1 + 5;
            if (ck >= 4) {ck -= 4; rk += 1;}
            unsigned r2 = rk;
            unsigned c2 = ck+1;
            if (c2 >= 4) {c2 -= 4; r2 += 1;}
            if (r2 >= valColK.size()) continue;
            auto flag = constraintMC(rk, ck, 1, rr+1, c1, r2, c2, mat, line1, line2, valX, valK, valColX, valColK, valColSR);
            if (!flag.first) return false;
            if (flag.second) {restart = true; break;}
          }
          if (pos1 != 5) {
            unsigned r2 = rr+1;
            unsigned c2 = c1 + 1;
            if (c2 >= 4) {c2 -= 4; r2 += 1;}
            if (r2 >= valColK.size()) continue;
            unsigned rk = r2 - 2;
            unsigned ck = c2 + 2;
            if (ck >= 4) {ck -= 4; rk += 1;}
            if (rk >= valColK.size() || valColK[rk][ck] == 5) continue;
            auto flag = constraintMC(rk, ck, 0, rr+1, c1, r2, c2, mat, line1, line2, valX, valK, valColX, valColK, valColSR);
            if (!flag.first) return false;
            if (flag.second) {restart = true; break;}
          }
        }
        if (restart) break;
      }
    }
  }

  if (mod_step == 2) {
    next_state_key = inv_updateSR(T[step], global_bound-current_bound, state_key%n_states, state_key/n_states, valX[r]);
    for (auto f : next_state_key) {
      for (unsigned c = 0; c < 4; ++c) {
        valColX[r][c] = (f/mypow[c])%5;
      }
      auto m = mat;
      auto vX = valX;
      auto vK = valK;
      auto l1 = line1;
      auto l2 = line2;
      bool isvalid = true;
      for (unsigned c = 0; (c < 4) && isvalid; ++c) isvalid = updateColX(r, c, vX, vK, l1, l2, m, valColX, valColSR, valColK);
      if (isvalid) {
        if (findBestTrail1(f, T, global_bound, current_bound, step-1, m, l1, l2, vX, vK, valColX, valColK, valColSR)) return true;
      }
    }
    for (unsigned c = 0; c < 4; ++c) {
      valColX[r][c] = 5;
    }
  }
  else {
    if (mod_step == 0) {

      if (r == 0) return findBestTrail1(state_key, T, global_bound, current_bound, step-1, mat, line1, line2, valX, valK, valColX, valColK, valColSR);

      next_state_key = inv_updateKey192((colk + 3)%6, (r == 1) ? 2 : 4, T[step], global_bound-current_bound, state_key%n_states, state_key/n_states);

      //if (r == 1 && !next_state_key.empty()) return findBestTrail(next_state_key[0], T, global_bound, current_bound, step-1, mat, line1, line2, valX, valK, valColX, valColK, valColSR);

      //cout << "state_key: "; for (unsigned c = 0; c < 6; ++c) {cout << ((state_key/n_states)/mypow[c])%5 << " | ";} cout << endl;
      for (auto f : next_state_key) {
        //cout << "f: "; for (unsigned c = 0; c < 6; ++c) {cout << ((f/n_states)/mypow[c])%5 << " | ";} cout << endl;
        if (r > 1) for (unsigned c = 0; c < 2; ++c) valColK[r-1][c+2] = ((f/n_states)/mypow[(colk + c)%6])%5;
        for (unsigned c = 0; c < 2; ++c) valColK[r][c] = ((f/n_states)/mypow[(colk + 2 + c)%6])%5;
        //cout << r << ": ";
        //for (unsigned c = 0; c < 4; ++c) cout << (unsigned) valColK[r][c] << "| ";
        //cout << endl;
        auto m = mat;
        auto vX = valX;
        auto vK = valK;
        auto l1 = line1;
        auto l2 = line2;
        bool isvalid = true;
        for (unsigned c = 0; (c < 2) && isvalid; ++c) isvalid = updateColK(r-1, c+2, vX, vK, l1, l2, m, valColX, valColSR, valColK);
        for (unsigned c = 0; (c < 2) && isvalid; ++c) isvalid = updateColK(r, c, vX, vK, l1, l2, m, valColX, valColSR, valColK);
        if (isvalid) if (findBestTrail1(f, T, global_bound, current_bound, step-1, m, l1, l2, vX, vK, valColX, valColK, valColSR)) return true;
      }
      for (unsigned c = 0; c < 2; ++c) valColK[r][c] = 5;
      for (unsigned c = 0; c < 2; ++c) valColK[r-1][c+2] = 5;
    }
    else {
      next_state_key = inv_updateMC_ARK(T[step], colk, global_bound-current_bound, state_key%n_states, state_key/n_states);
      // cout << "c: " << (unsigned) (c + dec_key) << endl;
      // cout << "nb sol: " << next_state_key.size() << endl;

      for (auto f : next_state_key) {
        for (unsigned c = 0; c < 4; ++c) valColSR[r][c] = ((f%n_states)/mypow[c])%5;

        auto m = mat;
        auto vX = valX;
        auto vK = valK;
        auto l1 = line1;
        auto l2 = line2;
        bool isvalid = true;
        for (unsigned c = 0; (c < 4) && isvalid; ++c) isvalid = updateColSR(r, c, vX, vK, l1, l2, m, valColX, valColSR, valColK);
        if (isvalid) {
          unsigned cost = current_bound + valColX[r+1][0] + valColX[r+1][1] + valColX[r+1][2] + valColX[r+1][3];
          if (colk >= 2) cost += ((f/n_states)/mypow[5])%5;
          if (findBestTrail1(f, T, global_bound, cost, step-1, m, l1, l2, vX, vK, valColX, valColK, valColSR)) return true;
        }
      }
      for (unsigned c = 0; c < 4; ++c) valColSR[r][c] = 5;
    }
  }

  return false;
}

bool findBestTrail2(vector<uint8_t> & myvec, unsigned state_key, vector<vector<uint8_t>> const & T, uint8_t & global_bound, uint8_t current_bound, int step, Matrix & mat, unsigned line1, unsigned line2, vector<vector<uint8_t>> & valX, vector<vector<uint8_t>> & valK, vector<vector<uint8_t>> & valColX, vector<vector<uint8_t>> & valColK, vector<vector<uint8_t>> & valColSR) {
  static unsigned const n_states = 5*5*5*5;
  static unsigned const n_keys = 5*5*5*5*5*5;
  static unsigned const mypow[8] = {1, 5, 5*5, 5*5*5, 5*5*5*5, 5*5*5*5*5, 5*5*5*5*5*5, 5*5*5*5*5*5*5};

  static unsigned cpt = 0;

  int mod_step = step%3;


  vector<unsigned> next_state_key;


  unsigned r = (step+1)/3;

  unsigned colk = (4*((r+1)%3))%6;

  if (mod_step == 2 && r+1 < valColK.size()) {
    bool restart = true;
    while (restart) {
      restart = false;
      for (unsigned rr = r; rr <= r+1; ++rr) {
        for (unsigned c1 = 0; c1 < 4; ++c1) {
          if (rr+1 >= valColK.size() || valColK[rr+1][c1] == 0) continue;
          unsigned pos1 = (4*(rr+1) + c1)%6;
          if (pos1 == 0) {
            unsigned rk = rr+1;
            unsigned ck = c1 + 5;
            if (ck >= 4) {ck -= 4; rk += 1;}
            unsigned r2 = rk;
            unsigned c2 = ck+1;
            if (c2 >= 4) {c2 -= 4; r2 += 1;}
            if (r2 >= valColK.size()) continue;
            auto flag = constraintMC(rk, ck, 1, rr+1, c1, r2, c2, mat, line1, line2, valX, valK, valColX, valColK, valColSR);
            if (!flag.first) return false;
            if (flag.second) {restart = true; break;}
          }
          if (pos1 != 5) {
            unsigned r2 = rr+1;
            unsigned c2 = c1 + 1;
            if (c2 >= 4) {c2 -= 4; r2 += 1;}
            if (r2 >= valColK.size()) continue;
            unsigned rk = r2 - 2;
            unsigned ck = c2 + 2;
            if (ck >= 4) {ck -= 4; rk += 1;}
            if (rk >= valColK.size() || valColK[rk][ck] == 5) continue;
            auto flag = constraintMC(rk, ck, 0, rr+1, c1, r2, c2, mat, line1, line2, valX, valK, valColX, valColK, valColSR);
            if (!flag.first) return false;
            if (flag.second) {restart = true; break;}
          }
        }
        if (restart) break;
      }
    }
  }

  if (mod_step == 2) {
    next_state_key = inv_updateSR(T[step], global_bound-current_bound, state_key%n_states, state_key/n_states, valX[r]);
    for (auto f : next_state_key) {
      for (unsigned c = 0; c < 4; ++c) {
        valColX[r][c] = (f/mypow[c])%5;
      }
      auto m = mat;
      auto vX = valX;
      auto vK = valK;
      auto l1 = line1;
      auto l2 = line2;
      bool isvalid = true;
      for (unsigned c = 0; (c < 4) && isvalid; ++c) isvalid = updateColX(r, c, vX, vK, l1, l2, m, valColX, valColSR, valColK);
      if (isvalid) {
        #pragma omp critical
        {
          if (myvec[f] < global_bound-current_bound) myvec[f] = global_bound-current_bound;
        }

      }
    }
    for (unsigned c = 0; c < 4; ++c) {
      valColX[r][c] = 5;
    }
  }
  else {
    if (mod_step == 0) {

      if (r == 0) return findBestTrail2(myvec, state_key, T, global_bound, current_bound, step-1, mat, line1, line2, valX, valK, valColX, valColK, valColSR);

      next_state_key = inv_updateKey192((colk + 3)%6, (r == 1) ? 2 : 4, T[step], global_bound-current_bound, state_key%n_states, state_key/n_states);

      //if (r == 1 && !next_state_key.empty()) return findBestTrail(next_state_key[0], T, global_bound, current_bound, step-1, mat, line1, line2, valX, valK, valColX, valColK, valColSR);

      //cout << "state_key: "; for (unsigned c = 0; c < 6; ++c) {cout << ((state_key/n_states)/mypow[c])%5 << " | ";} cout << endl;
      for (auto f : next_state_key) {
        //cout << "f: "; for (unsigned c = 0; c < 6; ++c) {cout << ((f/n_states)/mypow[c])%5 << " | ";} cout << endl;
        if (r > 1) for (unsigned c = 0; c < 2; ++c) valColK[r-1][c+2] = ((f/n_states)/mypow[(colk + c)%6])%5;
        for (unsigned c = 0; c < 2; ++c) valColK[r][c] = ((f/n_states)/mypow[(colk + 2 + c)%6])%5;
        //cout << r << ": ";
        //for (unsigned c = 0; c < 4; ++c) cout << (unsigned) valColK[r][c] << "| ";
        //cout << endl;
        auto m = mat;
        auto vX = valX;
        auto vK = valK;
        auto l1 = line1;
        auto l2 = line2;
        bool isvalid = true;
        for (unsigned c = 0; (c < 2) && isvalid; ++c) isvalid = updateColK(r-1, c+2, vX, vK, l1, l2, m, valColX, valColSR, valColK);
        for (unsigned c = 0; (c < 2) && isvalid; ++c) isvalid = updateColK(r, c, vX, vK, l1, l2, m, valColX, valColSR, valColK);
        if (isvalid) findBestTrail2(myvec, f, T, global_bound, current_bound, step-1, m, l1, l2, vX, vK, valColX, valColK, valColSR);
      }
      for (unsigned c = 0; c < 2; ++c) valColK[r][c] = 5;
      for (unsigned c = 0; c < 2; ++c) valColK[r-1][c+2] = 5;
    }
    else {
      next_state_key = inv_updateMC_ARK(T[step], colk, global_bound-current_bound, state_key%n_states, state_key/n_states);
      // cout << "c: " << (unsigned) (c + dec_key) << endl;
      // cout << "nb sol: " << next_state_key.size() << endl;

      for (auto f : next_state_key) {
        for (unsigned c = 0; c < 4; ++c) valColSR[r][c] = ((f%n_states)/mypow[c])%5;

        auto m = mat;
        auto vX = valX;
        auto vK = valK;
        auto l1 = line1;
        auto l2 = line2;
        bool isvalid = true;
        for (unsigned c = 0; (c < 4) && isvalid; ++c) isvalid = updateColSR(r, c, vX, vK, l1, l2, m, valColX, valColSR, valColK);
        if (isvalid) {
          unsigned cost = current_bound + valColX[r+1][0] + valColX[r+1][1] + valColX[r+1][2] + valColX[r+1][3];
          if (colk >= 2) cost += ((f/n_states)/mypow[5])%5;
          findBestTrail2(myvec, f, T, global_bound, cost, step-1, m, l1, l2, vX, vK, valColX, valColK, valColSR);
        }
      }
      for (unsigned c = 0; c < 4; ++c) valColSR[r][c] = 5;
    }
  }

  return false;
}

bool flag_solution_found = false;

void findBestTrail(unsigned state_key, vector<vector<uint8_t>> const & T, uint8_t & global_bound, uint8_t current_bound, int step, Matrix & mat, unsigned line1, unsigned line2, vector<vector<uint8_t>> & valX, vector<vector<uint8_t>> & valK, vector<vector<uint8_t>> & valColX, vector<vector<uint8_t>> & valColK, vector<vector<uint8_t>> & valColSR) {
  static unsigned const n_states = 5*5*5*5;
  static unsigned const n_keys = 5*5*5*5*5*5;
  static unsigned const mypow[8] = {1, 5, 5*5, 5*5*5, 5*5*5*5, 5*5*5*5*5, 5*5*5*5*5*5, 5*5*5*5*5*5*5};

  static unsigned cpt = 0;

  int mod_step = step%3;

  /*{
    cout << "bound0: " << (unsigned) current_bound << " (" << step << ")" << endl;
    cout << "    ";
    for (unsigned r = 0; r < valColX.size(); ++r) {
      if (r > 0) {
        cout << "|";
        for (unsigned c = 0; c < 4; ++c) cout << (unsigned) valColX[r][c] << "|";
        if (r < valColX.size()-1) cout << " --> ";
      }
      if (r < valColX.size()-1) {
        cout << "|";
        for (unsigned c = 0; c < 4; ++c) cout << (unsigned) valColSR[r][c] << "|";
        cout << " --> ";
      }
    }
    cout << endl;
    for (unsigned r = 0; r < valColK.size(); ++r) {
      cout << "|";
      for (unsigned c = 0; c < 4; ++c) cout << (unsigned) valColK[r][c] << "|";
      if (r < valColK.size()-1) cout << "     -->      ";
    }
    cout << endl;
    cout << " --------------- " << endl;
    for (unsigned r = 0; r < valColX.size(); ++r) {
      for (unsigned l = 0; l < 4; ++l) {
        for (unsigned c = 0; c < 4; ++c) {
          cout << (unsigned) valX[r][4*l + c] << " ";
        }
        cout << " | ";
        for (unsigned c = 0; c < 4; ++c) {
          cout << (unsigned) valK[r][4*l + c] << " ";
        }
        cout << endl;
      }
      cout << endl;
    }
    getchar();
  }*/

  vector<unsigned> next_state_key;

  if (step < -1) {
    for (unsigned r = 0; r < valColK.size(); ++r) {
      for (unsigned c = 0; c < 4; ++c) {
        if (valColK[r][c] > 0 && valColK[r][c] < 4) {
          for (unsigned l = 0; l < 4; ++l) {
            if (valK[r][4*l + c] == 2) {
              {
                auto m = mat;
                auto vX = valX;
                auto vK = valK;
                auto l1 = line1;
                auto l2 = line2;
                if (updateKeyVar(0, r, l, c, m, l1, l2, vX, vK, valColX, valColSR, valColK)) {
                  findBestTrail(state_key, T, global_bound, current_bound, step, m, l1, l2, vX, vK, valColX, valColK, valColSR);
                }
              }
              if (updateKeyVar(1, r, l, c, mat, line1, line2, valX, valK, valColX, valColSR, valColK)) {
                return findBestTrail(state_key, T, global_bound, current_bound, step, mat, line1, line2, valX, valK, valColX, valColK, valColSR);
              }
              else return;
            }
          }
        }
      }
    }
    for (unsigned r = 0; r < valColX.size(); ++r) {
      for (unsigned l = 0; l < 4; ++l) {
        for (unsigned c = 0; c < 4; ++c) {
          if (valX[r][4*l + c] == 2) {
            {
              auto m = mat;
              auto vX = valX;
              auto vK = valK;
              auto l1 = line1;
              auto l2 = line2;
              if (updateStateVar(0, r, l, c, m, l1, l2, vX, vK, valColX, valColSR, valColK)) {
                findBestTrail(state_key, T, global_bound, current_bound, step, m, l1, l2, vX, vK, valColX, valColK, valColSR);
              }
            }
            if (updateStateVar(1, r, l, c, mat, line1, line2, valX, valK, valColX, valColSR, valColK)) {
              return findBestTrail(state_key, T, global_bound, current_bound, step, mat, line1, line2, valX, valK, valColX, valColK, valColSR);
            }
            else return;
          }
        }
      }
    }
  }


  if (step == -1) {
    for (unsigned c = 0; c < 4; ++c) current_bound += valColSR[0][c];
    //current_bound += valColK[0][3];
    if (current_bound != global_bound) return;
    --step;

    for (unsigned i = 0; i < 16; ++i) {
      if (valK[0][i] != 0 || valX[0][i] != 0 || valK[1][i] != 0) {
        auto m = mat;
        auto vX = valX;
        auto vK = valK;
        auto l1 = line1;
        auto l2 = line2;
        return findBestTrail(state_key, T, global_bound, current_bound, step, m, l1, l2, vX, vK, valColX, valColK, valColSR);
      }
    }

    return;
  }


  if (step < 0) {
    #pragma omp critical
    {
      cout << "bound: " << (unsigned) current_bound << " (" << ++cpt << ")" << endl;
      cout << "    ";
      for (unsigned r = 0; r < valColX.size(); ++r) {
        if (r > 0) {
          cout << "|";
          for (unsigned c = 0; c < 4; ++c) cout << (unsigned) valColX[r][c] << "|";
          if (r < valColX.size()-1) cout << " --> ";
        }
        if (r < valColX.size()-1) {
          cout << "|";
          for (unsigned c = 0; c < 4; ++c) cout << (unsigned) valColSR[r][c] << "|";
          cout << " --> ";
        }
      }
      cout << endl;
      for (unsigned r = 0; r < valColK.size(); ++r) {
        cout << "|";
        for (unsigned c = 0; c < 4; ++c) cout << (unsigned) valColK[r][c] << "|";
        if (r < valColK.size()-1) cout << "     -->      ";
      }
      cout << endl;
      cout << " --------------- " << endl;
      for (unsigned r = 0; r < valColX.size(); ++r) {
        for (unsigned l = 0; l < 4; ++l) {
          for (unsigned c = 0; c < 4; ++c) {
            cout << (unsigned) valX[r][4*l + c] << " ";
          }
          cout << " | ";
          for (unsigned c = 0; c < 4; ++c) {
            cout << (unsigned) valK[r][4*l + c] << " ";
          }
          cout << endl;
        }
        cout << endl;
      }
      flag_solution_found = true;
    }
    return;
  }

  unsigned r = (step+1)/3;

  unsigned colk = (4*((r+1)%3))%6;

  if (mod_step == 2 && r+1 < valColK.size()) {
    bool restart = true;
    while (restart) {
      restart = false;
      for (unsigned rr = r; rr <= r+1; ++rr) {
        for (unsigned c1 = 0; c1 < 4; ++c1) {
          if (rr+1 >= valColK.size() || valColK[rr+1][c1] == 0) continue;
          unsigned pos1 = (4*(rr+1) + c1)%6;
          if (pos1 == 0) {
            unsigned rk = rr+1;
            unsigned ck = c1 + 5;
            if (ck >= 4) {ck -= 4; rk += 1;}
            unsigned r2 = rk;
            unsigned c2 = ck+1;
            if (c2 >= 4) {c2 -= 4; r2 += 1;}
            if (r2 >= valColK.size()) continue;
            auto flag = constraintMC(rk, ck, 1, rr+1, c1, r2, c2, mat, line1, line2, valX, valK, valColX, valColK, valColSR);
            if (!flag.first) return;
            if (flag.second) {restart = true; break;}
          }
          if (pos1 != 5) {
            unsigned r2 = rr+1;
            unsigned c2 = c1 + 1;
            if (c2 >= 4) {c2 -= 4; r2 += 1;}
            if (r2 >= valColK.size()) continue;
            unsigned rk = r2 - 2;
            unsigned ck = c2 + 2;
            if (ck >= 4) {ck -= 4; rk += 1;}
            if (rk >= valColK.size() || valColK[rk][ck] == 5) continue;
            auto flag = constraintMC(rk, ck, 0, rr+1, c1, r2, c2, mat, line1, line2, valX, valK, valColX, valColK, valColSR);
            if (!flag.first) return;
            if (flag.second) {restart = true; break;}
          }
        }
        if (restart) break;
      }
    }
  }

  if (mod_step == 2) {
    next_state_key = inv_updateSR(T[step], global_bound-current_bound, state_key%n_states, state_key/n_states, valX[r]);
    for (auto f : next_state_key) {
      for (unsigned c = 0; c < 4; ++c) {
        valColX[r][c] = (f/mypow[c])%5;
      }
      auto m = mat;
      auto vX = valX;
      auto vK = valK;
      auto l1 = line1;
      auto l2 = line2;
      bool isvalid = true;
      for (unsigned c = 0; (c < 4) && isvalid; ++c) isvalid = updateColX(r, c, vX, vK, l1, l2, m, valColX, valColSR, valColK);
      if (isvalid) findBestTrail(f, T, global_bound, current_bound, step-1, m, l1, l2, vX, vK, valColX, valColK, valColSR);
    }
    for (unsigned c = 0; c < 4; ++c) {
      valColX[r][c] = 5;
    }
  }
  else {
    if (mod_step == 0) {

      if (r == 0) return findBestTrail(state_key, T, global_bound, current_bound, step-1, mat, line1, line2, valX, valK, valColX, valColK, valColSR);

      next_state_key = inv_updateKey192((colk + 3)%6, (r == 1) ? 2 : 4, T[step], global_bound-current_bound, state_key%n_states, state_key/n_states);

      //if (r == 1 && !next_state_key.empty()) return findBestTrail(next_state_key[0], T, global_bound, current_bound, step-1, mat, line1, line2, valX, valK, valColX, valColK, valColSR);

      //cout << "state_key: "; for (unsigned c = 0; c < 6; ++c) {cout << ((state_key/n_states)/mypow[c])%5 << " | ";} cout << endl;
      for (auto f : next_state_key) {
        //cout << "f: "; for (unsigned c = 0; c < 6; ++c) {cout << ((f/n_states)/mypow[c])%5 << " | ";} cout << endl;
        if (r > 1) for (unsigned c = 0; c < 2; ++c) valColK[r-1][c+2] = ((f/n_states)/mypow[(colk + c)%6])%5;
        for (unsigned c = 0; c < 2; ++c) valColK[r][c] = ((f/n_states)/mypow[(colk + 2 + c)%6])%5;
        //cout << r << ": ";
        //for (unsigned c = 0; c < 4; ++c) cout << (unsigned) valColK[r][c] << "| ";
        //cout << endl;
        auto m = mat;
        auto vX = valX;
        auto vK = valK;
        auto l1 = line1;
        auto l2 = line2;
        bool isvalid = true;
        for (unsigned c = 0; (c < 2) && isvalid; ++c) isvalid = updateColK(r-1, c+2, vX, vK, l1, l2, m, valColX, valColSR, valColK);
        for (unsigned c = 0; (c < 2) && isvalid; ++c) isvalid = updateColK(r, c, vX, vK, l1, l2, m, valColX, valColSR, valColK);
        if (isvalid) findBestTrail(f, T, global_bound, current_bound, step-1, m, l1, l2, vX, vK, valColX, valColK, valColSR);
      }
      for (unsigned c = 0; c < 2; ++c) valColK[r][c] = 5;
      for (unsigned c = 0; c < 2; ++c) valColK[r-1][c+2] = 5;
    }
    else {
      next_state_key = inv_updateMC_ARK(T[step], colk, global_bound-current_bound, state_key%n_states, state_key/n_states);
      // cout << "c: " << (unsigned) (c + dec_key) << endl;
      // cout << "nb sol: " << next_state_key.size() << endl;

      bool greater[4];
      for (unsigned c = 0; c < 4; ++c) greater[c] = (valColSR[r+1][c] + valColX[r+2][c] >= ((valColK[r+2][c] == 0) ? 6 : 5));

      for (auto f : next_state_key) {
        bool isvalid = true;
        auto m = mat;
        auto vX = valX;
        auto vK = valK;
        auto l1 = line1;
        auto l2 = line2;
        for (unsigned c = 0; (c < 4) && isvalid; ++c) {
          valColSR[r][c] = ((f%n_states)/mypow[c])%5;
          if (valColSR[r][c] + valColX[r+1][c] >= ((valColK[r+1][c] == 0) ? 6 : 5)) {
            for (unsigned cc = 0; cc < 4; ++cc) {
              if (greater[cc]) {
                isvalid = isvalid && updateStateVar(0, r+1, (c-cc+4)%4, c, m, l1, l2, vX, vK, valColX, valColSR, valColK);
              }
            }
          }
        }

        for (unsigned c = 0; (c < 4) && isvalid; ++c) isvalid = updateColSR(r, c, vX, vK, l1, l2, m, valColX, valColSR, valColK);
        if (isvalid) {
          unsigned cost = current_bound + valColX[r+1][0] + valColX[r+1][1] + valColX[r+1][2] + valColX[r+1][3];
          if (colk >= 2) cost += ((f/n_states)/mypow[5])%5;
          findBestTrail(f, T, global_bound, cost, step-1, m, l1, l2, vX, vK, valColX, valColK, valColSR);
        }
      }
      for (unsigned c = 0; c < 4; ++c) valColSR[r][c] = 5;
    }
  }
}

void updateBounds(vector<vector<uint8_t>> & T, unsigned Round) {

  static unsigned const n_states = 5*5*5*5;
  static unsigned const n_keys = 5*5*5*5*5*5;
  static unsigned const mypow[8] = {1, 5, 5*5, 5*5*5, 5*5*5*5, 5*5*5*5*5, 5*5*5*5*5*5, 5*5*5*5*5*5*5};

  auto mat = AES192eqs(Round);

  unsigned pivot = 0;
  for (unsigned r = 0; r <= 4*Round+2; ++r) {
    for (unsigned i = 0; i < 16; ++i) {
      if ((r%4 == 2 && r < 4*Round+2) || (r%4 == 1 && r < 4*Round+1)) continue;
      if (mat.setAsPivot(16*r + i, pivot)) ++pivot;
    }
  }

  mat = mat.extract(pivot);

  uint8_t my_min = 255;
  for (auto x : T.back()) {
    if (x < my_min) my_min = x;
  }

  cout << "my min: " << (unsigned) my_min << endl;

  bool flag = false;

  for (unsigned b = my_min; b < 255; ++b) {
    bool found = false;
    #pragma omp parallel for schedule(dynamic)
    for (unsigned x = 0; x < n_states*n_keys; ++x) {
      if (T.back()[x] == b) {
        vector<vector<uint8_t>> valX (Round, vector<uint8_t> (16,2));
        vector<vector<uint8_t>> valK (Round, vector<uint8_t> (16,2));
        vector<vector<uint8_t>> valColK (Round, vector<uint8_t> (4,5));
        vector<vector<uint8_t>> valColX (Round, vector<uint8_t> (4,5));
        vector<vector<uint8_t>> valColSR (Round, vector<uint8_t> (4,5));

        uint8_t bb = b;
        auto y = x;
        for (unsigned c = 0; c < 4; ++c) {valColX[Round-1][c] = y%5; y = y/5;}
        unsigned colk = (4*((Round-1)%3))%6;
        for (unsigned c = 0; c < 4; ++c) {valColK[Round-1][c] = (y/mypow[(colk+c)%6])%5;}
        for (unsigned c = 0; c < 2; ++c) {valColK[Round-2][2+c] = (y/mypow[(colk+c+4)%6])%5;}


        auto m = mat;
        auto line1 = 0u;
        auto line2 = mat.nblines;
        for (unsigned c = 0; c < 4; ++c) {
          if (c >= 2) updateColK(Round-2, c, valX, valK, line1, line2, m, valColX, valColSR, valColK);
          updateColK(Round-1, c, valX, valK, line1, line2, m, valColX, valColSR, valColK);
          updateColX(Round-1, c, valX, valK, line1, line2, m, valColX, valColSR, valColK);
        }

        {
          auto step = T.size()-2;
          auto next_state_key = inv_updateMC_ARK(T[step], colk, bb, x%n_states, x/n_states);
          // cout << "c: " << (unsigned) (c + dec_key) << endl;
          // cout << "nb sol: " << next_state_key.size() << endl;
          bool found_local = false;
          for (auto f : next_state_key) {
            bool isvalid = true;
            for (unsigned c = 0; c < 4; ++c) {
              valColSR[Round-2][c] = ((f%n_states)/mypow[c])%5;
              if (valColK[Round-1][c] == 0 && valColX[Round-1][c] > 5-valColSR[Round-2][c]) isvalid = false;
              if (valColK[Round-1][c] > 0 && valColSR[Round-2][c] > 0 && valColX[Round-1][c] > valColK[Round-1][c]-1 && valColX[Round-1][c] > 4-valColSR[Round-2][c]) isvalid = false;
            }
            if (!isvalid) {
              continue;
            }
            //continue;
            auto mm = m;
            auto vX = valX;
            auto vK = valK;
            auto l1 = line1;
            auto l2 = line2;

            for (unsigned c = 0; (c < 4) && isvalid; ++c) isvalid = updateColSR(Round-2, c, vX, vK, l1, l2, mm, valColX, valColSR, valColK);
            if (isvalid) {
              unsigned cost = valColX[Round-1][0] + valColX[Round-1][1] + valColX[Round-1][2] + valColX[Round-1][3];
              if (colk >= 2) cost += ((f/n_states)/mypow[5])%5;
              if (findBestTrail1(f, T, bb, cost, step-1, mm, l1, l2, vX, vK, valColX, valColK, valColSR)) {
                found_local = true;
                #pragma omp critical
                {
                  found = true;
                }
                break;
              }
            }
          }
          if (!found_local) T.back()[x] += 1;
        }
      }
    }
    if (flag) break;
    if (found) {
      if (!flag) cout << "my min: " << (unsigned) b << endl;
      flag = true;
      //break;
    }
  }
}

vector<vector<uint8_t>> computeDynProg(uint8_t const global_bound, unsigned const Round) {
  static unsigned const n_states = 5*5*5*5;
  static unsigned const n_keys = 5*5*5*5*5*5;
  static unsigned const mypow[8] = {1, 5, 5*5, 5*5*5, 5*5*5*5, 5*5*5*5*5, 5*5*5*5*5*5, 5*5*5*5*5*5*5};
  static unsigned const sizeT = n_states*n_keys;
  static vector<uint8_t> const count = initPop5();
  static vector<unsigned> const shiftRows ({0, 1, 2, 3, 7, 4, 5, 6, 10, 11, 8, 9, 13, 14, 15, 12});

  static vector<uint8_t> best ({0,0,0,1,4,5,10,13,18});


  // intialisation of T
  vector<uint8_t> T (sizeT, global_bound);
  for (unsigned k = 0; k < n_keys; ++k) {
    uint8_t sboxes = 0;
    //sboxes += (k/mypow[7])%5;
    //sboxes += (k/mypow[3])%5;
    for (unsigned x = 0; x < n_states; ++x) {
      uint8_t sboxes2 = sboxes + count[x];
      if (sboxes2 >= global_bound) continue;
      T[k*n_states + x] = sboxes2;
    }
  }
  T[0] = global_bound;

  vector<vector<uint8_t>> res;

  res.emplace_back(T);

  unsigned colk = 4;

  for (unsigned r = 1; r < Round; ++r) {
  	cout << "here: " << r << "/" << Round << endl;

    if (r != 1) {
      updateBounds(res, r);
      cout << "updated" << endl;
      T = res.back();
      updateSR(T, global_bound); res.emplace_back(T);
      for (unsigned c = 0; c < 4; ++c) {
        updateKey192Column((c + colk)%6, T, global_bound);
      }
      res.emplace_back(T);
      if (colk >= 2) {
        for (unsigned k = 0; k < n_keys; ++k) {
          uint8_t sboxes = (k/mypow[5])%5;
          for (unsigned x = 0; x < n_states; ++x) {
            auto & src = T[k*n_states + x];
            if (src < global_bound) src += sboxes;
          }
        }
      }
    }
    else {
      for (unsigned c = 2; c < 4; ++c) {
        updateKey192Column((c + colk)%6, T, global_bound);
      }
      res.emplace_back(T);
      if (colk >= 2) {
        for (unsigned k = 0; k < n_keys; ++k) {
          uint8_t sboxes = (k/mypow[5])%5;
          for (unsigned x = 0; x < n_states; ++x) {
            auto & src = T[k*n_states + x];
            if (src < global_bound) src += sboxes;
          }
        }
      }
    }
    for (unsigned c = 0; c < 4; ++c) {
      updateMC_ARK(T, c, (c + colk)%6, global_bound);
    }
    colk = (colk + 4)%6;
    /*if (r+1 < Round && r+1 < best.size()) {
      for (auto & x : T) if (x < best[r+1]) x = best[r+1];
    }*/
    res.emplace_back(move(T));
  }
  
  cout << "here" << endl;
  
  return res;
}

int main(int argc, char const *argv[]) {
  unsigned Round = stoi(argv[1]);

  static unsigned const n_states = 5*5*5*5;
  static unsigned const n_keys = 5*5*5*5*5*5;
  static unsigned const mypow[8] = {1, 5, 5*5, 5*5*5, 5*5*5*5, 5*5*5*5*5, 5*5*5*5*5*5, 5*5*5*5*5*5*5};

  {
    auto mat = AES192eqs(Round);

    unsigned pivot = 0;
    for (unsigned r = 0; r <= 4*Round+2; ++r) {
      for (unsigned i = 0; i < 16; ++i) {
        if ((r%4 == 2 && r < 4*Round+2) || (r%4 == 1 && r < 4*Round+1)) continue;
        if (mat.setAsPivot(16*r + i, pivot)) ++pivot;
      }
    }

    mat = mat.extract(pivot);
    cout << mat << endl;
    //getchar();
    unsigned bound = 0;
    for (unsigned r = 0; r < Round; ++r) {
      if (r%4 == 0) bound += 1;
      else if (r%4 == 3) bound += 16;
      else bound += 4;
    }
    uint8_t global_bound = min(254u , bound);
    //global_bound = 14;
    cout << "global_bound: " << (unsigned) global_bound << endl;
    //getchar();


    {

      static vector<uint8_t> const count = initPop5();

      auto T = computeDynProg(global_bound, Round);
      uint8_t my_min = global_bound;
      for (auto x : T.back()) {
        if (x < my_min) my_min = x;
      }
      cout << "min bound: " << (unsigned) my_min << endl;
      my_min = 24;

      vector<uint8_t> myvec (n_states*n_keys,0);


      unsigned dec_Round = 3;
      for (unsigned b = my_min; b < global_bound; ++b) {
        /*Round -= dec_Round;
        auto mat2 = AES192eqs(Round);

        unsigned pivot = 0;
        for (unsigned r = 0; r <= 4*Round+2; ++r) {
          for (unsigned i = 0; i < 16; ++i) {
            if ((r%4 == 2 && r < 4*Round+2) || (r%4 == 1 && r < 4*Round+1)) continue;
            if (mat2.setAsPivot(16*r + i, pivot)) ++pivot;
          }
        }

        mat2 = mat2.extract(pivot);

        unsigned cpt_tmp = 0;
        for (unsigned x = 0; x < n_states*n_keys; ++x) {
          if (T[T.size()-1 - 3*dec_Round][x] <= (b*Round)/(Round+dec_Round)) cpt_tmp += 1;
        }
        unsigned cpt_tmp2 = 0;

        #pragma omp parallel for schedule(dynamic)
        for (unsigned x = 0; x < n_states*n_keys; ++x) {
          while (T[T.size()-1 -3*dec_Round][x] <= (b*Round)/(Round+dec_Round)) {
            vector<vector<uint8_t>> valX (Round, vector<uint8_t> (16,2));
            vector<vector<uint8_t>> valK (Round, vector<uint8_t> (16,2));
            vector<vector<uint8_t>> valColK (Round, vector<uint8_t> (4,5));
            vector<vector<uint8_t>> valColX (Round, vector<uint8_t> (4,5));
            vector<vector<uint8_t>> valColSR (Round, vector<uint8_t> (4,5));

            uint8_t bb = T[T.size()-1 - 3*dec_Round][x];
            auto y = x;
            for (unsigned c = 0; c < 4; ++c) {valColX[Round-1][c] = y%5; y = y/5;}
            unsigned colk = (4*((Round-1)%3))%6;
            for (unsigned c = 0; c < 4; ++c) {valColK[Round-1][c] = (y/mypow[(colk+c)%6])%5;}
            for (unsigned c = 0; c < 2; ++c) {valColK[Round-2][2+c] = (y/mypow[(colk+c+4)%6])%5;}


            auto m = mat2;
            auto l1 = 0u;
            auto l2 = mat2.nblines;
            for (unsigned c = 0; c < 4; ++c) {
              if (c >= 2) updateColK(Round-2, c, valX, valK, l1, l2, m, valColX, valColSR, valColK);
              updateColK(Round-1, c, valX, valK, l1, l2, m, valColX, valColSR, valColK);
              updateColX(Round-1, c, valX, valK, l1, l2, m, valColX, valColSR, valColK);
            }
            //cout << "here: " << (unsigned) count[x % (5*5*5*5)] << endl;
            //findBestTrail(x, T, bb, 0, T.size()-2, mat, 0, mat.nblines, valX, valK, valColX, valColK, valColSR);
            if (findBestTrail1(x, T, bb, 0, T.size()-2 - 3*dec_Round, m, l1, l2, valX, valK, valColX, valColK, valColSR)) {
              #pragma omp critical
              {
                cout << "\r" << ++cpt_tmp2 << "/" << cpt_tmp << " " << flush;
              }
              break;
            }
            else {
              T[T.size()-1 -3*dec_Round][x] += 1;
              if (T[T.size()-1 -3*dec_Round][x] > (b*Round)/(Round+dec_Round)) {
                #pragma omp critical
                {
                  cout << "\r" << ++cpt_tmp2 << "/" << cpt_tmp << " " << flush;
                }
              }
            }
          }
        }
        cout << "done" << endl;
        Round += dec_Round;*/
        map<unsigned, vector<unsigned>> mymap;
        unsigned totry = 0, cpt = 0, eliminated = 0;
        for (unsigned x = 0; x < n_states*n_keys; ++x) {
          if (T.back()[x] <= b) {totry += 1; mymap[T.back()[x]].emplace_back(x);}
        }


        #pragma omp parallel for schedule(dynamic)
        for (unsigned x = 0; x < n_states*n_keys; ++x) {
          if (T.back()[x] <= b) {
            vector<vector<uint8_t>> valX (Round, vector<uint8_t> (16,2));
            vector<vector<uint8_t>> valK (Round, vector<uint8_t> (16,2));
            vector<vector<uint8_t>> valColK (Round, vector<uint8_t> (4,5));
            vector<vector<uint8_t>> valColX (Round, vector<uint8_t> (4,5));
            vector<vector<uint8_t>> valColSR (Round, vector<uint8_t> (4,5));

            uint8_t bb = b;
            auto y = x;
            for (unsigned c = 0; c < 4; ++c) {valColX[Round-1][c] = y%5; y = y/5;}
            unsigned colk = (4*((Round-1)%3))%6;
            for (unsigned c = 0; c < 4; ++c) {valColK[Round-1][c] = (y/mypow[(colk+c)%6])%5;}
            for (unsigned c = 0; c < 2; ++c) {valColK[Round-2][2+c] = (y/mypow[(colk+c+4)%6])%5;}


            auto m = mat;
            auto line1 = 0u;
            auto line2 = mat.nblines;
            for (unsigned c = 0; c < 4; ++c) {
              if (c >= 2) updateColK(Round-2, c, valX, valK, line1, line2, m, valColX, valColSR, valColK);
              updateColK(Round-1, c, valX, valK, line1, line2, m, valColX, valColSR, valColK);
              updateColX(Round-1, c, valX, valK, line1, line2, m, valColX, valColSR, valColK);
            }

            {
              auto step = T.size()-2;
              auto next_state_key = inv_updateMC_ARK(T[step], colk, bb, x%n_states, x/n_states);
              // cout << "c: " << (unsigned) (c + dec_key) << endl;
              // cout << "nb sol: " << next_state_key.size() << endl;

              for (auto f : next_state_key) {
                bool isvalid = true;
                for (unsigned c = 0; c < 4; ++c) {
                  valColSR[Round-2][c] = ((f%n_states)/mypow[c])%5;
                  if (valColK[Round-1][c] == 0 && valColX[Round-1][c] > 5-valColSR[Round-2][c]) isvalid = false;
                  if (valColK[Round-1][c] > 0 && valColSR[Round-2][c] > 0 && valColX[Round-1][c] > valColK[Round-1][c]-1 && valColX[Round-1][c] > 4-valColSR[Round-2][c]) isvalid = false;
                }
                if (!isvalid) {
                  #pragma omp critical
                  {
                    ++eliminated;
                  }
                  continue;
                }
                //continue;
                auto mm = m;
                auto vX = valX;
                auto vK = valK;
                auto l1 = line1;
                auto l2 = line2;

                for (unsigned c = 0; (c < 4) && isvalid; ++c) isvalid = updateColSR(Round-2, c, vX, vK, l1, l2, mm, valColX, valColSR, valColK);
                if (isvalid) {
                  unsigned cost = valColX[Round-1][0] + valColX[Round-1][1] + valColX[Round-1][2] + valColX[Round-1][3];
                  if (colk >= 2) cost += ((f/n_states)/mypow[5])%5;
                  findBestTrail(f, T, bb, cost, step-1, mm, l1, l2, vX, vK, valColX, valColK, valColSR);
                }
              }
            }
            #pragma omp critical
            {
              cout << "\r" << ++cpt << "/" << totry << " (" << eliminated << ") " << flush;
            }

            //cout << "here: " << (unsigned) count[x % (5*5*5*5)] << endl;
            //findBestTrail(x, T, bb, 0, T.size()-2, mat, 0, mat.nblines, valX, valK, valColX, valColK, valColSR);
            //findBestTrail(x, T, bb, 0, T.size()-2, m, l1, l2, valX, valK, valColX, valColK, valColSR);
          }
        }
        cout << "cpt: " << cpt << endl;
        cout << "b : " << b << " - done" << endl;
        if (flag_solution_found) break;
        //getchar();
      }
    }


    // //global_bound = 17;
    // vector<pair<uint16_t, uint16_t>> mykey (Round, make_pair(0,0));
    // auto cost = searchBound5(global_bound, mykey, Round, 0);
    //cout << "cost : " << (unsigned) global_bound << endl;
    return 0;
  }


  return 0;
}
