#include <set>
#include <map>

#include "Matrix.hpp"

using namespace std;

void printVar(ostream &flux, int x) {
  if (x >= 0) {
    int i = x%16;
    int v = (x/16)%4;
    int r = (x/16)/4;
    if (v == 0) flux << "Z";
    else if (v == 1) flux << "K";
    else if (v == 2) flux << "X";
    else flux << "Y";
    flux << r << "[" << i << "]";
  }
  else {
    flux << "S(";
    printVar(flux, -x);
    flux << ")";
  }
}

Matrix::Matrix(vector<vector<pair<GFElement, int>>> const & sys) {
  set<int> variables;
  for (auto const & v : sys) {
    for (auto const & p : v) variables.emplace(p.second);
  }
  unsigned n = sys.size();
  unsigned m = 0;
  map<int, unsigned> mapvars;
  for (auto x : variables) {
    mapvars[x] = m++;
  }
  vector<int> v (variables.begin(), variables.end());
  vector<vector<GFElement>> mat (n, vector<GFElement> (m, 0));
  for (unsigned i = 0; i < n; ++i) {
    for (auto const & p : sys[i]) mat[i][mapvars[p.second]] = p.first;
  }

  // for (unsigned i = 0; i < n; ++i) {
  //   for (unsigned j = 0; j < m; ++j) {
  //     if (mat[i][j] == 0) continue;
  //     cout << " + " << mat[i][j] << ".";
  //     printVar(cout, v[j]);
  //   }
  //   cout << endl;
  // }
  // cout << " -------------- " << endl;
  // getchar();

  for (unsigned p = 0; p < n; ++p) {
    unsigned c = p;
    while (c < m && mat[p][c] == 0) ++c;
    if (c == m) swap(mat[p--], mat[--n]);
    else {
      if (c != p) {
        for (unsigned i = 0; i < n; ++i) swap(mat[i][p], mat[i][c]);
        swap(v[p], v[c]);
      }
      auto const coefInv = mat[p][p].getInverse();
      for (c = p; c < m; ++c) mat[p][c] *= coefInv;
      for (unsigned l = 0; l < n; ++l) {
        if (mat[l][p] == 0 || l == p) continue;
        auto const coef = mat[l][p];
        for (c = p; c < m; ++c) mat[l][c] -= coef*mat[p][c];
      }
    }
  }

  // for (unsigned i = 0; i < n; ++i) {
  //   for (unsigned j = 0; j < m; ++j) {
  //     if (mat[i][j] == 0) continue;
  //     cout << " + " << mat[i][j] << ".";
  //     printVar(cout, v[j]);
  //   }
  //   cout << endl;
  // }
  // cout << " -------------- " << endl;

  front = vector<int> (n);
  for (unsigned i = 0; i < n; ++i) front[i] = v[i];
  columns = vector<int> (m-n);
  for (unsigned i = n; i < m; ++i) columns[i-n] = v[i];
  space = vector<GFElement> (n*(m-n));
  for (unsigned i = 0; i < n; ++i) {
    for (unsigned j = n; j < m; ++j) space[i*(m-n) + j-n] = mat[i][j];
  }
  nbcols = m-n;
  nblines = n;
  lines = vector<GFElement *> (n);
  for (unsigned i = 0; i < n; ++i) lines[i] = space.data() + i*nbcols;
}

Matrix Matrix::extract(unsigned l) const {
  vector<vector<pair<GFElement, int>>> sys;
  for (; l < nblines; ++l) {
    vector<pair<GFElement, int>> v;
    v.emplace_back(1, front[l]);
    for (unsigned i = 0; i < nbcols; ++i) {
      if ((*this)(l,i) != 0) v.emplace_back((*this)(l,i), columns[i]);
    }
    sys.emplace_back(move(v));
  }
  return Matrix(sys);
}

void Matrix::swapLineColumn(unsigned l, unsigned c) {
  auto coef = (*this)(l,c).getInverse();
  for (unsigned i = 0; i < nbcols; ++i) (*this)(l,i) *= coef;
  (*this)(l,c) = 0;
  for (unsigned i = 0; i < nblines; ++i) {
    if ((*this)(i,c) == 0) continue;
    auto coef2 = (*this)(i,c);
    for (unsigned j = 0; j < nbcols; ++j) (*this)(i,j) += coef2*(*this)(l,j);
    (*this)(i,c) = coef*coef2;
  }
  (*this)(l,c) = coef;
  swap(front[l], columns[c]);
}

unsigned Matrix::checkZ(unsigned line1, vector<vector<uint8_t>> const & valX, vector<vector<uint8_t>> const & valK) {
  for (unsigned l = line1; l-- != 0; ) {
    unsigned c = 0;
    while (c < nbcols) {
      if ((*this)(l,c) != 0) {
        auto u = abs(columns[c]);
        int i = u%16;
        int v = (u/16)%4;
        int r = (u/16)/4;
        if (v == 2 && valX[r][i] == 1) break;
        else if (v == 1 && valK[r][i/4] == 1) break;
      }
      ++c;
    }
    if (c == nbcols) return l;
  }
  return line1;
}

// unsigned Matrix::checkZ(unsigned line1, vector<vector<uint8_t>> const & valX, vector<vector<uint8_t>> const & valK) {
//   unsigned opt = line1;
//   unsigned nopt = nbcols;
//   for (unsigned l = 0; l < line1; ++l) {
//     unsigned c = 0;
//     unsigned tmp = 0;
//     while (c < nbcols) {
//       if ((*this)(l,c) != 0) {
//         auto u = abs(columns[c]);
//         int i = u%16;
//         int v = (u/16)%4;
//         int r = (u/16)/4;
//         if (v == 2) {
//           if (valX[r][i] == 1) break;
//           else if (valX[r][i] == 2) tmp +=1;
//         }
//         else {
//           if (valK[r][i/4] == 1) break;
//           else if (valK[r][i/4] == 2) tmp += 1;
//         }
//       }
//       ++c;
//     }
//     if (c == nbcols && nopt > tmp) {
//       if (tmp == 1) cout << "un" << endl;
//       opt = l; nopt = tmp;
//     }
//   }
//   return opt;
// }

unsigned Matrix::checkZ(double * X) {
  for (unsigned l = 0; l < nblines; ++l) {
    auto const u = abs(front[l]);
    if (X[u] > 0.5) {
      unsigned c = 0;
      while (c < nbcols && ((*this)(l,c) == 0 || X[abs(columns[c])] > 0.5)) ++c;
      if (c < nbcols) swapLineColumn(l,c);
    }
  }
  for (unsigned l = 0; l < nblines; ++l) {
    auto const u = abs(front[l]);
    if (X[u] > 0.5 && (u/16)%4 == 1) {
      unsigned c = 0;
      while (c < nbcols && ((*this)(l,c) == 0 || (abs(columns[c])/16)%4 == 1)) ++c;
      if (c < nbcols) swapLineColumn(l,c);
    }
  }
  for (unsigned l = 0; l < nblines; ++l) {
    if (X[abs(front[l])] > 0.5 || (abs(front[l])/16)%4 != 1) continue;
    unsigned c = 0;
    while (c < nbcols && ((*this)(l,c) == 0 || (X[abs(columns[c])] > 0.5 && (abs(columns[c])/16)%4 == 1))) ++c;
    if (c == nbcols) return l;
  }
  for (unsigned l = 0; l < nblines; ++l) {
    if (X[abs(front[l])] > 0.5 /*|| (abs(front[l])/16)%4 == 1*/) continue;
    unsigned c = 0;
    while (c < nbcols && ((*this)(l,c) == 0 || X[abs(columns[c])] > 0.5)) ++c;
    if (c == nbcols) return l;
  }
  return nblines;
}

void Matrix::swapLines(unsigned l1, unsigned l2) {
  swap(front[l1], front[l2]);
  swap(lines[l1], lines[l2]);
}

bool Matrix::setAsPivot(int x, unsigned start) {
  for (unsigned l = start; l < nblines; ++l) {
    if (front[l] == x) {
      swap(front[l], front[start]);
      swap(lines[l], lines[start]);
      return true;
    }
  }
  for (unsigned c = 0; c < nbcols; ++c) {
    if (columns[c] == x) {
      for (unsigned l = start; l < nblines; ++l) {
        if ((*this)(l,c) != 0) {
          swapLineColumn(l,c);
          swap(front[l], front[start]);
          swap(lines[l], lines[start]);
          return true;
        }
      }
      return false;
    }
  }
  return false;
}

bool Matrix::setAsPivot(int x, unsigned start, unsigned end) {
  for (unsigned l = start; l < end; ++l) {
    if (front[l] == x) {
      swap(front[l], front[start]);
      swap(lines[l], lines[start]);
      return true;
    }
  }
  for (unsigned c = 0; c < nbcols; ++c) {
    if (columns[c] == x) {
      for (unsigned l = start; l < end; ++l) {
        if ((*this)(l,c) != 0) {
          swapLineColumn(l,c);
          swap(front[l], front[start]);
          swap(lines[l], lines[start]);
          return true;
        }
      }
      return false;
    }
  }
  return false;
}

bool Matrix::setColumnAsPivot(unsigned c, unsigned start, unsigned end) {
  for (unsigned l = start; l < end; ++l) {
    if ((*this)(l,c) != 0) {
      swapLineColumn(l,c);
      swap(front[l], front[start]);
      swap(lines[l], lines[start]);
      return true;
    }
  }
  return false;
}

void Matrix::eraseColumn(unsigned c) {
  nbcols -= 1;
  for (unsigned l = 0; l < nblines; ++l) {
    (*this)(l,c) = (*this)(l,nbcols);
  }
  columns[c] = columns[nbcols];
  columns.pop_back();
}

bool Matrix::isLinear(int x, unsigned start) {
  if (setAsPivot(x, start)) {
    return !setAsPivot(-x, start+1);
  }
  else return true;
}

vector<int> Matrix::checkK(double * X, double * SK, map<unsigned, unsigned> const & mapSK) {
  unsigned start = 0;
  vector<int> toprocess1;
  vector<int> toprocess2;
  for (auto x : front) {
    if (x >= 0 && X[x] < 0.5) {
      auto it = mapSK.find(x);
      if (it == mapSK.end()) { // unknown state variable
        if (setAsPivot(x, start)) ++start;
      }
      else {
        if (SK[it->second] < 0.5) toprocess2.emplace_back(x);
        else toprocess1.emplace_back(x);
      }
    }
    else if (x < 0 && X[-x] > 0.5 && SK[mapSK.at(-x)] < 0.5) toprocess1.emplace_back(x);
  }
  for (auto x : columns) {
    if (x >= 0 && X[x] < 0.5) {
      auto it = mapSK.find(x);
      if (it == mapSK.end()) { // unknown state variable
        if (setAsPivot(x, start)) ++start;
      }
      else {
        if (SK[it->second] < 0.5) toprocess2.emplace_back(x);
        else toprocess1.emplace_back(x);
      }
    }
    else if (x < 0 && X[-x] > 0.5 && SK[mapSK.at(-x)] < 0.5) toprocess1.emplace_back(x);
  }

  unsigned i = 0, n = toprocess2.size();
  while (i < n) {
    if (setAsPivot(toprocess2[i], start)) {
      if (setAsPivot(-toprocess2[i], start+1)) { // does not linearly appear
        ++i;
      }
      else { // linearly appear
        ++start;
        toprocess2[i] = toprocess2[--n];
        i = 0;
      }
    }
    else { // toprocess_v[i] does not appear in the system
      if (setAsPivot(-toprocess2[i], start)) ++start;
      toprocess2[i] = toprocess2[--n];
      i = 0;
    }
  }

  toprocess2.resize(n);
  for (auto x : toprocess1) if (setAsPivot(x, start)) toprocess2.emplace_back(abs(x));
  return toprocess2;
}

vector<int> Matrix::checkZ2(double * X, unsigned nX, set<unsigned> const & setSK) {
  unsigned start = 0;
  vector<int> toprocess2;
  for (auto x : front) {
    if (x >= 0 && X[x] < 0.5) {
      if (setSK.count(x) == 0) { // unknown state variable
        if (setAsPivot(x, start)) ++start;
      }
      else toprocess2.emplace_back(x);
    }
  }
  for (auto x : columns) {
    if (x >= 0 && X[x] < 0.5) {
      if (setSK.count(x) == 0) { // unknown state variable
        if (setAsPivot(x, start)) ++start;
      }
      else {
        toprocess2.emplace_back(x);
      }
    }
  }

  unsigned i = 0, n = toprocess2.size();
  while (i < n) {
    if (setAsPivot(toprocess2[i], start)) {
      if (setAsPivot(-toprocess2[i], start+1)) { // does not linearly appear
        ++i;
      }
      else { // linearly appear
        ++start;
        toprocess2[i] = toprocess2[--n];
        i = 0;
      }
    }
    else { // toprocess_v[i] does not appear in the system
      if (setAsPivot(-toprocess2[i], start)) ++start;
      toprocess2[i] = toprocess2[--n];
      i = 0;
    }
  }
  toprocess2.resize(n);
  if (n == 0) return toprocess2;

  for (unsigned j = 0; j < nX; ++j) {
    if (X[j] < 0.5) continue;
    auto back_start = start;
    auto res = toprocess2;
    if (setAsPivot(j, start) || (setSK.count(j) != 0 && setAsPivot(-j, start))) {
      i = 0; n = toprocess2.size();
      ++start;
      if (setSK.count(j) != 0 && setAsPivot(-j, start)) ++start;
      while (i < n) {
        if (setAsPivot(toprocess2[i], start)) {
          if (setAsPivot(-toprocess2[i], start+1)) { // does not linearly appear
            ++i;
          }
          else { // linearly appear
            ++start;
            toprocess2[i] = toprocess2[--n];
            i = 0;
          }
        }
        else { // toprocess_v[i] does not appear in the system
          if (setAsPivot(-toprocess2[i], start)) ++start;
          toprocess2[i] = toprocess2[--n];
          i = 0;
        }
      }
      toprocess2.resize(n);
    }
    if (toprocess2.empty()) {
      toprocess2 = move(res);
      start = back_start;
    }
    else X[j] = 0.0;
  }
  return toprocess2;
}

void Matrix::printLine(unsigned i) {
  printVar(cout, front[i]);
  for (unsigned j = 0; j < nbcols; ++j) {
    if ((*this)(i,j) != 0) {
      cout << " + ";
      if ((*this)(i,j) != 1) cout << (*this)(i,j) << ".";
      //flux << mat.columns[j];
      printVar(cout, columns[j]);
    }
  }
  cout << endl;
}

ostream& operator<<( ostream &flux, Matrix const& mat)
{
	for (unsigned i = 0; i < mat.nblines; ++i) {
    //flux << mat.front[i];
    printVar(flux, mat.front[i]);
    for (unsigned j = 0; j < mat.nbcols; ++j) {
      if (mat(i,j) != 0) {
        flux << " + ";
        if (mat(i,j) != 1) flux << mat(i,j) << ".";
        //flux << mat.columns[j];
        printVar(flux, mat.columns[j]);
      }
    }
    flux << endl;
  }
	return flux;
}
