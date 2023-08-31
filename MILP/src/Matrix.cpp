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

unsigned Matrix::checkZ(double * X) {
  for (unsigned l = 0; l < nblines; ++l) {
    auto const u = abs(front[l]);
    if (X[u] < 0.5) {
      unsigned c = 0;
      while (c < nbcols && ((*this)(l,c) == 0 || X[abs(columns[c])] < 0.5)) ++c;
      if (c < nbcols) swapLineColumn(l,c);
    }
  }

  for (unsigned l = 0; l < nblines; ++l) {
    if (X[abs(front[l])] < 0.5) continue;
    unsigned c = 0;
    while (c < nbcols && ((*this)(l,c) == 0 || X[abs(columns[c])] < 0.5)) ++c;
    if (c == nbcols) return l;
  }

  return nblines;
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

bool Matrix::isLinear(int x, unsigned start) {
  if (setAsPivot(x, start)) {
    return !setAsPivot(-x, start+1);
  }
  else return true;
}

bool isKeySboxAES128(int x) { // 16*(4r + 1) + 4*i + 3
  if (x < 0) x = -x;
  if ((x%4) != 3) return false;
  return ((x/16)%4 == 1);
}

bool isKeySboxAES192(int x) { // 16*(4r + 1) + 4*i + 3
  if (x < 0) x = -x;
  if ((x/16)%4 != 1) return false;
  int r = x/64;
  int c = x%4;
  return ((c+r)%6 == 5);
}

vector<int> Matrix::checkK128_256(double * X) {
  unsigned start = 0;
  vector<int> toprocess1;
  vector<int> toprocess2;
  for (auto x : front) {
    if (X[abs(x)] > 0.5) {
      if (!isKeySboxAES128(x)) { // unknown state variable
        if (setAsPivot(x, start)) ++start;
      }
      else if (x > 0) toprocess2.emplace_back(x);
    }
  }
  for (auto x : columns) {
    if (X[abs(x)] > 0.5) {
      if (!isKeySboxAES128(x)) { // unknown state variable
        if (setAsPivot(x, start)) ++start;
      }
      else if (x > 0) toprocess2.emplace_back(x);
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

  // if (n == 0) {
  //   for (unsigned i = 0; i < this->nblines; ++i) {
  //     if (((abs(front[i])/16)%4 != 1)) continue;
  //     //flux << mat.front[i];
  //     bool flag = false;
  //     if (X[abs(front[i])] > 0.5) {
  //       printVar(cout, this->front[i]);
  //       flag = true;
  //     }
  //     for (unsigned j = 0; j < nbcols; ++j) {
  //       if ((*this)(i,j) != 0 && X[abs(columns[j])] > 0.5) {
  //         if (flag) cout << " + ";
  //         else flag = true;
  //         if ((*this)(i,j) != 1) cout << (*this)(i,j) << ".";
  //         //flux << mat.columns[j];
  //         printVar(cout, columns[j]);
  //       }
  //     }
  //     if (flag) cout << endl;
  //   }
  // }

  return toprocess2;
}

vector<int> Matrix::checkK192(double * X) {
  unsigned start = 0;
  vector<int> toprocess1;
  vector<int> toprocess2;
  for (auto x : front) {
    if (X[abs(x)] > 0.5) {
      if (!isKeySboxAES192(x)) { // unknown state variable
        if (setAsPivot(x, start)) ++start;
      }
      else if (x > 0) toprocess2.emplace_back(x);
    }
  }
  for (auto x : columns) {
    if (X[abs(x)] > 0.5) {
      if (!isKeySboxAES192(x)) { // unknown state variable
        if (setAsPivot(x, start)) ++start;
      }
      else if (x > 0) toprocess2.emplace_back(x);
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

  // if (n == 0) {
  //   for (unsigned i = 0; i < this->nblines; ++i) {
  //     if (((abs(front[i])/16)%4 != 1)) continue;
  //     //flux << mat.front[i];
  //     bool flag = false;
  //     if (X[abs(front[i])] > 0.5) {
  //       printVar(cout, this->front[i]);
  //       flag = true;
  //     }
  //     for (unsigned j = 0; j < nbcols; ++j) {
  //       if ((*this)(i,j) != 0 && X[abs(columns[j])] > 0.5) {
  //         if (flag) cout << " + ";
  //         else flag = true;
  //         if ((*this)(i,j) != 1) cout << (*this)(i,j) << ".";
  //         //flux << mat.columns[j];
  //         printVar(cout, columns[j]);
  //       }
  //     }
  //     if (flag) cout << endl;
  //   }
  // }

  return toprocess2;
}

/*
vector<int> Matrix::checkK128_256(double * X) {
  unsigned start = 0;
  vector<int> toprocess1;
  vector<int> toprocess2;
  for (auto x : front) {
    if (x >= 0 && X[x] > 0.5) {
      if (!isKeySboxAES128(x)) { // unknown state variable
        if (setAsPivot(x, start)) ++start;
      }
      else toprocess2.emplace_back(x);
    }
  }
  for (auto x : columns) {
    if (x >= 0 && X[x] > 0.5) {
      if (!isKeySboxAES128(x)) { // unknown state variable
        if (setAsPivot(x, start)) ++start;
      }
      else toprocess2.emplace_back(x);
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
  return toprocess2;
}
*/

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

  // cout << " in ------------ " << endl;
  // for (unsigned i = start; i < nblines; ++i) {
  //   //flux << mat.front[i];
  //   printVar(cout, front[i]);
  //   for (unsigned j = 0; j < nbcols; ++j) {
  //     if ((*this)(i,j) != 0) {
  //       cout << " + ";
  //       if ((*this)(i,j) != 1) cout << (*this)(i,j) << ".";
  //       //flux << mat.columns[j];
  //       printVar(cout, columns[j]);
  //     }
  //   }
  //   cout << endl;
  // }
  // cout << " ------------ " << endl;

  //auto backstart = start;

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
  // for (unsigned i = backstart; i < start; ++i) {
  //   //flux << mat.front[i];
  //   printVar(cout, front[i]);
  //   for (unsigned j = 0; j < nbcols; ++j) {
  //     if ((*this)(i,j) != 0) {
  //       cout << " + ";
  //       if ((*this)(i,j) != 1) cout << (*this)(i,j) << ".";
  //       //flux << mat.columns[j];
  //       printVar(cout, columns[j]);
  //     }
  //   }
  //   cout << endl;
  // }
  //
  // cout << " out ------------ " << endl;
  // for (unsigned i = start; i < nblines; ++i) {
  //   //flux << mat.front[i];
  //   printVar(cout, front[i]);
  //   for (unsigned j = 0; j < nbcols; ++j) {
  //     if ((*this)(i,j) != 0) {
  //       cout << " + ";
  //       if ((*this)(i,j) != 1) cout << (*this)(i,j) << ".";
  //       //flux << mat.columns[j];
  //       printVar(cout, columns[j]);
  //     }
  //   }
  //   cout << endl;
  // }
  // cout << " ------------ " << endl;

  toprocess2.resize(n);
  for (auto x : toprocess1) if (setAsPivot(x, start)) toprocess2.emplace_back(abs(x));


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
