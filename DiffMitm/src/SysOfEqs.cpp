#include "SysOfEqs.hpp"

using namespace std;

// system of equations of AES-192
Matrix AES192eqs(unsigned R) {
  GFElement mc [4] = {2, 3, 1, 1};
  vector<unsigned> shiftRows ({0, 1, 2, 3, 7, 4, 5, 6, 10, 11, 8, 9, 13, 14, 15, 12});
  vector<unsigned> shiftRowsInv (16);
  for (unsigned i = 0; i < 16; ++i) shiftRowsInv[shiftRows[i]] = i;

  vector<vector<pair<GFElement, int>>> sys;
  // P + K0 = X0
  for (unsigned i = 0; i < 16; ++i) {
    vector<pair<GFElement, int>> v;
    v.emplace_back(1, i);
    v.emplace_back(1, i+16);
    v.emplace_back(1, i+32);
    sys.emplace_back(move(v));
  }
  for (unsigned r = 0; r < R-1; ++r) {
    for (unsigned c = 0; c < 4; ++c) {
      for (unsigned i = 0; i < 4; ++i) {
        vector<pair<GFElement, int>> v;
        for (unsigned j = 0; j < 4; ++j) v.emplace_back(mc[(j + 4 - i)%4], (4*r+3)*16 + shiftRowsInv[4*c + j]);
        v.emplace_back(1, (4*r+4)*16 + 4*i + c);
        sys.emplace_back(move(v));
      }
    }
  }
  for (unsigned i = 0; i < 16; ++i) {
    vector<pair<GFElement, int>> v;
    v.emplace_back(1, (4*R-1)*16 + i);
    v.emplace_back(1, (4*R)*16 + shiftRows[i]);
    sys.emplace_back(move(v));
  }
  for (unsigned r = 0; r < R; ++r) {
    for (unsigned i = 0; i < 16; ++i) {
      vector<pair<GFElement, int>> v;
      v.emplace_back(1, (4*r+4)*16 + i);
      v.emplace_back(1, (4*r+5)*16 + i);
      v.emplace_back(1, (4*r+6)*16 + i);
      sys.emplace_back(move(v));
    }
  }
    for (unsigned r = 1; r < R; ++r) {
        for (unsigned i = 0; i < 16; ++i) {
          if (4*r + (i%4) < 6) continue;
          if ((4*r + (i%4))%6 == 0) {
            if (i%4 == 0) {
              vector<pair<GFElement, int>> v;
              v.emplace_back(1, 16*(4*r + 1) + i);
              v.emplace_back(1, -(16*(4*(r-1) + 1) + ((i+7)%16)));
              v.emplace_back(1, 16*(4*(r-2) + 1) + (i+2));
              sys.emplace_back(move(v));
            }
            else {
              vector<pair<GFElement, int>> v;
              v.emplace_back(1, 16*(4*r + 1) + i);
              v.emplace_back(1, -(16*(4*r + 1) + ((i+3)%16)));
              v.emplace_back(1, 16*(4*(r-1) + 1) + (i-2));
              sys.emplace_back(move(v));
            }
          }
          else {
            vector<pair<GFElement, int>> v;
            v.emplace_back(1, 16*(4*r + 1) + i);
            if (i%4 == 0) v.emplace_back(1, 16*(4*(r-1) + 1) + i + 3);
            else v.emplace_back(1, 16*(4*r + 1) + i - 1);
            if (i%4 < 2) v.emplace_back(1, 16*(4*(r-2) + 1) + (i+2));
            else v.emplace_back(1, 16*(4*(r-1) + 1) + (i-2));
            sys.emplace_back(move(v));
          }
        }
    }
    return Matrix(sys);
}

 /* for (unsigned r = 1; r <= R; ++r) {
    for (unsigned c = 0; c < 4; ++c) {
      if (c + 4*r < 6) continue;
      unsigned u0 = (4*r + 1)*16 + 4*c;
      unsigned u1 = (c == 0) ? (4*(r-1) + 1)*16 + 4*3 : u0 - 4;
      unsigned u2 = (c >= 2) ? (4*(r-1) + 1)*16 + 4*(c-2) : (4*(r-2) + 1)*16 + 4*(c+2);
      if ((c + 4*r)%6 == 0) {
        for (unsigned i = 0; i < 4; ++i) {
          unsigned j = ((i+1)%4);
          vector<pair<GFElement, int>> v;
          v.emplace_back(1, -(u1 + j));
          v.emplace_back(1, u0 + i);
          v.emplace_back(1, u2 + i);
          sys.emplace_back(move(v));
        }
      }
      else {
        for (unsigned i = 0; i < 4; ++i) {
          vector<pair<GFElement, int>> v;
          v.emplace_back(1, u0 + i);
          v.emplace_back(1, u1 + i);
          v.emplace_back(1, u2 + i);
          sys.emplace_back(move(v));
        }
      }
    }
  }*/

// system of equations of AES-256
Matrix AES256eqs(unsigned R) {
  GFElement mc [4] = {2, 3, 1, 1};
  vector<unsigned> shiftRows ({0, 1, 2, 3, 7, 4, 5, 6, 10, 11, 8, 9, 13, 14, 15, 12});
  vector<unsigned> shiftRowsInv (16);
  for (unsigned i = 0; i < 16; ++i) shiftRowsInv[shiftRows[i]] = i;

  vector<vector<pair<GFElement, int>>> sys;
  // P + K0 = X0
  for (unsigned i = 0; i < 16; ++i) {
    vector<pair<GFElement, int>> v;
    v.emplace_back(1, i);
    v.emplace_back(1, i+16);
    v.emplace_back(1, i+32);
    sys.emplace_back(move(v));
  }
  for (unsigned r = 0; r < R-1; ++r) {
    for (unsigned c = 0; c < 4; ++c) {
      for (unsigned i = 0; i < 4; ++i) {
        vector<pair<GFElement, int>> v;
        for (unsigned j = 0; j < 4; ++j) v.emplace_back(mc[(j + 4 - i)%4], (4*r+3)*16 + shiftRowsInv[4*j + c]);
        v.emplace_back(1, (4*r+4)*16 + 4*i + c);
        sys.emplace_back(move(v));
      }
    }
  }
  for (unsigned i = 0; i < 16; ++i) {
    vector<pair<GFElement, int>> v;
    v.emplace_back(1, (4*R-1)*16 + i);
    v.emplace_back(1, (4*R)*16 + shiftRows[i]);
    sys.emplace_back(move(v));
  }
  for (unsigned r = 0; r < R; ++r) {
    for (unsigned i = 0; i < 16; ++i) {
      vector<pair<GFElement, int>> v;
      v.emplace_back(1, (4*r+4)*16 + i);
      v.emplace_back(1, (4*r+5)*16 + i);
      v.emplace_back(1, (4*r+6)*16 + i);
      sys.emplace_back(move(v));
    }
  }
  for (unsigned r = 2; r <= R; ++r) {
    for (unsigned i = 0; i < 4; ++i) {
      for (unsigned j = 1; j < 4; ++j) {
        // K_{r}[i] = K_r[i-4] + K_{r-2}[i]
        vector<pair<GFElement, int>> v;
        v.emplace_back(1, (4*r+1)*16 + 4*i + j);
        v.emplace_back(1, (4*r+1)*16 + 4*i + j - 1);
        v.emplace_back(1, (4*r-7)*16 + 4*i + j);
        sys.emplace_back(move(v));
      }

      unsigned j = 4*((r%2 == 0) ? (i+1)%4 : i) + 3;
      // K_{r}[i] = S(K_{r-1}[j]) + K_{r-2}[i]
      vector<pair<GFElement, int>> v;
      v.emplace_back(1, -((4*r-3)*16 + j));
      v.emplace_back(1, (4*r+1)*16 + 4*i);
      v.emplace_back(1, (4*r-7)*16 + 4*i);
      sys.emplace_back(move(v));

    }
  }
  return Matrix(sys);
}

// system of equations of AES-128
Matrix AES128eqs(unsigned R) {
  GFElement mc [4] = {2, 3, 1, 1};
  vector<unsigned> shiftRows ({0, 1, 2, 3, 7, 4, 5, 6, 10, 11, 8, 9, 13, 14, 15, 12});
  vector<unsigned> shiftRowsInv (16);
  for (unsigned i = 0; i < 16; ++i) shiftRowsInv[shiftRows[i]] = i;

  vector<vector<pair<GFElement, int>>> sys;
  // P + K0 = X0
  for (unsigned i = 0; i < 16; ++i) {
    vector<pair<GFElement, int>> v;
    v.emplace_back(1, i);
    v.emplace_back(1, i+16);
    v.emplace_back(1, i+32);
    sys.emplace_back(move(v));
  }
  for (unsigned r = 0; r < R-1; ++r) {
    for (unsigned c = 0; c < 4; ++c) {
      for (unsigned i = 0; i < 4; ++i) {
        vector<pair<GFElement, int>> v;
        for (unsigned j = 0; j < 4; ++j) v.emplace_back(mc[(j + 4 - i)%4], (4*r+3)*16 + shiftRowsInv[4*c + j]);
        v.emplace_back(1, (4*r+4)*16 + 4*i + c);
        sys.emplace_back(move(v));
      }
    }
  }
  for (unsigned i = 0; i < 16; ++i) {
    vector<pair<GFElement, int>> v;
    v.emplace_back(1, (4*R-1)*16 + i);
    v.emplace_back(1, (4*R)*16 + shiftRows[i]);
    sys.emplace_back(move(v));
  }
  for (unsigned r = 0; r < R; ++r) {
    for (unsigned i = 0; i < 16; ++i) {
      vector<pair<GFElement, int>> v;
      v.emplace_back(1, (4*r+4)*16 + i);
      v.emplace_back(1, (4*r+5)*16 + i);
      v.emplace_back(1, (4*r+6)*16 + i);
      sys.emplace_back(move(v));
    }
  }
    // Key schedule relations
  for (unsigned r = 1; r <= R; ++r) {
      for (unsigned i = 0; i < 16; ++i) {
        if (i%4 == 0) {
            vector<pair<GFElement, int>> v;
            v.emplace_back(1, 16*(4*r +1) + i);
            v.emplace_back(1, 16*(4*(r-1) +1) + i);
            v.emplace_back(1, -(16*(4*(r-1) + 1)  + (i+7)%16));
            sys.emplace_back(move(v));
        }
        else{
            vector<pair<GFElement, int>> v;
            v.emplace_back(1, 16*(4*r +1) + i);
            v.emplace_back(1, 16*(4*(r-1) +1) + i);
            v.emplace_back(1, 16*(4*r +1) + i-1);
            sys.emplace_back(move(v));
        }
      }
                           }
    return Matrix(sys);
}
