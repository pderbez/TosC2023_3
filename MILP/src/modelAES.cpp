#include <iostream>
#include <vector>
#include <fstream>
#include "modelAES.hpp"
#include "CustomCallback.hpp"
#include "SysOfEqs.hpp"
#include "/opt/gurobi/linux64/include/gurobi_c++.h"
//#include "/Library/gurobi952/macos_universal2/include/gurobi_c++.h"

using namespace std;


void addXorConstr(GRBModel & model, const GRBVar & a, const GRBVar & b, const GRBVar & c) {
    model.addConstr(1 - a + b + c >= 1);
    model.addConstr(a + 1 - b + c >= 1);
    model.addConstr(a + b + 1 - c >= 1);
}

void addXorConstr(GRBModel & model, const GRBVar & a, const GRBVar & b, const GRBVar & c, const GRBVar & d) {
    model.addConstr(1 - a + b + c + d >= 1);
    model.addConstr(a + 1 - b + c + d >= 1);
    model.addConstr(a + b + 1 - c + d >= 1);
    model.addConstr(a + b + c + 1 - d >= 1);
}

void addKSConstr(GRBModel & model, vector<GRBVar> & dX, int ck0, int ck1, int ck2, unsigned R) {
  if (ck0 < 0 || ck1 < 0 || ck2 < 0) return;

  int r0 = ck0 / 4;
  int r1 = ck1 / 4;
  int r2 = ck2 / 4;

  int c0 = ck0 % 4;
  int c1 = ck1 % 4;
  int c2 = ck2 % 4;

  if ( r0 <= R && r1 <= R && r2 <= R
      && r0 >= 0 && r1 >= 0 && r2 >= 0 ) {
        for (unsigned i = 0; i < 4; ++i) addXorConstr(model, dX[16*(4*r0+1) + c0 + 4*i], dX[16*(4*r1+1) + c1 + 4*i], dX[16*(4*r2+1) + c2 + 4*i]);
      }

}

void addMCConstr(GRBModel & model, vector<GRBVar> & dX, int ck0, int ck1, int ck2, unsigned R) {

    if (ck0 < 0 || ck1 < 0 || ck2 < 0) return;

    int r0 = ck0 / 4;
    int r1 = ck1 / 4;
    int r2 = ck2 / 4;

    if (abs(r0 - r1) > 1 || abs(r0 - r2) > 1 || abs(r2 - r1) > 1 ) return;

    int c0 = ck0 % 4;
    int c1 = ck1 % 4;
    int c2 = ck2 % 4;

    if ( r0 <= R && r1 <= R && r2 <= R
        && r0 >= 1 && r1 >= 1 && r2 >= 1 ) {
        vector<GRBVar> u(4);
        vector<GRBVar> v(4);

        for (unsigned i = 0; i < 4; ++i) {
            u[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
            v[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
        }
        //u[i] = X_{r0-1}[(c0+i)%4 + 4*i] + X_{r1-1}[(c1+i)%4 + 4*i]
        for (unsigned i = 0; i < 4; ++i) addXorConstr(model, u[i], dX[16*(4*(r0-1)+2) + (c0+i)%4 + 4*i], dX[16*(4*(r1-1)+2) + (c1+i)%4 + 4*i]);

        // v[i] = K_r2[c2 + 4i] + X_r0[c0 + 4i] + X_r1[c1 + 4i]
        for (unsigned i = 0; i < 4; ++i) addXorConstr(model, v[i], dX[16*(4*r2+1) + c2 + 4*i], dX[ 16*(4*r0+2)+ c0 + 4*i ], dX[16*(4*r1+2) + c1 + 4*i ]);

        GRBLinExpr e = 0;
        for (unsigned i = 0; i < 4; ++i) e += u[i] + v[i];
        GRBVar f = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
        model.addConstr(e <= 8*f);
        model.addConstr(e >= 5*f);
    }
}

void modelAES(unsigned R, int version) {

    // P --> + K  --> S --> L --> + K
    // P: 0..15
    // K0: 16..31
    // X0 : 32..47
    // Y0 = S(X0) : 48..63
    // Z0 = L(Y0) : 64..79
    // Kr[i] : (4*r+1)*16 + i
    // Xr[i] : (4*r+2)*16 + i
    // Yr[i] : (4*r+3)*16 + i
    // Zr[i] : (4*r+4)*16 + i
    // C = XR

    // definition of variables

    vector<GRBVar> dX ((4*R + 3)*16); // before ARK

    auto mat = (version == 128 ? AES128eqs(R) : (version == 192 ? AES192eqs(R) : AES256eqs(R)));

    GRBEnv env = GRBEnv(true);
    env.set("LogFile", "mip1.log");
    env.start();


    GRBModel model = GRBModel(env);

    for (unsigned i = 0; i < (4*R+3)*16; ++i) {
        dX[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
    }

    // S
    for (unsigned r = 0; r < R; ++r) {
      for (unsigned i = 0; i < 16; ++i) model.addConstr(dX[16*(4*r +2) + i] == dX[16*(4*r +3) + i]);
    }

    // At least one difference in the state or the key
    GRBLinExpr s = 0;
    for (unsigned i = 0; i < 16; ++i) s += dX[i] + dX[16 + i];
    if (R >= 1 && version != 128) for (unsigned i = 0; i < 16; ++i) s += dX[16*5 + i];
    model.addConstr(s >= 1);

    // Key update
    if (version == 128) {
        for (unsigned r = 1; r <= R; r++) {
          for (unsigned i = 0; i < 16; ++i) { //Key update
            if (i%4 == 0) {
                addXorConstr(model, dX[16*(4*r +1) + i], dX[16*(4*(r-1) +1) + i], dX[16*(4*(r-1) + 1)  + (i+7)%16]);
            }
            else {
                addXorConstr(model, dX[16*(4*r +1) + i], dX[16*(4*(r-1) +1) + i], dX[16*(4*r +1) + i-1] );
                if (i%4 > 1 && r >= 2) {
                    //addXorConstr(model, dX[16*(4*r +1) + i], dX[16*(4*(r-2) +1) + i], dX[16*(4*r +1) + i-2]);
                }
            }
          }
        }
    }
    if (version == 256) {
        for (unsigned r = 2; r <= R; ++r) {
            for (unsigned i = 0; i < 16; ++i) {
                if (i%4 == 0 && r%2 == 0) {
                    addXorConstr(model, dX[16*(4*r +1) + i], dX[16*(4*(r-2) +1) + i] , dX[16*(4*(r-1) + 1)  + ((i+7)%16)]) ;
                }
                else {
                    if (i%4 == 0)
                        addXorConstr(model, dX[16*(4*r +1) + i], dX[16*(4*(r-2) +1) + i], dX[16*(4*(r-1) + 1)  + i+3]);
                    else
                        addXorConstr(model, dX[16*(4*r +1) + i], dX[16*(4*(r-2) +1) + i], dX[16*(4*r +1) + i-1]);
                }
            }
        }
        for (unsigned ck = 0; ck <= 4*R; ++ck) {
          if (ck%4 >= 2) {
              addKSConstr(model, dX, ck, ck-2, ck-16, R);
          }
        }
    }
    if (version == 192) {
        for (unsigned r = 1; r <= R; ++r) {
            for (unsigned i = 0; i < 16; ++i) {
              if (4*r + (i%4) < 6) continue;
              if ((4*r + (i%4))%6 == 0) {
                GRBLinExpr e1 = (i%4 == 0) ? dX[16*(4*(r-1) + 1) + (i+7)%16] : dX[16*(4*r + 1) + (i+3)%16];
                GRBLinExpr e2 = (i%4 == 0) ? dX[16*(4*(r-2) +1) + (i+2)] : dX[16*(4*(r-1) + 1) + (i-2)];
                model.addConstr(1-dX[16*(4*r +1) + i] + e1 + e2 >= 1);
                model.addConstr(dX[16*(4*r +1) + i] + 1-e1 + e2 >= 1);
                model.addConstr(dX[16*(4*r +1) + i] + e1 + 1-e2 >= 1);

              }
              else {
                GRBLinExpr e1 = (i%4 == 0) ? dX[16*(4*(r-1) + 1) + i + 3] : dX[16*(4*r + 1) + i - 1];
                GRBLinExpr e2 = (i%4 < 2) ? dX[16*(4*(r-2) + 1) + (i+2)] : dX[16*(4*(r-1) + 1) + (i-2)];
                model.addConstr(1-dX[16*(4*r +1) + i] + e1 + e2 >= 1);
                model.addConstr(dX[16*(4*r +1) + i] + 1-e1 + e2 >= 1);
                model.addConstr(dX[16*(4*r +1) + i] + e1 + 1-e2 >= 1);
              }
            }
        }
        for (unsigned ck = 0; ck <= 4*R; ++ck) {
          if (ck%6 >= 2) {
              addKSConstr(model, dX, ck, ck-2, ck-12, R);
          }
          if (ck%6 >= 4) {
              addKSConstr(model, dX, ck, ck-4, ck-24, R);
          }
        }

    }


    // MC Constraints
    if (version == 128) {
        for (unsigned ck = 0; ck <= 4*R; ++ck) {
            if (ck%4 != 0) { // C'EST BON ?
                addMCConstr(model, dX, ck, ck-1, ck-4, R);
                addMCConstr(model, dX, ck-1, ck-4, ck, R);
                addMCConstr(model, dX, ck-4, ck, ck-1, R);
                // if (ck%4 >= 2) {
                //     addMCConstr(model, dX, ck, ck-8, ck-2, R);
                //     addMCConstr(model, dX, ck-8, ck-2, ck, R);
                //     addMCConstr(model, dX, ck-2, ck, ck-8, R);
                // }
            }
        }
    }
    if (version == 192) {
        for (unsigned ck = 0; ck <= 4*R; ++ck) {
            if (ck%6 != 0) {
                addMCConstr(model, dX, ck, ck-1, ck-6, R);
                addMCConstr(model, dX, ck-1, ck-6, ck, R);
                addMCConstr(model, dX, ck-6, ck, ck-1, R);
            }
            // if (ck%6 >= 2) {
            //     addMCConstr(model, dX, ck, ck-2, ck-12, R);
            //     addMCConstr(model, dX, ck-2, ck-12, ck, R);
            //     addMCConstr(model, dX, ck-12, ck, ck-2, R);
            // }
            // if (ck%6 >= 4) {
            //     addMCConstr(model, dX, ck, ck-4, ck-24, R);
            //     addMCConstr(model, dX, ck-4, ck-24, ck, R);
            //     addMCConstr(model, dX, ck-24, ck, ck-4, R);
            // }
        }

    }
    if (version == 256) {
        for (unsigned ck = 0; ck <= 4*R; ++ck) {
            if (ck%4 != 0) {
                addMCConstr(model, dX, ck, ck-1, ck-8, R);
                addMCConstr(model, dX, ck-1, ck-8, ck, R);
                addMCConstr(model, dX, ck-8, ck, ck-1, R);
            }
            // if (ck%4 >= 2) {
            //     addMCConstr(model, dX, ck, ck-2, ck-16, R);
            //     addMCConstr(model, dX, ck-2, ck-16, ck, R);
            //     addMCConstr(model, dX, ck-16, ck, ck-2, R);
            // }
        }
    }


    // SR and MC
    for (unsigned r = 0; r < R; ++r) {
      if (r == R-1) {
        for (unsigned c = 0; c < 4; ++c) {
          for (unsigned i = 0; i < 4; ++i) model.addConstr(dX[16*(4*(r+1))  + c  + 4*i] == dX[16*(4*(r) +2) + (c + i)%4 + 4*i]);
        }
      }
      else {
        for (unsigned c = 0; c < 4; ++c) {
            GRBLinExpr e = 0;
            for (unsigned i = 0; i < 4; ++i) e += dX[16*(4*(r+1))  + c  + 4*i];
            for (unsigned i = 0; i < 4; ++i) e += dX[16*(4*(r) +2) + (c + i)%4 + 4*i];
            GRBVar f = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
            model.addConstr(e <= 8*f);
            model.addConstr(e >= 5*f);
            f.set(GRB_IntAttr_BranchPriority, 100);
        }
      }
    }

    // ARK X_r[i] = K_r[i] + Z_{r-1}[i]
    for (unsigned r = 0; r <= R; ++r) {
        for (unsigned i = 0; i < 16; ++i) {
          model.addConstr(1-dX[16*(4*r +2)  + i] + dX[16*(4*r +1)  + i] + dX[16*(4*r)  + i] >= 1);
          model.addConstr(dX[16*(4*r +2)  + i] + 1-dX[16*(4*r +1)  + i] + dX[16*(4*r)  + i] >= 1);
          model.addConstr(dX[16*(4*r +2)  + i] + dX[16*(4*r +1)  + i] + 1-dX[16*(4*r)  + i] >= 1);
        }
    }

    // Count the number of active S-boxes
    GRBLinExpr obj = 0;
    for (unsigned r = 0; r < R; ++r) {
        for (unsigned i = 0; i < 16; ++i) obj += dX[16*(4*r + 2) + i];
    }
    if (version == 128) {
      for (unsigned r = 0; r < R; ++r) {
        for (unsigned i = 0; i < 4; ++i) obj += dX[16*(4*r + 1) + 3 + 4*i];
      }
    }
    if (version == 256) {
      for (unsigned r = 1; r < R; ++r) {
        for (unsigned i = 0; i < 4; ++i) obj += dX[16*(4*r + 1) + 3 + 4*i];
      }
    }

    if (version == 192) {
      for (unsigned r = 1; r < R; ++r) {
        if (r%3 == 1) {
          for (unsigned i = 0; i < 4; ++i) obj += dX[16*(4*r + 1) + 1 + 4*i];
        }
        if (r%3 == 2) {
          for (unsigned i = 0; i < 4; ++i) obj += dX[16*(4*r + 1) + 3 + 4*i];
        }
      }
      if (R%3 == 1) {
        for (unsigned i = 0; i < 4; ++i) obj += dX[16*(4*R + 1) + 1 + 4*i];
      }
    }



    model.setObjective(obj, GRB_MINIMIZE);

    mycallback cb ((4*R+3)*16, dX.data(), mat, version);
    model.setCallback(&cb);
    model.set(GRB_IntParam_LazyConstraints , 1);
    model.set(GRB_IntParam_Cuts , 2);
    model.set(GRB_IntParam_MIPFocus , 3);
    // Limit how many solutions to collect
    model.set(GRB_IntParam_PoolSolutions, 2000000);

    // Limit the search space by setting a gap for the worst possible solution that will be accepted
    model.set(GRB_DoubleParam_PoolGap, 0.001);

    // do a systematic search for the k-best solutions
    model.set(GRB_IntParam_PoolSearchMode, 2);
    //model.read("tune.prm");

    model.optimize();

    if (model.get(GRB_IntAttr_Status) != GRB_OPTIMAL) return;

    for (unsigned r = 0; r <= R; ++r) {
      for (unsigned i = 0; i < 16; ++i) {
        cout << (dX[16*4*r + i].get(GRB_DoubleAttr_X) < 0.5 ? 0 : 1) << " ";
        if (i%4 == 3) cout << endl;
      }
      cout << endl;
      for (unsigned i = 0; i < 16; ++i) {
        cout << (dX[16*(4*r +1) + i].get(GRB_DoubleAttr_X) < 0.5 ? 0 : 1) << " ";
        if (i%4 == 3) cout << endl;
      }
      cout << endl;
      for (unsigned i = 0; i < 16; ++i) {
        cout << (dX[16*(4*r +2) + i].get(GRB_DoubleAttr_X) < 0.5 ? 0 : 1) << " ";
        if (i%4 == 3) cout << endl;
      }
      cout << endl << " --------------- " << endl;
    }
    //getchar();
}
