#include <iostream>
#include <vector>
#include <fstream>
#include "modelAES.hpp"
#include "CustomCallback.hpp"
#include "SysOfEqs.hpp"
#include "/opt/gurobi/linux64/include/gurobi_c++.h"
//#include "/Library/gurobi952/macos_universal2/include/gurobi_c++.h"

using namespace std;

void propagateForward(vector<GRBVar> & dX, unsigned r, unsigned R, GRBModel & model) {
  // force propagation of proba 1 from Round r to Round r+1
  if (r < R-1) {
    for (unsigned c = 0; c < 4; ++c) {
      for (unsigned i = 1; i < 4; ++i) model.addConstr(dX[16*(4*r + 4) + c + 4*i] == dX[16*(4*r + 4) + c]);
      GRBLinExpr e = 0;
      for (unsigned i = 0; i < 4; ++i) e += dX[16*(4*r +2) + ((c + 5*i)%4) + 4*i];
      model.addConstr(dX[16*(4*r + 4) + c] <= e);
      model.addConstr(4*dX[16*(4*r + 4) + c] >= e);
      dX[16*(4*r + 4) + c].set(GRB_IntAttr_BranchPriority, 100);
      for (unsigned i = 0; i < 4; ++i) model.addConstr(dX[16*(4*r + 4) + c + 4*i] <= dX[16*(4*r + 6) + c + 4*i]);
    }
  }
  else {
    for (unsigned c = 0; c < 4; ++c) {
      for (unsigned i = 0; i < 4; ++i) {
        model.addConstr(dX[16*(4*r + 4) + c + 4*i] == dX[16*(4*r +2) + ((c + 5*i)%4) + 4*i]);
        model.addConstr(dX[16*(4*r + 4) + c + 4*i] <= dX[16*(4*r + 6) + c + 4*i]);
      }
    }
  }
}

void relationRoundDiff(vector<GRBVar> & dX, unsigned r, unsigned R, GRBModel & model) {
  // relation diff between Round r and Round r+1
  if (r < R - 1) {
    for(unsigned c = 0; c < 4; c++){
        GRBLinExpr e = 0;
        for (unsigned i = 0; i < 4; ++i) e += dX[16*(4*(r+1))  + c  + 4*i];
        for (unsigned i = 0; i < 4; ++i) e += dX[16*(4*(r) +2) + ((c + 5*i)%4) + 4*i];
        GRBVar f = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
        model.addConstr(e <= 8*f);
        model.addConstr(e >= 5*f);
        f.set(GRB_IntAttr_BranchPriority, 100);
    }
  }
 }

void propagateBackward(vector<GRBVar> & dX, unsigned r, unsigned R, GRBModel & model) {
  // force propagation of proba 1 from Round r to Round r-1
  if (r > 0 && r < R) {
    for(unsigned c = 0; c < 4; c++){
        GRBLinExpr e = 0;
        for (unsigned i = 0; i < 4; ++i) e += dX[16*(4*(r-1) + 4)  + c  + 4*i];
        for (unsigned i = 0; i < 4; ++i) e += dX[16*(4*(r-1) + 2) + ((c + 5*i)%4) + 4*i];
        GRBVar f = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
        model.addConstr(e <= 8*f);
        model.addConstr(e >= 5*f);
        f.set(GRB_IntAttr_BranchPriority, 100);
    }
    for(unsigned c = 0; c < 4; c++){
      GRBLinExpr e = 0;
      for (unsigned i = 0; i < 4; ++i) e += dX[16*(4*r + 2)  + c  + 4*i];
      for (unsigned i = 0; i < 4; ++i) {
        model.addConstr(4*dX[16*(4*(r-1) +2) + ((c + 5*i)%4) + 4*i] >= e);
      }
    }
  }
}






void modelAES(unsigned rin, unsigned rmiddle, unsigned rout, int version) {

  unsigned R = rin + rmiddle + rout;

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

    for (unsigned r = 0; r < R; ++r) {
      for (unsigned i = 0; i < 16; ++i) model.addConstr(dX[16*(4*r +2) + i] == dX[16*(4*r +3) + i]);
    }

    // At least one difference in the state or the key
    GRBLinExpr s = 0;
    for (unsigned i = 0; i < 16; ++i) s += dX[i];
    for (unsigned r = 0; r <= R; ++r) {
      for (unsigned i = 0; i < 16; ++i) s += dX[16*(4*r +2) + i] + dX[16*(4*r +1) + i];
    }
    model.addConstr(s >= 1);

    if(version == 256){
        for(unsigned r = 2; r <= R; r++){
            for (unsigned i = 0; i < 16; ++i) { //Key update
                if (i%4 == 0 && r%2 == 0) {
                    model.addConstr(1-dX[16*(4*r +1) + i] + dX[16*(4*(r-2) +1) + i] + dX[16*(4*(r-1) + 1)  + ((i+7)%16)] >= 1);
                    model.addConstr(dX[16*(4*r +1) + i] + 1-dX[16*(4*(r-2) +1) + i] + dX[16*(4*(r-1) + 1)  + ((i+7)%16)] >= 1);
                    model.addConstr(dX[16*(4*r +1) + i] + dX[16*(4*(r-2) +1) + i] + 1-dX[16*(4*(r-1) + 1)  + ((i+7)%16)] >= 1);
                }
                else {
                    if (i%4 == 0) {
                        model.addConstr(1-dX[16*(4*r +1) + i] + dX[16*(4*(r-2) +1) + i] + dX[16*(4*(r-1) + 1)  + i+3] >= 1);
                        model.addConstr(dX[16*(4*r +1) + i] + 1-dX[16*(4*(r-2) +1) + i] + dX[16*(4*(r-1) + 1)  + i+3] >= 1);
                        model.addConstr(dX[16*(4*r +1) + i] + dX[16*(4*(r-2) +1) + i] + 1-dX[16*(4*(r-1) + 1)  + i+3] >= 1);
                    }
                    else {
                        model.addConstr(1-dX[16*(4*r +1) + i] + dX[16*(4*(r-2) +1) + i] + dX[16*(4*r +1) + i-1] >= 1);
                        model.addConstr(dX[16*(4*r +1) + i] + 1-dX[16*(4*(r-2) +1) + i] + dX[16*(4*r +1) + i-1] >= 1);
                        model.addConstr(dX[16*(4*r +1) + i] + dX[16*(4*(r-2) +1) + i] + 1-dX[16*(4*r +1) + i-1] >= 1);
                    }
                }
            }
        }
    }
    // Version à modifier
    if(version == 192){
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
    }
    else if(version == 128){
        for(unsigned r = 1; r <= R; r++){
          for (unsigned i = 0; i < 16; ++i) { //Key update
            if (i%4 == 0) {
                model.addConstr(1-dX[16*(4*r +1) + i] + dX[16*(4*(r-1) +1) + i] + dX[16*(4*(r-1) + 1)  + (i+7)%16] >= 1);
                model.addConstr(dX[16*(4*r +1) + i] + 1-dX[16*(4*(r-1) +1) + i] + dX[16*(4*(r-1) + 1)  + (i+7)%16] >= 1);
                model.addConstr(dX[16*(4*r +1) + i] + dX[16*(4*(r-1) +1) + i] + 1-dX[16*(4*(r-1) + 1)  + (i+7)%16] >= 1);
            }
            else {
                model.addConstr(1-dX[16*(4*r +1) + i] + dX[16*(4*(r-1) +1) + i] + dX[16*(4*r +1) + i-1] >= 1);
                model.addConstr(dX[16*(4*r +1) + i] + 1-dX[16*(4*(r-1) +1) + i] + dX[16*(4*r +1) + i-1] >= 1);
                model.addConstr(dX[16*(4*r +1) + i] + dX[16*(4*(r-1) +1) + i] + 1-dX[16*(4*r +1) + i-1] >= 1);
                if (i%4 > 1 && r >= 2) {
                  model.addConstr(1-dX[16*(4*r +1) + i] + dX[16*(4*(r-2) +1) + i] + dX[16*(4*r +1) + i-2] >= 1);
                  model.addConstr(dX[16*(4*r +1) + i] + 1-dX[16*(4*(r-2) +1) + i] + dX[16*(4*r +1) + i-2] >= 1);
                  model.addConstr(dX[16*(4*r +1) + i] + dX[16*(4*(r-2) +1) + i] + 1-dX[16*(4*r +1) + i-2] >= 1);
                }
              }

          }
        }
    }

    for (unsigned r = 0; r <= R ; ++r) {
      //ARK (sX[r+1], kX[r+1] and zX[r+1])
      for (unsigned i = 0; i < 16; ++i) {
          model.addConstr(1-dX[16*(4*r +2)  + i] + dX[16*(4*r +1)  + i] + dX[16*(4*r)  + i] >= 1);
          model.addConstr(dX[16*(4*r +2)  + i] + 1-dX[16*(4*r +1)  + i] + dX[16*(4*r)  + i] >= 1);
          model.addConstr(dX[16*(4*r +2)  + i] + dX[16*(4*r +1)  + i] + 1-dX[16*(4*r)  + i] >= 1);
      }
    }

    // Force diff on P and C
    for (unsigned i = 0; i < 16; ++i) {
      model.addConstr(2*dX[i] >= dX[16 + i] + dX[32 + i]);
      model.addConstr(2*dX[16*(4*R + 2) + i] >= dX[16*4*R + i] + dX[16*(4*R + 1) + i]);
    }

    //Distinguisher part
    for(unsigned r = (rin == 0 ? rin : rin-1); r < rin + rmiddle; r++){
      relationRoundDiff(dX, r, R, model);
    }

    //rin
    for (unsigned r = 0; r < rin; ++r) propagateBackward(dX, r, R, model);

    //rout
    for (unsigned r = rin + rmiddle; r < R; ++r) propagateForward(dX, r, R, model);

    GRBLinExpr prob = 0;
    GRBLinExpr kin = 0;
    GRBLinExpr kout = 0;

    for (unsigned r = 0; r < rin; ++r) {
      for (unsigned i = 0; i < 16; ++i) kin += dX[16*(4*r + 2) + i];
    }

    for (unsigned r = rin; r < rin+rmiddle; ++r) {
      for (unsigned i = 0; i < 16; ++i) prob += dX[16*(4*r + 2) + i];
    }

    for (unsigned r = rin+rmiddle; r < R; ++r) {
      for (unsigned i = 0; i < 16; ++i) kout += dX[16*(4*r + 2) + i];
    }

    // force the last column of the round keys to 0 in the middle rounds
 /*   for (unsigned r = rin+3; r <= rin+rmiddle - 3; ++r) {
      for (unsigned i = 0; i < 4; ++i) model.addConstr(dX[16*(4*r + 1) + 4*i + 3] == 0);
    }*/

    if (version == 192) {
      for (unsigned r = rin+1; r <= rin+rmiddle - 1; ++r) {
           for (unsigned i = 0; i < 16; ++i)
           {
             if ((r+(i%4))%6 == 5) model.addConstr(dX[16*(4*r + 1) + i] == 0);
           }
         }
    }
    else {
      for (unsigned r = rin+1; r <= rin+rmiddle - 1; ++r) {
           for (unsigned i = 0; i < 16; ++i)
           {
             if (i%4 == 3) model.addConstr(dX[16*(4*r + 1) + i] == 0);
           }
         }
    }




    // for (unsigned i = 0; i < 16; ++i) {
    //   model.addConstr(dX[16*(4*3 + 2) + i] == 0);
    //   model.addConstr(dX[16*(4*5 + 2) + i] == 0);
    //   model.addConstr(dX[16*(4*7 + 2) + i] == 0);
    //   model.addConstr(dX[16*(4*9 + 2) + i] == 0);
    // }
    // //
    // {
    //   int mytab[] = {4,6};
    //   for (int r = 0; r < 2; ++r)
    //   {
    //     GRBLinExpr e = 0;
    //     for (unsigned i = 0; i < 16; ++i) e += dX[16*(4*mytab[r] + 2) + i];
    //     model.addConstr(e == 2);
    //   }
    // }
    //
    // {
    //   int mytab[] = {8};
    //   for (int r = 0; r < 1; ++r)
    //   {
    //     GRBLinExpr e = 0;
    //     for (unsigned i = 0; i < 16; ++i) e += dX[16*(4*mytab[r] + 2) + i];
    //     model.addConstr(e == 1);
    //   }
    // }
    //
    // {
    //   int mytab[] = {1};
    //   for (int r = 0; r < 1; ++r)
    //   {
    //     GRBLinExpr e = 0;
    //     for (unsigned i = 0; i < 16; ++i) e += dX[16*(4*mytab[r] + 2) + i];
    //     model.addConstr(e == 3);
    //   }
    // }
    //
    // {
    //   int mytab[] = {0};
    //   for (int r = 0; r < 1; ++r)
    //   {
    //     GRBLinExpr e = 0;
    //     for (unsigned i = 0; i < 16; ++i) e += dX[16*(4*mytab[r] + 2) + i];
    //     model.addConstr(e == 9);
    //   }
    // }
    //
    // {
    //   int mytab[] = {2, 10, 11};
    //   for (int r = 0; r < 3; ++r)
    //   {
    //     GRBLinExpr e = 0;
    //     for (unsigned i = 0; i < 16; ++i) e += dX[16*(4*mytab[r] + 2) + i];
    //     model.addConstr(e == 4);
    //   }
    // }



    // Constraint on the number of active Sboxes in the distinguisher
    model.addConstr(prob <= 21);
    GRBVar obj = model.addVar(0.0, 255.0, 0.0, GRB_INTEGER);

    model.addConstr(obj <= version - 1);
    model.addConstr(obj >= 8*kin + prob*6);
    model.addConstr(obj >= 8*kout + prob*6);

    GRBVar obj_e = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
    model.addConstr(obj <= 8*kin + prob*6 + 255*obj_e);
    model.addConstr(obj <= 8*kout + prob*6 + 255*(1 - obj_e));




    GRBLinExpr obj_l = obj;

    model.setObjective(obj_l, GRB_MINIMIZE);

    mycallback cb ((4*R+3)*16, dX.data(), mat, version);
    model.setCallback(&cb);
    model.set(GRB_IntParam_LazyConstraints , 1);
    model.read("tune.prm");

    model.optimize();

    //getchar();

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

    ofstream myfile;
    myfile.open ("attack" + to_string(obj.get(GRB_DoubleAttr_X)) + "_" + to_string(rin) + "_" + to_string(rmiddle) + "_" + to_string(rout)  + ".tex");
    myfile << "\\documentclass[10pt,a4paper]{article}" << endl;
    myfile << "\\usepackage[utf8]{inputenc}" << endl;
    myfile << "\\usepackage{amsmath}" << endl;
    myfile << "\\usepackage{amsfonts}" << endl;
    myfile << "\\usepackage{amssymb}" << endl;

    myfile << "\\usepackage{tikz}" << endl;
    myfile << "\\usepackage{subfigure}" << endl;
    myfile << "\\usetikzlibrary{fit}" << endl;
    myfile << "\\usetikzlibrary{shapes.geometric}" << endl;
    myfile << "\\usetikzlibrary{patterns}" << endl;

    myfile << "\\begin{document}" << endl;



    myfile << "\\begin{figure}" << endl;
    myfile << "\\centering" << endl;
    myfile << "\\begin{tikzpicture}[scale=0.26,>=stealth]" << endl;

    myfile << "\\begin{scope} " << endl;
    myfile << "\\begin{scope}[xshift=6cm]" << endl;
    myfile << "\\path (3,2) node[anchor=east] {$P$};" << endl;
    myfile << "\\draw[->] (3,2) --  +(3,0);" << endl;
    myfile << "\\end{scope}" << endl;

    myfile << "\\begin{scope}[xshift=8.5cm,yshift=-4.5cm]" << endl;
    myfile << "\\fill[white] (0,0) rectangle +(4,4);" << endl;
    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < 4; ++j) {
        if (dX[16 + 4*i + j].get(GRB_DoubleAttr_X) > 0.5)   myfile << "\\fill[gray!20] (" << j << "," << 3-i << ") rectangle +(1,1);" << endl;
      }
    }

    myfile << "\\draw (0,0) rectangle (4,4);" << endl;
    myfile << "\\draw (0,1) -- +(4,0);" << endl;
    myfile << "\\draw (0,2) -- +(4,0);" << endl;
    myfile << "\\draw (0,3) -- +(4,0);" << endl;
    myfile << "\\draw (1,0) -- +(0,4);" << endl;
    myfile << "\\draw (2,0) -- +(0,4);" << endl;
    myfile << "\\draw (3,0) -- +(0,4);" << endl;
    myfile << "\\node at (5,2) {$k_0$};" << endl;
    myfile << "\\draw (2,4) -- +(0,2.5);" << endl;
    myfile << "\\end{scope}" << endl;

    for (unsigned r = 0; r < R; ++r) {
      unsigned x = 6 + 17*(r%3);
      unsigned y = 12*(r/3);

      myfile << "\\begin{scope}[xshift=" << x << "cm, yshift=-" << y << "cm] % round " << r << endl;
      myfile << "\\begin{scope}[xshift=6cm]" << endl;
      myfile << "\\fill[white] (0,0) rectangle +(4,4);" << endl;

      for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
          if (dX[16*(4*r + 2) + 4*i + j].get(GRB_DoubleAttr_X) > 0.5)   myfile << "\\fill[gray!20] (" << j << "," << 3-i << ") rectangle +(1,1);" << endl;
        }
      }


      myfile << "\\draw (0,0) rectangle (4,4);" << endl;
      myfile << "\\draw (0,1) -- +(4,0);" << endl;
      myfile << "\\draw (0,2) -- +(4,0);" << endl;
      myfile << "\\draw (0,3) -- +(4,0);" << endl;
      myfile << "\\draw (1,0) -- +(0,4);" << endl;
      myfile << "\\draw (2,0) -- +(0,4);" << endl;
      myfile << "\\draw (3,0) -- +(0,4);" << endl;
      myfile << "\\node at (2,5) {$x_{" << r << "}$};" << endl;
      myfile << "\\draw[->] (4,2) --  +(1,0);" << endl;
      myfile << "\\end{scope}" << endl;

      myfile << "\\begin{scope}[xshift=11cm]" << endl;
      myfile << "\\fill[white] (0,0) rectangle +(4,4);" << endl;

      for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
          if (dX[16*(4*r + 2) + 4*i + j].get(GRB_DoubleAttr_X) > 0.5)   myfile << "\\fill[gray!20] (" << (j+i + 4)%4 << "," << 3-i << ") rectangle +(1,1);" << endl;
        }
      }


      myfile << "\\draw (0,0) rectangle (4,4);" << endl;
      myfile << "\\draw (0,1) -- +(4,0);" << endl;
      myfile << "\\draw (0,2) -- +(4,0);" << endl;
      myfile << "\\draw (0,3) -- +(4,0);" << endl;
      myfile << "\\draw (1,0) -- +(0,4);" << endl;
      myfile << "\\draw (2,0) -- +(0,4);" << endl;
      myfile << "\\draw (3,0) -- +(0,4);" << endl;
      myfile << "\\node at (2,5) {$z_{" << r << "}$};" << endl;
      myfile << "\\draw[->] (4,2) --  +(1,0);" << endl;
      myfile << "\\end{scope}" << endl;

      myfile << "\\begin{scope}[xshift=16cm]" << endl;
      myfile << "\\fill[white] (0,0) rectangle +(4,4);" << endl;

      for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
          if (dX[16*(4*r + 4) + 4*i + j].get(GRB_DoubleAttr_X) > 0.5)   myfile << "\\fill[gray!20] (" << j << "," << 3-i << ") rectangle +(1,1);" << endl;
        }
      }


      myfile << "\\draw (0,0) rectangle (4,4);" << endl;
      myfile << "\\draw (0,1) -- +(4,0);" << endl;
      myfile << "\\draw (0,2) -- +(4,0);" << endl;
      myfile << "\\draw (0,3) -- +(4,0);" << endl;
      myfile << "\\draw (1,0) -- +(0,4);" << endl;
      myfile << "\\draw (2,0) -- +(0,4);" << endl;
      myfile << "\\draw (3,0) -- +(0,4);" << endl;
      myfile << "\\node at (2,5) {$w_{" << r << "}$};" << endl;
      myfile << "\\draw[->] (4,2) --  +(3,0);" << endl;
      myfile << "\\end{scope}" << endl;

      myfile << "\\begin{scope}[xshift=19.5cm,yshift=-4.5cm]" << endl;
      myfile << "\\fill[white] (0,0) rectangle +(4,4);" << endl;

      for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
          if (dX[16*(4*(r+1) + 1) + 4*i + j].get(GRB_DoubleAttr_X) > 0.5)   myfile << "\\fill[gray!20] (" << j << "," << 3-i << ") rectangle +(1,1);" << endl;
        }
      }


      myfile << "\\draw (0,0) rectangle (4,4);" << endl;
      myfile << "\\draw (0,1) -- +(4,0);" << endl;
      myfile << "\\draw (0,2) -- +(4,0);" << endl;
      myfile << "\\draw (0,3) -- +(4,0);" << endl;
      myfile << "\\draw (1,0) -- +(0,4);" << endl;
      myfile << "\\draw (2,0) -- +(0,4);" << endl;
      myfile << "\\draw (3,0) -- +(0,4);" << endl;
      myfile << "\\node at (5,2) {$k_{" << r+1 << "}$};" << endl;
      myfile << "\\draw (2,4) --  +(0,2.5);" << endl;
      myfile << "\\end{scope}" << endl;


      myfile << "\\end{scope} % round " << r << endl;
    }

    myfile << "\\end{scope}" << endl;

    myfile << "\\end{tikzpicture}" << endl;
    myfile << "\\end{figure}" << endl;


    myfile << "\\end{document}" << endl;


    myfile.close();


}
