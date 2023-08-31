#include "CustomCallback.hpp"

using namespace std;

void mycallback::callback()
{
	static unsigned mylazies = 0;
	try {
		// if (where == GRB_CB_MIPNODE) {
		//   double * kXup = getNodeRel(varsKup, numvars);
		//   double * kXlo = getNodeRel(varsKlo, numvars);
		//   double * kKup = getNodeRel(varsSKup, mapSKup.size());
		//   double * kKlo = getNodeRel(varsSKlo, mapSKup.size());
		//
		//   double * zXup = getNodeRel(varsZup, numvars);
		//   double * dXup = getNodeRel(varsDup, numvars);
		//   double * sXup = getNodeRel(varsSup, numvars);
		//   double * zXlo = getNodeRel(varsZlo, numvars);
		//   double * dXlo = getNodeRel(varsDlo, numvars);
		//   double * sXlo = getNodeRel(varsSlo, numvars);
		//
		//   // for (unsigned r = 0; r < numvars/16; ++r) {
		//   //   cout << ((r%4 == 0) ? 4 : r%4)  << ": ";
		//   //   for (unsigned i = 0; i < 16; ++i) {
		//   //     cout << "(" << round(dXup[16*r + i]*10)/10.0 << "," << round(zXup[16*r + i]*10)/10.0 << "," << round(kXup[16*r + i]*10)/10.0 << "," << round(sXup[16*r + i]*10)/10.0 << ") ";
		//   //   }
		//   //   cout << endl;
		//   // }
		//   // cout << " ---- " << endl;
		//   // for (unsigned r = 0; r < numvars/16; ++r) {
		//   //   cout << ((r%4 == 0) ? 4 : r%4)  << ": ";
		//   //   for (unsigned i = 0; i < 16; ++i) {
		//   //     cout << "(" << round(dXlo[16*r + i]*10)/10.0 << "," << round(zXlo[16*r + i]*10)/10.0 << "," << round(kXlo[16*r + i]*10)/10.0 << "," << round(sXlo[16*r + i]*10)/10.0 << ") ";
		//   //   }
		//   //   cout << endl;
		//   // }
		//   for (unsigned r = 0; r < numvars/16; ++r) {
		//     cout << (r%4 == 0 ? 4 : r%4) << ": ";
		//     double cptd = 0, cpts = 0, cptk = 0, cptz = 0;
		//     for (unsigned i = 0; i < 16; ++i) {
		//       cptd += dXup[16*r + i];
		//       cptz += zXup[16*r + i];
		//       cptk += kXup[16*r + i];
		//       cpts += sXup[16*r + i];
		//     }
		//     cout << "(" << cptd << ", " << cptz << ", " << cptk << ", " << cpts << ")";
		//     cptd = 0; cpts = 0; cptk = 0; cptz = 0;
		//     for (unsigned i = 0; i < 16; ++i) {
		//       cptd += dXlo[16*r + i];
		//       cptz += zXlo[16*r + i];
		//       cptk += kXlo[16*r + i];
		//       cpts += sXlo[16*r + i];
		//     }
		//     cout << " | (" << cptd << ", " << cptz << ", " << cptk << ", " << cpts << ")";
		//     cout << endl;
		//   }
		//   getchar();
		//   delete[] zXup;
		//   delete[] dXup;
		//   delete[] sXup;
		//   delete[] zXlo;
		//   delete[] dXlo;
		//   delete[] sXlo;
		//
		//   // vector<double> kX (numvars), kK (mapSKup.size());
		//   //
		//   // for (unsigned i = 0; i < numvars; ++i) {
		//   //   if (kXup[i] > 0.99 || kXlo[i] > 0.99) kX[i] = 1.0;
		//   //   else kX[i] = 0;
		//   // }
		//   // for (unsigned i = 0; i < mapSKup.size(); ++i) {
		//   //   if (kKup[i] > 0.99 || kKlo[i] > 0.99) kK[i] = 1.0;
		//   //   else kX[i] = 0;
		//   // }
		//   // auto vup = mat.checkK(kX.data(), kK.data(), mapSKup);
		//   // if (!vup.empty()) {
		//   //   GRBLinExpr e = 0;
		//   //   for (unsigned i = 0; i < numvars; ++i) {
		//   //     if (kXup[i] > 0.99) e += 1 - varsKup[i];
		//   //     else if (kXlo[i] > 0.99) e += 1 - varsKlo[i];
		//   //   }
		//   //   for (unsigned i = 0; i < mapSKup.size(); ++i) {
		//   //     if (kKup[i] > 0.99) e += 1 - varsSKup[i];
		//   //     else if (kKlo[i] > 0.99) e += 1 - varsSKlo[i];
		//   //   }
		//   //   for (auto x : vup) {
		//   //     if (getNodeRel(varsZup[x]) + getNodeRel(varsZlo[x]) < 0.9) {
		//   //       addLazy(varsZup[x] + varsZlo[x] + e >= 1);
		//   //       cout << "\r" << ++mylazies;
		//   //     }
		//   //
		//   //   }
		//   // }
		//   //cout << "vup: " << vup.size() << endl;
		//   delete[] kXlo;
		//   delete[] kXup;
		//   delete[] kKup;
		//   delete[] kKlo;
		//
		// }
		if (where == GRB_CB_MIPSOL) {
			// {
			//   double * zXup = getSolution(varsZup, numvars);
			//   double * zXlo = getSolution(varsZlo, numvars);
			//   double * kXup = getSolution(varsKup, numvars);
			//   double * kXlo = getSolution(varsKlo, numvars);
			//   cout << "--- z ---" << endl;
			//   for (unsigned r = 0; r < numvars/16; ++r) {
			//     if (r%4 != 1) continue;
			//     cout << ((r%4 == 0) ? 4 : r%4)  << ": ";
			//     for (unsigned i = 0; i < 16; ++i) cout << ((zXup[16*r + i] < 0.5) ? 0 : 1) << " ";
			//     cout << "    |    ";
			//     for (unsigned i = 0; i < 16; ++i) cout << ((zXlo[16*r + i] < 0.5) ? 0 : 1) << " ";
			//     cout << endl;
			//   }
			//
			//   cout << "--- s ---" << endl;
			//   for (unsigned r = 0; r < numvars/16; ++r) {
			//     if (r%4 != 1) continue;
			//     cout << ((r%4 == 0) ? 4 : r%4)  << ": ";
			//     for (unsigned i = 0; i < 16; ++i) cout << ((kXup[16*r + i] < 0.5) ? 0 : 1) << " ";
			//     cout << "    |    ";
			//     for (unsigned i = 0; i < 16; ++i) cout << ((kXlo[16*r + i] < 0.5) ? 0 : 1) << " ";
			//     cout << endl;
			//   }
			//   delete[] zXup;
			//   delete[] kXup;
			//   delete[] zXlo;
			//   delete[] kXlo;
			// }
			bool flagAddLazy = false;
			double * X = getSolution(varsD, numvars);
			unsigned l = mat.checkZ(X);
			if (l < mat.nblines) {
				GRBLinExpr e = 1-varsD[abs(mat.getFront(l))];
				for (unsigned c = 0; c < mat.nbcols; ++c) {
					if (mat(l, c) != 0) e += varsD[abs(mat.getColumns(c))];
				}
				addLazy(e >= 1);
				//mat.printLine(l);
				//cout << "issue 1: " << (abs(mat.getFront(l))/16)%4 << endl;
				flagAddLazy = true;
			}
			else {
				// auto v = (version == 192 ? mat.checkK192(X) : mat.checkK128_256(X));
				// unsigned weak = 0;
				// if (v.size() > weak) {
				// 	GRBLinExpr e = 0;
				// 	for (auto x : v) e += 1-varsD[x];
				// 	for (unsigned i = 0; i < numvars; ++i) {
				// 		if (X[i] < 0.5) e += v.size() * varsD[i];
				// 	}
				// 	addLazy(e >= v.size() - weak);
				// 	//cout << "issue 3" << endl;
				// 	flagAddLazy = true;
				// }
			}
			delete[] X;
			if (flagAddLazy) cout << "\r" << ++mylazies;
			//getchar();
		}

	} catch (GRBException e) {
		cout << "Error number: " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	} catch (...) {
		cout << "Error during callback" << endl;
	}

}
