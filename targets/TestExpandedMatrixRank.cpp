#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include "TPZBndCondT.h"
#include "TPZCompElDisc.h"
#include "TPZCompMeshTools.h"
#include "TPZExtendGridDimension.h"
#include "TPZGenGrid2D.h"
#include "TPZGeoCube.h"
#include "TPZGeoLinear.h"
#include "TPZGeoMeshTools.h"
#include "TPZLinearAnalysis.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "TPZRefPattern.h"
#include "TPZRefPatternTools.h"
#include "TPZShapeH1.h"
#include "TPZShapeHDiv.h"
#include "TPZShapeHDivConstant.h"
#include "TPZShapeHDivOptimized.h"
#include "TPZVTKGeoMesh.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzfstrmatrix.h"
#include "pzgeoel.h"
#include "pzgeoelrefless.h"
#include "pzgeoelside.h"
#include "pzgeopoint.h"
#include "pzgeoquad.h"
#include "pzgeotetrahedra.h"
#include "pzintel.h"
#include "pzmanvector.h"
#include "pzmultiphysicscompel.h"
#include "pzmultiphysicselement.h"
#include "pzshapecube.h"
#include "pzshapelinear.h"
#include "pzshapepiram.h"
#include "pzshapepiramHdiv.h"
#include "pzshapeprism.h"
#include "pzshapequad.h"
#include "pzshapetetra.h"
#include "pzshapetriang.h"
#include "pzshtmat.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include "pztrnsform.h"
#include "pzvec_extras.h"
#include "tpzarc3d.h"
#include "tpzautopointer.h"
#include "tpzgeoblend.h"
#include "tpzgeoelrefpattern.h"
#include "tpzintpoints.h"
#include "tpzpermutation.h"
#include "tpztriangle.h"

template <class TSHAPE>
void TestHDivOptimized(int kFacet);

template <class TSHAPE>
void filterEquations(TPZVec<int> &filteredIndices, int kFacet);

int main(int argc, char **argv) {
  int kFacet = 3;
  TestHDivOptimized<pzshape::TPZShapeCube>(kFacet);
  // TestHDivOptimized<pzshape::TPZShapeQuad>(kFacet);
}

template <class TSHAPE>
void TestHDivOptimized(int kFacet) {
  TPZMaterialDataT<REAL> dataHDiv, dataHDivOptimized, dataHDivConst, dataH1;
  TPZManVector<int64_t, 27> ids(TSHAPE::NCornerNodes, 0);
  for (int i = 0; i < TSHAPE::NCornerNodes; i++) {
    ids[i] = i;
  }
  TPZManVector<int, 27> sideorient(TSHAPE::NFacets, 1);
  TPZManVector<int, 27> orders(TSHAPE::NFacets + 1, kFacet);
  TPZManVector<int, 27> h1orders(TSHAPE::NSides - TSHAPE::NCornerNodes, kFacet);
  TPZShapeHDiv<TSHAPE> hdiv_std;
  TPZShapeHDivOptimized<TSHAPE> hdiv_opt;
  TPZShapeHDivConstant<TSHAPE> hdiv_const;
  TPZShapeH1<TSHAPE> h1;

  hdiv_std.Initialize(ids, orders, sideorient, dataHDiv);
  hdiv_opt.Initialize(ids, orders, sideorient, dataHDivOptimized);
  hdiv_const.Initialize(ids, orders, sideorient, dataHDivConst);
  h1.Initialize(ids, h1orders, dataH1);
  int nshape_std = hdiv_std.NShapeF(dataHDiv);
  int nshapeInt_std = hdiv_std.NConnectShapeF(TSHAPE::NFacets, dataHDiv);
  int nshape_opt = hdiv_opt.NShapeF(dataHDivOptimized);
  int nshape_cust = nshape_std;
  int nshape_h1 = dataH1.fH1.fPhi.Rows();

  TPZVec<int> filteredIndices;
  filterEquations<TSHAPE>(filteredIndices, kFacet);
  int nFiltered = filteredIndices.size();
  if (nFiltered != nshape_h1-1) {
    std::cout << "Warning: Number of filtered indices (" << nFiltered << ") does not match expected (" << nshape_h1-1 << ")" << std::endl;
    DebugStop();
  }

  constexpr int dim = TSHAPE::Dimension;
  TPZFMatrix<REAL> phi_std(dim, nshape_std, 0.);
  TPZFMatrix<REAL> phi_cust(dim, nshape_cust, 0.);
  TPZFNMatrix<60, REAL> div_std(nshape_std, 1);
  TPZFNMatrix<60, REAL> div_filtered(nFiltered, 1);
  TPZFMatrix<REAL> phi_h1(nshape_h1, 1, 0.);
  TPZFMatrix<REAL> phi_h1_filtered(nshape_h1-1, 1, 0.);
  TPZFMatrix<REAL> dphi_h1(dim, nshape_h1, 0.);
  typename TSHAPE::IntruleType intrule(kFacet+6);
  int nintpoints = intrule.NPoints();
  TPZManVector<REAL, 3> point(dim, 0.);

  TPZFMatrix<REAL> Mstd(nshape_std, nshape_std, 0.);
  TPZFMatrix<REAL> Mcust(nshape_cust, nshape_cust, 0.);
  TPZFMatrix<REAL> Mstd_cust(nshape_std, nshape_cust, 0.);
  TPZFMatrix<REAL> Mdiv(nshape_cust, nshape_cust, 0.);
  TPZFMatrix<REAL> Mh1(nshape_h1, nshape_h1, 0.);
  TPZFMatrix<REAL> Mdiv_h1(nshape_cust, nshape_h1, 0.);

  TPZFMatrix<REAL> Mdiv_filtered(nFiltered, nFiltered, 0.);
  TPZFMatrix<REAL> Mdiv_filtered_h1(nFiltered, nshape_h1, 0.);
  TPZFMatrix<REAL> Mh1_filtered(nshape_h1-1, nshape_h1-1, 0.);

  for (int ip = 0; ip < nintpoints; ip++) {
    REAL weight;
    intrule.Point(ip, point, weight);
    hdiv_std.Shape(point, dataHDiv, phi_std, div_std);
    h1.Shape(point, dataH1, phi_h1, dphi_h1);

    // Filtered functions

    for (int i = 0; i < nFiltered; i++) {
      div_filtered(i, 0) = div_std(filteredIndices[i], 0);
    }

    for (int i =1; i < nshape_h1; i++) {
      phi_h1_filtered(i-1,0) = phi_h1(i,0);
    }

    // Matrices to check deRham compatibility
    Mdiv.AddContribution(0, 0, div_std, 0, div_std, 1, weight); // div x div
    Mh1.AddContribution(0, 0, phi_h1, 0, phi_h1, 1, weight); // L2 x L2
    Mdiv_h1.AddContribution(0, 0, div_std, 0, phi_h1, 1, weight); // div x L2

    // Matrices to check deRham compatibility of the filtered functions
    Mdiv_filtered.AddContribution(0, 0, div_filtered, 0, div_filtered, 1, weight); // div x div
    Mdiv_filtered_h1.AddContribution(0, 0, div_filtered, 0, phi_h1_filtered, 1, weight); // div x L2
    Mh1_filtered.AddContribution(0, 0, phi_h1_filtered, 0, phi_h1_filtered, 1, weight); // L2 x L2
  }

  {
    std::ofstream out("divergenceMatrixFull.txt");
    Mdiv.Print("GK = ", out, EMathematicaInput);
  }

  {
    std::ofstream out("divergenceMatrixFiltered.txt");
    Mdiv_filtered.Print("GKFilt = ", out, EMathematicaInput);
  }
  
  // TODO: Check with gi. Maybe cleaner way to do it
  TPZFMatrix<REAL> Mdiv_h1_t;
  Mh1.Decompose_LU();
  Mdiv_h1.Transpose(&Mdiv_h1_t);
  Mh1.Substitution(&Mdiv_h1_t); 
  TPZFMatrix<REAL> Res1 = Mdiv - Mdiv_h1 * Mdiv_h1_t;

  TPZFMatrix<REAL> Mdiv_filtered_h1_t;
  Mh1_filtered.Decompose_LU();
  Mdiv_filtered_h1.Transpose(&Mdiv_filtered_h1_t);
  Mh1_filtered.Substitution(&Mdiv_filtered_h1_t); 
  TPZFMatrix<REAL> Res2 = Mdiv_filtered - Mdiv_filtered_h1 * Mdiv_filtered_h1_t;

  std::cout << "\n------------------------------" << std::endl;
  std::cout << "Calling the test for " << typeid(TSHAPE).name() << " with order " << kFacet << std::endl;
  std::cout << "\nNumber of shape functions: " << nshape_std << std::endl;
  std::cout << "Number of internal shape functions: " << nshapeInt_std << std::endl;
  std::cout << "Number of non-null divergence shape functions: " << filteredIndices.size() << std::endl;
  std::cout << "\nIs de Rham compatible? " << Res1.MatrixNorm(1, 10, 1.e-10) << std::endl;
  std::cout << "Is de Rham compatible for the filtered functions? " << Res2.MatrixNorm(1, 10, 1.e-10) << std::endl;
  std::cout << "------------------------------" << std::endl;
}

template <class TSHAPE>
void filterEquations(TPZVec<int> &filteredIndices, int kFacet) {
  int nFiltered = 0;
  int count = 0;
  if (TSHAPE::Type() == ECube) {
    filteredIndices.Resize((kFacet+1)*(kFacet+1)*(kFacet+1)-1);
    filteredIndices.Fill(-1);
    count = TSHAPE::NFacets * (kFacet + 1) * (kFacet + 1); // Skip trace functions
    std::cout << "Count after trace functions: " << count << std::endl;

    // Pick edge functions of the first 8 edges, skipping  the 4th one.
    for (int i = 0; i < 8; i++) {
      if (i == 3) {
        count += kFacet;
        continue;
      }
      for (int j = 0; j < kFacet; j++) {
        filteredIndices[nFiltered] = count;
        nFiltered++;
        count++;
      }
    }

    count += 4 * kFacet; // Skip the four remaining edges

    std::cout << "Count after edges: " << count << "   nFiltered after edges: " << nFiltered << std::endl;

    // Pick the functions for the second direction of the faces 1,2,3, and 4
    count += 2 * kFacet * (kFacet - 1); // Skip the first face
    for (int i = 0; i < 4; i++) {
      count += kFacet * (kFacet - 1); // Skip the first direction of the face
      for (int j = 0; j < kFacet * (kFacet - 1); j++) {
        filteredIndices[nFiltered] = count;
        nFiltered++;
        count++;
      }
    }

    // Pick the functions of the first direction of the faces 6
    for (int j = 0; j < kFacet * (kFacet - 1); j++) {
      filteredIndices[nFiltered] = count;
      nFiltered++;
      count++;
    }

    count += kFacet * (kFacet - 1); // Skip the second direction of the face 6

    std::cout << "Count after faces: " << count << "  nFiltered after faces: " << nFiltered << std::endl;

    // Pick the internal functions in the third direction
    count += 2 * (kFacet * (kFacet - 1) * (kFacet - 1));
    for (int j = 0; j < kFacet * (kFacet - 1) * (kFacet - 1); j++) {
      filteredIndices[nFiltered] = count;
      nFiltered++;
      count++;
    }
  } else {
    std::cout << "Not implemented yet" << std::endl;
    return;
  }
  std::cout << "Filtered indices: ";
  for (int i = 0; i < nFiltered; i++) {
    std::cout << filteredIndices[i] << " ";
  }
  std::cout << std::endl;
}