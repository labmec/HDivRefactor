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

// Internal stuff
#include "TPZShapeHDivRefactor.h"

template <class TSHAPE>
void TestHDivOptimized(int kFacet);

template <class TSHAPE>
void filterEquations(TPZVec<int> &filteredIndices, int kFacet);

int main(int argc, char **argv) {
  for (int kFacet = 1; kFacet <= 3; kFacet++) {
    TestHDivOptimized<pzshape::TPZShapeQuad>(kFacet);
    TestHDivOptimized<pzshape::TPZShapeCube>(kFacet);
    TestHDivOptimized<pzshape::TPZShapeTriang>(kFacet);
    TestHDivOptimized<pzshape::TPZShapeTetra>(kFacet);
  }
}

template <class TSHAPE>
void TestHDivOptimized(int kFacet) {
  TPZMaterialDataT<REAL> dataHDiv, dataHDivOptimized, dataHDivRefactor, dataHDivConst, dataH1;
  TPZManVector<int64_t, 27> ids(TSHAPE::NCornerNodes, 0);
  for (int i = 0; i < TSHAPE::NCornerNodes; i++) {
    ids[i] = i;
  }
  TPZManVector<int, 27> sideorient(TSHAPE::NFacets, 1);
  TPZManVector<int, 27> orders(TSHAPE::NFacets + 1, kFacet);
  TPZManVector<int, 27> h1orders(TSHAPE::NSides - TSHAPE::NCornerNodes, kFacet);
  TPZShapeHDiv<TSHAPE> hdiv_std;
  TPZShapeHDivOptimized<TSHAPE> hdiv_opt;
  TPZShapeHDivRefactor<TSHAPE> hdiv_rft;
  TPZShapeHDivConstant<TSHAPE> hdiv_const;
  TPZShapeH1<TSHAPE> h1;

  hdiv_std.Initialize(ids, orders, sideorient, dataHDiv);
  hdiv_rft.Initialize(ids, orders, sideorient, dataHDivRefactor);
  hdiv_opt.Initialize(ids, orders, sideorient, dataHDivOptimized);
  hdiv_const.Initialize(ids, orders, sideorient, dataHDivConst);
  h1.Initialize(ids, h1orders, dataH1);
  int nshape_std = hdiv_std.NShapeF(dataHDiv);
  int nshape_rft = hdiv_rft.NShapeF(dataHDivRefactor);
  int nshape_opt = hdiv_opt.NShapeF(dataHDivOptimized);
  int nshapeInt_std = hdiv_std.NConnectShapeF(TSHAPE::NFacets, dataHDiv);
  int nshape_const = hdiv_const.NHDivShapeF(dataHDivConst);
  int nshape_h1 = dataH1.fH1.fPhi.Rows();

  constexpr int dim = TSHAPE::Dimension;
  int nFacets = TSHAPE::NFacets;
  TPZFMatrix<REAL> phi_std(dim, nshape_std, 0.);
  TPZFNMatrix<60, REAL> div_std(nshape_std, 1);
  TPZFMatrix<REAL> phi_rft(dim, nshape_rft, 0.);
  TPZFNMatrix<60, REAL> div_rft(nshape_rft, 1);
  TPZFMatrix<REAL> phi_h1(nshape_h1, 1, 0.);
  TPZFMatrix<REAL> dphi_h1(dim, nshape_h1, 0.);
  typename TSHAPE::IntruleType intrule(kFacet+6);
  int nintpoints = intrule.NPoints();
  TPZManVector<REAL, 3> point(dim, 0.);

  TPZFMatrix<REAL> Mdiv_rft(nshape_rft, nshape_rft, 0.);
  TPZFMatrix<REAL> Mdiv_rft_h1(nshape_rft, nshape_h1, 0.);
  TPZFMatrix<REAL> Mh1(nshape_h1, nshape_h1, 0.);

  TPZFMatrix<REAL> Mphi_std(nshape_std, nshape_std, 0.);
  TPZFMatrix<REAL> Mphi_rft(nshape_rft, nshape_rft, 0.);
  TPZFMatrix<REAL> Mphi_std_rft(nshape_std, nshape_rft, 0.);

  for (int ip = 0; ip < nintpoints; ip++) {
    REAL weight;
    intrule.Point(ip, point, weight);
    hdiv_std.Shape(point, dataHDiv, phi_std, div_std);
    hdiv_rft.Shape(point, dataHDivRefactor, phi_rft, div_rft);
    h1.Shape(point, dataH1, phi_h1, dphi_h1);

    {
      std::ofstream out("shapeFunctionsRefactor.txt");
      phi_rft.Print("GKRefactor = ", out, EMathematicaInput);
    }

    // Matrices to check deRham compatibility
    Mdiv_rft.AddContribution(0, 0, div_rft, 0, div_rft, 1, weight); // div x div
    Mh1.AddContribution(0, 0, phi_h1, 0, phi_h1, 1, weight); // L2 x L2
    Mdiv_rft_h1.AddContribution(0, 0, div_rft, 0, phi_h1, 1, weight); // div x L2

    // Matrices to check polynomial span of the refactor functions
    Mphi_std.AddContribution(0, 0, phi_std, 1, phi_std, 0, weight); // std x std
    Mphi_rft.AddContribution(0, 0, phi_rft, 1, phi_rft, 0, weight); // rft x rft
    Mphi_std_rft.AddContribution(0, 0, phi_std, 1, phi_rft, 0, weight); // std x rft
  }

  {
    std::ofstream out("divergenceMatrixRefactor.txt");
    Mdiv_rft.Print("GKRefactor = ", out, EMathematicaInput);
  }

  {
    std::ofstream out("L2MatrixRefactor.txt");
    Mphi_rft.Print("GK = ", out, EMathematicaInput);
  }
  
  // Span test
  TPZFMatrix<REAL> Mphi_std_rft_t;
  Mphi_std_rft.Transpose(&Mphi_std_rft_t);
  Mphi_rft.SolveDirect(Mphi_std_rft_t, ELU); 
  TPZFMatrix<REAL> Res1 = Mphi_std - Mphi_std_rft * Mphi_std_rft_t;

  // de Rham test
  TPZFMatrix<REAL> Mdiv_rft_h1_t;
  Mdiv_rft_h1.Transpose(&Mdiv_rft_h1_t);
  Mh1.SolveDirect(Mdiv_rft_h1_t, ELDLt); 
  TPZFMatrix<REAL> Res2 = Mdiv_rft - Mdiv_rft_h1 * Mdiv_rft_h1_t;

  std::cout << "\n------------------------------" << std::endl;
  std::cout << "Calling the test for " << typeid(TSHAPE).name() << " with order " << kFacet << std::endl;
  std::cout << "\nNumber of shape functions: " << nshape_std << std::endl;
  std::cout << "Number of internal shape functions: " << nshapeInt_std << std::endl;
  std::cout << "\nDoes it expands the same polynomial space? " << Res1.MatrixNorm(1, 10, 1.e-10) << std::endl;
  std::cout << "Is de Rham compatible? " << Res2.MatrixNorm(1, 10, 1.e-10) << std::endl;
  std::cout << "------------------------------" << std::endl;

  if (Res1.MatrixNorm(1, 10, 1.e-10) > 1.e-12) {
    std::cout << "\nWarning: Refactor HDiv does not seem to be expanding the same polynomial space! Norm of the residual: " << Res1.MatrixNorm(1, 10, 1.e-10) << std::endl;
    // DebugStop();
  }

  if (Res2.MatrixNorm(1, 10, 1.e-10) > 1.e-12) {
    std::cout << "\nWarning: Refactor HDiv does not seem to be de Rham compatible! Norm of the residual: " << Res2.MatrixNorm(1, 10, 1.e-10) << std::endl;
    // DebugStop();
  }
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

  } else if (TSHAPE::Type() == EQuadrilateral) {
    filteredIndices.Resize((kFacet+1)*(kFacet+1)-1);
    filteredIndices.Fill(-1);
    count = TSHAPE::NFacets * (kFacet + 1); // Skip trace functions
    std::cout << "Count after trace functions: " << count << std::endl;

    // Pick edge functions of the first 3 edges
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < kFacet; j++) {
        filteredIndices[nFiltered] = count;
        nFiltered++;
        count++;
      }
    }

    count += kFacet; // Skip the final edge

    std::cout << "Count after edges: " << count << "   nFiltered after edges: " << nFiltered << std::endl;

    // Pick the functions of the first direction of the face
    for (int j = 0; j < kFacet * (kFacet - 1); j++) {
      filteredIndices[nFiltered] = count;
      nFiltered++;
      count++;
    }

    std::cout << "Count after faces: " << count << "  nFiltered after faces: " << nFiltered << std::endl;

  } else if (TSHAPE::Type() == ETriangle) {
    filteredIndices.Resize((kFacet+1)*(kFacet+2)/2-1);
    filteredIndices.Fill(-1);
    count = TSHAPE::NFacets * (kFacet + 1); // Skip trace functions
    std::cout << "Count after trace functions: " << count << std::endl;

    // Pick edge functions of the first 2 edges
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < kFacet; j++) {
        filteredIndices[nFiltered] = count;
        nFiltered++;
        count++;
      }
    }

    count += kFacet; // Skip the final edge

    std::cout << "Count after edges: " << count << "   nFiltered after edges: " << nFiltered << std::endl;

    // Pick the functions of the first direction of the face
    for (int j = 0; j < kFacet * (kFacet - 1) / 2; j++) {
      filteredIndices[nFiltered] = count;
      nFiltered++;
      count++;
    }

    std::cout << "Count after faces: " << count << "  nFiltered after faces: " << nFiltered << std::endl;

  } else if (TSHAPE::Type() == ETetraedro) {
    filteredIndices.Resize((kFacet+3)*(kFacet+2)*(kFacet+1)/6-1);
    filteredIndices.Fill(-1);
    count = TSHAPE::NFacets * (kFacet + 2) * (kFacet + 1) / 2; // Skip trace functions
    std::cout << "Count after trace functions: " << count << std::endl;

    // Pick edge functions of the first 3 edges, skipping  the 2th one.
    for (int i = 0; i < 4; i++) {
      if (i == 2) {
        count += kFacet;
        continue;
      }
      for (int j = 0; j < kFacet; j++) {
        filteredIndices[nFiltered] = count;
        nFiltered++;
        count++;
      }
    }

    count += 2*kFacet; // Skip the two remaining edges

    std::cout << "Count after edges: " << count << "   nFiltered after edges: " << nFiltered << std::endl;

    // Pick the functions for the second direction of the faces 1,2,3
    count += kFacet * (kFacet - 1); // Skip the first face
    for (int i = 0; i < 3; i++) {
      count += kFacet * (kFacet - 1)/2; // Skip the first direction of the face
      for (int j = 0; j < kFacet * (kFacet - 1)/2; j++) {
        filteredIndices[nFiltered] = count;
        nFiltered++;
        count++;
      }
    }

    // Pick the functions of the first direction of the faces 6
    count += kFacet * (kFacet - 1) * (kFacet -2)/3; // Skipe the first two directions of the volume
    for (int j = 0; j < kFacet * (kFacet - 1) * (kFacet - 2) / 6; j++) {
      filteredIndices[nFiltered] = count;
      nFiltered++;
      count++;
    }

  } else {
    std::cout << "Not implemented yet" << std::endl;
    return;
  }

  // Print the filtered indices
  std::cout << "Filtered indices: ";
  for (int i = 0; i < nFiltered; i++) {
    std::cout << filteredIndices[i] << " ";
  }
  std::cout << std::endl;
}