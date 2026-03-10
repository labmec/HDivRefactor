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

int main(int argc, char **argv) {
  int kFacet = 2;
  TestHDivOptimized<pzshape::TPZShapeQuad>(kFacet);
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
  int nshape_opt = hdiv_opt.NShapeF(dataHDivOptimized);
  int nshape_const = hdiv_const.NHDivShapeF(dataHDivConst);
  int nshape_cust = nshape_opt;
  int nshape_h1 = 9; //h1.NShapeF(dataH1);

  int nfacet_std = 0, nfacet_opt = 0, nfacet_const = 0;
  for (int i = 0; i < TSHAPE::NFacets; i++) {
    nfacet_std += hdiv_std.NConnectShapeF(i, dataHDiv);
    nfacet_opt += hdiv_opt.NConnectShapeF(i, dataHDivOptimized);
    nfacet_const += hdiv_const.NConnectShapeF(i, dataHDivConst);
  }
  
  int nvol_std = hdiv_std.NConnectShapeF(TSHAPE::NFacets, dataHDiv);
  int nvol_opt = hdiv_opt.NConnectShapeF(TSHAPE::NFacets, dataHDivOptimized);
  int nvol_const = hdiv_const.NConnectShapeF(TSHAPE::NFacets, dataHDivConst);

  constexpr int dim = TSHAPE::Dimension;
  TPZFMatrix<REAL> phi_std(dim, nshape_std, 0.), phi_opt(dim, nshape_opt, 0.);
  TPZFMatrix<REAL> phi_const(dim, nshape_const, 0.), phi_cust(dim, nshape_cust, 0.);
  TPZFNMatrix<60, REAL> div_std(nshape_std, 1), div_opt(nshape_opt, 1);
  TPZFNMatrix<60, REAL> div_const(nshape_const, 1), div_cust(nshape_cust, 1);
  TPZFMatrix<REAL> phi_h1(nshape_h1, 1, 0.);
  TPZFMatrix<REAL> dphi_h1(dim, nshape_h1, 0.);
  typename TSHAPE::IntruleType intrule(3);
  int nintpoints = intrule.NPoints();
  TPZManVector<REAL, 3> point(dim, 0.);

  TPZFMatrix<REAL> Mstd(nshape_std, nshape_std, 0.);
  TPZFMatrix<REAL> Mcust(nshape_cust, nshape_cust, 0.);
  TPZFMatrix<REAL> Mstd_cust(nshape_std, nshape_cust, 0.);
  TPZFMatrix<REAL> Mdiv(nshape_cust, nshape_cust, 0.);
  TPZFMatrix<REAL> Mh1(nshape_h1, nshape_h1, 0.);
  TPZFMatrix<REAL> Mdiv_h1(nshape_cust, nshape_h1, 0.);

  for (int ip = 0; ip < nintpoints; ip++) {
    REAL weight;
    intrule.Point(ip, point, weight);
    hdiv_std.Shape(point, dataHDiv, phi_std, div_std);
    hdiv_opt.Shape(point, dataHDivOptimized, phi_opt, div_opt);
    hdiv_const.Shape(point, dataHDivConst, phi_const, div_const);
    h1.Shape(point, dataH1, phi_h1, dphi_h1);

    // TODO: Build phi_cust and div_cust here! Maybe using a function
    phi_cust = phi_opt;
    div_cust = div_opt;

    // Matrices to check deRham compatibility
    Mdiv.AddContribution(0, 0, div_cust, 0, div_cust, 1, weight); // div x div
    Mh1.AddContribution(0, 0, phi_h1, 0, phi_h1, 1, weight); // L2 x L2
    Mdiv_h1.AddContribution(0, 0, div_cust, 0, phi_h1, 1, weight); // div x L2

    // Matrices to check polynomial space
    Mstd.AddContribution(0, 0, phi_std, 1, phi_std, 0, weight); // std x std
    Mcust.AddContribution(0, 0, phi_cust, 1, phi_cust, 0, weight); // cust x cust
    Mstd_cust.AddContribution(0, 0, phi_std, 1, phi_cust, 0, weight); // std x cust
  }

  // Check with gi. Maybe cleaner way to do it
  TPZFMatrix<REAL> Mdiv_h1_t;
  Mh1.Decompose_LU();
  Mdiv_h1.Transpose(&Mdiv_h1_t);
  Mh1.Substitution(&Mdiv_h1_t); 
  TPZFMatrix<REAL> Res1 = Mdiv - Mdiv_h1 * Mdiv_h1_t;

  TPZFMatrix<REAL> Mstd_cust_t;
  Mcust.Decompose_LU();
  Mstd_cust.Transpose(&Mstd_cust_t);
  Mcust.Substitution(&Mstd_cust_t); 
  TPZFMatrix<REAL> Res2 = Mstd - Mstd_cust * Mstd_cust_t;

  std::cout << "\nIs deRham compatible? " << Res1.MatrixNorm(2, 10, 1.e-10) << std::endl;
  std::cout << "Does it span the same polynomial space? " << Res2.MatrixNorm(2, 10, 1.e-10) << "\n" << std::endl;
}