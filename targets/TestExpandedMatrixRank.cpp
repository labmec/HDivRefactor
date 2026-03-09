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

int main(int argc, char *argv) {
}

template <class TSHAPE>
void TestHDivOptimized(int kFacet) {
  TPZMaterialDataT<REAL> dataHDiv, dataHDivOptimized, dataHDivConst;
  TPZManVector<int64_t, 27> ids(TSHAPE::NCornerNodes, 0);
  for (int i = 0; i < TSHAPE::NCornerNodes; i++) {
    ids[i] = i;
  }
  TPZManVector<int, 27> sideorient(TSHAPE::NFacets, 1);
  TPZManVector<int, 27> orders(TSHAPE::NFacets + 1, kFacet);
  TPZShapeHDiv<TSHAPE> hdiv_std;
  TPZShapeHDivOptimized<TSHAPE> hdiv_opt;
  TPZShapeHDivConstant<TSHAPE> hdiv_const;
  hdiv_std.Initialize(ids, orders, sideorient, dataHDiv);
  hdiv_opt.Initialize(ids, orders, sideorient, dataHDivOptimized);
  hdiv_const.Initialize(ids, orders, sideorient, dataHDivConst);
  int nshape_std = hdiv_std.NShapeF(dataHDiv);
  int nshape_opt = hdiv_opt.NShapeF(dataHDivOptimized);
  int nshape_const = hdiv_const.NHDivShapeF(dataHDivConst);
  int nfacet_std = 0, nfacet_opt = 0, nfacet_const = 0;
  for (int i = 0; i < TSHAPE::NFacets; i++) {
    nfacet_std += hdiv_std.NConnectShapeF(i, dataHDiv);
    nfacet_opt += hdiv_opt.NConnectShapeF(i, dataHDivOptimized);
    nfacet_const += hdiv_const.NConnectShapeF(i, dataHDivConst);
  }
  int nvol_std = hdiv_std.NConnectShapeF(TSHAPE::NFacets, dataHDiv);
  int nvol_opt = hdiv_opt.NConnectShapeF(TSHAPE::NFacets, dataHDivOptimized);
  int nvol_const = hdiv_const.NConnectShapeF(TSHAPE::NFacets, dataHDivConst);

  auto hdiv_std_order = dataHDiv.fHDiv.fConnectOrders;
  auto hdiv_opt_order = dataHDivOptimized.fHDiv.fConnectOrders;
  auto hdiv_const_order = dataHDivConst.fHDiv.fConnectOrders;

  auto hdiv_std_shape = dataHDiv.fHDiv.fNumConnectShape;
  auto hdiv_opt_shape = dataHDivOptimized.fHDiv.fNumConnectShape;
  auto hdiv_const_shape = dataHDivConst.fHDiv.fNumConnectShape;

  CAPTURE(orders, hdiv_std_order, hdiv_opt_order, hdiv_const_order, hdiv_std_shape, hdiv_opt_shape, hdiv_const_shape);

  REQUIRE(nvol_std == nvol_opt);
  REQUIRE(nfacet_opt == nfacet_const);
  REQUIRE(nfacet_std + nvol_std == nshape_std);
  REQUIRE(nfacet_opt + nvol_opt == nshape_opt);

  for (int i = 0; i < TSHAPE::NFacets + 1; i++) {
    REQUIRE(hdiv_std_order[i] == orders[i]);
    REQUIRE(hdiv_opt_order[i] == orders[i]);
    REQUIRE(hdiv_std_shape[i] == hdiv_std.NConnectShapeF(i, dataHDiv));
    REQUIRE(hdiv_opt_shape[i] == hdiv_opt.NConnectShapeF(i, dataHDivOptimized));
    REQUIRE(hdiv_const_shape[i] == hdiv_const.NConnectShapeF(i, dataHDivConst));
  }

  // Evaluate the shape functions at the integration point.
  // The volume shape functions must be equal for hdiv_std and hdiv_opt at every integration point.
  // The facet shape functions must be equal for hdiv_opt and hdiv_const at every integration point.
  constexpr int dim = TSHAPE::Dimension;
  TPZFMatrix<REAL> phi_std(dim, nshape_std, 0.), phi_opt(dim, nshape_opt, 0.), phi_const(dim, nshape_const, 0.);
  TPZFNMatrix<60, REAL> div_std(nshape_std, 1), div_opt(nshape_opt, 1), div_const(nshape_const, 1);
  typename TSHAPE::IntruleType intrule(3);
  int nintpoints = intrule.NPoints();
  TPZManVector<REAL, 3> point(dim, 0.);
  for (int ip = 0; ip < nintpoints; ip++) {
    REAL weight;
    intrule.Point(ip, point, weight);
    hdiv_std.Shape(point, dataHDiv, phi_std, div_std);
    hdiv_opt.Shape(point, dataHDivOptimized, phi_opt, div_opt);
    hdiv_const.Shape(point, dataHDivConst, phi_const, div_const);

    for (int i = 0; i < nvol_std; i++) {
      for (int j = 0; j < dim; j++) {
        REAL val_std = phi_std(j, i + nfacet_std);
        REAL val_opt = phi_opt(j, i + nfacet_opt);
        if (abs(val_std - val_opt) > 1.e-10) {
          std::cout << "val_std " << val_std << " val_opt " << val_opt << std::endl;
        }
        REQUIRE((val_std - val_opt) == Catch::Approx(0.));
      }
    }

    for (int i = 0; i < nfacet_const; i++) {
      for (int j = 0; j < dim; j++) {
        REAL val_const = phi_const(j, i);
        REAL val_opt = phi_opt(j, i);
        if (abs(val_const - val_opt) > 1.e-10) {
          std::cout << "val_const " << val_const << " val_opt " << val_opt << std::endl;
        }
        REQUIRE((val_const - val_opt) == Catch::Approx(0.));
      }
    }
  }
}