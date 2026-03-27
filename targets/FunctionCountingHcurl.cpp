#include <iostream>

#include "TPZMaterialDataT.h"
#include "TPZShapeHDiv.h"
#include "TPZShapeHDivConstant.h"
#include "TPZShapeHDivOptimized.h"
#include "TPZShapeHCurl.h"
#include "TPZShapeHCurlNoGrads.h"

#include "pzmanvector.h"
#include "pzshapeprism.h"

// Internal stuff
#include "TPZShapeHDivRefactor.h"

enum EFamily {EStandard, EConstant, EOptimized, ERefactor, ECurl, ECurlNoGrads};

template <class TSHAPE>
void CountingFunctions(int kFacet, int Family);

int main(int argc, char **argv) {

  for (int kFacet = 1; kFacet <= 3; kFacet++) {
    CountingFunctions<pzshape::TPZShapePrism>(kFacet, ECurl);
  }
}

template <class TSHAPE>
void CountingFunctions(int kFacet, int Family) {
  TPZManVector<int, 27> sideorient(TSHAPE::NFacets, 1);
  TPZManVector<int, 27> ordersHdiv(TSHAPE::NFacets + 1, kFacet);
  TPZManVector<int, 27> ordersHcurl(TSHAPE::NSides-TSHAPE::NCornerNodes, kFacet);
  TPZManVector<int64_t, 27> ids(TSHAPE::NCornerNodes, 0);
  for (int i = 0; i < TSHAPE::NCornerNodes; i++) {
    ids[i] = i;
  }

  TPZMaterialDataT<REAL> data;

  // Initialize the shape with the appropriate family
  if (Family == EStandard) {
    TPZShapeHDiv<TSHAPE> shape;
    shape.Initialize(ids, ordersHdiv, sideorient, data);
    std::cout << "\nStandard Hdiv of order " << kFacet <<":"<< std::endl;
    std::cout << "Number of shape functions: " << shape.NShapeF(data) << std::endl;
    std::cout << "Number of shapes per connect: ";
    int nc = data.fHDiv.fNumConnectShape.size();
    for (int i = 0; i < nc; i++) {
      std::cout << shape.NConnectShapeF(i, data) << " ";
    }
    std::cout << std::endl;
  } 
  else if (Family == EConstant) {
    TPZShapeHDivConstant<TSHAPE> shape;
    shape.Initialize(ids, ordersHdiv, sideorient, data);
    std::cout << "\nConstant Hdiv of order " << kFacet <<":"<< std::endl;
    std::cout << "Number of shape functions: " << shape.NHDivShapeF(data) << std::endl;
    std::cout << "Number of shapes per connect: ";
    int nc = data.fHDiv.fNumConnectShape.size();
    for (int i = 0; i < nc; i++) {
      std::cout << shape.NConnectShapeF(i, data) << " ";
    }
    std::cout << std::endl;
  }
  else if (Family == EOptimized) {
    TPZShapeHDivOptimized<TSHAPE> shape;
    shape.Initialize(ids, ordersHdiv, sideorient, data);
    std::cout << "\nOptimized Hdiv of order " << kFacet <<":"<< std::endl;
    std::cout << "Number of shape functions: " << shape.NShapeF(data) << std::endl;
    std::cout << "Number of shapes per connect: ";
    int nc = data.fHDiv.fNumConnectShape.size();
    for (int i = 0; i < nc; i++) {
      std::cout << shape.NConnectShapeF(i, data) << " ";
    }
    std::cout << std::endl;
  }
  else if (Family == ERefactor) {
    TPZShapeHDivRefactor<TSHAPE> shape;
    shape.Initialize(ids, ordersHdiv, sideorient, data);
    std::cout << "\nRefactor Hdiv of order " << kFacet <<":"<< std::endl;
    std::cout << "Number of shape functions: " << shape.NShapeF(data) << std::endl;
    std::cout << "Number of shapes per connect: ";
    int nc = data.fHDiv.fNumConnectShape.size();
    for (int i = 0; i < nc; i++) {
      std::cout << shape.NConnectShapeF(i, data) << " ";
    }
    std::cout << std::endl;
  }
  else if (Family == ECurl) {
    TPZShapeHCurl<TSHAPE> shape;
    shape.Initialize(ids, ordersHcurl, data);
    // std::cout << data.fHCurl.fSDVecShapeIndex << std::endl;
    std::cout << "\nHCurl of order " << kFacet <<":"<< std::endl;
    std::cout << "Number of shape functions: " << shape.NHCurlShapeF(data) << std::endl;
    std::cout << "Number of shapes per connect: ";
    int nc = data.fHCurl.fNumConnectShape.size();
    for (int i = 0; i < nc; i++) {
      std::cout << data.fHCurl.fNumConnectShape[i] << " ";
    }
    std::cout << std::endl;
  }
  else if (Family == ECurlNoGrads) {
    TPZShapeHCurlNoGrads<TSHAPE> shape;
    shape.Initialize(ids, ordersHcurl, data);
    std::cout << "\nHCurl (No Grads) of order " << kFacet <<":"<< std::endl;
    std::cout << "Number of shape functions: " << shape.NHCurlShapeF(data) << std::endl;
    std::cout << "Number of shapes per connect: ";
    int nc = data.fHCurl.fNumConnectShape.size();
    for (int i = 0; i < nc; i++) {
      std::cout << data.fHCurl.fNumConnectShape[i] << " ";
    }
    std::cout << std::endl;
  }
}