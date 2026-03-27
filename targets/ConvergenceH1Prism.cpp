#include <iostream>
#include "TPZGenGrid3D.h"
#include "TPZRefPatternDataBase.h"
#include "TPZVTKGeoMesh.h"
#include "TPZVTKGenerator.h"
#include "TPZNullMaterial.h"
#include "DarcyFlow/TPZDarcyFlow.h"
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZLinearAnalysis.h"
#include "TPZSSpStructMatrix.h"
#include "pzstepsolver.h"
#include "TPZAnalyticSolution.h"
#include "pzlog.h"

// ----------------
// Global variables
// ----------------

TLaplaceExample1 gexact;

enum EnumMatIds {
  EDomain = 1,  // Material ID for the domain
  EBottom = 2, // Material ID bottom boundary
  EFront = 3,  // Material ID front boundary
  ERight = 4,  // Material ID right boundary
  EBack = 5,   // Material ID back boundary
  ELeft = 6,    // Material ID left boundary
  ETop = 7,    // Material ID top boundary
};

// Creates a geometric mesh using TPZGenGrid3D
TPZGeoMesh* createGeoMesh(const TPZManVector<int, 3> &nelDiv, 
  const TPZManVector<REAL, 3> &minX, const TPZManVector<REAL, 3> &maxX);

// Computes the diameter of a geometric element
REAL ElementDiameter(TPZGeoEl* gel);

// Computes the diameter of a mesh  
REAL MeshDiameter(TPZGeoMesh *gmesh);

// Creates a computational mesh for H1 approximation
TPZCompMesh* createCompMeshH1(TPZGeoMesh *gmesh, int order = 1);

// Computes and prints the convergence order of relevant errors
void ConvOrder(TPZFMatrix<REAL> &errorMat, TPZVec<REAL> &hvector, int nref);

// ----
// Main
// ----

int main (int argc, char * const argv[]) {

  // Initialize the logger (How it works?)
  #ifdef PZ_LOG
  TPZLogger::InitializePZLOG();
  #endif

  // Initializing uniform refinements for reference elements
  // (Used in TPZGenGrid2D?)
  gRefDBase.InitializeUniformRefPattern(EOned);
  gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
  gRefDBase.InitializeUniformRefPattern(ETriangle);
  gRefDBase.InitializeUniformRefPattern(EPrisma);

  const int nthreads = 0;
  const int refMax = 5;

  TPZGeoMesh* gmesh = nullptr;
  TPZCompMesh* cmesh = nullptr;

  TPZFMatrix<REAL> errorMat(refMax, 5, 0.); // To storage the errors
  TPZVec<REAL> hvector(refMax, 0.); // To storage the mesh size

  int order = 2; // Approximation order

  // Set a problem with analytic solution
  gexact.fExact = TLaplaceExample1::ESinSin; 

  // Initial geometric mesh
  gmesh = createGeoMesh({2, 2, 2}, {0., 0., 0.}, {1., 1., 1.});
  // TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);

  // Refinement loop
  for (int ref = 0; ref < refMax; ++ref) {

    // Mesh refinement
    if (gmesh && ref > 0) {
      TPZCheckGeom checkgeom(gmesh);
      checkgeom.UniformRefine(1);
    }

    // Get mesh diameter and store it
    hvector[ref] = MeshDiameter(gmesh);

    // ------ Create computational mesh ------

    if (cmesh) {
      delete cmesh; // Delete previous computational mesh (needed?)
    }

    cmesh = createCompMeshH1(gmesh, order);


    // ------ Assemble and solve ------

    TPZLinearAnalysis an(cmesh); // Analysis object
    TPZSSpStructMatrix<STATE> matsp(cmesh); // Sparse matrix structure for assembly
    matsp.SetNumThreads(nthreads); // Number of threads for assembly
    an.SetStructuralMatrix(matsp);

    // Set direct solver
    TPZStepSolver<STATE> step;
    auto directType = ECholesky;
    step.SetDirect(directType);
    an.SetSolver(step);
    an.Run();

    // an.Solution().Print("Solution");
    // an.Rhs().Print("Rhs")

    // ------ Compute errors ------

    TPZVec<REAL> errors(3, 0.);
    an.PostProcessError(errors, false, std::cout);

    int nerrors = 3; // Number of computed errors

    // Store errors in the matrix
    for(int j=0; j<nerrors; j++){
        errorMat(ref, j) = errors[j];
    }
  }

  // ------ Convergence orders ------
  ConvOrder(errorMat, hvector, refMax);

  // ------ Plotting ------
  // (Only the most refined solution)

  const std::string plotfile = "darcyplot";  // sem o .vtk no final
  constexpr int vtkRes{0};

  TPZManVector<std::string, 2> fields = {"Flux", "Pressure"};
  auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);

  vtk.Do();

  // ------ Clean up ------

  delete cmesh, 
         gmesh;
  
  return 0;
}

// ---------
// Functions
// ---------

TPZGeoMesh* createGeoMesh(
  const TPZManVector<int, 3> &nelDiv, 
  const TPZManVector<REAL, 3> &minX, 
  const TPZManVector<REAL, 3> &maxX) {

  TPZGenGrid3D generator(minX, maxX, nelDiv, MMeshType::EPrismatic);

  generator.BuildVolumetricElements(EDomain);
  TPZGeoMesh *gmesh = generator.BuildBoundaryElements(EBottom, EFront, ERight, EBack, ELeft, ETop);

  return gmesh;
}

TPZCompMesh* createCompMeshH1(TPZGeoMesh *gmesh, int order) {
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDimModel(gmesh->Dimension());
  cmesh->SetDefaultOrder(order); // Polynomial order
  cmesh->SetAllCreateFunctionsContinuous(); // H1 Elements

  // Add materials (weak formulation)
  TPZDarcyFlow *mat = new TPZDarcyFlow(EDomain, gmesh->Dimension());
  mat->SetConstantPermeability(1.0); // Set constant permeability
  mat->SetForcingFunction(gexact.ForceFunc(),4);
  mat->SetExactSol(gexact.ExactSolution(),4);
  cmesh->InsertMaterialObject(mat);

  // Add boundary conditions
  TPZManVector<REAL,1> val2(1,3.); // Part that goes to the RHS vector
  TPZFMatrix<REAL> val1(1,1,0.); // Part that goes to the Stiffnes matrix

  TPZBndCondT<REAL> *bcond = mat->CreateBC(mat, EBottom, 0, val1, val2);
  bcond->SetForcingFunctionBC(gexact.ExactSolution(),4);
  cmesh->InsertMaterialObject(bcond);

  bcond = mat->CreateBC(mat, ERight, 0, val1, val2);
  bcond->SetForcingFunctionBC(gexact.ExactSolution(),4);
  cmesh->InsertMaterialObject(bcond);
  
  bcond = mat->CreateBC(mat, EFront, 0, val1, val2);
  bcond->SetForcingFunctionBC(gexact.ExactSolution(),4);
  cmesh->InsertMaterialObject(bcond);

  bcond = mat->CreateBC(mat, EBack, 0, val1, val2);
  bcond->SetForcingFunctionBC(gexact.ExactSolution(),4);
  cmesh->InsertMaterialObject(bcond);

  bcond = mat->CreateBC(mat, ETop, 0, val1, val2);
  bcond->SetForcingFunctionBC(gexact.ExactSolution(),4);
  cmesh->InsertMaterialObject(bcond);
  
  bcond = mat->CreateBC(mat, ELeft, 0, val1, val2);
  bcond->SetForcingFunctionBC(gexact.ExactSolution(),4);
  cmesh->InsertMaterialObject(bcond);

  // Set up the computational mesh
  cmesh->AutoBuild();

  return cmesh;
}

REAL ElementDiameter(TPZGeoEl* gel) {
    REAL maxdist = 0.;
    int nnodes = gel->NNodes();
    for (int i = 0; i < nnodes; ++i) {
        TPZManVector<REAL,3> xi(3,0.), xj(3,0.);
        gel->Node(i).GetCoordinates(xi);
        for (int j = i+1; j < nnodes; ++j) {
            gel->Node(j).GetCoordinates(xj);
            REAL dist = 0.;
            for (int d = 0; d < gel->Dimension(); ++d) {
                dist += (xi[d] - xj[d]) * (xi[d] - xj[d]);
            }
            dist = sqrt(dist);
            if (dist > maxdist) maxdist = dist;
        }
    }
    return maxdist;
}

REAL MeshDiameter(TPZGeoMesh *gmesh) {
  REAL h = 0.;
  int64_t nel = gmesh->NElements();
  for (int64_t el = 0; el < nel; ++el) {
    TPZGeoEl *gel = gmesh->Element(el);
    if (gel && gel->Dimension() == gmesh->Dimension() && !gel->HasSubElement()) {
      REAL elSize = ElementDiameter(gel);
      if (elSize > h) h = elSize;
    }
  }
  return h;
}

void ConvOrder(TPZFMatrix<REAL> &errorMat, TPZVec<REAL> &hvector, int nref) {

  const char *errnames_h1[3] = {"H1", "L2", "Energy"};
  const char **errnames = nullptr;
  int nerrors = 3;
  errnames = errnames_h1;

  REAL order;

  std::cout << std::endl;
  for (int j = 0; j < nerrors; j++) {
    std::cout << errnames[j] << " error and order:" << std::endl;
    std::cout << "ref level 0 & " << std::scientific << std::setprecision(3) << errorMat(0, j)
              << " & - " << std::endl; // No order for the first refinement
    for (int i = 1; i < nref; i++) {
      order = log(errorMat(i, j) / errorMat(i - 1, j)) / log(hvector[i] / hvector[i - 1]);
      std::cout << "ref level " << i
                << " & " << std::scientific << std::setprecision(3) << errorMat(i, j)
                << " & " << std::fixed << std::setprecision(2) << order
                << std::endl;
    }
    std::cout << std::endl;
  }
}