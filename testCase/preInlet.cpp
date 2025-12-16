/*
This file is part of the HemoCell library

HemoCell is developed and maintained by the Computational Science Lab
in the University of Amsterdam. Any questions or remarks regarding this library
can be sent to: info@hemocell.eu

When using the HemoCell library in scientific work please cite the
corresponding paper: https://doi.org/10.3389/fphys.2017.00563

The HemoCell library is free software: you can redistribute it and/or
modify it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

The library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "hemocell.h"
#include <helper/voxelizeDomain.h>
#include "rbcHighOrderModel.h"
#include "pltSimpleModel.h"
#include "cellInfo.h"
#include "fluidInfo.h"
#include "particleInfo.h"
#include "writeCellInfoCSV.h"
#include <fenv.h>
#include "preInlet.h"

#include "palabos3D.h"
#include "palabos3D.hh"

using namespace hemo;

int main(int argc, char *argv[]) {
  if(argc < 2) {
    cout << "Usage: " << argv[0] << " <configuration.xml>" << endl;
    return -1;
  }

  HemoCell hemocell(argv[1], argc, argv);
  Config * cfg = hemocell.cfg;

  hlog << "(Stl preinlet) (Geometry) reading and voxelizing STL file " << (*cfg)["domain"]["geometry"].read<string>() << endl;
  MultiScalarField3D<int> *flagMatrix = 0;
  VoxelizedDomain3D<double> * voxelizedDomain = 0;
  getFlagMatrixFromSTL((*cfg)["domain"]["geometry"].read<string>(),
                       (*cfg)["domain"]["fluidEnvelope"].read<int>(),
                       (*cfg)["domain"]["refDirN"].read<int>(),
                       (*cfg)["domain"]["refDir"].read<int>(),
                       voxelizedDomain, flagMatrix,
                       (*cfg)["domain"]["blockSize"].read<int>());

  hlog << "(Stl preinlet) (Parameters) setting lbm parameters" << endl;
  param::lbm_base_parameters((*cfg));
  param::printParameters();

  hlog << "(PreInlets) creating preInlet" << endl;
  hemocell.preInlet = new hemo::PreInlet(&hemocell,flagMatrix);

  Box3D slice = flagMatrix->getBoundingBox();
  // You can put the inlet inside the domain, or any other side if needed
  //slice.x0 = slice.x1 - 10;

  // Now we just put it at the beginning of the domain
  slice.x0 = slice.x1 -10 ;

  // Direction:: -> define the inflow direction (preInlet is on the X negative side)
  hemocell.preInlet->preInletFromSlice(Direction::Xpos,slice);
  //hemocell.preInlet->preInletFromSlice(Direction::Xneg,slice) //Positive side

  hlog << "(Stl preinlet) (Fluid) Initializing Palabos Fluid Field" << endl;
  hemocell.initializeLattice(voxelizedDomain->getMultiBlockManagement());

  if (!hemocell.partOfpreInlet) {
    hemocell.lattice->periodicity().toggleAll(false);
  }

  // Setting Preinlet creation
  hemocell.preInlet->initializePreInlet();


  hlog << "(Stl preinlet) (Fluid) Setting up boundaries in Palabos Fluid Field" << endl;
  boundaryFromFlagMatrix(hemocell.lattice,flagMatrix,hemocell.partOfpreInlet);


  hemocell.preInlet->createBoundary();

  hemocell.lattice->toggleInternalStatistics(false);

  hemocell.latticeEquilibrium(1.,plb::Array<double, 3>(0.,0.,0.));

  // Driving Force
  hemocell.preInlet->calculateDrivingForce();

  hemocell.lattice->initialize();

  // Adding all the cells
  hemocell.initializeCellfield();
  string RDW =  (*cfg)["UDV"]["RDW"].read<string>();
  string PLT =  (*cfg)["UDV"]["PLT"].read<string>();

  hemocell.addCellType<RbcHighOrderModel>("RBC", RBC_FROM_SPHERE);
  hemocell.setMaterialTimeScaleSeparation("RBC", (*cfg)["ibm"]["stepMaterialEvery"].read<int>());
  hemocell.setInitialMinimumDistanceFromSolid("RBC", 0.2); //Micrometer! not LU

  if (RDW == "True")
  {
  hemocell.addCellType<RbcHighOrderModel>("LRBC", RBC_FROM_SPHERE);
  hemocell.setMaterialTimeScaleSeparation("LRBC", (*cfg)["ibm"]["stepMaterialEvery"].read<int>());
  hemocell.setInitialMinimumDistanceFromSolid("LRBC", 0.2); //Micrometer! not LU*/

  hemocell.addCellType<RbcHighOrderModel>("SRBC", RBC_FROM_SPHERE);
  hemocell.setMaterialTimeScaleSeparation("SRBC", (*cfg)["ibm"]["stepMaterialEvery"].read<int>());
  hemocell.setInitialMinimumDistanceFromSolid("SRBC", 0.2); //Micrometer! not LU
  }
  if (PLT == "True")
  {
  hemocell.addCellType<PltSimpleModel>("PLT", ELLIPSOID_FROM_SPHERE);
  hemocell.setMaterialTimeScaleSeparation("PLT", (*cfg)["ibm"]["stepMaterialEvery"].read<int>());
  }

  hemocell.ActivateRepulsion();

  vector<int> outputs = {OUTPUT_POSITION, OUTPUT_VELOCITY, OUTPUT_TRIANGLES, OUTPUT_FORCE,  OUTPUT_FORCE_VOLUME,OUTPUT_FORCE_BENDING,OUTPUT_FORCE_LINK,OUTPUT_FORCE_AREA, OUTPUT_FORCE_REPULSION};
  hemocell.setOutputs("RBC", outputs);
  if (RDW == "True")
  {
  hemocell.setOutputs("LRBC", outputs);
  hemocell.setOutputs("SRBC", outputs);
  }
  if (PLT == "True")
  {
  hemocell.setOutputs("PLT", outputs);
  }
  outputs = {OUTPUT_VELOCITY, OUTPUT_DENSITY, OUTPUT_FORCE, OUTPUT_BOUNDARY, OUTPUT_SHEAR_RATE, OUTPUT_STRAIN_RATE, OUTPUT_SHEAR_STRESS, OUTPUT_OMEGA, OUTPUT_CELL_DENSITY};
  hemocell.setFluidOutputs(outputs);

  double stepFunction = (*cfg)["UDV"]["stepFunction"].read<double>();
  unsigned int stepActive =  (*cfg)["UDV"]["stepActive"].read<unsigned int>();
  unsigned int stepDeActive =  (*cfg)["UDV"]["stepDeActive"].read<unsigned int>();
  
  string sinusoidal =  (*cfg)["UDV"]["sinusoidal"].read<string>();
  double frequency = (*cfg)["UDV"]["frequency"].read<double>(); //frequency with respect to tmax, not second!
  double amplitude = (*cfg)["UDV"]["amplitude"].read<double>(); //outlet pressure is 1





  // For the main simulation domain we have to define outlets
  if (!hemocell.partOfpreInlet) {
    Box3D bb = hemocell.lattice->getBoundingBox();
    Box3D outlet(bb.x0,bb.x0+2,bb.y0,bb.y1,bb.z0,bb.z1); //for negative
    //BOX3D outlet(bb.x1,bb.x1,bb.y0,bb.y1,bb.z0,bb.z1); //for positive 
    OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundary = new BoundaryConditionInstantiator3D
          < T, DESCRIPTOR, WrappedZouHeBoundaryManager3D<T,DESCRIPTOR> > ();
    boundary->addPressureBoundary0P(outlet,*hemocell.lattice,boundary::density);
    //boundary->addPressureBoundary0P(outlet,*hemocell.lattice,boundary::outflow);
    setBoundaryDensity(*hemocell.lattice,outlet, 1.0);
  }

  //loading the cellfield
  if (not cfg->checkpointed) {
    hemocell.loadParticles();
    hemocell.writeOutput();
  } else {
    hemocell.loadCheckPoint();
  }

  if (hemocell.iter == 0) {
    pcout << "Running fluid only simulation to get converged fluid field "  << (*cfg)["parameters"]["warmup"].read<plint>() << " iterations..." << endl;
    for (plint itrt = 0; itrt < (*cfg)["parameters"]["warmup"].read<plint>(); ++itrt) {
    //pcout << "\t Average energy.: " << computeAverageEnergy(*hemocell.lattice) << endl;
    hemocell.lattice->collideAndStream();
    }
  }

  unsigned int tmax = (*cfg)["sim"]["tmax"].read<unsigned int>();
  unsigned int tcsv = (*cfg)["sim"]["tcsv"].read<unsigned int>();
  unsigned int tmeas = (*cfg)["sim"]["tmeas"].read<unsigned int>();
  unsigned int tcheckpoint = (*cfg)["sim"]["tcheckpoint"].read<unsigned int>();

  //unsigned int tbalance = (*cfg)["sim"]["tbalance"].read<unsigned int>();
  if (RDW == "True")
  {
  pcout << "<<<<<<<<<<<<<<<<RDW is active>>>>>>>>>>>>>>>>>>>" << endl;
  }
  if (PLT == "True")
  {
  pcout << "<<<<<<<<<<<<<<<<PLT is active>>>>>>>>>>>>>>>>>>>" << endl;
  }
  if (stepFunction != 1)
  {
  pcout << "Step funciton will be applied between: " << stepActive << " and " << stepDeActive << endl;
  }

  if (sinusoidal == "True")
  {
  pcout << "Sinusoidal profile will be applied to outlet pressure with frequency of: " << frequency << " and the amplitude of: " << amplitude <<endl;
  }


  while (hemocell.iter < tmax ) {
  ParticleStatistics pinfo = ParticleInfo::calculateForceStatistics(&hemocell); double topN = param::df * 1.0e12;
  FluidStatistics finfo = FluidInfo::calculateVelocityStatistics(&hemocell); double toMpS = param::dx / param::dt;
    
    
    //T value = 0;
    Box3D bb = hemocell.lattice->getBoundingBox();
    Box3D outlet(bb.x0,bb.x0+2,bb.y0,bb.y1,bb.z0,bb.z1); //for negative
    if (stepFunction != 1)
    {

      if (hemocell.iter >= stepActive) 
      {
      setBoundaryDensity(*hemocell.lattice,outlet, stepFunction);
      }
      if (hemocell.iter >= stepDeActive)
      {
        stepFunction=1;
      }
    }

    if (sinusoidal == "True")
    {

      static T pi = std::acos((T)-1);
      double oscillation=1+amplitude*(std::sin(2*pi * (frequency*hemocell.iter) / (tmax)));
      setBoundaryDensity(*hemocell.lattice,outlet, oscillation);

    }

    //preinlet.update();
    hemocell.iterate();

    if (hemocell.partOfpreInlet) {
      //Set driving force as required after each iteration
      hemocell.preInlet->setDrivingForce();
    }
    hemocell.preInlet->applyPreInlet();

 if (hemocell.iter % tmeas == 0) {

      pcout << "(main) Stats. @ " <<  hemocell.iter << " (" << hemocell.iter * param::dt << " s):" << endl;
      pcout << "\t Average energy.: " << computeAverageEnergy(*hemocell.lattice) << endl;
      //pcout << "\t # of cells: " << CellInformationFunctionals::getTotalNumberOfCells(&hemocell);
      pcout << "\t num of RBC: " << CellInformationFunctionals::getNumberOfCellsFromType(&hemocell, "RBC") << endl;
      if (RDW == "True")
      {
      pcout << "\t num of SRBC: " << CellInformationFunctionals::getNumberOfCellsFromType(&hemocell, "SRBC") << endl;
      pcout << "\t num of LRBC: " << CellInformationFunctionals::getNumberOfCellsFromType(&hemocell, "LRBC") << endl;
      }
      if (PLT == "True")
      {
      pcout << "\t num of PLT: " << CellInformationFunctionals::getNumberOfCellsFromType(&hemocell, "PLT") << endl;
      }
      pcout << "\t Velocity(m/s) -  max.: " << finfo.max * toMpS << " mean: " << finfo.avg * toMpS <<  endl;
      pcout << "\t Velocity(LU)  -  max.: " << finfo.max << " mean: " << finfo.avg << endl;
      pcout << "\t Force(pN)     -  max.: " << pinfo.max * topN << " mean: " << pinfo.avg * topN << endl;

    hemocell.writeOutput();
    }
    if ( (pinfo.max * topN) > 100000) {
    pcout << "\t ! ! ! F O R C E E X C E E D ! ! ! " << pinfo.max *topN <<  " pN, at iter "  <<  hemocell.iter << endl;

    }

    if (hemocell.iter % tcsv == 0) {
      hlog << " /(>_>)/             " << hemocell.iter  <<   "             /(>_>)/ "  << endl;
      writeCellInfo_CSV(hemocell);
    }

    if (hemocell.iter % tcheckpoint == 0) {
      hemocell.saveCheckPoint();
    }
  }

  pcout << "End of the Simulation" << endl;

  return 0;
}
