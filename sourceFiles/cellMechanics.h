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
#ifndef HEMO_CELLMECHANICS
#define HEMO_CELLMECHANICS
namespace hemo {
  class CellMechanics;
}

#include "hemoCellParticleField.h"
#include "hemoCellParticle.h"
#include "commonCellConstants.h"
#include "meshMetrics.h"
#include "constantConversion.h"

namespace hemo {
class CellMechanics {
  public:
  const CommonCellConstants cellConstants;
  Config & cfg;

  CellMechanics(HemoCellField & cellfield, Config & modelCfg_) : cellConstants(CommonCellConstants::CommonCellConstantsConstructor(cellfield, modelCfg_)), cfg(modelCfg_) {}
  virtual ~CellMechanics() {};
  
  virtual void ParticleMechanics(std::map<int,std::vector<HemoCellParticle *>> &,const std::map<int,bool> &, pluint ctype) = 0 ;
  virtual void statistics() = 0;
  virtual void solidifyMechanics(const std::map<int,std::vector<int>>&,std::vector<HemoCellParticle>&,plb::BlockLattice3D<T,DESCRIPTOR> *,plb::BlockLattice3D<T,CEPAC_DESCRIPTOR> *, pluint ctype, HemoCellParticleField &) {};
  


  /*T calculate_kES(Config & cfg, plb::MeshMetrics<T> & meshmetric ){
    T eqLength = 5e-7/param::dx;
    return cfg["MaterialModel"]["kES"].read<T>() * param::kBT_lbm / eqLength;
  };

  T calculate_kAD(Config & cfg, plb::MeshMetrics<T> & meshmetric ){
  T eqLength = 5e-7/param::dx;
  return cfg["MaterialModel"]["kAD"].read<T>() * param::kBT_lbm / eqLength;
  };*/

  T calculate_kScale(Config & cfg){
    T kScale = cfg["MaterialModel"]["kScale"].read<T>();
    return  kScale;
  };

  T calculate_tauA(Config & cfg){
    T tauA = cfg["MaterialModel"]["tauA"].read<T>();
    return  tauA;
  };

  T calculate_tauV(Config & cfg){
  T tauV = cfg["MaterialModel"]["tauV"].read<T>();
    return  tauV;
  };

  T calculate_tauB(Config & cfg){
  T tauB = cfg["MaterialModel"]["tauB"].read<T>();
    return  tauB;
  };

  T calculate_tauL(Config & cfg){
  T tauL = cfg["MaterialModel"]["tauL"].read<T>();
    return  tauL;
  };

  /*T calculate_tauS(Config & cfg){
  T tauS = cfg["MaterialModel"]["tauS"].read<T>();
    return  tauS;
  };*/

  T calculate_kLink(Config & cfg, plb::MeshMetrics<T> & meshmetric){
    T kLink = cfg["MaterialModel"]["kLink"].read<T>();
    //T persistenceLength = cfg["MaterialModel"]["persistenceLength"].read<T>();
    //T plc = persistenceLength/param::dx; //* sqrt((meshmetric.getNumVertices()-2.0) / (23867-2.0)); // <- Constant from the Karniadakis group. We don't use this.

    T persistenceLengthFine = 7.5e-9; // In [m] -> this is the biological value.
    T plc = persistenceLengthFine/param::dx; //* sqrt((meshmetric.getNumVertices()-2.0) / (23867-2.0)); // <- Constant from the Karniadakis group. We don't use this.
    return  kLink * param::kBT_lbm/plc;
  };
  
  T calculate_kBend(Config & cfg, plb::MeshMetrics<T> & meshmetric ){
    T eqLength = 5e-7/param::dx;
    return cfg["MaterialModel"]["kBend"].read<T>() * param::kBT_lbm / eqLength;
  };

  T calculate_kVolume(Config & cfg, plb::MeshMetrics<T> & meshmetric){
    T kVolume =  cfg["MaterialModel"]["kVolume"].read<T>();
    T eqLength = 5e-7/param::dx;
    T NfacesScaling = 1280.0/cellConstants.triangle_list.size();
    return kVolume * NfacesScaling * param::kBT_lbm / eqLength;
  };

  T calculate_kArea(Config & cfg, plb::MeshMetrics<T> & meshmetric){
    T kArea =  cfg["MaterialModel"]["kArea"].read<T>();
    T eqLength = 5e-7/param::dx;
    T NfacesScaling = 1280.0/cellConstants.triangle_list.size();
    return kArea * NfacesScaling * param::kBT_lbm/(eqLength);
  };

  /*T calculate_kSurf(Config & cfg, plb::MeshMetrics<T> & meshmetric){
    T kSurf =  cfg["MaterialModel"]["kSurf"].read<T>();
    T eqLength = 5e-7/param::dx;
    T NfacesScaling = 1280.0/cellConstants.triangle_list.size();
    return kSurf * NfacesScaling * param::kBT_lbm/(eqLength);
  };*/

  T calculate_etaM(Config & cfg ){
    return cfg["MaterialModel"]["eta_m"].read<T>() * param::dx / param::dt / param::df;
  };
};
}
#endif
