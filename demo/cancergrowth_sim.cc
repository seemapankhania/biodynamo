#include "biology_module_util.h"
#include "cell.h"
#include "resource_manager.h"
#include "scheduler.h"
#include "simulation_object_util.h"
#include "bdm_interface.h"

namespace bdm {

// 1. Define growth behaviour
  struct GrowthModule {

    int curr_celltype;
    double rand_double;

  template <typename T>
    void Run(T* cell) {

      if (cell->GetCellTypeID()==0) {
       rand_double = ((double) rand() / (RAND_MAX));
       if (rand_double<0.7) {
        curr_celltype = 1;
      }
      else {
        curr_celltype = 2;
      }
    }

    if (cell->GetDiameter() <= 40) {
      cell->ChangeVolume(300000);
    } else {
      Divide(
        *cell,
        ResourceManager<Cell<Soa, variant<GrowthModule>>>::Get()->GetCells());
        //cell->SetCellProperty("celltype",curr_celltype);
    }
  }

  bool IsCopied(Event event) const { return true; }
};

// 2. Define biology modules that should be used in this simulation
typedef variant<GrowthModule> BiologyModules;

// 3. Use predefined cell class as is and pass biology module definitions
//    Hence the biology module template parameter can be ommitted later on
template <typename Backend = Scalar>
using MyCell = Cell<Backend, BiologyModules>;

void Simulate(size_t cells_per_dim = 16) {
  // 4. Get cell container
  auto cells = ResourceManager<MyCell<Soa>>::Get()->GetCells();
  cells->reserve(cells_per_dim * cells_per_dim * cells_per_dim);

  int curr_celltype, cellcounter = 1;
  double rand_double;

  int type1_counter=0;
  int type2_counter=0;
  double mass1_sum = 0.0;
  double mass2_sum = 0.0;

  // 5. Define initial model - in this case 3D grid of cells
  double space = 20;
  for (size_t i = 0; i < cells_per_dim; i++) {
    for (size_t j = 0; j < cells_per_dim; j++) {
      for (size_t k = 0; k < cells_per_dim; k++) {
        MyCell<Scalar> cell({i * space, j * space, k * space});
        cell.SetDiameter(30);
        cell.SetAdherence(0.4);

        cell.UpdateVolume();
        cell.AddBiologyModule(GrowthModule());

        rand_double = ((double) rand() / (RAND_MAX));
        if (rand_double<0.7) {

          curr_celltype = 1;
          cell.SetMass(1.2);
            //std::cout << "should imprint type 1 for normoxic cells" << std::endl;
        }
        else {
          curr_celltype = 2;
          cell.SetMass(0.8);
            //std::cout << "should imprint type hypoxic cells" << std::endl;
        }        
        cell.SetCellTypeID(curr_celltype);
        cells->push_back(cell);

        cellcounter++;
      }
    }
  }

  // 6. Run simulation for one timestep
  Scheduler scheduler;
  //double a_rand;


  for (int curr_time = 0; curr_time < 20; curr_time++) {

    type1_counter=0;
    type2_counter=0;
    mass1_sum = 0.0;
    double mass2_sum = 0.0;
    
    scheduler.Simulate<MyCell<Soa>>(1);

  // Readout relevant information

    for (size_t i = 0; i < cells->size(); i++) {
      auto &&cell = (*cells)[i];

      curr_celltype = cell.GetCellTypeID();
      if (curr_celltype<2) {
        type1_counter++;
        mass1_sum+=cell.GetMass();
      }
      else {
        type2_counter++;
        mass2_sum+=cell.GetMass();
      }

      //a_rand = ((double) rand() / (RAND_MAX));
      //if (a_rand<0.01) {
      //  std::cout << "a random cell diameter and type: " << cell.GetDiameter() << " / " << curr_celltype << std::endl;
      //}


    }

    std::cout << "Current time step: " << curr_time << ". Number of type 0 and 1 cells: " << type1_counter << " / " << type2_counter << std::endl;
    std::cout << "Mass of type 0 and 1 cells: " << mass1_sum << " / " << mass2_sum << std::endl;


  }

  DiscontinuousInterfaceData myDisc_fd(1);
  myDisc_fd.normoxic_cells_mass.assign(1,mass1_sum);
  myDisc_fd.hypoxic_cells_mass.assign(1,mass2_sum);
  myDisc_fd.normoxic_cells_population.assign(1,type1_counter);
  myDisc_fd.hypoxic_cells_population.assign(1,type2_counter);

}

}  // namespace bdm

int main() { bdm::Simulate(); }
