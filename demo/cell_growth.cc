#include "cell_growth.h"

#include <omp.h>
#include <cmath>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>

#include "backend.h"
#include "cell.h"
#include "compile_time_param.h"
#include "displacement_op.h"
#include "dividing_cell_op.h"
#include "exporter.h"
#include "resource_manager.h"
#include "scheduler.h"
#include "timing.h"
#include "timing_aggregator.h"

using bdm::Cell;
using bdm::SoaCell;
using bdm::Scalar;
using bdm::Soa;
using bdm::Timing;
using bdm::TimingAggregator;
using bdm::ExporterFactory;
using bdm::Exporter;
using bdm::Grid;

void Execute(size_t cells_per_dim, size_t iterations, size_t threads,
             size_t repititions, TimingAggregator *statistic,
             bool with_export) {
  const double space = 20;

  for (size_t r = 0; r < repititions; r++) {
    std::stringstream ss;
    ss << "measurement " << r << " - " << threads << " thread(s) - "
       << cells_per_dim << " cells per dim - " << iterations << " iteration(s)";
    statistic->AddDescription(ss.str());

    // bdm::TransactionalVector<Cell> cells;
    auto cells = Cell::NewEmptySoa();
    cells.reserve(cells_per_dim * cells_per_dim * cells_per_dim);
    {
      Timing timing("Setup", statistic);
      for (size_t i = 0; i < cells_per_dim; i++) {
        for (size_t j = 0; j < cells_per_dim; j++) {
          for (size_t k = 0; k < cells_per_dim; k++) {
            Cell cell({i * space, j * space, k * space});
            cell.SetDiameter(30);
            cell.SetAdherence(0.4);
            cell.SetMass(1.0);
            cell.UpdateVolume();
            cells.push_back(cell);
          }
        }
      }
    }

    // iterate for all (time) steps
    ExporterFactory exp_fac;

    auto exp_basic = exp_fac.GenerateExporter<SoaCell>("basic");
    auto exp_matlab = exp_fac.GenerateExporter<SoaCell>("matlab");
    auto exp_neuroml = exp_fac.GenerateExporter<SoaCell>("neuroml");
    auto exp_paraview = exp_fac.GenerateExporter<SoaCell>("paraview");

    if (with_export) {
      exp_paraview->CreatePVDFile("Results4Paraview", iterations, 1.0);
    }

    for (size_t i = 0; i < iterations; i++) {
      {
        Timing timing("Find Neighbors", statistic);
        Grid<>::GetInstance().UpdateGrid();
      }

      {
        Timing timing("Cell Growth", statistic);
        bdm::DividingCellOp biology;
        biology(&cells, 0);
      }

      {
        Timing timing("Displacement", statistic);
        bdm::DisplacementOp<> op;
        op(&cells, 0);
      }

      if (with_export) {
        Timing timing("Export", statistic);
        std::cout << "exporting now..." << std::endl;
        exp_basic->ToFile(cells, "FinalPositions.dat");
        exp_matlab->ToFile(cells, "FinalPositions.m");
        exp_neuroml->ToFile(cells, "FinalPositions.xml");
        exp_paraview->ToFile(cells, "Results4Paraview");
        exp_paraview->AddIteration();
      }
    }
  }
}

void Scaling(size_t cells_per_dim, size_t iterations, size_t repititions,
             TimingAggregator *statistic, bool with_export,
             const std::function<void(int &)> thread_inc =
                 [](int &i) {  // NOLINT(runtime/references)
                   i *= 2;
                 },
             const int max_threads = omp_get_max_threads()) {
  for (int i = 1; i <= max_threads; thread_inc(i)) {
    omp_set_num_threads(i);
    Execute(cells_per_dim, iterations, i, repititions, statistic, with_export);
  }
}

int main(int args, char **argv) {
  TimingAggregator statistic;
  size_t repititions = 1;
  bool do_export = true;
  if (args == 2 &&
      (std::string(argv[1]) == "help" || std::string(argv[1]) == "--help")) {
    // clang-format off
    std::cout << "SYNOPSIS\n"
              << "  ./cell_growth help | --help |\n"
              << "         [#repititions] | \n"
              << "         [#cells_per_dim #iterations #threads [#repititions] | \n"
              << "         --scaling [#repititions] | \n"
              << "         --detailed-scaling [#repititions] \n"
              << "\nDESCRIPTION\n"
              << "  Creates a three dimensional grid of cells, calculates neighbors, simulates \n"
              << "  cell growth and calculates displacement based on mechanical forces\n"
              << "  outputs runtime statistic for each operation\n"
              << "\nOPTIONS\n"
              << "  help | --help\n"
              << "    Explains usage of this binary and its command line options\n"
              << "\n  [#repititions]\n"
              << "     number of cells per dimension: 8 (-> total number of cells: 8^3 = 512)\n"
              << "     number of iterations:          1\n"
              << "     number of threads:             1\n"
              << "     number of repititions:         according to parameter - 1 if not specified\n"
              << "\n  --scaling [#repititions]\n"
              << "     executes the simulation several times with different number of threads\n"
              << "     number of cells per dimension: 128 (-> total number of cells: 128^3 =~ 2.1M)\n"
              << "     number of iterations:          1\n"
              << "     number of threads:             1 - logical CPUs on the system - incremented *= 2\n"
              << "     number of repititions:         according to parameter - 1 if not specified\n"
              << "\n  --detailed-scaling [#repititions]\n"
              << "     executes the simulation several times with different number of threads\n"
              << "     number of cells per dimension: 128 (-> total number of cells: 128^3 =~ 2.1M)\n"
              << "     number of iterations:          1\n"
              << "     number of threads:             1 - logical CPUs on the system - threads incremented += 1\n"
              << "     number of repititions:         according to parameter - 1 if not specified\n"
              << std::endl;
    // clang-format on
    return 0;
  } else if (args >= 4) {
    size_t cells;
    size_t iterations;
    size_t threads;
    std::istringstream(std::string(argv[1])) >> cells;
    std::istringstream(std::string(argv[2])) >> iterations;
    std::istringstream(std::string(argv[3])) >> threads;
    if (args == 5) {
      std::istringstream(std::string(argv[4])) >> repititions;
    }
    omp_set_num_threads(threads);
    Execute(cells, iterations, threads, repititions, &statistic, do_export);
  } else if (args >= 2 && std::string(argv[1]) == "--scaling") {
    if (args == 3) {
      std::istringstream(std::string(argv[2])) >> repititions;
    }
    Scaling(256, 1, repititions, &statistic, do_export);
  } else if (args >= 2 && std::string(argv[1]) == "--detailed-scaling") {
    if (args == 3) {
      std::istringstream(std::string(argv[2])) >> repititions;
    }
    Scaling(256, 1, repititions, &statistic, do_export, [](int &i) { i++; });
  } else {
    omp_set_num_threads(1);
    if (args == 2) {
      std::istringstream(std::string(argv[1])) >> repititions;
    }
    Execute(8, 1, 1, repititions, &statistic, do_export);
  }
  std::cout << statistic << std::endl;
  return 0;
}
