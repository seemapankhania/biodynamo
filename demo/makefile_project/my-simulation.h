// -----------------------------------------------------------------------------
//
// Copyright (C) The BioDynaMo Project.
// All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
//
// See the LICENSE file distributed with this work for details.
// See the NOTICE file distributed with this work for additional information
// regarding copyright ownership.
//
// -----------------------------------------------------------------------------

#ifndef DEMO_MAKEFILE_PROJECT_MY_SIMULATION_H_
#define DEMO_MAKEFILE_PROJECT_MY_SIMULATION_H_

#include "biodynamo.h"

namespace bdm {

// Define compile time parameter
template <typename Backend>
struct CompileTimeParam : public DefaultCompileTimeParam<Backend> {};

inline int Simulate(int argc, const char** argv) {
  InitializeBiodynamo(argc, argv);

  // Define initial model - in this example: single cell at origin
  Cell cell({0, 0, 0});
  cell.SetDiameter(30);
  ResourceManager<>::Get()->push_back(cell);

  // Run simulation for one timestep
  Scheduler<> scheduler;
  scheduler.Simulate(1);

  std::cout << "Simulation completed successfully!" << std::endl;
  return 0;
}

}  // namespace bdm

#endif  // DEMO_MAKEFILE_PROJECT_MY_SIMULATION_H_
