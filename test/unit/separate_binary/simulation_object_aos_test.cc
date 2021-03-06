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

#include "cell.h"
#include "compile_time_param.h"
#include "unit/simulation_object_test.h"

namespace bdm {

template <typename TBackend>
struct CompileTimeParam : public DefaultCompileTimeParam<TBackend> {
  using SimulationBackend = Scalar;
};

namespace simulation_object_aos_test_internal {

TEST(SimulationObjectTest, AosGetElementIndex) {
  Rm()->Clear();
  for (uint64_t i = 0; i < 10; i++) {
    Rm()->New<Cell>(1);
  }
  Rm()->Get<Cell>()->Commit();
  EXPECT_EQ(10u, Rm()->GetNumSimObjects());
  auto cells = Rm()->Get<Cell>();
  for (uint64_t i = 0; i < 10; i++) {
    EXPECT_EQ(i, (*cells)[i].GetElementIdx());
  }
}

}  // namespace simulation_object_aos_test_internal
}  // namespace bdm

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
