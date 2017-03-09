#include <iostream>
#include <time.h>
#include "geometry.h"
#include "simulator.h"
#include "force.h"
#include "engine.h"

int Simulator::run(int process_id)
{
  // std::cout << "test atan2: " << atan2(0.0, 0.1) << "\n";
  // getchar();

  double t1, t2;

  Geometry geometry;
  geometry.init();

  Force force;

  Engine engine;
  //engine.one_step_solution(geometry, force);
  engine.run(geometry, force, process_id);

  
  return 0;
}






