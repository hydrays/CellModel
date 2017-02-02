#include <iostream>
#include <time.h>
#include "geometry.h"
#include "simulator.h"
#include "force.h"
#include "engine.h"

int Simulator::run()
{
  double t1, t2;

  Geometry geometry;
  geometry.init();

  Force force;
  force.update_force_field(geometry);

  Engine engine;
  engine.one_step_solution(geometry, force);

  return 0;
}






