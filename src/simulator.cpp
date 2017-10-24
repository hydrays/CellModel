#include <iostream>
#include <time.h>
#include "geometry.h"
#include "simulator.h"
#include "force.h"
#include "engine.h"
#include "parameters.h"

int Simulator::run(int process_id)
{
  Parameters params;
  params.init();
  
  Geometry geometry;
  geometry.init(params);

  Force force;
  Engine engine;

  engine.run(geometry, force, params, process_id);
  
  return 0;
}






