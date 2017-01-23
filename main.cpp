#include <stdio.h>
#include <stdlib.h>
#include "working_data.h"
#include "simulator.h"
  
int main()
{
  WorkingData working_data;
  working_data.init();

  Simulator simulate;
  simulate.run(working_data);

  return 0;
}
