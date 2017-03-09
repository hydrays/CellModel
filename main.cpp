#include <stdio.h>
#include <stdlib.h>
#include "simulator.h"
  
int main()
{
    for ( int process_id=0; process_id < 100; process_id++ )
    {
	try
	{
	    Simulator simulate;
	    simulate.run(process_id);
	}
	catch (int e)
	{
	    std::cout << "An exception occurred. Exception Nr. " << e << '\n';
	}
    }
    return 0;
}
