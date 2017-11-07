#include <stdio.h>
#include <stdlib.h>
#include "simulator.h"
  
int main()
{
    std::ofstream flog("log.txt");
    if(!flog) 
    {
	std::cout << "file open error.\n";
	return -1; 
    }

    for ( int process_id=0; process_id < 50; process_id++ )
    {
	int status = 1;
	try
	{
	    Simulator simulate;
	    status = simulate.run(process_id);
	}
	catch (int e)
	{
	    std::cout << "An exception occurred. Exception Nr. " << e << "!\n";
	    flog << "Process " << process_id << "  exception: " << e << "!\n";
	}
	if ( status == 0 )
	{
	    flog << "Process " << process_id << "  terminated normally. \n";
	}
    }
    flog.close();	
    return 0;
}
