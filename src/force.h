#ifndef FORCE
#define FORCE

#include <algorithm>
#include <string>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <unordered_set>
#include <iomanip>
#include "ellipse.h"
#include "geometry.h"

class Force
{
public:
    std::vector<double> force_list;
    std::vector<double> pressure_list;
    double para_c = 2.0;

public:    
    int update_force_field(Geometry &geometry)
	{
	    compute_pressure(geometry);
	    compute_force_per_edge(geometry);
	}

    int compute_pressure(Geometry &geometry)
	{
	    pressure_list.resize(geometry.voronoi_cell_list.size());
	    for ( int k=0; k<geometry.voronoi_cell_list.size(); k++ )
	    {
		if ( geometry.voronoi_cell_list[k].area > 0.0 )
		{
		    pressure_list[k] = -(geometry.voronoi_cell_list[k].area - 4.0);
		}
		else
		{
		    std::cout << "compute_pressure: not suppose to happen.\n";
		    getchar();
		    pressure_list[k] = -1.0;
		}
	    }
	    return 0;
	}

    int compute_force_per_edge(Geometry &geometry)
    	{
    	    force_list.resize(geometry.voronoi_edge_list.size());
    	    for ( int k=0; k<geometry.voronoi_edge_list.size(); k++ )
    	    {
    		VoronoiEdge edge = geometry.voronoi_edge_list[k];
    		if ( edge.side_colors.size()==2 )
    		{
    		    double force = 0.0;
    		    // contact force
    		    int cell_i = geometry.map_cell_id[edge.side_colors[0]];
    		    int cell_j = geometry.map_cell_id[edge.side_colors[1]];
    		    force = force - 0.5*edge.length*(pressure_list[cell_i] + 
    			pressure_list[cell_j]);
    		    force_list[k] = force;
    		}
    		else
    		{
    		    force_list[k] = -1.0;
    		}
    	    }

	    // config force
    	    for ( int k=0; k<geometry.voronoi_cell_list.size(); k++ )
    	    {
    		VoronoiCell cell = geometry.voronoi_cell_list[k];
    		for ( int kk=0; kk<cell.edges.size(); kk++)
    		{
    		    int edge_id = cell.edges[kk].edge_id;
    		    force_list[edge_id] = force_list[edge_id] 
			+ para_c*cell.edges[kk].tri_sector_area/cell.area 
			- para_c*cell.edges[kk].ell_sector_area/cell.ellipse.area;
		}
    	    }
    	    return 0;
    	}

    int output_force_filed(Geometry &geometry, std::string file_index)
	{
	    std::ofstream foutput("out/force" + file_index + ".txt");
	    if(!foutput) 
	    {
		std::cout << "file open error.\n";
		return -1; 
	    }
	    std::ofstream foutput2("out/sum_force" + file_index + ".txt");
	    if(!foutput2) 
	    {
		std::cout << "file open error.\n";
		return -1; 
	    }
	    for ( int k=0; k<geometry.voronoi_cell_list.size(); k++ )
	    {
		if ( geometry.voronoi_cell_list[k].cell_id < 8000)
		{
		    double temp = sqrt(6.0*PI/sqrt(3));
		    double radius = sqrt(geometry.voronoi_cell_list[k].area/PI);
		    double a = geometry.voronoi_cell_list[k].ellipse.a;
		    double b = geometry.voronoi_cell_list[k].ellipse.b;
		    double sum_force_x = 0.0;
		    double sum_force_y = 0.0;
		    double oxi = geometry.voronoi_cell_list[k].ellipse.c1;
		    double oyi = geometry.voronoi_cell_list[k].ellipse.c2;
		    for ( int kk=0; kk<geometry.voronoi_cell_list[k].neighbor_id.size(); kk++ )
		    {
			int edge_id = geometry.voronoi_cell_list[k].edge_set[kk];
			int jcell_id = geometry.map_cell_id[geometry.voronoi_cell_list[k].neighbor_id[kk]];
			double cx = geometry.voronoi_edge_list[edge_id].c1;
			double cy = geometry.voronoi_edge_list[edge_id].c2;
			double oxj = geometry.voronoi_cell_list[jcell_id].ellipse.c1;
			double oyj = geometry.voronoi_cell_list[jcell_id].ellipse.c2;
			double theta = geometry.voronoi_cell_list[k].edges[kk].theta;
			foutput << geometry.voronoi_cell_list[k].cell_id << ", " << oxi << ", " << oyi << ", " 
				<< oxi + force_list[edge_id]*geometry.voronoi_cell_list[k].edges[kk].n1
				<< ", "
				<< oyi + force_list[edge_id]*geometry.voronoi_cell_list[k].edges[kk].n2
				<< ", " << force_list[edge_id] << "\n";
			sum_force_x = sum_force_x + force_list[edge_id]*geometry.voronoi_cell_list[k].edges[kk].n1;
			sum_force_y = sum_force_y + force_list[edge_id]*geometry.voronoi_cell_list[k].edges[kk].n2;
		    }
		    foutput2 << geometry.voronoi_cell_list[k].cell_id << ", " << oxi << ", " << oyi << ", " 
			     << oxi + sum_force_x << ", " 
			     << oyi + sum_force_y << ", "
			     << sum_force_x << ", "
			     << sum_force_y << "\n";
		    //foutput << "\n";
		}
	    }
	    foutput.close();
	    foutput2.close();
	    return 0;
	}
};



#endif //FORCE
