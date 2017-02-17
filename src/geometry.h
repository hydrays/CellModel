
#ifndef GEOMETRY
#define GEOMETRY

#include <algorithm>
#include <string>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <unordered_set>
#include <iomanip>
#include <map>
#include "ellipse.h"

#define PI 3.1415926535

struct PixelNode
{
    int irow;
    int icol;
    int degree;
    int node_id;

    PixelNode(int i, int j, int d, int node_id)
	: irow(i), icol(j), degree(d), node_id(node_id)
	{
	}
};

struct BoundaryData
{
    double L;
    double H;
    double margin_L;
    double margin_H;
    double res_x, res_y;
    int n, m;
    int nmin, nmax, mmin, mmax;    
    double xmin, xmax, ymin, ymax;    

    int set_boundary_data(double L_value, double H_value)
	{
	    L = L_value;
	    H = H_value;
	    int margin_n = 80;
	    int margin_m = 80;
	    int inner_n = 400;
	    int inner_m = 400;
	    n = inner_n + 2*margin_n;
	    m = inner_m + 2*margin_m;
	    res_x = 2.0*L/inner_m;
	    res_y = 2.0*H/inner_n;
	    margin_L = margin_m*res_x;
	    margin_H = margin_n*res_y;
	    // nmin = margin_n;
	    // nmax = margin_n + inner_n;
	    // mmin = margin_m;
	    // mmax = margin_m + inner_m;
	    nmin = 1;
	    nmax = n-1;
	    mmin = 1;
	    mmax = m-1;
	    xmin = -L-margin_L;
	    xmax = L + margin_L;
	    ymin = -H-margin_H;
	    ymax = H + margin_H;

	    if ( nmin <= 0 || nmax >= n || mmin <= 0 || mmax >= m )
	    {
		std::cout << "error: not enough marginal region. \n";
		getchar();
	    }
	    std::cout << "Status: boundary data initialized. \n";
	}
};

struct VoronoiNode
{
    double x;
    double y;
    std::unordered_set<int> color_set;
};

struct VoronoiEdge
{
    int edge_id;
    int node_id1;
    int node_id2;
    double x1, y1;
    double x2, y2;
    double c1, c2;
    double n1;
    double n2;
    double length;
    double theta;
    double theta1, theta2;
    double tri_sector_area;
    double ell_sector_area;
    std::vector<int> side_colors;
    VoronoiEdge(int edge_id, int id1, int id2)
	: edge_id(edge_id), node_id1(id1), node_id2(id2)
	{
	    n1 = 0.0;
	    n2 = 0.0;
	    length = -1.0;
	}
};

struct VoronoiCell
{
    bool is_valid;
    int cell_id;
    double c1, c2;
    Ellipse ellipse;
    std::vector<int> node_set;
    std::vector<int> edge_set;
    std::vector<VoronoiNode> nodes;
    std::vector<VoronoiEdge> edges;
    std::vector<VoronoiCell> neighbors;
    double area;
    std::vector<int> neighbor_id;
};

struct Hash {
    std::hash<int> int_hash;
    size_t operator() (const VoronoiEdge &edge) const {
	int temp;
	temp = (53 + int_hash(edge.node_id1)) * 53 + int_hash(edge.node_id2);
	return (temp);
    }
};

struct ColorHash {
    std::hash<int> int_hash;
    size_t operator() (const std::unordered_set<int> &node_color) const {
	std::vector<int> values;
	for (std::unordered_set<int>::const_iterator it = node_color.begin();
	     it != node_color.end(); it++ )
	{
	    values.push_back(*it);
	}
	sort(values.begin(), values.end());
	int temp = 0;
	for (std::vector<int>::iterator it=values.begin(); it!=values.end(); ++it)
	{
	    temp = (temp) * 53 + int_hash(*it);
	}
	return (temp);
    }
};

inline bool operator == (VoronoiEdge const& edge1, VoronoiEdge const& edge2)
{
    return ( (edge1.node_id1 == edge2.node_id1) && 
	     (edge2.node_id2 == edge2.node_id2) );
}

class Geometry
{
public:
    int MAX_NUM_ELLIPSE = 100000;
    int cell_number;
    BoundaryData boundary_data;
    std::vector<VoronoiCell> voronoi_cell_list;
    std::vector<VoronoiNode> voronoi_node_list;
    std::vector<VoronoiEdge> voronoi_edge_list;
    std::unordered_set<VoronoiEdge, Hash> primary_edge_set;
    std::vector<Ellipse> ellipse_list;
    std::vector<Ellipse> new_ellipse_list;
    std::map<int,int> map_cell_id;

    int init()
	{
	    //boundary_data.set_boundary_data(18.0, 11.9);
	    boundary_data.set_boundary_data(20.0, 17.5);
	    //boundary_data.set_boundary_data(15.0, 10.0);
	    //boundary_data.set_boundary_data(36.0, 36.0);

	    init_ellipse_list();
	    cell_number = ellipse_list.size();
	    std::cout << "Status: ellipse list initialized: (num: " << ellipse_list.size() << ") \n";
	    apply_boundary_condition_ellipse();
	    std::cout << "Status: ellipse extended: (num: " << ellipse_list.size() << ") \n";
	    //getchar();
	    ellipse_to_voronoi();
	    update_map();

	    //voronoi_to_ellipse();
	    //output_new_ellipse_list();
	}

    int update_map()
	{
	    map_cell_id.clear();
	    for ( int k=0; k<voronoi_cell_list.size(); k++ )
	    {
		map_cell_id.insert( std::pair<int,int>(voronoi_cell_list[k].cell_id, k) );
	    }
	    return 0;
	}
    
    int apply_boundary_condition_ellipse()
	{
	    double L = boundary_data.L;
	    double H = boundary_data.H;
	    double margin_L = boundary_data.margin_L;
	    double margin_H = boundary_data.margin_H;

	    for ( int k=0; k<ellipse_list.size(); k++ )
	    {
		if ( (ellipse_list[k].c1 < -L) || (ellipse_list[k].c1 > L) 
		     || (ellipse_list[k].c2 < -H) || (ellipse_list[k].c2 > H) )
		{
		    std::cout << "remove ellipse outside the bounding box..." 
			      << ellipse_list[k].ellipse_id << " "
			      << ellipse_list[k].c1 << " " 
			      << ellipse_list[k].c2 << "\n";
		    ellipse_list.erase(ellipse_list.begin()+k);
		    getchar();
		}
	    }
	    for ( int k=0; k<ellipse_list.size(); k++ )
	    {
		// left
		if ( (ellipse_list[k].c1 <= L) && (ellipse_list[k].c1 >= L - margin_L)
		     && (ellipse_list[k].c2 <= H) && (ellipse_list[k].c2 >= -H) )
		{
		    Ellipse ellipse = ellipse_list[k];
		    ellipse.c1 = ellipse.c1 - 2.0*L;
		    ellipse.ellipse_id = ellipse.ellipse_id + MAX_NUM_ELLIPSE;
		    ellipse_list.push_back(ellipse);		    
		}
		// right
		if ( (ellipse_list[k].c1 <= -L+margin_L) && (ellipse_list[k].c1 >= -L)
		     && (ellipse_list[k].c2 <= H) && (ellipse_list[k].c2 >= -H) )
		{
		    Ellipse ellipse = ellipse_list[k];
		    ellipse.c1 = ellipse.c1 + 2.0*L;
		    ellipse.ellipse_id = ellipse.ellipse_id + 2*MAX_NUM_ELLIPSE;
		    ellipse_list.push_back(ellipse);		    
		}
		// top
		if ( //(ellipse_list[k].c1 <= L) && (ellipse_list[k].c1 >= 0) )
		    (ellipse_list[k].c2 <= -H+margin_H) && (ellipse_list[k].c2 >= -H) )
		{
		    Ellipse ellipse = ellipse_list[k];
		    ellipse.c2 = ellipse.c2 + 2.0*H;
		    ellipse.ellipse_id = ellipse.ellipse_id + 3*MAX_NUM_ELLIPSE;
		    ellipse_list.push_back(ellipse);		    
		}
		// down
		if ( //(ellipse_list[k].c1 <= L) && (ellipse_list[k].c1 >= 0) )
		    (ellipse_list[k].c2 <= H) && (ellipse_list[k].c2 >= H - margin_H) )
		{
		    Ellipse ellipse = ellipse_list[k];
		    ellipse.c2 = ellipse.c2 - 2.0*H;
		    ellipse.ellipse_id = ellipse.ellipse_id + 4*MAX_NUM_ELLIPSE;
		    ellipse_list.push_back(ellipse);		    
		}
	    }
	    return 0;
	}


    int apply_boundary_condition_voronoi()
	{
	    double L = boundary_data.L;
	    double H = boundary_data.H;
	    double margin_L = boundary_data.margin_L;
	    double margin_H = boundary_data.margin_H;

	    for ( int k=0; k<voronoi_cell_list.size(); k++ )
	    {
		if ( voronoi_cell_list[k].cell_id < 8000 )
		{
		    if ( (voronoi_cell_list[k].ellipse.c1 < -L) )
		    {
			voronoi_cell_list[k].ellipse.c1 = voronoi_cell_list[k].ellipse.c1 + 2.0*L;
		    }
		    if ( (voronoi_cell_list[k].ellipse.c1 > L) )
		    {
			voronoi_cell_list[k].ellipse.c1 = voronoi_cell_list[k].ellipse.c1 - 2.0*L;
		    }
		    if ( (voronoi_cell_list[k].ellipse.c2 < -H) )
		    {
			voronoi_cell_list[k].ellipse.c2 = voronoi_cell_list[k].ellipse.c2 + 2.0*H;
		    }
		    if ( (voronoi_cell_list[k].ellipse.c2 > H) )
		    {
			voronoi_cell_list[k].ellipse.c2 = voronoi_cell_list[k].ellipse.c2 - 2.0*H;
		    }
		    // //if ( (voronoi_cell_list[k].ellipse.a > 1.0) )
		    // {
		    // 	voronoi_cell_list[k].ellipse.a = 0.5;
		    // }
		    // //if ( (voronoi_cell_list[k].ellipse.b > 1.0) )
		    // {
		    // 	voronoi_cell_list[k].ellipse.b = 0.5;
		    // }
		}
	    }
	    return 0;
	}

    int init_ellipse_list()
	{
	    //std::ifstream finit("initPos.txt");
	    std::ifstream finit("cellPos.txt");
	    if(!finit) 
	    {
		std::cout << "file open error.\n";
		return -1; 
	    }
	    
	    int id = 0;
	    double cx, cy;
	    //std::string temp_str;
	    char temp_str;
	    while (finit >> cx >> temp_str >>  cy)
	    {
		id = id + 1;
		//if ( id != 110 && id != 88 )
		{
		    double v1 = 1.5;
		    double v2 = 0.0;
		    double u1 = 0.0;
		    double u2 = 1.0;
		    //double a = sqrt(v1*v1 + v2*v2);
		    //double b = sqrt(u1*u1 + u2*u2);
		    //std::cout << " > ellipse: " << cx << " " << cy << " " << ellipse_list.size() << "\n";
		    Ellipse ellipse(id, cx, cy, v1, v2, u1, u2);
		    ellipse.area = PI*ellipse.a*ellipse.b;
		    ellipse.theta0 = 0.0;
		    ellipse_list.push_back(ellipse);
		}
	    }

	    // for ( int i=0; i<3; i++ )
	    // {
	    // 	for ( int j=0; j<3; j++ )
	    // 	{
	    // 	    id = id + 1;
	    // 	    double cx = 6.0 + 12.0*j;
	    // 	    double cy = 6.0 + 12.0*i;
	    // 	    double v1 = 4.0;
	    // 	    double v2 = 0.0;
	    // 	    double u1 = 0.0;
	    // 	    double u2 = 4.0;
	    // 	    if (i==0 && j==0)
	    // 	    {
	    // 		cx = cx + 4;
	    // 		cy = cy + 2;
	    // 	    }
	    // 	    if (i==0 && j==2)
	    // 	    {
	    // 		cx = cx - 1;
	    // 	    }
	    // 	    if (i==2 && j==0)
	    // 	    {
	    // 		cx = cx + 2;
	    // 	    }
	    // 	    std::cout << " > ellipse: " << cx << " " << cy << " " << ellipse_list.size() << "\n";
	    // 	    Ellipse ellipse(id, cx, cy, v1, v2, u1, u2);
	    // 	    ellipse_list.push_back(ellipse);
	    // 	}
	    // }

	    finit.close();
	    return 0;
	}

    int ellipse_to_voronoi()
	{
	    int n = boundary_data.n;
	    int m = boundary_data.m;
	    int *distance_table;
	    int *partition_table;
	    std::vector<PixelNode> node_list;
	    std::unordered_set<std::unordered_set<int>, ColorHash> node_color_set;

	    distance_table = (int *)malloc(n*m*sizeof(int));
	    partition_table = (int *)malloc(n*m*sizeof(int));

	    init_distance_table(ellipse_list, partition_table, distance_table);
	    //output_partition_table(n, m, distance_table, partition_table, node_list);
	    distance_mapping(partition_table, distance_table);
	    //output_partition_table(n, m, distance_table, partition_table, node_list);

	    //find_nodes(partition_table, node_list, raw_node_list, node_color_set);

	    voronoi_node_list.clear();
	    voronoi_edge_list.clear();
	    voronoi_cell_list.clear();
	    primary_edge_set.clear();

	    make_node_list(node_list, node_color_set, partition_table, distance_table);
	    refine_node_list();
	    make_cell_list(ellipse_list);
	    make_edge_list(); // must after "make_voronoi_cell_list"
	    refine_cell_list();

	    //write_ellipse_tofile();
	    //output_node_list();
	    //output_edge_list();
	    //output_cell_list();

	    delete(distance_table);
	    delete(partition_table);
	    return 0;
	}

    Ellipse voronoi_to_ellipse(VoronoiCell &cell)
	{
	    Ellipse new_ellipse = cell.ellipse;
	    double I1 = 0.0;
	    double I2 = 0.0;
	    double I3 = 0.0;
	    double c1 = 0.0;
	    double c2 = 0.0;
	    for ( int kk=0; kk<cell.node_set.size(); kk++ )
	    {
		c1 = c1 + voronoi_node_list[cell.node_set[kk]].x;
		c2 = c2 + voronoi_node_list[cell.node_set[kk]].y;
	    }
	    c1 = c1/cell.node_set.size();
	    c2 = c2/cell.node_set.size();

	    for ( int kk=1; kk<cell.node_set.size(); kk++ )
	    {
		double x1 = voronoi_node_list[cell.node_set[kk-1]].x - c1;
		double y1 = voronoi_node_list[cell.node_set[kk-1]].y - c2;
		double x2 = voronoi_node_list[cell.node_set[kk]].x - c1;
		double y2 = voronoi_node_list[cell.node_set[kk]].y - c2;
		I1 = I1 + (-(x2-x1)/12.0)*(y2*y2*y2 + y2*y2*y1 + y2*y1*y1 + y1*y1*y1);
		I2 = I2 + ((x2-x1)/24.0)*(x1*(3.0*y1*y1 + 2.0*y1*y2 + y2*y2) + 
					  x2*(3.0*y2*y2 + 2.0*y1*y2 + y1*y1));
		I3 = I3 + ((y2-y1)/12.0)*(x2*x2*x2 + x2*x2*x1 + x2*x1*x1 + x1*x1*x1); 
		//std::cout << x1 << ", " << y1 << ", ";
	    }		
	    int kk = cell.node_set.size()-1;
	    double x1 = voronoi_node_list[cell.node_set[kk]].x - c1;
	    double y1 = voronoi_node_list[cell.node_set[kk]].y - c2;
	    double x2 = voronoi_node_list[cell.node_set[0]].x - c1;
	    double y2 = voronoi_node_list[cell.node_set[0]].y - c2;
	    I1 = I1 + (-(x2-x1)/12.0)*(y2*y2*y2 + y2*y2*y1 + y2*y1*y1 + y1*y1*y1);
	    I2 = I2 + ((x2-x1)/24.0)*(x1*(3.0*y1*y1 + 2.0*y1*y2 + y2*y2) + 
				      x2*(3.0*y2*y2 + 2.0*y1*y2 + y1*y1));
	    I3 = I3 + ((y2-y1)/12.0)*(x2*x2*x2 + x2*x2*x1 + x2*x1*x1 + x1*x1*x1); 
	    double lambda1 = (I1 + I3 + sqrt(I1*I1 + I3*I3 + 4.0*I2*I2 - 2.0*I1*I3))/2.0;
	    double lambda2 = (I1 + I3 - sqrt(I1*I1 + I3*I3 + 4.0*I2*I2 - 2.0*I1*I3))/2.0;
	    // note that the order is different
	    double u1 = lambda1 - I3;
	    double u2 = I2;
	    double v1 = lambda2 - I3;
	    double v2 = I2;
	    if ( I2 == 0.0 )
	    {
		if ( fabs(lambda1 - I1) < 0.00001 )
		{
		    u1 = 1.0;
		    u2 = 0.0;
		    v1 = 0.0;
		    v2 = 1.0;
		}
		else if ( fabs(lambda2 - I1) < 0.00001 )
		{
		    u1 = 0.0;
		    u2 = 1.0;
		    v1 = 1.0;
		    v2 = 0.0;
		}
		else
		{
		    std::cout << "wrong eigenvalue\n";
		    std::cout << lambda1 << " " << lambda2 << " " << I1 << " " << I3 << "\n";
		    getchar();
		}
	    } 
	    if ( lambda1 <= 0 || lambda2 <= 0 )
	    {
		std::cout << "something wrong with ellipse \n";
		getchar();
	    }
	    double ratio_temp = lambda1 / lambda2;
	    double multi_temp = new_ellipse.a * new_ellipse.b;
	    // if ( ratio_temp > 2.0 )
	    // {
	    // 	ratio_temp = 2.0;
	    // }
	    // if ( ratio_temp < 1.0 )
	    // {
	    // 	ratio_temp = 1.0;
	    // }

	    lambda1 = sqrt(multi_temp*ratio_temp);
	    lambda2 = sqrt(multi_temp/ratio_temp);
	    // new_ellipse.c1 = c1;
	    // new_ellipse.c2 = c2;
	    new_ellipse.a = 1.0/sqrt(lambda2);
	    new_ellipse.b = 1.0/sqrt(lambda1);
	    new_ellipse.v1 = v1/sqrt(v1*v1 + v2*v2);
	    new_ellipse.v2 = v2/sqrt(v1*v1 + v2*v2);
	    new_ellipse.u1 = u1/sqrt(u1*u1 + u2*u2);
	    new_ellipse.u2 = u2/sqrt(u1*u1 + u2*u2);

	    //std::cout << cell.cell_id << ", " << lambda1 << ", " << lambda2 << "\n";
	    return new_ellipse;
	}

    int refine_cell_list()
	{
	    //for ( int k=0; k<voronoi_cell_list.size(); k++ )	    
	    //std::cout << "refining \n";
	    for ( std::vector<VoronoiCell>::iterator it=voronoi_cell_list.begin(); 
		  it!=voronoi_cell_list.end(); )	    
	    {
		std::vector<int> edges = it->edge_set;
		std::vector<int> refined_edges;
		for ( int kk=0; kk<edges.size(); kk++ )
		{
		    if ( voronoi_edge_list[edges[kk]].side_colors.size() == 2 )
		    {
		    	if ( voronoi_edge_list[edges[kk]].side_colors[0] == it->cell_id )
		    	{
		    	    it->neighbor_id.push_back(voronoi_edge_list[edges[kk]].side_colors[1]);
		    	}
		    	else
		    	{
		    	    it->neighbor_id.push_back(voronoi_edge_list[edges[kk]].side_colors[0]);
		    	}
			refined_edges.push_back(edges[kk]);
		    }
		    else
		    {
			//std::cout << "wrong edge here " << it->cell_id << "-" << kk << "\n";
			//voronoi_cell_list[k].edge_set.erase(voronoi_cell_list[k].edge_set.begin() + kk);
		    }
		}
		if ( (edges.size() != refined_edges.size()) || ( refined_edges.size() < 3 ) )
		{
		    //output_cell(*it);
		    if ( it->cell_id < 8000 )
		    {
			std::cout << "remove wrong cell? **** \n";
			output_cell(*it);
			for ( int kk=0; kk<edges.size(); kk++ )
			{
			    if ( voronoi_edge_list[edges[kk]].side_colors.size() != 2 )
			    {
				std::cout << "at edge " << edges[kk] << " side colors: "
					  << voronoi_edge_list[edges[kk]].side_colors.size() << "id: "
					  << voronoi_edge_list[edges[kk]].edge_id
					  << "( node1: " << voronoi_edge_list[edges[kk]].node_id1
					  << "(" << voronoi_edge_list[edges[kk]].x1
					  << ", " << voronoi_edge_list[edges[kk]].y1
					  << ") - node2: " << voronoi_edge_list[edges[kk]].node_id2
					  << "(" << voronoi_edge_list[edges[kk]].x2
					  << ", " << voronoi_edge_list[edges[kk]].y2
					  << ") )\n";
				output_node(voronoi_node_list[voronoi_edge_list[edges[kk]].node_id1]);
				output_node(voronoi_node_list[voronoi_edge_list[edges[kk]].node_id2]);
			    }
			}
			getchar();
		    }		    
		    voronoi_cell_list.erase(it);
		    //voronoi_cell_list[k].area = -1.0;
		    //it->is_valid = false;
		    //++it;
		}
		else
		{
		    it->is_valid = true;
		    double area = 0.0;
		    it->edges.clear();
		    it->nodes.clear();
		    it->nodes.push_back(voronoi_node_list[it->node_set[0]]);
		    for ( int kk=1; kk<it->node_set.size(); kk++ )
		    {
			int node_id1 = it->node_set[kk-1];
			int node_id2 = it->node_set[kk];
			double x1 = voronoi_node_list[node_id1].x;
			double y1 = voronoi_node_list[node_id1].y;
			double x2 = voronoi_node_list[node_id2].x;
			double y2 = voronoi_node_list[node_id2].y;
			it->nodes.push_back(voronoi_node_list[node_id2]);

			VoronoiEdge current_edge =  voronoi_edge_list[it->edge_set[kk-1]];
			double n1 = (y2-y1)/sqrt( (y2-y1)*(y2-y1) + (x2-x1)*(x2-x1) );
			double n2 = -(x2-x1)/sqrt( (y2-y1)*(y2-y1) + (x2-x1)*(x2-x1) );
			current_edge.n1 = n1;
			current_edge.n2 = n2;
			current_edge.theta = atan2(current_edge.c2 - it->c2, current_edge.c1 - it->c1);
			current_edge.theta1 = atan2(y1 - it->c2, x1 - it->c1);
			current_edge.theta2 = atan2(y2 - it->c2, x2 - it->c1);
			if ( current_edge.theta1 >= current_edge.theta2 )
			{
			    std::cout << "check theta1 and theta2 ? theta1:" <<
				current_edge.theta1 << " theta2: " << 
				current_edge.theta2 << "\n";
			    getchar();
			}
			double tri_sector_area = get_tri_sector_area(it->c1, it->c2, x1, y1, x2, y2);
			double ell_sector_area = get_ell_sector_area(it->ellipse, 
								     current_edge.theta1, 
								     current_edge.theta2);
			
			current_edge.tri_sector_area = tri_sector_area;
			current_edge.ell_sector_area = ell_sector_area;
			it->edges.push_back(current_edge);

			// check consistency
			if ( (x1 != current_edge.x1) && (x1 != current_edge.x2) )
			{
			    std::cout << "something wrong? " << x1 << " " <<
				x2 << "edge xy " << current_edge.x1 << " " <<
				current_edge.x2 << "size: " << it->edge_set.size() << "\n";
			    
			    output_cell(*it);
			    getchar();
			}

			area = area + (x1*y2 - x2*y1);
		    }
		    int node_id1 = it->node_set[it->node_set.size()-1];
		    int node_id2 = it->node_set[0];
		    double x1 = voronoi_node_list[node_id1].x;
		    double y1 = voronoi_node_list[node_id1].y;
		    double x2 = voronoi_node_list[node_id2].x;
		    double y2 = voronoi_node_list[node_id2].y;

		    VoronoiEdge current_edge =  voronoi_edge_list[it->edge_set[it->edge_set.size()-1]];
		    double n1 = (y2-y1)/sqrt( (y2-y1)*(y2-y1) + (x2-x1)*(x2-x1) );
		    double n2 = -(x2-x1)/sqrt( (y2-y1)*(y2-y1) + (x2-x1)*(x2-x1) );
		    current_edge.n1 = n1;
		    current_edge.n2 = n2;
		    current_edge.theta = atan2(current_edge.c2, current_edge.c1);
		    current_edge.theta1 = atan2(y1 - it->c2, x1 - it->c1);
		    current_edge.theta2 = 2.0*PI + atan2(y2 - it->c2, x2 - it->c1);
		    if ( current_edge.theta1 >= current_edge.theta2 )
		    {
			std::cout << "check theta1 and theta2 ? theta1:" <<
			    current_edge.theta1 << " theta2: " << 
			    current_edge.theta2 << "\n";
			getchar();
		    }
		    double tri_sector_area = get_tri_sector_area(it->c1, it->c2, x1, y1, x2, y2);
		    double ell_sector_area = get_ell_sector_area(it->ellipse, 
								 current_edge.theta1, 
								 current_edge.theta2);
		    
		    current_edge.tri_sector_area = tri_sector_area;
		    current_edge.ell_sector_area = ell_sector_area;
		    it->edges.push_back(current_edge);
		    // check consistency
		    if ( (x1 != current_edge.x1) && (x1 != current_edge.x2) )
		    {
			std::cout << "something wrong 2? " << x1 << " " <<
			    x2 << "edge xy " << current_edge.x1 << " " <<
			    current_edge.x2 << "size: " << it->edge_set.size() << "\n";
			
			output_cell(*it);
			getchar();
		    }

		    area = area + (x1*y2 - x2*y1);
		    it->area = 0.5*fabs(area);
		    // check consistency for area
		    double sum_area1 = 0.0;
		    double sum_area2 = 0.0;
		    //std::cout << "ellipse part area ";
		    for ( int kk=0; kk < it->edges.size(); kk++ )
		    {
			sum_area1 = sum_area1 + it->edges[kk].tri_sector_area;
			sum_area2 = sum_area2 + it->edges[kk].ell_sector_area;
			//std::cout << it->edges[kk].ell_sector_area << " ";
			//std::cout << it->edges[kk].theta1 << " " << it->edges[kk].theta2 << " ";
		    }
		    //std::cout << "\n";
		    if ( fabs(sum_area1 - it->area)/it->area > 0.25 )
		    {
			std::cout << "total area does not add up!\n";
			std::cout << "total area: " << it->area
				  << "sumed area: " << sum_area1 << "\n";
			getchar();
		    }
		    if ( fabs(sum_area2 - it->ellipse.area)/it->ellipse.area > 0.25 )
		    {
			std::cout << "ellipse area does not add up!\n";
			std::cout << "ellipse area: " << it->ellipse.area
				  << "sumed area: " << sum_area2
				  << it->ellipse.a << " " << it->ellipse.b << "\n";
			getchar();
		    }
		    ++it;
		}
	    }
	    return 0;
	}

    int make_node_list(const std::vector<PixelNode> &node_list, 
		       std::unordered_set<std::unordered_set<int>, ColorHash> &node_color_set,
		       int* partition_table,
		       const int* distance_table)
	{
	    int n = boundary_data.n;
	    int m = boundary_data.m;
	    int nmin = boundary_data.nmin;
	    int nmax = boundary_data.nmax;
	    int mmin = boundary_data.mmin;
	    int mmax = boundary_data.mmax;

	    double xmin = boundary_data.xmin;
	    double xmax = boundary_data.xmax;
	    double ymin = boundary_data.ymin;
	    double ymax = boundary_data.ymax;

	    double dx = boundary_data.res_x;
	    double dy = boundary_data.res_y;

	    if ( dx<=0 || dy<=0 )
	    {
		std::cout << "make_voronoi_node_list: parameter wrong \n";
		getchar();
	    }
	    voronoi_node_list.clear();
	    int nearby_values[8];

	    for ( int i=nmin; i<nmax; i++ )
	    {
		for ( int j=mmin; j<mmax; j++ )
		{
		    nearby_values[0] = partition_table[(i-1)*m+j];
		    nearby_values[1] = partition_table[(i-1)*m+j-1];
		    nearby_values[2] = partition_table[(i)*m+j-1];
		    nearby_values[3] = partition_table[(i+1)*m+j-1];
		    nearby_values[4] = partition_table[(i+1)*m+j];
		    nearby_values[5] = partition_table[(i+1)*m+j+1];
		    nearby_values[6] = partition_table[(i)*m+j+1];
		    nearby_values[7] = partition_table[(i-1)*m+j+1];

		    const size_t len = sizeof(nearby_values) / sizeof(nearby_values[0]);
		    std::unordered_set<int> node_color(nearby_values, nearby_values+len);
		    node_color.erase( -1 );
		    //if ( (node_color.size() >= 3) && 
		    if ( (node_color.size() > 5) &&
			 (node_color_set.find(node_color) == node_color_set.end()) )
		    {
			std::cout << "a high-fold node: " << node_color.size() << " colors: ";
			getchar();
		    }
		}
	    }

	    for ( int i=nmin; i<nmax; i++ )
	    {
		for ( int j=mmin; j<mmax; j++ )
		{
		    nearby_values[0] = partition_table[(i-1)*m+j];
		    nearby_values[1] = partition_table[(i-1)*m+j-1];
		    nearby_values[2] = partition_table[(i)*m+j-1];
		    nearby_values[3] = partition_table[(i+1)*m+j-1];
		    nearby_values[4] = partition_table[(i+1)*m+j];
		    nearby_values[5] = partition_table[(i+1)*m+j+1];
		    nearby_values[6] = partition_table[(i)*m+j+1];
		    nearby_values[7] = partition_table[(i-1)*m+j+1];

		    const size_t len = sizeof(nearby_values) / sizeof(nearby_values[0]);
		    std::unordered_set<int> node_color(nearby_values, nearby_values+len);
		    node_color.erase( -1 );
		    //if ( (node_color.size() >= 3) && 
		    if ( (node_color.size() == 5) &&
			 (node_color_set.find(node_color) == node_color_set.end()) )
		    {
			// std::cout << "new node_color: size " << node_color.size() << " colors: ";
			// for (std::unordered_set<int>::const_iterator it = node_color.begin();
			//      it != node_color.end(); it++ )
			// {
			//     std::cout << *it << " ";
			// }
			// std::cout << "\n";

			node_color_set.insert(node_color);

			VoronoiNode voronoi_node;
			//int node_id = i*m + j;
			//id = id + 1;
			//i = node_list[k].irow;
			//j = node_list[k].icol;

			//voronoi_node.node_id = node_id;
			voronoi_node.color_set = node_color;
			voronoi_node.x = xmin + (j+0.5)*dx;
			voronoi_node.y = ymin + (i+0.5)*dy;
		
			voronoi_node_list.push_back(voronoi_node);
			// partition_table[(i-1)*m+j] = 99799;
			// partition_table[(i-1)*m+j-1] = 99799;
			// partition_table[(i)*m+j-1] = 99799;
			// partition_table[(i+1)*m+j-1] = 99799;
			// partition_table[(i+1)*m+j] = 99799;
			// partition_table[(i+1)*m+j+1] = 99799;
			// partition_table[(i)*m+j+1] = 99799;
			// partition_table[(i-1)*m+j+1] = 99799;
			// partition_table[i*m+j] = 99799;
		    }
		}
	    }

	    for ( int i=nmin; i<nmax; i++ )
	    {
		for ( int j=mmin; j<mmax; j++ )
		{
		    nearby_values[0] = partition_table[(i-1)*m+j];
		    nearby_values[1] = partition_table[(i-1)*m+j-1];
		    nearby_values[2] = partition_table[(i)*m+j-1];
		    nearby_values[3] = partition_table[(i+1)*m+j-1];
		    nearby_values[4] = partition_table[(i+1)*m+j];
		    nearby_values[5] = partition_table[(i+1)*m+j+1];
		    nearby_values[6] = partition_table[(i)*m+j+1];
		    nearby_values[7] = partition_table[(i-1)*m+j+1];

		    const size_t len = sizeof(nearby_values) / sizeof(nearby_values[0]);
		    std::unordered_set<int> node_color(nearby_values, nearby_values+len);
		    node_color.erase( -1 );
		    //if ( (node_color.size() >= 3) && 
		    if ( (node_color.size() == 4) &&
			 (node_color_set.find(node_color) == node_color_set.end()) &&
			 !isSubset(node_color, node_color_set) )
		    {
			// std::cout << "new node_color: size " << node_color.size() << " colors: ";
			// for (std::unordered_set<int>::const_iterator it = node_color.begin();
			//      it != node_color.end(); it++ )
			// {
			//     std::cout << *it << " ";
			// }
			// std::cout << "\n";

			node_color_set.insert(node_color);

			VoronoiNode voronoi_node;
			//int node_id = i*m + j;
			//id = id + 1;
			//i = node_list[k].irow;
			//j = node_list[k].icol;

			//voronoi_node.node_id = node_id;
			voronoi_node.color_set = node_color;
			voronoi_node.x = xmin + (j+0.5)*dx;
			voronoi_node.y = ymin + (i+0.5)*dy;
		
			voronoi_node_list.push_back(voronoi_node);
			// partition_table[(i-1)*m+j] = 99799;
			// partition_table[(i-1)*m+j-1] = 99799;
			// partition_table[(i)*m+j-1] = 99799;
			// partition_table[(i+1)*m+j-1] = 99799;
			// partition_table[(i+1)*m+j] = 99799;
			// partition_table[(i+1)*m+j+1] = 99799;
			// partition_table[(i)*m+j+1] = 99799;
			// partition_table[(i-1)*m+j+1] = 99799;
			// partition_table[i*m+j] = 99799;
		    }
		}
	    }

	    for ( int i=nmin; i<nmax; i++ )
	    {
		for ( int j=mmin; j<mmax; j++ )
		{
		    nearby_values[0] = partition_table[(i-1)*m+j];
		    nearby_values[1] = partition_table[(i-1)*m+j-1];
		    nearby_values[2] = partition_table[(i)*m+j-1];
		    nearby_values[3] = partition_table[(i+1)*m+j-1];
		    nearby_values[4] = partition_table[(i+1)*m+j];
		    nearby_values[5] = partition_table[(i+1)*m+j+1];
		    nearby_values[6] = partition_table[(i)*m+j+1];
		    nearby_values[7] = partition_table[(i-1)*m+j+1];

		    const size_t len = sizeof(nearby_values) / sizeof(nearby_values[0]);
		    std::unordered_set<int> node_color(nearby_values, nearby_values+len);
		    node_color.erase( -1 );

		    if ( (node_color.size() == 3) && 
			 //isLocalMaximum(i, j, n, m, distance_table) &&
			 (node_color_set.find(node_color) == node_color_set.end()) &&
			 !isSubset(node_color, node_color_set) )
		    {
			// std::cout << "new node_color: ";
			// for (std::unordered_set<int>::const_iterator it = node_color.begin();
			//      it != node_color.end(); it++ )
			// {
			//     std::cout << *it << " ";
			// }
			// std::cout << "\n";

			node_color_set.insert(node_color);

			VoronoiNode voronoi_node;
			//int node_id = i*m + j;
			//id = id + 1;
			//i = node_list[k].irow;
			//j = node_list[k].icol;

			//voronoi_node.node_id = node_id;
			voronoi_node.color_set = node_color;
			voronoi_node.x = xmin + (j+0.5)*dx;
			voronoi_node.y = ymin + (i+0.5)*dy;
		
			voronoi_node_list.push_back(voronoi_node);
			// partition_table[(i-1)*m+j] = 99799;
			// partition_table[(i-1)*m+j-1] = 99799;
			// partition_table[(i)*m+j-1] = 99799;
			// partition_table[(i+1)*m+j-1] = 99799;
			// partition_table[(i+1)*m+j] = 99799;
			// partition_table[(i+1)*m+j+1] = 99799;
			// partition_table[(i)*m+j+1] = 99799;
			// partition_table[(i-1)*m+j+1] = 99799;
			// partition_table[i*m+j] = 99799;
		    }
		    // std::unordered_set<int> set1;
		    // set1.insert(3);
		    // set1.insert(6);
		    // std::unordered_set<int> set2;
		    // set2.insert(6);
		    // set2.insert(3);
		    // if ( set1 == set2 )
		    // {
		    // 	cout << "equal\n";
		    // }
		    // else
		    // {
		    // 	cout << "not equal!\n";			
		    // }
		}
	    }
	    // clean node

	    if ( node_color_set.size() != voronoi_node_list.size() )
	    {
		std::cout << "make_voronoi_node_list: unbalanced out come.\n";
		getchar();
	    }
	    

	    return 0;
	}

    int refine_node_list()
	{
	    double dx = boundary_data.res_x;
	    double dy = boundary_data.res_y;
	    if ( dx<=0 || dy<=0 )
	    {
		std::cout << "make_voronoi_node_list: parameter wrong \n";
		getchar();
	    }

	    double dist_up = 2.0*(dx+dy);
	    double dist_low = 0.1*(dx+dy);

	    std::vector<VoronoiNode> new_node_list;
	    for ( std::vector<VoronoiNode>::iterator it=voronoi_node_list.begin(); 
		  it!=voronoi_node_list.end(); ++it)	    
	    {
		int found = 0;
		for ( int k=0; k<new_node_list.size(); k++ )
		{
		    double dist = fabs(it->x - new_node_list[k].x) + fabs(it->y - new_node_list[k].y);
		    if ( (dist > dist_low) && (dist < dist_up) )
		    {
			// add the color_set of it to node k
			for ( std::unordered_set<int>::const_iterator it_color = it->color_set.begin();
			      it_color != it->color_set.end(); it_color++ )
			{
			    new_node_list[k].color_set.insert(*it_color);
			}
			found = 1;
			//std::cout << "remove nearby nodes and add colors\n";
			//std::cout << "old color_set \n";
			/* for (std::unordered_set<int>::const_iterator it2 = it->color_set.begin(); */
			/*      it2 != it->color_set.end(); it2++ ) */
			/* { */
			/*     std::cout << *it2 << " "; */
			/* } */
			//std::cout << "\n";
			//std::cout << "new color_set \n";
			/* for (std::unordered_set<int>::const_iterator it3 = new_node_list[k].color_set.begin(); */
			/*      it3 != new_node_list[k].color_set.end(); it3++ ) */
			/* { */
			/*     std::cout << *it3 << " "; */
			/* } */
			//getchar();
		    }
		}
		if ( found == 0 )
		{
		    new_node_list.push_back(*it);
		}
	    }
	    /* if (voronoi_node_list.size() != new_node_list.size() ) */
	    /* { */
	    /* 	std::cout << "old node list size: " << voronoi_node_list.size() << "\n"; */
	    /* 	std::cout << "new node list size: " << new_node_list.size() << "\n"; */
	    /* 	//getchar(); */
	    /* } */
	    voronoi_node_list = new_node_list;
	    return 0;
	}

    bool isLocalMaximum(const int i, const int j, 
			const int n, const int m, 
			const int * distance_table)
	{
	    return ( (distance_table[i*m+j] >= distance_table[(i-1)*m+j]) &&
		     (distance_table[i*m+j] >= distance_table[(i+1)*m+j]) &&
		     (distance_table[i*m+j] >= distance_table[i*m+j+1]) &&
		     (distance_table[i*m+j] >= distance_table[i*m+j-1]) );
	}

    bool isSubset(const std::unordered_set<int> test_set,
		  const std::unordered_set<std::unordered_set<int>, ColorHash> set_list)
	{
	    bool is_found = false;
	    for ( std::unordered_set<std::unordered_set<int>, ColorHash>::const_iterator it = set_list.begin();
		  it != set_list.end(); it++ )
	    {
		is_found = true;
		for ( std::unordered_set<int>::const_iterator it_test = test_set.begin();
		      it_test != test_set.end(); it_test++ )
		{
		    if ( it->find(*it_test) == it->end() )
		    {
			is_found = false;
			break;
		    }
		}
		if ( is_found )
		{
		    //std::cout << "found subset:" << test_set.size() << " in " << it->size() << "\n";
		    //getchar();
		    return is_found;
		}
	    }
	    return is_found;
	}

    int make_cell_list(std::vector<Ellipse> &ellipse_list)
	{
	    int cell_id;
	    for ( int k=0; k<ellipse_list.size(); k++ )	    
	    {
		VoronoiCell voronoi_cell;
		std::vector<int> node_set;
		std::vector<int> edge_set;
		node_set.clear();
		cell_id = ellipse_list[k].ellipse_id;
		voronoi_cell.ellipse = ellipse_list[k];
		voronoi_cell.cell_id = cell_id;
		voronoi_cell.c1 = ellipse_list[k].c1;
		voronoi_cell.c2 = ellipse_list[k].c2;

		//std::cout << "make_voronoi_cell_list: cell_id" << cell_id << ":";
		 
		for ( int kk=0; kk<voronoi_node_list.size(); kk++ )
		{
		    std::unordered_set<int>::iterator got = 
			voronoi_node_list[kk].color_set.find(cell_id);
		    if ( got != voronoi_node_list[kk].color_set.end() )
		    {	
			//std::cout << "find node\n"; 
			// we do not use node_id, simple use array index AS its id.
			//node_set.push_back(voronoi_node_list[kk].node_id);
			node_set.push_back(kk);
		    }
		}
		std::vector<int> sorted_node_set = 
		    voronoi_arrange_node(ellipse_list[k].c1, ellipse_list[k].c2, node_set);
		voronoi_cell.node_set = sorted_node_set;

		edge_set = prepare_edges(sorted_node_set);

		voronoi_cell.edge_set = edge_set;

		voronoi_cell_list.push_back(voronoi_cell);
	    }
	    return 0;
	}

    int make_edge_list()
	{
	    voronoi_edge_list.clear();
	    //make sure the edge_id is the same as the index in edge_list
	    VoronoiEdge dummy(0,0,0);
	    voronoi_edge_list.resize(primary_edge_set.size(), dummy);
	    for ( std::unordered_set<VoronoiEdge, Hash>::iterator it = primary_edge_set.begin();
		  it != primary_edge_set.end(); it++ )
	    {
		VoronoiEdge edge(it->edge_id, it->node_id1, it->node_id2);
		VoronoiNode node1 = voronoi_node_list[it->node_id1];
		VoronoiNode node2 = voronoi_node_list[it->node_id2];
		double length = sqrt((node1.x - node2.x)*(node1.x - node2.x) 
				     + (node1.y - node2.y)*(node1.y - node2.y));
		edge.length = length;
		edge.x1 = node1.x;
		edge.y1 = node1.y;
		edge.x2 = node2.x;
		edge.y2 = node2.y;
		edge.c1 = 0.5*(edge.x1 + edge.x2);
		edge.c2 = 0.5*(edge.y1 + edge.y2);
		edge.n1 = 0.0;
		edge.n2 = 0.0;

		std::unordered_set<int> color_set1 = node1.color_set;
		std::unordered_set<int> color_set2 = node2.color_set;

		std::vector<int> side_colors;
		find_common_colors(color_set1, color_set2, side_colors);

		if ( side_colors.size() == 1 )
		{
		    side_colors.clear();
		}
		else if ( side_colors.size() == 2 )
		{
		    // do nothing
		}
		else
		{
		    std::cout << "make_edge_list: something wrong.\n";
		    std::cout << side_colors.size() << "\n";
		    std::cout << edge.x1 << " " << edge.x2 << "\n";
		    std::cout << edge.y1 << " " << edge.y2 << "\n";
		    std::cout << it->node_id1 << " " << it->node_id2 << "\n";
		    for (std::unordered_set<int>::const_iterator it = color_set1.begin();
			 it != color_set1.end(); it++ )
		    {
			std::cout << *it << " ";
		    }
		    std::cout << "\n";
		    for (std::unordered_set<int>::const_iterator it = color_set2.begin();
			 it != color_set2.end(); it++ )
		    {
			std::cout << *it << " ";
		    }
		    std::cout << "\n";
		    side_colors.clear();
		    getchar();
		}
		edge.side_colors = side_colors;
		//make sure the edge_id is the same as the index in edge_list
		voronoi_edge_list[it->edge_id] = edge;
		//std::cout << it->edge_id << " " << it->node_id1 << " " << it->node_id2 << "\n";
	    }
	    //std::cout << "voronoi edge list size: " << voronoi_edge_list.size() << "\n";
	    for ( int k=0; k<voronoi_edge_list.size(); k++ )
	    {
		if ( voronoi_edge_list[k].edge_id != k )
		{
		    std::cout << "voronoi edge list wrong id: \n";
		    getchar();
		}
	    }
	    return 0;
	}

    int find_common_colors( const std::unordered_set<int> &color_set1, 
			    const std::unordered_set<int> &color_set2,
			    std::vector<int> &side_colors)
	{
	    int counter = 0;
	    side_colors.clear();
	    //std::cout << "\n";
	    for ( std::unordered_set<int>::const_iterator it=color_set1.begin(); 
		  it!=color_set1.end(); it++)
	    {
		//std::cout << "find " << *it << "in ";
		for ( std::unordered_set<int>::const_iterator it2=color_set2.begin(); 
		      it2!=color_set2.end(); it2++)
		{
		    //std::cout << *it2 << " ";
		}
		//std::cout << "\n";

		if ( color_set2.find(*it) != color_set2.end() )
		{
		    side_colors.push_back(*it);
		    counter = counter + 1;
		}
	    }
	    // if ( counter != 2 )
	    // {
	    // 	std::cout << "find_common_colors: wrong number of colors: " << counter << "\n";
	    // 	getchar();
	    // }
	    return 0;
	}

    std::vector<int> prepare_edges(const std::vector<int> &node_set)
	{
	    std::vector<int> edge_set;
	    edge_set.clear();
	    int start_node;
	    int end_node;
	    int node_index1;
	    int node_index2;
	    if ( node_set.size() > 1 )
	    {
		for ( int k=0; k<node_set.size(); k++ )
		{
		    if ( k == node_set.size()-1 )
		    {
		    	node_index1 = node_set[0];
		    	node_index2 = node_set[k];
		    }
		    else
		    {
		    	node_index1 = node_set[k+1];
		    	node_index2 = node_set[k];
		    }
		    if ( node_index1 > node_index2 )
		    {
			start_node = node_index2;
			end_node = node_index1;
		    }
		    else
		    {
			end_node = node_index2;
			start_node = node_index1;
		    }
		    VoronoiEdge edge(primary_edge_set.size(), start_node, end_node);
		    primary_edge_set.insert(edge);
		    std::unordered_set<VoronoiEdge, Hash>::iterator got = primary_edge_set.find(edge);
		    edge_set.push_back(got->edge_id);
		}
	    }
	    return edge_set;
	}

    std::vector<int> voronoi_arrange_node(const double cx, 
					  const double cy,
					  std::vector<int> &node_set)
	{
	    std::vector<double> angles;
	    angles.clear();
	    for ( int k=0; k<node_set.size(); k++ )
	    {
		double azimuth = atan2(voronoi_node_list[node_set[k]].y - cy, 
				       voronoi_node_list[node_set[k]].x - cx); 
		angles.push_back(azimuth);
		//std::cout << "voronoi_arrange_node: " << azimuth << ":";
	    }
	    //std::cout << "\n";
	    std::vector<size_t> idx = sort_indexes(angles);
	    std::vector<int> sorted_node_set(node_set.size());
	    for ( int k=0; k<node_set.size(); k++ )
	    {
		sorted_node_set[k] = node_set[idx[k]];
	    }	    
	    return sorted_node_set;
	}
    
    int update_ellipse_list()
	{
	    ellipse_list.clear();
	    for ( int k=0; k<voronoi_cell_list.size(); k++)
	    {
		if ( voronoi_cell_list[k].cell_id < 8000 )
		{
		    ellipse_list.push_back(voronoi_cell_list[k].ellipse);
		    //ellipse_list.push_back(voronoi_to_ellipse(voronoi_cell_list[k]));
		}
	    }
	    return 0;
	}

    int output_ellipse_list(std::string file_index)
	{
	    std::ofstream fellipse("out/ellipses" + file_index + ".txt");
	    if(!fellipse) 
	    {
		std::cout << "file open error.\n";
		return -1; 
	    }

	    for ( int k=0; k<voronoi_cell_list.size(); k++)
	    {
		if ( voronoi_cell_list[k].cell_id < 8000 )
		{
		    Ellipse converted_ellipse = voronoi_to_ellipse(voronoi_cell_list[k]);
		    //double aspect_ratio = converted_ellipse.a / converted_ellipse.b ;
		    double aspect_ratio = voronoi_cell_list[k].ellipse.accumu_a /
			voronoi_cell_list[k].ellipse.accumu_b;
		    fellipse << voronoi_cell_list[k].ellipse.ellipse_id << ", " <<
			voronoi_cell_list[k].ellipse.c1 << ", " << 
			voronoi_cell_list[k].ellipse.c2 << ", " << 
			voronoi_cell_list[k].ellipse.accumu_a << ", " << 
			voronoi_cell_list[k].ellipse.accumu_b << ", " << 
			converted_ellipse.a << ", " << 
			converted_ellipse.b << ", " << 
			converted_ellipse.v1 << ", " << 
			converted_ellipse.v2 << ", " << 
			voronoi_cell_list[k].ellipse.c1 + aspect_ratio*voronoi_cell_list[k].ellipse.v1 << ", " <<
		        voronoi_cell_list[k].ellipse.c2 + aspect_ratio*voronoi_cell_list[k].ellipse.v2 << "\n";
		}
	    }

	    fellipse.close();
	    return 0;
	}

    int output_cell_position()
	{
	    std::ofstream fout("out/cellPos.txt");
	    if(!fout) 
	    {
		std::cout << "file open error.\n";
		return -1; 
	    }

	    for ( int k=0; k<voronoi_cell_list.size(); k++)
	    {
		if ( voronoi_cell_list[k].cell_id < 8000 )
		{
		    fout << voronoi_cell_list[k].ellipse.c1 << "," << 
			voronoi_cell_list[k].ellipse.c2 << "\n";
		}
	    }

	    fout.close();
	    return 0;
	}

    int write_ellipse_tofile()
	{
	    std::ofstream fellipse("out/ellipses.txt");
	    if(!fellipse) 
	    {
		std::cout << "file open error.\n";
		return -1; 
	    }

	    for ( int k=0; k<ellipse_list.size(); k++)
	    {
		fellipse << ellipse_list[k].ellipse_id << ", " <<
		    ellipse_list[k].c1 << "," << 
		    ellipse_list[k].c2 << "," <<
		    ellipse_list[k].v1 << "," << 
		    ellipse_list[k].v2 << "," <<
		    ellipse_list[k].u1 << "," << 
		    ellipse_list[k].u2 << "," <<
		    ellipse_list[k].a << "," << 
		    ellipse_list[k].b << "\n";
	    }

	    fellipse.close();
	    return 0;
	}

    int output_new_ellipse_list()
	{
	    std::ofstream fellipse("out/new_ellipses.txt");
	    if(!fellipse) 
	    {
		std::cout << "file open error.\n";
		return -1; 
	    }

	    for ( int k=0; k<new_ellipse_list.size(); k++)
	    {
		fellipse << new_ellipse_list[k].ellipse_id << ", " <<
		    new_ellipse_list[k].c1 << "," << 
		    new_ellipse_list[k].c2 << "," <<
		    new_ellipse_list[k].v1 << "," << 
		    new_ellipse_list[k].v2 << "," <<
		    new_ellipse_list[k].u1 << "," << 
		    new_ellipse_list[k].u2 << "," <<
		    new_ellipse_list[k].a << "," << 
		    new_ellipse_list[k].b << "\n";
	    }

	    fellipse.close();
	    return 0;
	}

    int output_partition_table(const int n, const int m, 
			       const int* distance_table, 
			       const int* partition_table, 
			       const std::vector<PixelNode> &node_list)
	{
	    std::ofstream fout("out/partition_table.txt");
	    std::ofstream fout2("out/distance_table.txt");
	    if(!fout || !fout2) 
	    {
		std::cout << "file open error.\n";
		return -1; 
	    }

	    // for ( int i=0; i<n; i++ )
	    // {
	    // 	for ( int j=0; j<m; j++ )
	    // 	{
	    // 	    fout << distance_table[i*m + j] << "," ;
	    // 	}
	    // 	fout << "\n";
	    // }

	    std::cout << "\n output pixels \n";
	    for ( int i=0; i<n; i++ )
	    {
		for ( int j=0; j<m; j++ )
		{
		    bool found_flag;
		    found_flag = false;
		    fout2 << distance_table[i*m+j] << ",";
		    for ( int k=0; k<node_list.size(); k++ )
		    {
			if ( node_list[k].irow == i && node_list[k].icol == j )
			{
			    fout << std::setw(2) << std::right << "-1" << ",";
			    found_flag = true;
			    break;
			}
		    }
		    if ( found_flag == false )
		    {
			fout << std::setw(2) << std::right << partition_table[i*m+j] << ",";
		    }
		}
		fout << "\n";
		fout2 << "\n";
	    }

	    fout.close();
	    fout2.close();
	    return 0;
	}

    int output_node_list()
	{
	    std::ofstream fnode("out/nodes.txt");
	    if(!fnode) 
	    {
		std::cout << "file open error.\n";
		return -1; 
	    }

	    for ( int k=0; k<voronoi_node_list.size(); k++)
	    {
		std::cout << "Node " << k << "(" 
			  << voronoi_node_list[k].x << ","
			  << voronoi_node_list[k].y << ") colors: ";
		for ( std::unordered_set<int>::iterator it=voronoi_node_list[k].color_set.begin(); 
		      it!=voronoi_node_list[k].color_set.end(); it++)
		{
		    std::cout << std::setw(2) << std::right << *it;
		}
		std::cout << "\n";		

		fnode << k << ", " <<
		    voronoi_node_list[k].x << "," << 
		    voronoi_node_list[k].y << ",";
		for ( std::unordered_set<int>::iterator it=voronoi_node_list[k].color_set.begin(); 
		      it!=voronoi_node_list[k].color_set.end(); it++)
		{
		    fnode << *it << ",";
		}
		fnode << "\n";		
	    }

	    fnode.close();
	    return 0;
	}

    int output_node(VoronoiNode &node)
	{
	    std::cout << "output node: "
		      << node.x << ","
		      << node.y << ") colors: ";
	    for ( std::unordered_set<int>::iterator it=node.color_set.begin(); 
		  it!=node.color_set.end(); it++)
	    {
		std::cout << *it << " ";
	    }
	    std::cout << "\n";		
	    return 0;
	}

    int output_cell_list(std::string file_index)
	{
	    std::ofstream fnode("out/cell_nodes" + file_index + ".txt");
	    std::ofstream fedge("out/cell_edges" + file_index + ".txt");
	    std::ofstream farea("out/cell_area" + file_index + ".txt");
	    if(!fnode || !fedge || !farea) 
	    {
		std::cout << "file open error.\n";
		return -1; 
	    }

	    for ( int k=0; k<voronoi_cell_list.size(); k++)
	    {
		// std::cout << "Cell " << voronoi_cell_list[k].cell_id << "\n node_set: ";
		// for ( int kk=0; kk<voronoi_cell_list[k].node_set.size(); kk++)
		// {
		//     std::cout << voronoi_cell_list[k].node_set[kk] << "(" <<
		// 	voronoi_node_list[voronoi_cell_list[k].node_set[kk]].x << ", " <<
		// 	voronoi_node_list[voronoi_cell_list[k].node_set[kk]].y << ") ";
		// }
		// std::cout << "\n edge_set: \n";		
		// for ( int kk=0; kk<voronoi_cell_list[k].edge_set.size(); kk++)
		// {
		//     std::cout << voronoi_edge_list[voronoi_cell_list[k].edge_set[kk]].edge_id 
		// 	      << "( node1: " << voronoi_edge_list[voronoi_cell_list[k].edge_set[kk]].node_id1
		// 	      << "(" << voronoi_edge_list[voronoi_cell_list[k].edge_set[kk]].x1 
		// 	      << ", " << voronoi_edge_list[voronoi_cell_list[k].edge_set[kk]].y1
		// 	      << ") - node2: " << voronoi_edge_list[voronoi_cell_list[k].edge_set[kk]].node_id2
		// 	      << "(" << voronoi_edge_list[voronoi_cell_list[k].edge_set[kk]].x2
		// 	      << ", " << voronoi_edge_list[voronoi_cell_list[k].edge_set[kk]].y2
		// 	      << ") )\n";
		// }
		// std::cout << "\n area: " << voronoi_cell_list[k].area;
		// std::cout << "\n neighbor_id: ";		
		// for ( int kk=0; kk<voronoi_cell_list[k].neighbor_id.size(); kk++)
		// {
		//     std::cout << voronoi_cell_list[k].neighbor_id[kk] << " ";
		// }
		// std::cout << "\n";		

    		VoronoiCell cell = voronoi_cell_list[k];
		std::vector<int> node_set = cell.node_set;
		std::vector<int> edge_set = cell.edge_set;
		
		for ( int kk=0; kk<node_set.size(); kk++)
		{
		    fnode << node_set[kk] << ", " <<
			voronoi_node_list[node_set[kk]].x << "," << 
			voronoi_node_list[node_set[kk]].y << ",";
		    fnode << "\n";		
		}

		for ( int kk=0; kk<edge_set.size(); kk++)
		{
		    fedge << voronoi_edge_list[edge_set[kk]].edge_id << ", " <<  
			voronoi_edge_list[edge_set[kk]].node_id1 << ", " <<
			voronoi_edge_list[edge_set[kk]].node_id2 << ", " <<
			voronoi_edge_list[edge_set[kk]].x1 << ", " <<
			voronoi_edge_list[edge_set[kk]].y1 << ", " <<
			voronoi_edge_list[edge_set[kk]].x2 << ", " <<
			voronoi_edge_list[edge_set[kk]].y2 << "\n";
		}
		
		farea << cell.cell_id << " ";  
		for ( int kk=0; kk<edge_set.size(); kk++)
		{
		    farea << cell.edges[kk].tri_sector_area/cell.area 
			  << " " << cell.edges[kk].ell_sector_area/cell.ellipse.area << " ";
		}
		farea << "\n";

	    }
	    fnode.close();
	    fedge.close();
	    farea.close();
	    return 0;
	}

    int output_cell(VoronoiCell &cell)
	{
	    std::cout << "Cell " << cell.cell_id
		      << " " << cell.c1 << " " << cell.c2 << "\n node_set: ";
	    for ( int kk=0; kk<cell.node_set.size(); kk++)
	    {
		std::cout << cell.node_set[kk] << "(" <<
		    voronoi_node_list[cell.node_set[kk]].x << ", " <<
		    voronoi_node_list[cell.node_set[kk]].y << ") ";
	    }
	    std::cout << "\n edge_set: \n";		
	    for ( int kk=0; kk<cell.edge_set.size(); kk++)
	    {
		std::cout << voronoi_edge_list[cell.edge_set[kk]].edge_id 
			  << "( node1: " << voronoi_edge_list[cell.edge_set[kk]].node_id1
			  << "(" << voronoi_edge_list[cell.edge_set[kk]].x1 
			  << ", " << voronoi_edge_list[cell.edge_set[kk]].y1
			  << ") - node2: " << voronoi_edge_list[cell.edge_set[kk]].node_id2
			  << "(" << voronoi_edge_list[cell.edge_set[kk]].x2
			  << ", " << voronoi_edge_list[cell.edge_set[kk]].y2
			  << ") )\n";
	    }
	    std::cout << "\n area: " << cell.area;
	    std::cout << "\n neighbor_id: ";		
	    for ( int kk=0; kk<cell.neighbor_id.size(); kk++)
	    {
		std::cout << cell.neighbor_id[kk] << " ";
	    }
	    std::cout << "\n";		
	    
	    return 0;
	}

    int output_edge_list()
	{
	    std::ofstream fedge("out/edges.txt");
	    if(!fedge) 
	    {
		std::cout << "file open error.\n";
		return -1; 
	    }
	    
	    std::cout << "output edges: " << voronoi_edge_list.size(); 
	    for ( int k=0; k<voronoi_edge_list.size(); k++)
	    {
		if ( voronoi_edge_list[k].side_colors.size() > 0 )
		{
		    std::cout << "Edge (" << k << "): " << voronoi_edge_list[k].edge_id << "(" 
			      << voronoi_edge_list[k].node_id1 << ","
			      << voronoi_edge_list[k].node_id2 << ") \n";

		    fedge <<  voronoi_edge_list[k].edge_id << ", " <<  
			voronoi_edge_list[k].node_id1 << ", " <<
			voronoi_edge_list[k].node_id2 << ", " <<
			voronoi_edge_list[k].x1 << ", " <<
			voronoi_edge_list[k].y1 << ", " <<
			voronoi_edge_list[k].x2 << ", " <<
			voronoi_edge_list[k].y2 << "\n";
		}
	    }

	    fedge.close();
	    return 0;
	}

    int distance_mapping(//BoundaryData boundary_data, 
	int* partition_table,
	int* distance_table)
	{
	    int is_done = 0;
	    int filling_value = 0;
	    int has_unfilled_pixel;
	    int min_nearby_dist;
	    int partition;

	    int n = boundary_data.n;
	    int m = boundary_data.m;
	    int nmin = boundary_data.nmin;
	    int nmax = boundary_data.nmax;
	    int mmin = boundary_data.mmin;
	    int mmax = boundary_data.mmax;
	
	    while ( is_done==0 )
	    {
		filling_value = filling_value + 1;

// output
/* std::cout << "\n output distance_table \n"; */
/* for ( int i=nmin; i<nmax; i++ ) */
/* { */
/* 	for ( int j=mmin; j<mmax; j++ ) */
/* 	{ */
/* 	    std::cout << std::setw(2) << std::right << distance_table[i*m+j] << " "; */
/* 	} */
/* 	std::cout << "\n"; */
/* } */
/* std::cout << "\n output partition_table \n"; */
/* for ( int i=nmin; i<nmax; i++ ) */
/* { */
/* 	for ( int j=mmin; j<mmax; j++ ) */
/* 	{ */
/* 	    std::cout << std::setw(2) << std::right << partition_table[i*m+j] << " "; */
/* 	} */
/* 	std::cout << "\n"; */
/* } */

		has_unfilled_pixel = 0;
		for ( int i=nmin; i<nmax; i++ )
		{
		    for ( int j=mmin; j<mmax; j++ )
		    {
			if ( distance_table[i*m+j] > filling_value )
			{
			    std::cout << "error filling values\n";
			    getchar();
			}
			if ( distance_table[i*m+j] == filling_value )
			{
			    min_nearby_dist = 100000;
			    partition = -1;
			    if ( distance_table[(i+1)*m+j] < min_nearby_dist )
			    {
				min_nearby_dist = distance_table[(i+1)*m+j];
				partition = partition_table[(i+1)*m+j];
			    }
			    if ( distance_table[(i-1)*m+j] < min_nearby_dist )
			    {
				min_nearby_dist = distance_table[(i-1)*m+j];
				partition = partition_table[(i-1)*m+j];
			    }
			    if ( distance_table[i*m+j+1] < min_nearby_dist )
			    {
				min_nearby_dist = distance_table[i*m+j+1];
				partition = partition_table[i*m+j+1];
			    }
			    if ( distance_table[i*m+j-1] < min_nearby_dist )
			    {
				min_nearby_dist = distance_table[i*m+j-1];
				partition = partition_table[i*m+j-1];
			    }
			    if ( min_nearby_dist == filling_value)
			    {
				distance_table[i*m+j] = filling_value + 1;
//std::cout << "update " << i << " " << j << " " << partition_table[i*m+j] << "\n";
			    }
			    else if ( min_nearby_dist == filling_value - 1)
			    {
				partition_table[i*m+j] = partition;
				if ( has_unfilled_pixel==0 )
				{
				    has_unfilled_pixel = 1;
				}
			    }
			}
		    }
		}
		if ( has_unfilled_pixel==0 )
		{
		    is_done = 1;
		}
		//std::cout << "distance " << filling_value << "\n";
	    }
	    return 0;
	}

    int init_distance_table(std::vector<Ellipse> &ellipse_list,
			    //BoundaryData boundary_data, 
			    int* partition_table,
			    int* distance_table)
	{
	    int n = boundary_data.n;
	    int m = boundary_data.m;
	    int nmin = boundary_data.nmin;
	    int nmax = boundary_data.nmax;
	    int mmin = boundary_data.mmin;
	    int mmax = boundary_data.mmax;

	    double xmin = boundary_data.xmin;
	    double xmax = boundary_data.xmax;
	    double ymin = boundary_data.ymin;
	    double ymax = boundary_data.ymax;

	    double dx = boundary_data.res_x;
	    double dy = boundary_data.res_y;

	    if ( dx<=0 || dy<=0 )
	    {
		std::cout << "make_distance_table: parameter wrong \n";
		getchar();
	    }
	
	    double pixel_x;
	    double pixel_y;
	    for ( int i=0; i<n; i++ )
	    {
		pixel_y = ymin + (i+0.5)*dy;
		for ( int j=0; j<m; j++ )
		{
		    pixel_x = xmin + (j+0.5)*dx;
		    distance_table[i*m+j] = -1;
		    partition_table[i*m+j] = -1;
		    for ( int k=0; k<ellipse_list.size(); k++ )
		    {
			if ( ellipse_list[k].isInside(pixel_x, pixel_y) )
			{
//std::cout << "pixel " << i << ", "<< j << "inside ellipse" << k << "\n";
			    if (distance_table[i*m+j] == -1)
			    {
				distance_table[i*m+j] = 0;
				partition_table[i*m+j] = ellipse_list[k].ellipse_id;
			    }
			    else
			    {
				distance_table[i*m+j] = 1;
				partition_table[i*m+j] = -1;
			    }
			}
		    }
		    if (distance_table[i*m+j] == -1)
		    {
			distance_table[i*m+j] = 1;
		    }
		}
	    }
	    return 0;
	}

    template <typename T>
    std::vector<size_t> sort_indexes(const std::vector<T> &v) 
    	{

    	    // initialize original index locations
    	    std::vector<size_t> idx(v.size());
    	    iota(idx.begin(), idx.end(), 0);
	    
    	    // sort indexes based on comparing values in v
	    std::sort(idx.begin(), idx.end(),
		      [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

    	    return idx;
    	}

    double get_tri_sector_area(const double c1, const double c2,
			       const double x1, const double y1, 
			       const double x2, const double y2)
	{
	    double tri_area;
	    tri_area = 0.5*fabs(c1*y1 - c2*x1 + x1*y2 - y1*x2 + x2*c2 - c1*y2);
	    return tri_area;
	}

    double get_ell_sector_area(const Ellipse &ellipse, 
			       const double theta1,
			       const double theta2)
	{
	    double ell_area;
	    double new_theta1 = atan2(ellipse.a*sin(theta1), ellipse.b*cos(theta1));
	    double new_theta2 = atan2(ellipse.a*sin(theta2), ellipse.b*cos(theta2));
	    if ( theta2 > PI )
	    {
		new_theta2 = new_theta2 + 2*PI;
	    }
	    //std::cout << "theta1 " << theta1 << "new theta1 " << new_theta1 << "\n";
	    //std::cout << "theta2 " << theta2 << "new theta2 " << new_theta2 << "\n";
	    ell_area = 0.5*ellipse.a*ellipse.b*fabs(new_theta1 - new_theta2);
	    return ell_area;
	}

};


#endif //GEOMETRY
