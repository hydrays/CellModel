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
	    int margin_n = 40;
	    int margin_m = 40;
	    int inner_n = 200;
	    int inner_m = 200;
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

struct VoronoiCell
{
    bool is_valid;
    int cell_id;
    std::vector<int> node_set;
    std::vector<int> edge_set;
    double area;
    std::vector<int> neighbor_id;
};

struct VoronoiNode
{
    int node_id;
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
    std::vector<int> side_colors;
    VoronoiEdge(int edge_id, int id1, int id2)
	: edge_id(edge_id), node_id1(id1), node_id2(id2)
	{
	    n1 = 0.0;
	    n2 = 0.0;
	    length = -1.0;
	}
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
	    boundary_data.set_boundary_data(15.0, 10.0);
	    //boundary_data.set_boundary_data(36.0, 36.0);

	    init_ellipse_list();
	    std::cout << "Status: ellipse list initialized. \n";

	    apply_boundary_condition_ellipse();
	    ellipse_to_voronoi();
	    update_map();

	    voronoi_to_ellipse();
	    output_new_ellipse_list();
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

    int init_ellipse_list()
	{
	    std::ifstream finit("initPos.txt");
	    if(!finit) 
	    {
		std::cout << "file open error.\n";
		return -1; 
	    }
	    
	    int id = 0;
	    double cx, cy;
	    while (finit >> cx >>  cy)
	    {
		id = id + 1;
		double v1 = 0.5;
		double v2 = 0.0;
		double u1 = 0.0;
		double u2 = 0.5;
		double a = sqrt(v1*v1 + v2*v2);
		double b = sqrt(u1*u1 + u2*u2);
		double area = PI*a*b;
		std::cout << " > ellipse: " << cx << " " << cy << " " << ellipse_list.size() << "\n";
		Ellipse ellipse(id, cx, cy, v1, v2, u1, u2);
		ellipse.a = a;
		ellipse.b = b;
		ellipse.area = area;
		ellipse_list.push_back(ellipse);
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
	    distance_mapping(partition_table, distance_table);

	    //find_nodes(partition_table, node_list, raw_node_list, node_color_set);

	    make_node_list(node_list, node_color_set, partition_table, distance_table);
	    make_cell_list(ellipse_list);
	    make_edge_list(); // must after "make_voronoi_cell_list"
	    refine_cell_list();

	    output_ellipse_list();
	    output_node_list();
	    output_edge_list();
	    output_cell_list();
	    output_partition_table(n, m, partition_table, node_list);

	    delete(distance_table);
	    delete(partition_table);
	    return 0;
	}

    int voronoi_to_ellipse()
	{
	    //std::vector<Ellipse> new_ellipse_list;
	    new_ellipse_list.clear();
	    for ( int k=0; k<voronoi_cell_list.size(); k++ )	    
	    {
		Ellipse new_ellipse = ellipse_list[k];
		if ( voronoi_cell_list[k].node_set.size() > 3 )
		{
		    double I1 = 0.0;
		    double I2 = 0.0;
		    double I3 = 0.0;
		    double c1 = 0.0;
		    double c2 = 0.0;
		    for ( int kk=0; kk<voronoi_cell_list[k].node_set.size(); kk++ )
		    {
			c1 = c1 + voronoi_node_list[voronoi_cell_list[k].node_set[kk]].x;
			c2 = c2 + voronoi_node_list[voronoi_cell_list[k].node_set[kk]].y;
		    }
		    c1 = c1/voronoi_cell_list[k].node_set.size();
		    c2 = c2/voronoi_cell_list[k].node_set.size();

		    for ( int kk=1; kk<voronoi_cell_list[k].node_set.size(); kk++ )
		    {
			double x1 = voronoi_node_list[voronoi_cell_list[k].node_set[kk-1]].x - c1;
			double y1 = voronoi_node_list[voronoi_cell_list[k].node_set[kk-1]].y - c2;
			double x2 = voronoi_node_list[voronoi_cell_list[k].node_set[kk]].x - c1;
			double y2 = voronoi_node_list[voronoi_cell_list[k].node_set[kk]].y - c2;
			I1 = I1 + (-(x2-x1)/12.0)*(y2*y2*y2 + y2*y2*y1 + y2*y1*y1 + y1*y1*y1);
			I2 = I2 + ((x2-x1)/24.0)*(x1*(3.0*y1*y1 + 2.0*y1*y2 + y2*y2) + 
						  x2*(3.0*y2*y2 + 2.0*y1*y2 + y1*y1));
			I3 = I3 + ((y2-y1)/12.0)*(x2*x2*x2 + x2*x2*x1 + x2*x1*x1 + x1*x1*x1); 
			std::cout << x1 << ", " << y1 << ", ";
		    }		
		    int kk = voronoi_cell_list[k].node_set.size()-1;
		    double x1 = voronoi_node_list[voronoi_cell_list[k].node_set[kk]].x - c1;
		    double y1 = voronoi_node_list[voronoi_cell_list[k].node_set[kk]].y - c2;
		    double x2 = voronoi_node_list[voronoi_cell_list[k].node_set[0]].x - c1;
		    double y2 = voronoi_node_list[voronoi_cell_list[k].node_set[0]].y - c2;
		    I1 = I1 + (-(x2-x1)/12.0)*(y2*y2*y2 + y2*y2*y1 + y2*y1*y1 + y1*y1*y1);
		    I2 = I2 + ((x2-x1)/24.0)*(x1*(3.0*y1*y1 + 2.0*y1*y2 + y2*y2) + 
					      x2*(3.0*y2*y2 + 2.0*y1*y2 + y1*y1));
		    I3 = I3 + ((y2-y1)/12.0)*(x2*x2*x2 + x2*x2*x1 + x2*x1*x1 + x1*x1*x1); 
		    std::cout << x1 << ", " << y1 << "\n";
		    double lambda1 = (I1 + I3 + sqrt(I1*I1 + I3*I3 + 4.0*I2*I2 - 2.0*I1*I3))/2.0;
		    double lambda2 = (I1 + I3 - sqrt(I1*I1 + I3*I3 + 4.0*I2*I2 - 2.0*I1*I3))/2.0;
		    double v1 = lambda1 - I3;
		    double v2 = I2;
		    double u1 = lambda2 - I3;
		    double u2 = I2;
		    new_ellipse.c1 = c1;
		    new_ellipse.c2 = c2;
		    new_ellipse.v1 = v1;
		    new_ellipse.v2 = v2;
		    new_ellipse.u1 = u1;
		    new_ellipse.u2 = u2;
		    new_ellipse.a = lambda1;
		    new_ellipse.b = lambda2;

		    std::cout << "ellipse " << k << ":[" << I1 << " " <<
			I2 << " " << I3 << "]\n";
		}
		else
		{
		    new_ellipse.v1 = 0.0;
		    new_ellipse.v2 = 0.0;
		    new_ellipse.u1 = 0.0;
		    new_ellipse.u2 = 0.0;
		    new_ellipse.a = 0.0;
		    new_ellipse.b = 0.0;
		}
		new_ellipse_list.push_back(new_ellipse);
	    }
	    return 0;
	}

    int refine_cell_list()
	{
	    //for ( int k=0; k<voronoi_cell_list.size(); k++ )	    
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
			//std::cout << "erase " << k << "-" << kk << "\n";
			//voronoi_cell_list[k].edge_set.erase(voronoi_cell_list[k].edge_set.begin() + kk);
		    }
		}
		if ( edges.size() != refined_edges.size() )
		{
		    it->edge_set.clear();
		    it->edge_set = refined_edges;
		}
		
		if ( refined_edges.size() < 3 )
		{
		    voronoi_cell_list.erase(it);
		    //voronoi_cell_list[k].area = -1.0;
		}
		else
		{
		    it->is_valid = true;
		    double area = 0.0;
		    for ( int kk=1; kk<it->node_set.size(); kk++ )
		    {
			int node_id1 = it->node_set[kk-1];
			int node_id2 = it->node_set[kk];
			double x1 = voronoi_node_list[node_id1].x;
			double y1 = voronoi_node_list[node_id1].y;
			double x2 = voronoi_node_list[node_id2].x;
			double y2 = voronoi_node_list[node_id2].y;
			area = area + (x1*y2 - x2*y1);
		    }
		    int node_id1 = it->node_set[it->node_set.size()-1];
		    int node_id2 = it->node_set[0];
		    double x1 = voronoi_node_list[node_id1].x;
		    double y1 = voronoi_node_list[node_id1].y;
		    double x2 = voronoi_node_list[node_id2].x;
		    double y2 = voronoi_node_list[node_id2].y;
		    area = area + (x1*y2 - x2*y1);
		    it->area = 0.5*fabs(area);
		    ++it;
		}
	    }
	    return 0;
	}

    int make_node_list(const std::vector<PixelNode> &node_list, 
		       std::unordered_set<std::unordered_set<int>, ColorHash> &node_color_set,
		       const int* partition_table,
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
		    if ( (node_color.size() == 4) &&
			 (node_color_set.find(node_color) == node_color_set.end()) )
		    {
			std::cout << "new node_color: size " << node_color.size() << " colors: ";
			for (std::unordered_set<int>::const_iterator it = node_color.begin();
			     it != node_color.end(); it++ )
			{
			    std::cout << *it << " ";
			}
			std::cout << "\n";

			node_color_set.insert(node_color);

			VoronoiNode voronoi_node;
			int node_id = i*m + j;
			//id = id + 1;
			//i = node_list[k].irow;
			//j = node_list[k].icol;

			voronoi_node.node_id = node_id;
			voronoi_node.color_set = node_color;
			voronoi_node.x = xmin + (j+0.5)*dx;
			voronoi_node.y = ymin + (i+0.5)*dy;
		
			voronoi_node_list.push_back(voronoi_node);
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
			 isLocalMaximum(i, j, n, m, distance_table) &&
			 (node_color_set.find(node_color) == node_color_set.end()) &&
			 !isSubset(node_color, node_color_set) )
		    {
			std::cout << "new node_color: ";
			for (std::unordered_set<int>::const_iterator it = node_color.begin();
			     it != node_color.end(); it++ )
			{
			    std::cout << *it << " ";
			}
			std::cout << "\n";

			node_color_set.insert(node_color);

			VoronoiNode voronoi_node;
			int node_id = i*m + j;
			//id = id + 1;
			//i = node_list[k].irow;
			//j = node_list[k].icol;

			voronoi_node.node_id = node_id;
			voronoi_node.color_set = node_color;
			voronoi_node.x = xmin + (j+0.5)*dx;
			voronoi_node.y = ymin + (i+0.5)*dy;
		
			voronoi_node_list.push_back(voronoi_node);
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
		voronoi_cell.cell_id = cell_id;

		//std::cout << "make_voronoi_cell_list: cell_id" << cell_id << ":";
		 
		for ( int kk=0; kk<voronoi_node_list.size(); kk++ )
		{
		    std::unordered_set<int>::iterator got = 
			voronoi_node_list[kk].color_set.find(cell_id);
		    if ( got != voronoi_node_list[kk].color_set.end() )
		    {	
			//std::cout << "find node\n"; 
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
		    getchar();
		}
		edge.side_colors = side_colors;
		voronoi_edge_list.push_back(edge);
		//std::cout << it->node_id1 << " " << it->node_id2 << "\n";
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
	    int start_node_id;
	    int end_node_id;
	    for ( int k=1; k<node_set.size(); k++ )
	    {
		int node_id1 = node_set[k];
		int node_id2 = node_set[k-1];
		if ( node_id1 > node_id2 )
		{
		    start_node_id = node_id2;
		    end_node_id = node_id1;
		}
		else
		{
		    end_node_id = node_id2;
		    start_node_id = node_id1;
		}
		VoronoiEdge edge(primary_edge_set.size(), start_node_id, end_node_id);
		primary_edge_set.insert(edge);
		std::unordered_set<VoronoiEdge, Hash>::iterator got = primary_edge_set.find(edge);
		edge_set.push_back(got->edge_id);
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
    
    int output_ellipse_list()
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
		    ellipse_list[k].c2 << "\n";
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
			       const int* partition_table, 
			       const std::vector<PixelNode> &node_list)
	{
	    std::ofstream fout("out/partition_table.txt");
	    if(!fout) 
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
	    }

	    fout.close();
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
		std::cout << "Node " << voronoi_node_list[k].node_id << "(" 
			  << voronoi_node_list[k].x << ","
			  << voronoi_node_list[k].y << ") colors: ";
		for ( std::unordered_set<int>::iterator it=voronoi_node_list[k].color_set.begin(); 
		      it!=voronoi_node_list[k].color_set.end(); it++)
		{
		    std::cout << std::setw(2) << std::right << *it;
		}
		std::cout << "\n";		

		fnode << voronoi_node_list[k].node_id << ", " <<
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

    int output_cell_list()
	{
	    // output nodes
	    std::ofstream fnode("out/cell_nodes.txt");
	    std::ofstream fedge("out/cell_edges.txt");
	    if(!fnode || !fedge) 
	    {
		std::cout << "file open error.\n";
		return -1; 
	    }

	    // output edges
	    for ( int k=0; k<voronoi_cell_list.size(); k++)
	    {
		std::cout << "Cell " << voronoi_cell_list[k].cell_id << "\n node_set: ";
		for ( int kk=0; kk<voronoi_cell_list[k].node_set.size(); kk++)
		{
		    std::cout << voronoi_cell_list[k].node_set[kk] << "(" <<
			voronoi_node_list[voronoi_cell_list[k].node_set[kk]].x << ", " <<
			voronoi_node_list[voronoi_cell_list[k].node_set[kk]].y << ") ";
		}
		std::cout << "\n edge_set: ";		
		for ( int kk=0; kk<voronoi_cell_list[k].edge_set.size(); kk++)
		{
		    std::cout << voronoi_cell_list[k].edge_set[kk] << " ";
		}
		std::cout << "\n area: " << voronoi_cell_list[k].area;
		std::cout << "\n neighbor_id: ";		
		for ( int kk=0; kk<voronoi_cell_list[k].neighbor_id.size(); kk++)
		{
		    std::cout << voronoi_cell_list[k].neighbor_id[kk] << " ";
		}
		std::cout << "\n";		
		std::vector<int> node_set = voronoi_cell_list[k].node_set;
		std::vector<int> edge_set = voronoi_cell_list[k].edge_set;

		for ( int kk=0; kk<node_set.size(); kk++)
		{
		    fnode << voronoi_node_list[node_set[kk]].node_id << ", " <<
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
	    }
	    fnode.close();
	    fedge.close();
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
	    
	    for ( int k=0; k<voronoi_edge_list.size(); k++)
	    {
		if ( voronoi_edge_list[k].side_colors.size() > 0 )
		{
		    std::cout << "Edge " << voronoi_edge_list[k].edge_id << "(" 
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
		std::cout << "distance " << filling_value << "\n";
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
};

#endif //GEOMETRY
