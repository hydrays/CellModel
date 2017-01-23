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
#include "ellipse.h"
#include "boundary_data.h"

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
	    res_x = L/inner_m;
	    res_y = H/inner_n;
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
    double n1;
    double n2;
    double length;
    int side_color[2];
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

    int init()
	{
	    boundary_data.set_boundary_data(15.0, 10.0);
	    //boundary_data.set_boundary_data(36.0, 36.0);

	    init_ellipse_list();
	    std::cout << "Status: ellipse list initialized. \n";

	    apply_boundary_condition_ellipse();
	    ellipse_to_voronoi();
	}

    int apply_boundary_condition_ellipse()
	{
	    double L = boundary_data.L;
	    double H = boundary_data.H;
	    double margin_L = boundary_data.margin_L;
	    double margin_H = boundary_data.margin_H;

	    for ( int k=0; k<ellipse_list.size(); k++ )
	    {
		if ( (ellipse_list[k].c1 < 0.0) || (ellipse_list[k].c1 > L) 
		    || (ellipse_list[k].c2 < 0.0) || (ellipse_list[k].c2 > H) )
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
		     && (ellipse_list[k].c2 <= H) && (ellipse_list[k].c2 >= 0) )
		{
		    Ellipse ellipse = ellipse_list[k];
		    ellipse.c1 = ellipse.c1 - L;
		    ellipse.ellipse_id = ellipse.ellipse_id + MAX_NUM_ELLIPSE;
		    ellipse_list.push_back(ellipse);		    
		}
		// right
		if ( (ellipse_list[k].c1 <= margin_L) && (ellipse_list[k].c1 >= 0)
		     && (ellipse_list[k].c2 <= H) && (ellipse_list[k].c2 >= 0) )
		{
		    Ellipse ellipse = ellipse_list[k];
		    ellipse.c1 = ellipse.c1 + L;
		    ellipse.ellipse_id = ellipse.ellipse_id + 2*MAX_NUM_ELLIPSE;
		    ellipse_list.push_back(ellipse);		    
		}
		// top
		if ( //(ellipse_list[k].c1 <= L) && (ellipse_list[k].c1 >= 0) )
		     (ellipse_list[k].c2 <= margin_H) && (ellipse_list[k].c2 >= 0) )
		{
		    Ellipse ellipse = ellipse_list[k];
		    ellipse.c2 = ellipse.c2 + H;
		    ellipse.ellipse_id = ellipse.ellipse_id + 3*MAX_NUM_ELLIPSE;
		    ellipse_list.push_back(ellipse);		    
		}
		// down
		if ( //(ellipse_list[k].c1 <= L) && (ellipse_list[k].c1 >= 0) )
		     (ellipse_list[k].c2 <= H) && (ellipse_list[k].c2 >= H - margin_H) )
		{
		    Ellipse ellipse = ellipse_list[k];
		    ellipse.c2 = ellipse.c2 - H;
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
		std::cout << " > ellipse: " << cx << " " << cy << " " << ellipse_list.size() << "\n";
		Ellipse ellipse(id, cx, cy, v1, v2, u1, u2);
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

	    make_voronoi_node_list(node_list, node_color_set, partition_table, distance_table);
	    make_voronoi_cell_list(ellipse_list);
	    make_voronoi_edge_list(); // must after "make_voronoi_cell_list"

	    output_ellipse_list();
	    output_voronoi_node_list();
	    output_voronoi_edge_list();
	    output_voronoi_cell_list();
	    output_partition_table(n, m, partition_table, node_list);

	    delete(distance_table);
	    delete(partition_table);
	    return 0;
	}

    int voronoi_to_ellipse()
	{
	    return 0;
	}

    int make_voronoi_node_list(const std::vector<PixelNode> &node_list, 
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

	    double xmin = -boundary_data.margin_L;
	    double xmax = boundary_data.L + boundary_data.margin_L;
	    double ymin = -boundary_data.margin_H;
	    double ymax = boundary_data.H + boundary_data.margin_H;

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
			voronoi_node.x = ymin + (j+0.5)*dx;
			voronoi_node.y = xmin + (i+0.5)*dy;
		
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
			voronoi_node.x = ymin + (j+0.5)*dx;
			voronoi_node.y = xmin + (i+0.5)*dy;
		
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

    int make_voronoi_cell_list(std::vector<Ellipse> &ellipse_list)
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
		std::vector<int> sorted_node_set = voronoi_arrange_node(ellipse_list[k].c1, ellipse_list[k].c2, node_set);
		voronoi_cell.node_set = sorted_node_set;

		edge_set = prepare_edges(sorted_node_set);

		voronoi_cell.edge_set = edge_set;

		voronoi_cell_list.push_back(voronoi_cell);
	    }
	    return 0;
	}

    int make_voronoi_edge_list()
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

		std::unordered_set<int> color_set1 = node1.color_set;
		std::unordered_set<int> color_set2 = node2.color_set;

		voronoi_edge_list.push_back(edge);
		//std::cout << it->node_id1 << " " << it->node_id2 << "\n";
	    }
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

    int output_partition_table(const int n, const int m, 
			       const int* partition_table, 
			       const std::vector<PixelNode> node_list)
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

    int output_voronoi_node_list()
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

    int output_voronoi_cell_list()
	{
	    for ( int k=0; k<voronoi_cell_list.size(); k++)
	    {
		std::cout << "Cell " << voronoi_cell_list[k].cell_id << " node_set: ";
		for ( int kk=0; kk<voronoi_cell_list[k].node_set.size(); kk++)
		{
		    std::cout << std::setw(4) << std::right << voronoi_cell_list[k].node_set[kk];
		}
		std::cout << "\n";		
		for ( int kk=0; kk<voronoi_cell_list[k].edge_set.size(); kk++)
		{
		    std::cout << std::setw(4) << std::right << voronoi_cell_list[k].edge_set[kk];
		}
		std::cout << "\n";		
	    }
	    return 0;
	}

    int output_voronoi_edge_list()
	{
	    std::ofstream fedge("out/edges.txt");
	    if(!fedge) 
	    {
		std::cout << "file open error.\n";
		return -1; 
	    }
	    
	    for ( int k=0; k<voronoi_edge_list.size(); k++)
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

	    double xmin = -boundary_data.margin_L;
	    double xmax = boundary_data.L + boundary_data.margin_L;
	    double ymin = -boundary_data.margin_H;
	    double ymax = boundary_data.H + boundary_data.margin_H;

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
		pixel_y = xmin + (i+0.5)*dy;
		for ( int j=0; j<m; j++ )
		{
		    pixel_x = ymin + (j+0.5)*dx;
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

//     int find_nodes(int* partition_table, std::vector<PixelNode> &node_list, 
// 		   std::vector<PixelNode> &raw_node_list,
// 		   unordered_set<unordered_set<int>, ColorHash> &node_color_set)
// 	{
// 	    int n = boundary_data.n;
// 	    int m = boundary_data.m;
// 	    int nmin = boundary_data.nmin;
// 	    int nmax = boundary_data.nmax;
// 	    int mmin = boundary_data.mmin;
// 	    int mmax = boundary_data.mmax;

// 	    int nearby_values[8];
// 	    int num_colors;
// 	    int color_value;
// 	    int first_value;
// 	    int node_id;
// 	    for ( int i=nmin; i<nmax; i++ )
// 	    {
// 		for ( int j=mmin; j<mmax; j++ )
// 		{
// 		    nearby_values[0] = partition_table[(i-1)*m+j];
// 		    nearby_values[1] = partition_table[(i-1)*m+j-1];
// 		    nearby_values[2] = partition_table[(i)*m+j-1];
// 		    nearby_values[3] = partition_table[(i+1)*m+j-1];
// 		    nearby_values[4] = partition_table[(i+1)*m+j];
// 		    nearby_values[5] = partition_table[(i+1)*m+j+1];
// 		    nearby_values[6] = partition_table[(i)*m+j+1];
// 		    nearby_values[7] = partition_table[(i-1)*m+j+1];

// 		    const size_t len = sizeof(nearby_values) / sizeof(nearby_values[0]);
// 		    std::unordered_set<int> node_color(nearby_values, nearby_values+len);
// 		    node_color.erase( -1 );
// 		    if ( (node_color.size() >= 3) && 
// 			 (node_color_set.find(node_color) == node_color_set.end()) )
// 		    {
// 			cout << "new node_color: ";
// 			for (std::unordered_set<int>::const_iterator it = node_color.begin();
// 			     it != node_color.end(); it++ )
// 			{
// 			    cout << *it << " ";
// 			}
// 			cout << "\n";

// 			node_color_set.insert(node_color);
// 			node_id = i*m + j;
// 			node_list.push_back(PixelNode(i,j, node_color.size(), node_id));
// 		    }

// 		    // std::unordered_set<int> set1;
// 		    // set1.insert(3);
// 		    // set1.insert(6);
// 		    // std::unordered_set<int> set2;
// 		    // set2.insert(6);
// 		    // set2.insert(3);
// 		    // if ( set1 == set2 )
// 		    // {
// 		    // 	cout << "equal\n";
// 		    // }
// 		    // else
// 		    // {
// 		    // 	cout << "not equal!\n";			
// 		    // }
// 		}
// 	    }

// /* std::cout << "\n output pixels \n"; */
// /* for ( int i=0; i<n; i++ ) */
// /* { */
// /*     for ( int j=0; j<m; j++ ) */
// /*     { */
// /* 	bool found_flag; */
// /* 	found_flag = false; */
// /* 	for ( int k=0; k<raw_node_list.size(); k++ ) */
// /* 	{ */
// /* 	    if ( raw_node_list[k].irow == i && raw_node_list[k].icol == j ) */
// /* 	    { */
// /* 		std::cout << std::setw(2) << std::right << "*" << " "; */
// /* 		found_flag = true; */
// /* 		break; */
// /* 	    } */
// /* 	} */
// /* 	if ( found_flag == false ) */
// /* 	{ */
// /* 	    std::cout << std::setw(2) << std::right << partition_table[i*m+j] << " "; */
// /* 	} */
// /*     } */
// /*     std::cout << "\n"; */
// /* } */

// 	    // clean_node(n, m, raw_node_list, node_list);

// 	    // for ( int k=0; k<node_list.size(); k++ )
// 	    // {
// 	    // 	std::cout << node_list[k].irow << " " << node_list[k].icol << "\n";
// 	    // }

// /* std::cout << "\n output pixels size 2: " << node_list.size() << "\n"; */
// /* for ( int i=0; i<n; i++ ) */
// /* { */
// /*     for ( int j=0; j<m; j++ ) */
// /*     { */
// /* 	bool found_flag; */
// /* 	found_flag = false; */
// /* 	for ( int k=0; k<node_list.size(); k++ ) */
// /* 	{ */
// /* 	    if ( node_list[k].irow == i && node_list[k].icol == j ) */
// /* 	    { */
// /* 		std::cout << std::setw(2) << std::right << "*" << " "; */
// /* 		found_flag = true; */
// /* 		break; */
// /* 	    } */
// /* 	} */
// /* 	if ( found_flag == false ) */
// /* 	{ */
// /* 	    std::cout << std::setw(2) << std::right << partition_table[i*m+j] << " "; */
// /* 	} */
// /*     } */
// /*     std::cout << "\n"; */
// /* } */

// 	    return 0;
// 	}

    // int clean_node(int n, int m, 
    // 		   std::vector<PixelNode> &raw_node_list,
    // 		   std::vector<PixelNode> &node_list)
    // 	{
    // 	    std::unordered_set<int> raw_node_set;
    // 	    int node_id;

    // 	    for ( int k=0; k<raw_node_list.size(); k++ )
    // 	    {
    // 		raw_node_set.insert(raw_node_list[k].node_id);
    // 	    }

    // 	    if ( raw_node_list.size() != raw_node_set.size() )
    // 	    {
    // 		std::cout << "clearn_node: something wrong with node id assignment.\n";
    // 		getchar();
    // 	    }
	
    // 	    for ( int k=0; k<raw_node_list.size(); k++ )
    // 	    {
    // 		std::unordered_set<int>::iterator got = raw_node_set.find(raw_node_list[k].node_id);
    // 		if ( got != raw_node_set.end() )
    // 		{
    // 		    if ( raw_node_list[k].degree > 4 )
    // 		    {
    // 			std::cout << "find node: something wrong with node degree.\n";
    // 			getchar();
    // 		    }
    // 		    else if ( raw_node_list[k].degree == 4 )
    // 		    {
    // 			node_id = raw_node_list[k].node_id;
    // 			raw_node_set.erase(node_id-m);
    // 			raw_node_set.erase(node_id-m-1);
    // 			raw_node_set.erase(node_id-1);
    // 			raw_node_set.erase(node_id+m-1);
    // 			raw_node_set.erase(node_id+m);
    // 			raw_node_set.erase(node_id+m+1);
    // 			raw_node_set.erase(node_id+1);
    // 			raw_node_set.erase(node_id-m+1);
    // 		    }
    // 		}
    // 	    }
    // 	    std::cout << "size:" << raw_node_set.size() << "\n";

    // 	    for ( int k=0; k<raw_node_list.size(); k++ )
    // 	    {
    // 		std::unordered_set<int>::iterator got = raw_node_set.find(raw_node_list[k].node_id);
    // 		if ( got != raw_node_set.end() )
    // 		{
    // 		    if ( raw_node_list[k].degree == 3 )
    // 		    {
    // 			node_id = raw_node_list[k].node_id;
    // 			raw_node_set.erase(node_id-m);
    // 			raw_node_set.erase(node_id-m-1);
    // 			raw_node_set.erase(node_id-1);
    // 			raw_node_set.erase(node_id+m-1);
    // 			raw_node_set.erase(node_id+m);
    // 			raw_node_set.erase(node_id+m+1);
    // 			raw_node_set.erase(node_id+1);
    // 			raw_node_set.erase(node_id-m+1);
    // 		    }
    // 		}
    // 	    }
    // 	    std::cout << "size:" << raw_node_set.size() << "\n";
		
    // 	    node_list.clear();
    // 	    std::unordered_set<int>::iterator got;
    // 	    for ( int k=0; k<raw_node_list.size(); k++ )
    // 	    {
    // 		got = raw_node_set.find(raw_node_list[k].node_id);
    // 		if ( got != raw_node_set.end() )
    // 		{
    // 		    node_list.push_back(raw_node_list[k]);
    // 		}
    // 	    }
    // 	    std::cout << "size2:" << node_list.size() << "\n";
    // 	}

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
