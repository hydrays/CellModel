#ifndef ENGINE
#define ENGINE

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
#include <time.h>
//#include "SuiteSparseQR.hpp"

struct DivInfo
{
    double x;
    double y;
    int type;
    double ratio;
};

class Engine
{
    double eta = 1.0;
    double dt = 0.2;
    double q0 = 1.56;
    double scaling_beta1 = 0.5;
    double scaling_beta2 = 0.05;
    double A1 = 1.0;
    double A2 = 0.1;
    double F_factor = 0.0;
    double p0 = 0.63;
    int div_counter;
	
public:    
    int run(Geometry &geometry, Force &force)
	{
	    std::ofstream foutput("out/statistics.txt");
	    if(!foutput) 
	    {
		std::cout << "file open error.\n";
		return -1; 
	    }

	    srand (time(NULL));
	    div_counter = 0;
	    F_factor = 0.0;
	    for ( int istep=100; istep<1600; istep++ )
	    //int istep = 100;
	    {
		std::cout << "ellipse num: " << geometry.ellipse_list.size() << "\n";
		//for ( int jstep=0; jstep<10; jstep++ )
		//{
		geometry.output_cell_list(std::to_string(istep));
		geometry.output_ellipse_list(std::to_string(istep));
		force.update_force_field(geometry);
		force.output_force_filed(geometry, std::to_string(istep));
		//force.output_force_filed(geometry, std::to_string(istep));
		one_step_solution_modify(geometry, force, std::to_string(istep));
		//}

		if ( istep % 10 == 0)
		//if ( istep == 110)
		//if ( 0 )
		{
		    DivInfo div_info = cell_division(geometry, force);
		    div_counter = div_counter + 1;
		    if ( div_info.type==1 )
		    {
			F_factor = F_factor + p0;
		    }
		    else
		    {
			F_factor = F_factor - (1-p0);
		    }
			
		    foutput << div_counter << ", "
			    << geometry.boundary_data.L << ", "
			    << geometry.boundary_data.H << ", "
			    << div_info.x << ", "
			    << div_info.y << ", "
			    << div_info.type << ", "
			    << F_factor << ", "
			    << div_info.ratio << "\n";
		    foutput.flush();
		}
		
		geometry.apply_boundary_condition_voronoi();

		//getchar();
		geometry.update_ellipse_list();

		//geometry.output_ellipse_list(std::to_string(100));
		//geometry.write_ellipse_tofile();

		std::cout << "istep: " << istep << 
		    "ellipse num: " << geometry.ellipse_list.size() << "\n";
		//if ( geometry.ellipse_list.size() < 140 ) getchar();
		//std::cout << "here_out2 \n";
		geometry.apply_boundary_condition_ellipse();
		//std::cout << "here_out3 \n";
		//std::cout << "ellipse num: " << geometry.ellipse_list.size() << "\n";
		//geometry.write_ellipse_tofile();
		//std::cout << "here_out0 \n";
		geometry.ellipse_to_voronoi();
		//std::cout << "here_out1 \n";
		geometry.update_map();
	    }
	    geometry.output_cell_position();

	    foutput.close();
	    return 0;
	}

    DivInfo cell_division(Geometry &geometry, Force &force)
	{
	    // divided cell
	    int div_index;
	    DivInfo div_info;
	    
	    while (true)
	    {
		div_index = rand() % geometry.voronoi_cell_list.size();
		if ( geometry.voronoi_cell_list[div_index].cell_id < 8000 )
		{
		    std::cout << "dividing cell: " << div_index << "==";
		    break;
		}
	    }
	    std::cout << div_index << "\n";
	    // ***************
	    //div_index = 100;
	    // ***************
	    //getchar();
	    Ellipse dividing_ellipse = geometry.voronoi_to_ellipse(geometry.voronoi_cell_list[div_index]);
	    VoronoiCell new_cell = geometry.voronoi_cell_list[div_index];
	    Ellipse new_ellipse = geometry.voronoi_cell_list[div_index].ellipse;
	    double origin_c1 = new_ellipse.c1;
	    double origin_c2 = new_ellipse.c2;
	    double delta_H = 0.0;
	    double delta_L = 0.0;
	    double dividing_angle;
	    if ( dividing_ellipse.a / dividing_ellipse.b > q0 ) // fix
	    {
		dividing_angle = atan2(dividing_ellipse.v2, dividing_ellipse.v1);
		div_info.type = 1;
	    }
	    else
	    {
		dividing_angle = fmod( double(rand() % 100000), PI );
		div_info.type = 0;
	    }
	    div_info.x = origin_c1;
	    div_info.y = origin_c2;
	    new_ellipse.c1 = origin_c1 + 0.5*cos(dividing_angle); 
	    new_ellipse.c2 = origin_c2 + 0.5*sin(dividing_angle); 
	    geometry.voronoi_cell_list[div_index].ellipse.c1 = origin_c1 - 0.5*cos(dividing_angle);
	    geometry.voronoi_cell_list[div_index].ellipse.c2 = origin_c2 - 0.5*sin(dividing_angle);

	    new_cell.c1 = new_ellipse.c1;
	    new_cell.c2 = new_ellipse.c2;
	    geometry.voronoi_cell_list[div_index].c1 = geometry.voronoi_cell_list[div_index].ellipse.c1;
	    geometry.voronoi_cell_list[div_index].c2 = geometry.voronoi_cell_list[div_index].ellipse.c2;
	    
	    delta_H = fabs(cos(dividing_angle)*scaling_beta1/(scaling_beta1 + scaling_beta2))
		+ fabs(sin(dividing_angle)*scaling_beta2/(scaling_beta1 + scaling_beta2));
	    delta_L = fabs(sin(dividing_angle)*scaling_beta1/(scaling_beta1 + scaling_beta2))
		+ fabs(cos(dividing_angle)*scaling_beta2/(scaling_beta1 + scaling_beta2));
	    delta_H = 4.0*0.25*delta_H / geometry.boundary_data.L;
	    delta_L = 4.0*0.25*delta_L / geometry.boundary_data.H;
	    geometry.boundary_data.set_boundary_data(geometry.boundary_data.L + delta_L, 
						     geometry.boundary_data.H + delta_H);

	    // naming the cell
	    geometry.cell_number = geometry.cell_number + 1;
	    new_ellipse.ellipse_id = geometry.cell_number;
	    new_cell.ellipse = new_ellipse;
	    new_cell.cell_id = new_ellipse.ellipse_id;
	    geometry.voronoi_cell_list.push_back(new_cell);

	    double ratio = 0.0;
	    double total_counter = 0.0;
	    // update cell positions
	    for ( int k=0; k<geometry.voronoi_cell_list.size(); k++ )
	    {
	    	//geometry.output_cell(geometry.voronoi_cell_list[k]);
		if ( geometry.voronoi_cell_list[k].cell_id < 8000)
		//if (geometry.voronoi_cell_list[k].cell_id == 1)
		{
		    total_counter = total_counter + 1.0;
		    Ellipse converted_ellipse = geometry.voronoi_to_ellipse(geometry.voronoi_cell_list[k]);
		    if ( converted_ellipse.a / converted_ellipse.b > q0 )
		    {
			ratio = ratio + 1.0;
		    }
		    double x_coord = geometry.voronoi_cell_list[k].ellipse.c1 - origin_c1;
		    double y_coord = geometry.voronoi_cell_list[k].ellipse.c2 - origin_c2;
		    double x_coord1 = x_coord*cos(dividing_angle) + y_coord*sin(dividing_angle);
		    double y_coord1 = -x_coord*sin(dividing_angle) + y_coord*cos(dividing_angle);
		    /* double x_coord_new1 = x_coord1*( 1.0 +  */
		    /* 				     A1*exp( -( scaling_beta1*fabs(x_coord1) +  */
		    /* 						scaling_beta2*fabs(y_coord1) ) ) ); */
		    /* double y_coord_new1 = y_coord1*( 1.0 +  */
		    /* 				     A2*exp( -( scaling_beta1*fabs(x_coord1) +  */
		    /* 						scaling_beta2*fabs(y_coord1) ) ) ); */
		    double x_coord_new1 = x_coord1*( 1.0 +
		    				     exp( -( scaling_beta1*fabs(x_coord1) +
		    						scaling_beta2*y_coord1*y_coord1 ) ) );
		    double y_coord_new1 = y_coord1;
		    double x_coord_new = x_coord_new1*cos(dividing_angle) - y_coord_new1*sin(dividing_angle);
		    double y_coord_new = x_coord_new1*sin(dividing_angle) + y_coord_new1*cos(dividing_angle);
		    geometry.voronoi_cell_list[k].ellipse.c1 = 
			geometry.voronoi_cell_list[k].ellipse.c1 + x_coord_new - x_coord;
		    geometry.voronoi_cell_list[k].ellipse.c2 = 
			geometry.voronoi_cell_list[k].ellipse.c2 + y_coord_new - y_coord;
		}
	    }
	    div_info.ratio = ratio / total_counter;
	    return div_info;
	}

    int one_step_solution_modify(Geometry &geometry, Force &force, std::string file_index)
	{
	    std::ofstream foutput("out/eigen" + file_index + ".txt");
	    if(!foutput) 
	    {
		std::cout << "file open error.\n";
		return -1; 
	    }

	    for ( int k=0; k<geometry.voronoi_cell_list.size(); k++ )
	    {
	    	//geometry.output_cell(geometry.voronoi_cell_list[k]);
		if ( geometry.voronoi_cell_list[k].cell_id < 8000)
		//if (geometry.voronoi_cell_list[k].cell_id == 1)
		{
		    std::cout << "k: " << k << "cell id: " 
			      << geometry.voronoi_cell_list[k].cell_id << "\n";
		    double temp = sqrt(6.0*PI/sqrt(3));
		    double radius = sqrt(geometry.voronoi_cell_list[k].area/PI);
		    double a = geometry.voronoi_cell_list[k].ellipse.a;
		    double b = geometry.voronoi_cell_list[k].ellipse.b;
		    double sum_force_x = 0.0;
		    double sum_force_y = 0.0;
		    double sum_force2_x = 0.0;
		    double sum_force2_y = 0.0;
		    double sigma11 = 0.0;
		    double sigma12 = 0.0;
		    double sigma21 = 0.0;
		    double sigma22 = 0.0;
		    double oxi = geometry.voronoi_cell_list[k].ellipse.c1;
		    double oyi = geometry.voronoi_cell_list[k].ellipse.c2;
		    for ( int kk=0; kk<geometry.voronoi_cell_list[k].neighbor_id.size(); kk++ )
		    {
			int edge_id = geometry.voronoi_cell_list[k].edge_set[kk];
			int jcell_id = geometry.map_cell_id[geometry.voronoi_cell_list[k].neighbor_id[kk]];
			double eta_hat = eta*geometry.voronoi_edge_list[edge_id].length/dt;
			double cx = geometry.voronoi_edge_list[edge_id].c1;
			double cy = geometry.voronoi_edge_list[edge_id].c2;
			double oxj = geometry.voronoi_cell_list[jcell_id].ellipse.c1;
			double oyj = geometry.voronoi_cell_list[jcell_id].ellipse.c2;
			double theta = geometry.voronoi_cell_list[k].edges[kk].theta;
			sum_force_x = sum_force_x + 
			    force.force_list[edge_id]*geometry.voronoi_cell_list[k].edges[kk].n1;
			sum_force_y = sum_force_y + 
			    force.force_list[edge_id]*geometry.voronoi_cell_list[k].edges[kk].n2;

			sum_force2_x = sum_force2_x + 
			    force.force_list2[edge_id]*geometry.voronoi_cell_list[k].edges[kk].n1;
			sum_force2_y = sum_force2_y + 
			    force.force_list2[edge_id]*geometry.voronoi_cell_list[k].edges[kk].n2;

			sigma11 = sigma11 + force.force_list[edge_id]*cos(theta)*cos(theta);
			sigma12 = sigma12 + force.force_list[edge_id]*cos(theta)*sin(theta);
			sigma21 = sigma21 + force.force_list[edge_id]*sin(theta)*cos(theta);
			sigma22 = sigma22 + force.force_list[edge_id]*sin(theta)*sin(theta);
		    }
		    // geometry.voronoi_cell_list[k].ellipse.c1 = 
		    // 	geometry.voronoi_cell_list[k].ellipse.c1 + sum_force_x*dt;

		    geometry.voronoi_cell_list[k].ellipse.c1 = 
		    	geometry.voronoi_cell_list[k].ellipse.c1 + (sum_force_x + sum_force2_x)*dt;
		    geometry.voronoi_cell_list[k].ellipse.c2 = 
		    	geometry.voronoi_cell_list[k].ellipse.c2 + (sum_force_y + sum_force2_y)*dt;

		    // geometry.voronoi_cell_list[k].ellipse.c1 = 
		    // 	geometry.voronoi_cell_list[k].ellipse.c1 + sum_force_x*dt;
		    // geometry.voronoi_cell_list[k].ellipse.c2 = 
		    // 	geometry.voronoi_cell_list[k].ellipse.c2 + sum_force_y*dt;

		    double trace = sigma11 + sigma22;
		    double det = sigma11*sigma22 - sigma12*sigma21;
		    double lambda1 = 0.5*trace + sqrt(0.25*trace*trace - det);
		    double lambda2 = 0.5*trace - sqrt(0.25*trace*trace - det);
		    double v1, v2;
		    double u1, u2;
		    if ( sigma12 !=0 )
		    {
			v1 = lambda1 - sigma22;
			v2 = sigma12;
			u1 = lambda2 - sigma22;
			u2 = sigma12;
		    }
		    else
		    {
			v1 = 1;
			v2 = 0;
			u1 = 0;
			u2 = 1;
		    }
		    temp = sqrt(v1*v1 + v2*v2);
		    v1 = v1/temp;
		    v2 = v2/temp;
		    temp = sqrt(u1*u1 + u2*u2);
		    u1 = u1/temp;
		    u2 = u2/temp;
		    // geometry.voronoi_cell_list[k].ellipse.c1 = 
		    // 	geometry.voronoi_cell_list[k].ellipse.c1 + sum_force_x*dt;
		    
		    // std::cout << "sigma: " << sigma11 << ", " << sigma12 << ", " << sigma22 << "\n";
		    // std::cout << "lambda: " << lambda1 << ", " << lambda2 << "\n";
		    // std::cout << "v: " << v1 << ", " << v2 << "\n";
		    // std::cout << "u: " << u1 << ", " << u2 << "\n";
		    // getchar();

		    foutput << geometry.voronoi_cell_list[k].cell_id << ", " 
			    << oxi << ", " << oyi << ", " 
			    << lambda1 << ", " << lambda2 << ", " 
			    << oxi+5*lambda1*v1 << ", " << oyi+5*lambda1*v2 << ", " 
			    << oxi+5*lambda2*u1 << ", " << oyi+5*lambda2*u2 << "\n";
		}
	    }
	    foutput.close();
	    return 0;
	}

    // int one_step_solution(Geometry &geometry, Force &force)
    // 	{
    // 	    cholmod_common common;
    // 	    cholmod_sparse *A;
    // 	    cholmod_dense *x, *b, *residual = NULL ;
    // 	    double residual_norm, one [2] = {1,0}, minusone [2] = {-1,0} ;

    // 	    cholmod_l_start(&common);
    // 	    size_t n = geometry.voronoi_cell_list.size();
    // 	    cholmod_dense *A_dense = cholmod_l_zeros(6*n+2, 4*n, CHOLMOD_REAL, &common);
    // 	    if (A_dense == nullptr) {
    // 	    	std::cout << "failed! " << "\n";
    // 	    	return EXIT_FAILURE;
    // 	    }
    // 	    double* v = (double*) A_dense->x;
    // 	    int d = A_dense->d;
    // 	    std::cout << "leading dimension " << d << "rows: " << 6*n+2 << "\n";
    // 	    b = cholmod_l_zeros(A_dense->nrow, 1, A_dense->xtype, &common);
    // 	    double* bv = (double*) b->x;
	    
    // 	    int irow, icol;
    // 	    double mvalue;
    // 	    for ( int k=0; k<geometry.voronoi_cell_list.size(); k++ )
    // 	    {
    // 	    	//geometry.output_cell(geometry.voronoi_cell_list[k]);
    // 		if ( geometry.voronoi_cell_list[k].cell_id < 8000)
    // 		{
    // 		    double temp = sqrt(6.0*PI/sqrt(3));
    // 		    double radius = sqrt(geometry.voronoi_cell_list[k].area/PI);
    // 		    double a = geometry.voronoi_cell_list[k].ellipse.a;
    // 		    double b = geometry.voronoi_cell_list[k].ellipse.b;
    // 		    for ( int kk=0; kk<geometry.voronoi_cell_list[k].neighbor_id.size(); kk++ )
    // 		    {
    // 			int edge_id = geometry.voronoi_cell_list[k].edge_set[kk];
    // 			int jcell_id = geometry.map_cell_id[geometry.voronoi_cell_list[k].neighbor_id[kk]];
    // 			double eta_hat = eta*geometry.voronoi_edge_list[edge_id].length/dt;
    // 			double cx = geometry.voronoi_edge_list[edge_id].c1;
    // 			double cy = geometry.voronoi_edge_list[edge_id].c2;
    // 			double oxi = geometry.voronoi_cell_list[k].ellipse.c1;
    // 			double oyi = geometry.voronoi_cell_list[k].ellipse.c2;
    // 			double oxj = geometry.voronoi_cell_list[jcell_id].ellipse.c1;
    // 			double oyj = geometry.voronoi_cell_list[jcell_id].ellipse.c2;
    // 			double theta = geometry.voronoi_cell_list[k].edges[kk].theta;

    // 			//std::cout << theta << " ";
    // 			// 6i
    // 			irow = 6*k;
    // 			bv[irow] = bv[irow] + force.force_list[edge_id]*
    // 			    geometry.voronoi_cell_list[k].edges[kk].n1;

    // 			icol = 4*k;
    // 			mvalue = -eta_hat;
    // 			v[irow+icol*d] = v[irow+icol*d] + mvalue;

    // 			icol = 4*k+2;
    // 			mvalue = -eta_hat*(cx - oxi);
    // 			v[irow+icol*d] = v[irow+icol*d] + mvalue;

    // 			icol = 4*k+3;
    // 			mvalue = -eta_hat*(cy - oyi);
    // 			v[irow+icol*d] = v[irow+icol*d] + mvalue;

    // 			// 6i+2
    // 			irow = 6*k+2;
    // 			bv[irow] = bv[irow] + force.force_list[edge_id]*
    // 			    geometry.voronoi_cell_list[k].edges[kk].n1*cos(theta)
    // 			    - temp*radius*a;

    // 			icol = 4*k;
    // 			mvalue = -eta_hat*cos(theta);
    // 			v[irow+icol*d] = v[irow+icol*d] + mvalue;

    // 			icol = 4*k+2;
    // 			//mvalue = -eta_hat*(cx - oxi)*cos(theta) + temp*radius;
    // 			mvalue = -eta_hat*(cx - oxi)*cos(theta);
    // 			v[irow+icol*d] = v[irow+icol*d] + mvalue;

    // 			icol = 4*k+3;
    // 			mvalue = -eta_hat*(cy - oyi)*cos(theta);
    // 			v[irow+icol*d] = v[irow+icol*d] + mvalue;

    // 			// 6i+3
    // 			irow = 6*k+3;
    // 			bv[irow] = bv[irow] + force.force_list[edge_id]*
    // 			    geometry.voronoi_cell_list[k].edges[kk].n1*sin(theta)
    // 			    -temp*radius*b;

    // 			icol = 4*k;
    // 			mvalue = -eta_hat*sin(theta);
    // 			v[irow+icol*d] = v[irow+icol*d] + mvalue;

    // 			icol = 4*k+2;
    // 			mvalue = -eta_hat*(cx - oxi)*sin(theta);
    // 			v[irow+icol*d] = v[irow+icol*d] + mvalue;

    // 			icol = 4*k+3;
    // 			//mvalue = -eta_hat*(cy - oyi)*sin(theta) + temp*radius;
    // 			mvalue = -eta_hat*(cy - oyi)*sin(theta);
    // 			v[irow+icol*d] = v[irow+icol*d] + mvalue;

    // 			// 6i+1
    // 			irow = 6*k+1;
    // 			bv[irow] = bv[irow] + force.force_list[edge_id]*
    // 			    geometry.voronoi_cell_list[k].edges[kk].n2;

    // 			icol = 4*k+1;
    // 			mvalue = -eta_hat;
    // 			v[irow+icol*d] = v[irow+icol*d] + mvalue;

    // 			icol = 4*k+3;
    // 			mvalue = -eta_hat*(cx - oxi);
    // 			v[irow+icol*d] = v[irow+icol*d] + mvalue;

    // 			icol = 4*k+2;
    // 			mvalue = eta_hat*(cy - oyi);
    // 			v[irow+icol*d] = v[irow+icol*d] + mvalue;
    // 			// 6i+4
    // 			irow = 6*k+4;
    // 			bv[irow] = bv[irow] + force.force_list[edge_id]*
    // 			    geometry.voronoi_cell_list[k].edges[kk].n2*cos(theta)
    // 			    - temp*radius*b;

    // 			icol = 4*k+1;
    // 			mvalue = -eta_hat*cos(theta);
    // 			v[irow+icol*d] = v[irow+icol*d] + mvalue;

    // 			icol = 4*k+3;
    // 			//mvalue = -eta_hat*(cx - oxi)*cos(theta) + temp*radius;
    // 			mvalue = -eta_hat*(cx - oxi)*cos(theta);
    // 			v[irow+icol*d] = v[irow+icol*d] + mvalue;

    // 			icol = 4*k+2;
    // 			mvalue = eta_hat*(cy - oyi)*cos(theta);
    // 			v[irow+icol*d] = v[irow+icol*d] + mvalue;

    // 			// 6i+5
    // 			irow = 6*k+5;
    // 			bv[irow] = bv[irow] + force.force_list[edge_id]*
    // 			    geometry.voronoi_cell_list[k].edges[kk].n2*sin(theta)
    // 			    + temp*radius*a;

    // 			icol = 4*k+1;
    // 			mvalue = -eta_hat*sin(theta);
    // 			v[irow+icol*d] = v[irow+icol*d] + mvalue;

    // 			icol = 4*k+3;
    // 			mvalue = -eta_hat*(cx - oxi)*sin(theta);
    // 			v[irow+icol*d] = v[irow+icol*d] + mvalue;

    // 			icol = 4*k+2;
    // 			//mvalue = eta_hat*(cy - oyi)*sin(theta) - temp*radius;
    // 			mvalue = eta_hat*(cy - oyi)*sin(theta);
    // 			v[irow+icol*d] = v[irow+icol*d] + mvalue;

    // 			if (geometry.voronoi_cell_list[k].neighbor_id[kk] < 8000)
    // 			{
    // 			// 6i
    // 			irow = 6*k;
    // 			icol = 4*jcell_id;
    // 			mvalue = eta_hat;
    // 			v[irow+icol*d] = v[irow+icol*d] + mvalue;

    // 			icol = 4*jcell_id+2;
    // 			mvalue = eta_hat*(cx - oxj);
    // 			v[irow+icol*d] = v[irow+icol*d] + mvalue;
		    
    // 			icol = 4*jcell_id+3;
    // 			mvalue = eta_hat*(cy - oyj);
    // 			v[irow+icol*d] = v[irow+icol*d] + mvalue;
		    
    // 			// 6i+2
    // 			irow = 6*k+2;
    // 			icol = 4*jcell_id;
    // 			mvalue = eta_hat*cos(theta);
    // 			v[irow+icol*d] = v[irow+icol*d] + mvalue;

    // 			icol = 4*jcell_id+2;
    // 			mvalue = eta_hat*(cx - oxj)*cos(theta);
    // 			v[irow+icol*d] = v[irow+icol*d] + mvalue;
		    
    // 			icol = 4*jcell_id+3;
    // 			mvalue = eta_hat*(cy - oyj)*cos(theta);
    // 			v[irow+icol*d] = v[irow+icol*d] + mvalue;

    // 			// 6i+3
    // 			irow = 6*k+3;
    // 			icol = 4*jcell_id;
    // 			mvalue = eta_hat*sin(theta);
    // 			v[irow+icol*d] = v[irow+icol*d] + mvalue;

    // 			icol = 4*jcell_id+2;
    // 			mvalue = eta_hat*(cx - oxj)*sin(theta);
    // 			v[irow+icol*d] = v[irow+icol*d] + mvalue;
		    
    // 			icol = 4*jcell_id+3;
    // 			mvalue = eta_hat*(cy - oyj)*sin(theta);
    // 			v[irow+icol*d] = v[irow+icol*d] + mvalue;

    // 			// 6i+1
    // 			irow = 6*k+1;
    // 			icol = 4*jcell_id+1;
    // 			mvalue = eta_hat;
    // 			v[irow+icol*d] = v[irow+icol*d] + mvalue;

    // 			icol = 4*jcell_id+3;
    // 			mvalue = eta_hat*(cx - oxj);
    // 			v[irow+icol*d] = v[irow+icol*d] + mvalue;
		    
    // 			icol = 4*jcell_id+2;
    // 			mvalue = -eta_hat*(cy - oyj);
    // 			v[irow+icol*d] = v[irow+icol*d] + mvalue;
		    
    // 			// 6i+4
    // 			irow = 6*k+4;
    // 			icol = 4*jcell_id+1;
    // 			mvalue = eta_hat*cos(theta);
    // 			v[irow+icol*d] = v[irow+icol*d] + mvalue;

    // 			icol = 4*jcell_id+3;
    // 			mvalue = eta_hat*(cx - oxj)*cos(theta);
    // 			v[irow+icol*d] = v[irow+icol*d] + mvalue;
		    
    // 			icol = 4*jcell_id+2;
    // 			mvalue = -eta_hat*(cy - oyj)*cos(theta);
    // 			v[irow+icol*d] = v[irow+icol*d] + mvalue;

    // 			// 6i+5
    // 			irow = 6*k+5;
    // 			icol = 4*jcell_id+1;
    // 			mvalue = eta_hat*sin(theta);
    // 			v[irow+icol*d] = v[irow+icol*d] + mvalue;

    // 			icol = 4*jcell_id+3;
    // 			mvalue = eta_hat*(cx - oxj)*sin(theta);
    // 			v[irow+icol*d] = v[irow+icol*d] + mvalue;
		    
    // 			icol = 4*jcell_id+2;
    // 			mvalue = -eta_hat*(cy - oyj)*sin(theta);
    // 			v[irow+icol*d] = v[irow+icol*d] + mvalue;
    // 			}
    // 		    }
    // 		}
    // 	    }
    // 	    A = cholmod_l_dense_to_sparse(A_dense, true, &common);
    // 	    cholmod_l_print_sparse(A, "sparse", &common);
    // 	    // FILE *pf_solution;
    // 	    // pf_solution = fopen ("solution.txt","w");
    // 	    // if (pf_solution == NULL) 
    // 	    // {
    // 	    // 	std::cout << " open file failed \n";
    // 	    // 	getchar();
    // 	    // }
    // 	    // cholmod_l_write_dense(pf_solution, A_dense, "dense", &common);
	    
    // 	    x = SuiteSparseQR<double>(A, b, &common);
    // 	    double* xv = (double*) x->x;
    // 	    for ( int k=0; k < geometry.voronoi_cell_list.size(); k++ )
    // 	    {
    // 		geometry.voronoi_cell_list[k].ellipse.c1 = geometry.voronoi_cell_list[k].ellipse.c1 + xv[4*k];
    // 		geometry.voronoi_cell_list[k].ellipse.c2 = geometry.voronoi_cell_list[k].ellipse.c2 + xv[4*k+1];
    // 		geometry.voronoi_cell_list[k].ellipse.a = geometry.voronoi_cell_list[k].ellipse.a + xv[4*k+2];
    // 		geometry.voronoi_cell_list[k].ellipse.b = geometry.voronoi_cell_list[k].ellipse.b + xv[4*k+3];
    // 		if  (xv[4*k+1] > 50.0)
    // 		{
    // 		    std::cout << "k: " << k << " - " << 
    // 			geometry.voronoi_cell_list[k].cell_id << "\n";
    // 		    std::cout << "ellipse k: " << k << " - " << 
    // 			geometry.voronoi_cell_list[k].ellipse.ellipse_id << "\n";
    // 		    geometry.output_cell(geometry.voronoi_cell_list[k]);		    
    // 		    getchar();
    // 		}
    // 	    }

    // 	    // cholmod_l_print_dense(x, "dense", &common);
    // 	    // FILE *pf_solution;
    // 	    // pf_solution = fopen ("solution.txt","w");
    // 	    // if (pf_solution == NULL) 
    // 	    // {
    // 	    // 	std::cout << " open file failed \n";
    // 	    // 	getchar();
    // 	    // }
    // 	    // cholmod_l_write_dense(pf_solution, b, "dense", &common);
	    
    // 	    residual = cholmod_l_copy_dense(b, &common);
    // 	    cholmod_l_sdmult(A, 0, minusone, one, x, residual, &common);
    // 	    residual_norm = cholmod_l_norm_dense(residual, 2, &common) ;
	    
    // 	    std::cout << "|| A x - b ||_2 = " << residual_norm << "\n";

    // 	    cholmod_l_free_dense(&residual, &common);
    // 	    cholmod_l_free_dense(&A_dense, &common);
    // 	    cholmod_l_free_sparse(&A, &common);
    // 	    cholmod_l_free_dense(&x, &common);
    // 	    cholmod_l_free_dense(&b, &common);
    // 	    cholmod_l_finish(&common);
    // 	}

    // int one_step_solution_simple(Geometry &geometry, Force &force)
    // 	{
    // 	    cholmod_common common;
    // 	    cholmod_sparse *A;
    // 	    cholmod_dense *x, *b, *residual = NULL ;
    // 	    double residual_norm, one [2] = {1,0}, minusone [2] = {-1,0} ;

    // 	    cholmod_l_start(&common);
    // 	    size_t n = geometry.voronoi_cell_list.size();
    // 	    cholmod_dense *A_dense = cholmod_l_zeros(2*n+2, 2*n, CHOLMOD_REAL, &common);
    // 	    if (A_dense == nullptr) {
    // 	    	std::cout << "failed! " << "\n";
    // 	    	return EXIT_FAILURE;
    // 	    }
    // 	    double* v = (double*) A_dense->x;
    // 	    int d = A_dense->d;
    // 	    std::cout << "leading dimension " << d << "rows: " << 2*n+2 << "\n";
    // 	    b = cholmod_l_zeros(A_dense->nrow, 1, A_dense->xtype, &common);
    // 	    double* bv = (double*) b->x;
	    
    // 	    int irow, icol;
    // 	    double mvalue;
    // 	    for ( int k=0; k<geometry.voronoi_cell_list.size(); k++ )
    // 	    {
    // 	    	//geometry.output_cell(geometry.voronoi_cell_list[k]);
    // 		if ( geometry.voronoi_cell_list[k].cell_id < 8000)
    // 		//if (geometry.voronoi_cell_list[k].cell_id == 1)
    // 		{
    // 		    std::cout << "k" << k << "\n";
    // 		    double temp = sqrt(6.0*PI/sqrt(3));
    // 		    double radius = sqrt(geometry.voronoi_cell_list[k].area/PI);
    // 		    double a = geometry.voronoi_cell_list[k].ellipse.a;
    // 		    double b = geometry.voronoi_cell_list[k].ellipse.b;
    // 		    for ( int kk=0; kk<geometry.voronoi_cell_list[k].neighbor_id.size(); kk++ )
    // 		    {
    // 			int edge_id = geometry.voronoi_cell_list[k].edge_set[kk];
    // 			int jcell_id = geometry.map_cell_id[geometry.voronoi_cell_list[k].neighbor_id[kk]];
    // 			double eta_hat = eta*geometry.voronoi_edge_list[edge_id].length/dt;
    // 			double cx = geometry.voronoi_edge_list[edge_id].c1;
    // 			double cy = geometry.voronoi_edge_list[edge_id].c2;
    // 			double oxi = geometry.voronoi_cell_list[k].ellipse.c1;
    // 			double oyi = geometry.voronoi_cell_list[k].ellipse.c2;
    // 			double oxj = geometry.voronoi_cell_list[jcell_id].ellipse.c1;
    // 			double oyj = geometry.voronoi_cell_list[jcell_id].ellipse.c2;
    // 			double theta = geometry.voronoi_cell_list[k].edges[kk].theta;

    // 			//std::cout << theta << " ";
    // 			// 2i
    // 			irow = 2*k;
    // 			bv[irow] = bv[irow] + force.force_list[edge_id]*
    // 			    geometry.voronoi_cell_list[k].edges[kk].n1;
    // 			std::cout << "b: " << irow << " -> " << bv[irow] << "\n";

    // 			icol = 2*k;
    // 			//mvalue = -eta_hat;
    // 			mvalue = eta_hat;
    // 			v[irow+icol*d] = v[irow+icol*d] + mvalue;
    // 			//std::cout << "v: " << irow << "-" << icol << "->" 
    // 			//	  << v[irow+icol*d] << "\n";

    // 			// 2i+1
    // 			irow = 2*k+1;
    // 			bv[irow] = bv[irow] + force.force_list[edge_id]*
    // 			    geometry.voronoi_cell_list[k].edges[kk].n2;
    // 			std::cout << "b: " << irow << " -> " << bv[irow] << "\n";

    // 			icol = 2*k+1;
    // 			//mvalue = -eta_hat;
    // 			mvalue = eta_hat;
    // 			v[irow+icol*d] = v[irow+icol*d] + mvalue;
    // 			//std::cout << "v: " << irow << "-" << icol << "->" 
    // 			//	  << v[irow+icol*d] << "\n";

    // 			if ( geometry.voronoi_cell_list[k].cell_id < 8000)
    // 			    //if (geometry.voronoi_cell_list[k].neighbor_id[kk] == 1)
    // 			    //if (geometry.voronoi_cell_list[k].cell_id == 1)
    // 			{
    // 			// 2i
    // 			irow = 2*k;
    // 			icol = 2*jcell_id;
    // 			//mvalue = eta_hat;
    // 			mvalue = -eta_hat;
    // 			v[irow+icol*d] = v[irow+icol*d] + mvalue;
    // 			std::cout << "v: " << irow << "-" << icol << "->" 
    // 				  << v[irow+icol*d] << "\n";

    // 			// 2i+1
    // 			irow = 2*k+1;
    // 			icol = 2*jcell_id+1;
    // 			//mvalue = eta_hat;
    // 			mvalue = -eta_hat;
    // 			v[irow+icol*d] = v[irow+icol*d] + mvalue;
    // 			std::cout << "v: " << irow << "-" << icol << "->" 
    // 				  << v[irow+icol*d] << "\n";
    // 			}
    // 		    }
    // 		}
    // 		else
    // 		{
    // 		    std::cout << "not k " << k << "\n";
    // 		}
    // 	    }
    // 	    getchar();
    // 	    A = cholmod_l_dense_to_sparse(A_dense, true, &common);
    // 	    cholmod_l_print_sparse(A, "sparse", &common);
    // 	    // FILE *pf_solution;
    // 	    // pf_solution = fopen ("solution.txt","w");
    // 	    // if (pf_solution == NULL) 
    // 	    // {
    // 	    // 	std::cout << " open file failed \n";
    // 	    // 	getchar();
    // 	    // }
    // 	    // cholmod_l_write_dense(pf_solution, A_dense, "dense", &common);
	    
    // 	    x = SuiteSparseQR<double>(A, b, &common);
    // 	    double* xv = (double*) x->x;
    // 	    double sum1 = 0.0;
    // 	    double sum2 = 0.0;
    // 	    for ( int k=0; k<x->nrow/2; k++)
    // 	    {
    // 		std::cout << "k: " << k << "xv->" << xv[2*k] << ", " << xv[2*k+1] << "\n";
    // 		sum1 = sum1 + xv[2*k];
    // 		sum2 = sum2 + xv[2*k+1];
    // 	    }
    // 	    std::cout << "sum" << sum1 << ", " << sum2 << "\n";	    
    // 	    getchar();
    // 	    for ( int k=0; k < geometry.voronoi_cell_list.size(); k++ )
    // 	    {
    // 		geometry.voronoi_cell_list[k].ellipse.c1 = geometry.voronoi_cell_list[k].ellipse.c1 + xv[2*k];
    // 		geometry.voronoi_cell_list[k].ellipse.c2 = geometry.voronoi_cell_list[k].ellipse.c2 + xv[2*k+1];
    // 		if  (xv[2*k+1] > 50.0)
    // 		{
    // 		    std::cout << "k: " << k << " - " << 
    // 			geometry.voronoi_cell_list[k].cell_id << "\n";
    // 		    std::cout << "ellipse k: " << k << " - " << 
    // 			geometry.voronoi_cell_list[k].ellipse.ellipse_id << "\n";
    // 		    geometry.output_cell(geometry.voronoi_cell_list[k]);		    
    // 		    getchar();
    // 		}
    // 	    }

    // 	    residual = cholmod_l_copy_dense(b, &common);
    // 	    cholmod_l_sdmult(A, 0, minusone, one, x, residual, &common);
    // 	    residual_norm = cholmod_l_norm_dense(residual, 2, &common) ;
	    
    // 	    std::cout << "|| A x - b ||_2 = " << residual_norm << "\n";

    // 	    cholmod_l_free_dense(&residual, &common);
    // 	    cholmod_l_free_dense(&A_dense, &common);
    // 	    cholmod_l_free_sparse(&A, &common);
    // 	    cholmod_l_free_dense(&x, &common);
    // 	    cholmod_l_free_dense(&b, &common);
    // 	    cholmod_l_finish(&common);
    // 	}

    // int test()
    // 	{
    // 	    cholmod_common common, *cc ;
    // 	    cholmod_sparse *A ;
    // 	    cholmod_dense *x, *b, *residual = NULL ;
    // 	    double residual_norm, one [2] = {1,0}, minusone [2] = {-1,0} ;

    // 	    cholmod_l_start(&common);
    // 	    {
    // 		const size_t n = 1000;
    // 		const size_t nnz = 2 + 2 * (n - 2);
    // 		cholmod_triplet *T = cholmod_l_allocate_triplet(n, n, nnz, 0, CHOLMOD_REAL, &common);
    // 		if (T == nullptr) {
    // 		    perror("cholmod_l_allocate_triplet");
    // 		    return EXIT_FAILURE;
    // 		}
		
    // 		size_t k = 0;
    // 		long* i = (long*) T->i;
    // 		long* j = (long*) T->j;
    // 		double* v = (double*) T->x;

    // 		i[k] = 0; j[k] = 1; v[k] = 1.0; ++k;
    // 		for (size_t row = 1; row < n-1; ++row) {
    // 		    i[k] = row; j[k] = row - 1; v[k] = 0.5; ++k;
    // 		    i[k] = row; j[k] = row + 1; v[k] = 0.5; ++k;
    // 		}
    // 		i[k] = n - 1; j[k] = n - 2; v[k] = 1.0; ++k;

    // 		T->nnz = k;

    // 		cholmod_l_print_triplet(T, "triplet", &common);

    // 		A = cholmod_l_triplet_to_sparse(T, T->nnz, &common);
    // 		if (A == nullptr) {
    // 		    perror("cholmod_l_triplet_to_sparse");
    // 		    return EXIT_FAILURE;
    // 		}

    // 		cholmod_l_print_sparse(A, "sparse", &common);
    // 		cholmod_l_free_triplet(&T, &common);
    // 	    }
    // 	    b = cholmod_l_ones(A->nrow, 1, A->xtype, &common);

    // 	    x = SuiteSparseQR<double>(A, b, &common);
	    
    // 	    residual = cholmod_l_copy_dense(b, &common);
    // 	    cholmod_l_sdmult(A, 0, minusone, one, x, residual, &common);
    // 	    residual_norm = cholmod_l_norm_dense(residual, 2, &common) ;

    // 	    std::cout << "|| A x - b ||_2 = " << residual_norm << "\n";

    // 	    cholmod_l_free_dense(&residual, &common);
    // 	    cholmod_l_free_sparse(&A, &common);
    // 	    cholmod_l_free_dense(&x, &common);
    // 	    cholmod_l_free_dense(&b, &common);
    // 	    cholmod_l_finish(&common);
    // 	}

};

#endif //ENGINE
