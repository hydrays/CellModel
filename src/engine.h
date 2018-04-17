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
#include "parameters.h"
#include <time.h>
#include <random>

struct DivInfo
{
    double x;
    double y;
    int type;
    double ratio;
    double dividing_angle;
    int counter = 0;
    int div_index;
};

struct SysInfo
{
    double ratio;
    double A;
};

class Engine
{
    std::random_device rd;
    std::mt19937 gen{rd()};
    std::uniform_real_distribution<> runif{0.0, 1.0};
    std::normal_distribution<> rnorm{0.0, 1.0};
    std::uniform_int_distribution<> rid{1, 32767};
    double p0 = 0.6;
    double sig = 0.0;
    double s1 = 0.0;
    double s2 = 0.0;
    
public:    
    int run(Geometry &geometry, Force &force, Parameters &params, int process_id)
    {
	std::ofstream foutput(params.output_dir + "/statistics" + std::to_string(process_id) +".txt");
	if(!foutput) 
	{
	    std::cout << "file open error.\n";
	    return -1; 
	}

	DivInfo div_info;
	SysInfo sys_info;

	for ( int istep=0; istep<params.Nstep; istep++ )
	{
	    if ( params.flag_output_cell )
	    {
		geometry.output_cell_list(std::to_string(istep), params);
		geometry.output_ellipse_list(std::to_string(istep), params);
	    }
	    force.update_force_field(geometry, params);
	    if ( params.flag_output_cell )
	    {
		force.output_force_filed(geometry, std::to_string(istep), params);
	    }
	    update_cell_position(geometry, force, params, std::to_string(istep));

	    cell_division(geometry, force, params, div_info);
	    collect_info(geometry, force, params, sys_info);

	    foutput << istep << ", " << div_info.counter << ", "
		    << div_info.div_index << ", "
		    << geometry.boundary_data.L_min << ", "
		    << geometry.boundary_data.L_max << ", "
		    << geometry.boundary_data.H << ", "
		    << div_info.x << ", "
		    << div_info.y << ", "
		    << div_info.type << ", "
		    << sys_info.ratio << ", "
		    << sys_info.A << ", "
		    << div_info.dividing_angle << ", "
		    << sig << "\n";
	    foutput.flush();
	    
	    //geometry.apply_periodic_boundary_condition_voronoi();
	    geometry.apply_mixed_boundary_condition_voronoi();
	    geometry.update_ellipse_list();

	    //geometry.output_ellipse_list(std::to_string(100));
	    //geometry.write_ellipse_tofile();

	    std::cout << "step: " << istep
		      << " ellipse num: "
		      << geometry.ellipse_list.size() << "\n";
	    geometry.apply_mixed_boundary_condition_ellipse();
	    geometry.ellipse_to_voronoi(params);
	    geometry.update_map();
	}
	geometry.output_cell_position(params);

	foutput.close();
	return 0;
    }

    int cell_division(Geometry &geometry, Force &force, Parameters &params, DivInfo &div_info)
    {
	// divided cell
	div_info.div_index = -1;

	double u = runif(gen);
	if ( u < params.mu*params.dt*geometry.voronoi_cell_list.size() )
	{
	    div_info.div_index = rid(gen) % geometry.voronoi_cell_list.size();
	    if ( geometry.voronoi_cell_list[div_info.div_index].cell_id < 8000 &&
		 geometry.voronoi_cell_list[div_info.div_index].ellipse.type < 10)
	    {
		std::cout << "dividing cell: " << div_info.div_index << "\n";
		div_info.counter = div_info.counter + 1;

		VoronoiCell new_cell = geometry.voronoi_cell_list[div_info.div_index];
		Ellipse new_ellipse = geometry.voronoi_cell_list[div_info.div_index].ellipse;
		double origin_c1 = new_ellipse.c1;
		double origin_c2 = new_ellipse.c2;
		double delta_H = 0.0;
		double delta_L = 0.0;
		double dividing_angle;

		Ellipse converted_ellipse =
		    geometry.voronoi_to_ellipse(geometry.voronoi_cell_list[div_info.div_index]);	
		double s = converted_ellipse.a/converted_ellipse.b - params.q0;

		if ( geometry.voronoi_cell_list[div_info.div_index].ellipse.type == 2 ) // rotate
		{
		    dividing_angle = runif(gen)*PI;
		    div_info.type = 0;
		}
		else
		{
		    if ( params.two_population_model==1 && s<=0.0 ) // rotate
			//if ( 0 )
		    {
			dividing_angle = runif(gen)*PI;
			div_info.type = 0;
		    }
		    else
		    {
			dividing_angle = atan2(converted_ellipse.v2, converted_ellipse.v1);
			dividing_angle = dividing_angle;// + rnorm(gen)*PI*15.37/180;
			div_info.type = 1;
		    }
		}
		div_info.dividing_angle = dividing_angle;
		div_info.x = origin_c1;
		div_info.y = origin_c2;
		if ( div_info.type == 0 )
		{
		    new_ellipse.c1 = origin_c1 + 0.5*cos(dividing_angle); 
		    new_ellipse.c2 = origin_c2 + 0.5*sin(dividing_angle); 
		    geometry.voronoi_cell_list[div_info.div_index].ellipse.c1 = origin_c1 - 0.5*cos(dividing_angle);
		    geometry.voronoi_cell_list[div_info.div_index].ellipse.c2 = origin_c2 - 0.5*sin(dividing_angle);
		}
		else
		{
		    new_ellipse.c1 = origin_c1 + 0.5*cos(dividing_angle); 
		    new_ellipse.c2 = origin_c2 + 0.5*sin(dividing_angle); 
		    geometry.voronoi_cell_list[div_info.div_index].ellipse.c1 = origin_c1 - 0.5*cos(dividing_angle);
		    geometry.voronoi_cell_list[div_info.div_index].ellipse.c2 = origin_c2 - 0.5*sin(dividing_angle);
		}

		// naming the cell
		geometry.cell_number = geometry.cell_number + 1;
		new_ellipse.ellipse_id = geometry.cell_number;
		new_cell.ellipse = new_ellipse;
		new_cell.cell_id = new_ellipse.ellipse_id;
		geometry.voronoi_cell_list.push_back(new_cell);

		// compute feedback signal
		s1 = s1 + fabs(cos(dividing_angle));
		s2 = s2 + fabs(sin(dividing_angle));
		sig = (1-p0)*s1 - p0*s2;
	    }
	    else
	    {
		div_info.div_index = -1;
	    }
	}
	else
	{
	    div_info.div_index = -1;
	}
	return 0;
    }

    int boundary_cell_division(Geometry &geometry, Force &force, int bc_type)
    {
	int div_index;
	while (true)
	{
	    div_index = rid(gen) % geometry.voronoi_cell_list.size();
	    if ( geometry.voronoi_cell_list[div_index].cell_id < 8000 &&
		 geometry.voronoi_cell_list[div_index].ellipse.type == bc_type)
	    {
		std::cout << "[adding boundary cell: " << div_index << "]\n";
		//getchar();
		break;
	    }
	}
	VoronoiCell new_cell = geometry.voronoi_cell_list[div_index];
	Ellipse new_ellipse = geometry.voronoi_cell_list[div_index].ellipse;
	double origin_c1 = new_ellipse.c1;
	double origin_c2 = new_ellipse.c2;
	double delta_H = 0.0;
	double delta_L = 0.0;
	double dividing_angle;

	dividing_angle = PI/2;
	new_ellipse.c1 = origin_c1 + 0.5*cos(dividing_angle); 
	new_ellipse.c2 = origin_c2 + 0.5*sin(dividing_angle); 
	geometry.voronoi_cell_list[div_index].ellipse.c1 = origin_c1 - 0.5*cos(dividing_angle);
	geometry.voronoi_cell_list[div_index].ellipse.c2 = origin_c2 - 0.5*sin(dividing_angle);

	geometry.cell_number = geometry.cell_number + 1;
	new_ellipse.ellipse_id = geometry.cell_number;
	new_cell.ellipse = new_ellipse;
	new_cell.cell_id = new_ellipse.ellipse_id;
	geometry.voronoi_cell_list.push_back(new_cell);

	return 0;
    }

    int collect_info(Geometry &geometry, Force &force, Parameters &params, SysInfo &sys_info)
    {
	double ratio = 0.0;
	double A = 0.0;
	double total_counter = 0.0;
	// update cell positions
	for ( int k=0; k<geometry.voronoi_cell_list.size(); k++ )
	{
	    //geometry.output_cell(geometry.voronoi_cell_list[k]);
	    /* if ( geometry.voronoi_cell_list[k].cell_id < 8000 && */
	    /* 	 geometry.voronoi_cell_list[k].ellipse.type < 10 ) */
	    if ( geometry.voronoi_cell_list[k].cell_id < 8000 &&
		 geometry.voronoi_cell_list[k].ellipse.type == 1)
	    {
		total_counter = total_counter + 1.0;
		Ellipse converted_ellipse = geometry.voronoi_to_ellipse(geometry.voronoi_cell_list[k]);
		double phi;
		//double phi = atan(fabs(geometry.voronoi_cell_list[k].ellipse.v2/
		//		       geometry.voronoi_cell_list[k].ellipse.v1));
		if ( params.two_population_model==1 && converted_ellipse.a/converted_ellipse.b < params.q0 )
		{
		    phi = runif(gen)*PI*0.5;
		}
		else
		{
		    phi = atan(fabs(converted_ellipse.v2/converted_ellipse.v1));
		}
		A = A + sin(phi) - cos(phi);
		//printf("v1: %f; v2: %f, phi %f \n", converted_ellipse.v1, converted_ellipse.v2, phi);
		if ( converted_ellipse.a / converted_ellipse.b > params.q0 )
		{
		    ratio = ratio + 1.0;
		}		    
	    }
	}
	ratio = ratio / total_counter;
	A = A / total_counter;
	//printf("total_counter: %f; A: %f \n", total_counter, A);
	sys_info.ratio = ratio;
	sys_info.A = A;
	return 0;
    }

    int update_cell_position(Geometry &geometry, Force &force, Parameters &params, std::string file_index)
    {
	double L_min_old = geometry.boundary_data.L_min;
	double L_max_old = geometry.boundary_data.L_max;
	double L_min = 0.0;
	double L_max = 0.0;
	double delta_H = 0.0;
	int L_min_counter = 0;
	int L_max_counter = 0;
	double external_force_y = 0.0;
	double boundary_force1_y = 0.0;
	double boundary_force2_y = 0.0;
	for ( int k=0; k<geometry.voronoi_cell_list.size(); k++ )
	{
	    if ( geometry.voronoi_cell_list[k].cell_id < 8000)
	    {
		/* std::cout << "k: " << k << "cell id: "  */
		/* 	      << geometry.voronoi_cell_list[k].cell_id << "\n"; */
		/* double temp = sqrt(6.0*PI/sqrt(3)); */
		/* double radius = sqrt(geometry.voronoi_cell_list[k].area/PI); */
		/* double a = geometry.voronoi_cell_list[k].ellipse.a; */
		/* double b = geometry.voronoi_cell_list[k].ellipse.b; */
		double sum_force_x = 0.0;
		double sum_force_y = 0.0;
		double sum_force2_x = 0.0;
		double sum_force2_y = 0.0;
		double sum_force3_x = 0.0;
		double sum_force3_y = 0.0;
		double torque = 0.0;
		/* double sigma11 = 0.0; */
		/* double sigma12 = 0.0; */
		/* double sigma21 = 0.0; */
		/* double sigma22 = 0.0; */
		double oxi = geometry.voronoi_cell_list[k].ellipse.c1;
		double oyi = geometry.voronoi_cell_list[k].ellipse.c2;
		for ( int kk=0; kk<geometry.voronoi_cell_list[k].neighbor_id.size(); kk++ )
		{
		    int edge_id = geometry.voronoi_cell_list[k].edge_set[kk];
		    //int jcell_id = geometry.map_cell_id[geometry.voronoi_cell_list[k].neighbor_id[kk]];
		    //double eta_hat = eta*geometry.voronoi_edge_list[edge_id].length/dt;
		    double cx = geometry.voronoi_edge_list[edge_id].c1;
		    double cy = geometry.voronoi_edge_list[edge_id].c2;
		    //double oxj = geometry.voronoi_cell_list[jcell_id].ellipse.c1;
		    //double oyj = geometry.voronoi_cell_list[jcell_id].ellipse.c2;
		    //double theta = geometry.voronoi_cell_list[k].edges[kk].theta;
		    torque = torque +
		    	(force.force_list[edge_id]+force.force_list2[edge_id])
		    	*geometry.voronoi_cell_list[k].edges[kk].n2*(cx-oxi) -
		    	(force.force_list[edge_id]+force.force_list2[edge_id])
		    	*geometry.voronoi_cell_list[k].edges[kk].n1*(cy-oyi);
		    /* torque = torque + */
		    /* 	force.force_list[edge_id] */
		    /* 	*geometry.voronoi_cell_list[k].edges[kk].n2*(cx-oxi) - */
		    /* 	force.force_list[edge_id] */
		    /* 	*geometry.voronoi_cell_list[k].edges[kk].n1*(cy-oyi); */

		    sum_force_x = sum_force_x + 
			force.force_list[edge_id]*geometry.voronoi_cell_list[k].edges[kk].n1;
		    sum_force_y = sum_force_y + 
			force.force_list[edge_id]*geometry.voronoi_cell_list[k].edges[kk].n2;

		    sum_force2_x = sum_force2_x + 
			force.force_list2[edge_id]*geometry.voronoi_cell_list[k].edges[kk].n1;
		    sum_force2_y = sum_force2_y + 
			force.force_list2[edge_id]*geometry.voronoi_cell_list[k].edges[kk].n2;

		    // compute vertical enternal force
		    if ( geometry.voronoi_cell_list[k].ellipse.type == 10 )
		    {
			boundary_force1_y = boundary_force1_y +
			    force.force_list[edge_id]*fabs(geometry.voronoi_cell_list[k].edges[kk].n2) +
			    force.force_list2[edge_id]*fabs(geometry.voronoi_cell_list[k].edges[kk].n2);
			/* sigma11 = sigma11 + force.force_list[edge_id]*cos(theta)*cos(theta); */
			/* sigma12 = sigma12 + force.force_list[edge_id]*cos(theta)*sin(theta); */
			/* sigma21 = sigma21 + force.force_list[edge_id]*sin(theta)*cos(theta); */
			/* sigma22 = sigma22 + force.force_list[edge_id]*sin(theta)*sin(theta); */
		    }
		    else if ( geometry.voronoi_cell_list[k].ellipse.type == 20 )
		    {
			boundary_force2_y = boundary_force2_y +
			    force.force_list[edge_id]*fabs(geometry.voronoi_cell_list[k].edges[kk].n2) +
			    force.force_list2[edge_id]*fabs(geometry.voronoi_cell_list[k].edges[kk].n2);
		    }
		    else
		    {
			external_force_y = external_force_y +
			    //force.force_list[edge_id]*fabs(geometry.voronoi_cell_list[k].edges[kk].n2);
			    force.force_list[edge_id]*fabs(geometry.voronoi_cell_list[k].edges[kk].n2) +
			    force.force_list2[edge_id]*fabs(geometry.voronoi_cell_list[k].edges[kk].n2);
		    }
		}
		// boundary force
		double F;
		if (params.feedback == 0)
		{
		    F = params.stretching_force;
		}
		else if (params.feedback == 1)
		{
		    F = params.stretching_force - params.alpha*sig/(1.0 + fabs(sig));
		    if ( F < 0.0 )
		    {
			F = 0.0;
		    }
		}
		else
		{
		    printf("feedback value can only be 1 and 0...exit!\n");
		    getchar();
		}
		if ( geometry.voronoi_cell_list[k].ellipse.type == 10 )
		{
		    sum_force3_x = -F -
			2.0*params.eta*(geometry.voronoi_cell_list[k].ellipse.c1 - L_min_old);
		    L_min = L_min + geometry.voronoi_cell_list[k].ellipse.c1;
		    L_min_counter = L_min_counter + 1;
		}
		else if ( geometry.voronoi_cell_list[k].ellipse.type == 20 )
		{
		    sum_force3_x = F -
			2.0*params.eta*(geometry.voronoi_cell_list[k].ellipse.c1 - L_max_old);
		    L_max = L_max + geometry.voronoi_cell_list[k].ellipse.c1;
		    L_max_counter = L_max_counter + 1;
		}
		else if ( geometry.voronoi_cell_list[k].ellipse.type < 10 )
		{
		    double dist_right = L_max_old - geometry.voronoi_cell_list[k].ellipse.c1;
		    double dist_left = geometry.voronoi_cell_list[k].ellipse.c1 - L_min_old;
		    if ( dist_right < 1e-4 )
		    {
			std::cout << "cell pertruding right boundary ... pause \n";
			getchar();
		    }
		    else if ( dist_right < 0.5 )
		    {
			sum_force3_x = -1.0*params.eta/dist_right;
		    }
		    else
		    {
			// do nothing
		    }
		    if ( dist_left < 1e-4 )
		    {
			std::cout << "cell pertruding left boundary ... pause \n";
			getchar();
		    }
		    if ( dist_left < 0.25 )
		    {
			sum_force3_x = 1.0*params.eta/dist_left;
		    }
		    else
		    {
			// do nothing
		    }
		}
		else
		{
		    std::cout << "in force calculation... wrong option.\n";
		}
		
		geometry.voronoi_cell_list[k].ellipse.c1 = 
		    geometry.voronoi_cell_list[k].ellipse.c1 +
		    (sum_force_x + sum_force2_x + sum_force3_x)*params.dt;
		geometry.voronoi_cell_list[k].ellipse.c2 = 
		    geometry.voronoi_cell_list[k].ellipse.c2 +
		    (sum_force_y + sum_force2_y + sum_force3_x)*params.dt;

		/* if ( geometry.voronoi_cell_list[k].ellipse.type == 1 ) */
		/* { */
		/*     double v1, v2; */
		/*     double u1, u2; */
		/*     double new_v1, new_v2; */
		/*     double new_u1, new_u2; */
		/*     double a, b, new_a, new_b; */
		/*     double aspect_ratio, new_aspect_ratio; */
		/*     double F_signal, e0; */
		
		/*     v1 = geometry.voronoi_cell_list[k].ellipse.v1; */
		/*     v2 = geometry.voronoi_cell_list[k].ellipse.v2; */
		/*     u1 = geometry.voronoi_cell_list[k].ellipse.u1; */
		/*     u2 = geometry.voronoi_cell_list[k].ellipse.u2; */

		/*     //compression */
		/*     a = sqrt(v1*v1 + v2*v2); */
		/*     b = sqrt(u1*u1 + u2*u2); */
		/*     aspect_ratio = a/b; */
		/*     F_signal = params.stretching_force*0.1; */
		/*     e0 = 1.0; */
		/*     new_aspect_ratio = aspect_ratio + */
		/* 	params.dt*(F_signal * fabs(v1)/sqrt(v1*v1 + v2*v2) - */
		/* 		   (aspect_ratio - e0)); */
		/*     if ( new_aspect_ratio < 1 ) */
		/*     { */
		/*     	new_aspect_ratio = 1.0; */
		/*     } */
		/*     new_a = sqrt(a*b*new_aspect_ratio); */
		/*     new_b = sqrt(a*b/new_aspect_ratio); */
		/*     new_v1 = v1*new_a/sqrt(v1*v1 + v2*v2); */
		/*     new_v2 = v2*new_a/sqrt(v1*v1 + v2*v2); */
		/*     new_u1 = u1*new_b/sqrt(u1*u1 + u2*u2); */
		/*     new_u2 = u2*new_b/sqrt(u1*u1 + u2*u2); */

		/*     /\* std::cout << k << " " << new_a << " " << new_aspect_ratio << "\n"; *\/ */
		/*     //rotation */
		/*     double d_phi = 25.0*torque*params.dt + 0.5*rnorm(gen); */
		/*     v1 = new_v1; */
		/*     v2 = new_v2; */
		/*     u1 = new_u1; */
		/*     u2 = new_u2; */
		/*     new_v1 = v1*cos(d_phi) - v2*sin(d_phi); */
		/*     new_v2 = v1*sin(d_phi) + v2*cos(d_phi); */
		/*     new_u1 = u1*cos(d_phi) - u2*sin(d_phi); */
		/*     new_u2 = u1*sin(d_phi) + u2*cos(d_phi); */

		/*     geometry.voronoi_cell_list[k].ellipse.v1 = new_v1; */
		/*     geometry.voronoi_cell_list[k].ellipse.v2 = new_v2; */
		/*     geometry.voronoi_cell_list[k].ellipse.u1 = new_u1; */
		/*     geometry.voronoi_cell_list[k].ellipse.u2 = new_u2; */
		/* } */
		
		if ( geometry.voronoi_cell_list[k].ellipse.type == 1 )
		{
		    double d_phi = 1.0*torque*params.dt;// + 0.1*rnorm(gen);
		    //double d_phi = 0.0;
		    double v1, v2;
		    double u1, u2;
		    v1 = geometry.voronoi_cell_list[k].ellipse.v1;
		    v2 = geometry.voronoi_cell_list[k].ellipse.v2;
		    u1 = geometry.voronoi_cell_list[k].ellipse.u1;
		    u2 = geometry.voronoi_cell_list[k].ellipse.u2;
		    geometry.voronoi_cell_list[k].ellipse.v1 = v1*cos(d_phi) - v2*sin(d_phi);
		    geometry.voronoi_cell_list[k].ellipse.v2 = v1*sin(d_phi) + v2*cos(d_phi);
		    geometry.voronoi_cell_list[k].ellipse.u1 = u1*cos(d_phi) - u2*sin(d_phi);
		    geometry.voronoi_cell_list[k].ellipse.u2 = u1*sin(d_phi) + u2*cos(d_phi);
		}
	    }
	}
	L_min = L_min/L_min_counter;
	L_max = L_max/L_max_counter;
	delta_H = -params.H_rate*external_force_y/(L_max - L_min);
	/* std::cout << "test: " << L_min << ", " << L_max << ", " */
	/* 	  << geometry.boundary_data.H << " external_force: " */
	/* 	  << external_force_y */
	/* 	  << " delta_H: " << delta_H << "\n"; */
	geometry.boundary_data.set_boundary_data(L_min, L_max,
						 geometry.boundary_data.H + delta_H);

	// extend the boundary wall if needed
	if ( boundary_force1_y/(2.0*geometry.boundary_data.H) > 0.1 )
	{
	    boundary_cell_division(geometry, force, 10);
	}
	if ( boundary_force2_y/(2.0*geometry.boundary_data.H) > 0.1 )
	{
	    boundary_cell_division(geometry, force, 20);
	}
	//foutput.close();
	return 0;
    }

};

#endif //ENGINE
