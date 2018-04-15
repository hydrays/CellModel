#ifndef PARAMETERS
#define PARAMETERS

#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <vector>

class Parameters
{
public:
    std::string output_dir;
    std::string cell_position_file;
    int Nstep;
    double dt;
    double eta;    
    int two_population_model;
    double q0;
    double a, b;
    double a_braf, b_braf;
    double braf_mosaic_percentage;
    double mu;
    double stretching_force;
    double H_rate;
    int flag_output_cell;
    int flag_output_force;
    int flag_record_final_state;    
	
public:
    bool init()
    {
	std::cout << "initializing parameters" << std::endl;

	// set parameters
	boost::property_tree::ptree pTree;

	try {
	    read_xml("config.xml", pTree);
	    std::cout << "reading config file: " << "config.xml" << std::endl;
	}
	catch (boost::property_tree::xml_parser_error e) {
	    std::cout << "error: no config file found." << std::endl;
	    getchar();
	}

	try {
	    output_dir = pTree.get<std::string>("main.output_dir");
	    std::cout << "output_dir: " << output_dir << std::endl;

	    cell_position_file = pTree.get<std::string>("main.cell_position_file");
	    std::cout << "cell_position_file: " << cell_position_file << std::endl;
	    
	    Nstep = pTree.get<int>("main.Nstep");
	    std::cout << "Nstep: " << Nstep << std::endl;

	    dt = pTree.get<double>("main.dt");
	    std::cout << "dt: " << dt << std::endl;

	    eta = pTree.get<double>("main.eta");
	    std::cout << "eta: " << eta << std::endl;
	    
	    two_population_model = pTree.get<int>("main.two_population_model");
	    std::cout << "two_population_model: " << two_population_model << std::endl;

	    q0 = pTree.get<double>("main.q0");
	    std::cout << "q0: " << q0 << std::endl;

	    a = pTree.get<double>("main.a");
	    std::cout << "a: " << a << std::endl;

	    b = pTree.get<double>("main.b");
	    std::cout << "b: " << b << std::endl;

	    a_braf = pTree.get<double>("main.a_braf");
	    std::cout << "a_braf: " << a_braf << std::endl;

	    b_braf = pTree.get<double>("main.b_braf");
	    std::cout << "b_braf: " << b_braf << std::endl;

	    braf_mosaic_percentage = pTree.get<double>("main.braf_mosaic_percentage");
	    std::cout << "braf_mosaic_percentage: " << braf_mosaic_percentage << std::endl;
		
	    mu = pTree.get<double>("main.mu");
	    std::cout << "mu: " << mu << std::endl;
	    
	    stretching_force = pTree.get<double>("main.stretching_force");
	    std::cout << "stretching_force: " << stretching_force << std::endl;

	    H_rate = pTree.get<double>("main.H_rate");
	    std::cout << "H_rate: " << H_rate << std::endl;

	    flag_output_cell = pTree.get<int>("main.flag_output_cell");
	    std::cout << "flag_output_cell: " << flag_output_cell << std::endl;

	    flag_output_force = pTree.get<int>("main.flag_output_force");
	    std::cout << "flag_output_force: " << flag_output_force << std::endl;
	    
	    flag_record_final_state = pTree.get<int>("main.flag_record_final_state");
	    std::cout << "flag_record_final_state: " << flag_record_final_state << std::endl;
	    
	    std::cout << "Done parameter initialize: Success!" << std::endl;	    
	}
	catch(boost::property_tree::ptree_bad_path e) {
	    std::cout << "error: parameter not set correctly." << std::endl;
	    getchar();
	}
    }
    
    ~Parameters()
    {
    }    
};

#endif //PARAMETERS
