#ifndef PREPROCESS_HPP
#define PREPROCESS_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#define pi acos(-1)
#define eps 2.2204e-16
#include "/Users/chongran/lib/yaml-shared/include/yaml-cpp/yaml.h"

class Preprocess
{
public:
	struct Types
	{
		int num_nodes;
		int dimension;
		std::string name;
		Types(int num, int dim, const std::string& type_name) :
		num_nodes(num), dimension(dim), name(type_name)
		{
			std::cout << "type: " << name << std::setw(8)<< "\t" << "dimension: " 
			<< dimension << "\t" << "number of nodes: " << num_nodes << std::endl;
		}
	};

	struct Nodes
	{
		int ID_nodes;
		std::vector<double> coor_nodes;

		Nodes(const int id, std::vector<double> coor) :
		ID_nodes(id), coor_nodes(coor)
		{
			std::cout << "Nodes ID:" << ID_nodes << std::setw(8) << "\t" << "coordinates: ( "
			<< std::setw(8) << coor_nodes[0] << ", " << std::setw(8) << coor_nodes[1]
			<< ", " << std::setw(8) << coor_nodes[2] << " )" << std::endl;
		}
	};

	struct Lines
	{
		int id_line_elem;
		int phy_grp; // physical group of line is 1
		int type_elem;
		int bc_tag;
		int id_line;
		std::vector<int> node_line_elem;

		Lines(int id, int grp, int type, int tag, int line, std::vector<int> node) : 
		id_line_elem(id),  phy_grp(grp), type_elem(type), bc_tag(tag), id_line(line),
		node_line_elem(node)
		{
			std::cout << "\nNo. " << id << std::setw(2) << " line element"
			<< " physical group: " << grp << std::setw(2)
			<< " line element type: " << type << std::setw(2)
			<< " boundary condition type: " << tag << std::setw(2)
			<< " at No. " << line << " line, using node "
			<< std::setw(2) << node[0] << " and node " << std::setw(2)
			<< node[1] << std::endl; 
		}
	};

	struct Squares
	{
		int id_square_elem;
		int phy_grp; // physical group of square is 2
		int type_elem;
		int bc_tag;
		int id_square;
		std::vector<int> node_square_elem;

		Squares(int id, int grp, int type, int tag, int square, std::vector<int> node) : 
		id_square_elem(id),  phy_grp(grp), type_elem(type), bc_tag(tag), id_square(square),
		node_square_elem(node)
		{
			std::cout << "\nNo. " << id << std::setw(2) << " square element"
			<< " physical group: " << grp << std::setw(2)
			<< " square element type: " << type << std::setw(2)
			<< " boundary condition type: " << tag << std::setw(2)
			<< " at No. " << square << " square, using node "
			<< std::setw(2) << node[0] << " , node " << std::setw(2)
			<< node[1] << " and node " << std::setw(2) << node[2]
			<< std::endl; 
		}
	};

	struct PhysicalNames
	{
		int num_phy_grp;
		int degree_grp;
		int id_grp;
		std::string region;

		PhysicalNames(int num, int degree, int grp, std::string reg) :
		num_phy_grp(num), degree_grp(degree), id_grp(grp), region(reg)
		{
			std::cout << "There are " << num << " physical groups."
			<< " This is No. " << id_grp << ". It's degree "
			<< degree << ", with BC: " << reg << std::endl;
		}

	};

	struct Mesh
	{
		std::vector<Types> Types;
		int num_nodes;
		std::vector<Nodes> Nodes;
		int num_elem;
		std::vector<Lines> Lines;
		std::vector<Squares> Squares;
		std::vector<PhysicalNames> PhysicalNames;
	};

	struct Info
	{
		int prob_type;
		int num_dir;
		int quad_deg;
		int elem_deg;
		double EE; // Young's modulus
		double nu; // Poisson's modulus 
		int num_DB; // number of Dirichlet boundary
		int num_NB;	// number of Neumann boundary

		Info(int pt, int ndir, int qd, int ed, double ee, double nn, int nD, int nN) :
		prob_type(pt), num_dir(ndir), quad_deg(qd), elem_deg(ed), EE(ee), nu(nn),
		num_DB(nD), num_NB(nN)
		{
			if (pt == 0)
			{
				std::cout << "=====Plane strain problem=====" << std::endl;
			}
			else if (pt == 1)
			{
				std::cout << "=====Plane stress problem====" << std::endl;
			}
			std::cout << "Number of direction: " << ndir << "\n";
			std::cout << "Quadrature degree: " << qd << "\n";
			std::cout << "Element degree: " << ed << "\n";
			std::cout << "Young's modulus: " << ee << "\n";
			std::cout << "Poisson's modulus: " << nn << "\n" << std::endl;
		}
	};

	struct Tri_quad
	{
    	std::vector<std::vector<double>> qp; // Quadrature points
    	std::vector<double> wq; // Quadrature weights
    	int nqp; // Number of quadrature points

    	Tri_quad() : nqp(0) {}

    	Tri_quad(std::vector<std::vector<double>> aa, std::vector<double> bb, int cc) :
    	qp(aa), wq(bb), nqp(cc)
    	{
    		std::cout << "\nUse " << cc << " quadrature nodes in triangle elements" << "\n";
    		std::cout << "weight\t\t coordinates\n";
    		for (int ii=0; ii<bb.size(); ++ii)
    		{
    			std::cout << std::setprecision(4) << std::fixed;
    			std::cout << std::setw(6) << bb[ii] << "\t" << std::setw(15) << aa[0][ii] << " ";
    			std::cout<< std::setw(6) << aa[1][ii] << " " << std::setw(6) << aa[2][ii] << std::endl;
    		}
    	}
    };

    struct Line_quad
	{
    	std::vector<double> qp; // Quadrature points
    	std::vector<double> wq; // Quadrature weights

    	Line_quad() {}

    	Line_quad(std::vector<double> aa, std::vector<double> bb) :
    	qp(aa), wq(bb)
    	{
    		std::cout << "\nUse " << aa.size() << " quadrature nodes in line elements" << "\n";
    		std::cout << "weight\t\t coordinates\n";
    		for (int ii=0; ii<aa.size(); ++ii)
    		{
    			std::cout << std::setprecision(4) << std::fixed;
    			std::cout << std::setw(6) << bb[ii] << "\t" << std::setw(15) << aa[ii] <<std::endl;
    		}
    	}
    };

	Preprocess() noexcept;
	~Preprocess() noexcept;
	Mesh read_gmsh(const std::string& file_name) const;
	Info read_info(const std::string& file_name) const;
	std::vector<std::vector<int>> get_ID(const Preprocess::Mesh &mesh, const Preprocess::Info &info) const;
	std::vector<std::vector<int>> get_IEN(const Preprocess::Mesh &mesh, const Preprocess::Info &info) const;
	std::vector<std::vector<int>> get_LM(const std::vector<std::vector<int>> &ID, const std::vector<std::vector<int>> &IEN) const;
	void print(const std::vector<std::vector<int>> &mat) const;
	Tri_quad TriQuad(const int &degree) const;
	Line_quad Gauss(const int &NN, const double &a, const double &b) const;
};

#endif