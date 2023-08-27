#include "preprocess.hpp"

Preprocess::Preprocess() noexcept
{}

Preprocess::~Preprocess() noexcept
{}

// Load gmsh file
Preprocess::Mesh Preprocess::read_gmsh(const std::string& file_name) const
{
	Mesh mesh;
	std::vector<Types> typesList = {
		{2, 1, "LINES"},
		{3, 2, "TRIANGLES"},
		{ 4,  2, "QUADS"},
		{ 4,  3, "TETS"},
		{ 8,  3, "HEXAS"},
		{ 6,  3, "PRISMS"},
		{ 5,  3, "PYRAMIDS"},
		{ 3,  1, "LINES3"},
		{ 6,  2, "TRIANGLES6"},
		{ 9,  2, "QUADS9"},
		{ 10,  3, "TETS10"},
		{ 27,  3, "HEXAS27"},
		{ 18,  3, "PRISMS18"},
		{ 14,  3, "PYRAMIDS14"},
		{ 1,  0, "POINTS"},
		{ 8,  3, "QUADS8"},
		{ 20,  3, "HEXAS20"},
		{ 15,  3, "PRISMS15"},
		{ 13,  3, "PYRAMIDS13"},
	};
	for (const auto& types : typesList) {
		mesh.Types.emplace_back(types);
	}

	std::ifstream meshFile(file_name);

	if (!meshFile.is_open())
	{
		std::cerr << "Error: Cannot open the file, check file exist or not." 
		<< std::endl;
		return mesh;
	}

	int num_phy_grp;
	std::string sentence;
	while (std::getline(meshFile, sentence))
	{
		if (sentence == "$PhysicalNames")
		{
			getline(meshFile, sentence);
			std::istringstream iss(sentence);
			iss >> num_phy_grp;
			break;
		}
	}

	while (std::getline(meshFile, sentence))
	{
		if (sentence == "$EndPhysicalNames") break;
		int degree_grp;
		int id_grp;
		std::string region;
		std::istringstream iss(sentence);
		iss >> degree_grp >> id_grp >> region;
		mesh.PhysicalNames.emplace_back(num_phy_grp, degree_grp, id_grp, region);
	}

	while (std::getline(meshFile, sentence))
	{
		if (sentence == "$Nodes")
		{
			getline(meshFile, sentence);
			std::istringstream iss(sentence);
			int num_nodes;
			if (iss >> num_nodes)
			{
				mesh.num_nodes = num_nodes;
				std::cout << "number of nodes: " << num_nodes << std::endl;
				break;
			}
		}
	}

	while (std::getline(meshFile, sentence))
	{
		std::istringstream iss(sentence);
		int id;
		std::vector<double> coor(3);
		if (iss >> id >> coor[0] >> coor[1] >> coor[2])
		{
			mesh.Nodes.emplace_back(id, coor);
		}
		if (sentence == "$EndNodes") break;
	}

	while (std::getline(meshFile, sentence))
	{
		if (sentence == "$Elements")
		{
			getline(meshFile, sentence);
			std::istringstream iss(sentence);
			int num_elem;
			if (iss >> num_elem)
			{
				mesh.num_elem = num_elem;
				std::cout << "number of elements: " << num_elem << std::endl;
				break;
			}
		}
	}
	while (std::getline(meshFile, sentence)) 
	{
		std::istringstream iss(sentence);
		int info[8] = {0, 0, 0, 0, 0, 0, 0, 0};
		if (sentence == "$EndElements") break;
		if (iss >> info[0] >> info[1] >> info[2] >> info[3] >> info[4] >> info[5] >> info[6] >> info[7])
		{
			std::vector<int> node_square_elem = {info[5], info[6], info[7]};
			mesh.Squares.emplace_back(info[0], info[1], info[2], info[3], info[4], node_square_elem);
		}
		else 
		{
			iss >> info[0] >> info[1] >> info[2] >> info[3] >> info[4] >> info[5] >> info[6];
			std::vector<int> node_line_elem = {info[5], info[6]};
			mesh.Lines.emplace_back(info[0], info[1], info[2], info[3], info[4], node_line_elem);

		}
	}

	meshFile.close();

	return mesh;
}

// Load info
Preprocess::Info Preprocess::read_info(const std::string& file_name) const
{
	YAML::Node input = YAML::LoadFile(file_name);
	Info info( input["Problem Type"]["Type"].as<int>(),
			   input["Number of Direction"].as<int>(),
			   input["Quadrature Degree"]["Degree"].as<int>(),
			   input["Element Degree"]["Degree"].as<int>(),
			   input["Material Property"]["YoungsModulus"].as<double>(),
			   input["Material Property"]["PoissonsRatio"].as<double>(),
			   input["Dirichlet Boundary"]["Count"].as<int>(),
			   input["Neumann Boundary"]["Count"].as<int>() );
	return info;
}

std::vector<std::vector<int>> Preprocess::get_ID(const Preprocess::Mesh &mesh, const Preprocess::Info &info) const
{
	std::vector<std::vector<int>> mat( info.num_dir, std::vector<int>(mesh.num_nodes, 1) );

    for (const Lines &line : mesh.Lines)
    {
    	if (line.bc_tag == 2)
    	{
    		mat[0][line.node_line_elem[0]] = 0;
    		mat[0][line.node_line_elem[1]] = 0;
    	}
    	if (line.bc_tag == 3)
    	{
    		mat[1][line.node_line_elem[0]] = 0;
    		mat[1][line.node_line_elem[1]] = 0;
    	}
    }
    for (const Squares &square : mesh.Squares)
    {
    	if (square.bc_tag == 2)
    	{
    		mat[0][square.node_square_elem[0]] = 0;
    		mat[0][square.node_square_elem[1]] = 0;
    		mat[0][square.node_square_elem[2]] = 0;    	
    	}
    	if (square.bc_tag == 3)
    	{
    		mat[1][square.node_square_elem[0]] = 0;
    		mat[1][square.node_square_elem[1]] = 0;
    		mat[1][square.node_square_elem[2]] = 0;    	
    	}
    }
    int counter = 0;
    for (int ii=0; ii<mesh.num_nodes; ++ii) 
    {
    	for (int jj=0; jj<info.num_dir; ++jj) 
    	{
    		if (mat[jj][ii] != 0)
    		{
    			counter += 1;
    			mat[jj][ii] = counter;
    		}
    	}
    }
    print(mat);
    return mat;
}

std::vector<std::vector<int>> Preprocess::get_IEN(const Preprocess::Mesh &mesh, const Preprocess::Info &info) const
{
	std::vector<std::vector<int>> mat( mesh.Squares[0].node_square_elem.size(),
		std::vector<int>(mesh.Lines.size()+mesh.Squares.size(), 0) );

	for (const Lines &line : mesh.Lines)
	{
		mat[0][line.id_line_elem-1] = line.node_line_elem[0];
		mat[1][line.id_line_elem-1] = line.node_line_elem[1];
	}

    for (const Squares &square : mesh.Squares)
	{
		mat[0][square.id_square_elem-1] = square.node_square_elem[0];
		mat[1][square.id_square_elem-1] = square.node_square_elem[1];
		mat[2][square.id_square_elem-1] = square.node_square_elem[2];
	}

	print(mat);
	return mat;
}

std::vector<std::vector<int>> Preprocess::get_IEN(const std::vector<std::vector<int>> &ID, const std::vector<std::vector<int>> &IEN) const
{
	std::vector<std::vector<int>> mat(2*IEN.size(), std::vector<int>(IEN[0].size(), 0));
	for (int ii=0; ii<IEN.size(); ++ii) 
	{
		for (int jj=0; jj<IEN[0].size(); ++jj) 
		{
			mat[ii][jj] = ID[0][IEN[ii][jj]];
		}
	}

	for (int ii=IEN.size(); ii<2*IEN.size(); ++ii) 
	{
		for (int jj=0; jj<IEN[0].size(); ++jj) 
		{
			mat[ii][jj] = ID[1][IEN[ii-IEN.size()][jj]];
		}
	}

	print(mat);
	return mat;
}

Preprocess::Tri_quad Preprocess::TriQuad(const int &degree) const
{
	if (degree == 2)
	{
		Tri_quad result(
			{ {0.666666666666667, 0.166666666666667, 0.166666666666667},
			{0.166666666666667, 0.666666666666667, 0.166666666666667},
			{0.166666666666667, 0.166666666666667, 0.666666666666667} }, 
			{0.333333333333333, 0.333333333333333, 0.333333333333333},
			3 );
		for (double &w : result.wq) 
		{
			w *= 0.5;
		}
		return result;
	}
	else if (degree == 3) 
	{
		Tri_quad result(
			{ {0.333333333333333, 0.6, 0.2, 0.2},
			{0.333333333333333, 0.2, 0.6, 0.2},
			{0.333333333333333, 0.2, 0.2, 0.6} },
			{0.520833333333333, 0.520833333333333, 0.520833333333333, -0.56250},			
			4 );
		for (double &w : result.wq) 
		{
			w *= 0.5;
		}
		return result;
	} 
	else if (degree == 4)
	{
		Tri_quad result(
			{ {0.816847572980459, 0.091576213509771, 0.091576213509771, 0.091576213509771, 0.108103018168070, 0.445948490915965},
			{0.091576213509771, 0.816847572980459, 0.091576213509771, 0.445948490915965, 0.108103018168070, 0.091576213509771},
			{0.091576213509771, 0.091576213509771, 0.816847572980459, 0.445948490915965, 0.445948490915965, 0.108103018168070} },
			{0.109951743655322, 0.109951743655322, 0.109951743655322, 0.223381589678011, 0.223381589678011, 0.223381589678011},
			6 );
		    // Scale the weights for triangles in a plane
		for (double &w : result.wq) 
		{
			w *= 0.5;
		}
		return result;
	}
	else 
	{
		Tri_quad result;
		std::cout << "TriangularQuad: Please check the data input." << std::endl;
		std::cout << "The degree input should be an integer from 2 to 7." << std::endl;
		return result;
	}
}

Preprocess::Line_quad Preprocess::Gauss(const int &NN, const double &a, const double &b) const
{
	int N = NN - 1;
	const int N1 = N + 1;
	const int N2 = N + 2;
	double xu[N1];
	for (int ii = 0; ii < N1; ii++)
	{
		xu[ii] = -1.0 + double(ii) * 2.0/(double(N1)-1.0);
	}
	double y[N1];
	for (int ii = 0; ii < N1; ii++)
	{
		y[ii] = cos((2.0*double(ii)+1)*pi/(2.0*double(N)+2.0))+(0.27/double(N1))*sin(pi*xu[ii]*double(N)/double(N2));
	} 
	double L[N1][N2];
	double * Lp = new double[N1];
	double y0[N1];
	double error = 1.0;

	while (error > eps)
	{
		for (int ii = 0; ii < N1; ii++)
		{
			L[ii][0] = 1.0;
			L[ii][1] = y[ii]; 
		}
		for (int ii = 1; ii < N1; ii++)
		{
			for (int jj = 0; jj < N1; jj++)
			{
				L[jj][ii+1] = ( (2.0*double(ii)+1.0) * y[jj] * L[jj][ii] - double(ii)*L[jj][ii-1] ) / double(ii+1);
			}
		}
		for (int ii = 0; ii < N1; ii++)
		{
			Lp[ii] = double(N2) * (L[ii][N1-1] - y[ii]*L[ii][N2-1] ) / (1.0 - y[ii]*y[ii]);
		}
		for (int ii = 0; ii < N1; ii++)
		{
			y0[ii] = y[ii];
			y[ii]  =y0[ii] - L[ii][N2-1] / Lp[ii];
		}

		double error0 = 0.0;
		for (int ii = 0; ii < N1; ii++)
		{
			error = (error0 > abs( y[ii]-y0[ii] )) ? error : abs(y[ii]-y0[ii]);
		}
		error0 = error;
	}

	std::vector<double> x(N1);
	std::vector<double> w(N1);

	for (int ii = 0; ii < N1; ++ii) {
		x[ii] = (a*(1.0-y[ii])+b*(1.0+y[ii])) / 2.0;
		w[ii] = (b-a) / ((1-y[ii]*y[ii]) * Lp[ii]*Lp[ii]) * (double(N2)/double(N1)) * (double(N2)/double(N1));
	}

    Line_quad result(x, w);

    return result;
}



void Preprocess::print(const std::vector<std::vector<int>> &mat) const
{
	std::cout << std::endl;
	for (const std::vector<int> &row : mat) 
	{
		for (int value : row) 
		{
            std::cout << std::setw(2) << value << " "; // Adjust the setw value as needed
        }
        std::cout << std::endl;
    }
}
