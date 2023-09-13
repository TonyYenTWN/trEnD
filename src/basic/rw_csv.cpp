// Read and write csv files
#include "rw_csv.h"
//#include <iostream> // just for checking

std::vector <int> basic::get_file_dim(std::string filename, bool row_name, bool col_name){
	std::ifstream in(filename);
	std::vector <int> dim;
	dim.reserve(2);

	if(in){
		std::string line;
		int row_ID = 0;
		int col_ID = 0;

		while(std::getline(in, line)){
			row_ID += 1;

			if(row_ID == 1){
		  	  	std::stringstream sep(line);
		  	  	std::string field;
		  	  	col_ID = 0;

		  	  	while(std::getline(sep, field, ',')){
		  	  		col_ID += 1;
			  	}

				dim.push_back(col_ID - row_name);
			}
		}

		dim.push_back(row_ID - col_name);
	}

	std::reverse(dim.begin(), dim.end());
	in.close();
	return dim;
}

std::vector <std::string> basic::get_col_name(std::string filename, int col_num){
	std::vector <std::string> col_names;
	col_names.reserve(col_num);

	// Read header row
	std::ifstream in(filename);
	std::string line;
	std::getline(in, line);

	// Seperate and read each column in the header row
	std::stringstream sep(line);
	std::string field;

	while(std::getline(sep, field, ',')){
		col_names.push_back(field);
	}

	// Close file
	in.close();
	return col_names;
}

std::vector <std::string> basic::get_row_name(std::string filename, int row_num){
	std::vector <std::string> row_names;
	row_names.reserve(row_num);

	std::ifstream in(filename);

	if(in){
	  	std::string line;
	  	std::getline(in, line); // skip the first line
	  	//int row_ID = 0;

	  	while(getline(in, line)){
	  	  	std::stringstream sep(line);
	  	  	std::string field;
	  	  	std::getline(sep, field, ',');

			row_names.push_back(field);
	  	}
	}

  	// Close file
	in.close();
	return row_names;
}

std::map<std::string, std::vector <std::string>> basic::read_config_file(std::string filename){
    std::ifstream in(filename);
    std::map<std::string, std::vector <std::string>> key_value;

    if(in){
	  	std::string line;

	  	// the first line
	  	// find number of values for each type
	  	int num_col = 0;
	  	std::getline(in, line);
	  	{
	  	  	std::stringstream sep(line);
	  	  	std::string field;

	  	  	// skip col for row name
            std::getline(sep, field, ',');

            while(std::getline(sep, field, ',')){
                num_col += 1;
            }
	  	}
//	  	std::cout << num_col << "\n";

	  	while(getline(in, line)){
	  	  	std::stringstream sep(line);
	  	  	std::string field;
	  	  	std::vector <std::string> value_vector;
            value_vector.reserve(num_col);

	  	  	// read row name
	  	  	std::getline(sep, field, ',');
            std::string key = field;
//            std::cout << field << ":\t";

            // read value names
            while(std::getline(sep, field, ',')){
                value_vector.push_back(field);
//                std::cout << field << "\t";
            }
//            std::cout << "\n";
            key_value.insert(std::pair<std::string, std::vector <std::string>>(key, value_vector));
	  	}
    }

    // Close file
	in.close();
	return key_value;
};

Eigen::MatrixXd basic::read_file(int num_row, int num_col, std::string filename, bool row_name){
	std::ifstream in(filename);
	Eigen::MatrixXd data(num_row, num_col - row_name);

	if(in){
	  	std::string line;
	  	std::getline(in, line); // skip the first line
	  	int row_ID = 0;
	  	int col_ID;

	  	while(getline(in, line)){
	  	  	std::stringstream sep(line);
	  	  	std::string field;
	  	  	col_ID = 0;
	  	  	if(row_name){
				std::getline(sep, field, ',');
				col_ID += 1;
	  	  	}

	  	  	while(std::getline(sep, field, ',')){
	  	  		data(row_ID, col_ID - row_name) = std::stod(field);
	  	  		col_ID += 1;
		  	}

	  	  row_ID += 1;
		}
  	}

  	// Close file
	in.close();
	return data;
}

void basic::write_file(Eigen::MatrixXd data, std::string filename, std::vector<std::string> col_name, int precision){

	int num_row = data.rows();
	int num_col = data.cols();

	std::fstream fout;

	// opening an existing csv file or creating a new csv file
    fout.open(filename, std::ios::out);
    fout << std::fixed << std::setprecision(precision);

    // Write the column names of the csv file
    for(int col_ID = 0; col_ID < num_col; ++ col_ID){
    	fout << col_name[col_ID];
    	if(col_ID < num_col - 1){
    		fout << " ";
		}
	}
	fout << "\n";

    // Write the content of the csv file
    for(int row_ID = 0; row_ID < num_row; ++ row_ID){
    	fout << data.row(row_ID);
    	fout << "\n";
	}

    fout.close();   // closing csv file
}
