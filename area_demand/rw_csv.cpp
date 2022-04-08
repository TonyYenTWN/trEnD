// Read and write csv files
//#pragma once

#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include "Basic_Definitions.h"

using namespace std;
using namespace Eigen;

MatrixXd read_file(int num_row, int num_col, string filename){
	ifstream in(filename);
	MatrixXd data(num_row, num_col);
	
	if(in){
  	string line;
  	getline(in, line); // skip the first line
  	int row_ID = 0;
  	int col_ID;
  	
  	while(getline(in, line)){
  	  	stringstream sep(line);
  	  	string field;
  	  	col_ID = 0; 
  	  
  	  	while(getline(sep, field, ',')){
  	  		data(row_ID, col_ID) = stod(field);
  	  		col_ID += 1;
	  	}
  	  
  	  row_ID += 1;
	}
  }
  
  return data;
}

void write_file(MatrixXd data, string filename, vector<string> col_name){
	
	int num_row = data.rows();
	int num_col = data.cols();
	
	fstream fout;
	
	// opening an existing csv file or creating a new csv file
    fout.open(filename, ios::out);
    fout << fixed << setprecision(20);
    
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