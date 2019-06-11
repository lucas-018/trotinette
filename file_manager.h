#ifndef DEF_FILE_MANAGER
#define DEF_FILE_MANAGER
#pragma once

#include <iostream>
#include "field.h"
#include <string>
#include <fstream>


std::vector<std::string> split(std::string line, char separator);
std::vector<int> getNumbersInteger(std::string line, char separator);
std::vector<float> getNumbersFloat(std::string line, char separator);


template<class T>
T getFromLine(std::string line, char separator)
{
	T res(0);
	return res;
}

template<>
Eigen::Vector3d getFromLine<Eigen::Vector3d>(std::string line, char separator);

template<>
float getFromLine<float>(std::string line, char separator);

template<class T>
std::vector<T> getData(std::string file_name, int& size_x, int& size_y, int& size_z, int& dim)
{
	std::vector<T> unfold_field;
	std::ifstream file(file_name, ios::in);
	if (file)
	{
		std::string first_line;
		std::getline(file, first_line);
		std::vector<int> dimensions = getNumbersInteger(first_line, ' ');
		size_x = dimensions[0];
		size_y = dimensions[1];
		size_z = dimensions[2];
		dim = dimensions[3];
		std::string current_line;
		while (std::getline(file, current_line))
		{
			unfold_field.push_back(getFromLine<T>(current_line, ' '));
		}
		file.close();
	}
	else
	{
		cerr << "Could not open " << file_name << "..." << endl;
	}
	return unfold_field;
}

template<class T>
Field3D<T> getField(std::string file_name)
{
	int size_x = 0;
	int size_y = 0;
	int size_z = 0;
	int dim;
	std::vector<T> unfold_field = getData<T>(file_name, size_x, size_y, size_z, dim);
	return Field3D<T>::fromArray(unfold_field, size_x, size_y, size_z);
}





#endif

















