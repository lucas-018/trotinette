#include "file_manager.h"


using namespace std;
using namespace Eigen;


vector<string> split(string line, char separator)
{
	vector<string> table;
	string current_word = "";
	for (int i = 0; i < (int)line.size(); ++i)
	{
		if (line[i] == separator)
		{
			table.push_back(current_word);
			current_word = "";
		}
		else
		{
			current_word += line[i];
		}
	}
	table.push_back(current_word);
	return table;
}


vector<int> getNumbersInteger(string line, char separator)
{
	vector<string> str_table = split(line, separator);
	vector<int> d_table;
	for (int i = 0; i < (int)str_table.size(); ++i)
	{
		d_table.push_back(stoi(str_table[i]));
	}
	return d_table;
}

vector<float> getNumbersFloat(string line, char separator)
{
	vector<string> str_table = split(line, separator);
	vector<float> d_table;
	for (int i = 0; i < (int)str_table.size(); ++i)
	{
		d_table.push_back(stof(str_table[i]));
	}
	return d_table;
}


template<>
Vector3d getFromLine<Vector3d>(string line, char separator)
{
	vector<string> str_table = split(line, separator);
	Vector3d vec(str_table.size());
	for (int i = 0; i < (int)str_table.size(); ++i)
	{
		vec[i] = stod(str_table[i]);
	}
	return vec;
}

template<>
float getFromLine<float>(string line, char separator)
{
	vector<string> str_table = split(line, separator);
	return stof(str_table[0]);
}


