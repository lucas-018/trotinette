#ifndef DEF_GENERATOR
#define DEF_GENERATOR
#pragma once


#include <iostream>
#include <vector>
#include <map>
#include "Eigen/Dense"
#include "field.h"

class Generator
{
protected:
public:
	virtual std::vector<float> operator()(std::vector<float*> data, int ix, int iy, int iz)=0;
};


class GrayScottGenerator : public Generator
{
private:
	Kernel3D<float> m_laplacian;
	int m_size;
	float m_diff_a;
	float m_diff_b;
	float m_kill;
	Field3D<float> m_feed;
public:
	GrayScottGenerator(float t_diff_a, float t_diff_b, float t_kill, const Field3D<float>& t_feed);
	std::vector<float> convolve(std::vector<float*> data, int ix, int iy, int iz);
	virtual std::vector<float> operator()(std::vector<float*> data, int ix, int iy, int iz);
};




#endif