#ifndef DEF_GENERATOR
#define DEF_GENERATOR
#pragma once

#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>


#include <iostream>
#include <vector>
#include <map>
#include "Eigen/Dense"
#include "field.h"

typedef Eigen::SparseMatrix<float, Eigen::ColMajor, long> SpMat;


class Generator
{
protected:
	int m_size;
public:
	int size()const;
	virtual std::vector<float> operator()(const std::vector<float*>& data, int ix, int iy, int iz)const=0;
};


class GrayScottGenerator : public Generator
{
private:
	Kernel3D<float> m_laplacian;
	float m_diff_a;
	float m_diff_b;
	float m_kill;
	Field3D<float> m_feed;
public:
	GrayScottGenerator(float t_diff_a, float t_diff_b, float t_kill, const Field3D<float>& t_feed);
	~GrayScottGenerator();
	float diffA()const;
	float diffB()const;
	float kill()const;
	Field3D<float> feed()const;
	void setKernel(const Kernel3D<float>& kernel);
	Kernel3D<float> getKernel()const;
	std::vector<float> convolve(const std::vector<float*>& data, int ix, int iy, int iz)const;
	virtual std::vector<float> operator()(const std::vector<float*>& data, int ix, int iy, int iz)const;
};



class SemiGroup
{
protected:
	Generator* m_generator;
public:
	Generator* getGenerator();
	virtual void operate(const std::vector<float*>& former_state, std::vector<float*>& new_state)const=0;
};

class GrayScottSG : public SemiGroup
{
private:
	bool m_is_explicit;
	float m_diff_a;
	float m_diff_b;
	float m_kill;
	float m_delta_t;
	Field3D<float> m_feed;
	SpMat m_left_matrix_a;
	SpMat m_left_matrix_b;
	Eigen::ConjugateGradient<SpMat, Eigen::Lower | Eigen::Upper> m_solver_a;
	Eigen::ConjugateGradient<SpMat, Eigen::Lower | Eigen::Upper> m_solver_b;
	int m_max_iter;
public:
	GrayScottSG(float t_diff_a, float t_diff_b, float t_kill, const Field3D<float>& t_feed, bool t_is_explicit);
	~GrayScottSG();
	void setTimeStep(float delta_t);
	void setMaxIterations(int t_max_iter);
	void setKernel(const Kernel3D<float>& kernel);
	Kernel3D<float> getKernel()const;
	std::vector<Eigen::VectorXf> getRightMembers(const std::vector<Eigen::VectorXf>& former_state_vecs)const;
	void buildMatrices();
	void explicitScheme(const std::vector<float*>& former_state, std::vector<float*>& new_state)const;
	void crankNicolsonScheme(const std::vector<float*>& former_state, std::vector<float*>& new_state)const;
	void operate(const std::vector<float*>& former_state, std::vector<float*>& new_state)const;
	
};



SpMat laplacian3DMatrix(int size, const Kernel3D<float>& kernel);

#endif