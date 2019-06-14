#include "generator.h"

using namespace std;
using namespace Eigen;


int Generator::size()const
{
	return m_size;
}


GrayScottGenerator::GrayScottGenerator(float t_diff_a, float t_diff_b, float t_kill, const Field3D<float>& t_feed)
	:m_diff_a(t_diff_a),
	m_diff_b(t_diff_b),
	m_kill(t_kill),
	m_feed(t_feed)
{
	m_size = t_feed.sizeX();
}


float GrayScottGenerator::diffA()const
{
	return m_diff_a;
}

float GrayScottGenerator::diffB()const
{
	return m_diff_b;
}

float GrayScottGenerator::kill()const
{
	return m_kill;
}

Field3D<float> GrayScottGenerator::feed()const
{
	return m_feed;
}

void GrayScottGenerator::setKernel(const Kernel3D<float>& kernel)
{
	m_laplacian = kernel;
}

Kernel3D<float> GrayScottGenerator::getKernel()const
{
	return m_laplacian;
}

vector<float> GrayScottGenerator::convolve(const vector<float*>& data, int ix, int iy, int iz)const
{
	vector<float> dd;
	for (int i = 0; i < 2; ++i)
	{
		dd.push_back(0);
	}
	int start = m_laplacian.startIndex();
	int stop = m_laplacian.stopIndex();
	for (int relx = start; relx < stop; ++relx)
	{
		for (int rely = start; rely < stop; ++rely)
		{
			for (int relz = start; relz < stop; ++relz)
			{
				for (int i = 0; i < 2; ++i)
				{
					dd[i] += *element_at(data[i], ix + relx, iy + rely, iz + relz, m_size)*m_laplacian.getRelValue(relx, rely, relz)*6;
				}
			}
		}
	}
	return dd;
}

vector<float> GrayScottGenerator::operator()(const vector<float*>& data, int ix, int iy, int iz)const
{

	float value_a = *element_at(data[0], ix, iy, iz, m_size);
	float value_b = *element_at(data[1], ix, iy, iz, m_size);
	vector<float> dd = convolve(data, ix, iy, iz);
	vector<float> variation;
	float da = m_diff_a * dd[0] - value_a*value_b*value_b + m_feed.getValue(ix, iy, iz) * (1.0f - value_a);
	float db = m_diff_b * dd[1] + value_a*value_b*value_b - (m_feed.getValue(ix, iy, iz) + m_kill)*value_b;
	variation.push_back(da);
	variation.push_back(db);
	return variation;
}



Generator* SemiGroup::getGenerator()
{
	return m_generator;
}


GrayScottSG::GrayScottSG(float t_diff_a, float t_diff_b, float t_kill, const Field3D<float>& t_feed, bool t_is_explicit)
	:m_diff_a(t_diff_a),
	m_diff_b(t_diff_b),
	m_kill(t_kill),
	m_feed(t_feed),
	m_delta_t(1.0f),
	m_is_explicit(t_is_explicit),
	m_max_iter(100)
{
	m_generator = new GrayScottGenerator(t_diff_a, t_diff_b, t_kill, t_feed);
	this->setKernel(laplacianKernel(0.6f));
}

GrayScottSG::~GrayScottSG()
{
	delete m_generator;
}


void GrayScottSG::setTimeStep(float delta_t)
{
	m_delta_t = delta_t;
}


void GrayScottSG::setMaxIterations(int t_max_iter)
{
	m_max_iter = t_max_iter;
}


void GrayScottSG::setKernel(const Kernel3D<float>& kernel)
{
	(static_cast<GrayScottGenerator*>(m_generator))->setKernel(kernel);
}


Kernel3D<float> GrayScottSG::getKernel()const
{
	return (static_cast<GrayScottGenerator*>(m_generator))->getKernel();
}


void GrayScottSG::buildMatrices()
{
	int size = m_generator->size();
	SparseMatrix<float> laplacian_matrix = laplacian3DMatrix(size, this->getKernel());

	m_left_matrix_a = -0.5*m_delta_t*m_diff_a*laplacian_matrix;
	m_left_matrix_b = -0.5*m_delta_t*m_diff_b*laplacian_matrix;

	m_right_matrix_a = 0.5*m_delta_t*m_diff_a*laplacian_matrix;
	m_right_matrix_b = 0.5*m_delta_t*m_diff_b*laplacian_matrix;

	for (int k = 0; k < size*size*size; ++k)
	{
		int x_pos = (int)(k / size / size);
		int y_pos = (int)((k - x_pos*size*size) / size);
		int z_pos = (int)(k - x_pos*size*size - y_pos*size);

		m_left_matrix_a.coeffRef(k, k) += 0.5f*m_delta_t*m_feed(x_pos, y_pos, z_pos) + 1.0f;
		m_left_matrix_b.coeffRef(k, k) += 0.5f*m_delta_t*(m_feed(x_pos, y_pos, z_pos) + m_kill) + 1.0f;

		m_right_matrix_a.coeffRef(k, k) += -0.5f*m_delta_t*m_feed(x_pos, y_pos, z_pos) + 1.0f;
		m_right_matrix_b.coeffRef(k, k) += -0.5f*m_delta_t*(m_feed(x_pos, y_pos, z_pos) + m_kill) + 1.0f;
	}

	m_solver_a.compute(m_left_matrix_a);
	m_solver_b.compute(m_left_matrix_b);

	m_solver_a.setMaxIterations(m_max_iter);
	m_solver_b.setMaxIterations(m_max_iter);
}


void GrayScottSG::explicitScheme(const vector<float*>& former_state, vector<float*>& new_state)const
{
	int size = m_generator->size();
	for (int ix = 0; ix < size; ++ix)
	{
		for (int iy = 0; iy < size; ++iy)
		{
			for (int iz = 0; iz < size; ++iz)
			{
				vector<float> values = (*m_generator)(former_state, ix, iy, iz);
				for (int k = 0; k < (int)values.size(); ++k)
				{
					*element_at(new_state[k], ix, iy, iz, size) = *element_at(former_state[k], ix, iy, iz, size) + m_delta_t*values[k];
				}
			}
		}
	}
}


void GrayScottSG::crankNicolsonScheme(const vector<float*>& former_state, vector<float*>& new_state)const
{
	int size = m_generator->size();
	VectorXf former_a(size*size*size);
	VectorXf offset_a(size*size*size);
	VectorXf former_b(size*size*size);
	for (int ix = 0; ix < size; ++ix)
	{
		for (int iy = 0; iy < size; ++iy)
		{
			for (int iz = 0; iz < size; ++iz)
			{
				int index = ix*size*size + iy*size + iz;
				former_a[index] = *element_at(former_state[0], ix, iy, iz, size);
				former_b[index] = *element_at(former_state[1], ix, iy, iz, size);
				offset_a[index] = m_delta_t*m_feed.getValue(ix, iy, iz);
			}
		}
	}
	VectorXf interaction = m_delta_t*former_a.cwiseProduct(former_b).cwiseProduct(former_b);

	VectorXf right_member_a = m_right_matrix_a*former_a - interaction + offset_a;
	VectorXf right_member_b = m_right_matrix_b*former_b + interaction;


	VectorXf new_a = m_solver_a.solveWithGuess(right_member_a, former_a);//TO DO : Try "right_member_a" as a guess

	VectorXf new_b = m_solver_b.solveWithGuess(right_member_b, former_b);//TO DO : Try "right_member_b" as a guess

	for (int ix = 0; ix < size; ++ix)
	{
		for (int iy = 0; iy < size; ++iy)
		{
			for (int iz = 0; iz < size; ++iz)
			{
				int index = ix*size*size + iy*size + iz;
				*element_at(new_state[0], ix, iy, iz, size) = new_a[index];
				*element_at(new_state[1], ix, iy, iz, size) = new_b[index];
			}
		}
	}
}



void GrayScottSG::operate(const vector<float*>& former_state, vector<float*>& new_state)const
{
	if (m_is_explicit)
	{
		explicitScheme(former_state, new_state);
	}
	else
	{
		crankNicolsonScheme(former_state, new_state);
	}
}


SparseMatrix<float> laplacian3DMatrix(int size, const Kernel3D<float>& kernel)
{
	int kernel_size = kernel.size();
	VectorXi reserve_vector = (kernel_size*kernel_size*kernel_size)*VectorXi::Ones(size*size*size);
	SparseMatrix<float> laplacian(size*size*size, size*size*size);
	int start = kernel.startIndex();
	int stop = kernel.stopIndex();
	laplacian.reserve(reserve_vector);
	for (int k = 0; k < size; ++k)
	{
		int x_pos = (int)(k / size / size);
		int y_pos = (int)((k - x_pos*size*size) / size);
		int z_pos = (int)(k - x_pos*size*size - y_pos*size);
		for (int relx = start; relx < stop; ++relx)
		{
			for (int rely = start; rely < stop; ++rely)
			{
				for (int relz = start; relz < stop; ++relz)
				{
					int index = min(max(x_pos + relx, 0), size - 1)*size*size + min(max(y_pos + rely, 0), size - 1)*size + min(max(z_pos + relz, 0), size - 1);
					laplacian.coeffRef(k, index) += kernel.getRelValue(relx, rely, relz);
					laplacian.
				}
			}
		}
	}
	return laplacian;
}