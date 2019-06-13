#include "generator.h"

using namespace std;
using namespace Eigen;


GrayScottGenerator::GrayScottGenerator(float t_diff_a, float t_diff_b, float t_kill, const Field3D<float>& t_feed)
	:m_diff_a(t_diff_a),
	m_diff_b(t_diff_b),
	m_kill(t_kill),
	m_feed(t_feed),
	m_size(t_feed.sizeX())
{
	
}

vector<float> GrayScottGenerator::convolve(vector<float*> data, int ix, int iy, int iz)
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
					dd[i] += *element_at(data[i], ix + relx, iy + rely, iz + relz, m_size)*m_laplacian(relx, rely, relz);
				}
			}
		}
	}
	return dd;
}

vector<float> GrayScottGenerator::operator()(vector<float*> data, int ix, int iy, int iz)
{

	float value_a = *element_at(data[0], ix, iy, iz, m_size);
	float value_b = *element_at(data[1], ix, iy, iz, m_size);
	vector<float> dd = convolve(data, ix, iy, iz);
	vector<float> variation;
	float da = m_diff_a * dd[0] - value_a*value_b*value_b + m_feed(ix, iy, iz) * (1.0 - value_a);
	float db = m_diff_b * dd[1] + value_a*value_b*value_b - (m_feed(ix, iy, iz) + m_kill)*value_b;
	variation.push_back(da);
	variation.push_back(db);
	return variation;
}