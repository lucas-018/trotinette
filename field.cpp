#include "field.h"

using namespace std;
using namespace Eigen;



Field3D<float> getZeroField(int size_x, int size_y, int size_z)
{
	Field3D<float> field(size_x, size_y, size_z);
	field.initElements(0);
	return field;
}


Field3D<float> getZeroField(int size)
{
	return getZeroField(size, size, size);
}



Field3D<float> getRandomField(int size_x, int size_y, int size_z, float min_value, float max_value)
{
	Field3D<float> field(size_x, size_y, size_z);
	VectorXf rand_numbers = VectorXf::Random(size_x*size_y*size_z);
	for (int ix = 0; ix < size_x; ++ix)
	{
		for (int iy = 0; iy < size_x; ++iy)
		{
			for (int iz = 0; iz < size_x; ++iz)
			{
				field(ix, iy, iz) = (float)(min_value + 0.5*(1.0 + rand_numbers(ix*size_y*size_z + iy*size_z + iz))*(max_value - min_value));
			}
		}
	}
	return field;
}

Field3D<float> getRandomField(int size, float min_value, float max_value)
{
	return getRandomField(size, size, size, min_value, max_value);
}


void setBallDensity(Field3D<float>* field, Vector3f center, float radius, float level)
{
	float dist;
	for (int ix = 0; ix < field->sizeX(); ++ix)
	{
		for (int iy = 0; iy < field->sizeY(); ++iy)
		{
			for (int iz = 0; iz < field->sizeZ(); ++iz)
			{
				Vector3f point((float)ix, (float)iy, (float)iz);
				dist = (point - center).norm();
				if (dist <= radius)
				{
					field->setValue(ix, iy, iz, level);
				}
			}
		}
	}
}


void setRandomBallsDensity(Field3D<float>* field, int number, float radius, float level)
{
	for (int i = 0; i < number; ++i)
	{
		Vector3f rand = Vector3f::Random();
		Vector3f center(0, 0, 0);
		center(0) = (float)(0.5*(1.0 + rand(0))*field->sizeX());
		center(1) = (float)(0.5*(1.0 + rand(1))*field->sizeY());
		center(2) = (float)(0.5*(1.0 + rand(2))*field->sizeZ());

		setBallDensity(field, center, radius, level);
	}
}

MatrixXf getRandomCenters(int number, int size_x, int size_y, int size_z)
{
	MatrixXf mat(number, 3);
	for (int i = 0; i < number; ++i)
	{
		Vector3f rand = Vector3f::Random();
		mat(i, 0) = (float)(0.5*(1.0 + rand(0))*size_x);
		mat(i, 1) = (float)(0.5*(1.0 + rand(1))*size_y);
		mat(i, 2) = (float)(0.5*(1.0 + rand(2))*size_z);
	}
	return mat;
}


MatrixXf getRandomCenters(int number, int size)
{
	return getRandomCenters(number, size, size, size);
}


void excludeBallsDensity(Field3D<float>* field, Eigen::MatrixXf centers, float radius, float exclude_value)
{
	for (int ix = 0; ix < field->sizeX(); ++ix)
	{
		for (int iy = 0; iy < field->sizeX(); ++iy)
		{
			for (int iz = 0; iz < field->sizeX(); ++iz)
			{
				Vector3f point((float)ix, (float)iy, (float)iz);
				bool is_outside = true;
				for (int k = 0; k < centers.rows(); ++k)
				{
					Vector3f center = centers.row(k);
					float dist = (point - center).norm();
					if (dist <= radius)
					{
						is_outside = false;
						break;
					}
				}
				if (is_outside)
				{
					field->setValue(ix, iy, iz, exclude_value);
				}
			}
		}
	}
}


Kernel3D<float> laplacianKernel(float sigma)
{
	Kernel3D<float> laplacian(3);
	float denom = 0;
	for (int relx = -1; relx < 2; ++relx)
	{
		for (int rely = -1; rely < 2; ++rely)
		{
			for (int relz = -1; relz < 2; ++relz)
			{
				if (relx == 0 && rely == 0 && relz == 0)
				{
					laplacian(relx, rely, relz) = -1.0;
				}
				else
				{
					Vector3f point((float)relx, (float)rely, (float)relz);
					float dist = point.norm();
					float value = exp(-dist*dist / (2 * sigma*sigma));
					laplacian(relx, rely, relz) = value;
					denom += value;
				}
			}
		}
	}

	for (int relx = -1; relx < 2; ++relx)
	{
		for (int rely = -1; rely < 2; ++rely)
		{
			for (int relz = -1; relz < 2; ++relz)
			{
				if (relx != 0 || rely != 0 || relz != 0)
				{
					laplacian(relx, rely, relz) /= denom;
				}
			}
		}
	}
	return laplacian;
}


Kernel3D<float> laplacianLikeKernel(Vector3f sigmas, Vector3f orientation)
{
	Vector3f ref_vector(1, 0, 0);
	Quaternion<float> q_rotation = Quaternion<float>::FromTwoVectors(ref_vector, orientation);
	Kernel3D<float> oriented_laplacian(3);
	float denom = 0;
	for (int relx = -1; relx < 2; ++relx)
	{
		for (int rely = -1; rely < 2; ++rely)
		{
			for (int relz = -1; relz < 2; ++relz)
			{
				if (relx != 0 || rely || 0 && relz || 0)
				{
					Vector3f eval_point((float)relx, (float)rely, (float)relz);
					Vector3f turned = q_rotation._transformVector(eval_point);
					float distance_2 = turned.dot(sigmas.cwiseInverse()) / 2.0f;
					float value = exp(-distance_2);
					denom += value;
					oriented_laplacian(relx, rely, relz) = value;
				}
				else
				{
					oriented_laplacian(relx, rely, relz) = -1.0;
				}
			}
		}
	}
	for (int relx = -1; relx < 2; ++relx)
	{
		for (int rely = -1; rely < 2; ++rely)
		{
			for (int relz = -1; relz < 2; ++relz)
			{
				if (relx != 0 || rely || 0 && relz || 0)
				{
					oriented_laplacian(relx, rely, relz) /= denom;
				}
			}
		}
	}
	return oriented_laplacian;
}



/*
float* element_at(float* origin, int ix, int iy, int iz, int size_x, int size_y, int size_z)
{
int ind_x = ix;
if (ind_x < 0)
{
ind_x = size_x + ind_x;
}
else if (ind_x >= size_x)
{
ind_x = ind_x - size_x;
}
int ind_y = iy;
if (ind_y < 0)
{
ind_y = size_y + ind_y;
}
else if (ind_y >= size_y)
{
ind_y = ind_y - size_y;
}
int ind_z = iz;
if (ind_z < 0)
{
ind_z = size_z + ind_z;
}
else if (ind_z >= size_z)
{
ind_z = ind_z - size_z;
}
return origin + ind_x*size_y*size_z + ind_y*size_z + ind_z;
}*/


float* element_at(float* origin, int ix, int iy, int iz, int size_x, int size_y, int size_z)
{
	return origin + min(max(0, ix), size_x - 1)*size_y*size_z + min(max(0, iy), size_y - 1)*size_z + min(max(0, iz), size_z - 1);
}

float* element_at(float* origin, int ix, int iy, int iz, int size)
{
	return element_at(origin, ix, iy, iz, size, size, size);
}