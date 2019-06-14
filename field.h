#ifndef DEF_FIELD
#define DEF_FIELD
#pragma once
#include <iostream>
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <vector>


template<class T>
class Field3D
{
protected:
	int m_size_x;
	int m_size_y;
	int m_size_z;
	T* m_data;
	int* m_copies;
public:
	Field3D() { m_data = 0; }
	
	Field3D(int t_size)
		:m_size_x(t_size),
		m_size_y(t_size),
		m_size_z(t_size)
	{
		m_copies = new int(1);
		m_data = new T[t_size*t_size*t_size];
	}

	Field3D(int t_size_x, int t_size_y, int t_size_z)
		:m_size_x(t_size_x),
		m_size_y(t_size_y),
		m_size_z(t_size_z)
	{
		m_copies = new int(1);
		m_data = new T[t_size_x*t_size_y*t_size_z];
	}

	~Field3D()
	{
		if (*m_copies==1)
		{
			delete[] m_data;
			delete m_copies;
		}
		else
		{
			*m_copies -= 1;
		}
	}

	int sizeX()const 
	{
		return m_size_x;
	}

	int sizeY()const
	{
		return m_size_y;
	}

	int sizeZ()const
	{
		return m_size_z;
	}

	void initElements(T value)
	{
		for (int ix = 0; ix < m_size_x; ++ix)
		{
			for (int iy = 0; iy < m_size_y; ++iy)
			{
				for (int iz = 0; iz < m_size_z; ++iz)
				{
					this->setValue(ix, iy, iz, value);
				}
			}
		}
	}

	void shallowCopy(Field3D<T>& field)
	{
		m_size_x = field.sizeX();
		m_size_y = field.sizeY();
		m_size_z = field.sizeZ();
		m_data = field.getPointer();
		m_copies = field.getCopiesPointer();
		*m_copies += 1;
	}

	void deepCopy(const Field3D<T>& field)
	{
		m_size_x = field.sizeX();
		m_size_y = field.sizeY();
		m_size_z = field.sizeZ();
		m_data = new T[m_size_x*m_size_y*m_size_z];
		for (int ix = 0; ix < m_size_x; ++ix)
		{
			for (int iy = 0; iy < m_size_y; ++iy)
			{
				for (int iz = 0; iz < m_size_z; ++iz)
				{
					this->setValue(ix, iy, iz, field.getValue(ix, iy, iz));
				}
			}
		}
		m_copies = new int(1);
	}

	
	Field3D(const Field3D<T>& field)
	{
		this->deepCopy(field);
	}

	void setValue(int ix, int iy, int iz, T value)
	{
		int index = ix*m_size_y*m_size_z + iy*m_size_z + iz;
		assert(index >= 0 && index < m_size_x*m_size_y*m_size_z);
		m_data[index] = value;
	}

	T getValue(int ix, int iy, int iz)const
	{
		int index = ix*m_size_y*m_size_z + iy*m_size_z + iz;
		assert(index >= 0 && index < m_size_x*m_size_y*m_size_z);
		return m_data[index];
	}

	virtual T& operator()(int ix, int iy, int iz)
	{
		int index = ix*m_size_y*m_size_z + iy*m_size_z + iz;
		assert(index >= 0 && index < m_size_x*m_size_x*m_size_z);
		return m_data[index];
	}

	static Field3D<T> fromArray(std::vector<T> v_data, int t_size_x, int t_size_y, int t_size_z)
	{
		assert(t_size_x*t_size_y*t_size_z == v_data.size());
		Field3D<T> field(t_size_x, t_size_y, t_size_z);
		for (int ix = 0; ix < t_size_x; ++ix)
		{
			for (int iy = 0; iy < t_size_y; ++iy)
			{
				for (int iz = 0; iz < t_size_z; ++iz)
				{
					field.setValue(ix, iy, iz, v_data[ix*t_size_y*t_size_z + iy*t_size_z + iz]);
				}
			}
		}
		return field;
	}
	
	T* getPointer()
	{
		return m_data;
	}

	int* getCopiesPointer()
	{
		return m_copies;
	}

	void checkSizes(const Field3D<T>& field)const
	{
		assert(this->sizeX() == field.sizeX() && this->sizeY() == field.sizeY() && this->sizeZ() == field.sizeZ());
	}

	/*
	Classic element-wise addition when it is defined for type T.
	*/
	Field3D<T> operator+(const Field3D<T>& field)const
	{
		checkSizes(field);
		Field3D result(this->sizeX(), this->sizeY(), this->sizeZ());
		for (int ix = 0; ix < this->sizeX(); ++ix)
		{
			for (int iy = 0; iy < this->sizeY(); ++iy)
			{
				for (int iz = 0; iz < this->sizeZ(); ++iz)
				{
					result(ix, iy, iz) = this->getValue(ix, iy, iz) + field.getValue(ix, iy, iz);
				}
			}
		}
		return result;
	}

	Field3D<T> operator*(float scalar)const
	{
		Field3D result(this->sizeX(), this->sizeY(), this->sizeZ());
		for (int ix = 0; ix < this->sizeX(); ++ix)
		{
			for (int iy = 0; iy < this->sizeY(); ++iy)
			{
				for (int iz = 0; iz < this->sizeZ(); ++iz)
				{
					result(ix, iy, iz) = scalar*this->getValue(ix, iy, iz);
				}
			}
		}
		return result;
	}

	void operator=(const Field3D<T>& field)
	{
		this->deepCopy(field);
	}
};


//Does not work properly...
template<class T>
Field3D<T> operator*(float scalar, const Field3D<T>& field)
{
	return field*scalar;
}


template<class T>
class Kernel3D : public Field3D<T>
{
public:
	Kernel3D() {}

	Kernel3D(int t_size)
	{
		m_size_x = t_size;
		m_size_y = t_size;
		m_size_z = t_size;
		m_copies = new int(1);
		m_data = new T[t_size*t_size*t_size];
	}

	~Kernel3D()
	{
		if (*m_copies == 1)
		{
			delete[] m_data;
			delete m_copies;
		}
		else
		{
			m_copies -= 1;
		}
	}

	int size()const
	{
		return m_size_x;
	}


	int startIndex()const
	{
		return -(int)(m_size_x / 2);
	}

	int stopIndex()const
	{
		return m_size_x - (int)(m_size_x / 2);
	}
	
	Kernel3D(const Kernel3D<T>& kernel)
	{
		this->deepCopy(kernel);
	}

	/*
	define accessors with signed indices
	*/
	T getRelValue(int relx, int rely, int relz)const
	{
		int shift = (int)(m_size_x / 2);
		int index = (relx + shift)*m_size_y*m_size_z + (rely + shift)*m_size_z + (relz + shift);
		assert(index >= 0 && index < m_size_x*m_size_x*m_size_z);
		return m_data[index];
	}
	void setRelValue(int relx, int rely, int relz, float value)
	{
		int shift = (int)(m_size_x / 2);
		int index = (relx + shift)*m_size_y*m_size_z + (rely + shift)*m_size_z + (relz + shift);
		assert(index >= 0 && index < m_size_x*m_size_x*m_size_z);
		m_data[index] = value;
	}
	virtual T& operator()(int relx, int rely, int relz)
	{
		int shift = (int)(m_size_x / 2);
		int index = (relx+shift)*m_size_y*m_size_z + (rely+shift)*m_size_z + (relz+shift);
		assert(index >= 0 && index < m_size_x*m_size_x*m_size_z);
		return m_data[index];
	}
	
};


/*
float 3D field with all elements initialized to 0
*/
Field3D<float> getZeroField(int size_x, int size_y, int size_z);
Field3D<float> getZeroField(int size);

/*
float 3D field with random values between min_value and max_value
*/
Field3D<float> getRandomField(int size_x, int size_y, int size_z, float min_value, float max_value);
Field3D<float> getRandomField(int size, float min_value, float max_value);

/*
sets field values to "level" when corresponding points are inside a ball(center, radius) 
*/
void setBallDensity(Field3D<float>* field, Eigen::Vector3f center, float radius, float level);

/*
sets field values to "level" when corresponding points are inside one of the balls with random centers
*/
void setRandomBallsDensity(Field3D<float>* field, int number, float radius, float level);

/*
get random positions inside a box of specified size
*/
Eigen::MatrixXf getRandomCenters(int number, int size_x, int size_y, int size_z);
Eigen::MatrixXf getRandomCenters(int number, int size);

/*
everything outside of balls with specified centers and radius is set to zero
*/
void excludeBallsDensity(Field3D<float>* field, Eigen::MatrixXf centers, float radius, float exclude_value);

/*
3*3*3 discrete kernel for laplacian operator
*/
Kernel3D<float> laplacianKernel(float sigma);

/*
oriented laplacian kernel to increase the diffusion accross one axis in particular
*/
Kernel3D<float> laplacianLikeKernel(Eigen::Vector3f sigmas, Eigen::Vector3f orientation);



/*
very unfriendly accessors for vtkImageData elements in an image with size (size_x*size_y*size_z):
usage:
to access element at (ix, iy, iz),
replace "origin" by :  static_cast<float*>(vtk_pointer->GetScalarPointer())

where vtk_pointer is a vtkImageData*.
*/
float* element_at(float* origin, int ix, int iy, int iz, int size_x, int size_y, int size_z);
float* element_at(float* origin, int ix, int iy, int iz, int size);



#endif