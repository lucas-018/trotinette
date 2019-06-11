#include "evolution.h"

using namespace std;
using namespace Eigen;

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

ChemicalSystem::ChemicalSystem(int t_size)
	:m_size(t_size)
{
	for (int i = 0; i < 2; ++i)
	{
		vtkImageData* image = vtkImageData::New();
		image->SetDimensions(t_size, t_size, t_size);
		image->AllocateScalars(VTK_FLOAT, 1);
		m_images.push_back(image);

		vtkImageData* buffer = vtkImageData::New();
		buffer->SetDimensions(t_size, t_size, t_size);
		buffer->AllocateScalars(VTK_FLOAT, 1);
		m_buffers.push_back(buffer);

		m_masses.push_back(0);
	}
	m_surface = 0;
	m_mapper = 0;
	m_actor = 0;
	m_has_kernel = false;
}

ChemicalSystem::~ChemicalSystem()
{
	for (int i = 0; i < 2; ++i)
	{
		if (m_images[i])
		{
			m_images[i]->Delete();
		}
		if (m_buffers[i])
		{
			m_buffers[i]->Delete();
		}
	}
	cleanUp();
}

void ChemicalSystem::cleanUp()
{
	if (m_surface)
	{
		m_surface->Delete();
		m_surface = 0;
	}
	if (m_mapper)
	{
		m_mapper->Delete();
		m_mapper = 0;
	}
	if (m_actor)
	{
		m_actor->Delete();
		m_actor = 0;
	}
}

void ChemicalSystem::initialize(float margin)
{
	int start = (int)(m_size / margin);
	int end = (int)(m_size - start);
	float* origin_a = static_cast<float*>(m_images[0]->GetScalarPointer());
	float* origin_b = static_cast<float*>(m_images[1]->GetScalarPointer());
	for (int ix = 0; ix < m_size; ++ix)
	{
		for (int iy = 0; iy < m_size; ++iy)
		{
			for (int iz = 0; iz < m_size; ++iz)
			{
				if (ix >= start && ix <= end && iy >= start && iy <= end && iz >= start && iz <= end)
				{
					*element_at(origin_a, ix, iy, iz, m_size) = 0.6;
					*element_at(origin_b, ix, iy, iz, m_size) = 0.4;
				}
				else
				{
					*element_at(origin_a, ix, iy, iz, m_size) = 0.0;
					*element_at(origin_b, ix, iy, iz, m_size) = 0.0;
				}
			}
		}
	}
}


void ChemicalSystem::setKernel(const Kernel3D<float>& kernel)
{
	m_kernel.deepCopy(kernel);
	m_has_kernel = true;
}


void ChemicalSystem::initialize(const Field3D<float>& field_A, const Field3D<float>& field_B)
{
	float* origin_a = static_cast<float*>(m_images[0]->GetScalarPointer());
	float* origin_b = static_cast<float*>(m_images[1]->GetScalarPointer());
	for (int ix = 0; ix < m_size; ++ix)
	{
		for (int iy = 0; iy < m_size; ++iy)
		{
			for (int iz = 0; iz < m_size; ++iz)
			{
				
				float value_a = field_A.getValue(ix, iy, iz);
				float value_b = field_B.getValue(ix, iy, iz);

				*element_at(origin_a, ix, iy, iz, m_size) = value_a;
				*element_at(origin_b, ix, iy, iz, m_size) = value_b;
				
				m_masses[0] += value_a;
				m_masses[1] += value_b;
			}
		}
	}
}


void ChemicalSystem::connect(vtkRenderer* p_renderer, float level)
{
	vtkSmartPointer<vtkNamedColors> colors = vtkSmartPointer<vtkNamedColors>::New();

	m_surface = vtkContourFilter::New();
	m_surface->SetInputData(m_images[0]);
	m_surface->SetValue(0, level);
	m_surface->Update();

	m_mapper = vtkPolyDataMapper::New();
	m_mapper->SetInputConnection(m_surface->GetOutputPort());
	m_mapper->ScalarVisibilityOff();

	m_actor = vtkActor::New();
	m_actor->SetMapper(m_mapper);
	m_actor->GetProperty()->SetColor(colors->GetColor3d("PaleTurquoise").GetData());

	p_renderer->AddActor(m_actor);
	p_renderer->SetBackground(colors->GetColor3d("PeachPuff").GetData());
	p_renderer->ResetCamera();
	p_renderer->Render();
}


void ChemicalSystem::update(float delta_t, int num_steps, Parameters params, vtkRenderer* p_renderer)
{
	for (int n = 0; n < num_steps; ++n)
	{
		float* old_a;
		float* old_b;
		float* new_a;
		float* new_b;

		switch (n % 2)
		{
		case 0:
			old_a = static_cast<float*>(m_images[0]->GetScalarPointer());
			old_b = static_cast<float*>(m_images[1]->GetScalarPointer());
			new_a = static_cast<float*>(m_buffers[0]->GetScalarPointer());
			new_b = static_cast<float*>(m_buffers[1]->GetScalarPointer());
			break;
		case 1:
			old_a = static_cast<float*>(m_buffers[0]->GetScalarPointer());
			old_b = static_cast<float*>(m_buffers[1]->GetScalarPointer());
			new_a = static_cast<float*>(m_images[0]->GetScalarPointer());
			new_b = static_cast<float*>(m_images[1]->GetScalarPointer());
			break;
		}

		for (int ix = 0; ix < m_size; ++ix)
		{
			for (int iy = 0; iy < m_size; ++iy)
			{
				for (int iz = 0; iz < m_size; ++iz)
				{
					float value_a = *element_at(old_a, ix, iy, iz, m_size);
					float value_b = *element_at(old_b, ix, iy, iz, m_size);

					float dda = 0;
					float ddb = 0;
					if (m_has_kernel)
					{
						for (int relx = -1; relx < 2; ++relx)
						{
							for (int rely = -1; rely < 2; ++rely)
							{
								for (int relz = -1; relz < 2; ++relz)
								{
									dda += *element_at(old_a, ix + relx, iy + rely, iz + relz, m_size)*m_kernel(relx, rely, relz);
								}
							}
						}
					}
					else
					{
						dda += *element_at(old_a, ix - 1, iy, iz, m_size);
						dda += *element_at(old_a, ix + 1, iy, iz, m_size);
						dda += *element_at(old_a, ix, iy - 1, iz, m_size);
						dda += *element_at(old_a, ix, iy + 1, iz, m_size);
						dda += *element_at(old_a, ix, iy, iz - 1, m_size);
						dda += *element_at(old_a, ix, iy, iz + 1, m_size);
						dda -= 6 * value_a;


						ddb += *element_at(old_b, ix - 1, iy, iz, m_size);
						ddb += *element_at(old_b, ix + 1, iy, iz, m_size);
						ddb += *element_at(old_b, ix, iy - 1, iz, m_size);
						ddb += *element_at(old_b, ix, iy + 1, iz, m_size);
						ddb += *element_at(old_b, ix, iy, iz - 1, m_size);
						ddb += *element_at(old_b, ix, iy, iz + 1, m_size);
						ddb -= 6 * value_b;
					}
					

					float da = params.diff_a*dda - value_a*value_b*value_b + params.feed(ix, iy, iz)*(1.0 - value_a);
					float db = params.diff_b*ddb + value_a*value_b*value_b - (params.feed(ix, iy, iz) + params.kill)*value_b;

					da += 1e-10f;
					db += 1e-10f;

					*element_at(new_a, ix, iy, iz, m_size) = value_a + delta_t*da;
					*element_at(new_b, ix, iy, iz, m_size) = value_b + delta_t*db;

					m_masses[0] += delta_t*da;
					m_masses[1] += delta_t*db;
				}
			}
		}
	}
	if (num_steps % 2 == 1)
	{
		m_images[0]->DeepCopy(m_buffers[0]);
		m_images[1]->DeepCopy(m_buffers[1]);
	}

	float level_a = m_masses[0] / (m_size*m_size*m_size);
	float level_b = m_masses[1] / (m_size*m_size*m_size);

	cout << "Mass A = " << level_a << endl;
	cout << "Mass B = " << level_b << endl << endl;

	m_surface->SetValue(0, level_a);
	m_surface->Update();
	m_surface->Modified();
	p_renderer->Render();
}


vtkAnnimationCueObserver::vtkAnnimationCueObserver()
{
	m_renderer = 0;
	m_ren_win = 0;
	m_chem_sys = 0;
	m_counter = 0;
}


vtkAnnimationCueObserver::~vtkAnnimationCueObserver()
{
	
}

vtkAnnimationCueObserver* vtkAnnimationCueObserver::New()
{
	return  new vtkAnnimationCueObserver;
}

void vtkAnnimationCueObserver::SetParams(float delta_t, int num_steps, Parameters params)
{
	m_delta_t = delta_t;
	m_num_steps = num_steps;
	m_params = params;
}


void vtkAnnimationCueObserver::Execute(vtkObject *vtkNotUsed(caller), unsigned long eventId, void * callData)
{
	switch (eventId)
	{
	case vtkCommand::StartAnimationCueEvent:
		m_chem_sys->connect(m_renderer, 0.05);
		break;
	case vtkCommand::EndAnimationCueEvent:
		(void)m_renderer;
		m_chem_sys->cleanUp();
		break;
	case vtkCommand::AnimationCueTickEvent:
		m_chem_sys->update(m_delta_t, m_num_steps, m_params, m_renderer);
		break;
	}
	m_ren_win->Render();
}