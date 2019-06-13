#include "evolution.h"

using namespace std;
using namespace Eigen;



ChemicalSystem::ChemicalSystem(int t_size)
	:m_size(t_size),
	m_level(0.5)
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
	m_image_to_show = vtkImageData::New();
	m_image_to_show->SetDimensions(t_size, t_size, t_size);
	m_image_to_show->AllocateScalars(VTK_FLOAT, 1);


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
	if (m_image_to_show)
	{
		m_image_to_show->Delete();
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


float ChemicalSystem::getLevel()const
{
	return m_level;
}

void ChemicalSystem::setLevel(float value)
{
	m_level = value;
	cout << m_level << endl;
	m_surface->SetValue(0, value);
	m_surface->Update();
	m_surface->Modified();
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
				float value_a = 0.0;
				float value_b = 0.0;
				if (ix >= start && ix <= end && iy >= start && iy <= end && iz >= start && iz <= end)
				{
					value_a = 0.6;
					value_b = 0.4;
				}
				*element_at(origin_a, ix, iy, iz, m_size) = value_a;
				*element_at(origin_b, ix, iy, iz, m_size) = value_b;
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
	float* origin_show = static_cast<float*>(m_image_to_show->GetScalarPointer());
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

				float value_show = 0;
				if (value_a + value_b > 0)
				{
					value_show = value_b / (value_a + value_b);
				}
				*element_at(origin_show, ix, iy, iz, m_size) = value_show;
			}
		}
	}
}


void ChemicalSystem::connect(vtkRenderer* p_renderer, float level)
{
	vtkSmartPointer<vtkNamedColors> colors = vtkSmartPointer<vtkNamedColors>::New();

	m_surface = vtkContourFilter::New();
	//m_surface->SetInputData(m_images[0]);
	m_surface->SetInputData(m_image_to_show);
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


	float* origin_show = static_cast<float*>(m_image_to_show->GetScalarPointer());
	float* origin_a = static_cast<float*>(m_images[0]->GetScalarPointer());
	float* origin_b = static_cast<float*>(m_images[1]->GetScalarPointer());
	for (int ix = 0; ix < m_size; ++ix)
	{
		for (int iy = 0; iy < m_size; ++iy)
		{
			for (int iz = 0; iz < m_size; ++iz)
			{
				float value_a = *element_at(origin_a, ix, iy, iz, m_size);
				float value_b = *element_at(origin_b, ix, iy, iz, m_size);
				float value_show = 0;
				if (value_a + value_b > 0)
				{
					value_show = value_b / (value_a + value_b);
				}
				*element_at(origin_show, ix, iy, iz, m_size) = value_show;
			}
		}
	}

	float level_a = m_masses[0] / (m_size*m_size*m_size);
	float level_b = m_masses[1] / (m_size*m_size*m_size);

	cout << "Mass A = " << level_a << endl;
	cout << "Mass B = " << level_b << endl << endl;

	//m_surface->SetValue(0, min(0.95, 1.5*level_a));
	if (level_b > 0)
	{
		m_level = level_b / (level_a + level_b);
	}
	else
	{
		m_level = 0;
	}
	m_surface->SetValue(0, m_level);
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


void vtkAnnimationCueObserver::Execute(vtkObject *caller, unsigned long eventId, void * callData)
{
	float level = m_chem_sys->getLevel();
	vtkRenderWindowInteractor* iren = static_cast<vtkRenderWindowInteractor*>(caller);
	string symbol;
	switch (eventId)
	{
	case vtkCommand::StartAnimationCueEvent:
		m_chem_sys->connect(m_renderer, 0.05);
		break;
	case vtkCommand::EndAnimationCueEvent:
		//(void)m_renderer;
		//m_chem_sys->cleanUp();
		break;
	case vtkCommand::KeyPressEvent:
		symbol = iren->GetKeySym();
		if (symbol == "Down")
		{
			level = max(0.0, level - 0.01);
			m_chem_sys->setLevel(level);
			m_renderer->Render();
		}
		else if (symbol == "Up")
		{
			level = min(1.0, level + 0.01);
			m_chem_sys->setLevel(level);
			m_renderer->Render();
		}
		break;
	case vtkCommand::AnimationCueTickEvent:
		m_chem_sys->update(m_delta_t, m_num_steps, m_params, m_renderer);
		break;
	}
	m_ren_win->Render();
}

