#ifndef DEF_EVOLUTION
#define DEF_EVOLUTION
#pragma once

#include <iostream>
#include "Eigen/Dense"
#include <math.h>
#include <vector>

#include "field.h"
#include "generator.h"

#include <vtkActor.h>
#include <vtkContourFilter.h>
#include <vtkMath.h>
#include <vtkNamedColors.h>
#include <vtkPointData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkFloatArray.h>
#include <vtkStructuredPoints.h>
#include <vtkCommand.h>
#include <vtkAnimationScene.h>


struct Parameters
{
	float diff_a;
	float diff_b;
	Field3D<float> feed;
	float kill;
};




class ChemicalSystem
{
private:
	std::vector<float> m_masses;//total quantities of the two species accross the volume
	int m_size;//cubical box for the moment...
	bool m_has_kernel;//false if the usual trivial laplacian kernel is used
	bool m_has_semi_group;
	bool m_vector_update;
	Kernel3D<float> m_kernel;//laplacian kernel if defined
	SemiGroup* m_evolution;
	vtkImageData* m_image_to_show;
	float m_level;
public:
	std::vector<vtkImageData*> m_images;//3D scalar fields for A and B species  -  TO DO: put in private
	std::vector<vtkImageData*> m_buffers;// buffers for A and B species used during the updating steps  - TO DO: put in private
	
	

	vtkContourFilter* m_surface;//isosurface
	
	vtkPolyDataMapper* m_mapper;//part of vtk pipeline - not very interesting
	vtkActor* m_actor;//part of vtk pipeline - not very interesting
	
	ChemicalSystem(int t_size);
	~ChemicalSystem();
	void cleanUp();//
	
	float getLevel()const;
	void setLevel(float value);

	void initialize(float margin);
	void initialize(const Field3D<float>& field_A, const Field3D<float>& field_B);//copies field_A and field_B elements into m_images[0] and m_images[1]
	void setKernel(const Kernel3D<float>& kernel);//defines the laplacian kernel to use if different from the trivial one
	void setSemiGroup(SemiGroup* t_evolution);

	void setVectorUpdateOn();
	void setVectorUpdateOff();

	void computeMasses();

	void connect(vtkRenderer* p_renderer, float level);//connects evolution to the vtk rendering pipeline
	
	void update(float delta_t, int num_steps, Parameters params, vtkRenderer *p_renderer);
};





/*
class that manages vtk animation - not very interesting
*/
class vtkAnnimationCueObserver : public vtkCommand
{
private:
	
	float m_delta_t;
	int m_num_steps;
	Parameters m_params;
	int m_counter;
public:
	vtkRenderer* m_renderer;
	vtkRenderWindow* m_ren_win;
	ChemicalSystem* m_chem_sys;
	vtkAnnimationCueObserver();
	~vtkAnnimationCueObserver();
	static vtkAnnimationCueObserver *New();
	//void SetRenderWindow(vtkRenderWindow* ren_win);
	//void SetRenderer(vtkRenderer* t_renderer);
	void SetParams(float delta_t, int num_steps, Parameters params);
	virtual void Execute(vtkObject *vtkNotUsed(caller), unsigned long eventId, void * callData);
};





#endif