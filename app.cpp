#include "evolution.h"
#include "file_manager.h"

using namespace std;
using namespace Eigen;






int main(int argc, char** argv)
{

/*----------------------------------------------------------------------------------
							PARAMETERS SETTING
-----------------------------------------------------------------------------------*/

	Parameters params;

	float delta_t;//time step
	int num_steps = 2; //number of steps between two frames (I think a bit more efficient if it is even)
	float diff; //reference diffusion coefficoent
	float feed_default_value = 0.05;
	

	if (argc >= 5)//quickly test different parameters in command line
	//i.e :  prog.exe 1.0 0.041 0.055 0.062
	{
		delta_t = atof(argv[1]);
		diff = atof(argv[2]);
		feed_default_value = atof(argv[3]);
		params.kill = atof(argv[4]);
	}
	else
	{
		delta_t = 0.1;
		float diff = 0.041;
		feed_default_value = 0.048;
		params.kill = 0.064;
	}

	params.diff_a = 2 * diff; 
	params.diff_b = 0.5*diff;

	


/*----------------------------------------------------------------------------------
						         INITIAL FIELDS
-----------------------------------------------------------------------------------*/

	//string file_name_density = "E:/users/lbt10/code/test_density.fld";
	//Field3D<float> density_field = getField<float>(file_name_density);
	
	//int box_size = density_field.sizeX();
	int box_size = 100;
	cout << endl << "BOX SIZE : " << box_size << endl<<endl;

	/*---------------------------
	A and B fields initialization
	-----------------------------*/

	int n_centers = 20;//number of balls
	float radius = ((float)box_size) / 8.0;//radius of the balls
	float init_level_a = 0.5;//initial level of A where it exists
	float init_level_b = 0.5; //initial level of B where it exists
	float epsilon = 0.1; //noise level

	
	Field3D<float> init_field_a(box_size);
	init_field_a.initElements(init_level_a);

	Field3D<float> init_field_b(box_size);
	init_field_b.initElements(init_level_b);

	Field3D<float> noise_a = getRandomField(box_size, -epsilon, epsilon);
	Field3D<float> noise_b = getRandomField(box_size, -epsilon, epsilon);

	init_field_a = init_field_a + noise_a;
	init_field_b = init_field_b + noise_b;

	MatrixXf centers = getRandomCenters(n_centers, box_size);

	excludeBallsDensity(&init_field_a, centers, radius, 0.0);
	excludeBallsDensity(&init_field_b, centers, radius, 0.0);


	/*-----------------------
	Feed field initialization
	------------------------*/

	
	params.feed = getZeroField(box_size);
	params.feed.initElements(feed_default_value);//uniform
	
	/*
	params.feed = feed_default_value*density_field;
	*/

/*----------------------------------------------------------------------------------
							VTK RENDRING PIPELINE
-----------------------------------------------------------------------------------*/

	int T = 7200; //total period of animation / evolution

	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> ren_win = vtkSmartPointer<vtkRenderWindow>::New();
	ren_win->AddRenderer(renderer);
	ren_win->SetMultiSamples(0);
	vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	iren->SetRenderWindow(ren_win);
	
	ren_win->SetSize(1000, 700);

	ren_win->Render();

	vtkSmartPointer<vtkAnimationScene> scene = vtkSmartPointer<vtkAnimationScene>::New();
	scene->SetModeToRealTime();
	scene->SetLoop(0);
	scene->SetFrameRate(10);
	scene->SetStartTime(1);
	scene->SetEndTime(T);

	vtkSmartPointer<vtkAnimationCue> cue = vtkSmartPointer<vtkAnimationCue>::New();
	cue->SetStartTime(1);
	cue->SetEndTime(T);
	scene->AddCue(cue);

	ChemicalSystem chem_sys(box_size);
	chem_sys.initialize(init_field_a, init_field_b);
	chem_sys.setKernel(laplacianKernel(0.6));

	vtkSmartPointer<vtkAnnimationCueObserver> observer = vtkSmartPointer<vtkAnnimationCueObserver>::New();
	observer->m_ren_win = ren_win;
	observer->m_renderer = renderer;
	observer->m_chem_sys = &chem_sys;
	observer->SetParams(delta_t, num_steps, params);

	cue->AddObserver(vtkCommand::StartAnimationCueEvent, observer);
	cue->AddObserver(vtkCommand::AnimationCueTickEvent, observer);
	cue->AddObserver(vtkCommand::EndAnimationCueEvent, observer);
	cue->AddObserver(vtkCommand::KeyPressEvent, observer);

	
	scene->Play();

	scene->Stop();

	iren->AddObserver(vtkCommand::KeyPressEvent, observer);
	iren->Start();
	


	return EXIT_SUCCESS;
}
