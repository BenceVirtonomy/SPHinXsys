/**
 * @file 	electrophysiology.cpp
 * @brief 	This is the first test to validate our PED-ODE solver for solving
 * 			monodomain model closed by a physiology reaction.
 * @author 	Chi Zhang and Xiangyu Hu
 * @version 0.2.1
 * 			From here, I will denote version a beta, e.g. 0.2.1, other than 0.1 as
 * 			we will introduce cardiac electrophysiology and cardaic mechanics herein.
 * 			Chi Zhang
 */

/** SPHinXsys Library. */
#include "sphinxsys.h"
/** Namespace cite here. */
using namespace SPH;
/** Geometry parameter. */
Real L = 1.0; 	
Real H = 1.0;
/** Particle spacing and boudary dummy particles. */
Real particle_spacing_ref = H / 50.0;
Real BW = 4.0 * particle_spacing_ref;
/** Material properties. */
Real rho_0 = 1.0; 	
Real a_0[4] = {1.0, 0.0, 0.0, 0.0};
Real b_0[4] = {1.0, 0.0, 0.0, 0.0};
Vec2d d_0(1.0, 0.0);
/** Electrophysiology prameters. */
Real c_m = 1.0;
Real k = 8.0;
Real a = 0.15;
Real mu_1 = 0.2;
Real mu_2 = 0.3;
Real epsilon = 0.04;

std::vector<Point> CreatShape()
{
	std::vector<Point> shape;
	shape.push_back(Point(0.0, 0.0));
	shape.push_back(Point(0.0,  H));
	shape.push_back(Point( L,  H));
	shape.push_back(Point( L, 0.0));
	shape.push_back(Point(0.0, 0.0));
	return shape;
}
/** 
 * Define geometry and initial conditions of SPH bodies. 
 * */
class MuscleBody : public SolidBody
{
public:
	MuscleBody(SPHSystem &system, string body_name, int refinement_level, ParticlesGeneratorOps op)
		: SolidBody(system, body_name, refinement_level, op)
	{
		std::vector<Point> body_shape = CreatShape();
		body_region_.add_geometry(new Geometry(body_shape), RegionBooleanOps::add);
		body_region_.done_modeling();
	}
};
/**
 * Voltage observer body definition.
 */
class VoltageObserver : public ObserverLagrangianBody
{
public:
	VoltageObserver(SPHSystem &system, string body_name, int refinement_level, ParticlesGeneratorOps op)
		: ObserverLagrangianBody(system, body_name, refinement_level, op)
	{
		/** postion and volume. */
		body_input_points_volumes_.push_back(make_pair(Point(0.3, 0.7), 0.0));
	}
};

/** 
 * The main program. 
 */
int main()
{
	/** 
	 * Build up context -- a SPHSystem. 
	 */
	SPHSystem system(Vec2d(- BW, - BW), Vec2d(L + BW, H + BW), particle_spacing_ref);
		GlobalStaticVariables::physical_time_ = 0.0;
	/** 
	 * Configuration of materials, crate particle container and muscle body. 
	 */
	MuscleBody *muscle_body  =  new MuscleBody(system, "MuscleBody", 0, ParticlesGeneratorOps::lattice);
	Muscle 						material("Muscle", muscle_body, a_0, b_0,d_0, rho_0, 1.0);
	MuscleParticles 			particles(muscle_body);
	ElectrophysiologyReaction 	reaction_model("TwoVariableModel",muscle_body, c_m, k, a, mu_1, mu_2, epsilon);
	/**
	 * Particle and body creation of fluid observer.
	 */
	VoltageObserver *voltage_observer = new VoltageObserver(system, "VoltageObserver", 0, ParticlesGeneratorOps::direct);
	ObserverParticles 					observer_particles(voltage_observer);
	/** 
	 * Set body contact map. 
	 */
	SPHBodyTopology body_topology = { { muscle_body, {  } }, {voltage_observer, {muscle_body}} };
	system.SetBodyTopology(&body_topology);
	/**
	 * Simulation data structure set up.
	 */
	system.SetupSPHSimulation();
	/**
	 * The main dynamics algorithm is defined start here.
	 */
	/** 
	 * Initialization for electrophysiology computation. 
	 */	
	electro_physiology::ElectroPhysiologyInitialCondition 		initialization(muscle_body);
	/** 
	 * Iinitialization of specific test case.
	 */
	electro_physiology::TransmembranePotentialInitialization 	potential_initialization(muscle_body);
	/** 
	 * Corrected strong configuration. 
	 */	
	electro_physiology::CorrectConfiguration 					correct_configuration(muscle_body);
	/** 
	 * Time step size caclutation. 
	 */
	electro_physiology::getDiffusionTimeStepSize 				get_time_step_size(muscle_body);
	/** 
	 * Diffusion process for diffusion body. 
	 */
	electro_physiology::DiffusionRelaxation 					diffusion_relaxation(muscle_body);
	/** 
	 * Sovlers for ODE system 
	 */
	electro_physiology::TransmembranePotentialReaction 			voltage_reaction(muscle_body);
	electro_physiology::GateVariableReaction  					gate_variable_reaction(muscle_body);
	/** Outputs. */
	/**
	 * Simple input and outputs.
	 */
	In_Output 							in_output(system);
	WriteBodyStatesToPlt 				write_states(in_output, system.real_bodies_);
	WriteObservedVoltage	write_recorded_voltage(in_output, voltage_observer, {muscle_body});
	/** 
	 * Pre-simultion. 
	 */
	initialization.parallel_exec();
	potential_initialization.parallel_exec();
	correct_configuration.parallel_exec();
	/** 
	 * Output global basic parameters. 
	 */
	write_states.WriteToFile(GlobalStaticVariables::physical_time_);
	write_recorded_voltage.WriteToFile(GlobalStaticVariables::physical_time_);

	int ite 		= 0;
	Real T0 		= 16.0;
	Real End_Time 	= T0;
	Real D_Time 	= 0.5; 				/**< Time period for output */
	Real Dt 		= 0.001 * D_Time;	/**< Time period for data observing */
	Real dt		 	= 0.0;
	/** Statistics for computing time. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	/** Main loop starts here. */ 
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integeral_time = 0.0;
		while (integeral_time < D_Time) 
		{
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt) 
			{
				if (ite % 1000 == 0) 
				{
					cout << "N=" << ite << " Time: "
						<< GlobalStaticVariables::physical_time_ << "	dt: "
						<< dt << "\n";
				}
				/**Strang's splitting method. */
				voltage_reaction.parallel_exec(0.5 * dt);
				gate_variable_reaction.parallel_exec(0.5 * dt);

				diffusion_relaxation.parallel_exec(dt);

				gate_variable_reaction.parallel_exec(0.5 * dt);
				voltage_reaction.parallel_exec(0.5 * dt);

				ite++;
				dt = get_time_step_size.parallel_exec();
				relaxation_time += dt;
				integeral_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}
			write_recorded_voltage.WriteToFile(GlobalStaticVariables::physical_time_);
		}

		tick_count t2 = tick_count::now();
		write_states.WriteToFile(GlobalStaticVariables::physical_time_);
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	cout << "Total wall time for computation: " << tt.seconds() << " seconds." << endl;

	return 0;
}