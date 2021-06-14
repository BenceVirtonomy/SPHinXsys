#include "solid_structural_simulation_class.h"

////////////////////////////////////////////////////
/* global functions in StructuralSimulation  */
////////////////////////////////////////////////////

BodyPartByParticleTriMesh::BodyPartByParticleTriMesh(SPHBody* body, string body_part_name, TriangleMeshShape* triangle_mesh_shape)
: BodyPartByParticle(body, body_part_name)
{	
	body_part_shape_ = new ComplexShape(body_part_name);
	body_part_shape_->addTriangleMeshShape(triangle_mesh_shape, ShapeBooleanOps::add);
	tagBodyPart();
}

ImportedModel::ImportedModel(SPHSystem &system, string body_name, TriangleMeshShape* triangle_mesh_shape, Real resolution)
	: SolidBody(system, body_name, resolution)
{
	ComplexShape original_body_shape;
	original_body_shape.addTriangleMeshShape(triangle_mesh_shape, ShapeBooleanOps::add);
	body_shape_ = new LevelSetComplexShape(this, original_body_shape, true);
}

SolidBodyForSimulation::SolidBodyForSimulation(SPHSystem &system, string body_name, TriangleMeshShape& triangle_mesh_shape, Real resolution, Real physical_viscosity, LinearElasticSolid& material_model):
	imported_model_(ImportedModel(system, body_name, &triangle_mesh_shape, resolution)),
	//material_model_(material_model),
	elastic_solid_particles_(ElasticSolidParticles(&imported_model_, &material_model)),
	inner_body_relation_(InnerBodyRelation(&imported_model_)),

	correct_configuration_(solid_dynamics::CorrectConfiguration(&inner_body_relation_)),
	stress_relaxation_first_half_(solid_dynamics::StressRelaxationFirstHalf(&inner_body_relation_)),
	stress_relaxation_second_half_(solid_dynamics::StressRelaxationSecondHalf(&inner_body_relation_)),
	damping_random_(DampingWithRandomChoice<DampingPairwiseInner<indexVector, Vec3d>>(&inner_body_relation_, 0.1, "Velocity", physical_viscosity))
{}

void ExpandBoundingBox(BoundingBox* original, BoundingBox* additional)
{
	for(int i = 0; i < original->first.size(); i++)
	{
		if ( additional->first[i] < original->first[i] )
		{
			original->first[i] = additional->first[i];
		}
		if ( additional->second[i] > original->second[i] )
		{
			original->second[i] = additional->second[i];
		}
	}
}

void RelaxParticlesSingleResolution(In_Output* in_output,
									bool write_particles_to_file,
									ImportedModel* imported_model,
									ElasticSolidParticles* imported_model_particles,
									InnerBodyRelation* imported_model_inner)
{	

	WriteBodyStatesToVtu write_imported_model_to_vtu(*in_output, { imported_model });
	WriteMeshToPlt write_mesh_cell_linked_list(*in_output, imported_model, imported_model->mesh_cell_linked_list_);

	//----------------------------------------------------------------------
	//	Methods used for particle relaxation.
	//----------------------------------------------------------------------
	RandomizePartilePosition  random_imported_model_particles(imported_model);
	/** A  Physics relaxation step. */
	relax_dynamics::RelaxationStepInner relaxation_step_inner(imported_model_inner, true);
	//----------------------------------------------------------------------
	//	Particle relaxation starts here.
	//----------------------------------------------------------------------
	random_imported_model_particles.parallel_exec(0.25);
	relaxation_step_inner.surface_bounding_.parallel_exec();
	if (write_particles_to_file)
	{
		write_imported_model_to_vtu.WriteToFile(0.0);
	} 
	imported_model->updateCellLinkedList();
	if (write_particles_to_file)
	{
		write_mesh_cell_linked_list.WriteToFile(0.0);
	}
	//----------------------------------------------------------------------
	//	Particle relaxation time stepping start here.
	//----------------------------------------------------------------------
	int ite_p = 0;
	while (ite_p < 500)
	{
		relaxation_step_inner.parallel_exec();
		ite_p += 1;
		if (ite_p % 100 == 0)
		{
			cout << fixed << setprecision(9) << "Relaxation steps for the imported model N = " << ite_p << "\n";
			if (write_particles_to_file)
			{
				write_imported_model_to_vtu.WriteToFile(Real(ite_p) * 1.0e-4);
			}
		}
	}
	cout << "The physics relaxation process of imported model finish !" << endl;
}

///////////////////////////////////////
/* StructuralSimulation members */
///////////////////////////////////////

StructuralSimulation::StructuralSimulation(StructuralSimulationInput* input) :
	relative_input_path_(input->relative_input_path),
	imported_stl_list_(input->imported_stl_list),
	scale_stl_(input->scale_stl),
	translation_list_(input->translation_list),
	default_resolution_(input->resolution_list[0]),
	resolution_list_(input->resolution_list),
	material_model_list_(input->material_model_list),
	physical_viscosity_(input->physical_viscosity),
	system_(SPHSystem(BoundingBox(Vec3d(0), Vec3d(0)), default_resolution_)),
	in_output_(In_Output(system_)),
	contacting_bodies_list_(input->contacting_bodies_list)
{
	// scaling of translation and resolution
	ScaleTranslationAndResolution();
	// creating body mesh list for triangular mesh shapes storage
	CreateBodyMeshList();
	// set up the system
	CalculateSystemBoundaries();
	system_.run_particle_relaxation_ = true;
	// initialize solid bodies with their properties
	InitializeElasticSolidBodies();
	// contacts
	InitializeAllContacts();
}

StructuralSimulation::~StructuralSimulation()
{
	// solid bodies
	for (auto sb: solid_body_list_)
	{
		delete sb;
	}
	// contact
	for (auto cl : contact_list_)
	{
		delete cl;
	}
	for (auto cd : contact_density_list_)
	{
		delete cd;
	}
	for (auto cf : contact_force_list_)
	{
		delete cf;
	}
}

//void StructuralSimulation::AddPrimitiveCuboid(Vec3d halfsize_cuboid, Vec3d translation, Real resolution, LinearElasticSolid& material)
//{
//    primitive_shape_list_.push_back( TriangleMeshShape(halfsize_cuboid, 20, translation) );
//	resolution_list_.push_back(resolution);
//	material_model_list_.push_back(material);
//}

void StructuralSimulation::ScaleTranslationAndResolution()
{
	for (int i = 0; i < translation_list_.size(); i++)
	{
		translation_list_[i] *= scale_stl_;
	}
	default_resolution_ *= scale_stl_;
	for (int i = 0; i < resolution_list_.size(); i++)
	{
		resolution_list_[i] *= scale_stl_;
	}
}

void StructuralSimulation::CalculateSystemBoundaries()
{	
	for (auto body_mesh: body_mesh_list_)
	{
		BoundingBox additional = body_mesh.findBounds();
		ExpandBoundingBox(&system_.system_domain_bounds_, &additional);
	}
}

void StructuralSimulation::CreateBodyMeshList()
{
	for (int i = 0; i < imported_stl_list_.size(); i++)
	{
		string relative_input_path_copy = relative_input_path_;
		body_mesh_list_.push_back(TriangleMeshShape(relative_input_path_copy.append(imported_stl_list_[i]), translation_list_[i], scale_stl_));
	}
}

void StructuralSimulation::InitializeElasticSolidBodies()
{
	for (int i = 0; i < body_mesh_list_.size(); i++)
	{
		SolidBodyForSimulation* sb = new SolidBodyForSimulation(system_, imported_stl_list_[i], body_mesh_list_[i], resolution_list_[i], physical_viscosity_, material_model_list_[i]);
		solid_body_list_.push_back(sb);
		RelaxParticlesSingleResolution(&in_output_, false, sb->GetImportedModel(), sb->GetElasticSolidParticles(), sb->GetInnerBodyRelation());
	}
}

void StructuralSimulation::InitializeContactBetweenTwoBodies(int first, int second)
{	
	ImportedModel* first_body = solid_body_list_[first]->GetImportedModel();
	ImportedModel* second_body = solid_body_list_[second]->GetImportedModel();
	SolidContactBodyRelation* first_contact = new SolidContactBodyRelation(first_body, {second_body});
	SolidContactBodyRelation* second_contact = new SolidContactBodyRelation(second_body, {first_body});
	contact_list_.push_back(first_contact);
	contact_list_.push_back(second_contact);

	contact_density_list_.push_back(new solid_dynamics::ContactDensitySummation (first_contact));
	contact_density_list_.push_back(new solid_dynamics::ContactDensitySummation (second_contact));

	contact_force_list_.push_back(new solid_dynamics::ContactForce (first_contact));
	contact_force_list_.push_back(new solid_dynamics::ContactForce (second_contact));
}

void StructuralSimulation::InitializeAllContacts()
{
	for (auto pair: contacting_bodies_list_)
	{
		InitializeContactBetweenTwoBodies(pair.first, pair.second);
	}
}

void StructuralSimulation::InitializeGravity()
{
	int i = 0;
	for (auto solid_body: solid_body_list_)
	{	
		if ( find(body_indeces_gravity_.begin(), body_indeces_gravity_.end(), i) != body_indeces_gravity_.end() )
		{	
			initialize_gravity_.push_back(new InitializeATimeStep(solid_body->GetImportedModel(), new Gravity(*gravity_[i])));
		}
		else
		{
			initialize_gravity_.push_back(new InitializeATimeStep(solid_body->GetImportedModel()));
		}
		i++;
	}
}

void StructuralSimulation::AddGravity(int body_index, Vec3d* gravity)
{
	body_indeces_gravity_.push_back(body_index);
	gravity_.push_back(gravity);
}

void StructuralSimulation::InitializeAccelerationForBodyPartInBoundingBox()
{	
	int i = 0;
	for (auto body_index: body_indeces_accelerations_)
	{
		acceleration_for_body_part_.push_back(new solid_dynamics::AccelerationForBodyPartInBoundingBox(solid_body_list_[body_index]->GetImportedModel(), bounding_boxes_[i], accelerations_[i]));
        i++;
    }
	
}

void StructuralSimulation::AddAccelerationForBodyPartInBoundingBox(int body_index, BoundingBox* bounding_box, Vec3d acceleration)
{
	body_indeces_accelerations_.push_back(body_index);
	bounding_boxes_.push_back(bounding_box);
	accelerations_.push_back(acceleration);
}

void StructuralSimulation::InitializeSpringDamperConstraintParticleWise()
{	
	int i = 0;
	for (auto body_index: body_indeces_spring_damper_)
	{
		spring_damper_contraint_.push_back(new solid_dynamics::SpringDamperConstraintParticleWise(solid_body_list_[body_index]->GetImportedModel(), stiffnesses_[i], damping_ratios_[i]));
        i++;
    }
}

void StructuralSimulation::AddSpringDamperConstraintParticleWise(int body_index, Vec3d stiffness, Real damping_ratio)
{
	body_indeces_spring_damper_.push_back(body_index);
	stiffnesses_.push_back(stiffness);
	damping_ratios_.push_back(damping_ratio);
}

void StructuralSimulation::InitializeConstrainSolidBodyRegion()
{	
	for (auto body_index: body_indeces_fixed_contraint_)
	{
		BodyPartByParticleTriMesh* bp = new BodyPartByParticleTriMesh(solid_body_list_[body_index]->GetImportedModel(), imported_stl_list_[body_index], &body_mesh_list_[body_index]);
		fixed_contraint_.push_back(new solid_dynamics::ConstrainSolidBodyRegion(solid_body_list_[body_index]->GetImportedModel(), bp));
	}
}

void StructuralSimulation::AddConstrainSolidBodyRegion(int body_index)
{
	body_indeces_fixed_contraint_.push_back(body_index);
}

void StructuralSimulation::ExecuteCorrectConfiguration()
{
	for (auto solid_body: solid_body_list_)
	{
		solid_body->GetCorrectConfiguration()->parallel_exec();
	}
}

void StructuralSimulation::ExecuteInitializeATimeStep()
{
	for (auto ig: initialize_gravity_)
	{
		ig->parallel_exec();
	}
}

void StructuralSimulation::ExecuteAccelerationForBodyPartInBoundingBox()
{
	for (auto acc: acceleration_for_body_part_)
	{
		acc->parallel_exec();
	}
}

void StructuralSimulation::ExecuteSpringDamperConstraintParticleWise()
{
	for (auto sd: spring_damper_contraint_)
	{
		sd->parallel_exec();
	}
}

void StructuralSimulation::ExecuteContactDensitySummation()
{
	for (auto cd: contact_density_list_)
	{
		cd->parallel_exec();
	}
}

void StructuralSimulation::ExecuteContactForce()
{
	for (auto cf: contact_force_list_)
	{
		cf->parallel_exec();
	}
}

void StructuralSimulation::ExecuteStressRelaxationFirstHalf(Real dt)
{
	for (auto solid_body : solid_body_list_)
	{
		solid_body->GetStressRelaxationFirstHalf()->parallel_exec(dt);
	}
}

void StructuralSimulation::ExecuteConstrainSolidBodyRegion()
{
	for (auto fc: fixed_contraint_)
	{
		fc->parallel_exec();
	}
}

void StructuralSimulation::ExecuteDamping(Real dt)
{
	for (auto solid_body : solid_body_list_)
	{
		solid_body->GetDampingWithRandomChoice()->parallel_exec(dt);
	}
}

void StructuralSimulation::ExecuteStressRelaxationSecondHalf(Real dt)
{
	for (auto solid_body : solid_body_list_)
	{
		solid_body->GetStressRelaxationSecondHalf()->parallel_exec(dt);
	}
}

void StructuralSimulation::ExecuteUpdateCellLinkedList()
{
	for (auto solid_body : solid_body_list_)
	{
		solid_body->GetImportedModel()->updateCellLinkedList();
	}
}

void StructuralSimulation::ExecuteContactUpdateConfiguration()
{
	for (auto cl: contact_list_)
	{
		cl->updateConfiguration();
	}
}

void StructuralSimulation::RunSimulationStep(int &ite, Real &dt, Real &integration_time)
{
	if (ite % 100 == 0) cout << "N=" << ite << " Time: " << GlobalStaticVariables::physical_time_ << "	dt: " << dt << "\n";

	/** ACTIVE BOUNDARY CONDITIONS */
	ExecuteInitializeATimeStep();
	ExecuteAccelerationForBodyPartInBoundingBox();
	ExecuteSpringDamperConstraintParticleWise();

	/** CONTACT */
	ExecuteContactDensitySummation();
	ExecuteContactForce();

	/** STRESS RELAXATOIN, DAMPING, PASSIVE CONSTRAINTS */
	ExecuteStressRelaxationFirstHalf(dt);
	ExecuteConstrainSolidBodyRegion();
	ExecuteDamping(dt);
	ExecuteConstrainSolidBodyRegion();
	ExecuteStressRelaxationSecondHalf(dt);
	
	/** UPDATE TIME STEP SIZE, INCREMENT */
	ite++;
	dt = system_.getSmallestTimeStepAmongSolidBodies();
	integration_time += dt;
	GlobalStaticVariables::physical_time_ += dt;
	
	/** UPDATE BODIES CELL LINKED LISTS */
	ExecuteUpdateCellLinkedList();
	
	/** UPDATE CONTACT CONFIGURATION */
	ExecuteContactUpdateConfiguration();
}

void StructuralSimulation::RunSimulation(Real end_time)
{
	WriteBodyStatesToVtu write_states(in_output_, system_.real_bodies_);
	GlobalStaticVariables::physical_time_ = 0.0;
	
	/** INITIALALIZE SYSTEM */
	system_.initializeSystemCellLinkedLists();
	system_.initializeSystemConfigurations();

	/** INITIAL CONDITION */
	ExecuteCorrectConfiguration();
	
	/** Statistics for computing time. */
	write_states.WriteToFile(GlobalStaticVariables::physical_time_);
	int ite = 0;
	Real output_period = 0.1 / 100.0;
	Real dt = 0.0;
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	/** Main loop */
	while (GlobalStaticVariables::physical_time_ < end_time)
	{
		Real integration_time = 0.0;
		while (integration_time < output_period) 
		{
			RunSimulationStep(ite, dt, integration_time);
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
}