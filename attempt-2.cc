/*
First attempt at coding this project up
Date Started: 4-10-19
Date Milestone: 6-26-19 (the Neumann BCs)
Date Milestone: 7-10-19 (Gradient/Heat Flux of solution)
Author: Wolfgang Bangerth
Author: John Paul Morgan
Coding references: step-1, step-4, step-7, step-22, step-26, step-59
Website references: https://www.dealii.org/developer/doxygen/deal.II/step_26.html
*/

/*
From Lukas Bystrickey
https://people.sc.fsu.edu/~lb13f/bystricky_thesis.pdf
*/

/*
Here are the basic includes for all the things we will be doing
*/
#include <deal.II/base/utilities.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/smartpointer.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/sparse_direct.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_in.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/derivative_approximation.h>

#include <fstream>
#include <iostream>

#include <cmath>
#include <map>
#include <deal.II/grid/grid_tools.h>
#include <array>


// NEUMANN BOUNDARY CONDITIONS(step-7)
// Velocity Vectors(step-22)

// Note: Anywhere there is "Point<dim> &p", please remember that the points p will refer the the spatial axis. So p[0] = x, p[1] = y and p[2] = z

// Note: Within the code there were pieces that read "fe.degree+1" that were changed to 
// "fe->degree+1". This is to keep inline with step-7, which uses the package "smartpointer"

// The linear system we are solving is:
/*
(M + dt*theta*A)*u^{n} =
 M*u^{n-1} - dt*(1 - theta)*A*u^{n-1} + dt*[(1 - theta)*F^{n-1} + theta*F^{n}

With:
dt = time step size
M = mass matrix
A = stiffness matrix
theta = 0, 1/2, 1
*/

namespace Attempt2
{

	using namespace dealii;


	template <int dim>
	class HeatEquation // Note, this is a class, and NOT a void, so stuff stays here
{

public:// Here, we need to make the HeatEquation public so it can be accessed by every function
	//HeatEquation(); // Old, from step-26
	HeatEquation(const FiniteElement<dim> &fe); // New, from step-7

// Functions declared
void run();
void setup_system();
void solve_time_step();
void output_results() const;
void post_process() const;

Triangulation<dim> triangulation; // Defines object for triangulation
//FE_Q<dim> fe; // Defines the dimension for FEM. Old, from step-26
DoFHandler<dim> dof_handler; // Defines the DoF for which dimension we're in

SmartPointer<const FiniteElement<dim>> fe; // Needed for the new HeatEquation() with argument

ConstraintMatrix constraints;

SparsityPattern sparsity_pattern;
SparseMatrix<double> mass_matrix;
SparseMatrix<double> laplace_matrix;
SparseMatrix<double> system_matrix;

Vector<double> solution; // Solution for current timestep
Vector<double> gradient_vectorz; // For the gradiant of solution
Vector<double> old_solution; // Solution from previous timestep
Vector<double> system_rhs; // The RHS of A*x = b

double time;
double time_step;
unsigned int timestep_number;

const double theta; // Dictates which solver we will use(1/2 = Crank-Nicolson, 1 = Implicit Euler, 0 = Explicit Euler). 
//We LOVE Implicit Euler. This only uses the previous time step to solve for the current timestep
//We run away from Explicit Euler.

double k; // Our constant for the value of the Neumann BC on Gamma_3
double D;

};





	template <int dim>
	class InitialValues : public Function<dim> // IVs were copied from step-23
{

public:
InitialValues () : Function<dim>() {}

virtual double value (const Point<dim> &p, const unsigned int component = 0) const;

};





	template <int dim>
	double InitialValues<dim>::value (const Point<dim> &p, const unsigned int component) const
{

// We return our function at time t = 0
//return (std::sin(numbers::PI * p[0]) * std::sin(numbers::PI * p[1]) * std::cos(numbers::PI * p[2]));

return 0;

}





	template <int dim>
	class RightHandSide : public Function<dim> 
// This is where we will construct the RHS of the PDE. i.e.: u_t - (u_xx + u_yy + u_zz) = f
{

public:
RightHandSide ()
:
Function<dim>()
{
}

virtual double value (const Point<dim> &p, const unsigned int /*component = 0*/) const;

};





	template <int dim>
	double RightHandSide<dim>::value (const Point<dim> &/*p*/, const unsigned int /*component*/) const
{

return 0; // Here, f = 0 in the PDE u_t - (u_xx + u_yy + u_zz) = f

}





	template <int dim>
	class BoundaryValues : public Function<dim>// Grabs the BCs constructed in the next function
{

public:
virtual double value (const Point<dim> &p, const unsigned int /*component = 0*/) const;
// Here, override takes care of the "unused parameter 'component'". Use override for homo BC

};





	template <int dim>
	double BoundaryValues<dim>::value(const Point<dim> &/*p*/, const unsigned int /*component*/) const
{

// For nonhomogeneous Dirichlet BC, edit this part

// Also, for nonhomo BC, uncomment the /*p*/
/*
if ((p[0] <= 1) && (p[0] >= 0) && (p[1] <= 1) && (p[1] >= 0) && (p[2] == 0)) // Bottom of cube mesh
return (std::sin(p[0] * numbers::PI) * std::sin(p[1] * numbers::PI) * 1 * std::exp(-1 * 3 * numbers::PI * numbers::PI * this->get_time()));


else if ((p[0] <= 1) && (p[0] >= 0) && (p[1] <= 1) && (p[1] >= 0) && (p[2] == 1)) // Top of cube mesh
return (std::sin(p[0] * numbers::PI) * std::sin(p[1] * numbers::PI) * -1 * std::exp(-1 * 3 * numbers::PI * numbers::PI * this->get_time()));

else
return 0;
*/

// For homogeneous Dirichlet BC, use this part

// For homo BC, change &p to &/ * p * /
return 0;

}





	template <int dim> // Declare constants
	HeatEquation<dim>::HeatEquation (const FiniteElement<dim> &fe) // Note: This function first looks into HeatEquation
:
D(50), // Diffusion coeff
k(1), // Influx 
fe(&fe),
dof_handler(triangulation),
time (0.0),
time_step(1. / 10),
timestep_number (0),
theta(1) // Scheme (1/2 = Crank-Nicolson, 1 = Implicit Euler, 0 = Explicit Euler) 
{
}





	template <int dim>
	void HeatEquation<dim>::setup_system()
{

dof_handler.distribute_dofs(*fe);

std::cout << "The system is being set up." << std::endl;

constraints.clear ();

DoFTools::make_hanging_node_constraints (dof_handler,constraints);
// This computes constraints from hanging nodes. In our mesh, we will NOT have hanging nodes

constraints.close();

DynamicSparsityPattern dsp(dof_handler.n_dofs());

DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, /* keep_constrained_dofs = */ true);
sparsity_pattern.copy_from(dsp);

mass_matrix.reinit(sparsity_pattern);
laplace_matrix.reinit(sparsity_pattern);
system_matrix.reinit(sparsity_pattern);

MatrixCreator::create_mass_matrix(dof_handler, QGauss<dim>(fe->degree+1), mass_matrix);
MatrixCreator::create_laplace_matrix(dof_handler, QGauss<dim>(fe->degree+1), laplace_matrix);

solution.reinit(dof_handler.n_dofs());
gradient_vectorz.reinit(dof_handler.n_dofs());
old_solution.reinit(dof_handler.n_dofs());
system_rhs.reinit(dof_handler.n_dofs());

std::cout << "The system has been set up!" << std::endl;

}





	template <int dim>
	void HeatEquation<dim>::solve_time_step()
{

SolverControl solver_control(1000,1e-8 * system_rhs.l2_norm());
SolverCG<> cg(solver_control);

PreconditionSSOR<> preconditioner;
preconditioner.initialize(system_matrix, 1.0);

cg.solve(system_matrix, solution, system_rhs, preconditioner);
// cg.solve(A, x, b, precond)

constraints.distribute(solution);

//std::cout << "Solving for time step." << std::endl;

std::cout << "     " << solver_control.last_step()
          << " CG iterations." << std::endl; // This line will tell us how many CG iterations it took to solve the linear system

}






template <int dim>
class HeatFluxPostprocessor : public DataPostprocessorVector<dim>
{
public:
  HeatFluxPostprocessor ()
    :
    DataPostprocessorVector<dim> ("heat_flux",
                                  update_gradients | update_quadrature_points)
  {}
  virtual
  void
  evaluate_scalar_field (const DataPostprocessorInputs::Scalar<dim> &input_data,
                         std::vector<Vector<double> >               &computed_quantities) const
  {
    AssertDimension (input_data.solution_gradients.size(),
                     computed_quantities.size());
    for (unsigned int p=0; p<input_data.solution_gradients.size(); ++p)
      {
        AssertDimension (computed_quantities[p].size(), dim);
        for (unsigned int d=0; d<dim; ++d)
          // like above, but also multiply the gradients with
          // the coefficient evaluated at the current point:
        /*  
	computed_quantities[p][d]
            = coefficient (input_data.evaluation_points[p]) *
              input_data.solution_gradients[p][d];
	*/
		computed_quantities[p][d]
            = input_data.solution_gradients[p][d];      
	}
  }
};






template <int dim>
class GradientPostprocessor : public DataPostprocessorVector<dim>
{
public:
  GradientPostprocessor ()
    :
    DataPostprocessorVector<dim> ("grad_u",
                                  update_gradients)
  {}
  virtual
  void
  evaluate_scalar_field (const DataPostprocessorInputs::Scalar<dim> &input_data,
                         std::vector<Vector<double> >               &computed_quantities) const
  {
    AssertDimension (input_data.solution_gradients.size(),
                     computed_quantities.size());
    for (unsigned int p=0; p<input_data.solution_gradients.size(); ++p)
      {
        AssertDimension (computed_quantities[p].size(), dim);
        for (unsigned int d=0; d<dim; ++d)
          computed_quantities[p][d]
            = input_data.solution_gradients[p][d];
      }
  }
};









	template <int dim>
	void HeatEquation<dim>::output_results() const
{

GradientPostprocessor<dim> gradient_postprocessor;
HeatFluxPostprocessor<dim> heat_flux_postprocessor;

DataOut<dim> data_out;
data_out.attach_dof_handler(dof_handler);

data_out.add_data_vector(solution, "U");
data_out.add_data_vector(solution, gradient_postprocessor);
data_out.add_data_vector(solution, heat_flux_postprocessor);

data_out.build_patches();

const std::string filename = "/location/to/where/you/will/put/file/FILENAME-"
                                 + Utilities::int_to_string(timestep_number, 3) +
                                 ".vtk";
std::ofstream output(filename.c_str()); 

data_out.write_vtk(output);

}





	template <int dim>
	void HeatEquation<dim>::run()
{

GridIn<dim> gridin;
gridin.attach_triangulation(triangulation);
std::ifstream f("FILENAME.msh");
gridin.read_msh(f);

setup_system();

Vector<double> tmp;
Vector<double> forcing_terms;

tmp.reinit (solution.size());
forcing_terms.reinit (solution.size());

// For zeroes, use Functions::ZeroFunction<dim>()
// For custom IVs/ICs, use InitialValues<dim>()
VectorTools::interpolate(dof_handler, Functions::ZeroFunction<dim>(), old_solution);

solution = old_solution;

output_results();

//Begins boundary_vector section
const QGauss<dim-1> face_quadrature_formula(3); //Old: 3
FEFaceValues<dim> fe_face_values (*fe, face_quadrature_formula, 
					update_values | 
					update_quadrature_points | 
					update_JxW_values);

const unsigned int n_face_q_points = face_quadrature_formula.size();

const unsigned int dofs_per_cell = fe->dofs_per_cell;

Vector<double> boundary_vector(dofs_per_cell); // From step-7
std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);


while (time <= 10.00)
{

time += time_step;
++timestep_number;

std::cout << "Our time is " << time << std::endl;

mass_matrix.vmult(system_rhs, old_solution);
// system_rhs = mass_matrix * old_soln = M * u^{n-1}

laplace_matrix.vmult(tmp, old_solution);
// tmp = laplace_matrix * old_soln = A*u^{n-1}

system_rhs.add(-(1 - theta) * time_step * D, tmp);
// system_rhs = system_rhs + (-(1-theta)*time_step))*tmp = M*u^{n-1} - (1-theta)*dt*A*u^{n-1}


RightHandSide<dim> rhs_function;
rhs_function.set_time(time);
VectorTools::create_right_hand_side(dof_handler, QGauss<dim>(fe->degree+1), rhs_function, tmp);
forcing_terms = tmp;
forcing_terms *= time_step * theta; // This is the theta*F^{n}

rhs_function.set_time(time - time_step);
VectorTools::create_right_hand_side(dof_handler, QGauss<dim>(fe->degree+1), rhs_function, tmp);

forcing_terms.add(time_step * (1 - theta), tmp);
// The 8 previous lines contruct the term: dt*[(1 - theta)*F^{n-1} + theta*F^{n}]


system_rhs += forcing_terms;
// syrtem_rhs = system_rhs + forcing_terms
//This: = M*u^{n-1} - (1 - theta)*dt*A*u^{n-1} + dt*[(1 - theta)*F^{n-1} + theta*F^{n}]


//Here, we construct the Neumann BC vector
for (const auto &cell : dof_handler.active_cell_iterators())
 {
  boundary_vector = 0.;
  for (unsigned int face_num = 0; face_num < GeometryInfo<dim>::faces_per_cell; ++face_num)
   {
  
   if (cell->face(face_num)->at_boundary() && 
       (cell->face(face_num)->boundary_id() == 2)) // Represents Gamma_3
    {
     fe_face_values.reinit(cell, face_num);
      for (unsigned int q_point = 0; q_point < n_face_q_points; ++q_point)
       {
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
         boundary_vector(i) += (time_step *                             // delta_t
				k * 					// g(x_q)  
	 		      fe_face_values.shape_value(i, q_point) *  // phi_i(x_q)
			      fe_face_values.JxW(q_point));             // dx
       }
    }

/*
   if (cell->face(face_num)->at_boundary() && 
       (cell->face(face_num)->boundary_id() == 1)) // Represents Gamma_2, the sides
    {
     fe_face_values.reinit(cell, face_num);
      for (unsigned int q_point = 0; q_point < n_face_q_points; ++q_point)
       {
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
         boundary_vector(i) += (time_step *                             // delta_t
				10 * 					// g(x_q)  
	 		      fe_face_values.shape_value(i, q_point) *  // phi_i(x_q)
			      fe_face_values.JxW(q_point));             // dx
       }
    }
*/

  }


 cell->get_dof_indices (local_dof_indices);
  for (unsigned int i = 0; i < dofs_per_cell; ++i)
   {
    for (unsigned int j = 0; j < dofs_per_cell; ++j)
     system_rhs(local_dof_indices[i]) += boundary_vector(i);
   }

 }//Ends boundary_vector section


system_matrix.copy_from(mass_matrix);
system_matrix.add(D * theta * time_step, laplace_matrix);
// The 2 previous lines construct the term: M + dt*theta*A

constraints.condense (system_matrix, system_rhs);


// Here, we now apply our BCs to the remiaing boundary indicators. Any part of the system that had Dirichlet BCs are applied here.
{
BoundaryValues<dim> boundary_values_function;
boundary_values_function.set_time(time);

std::map<types::global_dof_index, double> boundary_values;
VectorTools::interpolate_boundary_values(dof_handler, 0, boundary_values_function, boundary_values); 
// Here, 0 -> Gamma_1

MatrixTools::apply_boundary_values(boundary_values, system_matrix, solution, system_rhs);
}


solve_time_step();

output_results();

old_solution = solution;


}

} // This ends the void HeatEquation<dim>::run()

}// This ends the namespace Attempt2





/*
Here, we will code everything inside the int main
*/
int main()
{
const unsigned int dim = 3;
try
{
using namespace dealii;
using namespace Attempt2;

std::cout << "The program has begun." << std::endl;

FE_Q<dim> fe(1); //1 = Linear, 2 = Quadratic, etc
HeatEquation<dim> heat_equation_solver(fe);
heat_equation_solver.run();
}

  catch (std::exception &exc) // Catches errors that have occured that deal.ii can understand
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl << exc.what()
                << std::endl << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...) // (...) must be the declared last in all catch ()'s. Catches all other error's deal.ii cannot understand
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl << "Aborting!"
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1; // return 1 means an error occured, and terminates the program prematurely
    }





std::cout << "The program has ended." << std::endl;

return 0;
}
