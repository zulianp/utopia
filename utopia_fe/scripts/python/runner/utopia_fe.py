import os
import yaml
import subprocess
import rich
import sys, getopt
import glob
import time
import copy
import platform

def check_python_version():
    v = platform.python_version_tuple()
    if(int(v[0]) < 3 or int(v[1]) < 7):
        sys.exit("Python version >= 3.7 is required!")

check_python_version()        

console = rich.get_console()

try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

class Newton(yaml.YAMLObject):
    yaml_tag = u'!Newton'
    def __init__(self, linear_solver, max_it = 40, atol = 1e-8, verbose = True):
        self.linear_solver = linear_solver
        self.verbose = verbose
        self.max_it = max_it
        self.atol = atol
        self.dumping = 1
        self.inverse_diagonal_scaling = False

class LinearSolver(yaml.YAMLObject):
    yaml_tag = u'!LinearSolver'
    def __init__(self, backend, type, verbose = False):
        self.backend = backend
        self.type = type
        self.verbose = verbose
        # self.max_it =
    def set_pc_type(self, pc_type):
        self.pc_type = pc_type

    def set_ksp_type(self, ksp_type):
        self.ksp_type = ksp_type

    def set_max_it(self, max_it):
        self.max_it = max_it

class AMG(LinearSolver):
    yaml_tag = u'!AMG'
    def __init__(self, verbose = False):
        super().__init__("petsc", "ksp", verbose)
        self.set_pc_type("hypre")
        self.set_ksp_type("gmres")

class LU(LinearSolver):
    yaml_tag = u'!LU'
    def __init__(self, verbose = False):
        super().__init__("petsc", "ksp", verbose)
        self.set_pc_type("lu")
        self.set_ksp_type("preonly")

class BDDLinearSolver(LinearSolver):
    yaml_tag = u'!LinearSolver'

    def __init__(self, max_it = 5000, verbose = False):
        super().__init__("petsc", "bdd", verbose)
        self.preconditioner_type = "amg"
        # self.inner_solver = LinearSolver("any", "cg", verbose)
        self.inner_solver = LinearSolver("petsc", "ksp", verbose)
        self.inner_solver.set_pc_type("hypre")
        self.inner_solver.set_ksp_type("gmres")
        self.inner_solver.max_it = max_it
        self.max_it = max_it

class QPSolver(yaml.YAMLObject):
    yaml_tag = u'!QPSolver'
    def __init__(self, type):
        self.type = type
        self.backend = "any"
        self.verbose = True

    def set_linear_solver(self, linear_solver):
        self.linear_solver = linear_solver

class MPRGP(QPSolver):
    yaml_tag = u'!MPRGP'
    def __init__(self,atol=1e-11):
        super().__init__("mprgp")
        self.backend = 'any'
        self.atol = atol

class InteriorPointSolver(QPSolver):
    yaml_tag = u'!InteriorPointSolver'
    def __init__(self,atol=1e-11):
        super().__init__("ipm")
        self.backend = 'any'
        self.atol = atol
        self.stol = 1e-10
        self.dumping_parameter = 0.95
        self.min_val = 1e-20
        # self.linear_solver = BDDLinearSolver(3000, False)
        self.linear_solver = AMG(False)
        self.max_it = 40
        # self.linear_solver = AMG(False)
        self.linear_solver.rtol = 1e-8

class BarrierQPSolver(QPSolver):
    yaml_tag = u'!BarrierQPSolver'
    def __init__(self, atol=1e-11):
        super().__init__("logbarrier")
        self.backend = 'any'
        self.function_type = 'BoundedLogBarrier'
        self.atol = atol
        self.infinity = 5.e-3
        self.trivial_obstacle = False
        self.barrier_parameter = 1e-9
        self.barrier_parameter_shrinking_factor = 0.5
        self.min_barrier_parameter = 1e-9
        self.barrier_parameter_thickness = 1e-4
        self.soft_boundary = 1e-16
        self.zero = 1e-20
        self.allow_projection = True
        self.enable_line_search = True
        self.verbose = True

class ObstacleSolver(yaml.YAMLObject):
    yaml_tag = u'!ObstacleSolver'

    def __init__(self, qp_solver, max_it = 20, max_constraints_iterations = 30, atol = 0.0001, verbose = True):
        self.qp_solver = qp_solver
        self.max_it = max_it
        self.max_constraints_iterations = max_constraints_iterations
        self.atol = atol
        self.verbose = verbose

class Executable(yaml.YAMLObject):
    yaml_tag = u'!Executable'
    def __init__(self, name, mpi_processes = 1, launch_cmd = "mpiexec"):
        self.name = name
        self.path = self.find_path()
        self.check_path()
        self.verbose = True
        self.mpi_processes = mpi_processes
        self.launch_cmd = launch_cmd

    def find_path(self):
        path = os.environ.get(self.name)
        return path

    def check_path(self):
        if(not self.path):
            print("ERROR: Environment variable " + self.name + " must be defined! For instance with \'export " + self.name + "=<path_to_executable>")
            exit()

    def mpi_wrap(self, exec_str):
        if self.mpi_processes > 1:
            return f"{self.launch_cmd} -np {self.mpi_processes} {exec_str}"
        else:
            return exec_str


    def run(self, current_work_dir = './'):
        exec_str = self.mpi_wrap(f"{self.path} ")

        console.print(locals())
        console.print(exec_str, style='bold blue')
        return subprocess.run(exec_str,  shell=True, check=True, capture_output=(not self.verbose))

    def run_with_args(self, args, current_work_dir = './'):
        exec_str = self.mpi_wrap(f"{self.path} {args}")
        console.print(locals())
        console.print(exec_str, style='bold blue')
        return subprocess.run(exec_str,  shell=True, check=True, capture_output=(not self.verbose), cwd = current_work_dir)

    def describe(self):
        console.print("Using executable: %s=%s" % (self.name, self.path))

class Env(yaml.YAMLObject):
    """Env"""
    yaml_tag = u'!Env'
    step_id = 0

    @classmethod
    def to_yaml(cls,dumper,data):
        new_data = copy.deepcopy(data)
        new_data.__dict__.clear()
        # for item in cls.hidden_fields:
        #     del new_data.__dict__[item]
        return dumper.represent_yaml_object(cls.yaml_tag, new_data, cls,
                                        flow_style=cls.yaml_flow_style)

    def next_step_id(self):
        ret = self.step_id
        self.step_id += 1
        return ret

    def next_step_id_str(self):
        return str(self.next_step_id())

    def __init__(self, mpi_processes=1):

        launch_cmd = "mpiexec"
        env_launch_cmd = os.environ.get('UTOPIA_LAUNCH_CMD')
        if env_launch_cmd:
            console.print(f"Overriding launch_cmd: {env_launch_cmd}")
            launch_cmd = env_launch_cmd

        utopia_fe_exec = Executable('UTOPIA_FE_EXEC', mpi_processes, launch_cmd)
        trilinos_decompo_exec = Executable('TRILINOS_DECOMP', 1, launch_cmd)
        trilinos_epu_exec = Executable('TRILINOS_EPU', 1, launch_cmd)


        utopia_fe_exec.describe()
        trilinos_decompo_exec.describe()
        trilinos_epu_exec.describe()

        self.utopia_fe_exec = utopia_fe_exec
        self.trilinos_decompo_exec = trilinos_decompo_exec
        self.trilinos_epu_exec = trilinos_epu_exec
        self.mpi_processes = mpi_processes


        console.print("MPI processes: " + str(mpi_processes))

class Mesh(yaml.YAMLObject):
    yaml_tag = u'!Mesh'

    def __init__(self, env, path, auto_decompose = True):
        self.type = 'file'
        self.path = path
        self.auto_aura = True
        self.env = env
        self.auto_decompose = auto_decompose

    def clone(self):
        c = Mesh(self.env, self.path, self.auto_decompose)
        c.type = self.type;
        c.auto_aura = self.auto_aura;
        c.env = self.env;
        return c

    def __repr__(self):
        return "%s(type=%r,path=%r,auto_aura=%r)" % (
            self.__class__.__name__, self.type, self.path, self.auto_aura)

    def find_decomposed_files(self):
        pattern = self.path + "." + str(self.env.mpi_processes) + ".*"
        decomposed_files = glob.glob(pattern)
        return sorted(decomposed_files)

    def ensure(self):
        if self.auto_decompose and self.env.mpi_processes > 1:
            decomposed_files = self.find_decomposed_files()
            if not decomposed_files:
                print("Decomposed file not found! Generating decomposition...")
                # print(os.path.dirname(path))
                self.partition(self.env.mpi_processes, os.path.dirname(self.path))
            else:
                time_decomp = os.path.getmtime(decomposed_files[0])

                if os.path.exists(self.path):
                    time_original = os.path.getmtime(self.path)

                    if time_decomp < time_original:
                        print("Decomposed file are older than original mesh file!\n"
                            "Original:\t" + time.ctime(time_original) + "\nDecomposition:\t" + time.ctime(time_decomp) +
                            ".\nGenerating decomposition...")
                        self.partition(self.env.mpi_processes, os.path.dirname(self.path))

                print("found files:\n")
                for f in decomposed_files:
                    print("\t" + f)

    def partition(self, num_paritions, output_dir):
        args = f"-p {num_paritions} {os.path.basename(self.path)}"
        self.env.trilinos_decompo_exec.run_with_args(args, output_dir)

    def unify(self):
        decomposed_files = self.find_decomposed_files()

        if decomposed_files:
            file0 = decomposed_files[0]
            output_dir = os.path.dirname(file0)
            base_name = os.path.basename(file0)
            args = f"-auto {base_name}"
            self.env.trilinos_epu_exec.run_with_args(args, output_dir)

class DirichletBC(yaml.YAMLObject):
    yaml_tag = u'!DirichletBC'

    def __init__(self, name, value, var = 0):
        self.name = name
        self.value= value
        self.var = var

    def clone(self):
        c = DirichletBC(self.name, self.value, self.var)
        return c

class Var(yaml.YAMLObject):
    yaml_tag = u'!Var'

    def __init__(self, name, n_components = 1, order = 'FIRST'):
        self.name = name
        self.n_components = n_components
        self.order = order

    def clone(self):
        c = Var(self.name, self.n_components, self.order)
        return c

class FunctionSpace(yaml.YAMLObject):
    yaml_tag = u'!FunctionSpace'

    def __init__(self, mesh, variables, boundary_conditions):
        self.mesh = mesh
        self.variables = variables
        self.boundary_conditions = boundary_conditions
        self.read_state = False

    def clone(self):
        c = FunctionSpace(self.mesh.clone(), copy.deepcopy(self.variables), copy.deepcopy(self.boundary_conditions))
        c.read_state = self.read_state
        return c

    def ensure(self):
        self.mesh.ensure()

class FEProblem(yaml.YAMLObject):
    yaml_tag = u'!FEProblem'


class NLSolve(yaml.YAMLObject):
    yaml_tag = u'!ObstacleSimulation'

    def __init__(self, env, space, problem, solver):
        self.env = env
        self.app = 'stk_nlsolve'
        self.space = space
        self.problem = problem
        self.solver = solver
        self.executable = env.utopia_fe_exec
        self.input_file = os.environ.get('PWD') + '/' + env.next_step_id_str() + '_' + self.app + '.yaml'

    def generate_YAML(self):
        out = open(self.input_file, "w")
        n = out.write(yaml.dump(self))
        out.close()

    def run(self):
        self.space.ensure()
        self.generate_YAML()
        args = f"@file {self.input_file}"
        self.executable.run_with_args(args)

        out_db = Mesh(self.space.mesh.env, self.problem.output_path)
        out_db.unify()


class ObstacleSimulation(yaml.YAMLObject):
    yaml_tag = u'!ObstacleSimulation'

    def __init__(self, env, space,  material, forcing_functions, obstacle, solver, output_path):
        self.env = env
        self.app = 'stk_obs'
        self.executable = env.utopia_fe_exec
        self.space = space
        self.obstacle = obstacle
        self.input_file = os.environ.get('PWD') + '/'  + env.next_step_id_str() + '_' + self.app + '.yaml'
        self.solver = solver
        self.output_path = output_path
        self.enable_line_search = True

        # problem = Problem(assembly, output_path)
        self.assembly = assembly = Assembly(material, forcing_functions)

    def generate_YAML(self):
        out = open(self.input_file, "w")
        n = out.write(yaml.dump(self))
        out.close()

    def run(self):
        self.space.ensure()
        self.generate_YAML()
        args = f"@file {self.input_file}"
        self.executable.run_with_args(args)

        out_db = Mesh(self.space.mesh.env, self.output_path)
        out_db.unify()

class Material(yaml.YAMLObject):
    yaml_tag = u'!Material'

    def __init__(self, type, quadrature_order = 0):
        self.type = type
        self.quadrature_order = quadrature_order

class Mass(Material):
    yaml_tag = u'!Mass'
    def __init__(self, density, n_components, quadrature_order = 2):
        super().__init__("Mass", quadrature_order)
        self.density = density
        self.n_components = n_components
        self.lumped = False
        self.verbose = True


class LaplaceOperator(Material):
    yaml_tag = u'!LaplaceOperator'

    def __init__(self, quadrature_order = 0):
        super().__init__("LaplaceOperator", quadrature_order)

class VectorLaplaceOperator(Material):
    yaml_tag = u'!VectorLaplaceOperator'

    def __init__(self, quadrature_order = 0):
        super().__init__("VectorLaplaceOperator", quadrature_order)

class ElasticMaterial(Material):
    yaml_tag = u'!ElasticMaterial'

    def __init__(self, type,  quadrature_order = 0):
        super().__init__(type, quadrature_order)

    def set_shear_modulus(self, shear_modulus):
        self.shear_modulus = shear_modulus

    def set_bulk_modulus(self, bulk_modulus):
        self.bulk_modulus = bulk_modulus

    def set_first_lame_parameter(self, first_lame_parameter):
        self.first_lame_parameter = first_lame_parameter

    def set_poisson_ratio(self, poisson_ratio):
        self.poisson_ratio = poisson_ratio

    def set_young_modulus(self, young_modulus):
        self.young_modulus = young_modulus

class LinearElasticity(ElasticMaterial):
    yaml_tag = u'!LinearElasticity'
    def __init__(self, quadrature_order = 0):
        super().__init__("LinearElasticity", quadrature_order)
        self.verbose = True

class CompressibleNeoHookean(ElasticMaterial):
    yaml_tag = u'!CompressibleNeoHookean'
    def __init__(self, quadrature_order = 0):
        super().__init__("NeoHookeanOgden", quadrature_order)
        self.verbose = True

class Assembly(yaml.YAMLObject):
    yaml_tag = u'!Assembly'

    def __init__(self, material, forcing_functions: [], debug = True):
        self.material = material
        self.debug = debug

        if(len(forcing_functions) != 0):
            self.forcing_functions = forcing_functions

class DomainForcingFunction(yaml.YAMLObject):
    yaml_tag = u'!DomainForcingFunction'

    def __init__(self, value, n_components, component: 0):
        self.value = value
        self.n_components = n_components
        self.component = component
        self.type = 'value'
        self.density = 1
        self.verbose = True


class Problem(yaml.YAMLObject):
    yaml_tag = u'!Problem'

    def __init__(self, assembly, output_path):
        self.assembly = assembly
        self.output_path = output_path


class GenerateDistanceFunction(NLSolve):
    yaml_tag = u'!DistanceProblem'

    def __init__(self, env, space, rhs_value, output_path):
        material = LaplaceOperator()
        forcing_functions = [DomainForcingFunction(rhs_value, 1, 0)]
        assembly = Assembly(material, forcing_functions)
        problem = Problem(assembly, output_path)
        solver = Newton(LinearSolver('petsc', 'bdd', True))

        super().__init__(env, space, problem, solver)

class Time(yaml.YAMLObject):
    yaml_tag = u'!Time'

    def __init__(self, delta, n_steps):
        self.delta = delta
        self.n_steps = n_steps

class BarrierProblem(yaml.YAMLObject):
    yaml_tag = u'!BarrierProblem'

    def __init__(self, space,  material, mass, forcing_functions, obstacle, time, output_path):
        self.obstacle = obstacle
        self.output_path = output_path
        self.time = time

        # self.function_type = 'LogBarrier'
        self.function_type = 'BoundedLogBarrier'
        self.assembly = Assembly(material, forcing_functions)
        self.mass = Assembly(mass, [])

        # Barrier numerics
        self.infinity = 5.e-3
        self.barrier_parameter = 1
        self.barrier_parameter_shrinking_factor = 0.5
        self.barrier_parameter_thickness = 1e-4
        self.min_barrier_parameter = 1e-6
        self.soft_boundary = 1e-16
        self.zero = 1e-20
        self.enable_line_search = True

        # Algorithmic branches
        self.trivial_obstacle = False

        self.allow_projection = True
        self.use_barrier_mass_scaling = True
        self.zero_initial_guess = False

        # self.allow_projection = False
        # self.use_barrier_mass_scaling = False
        # self.zero_initial_guess = True

class BarrierObstacleSimulation(NLSolve):
    yaml_tag = u'!BarrierObstacleSimulation'

    def __init__(self, env, space,  material, mass, forcing_functions, obstacle, time, solver, output_path):
        # self.integrator = 'ObstacleVelocityNewmark'
        self.integrator = 'ObstacleStabilizedVelocityNewmark'
        # self.integrator = 'Newmark'
        # self.integrator = 'VelocityNewmark'
        # self.integrator = 'ObstacleNewmark'
        # self.integrator = 'ObstacleStabilizedNewmark'
        problem =  BarrierProblem(space, material, mass, forcing_functions, obstacle, time, output_path)
        super().__init__(env, space, problem, solver)

class MeshObstacle(yaml.YAMLObject):
    yaml_tag = u'!MeshObstacle'

    def __init__(self, path):
        self.type = 'file'
        self.path = path
        self.gap_negative_bound = -1.0e-04
        self.gap_positive_bound = 1.e-4
        self.invert_face_orientation = True
        self.margin = 1e-6

class DistanceField(yaml.YAMLObject):
    yaml_tag = u'!DistanceField'

    def __init__(self, db, fields):
        self.mesh = db
        self.fields = fields

class ImplicitObstacle(yaml.YAMLObject):
    yaml_tag = u'!ImplicitObstacle'

    def __init__(self, domain):
        self.type = 'implicit'
        self.domain = domain
        self.shift_field = True
        self.field_rescale = 1
        self.field_offset = -2e-5
        self.export_tensor = False

class FSIInit:
    def __init__(self, env, solid_mesh, fluid_mesh, workspace_dir, containement_iterations = 30):
        if not os.path.exists(workspace_dir):
            os.mkdir(workspace_dir)

        implitic_obstacle_db = workspace_dir + '/distance.e'
        implitic_obstacle_result_db = workspace_dir + '/contained_solid.e'
        obstacle_at_rest_result_db = workspace_dir + '/obstacle_at_rest.e'
        nonlinear_obstacle_at_rest_result_db = workspace_dir + '/nonlinear_obstacle_at_rest.e'
        fsi_ready_db = workspace_dir + '/solid.e'
        solid_test_db = workspace_dir + '/solid_test.e'

        obstacle_mesh = Mesh(env, implitic_obstacle_db)

        #################################################################
        ### Step 0

        distance_var =  Var("distance")

        distance_fs = FunctionSpace(
            Mesh(env, fluid_mesh),
            [distance_var],
            [
                DirichletBC('INLET', 0.),
                DirichletBC('OUTLET', 0.),
                DirichletBC('WALLS', 0.),
                DirichletBC('SYMMETRY', 0)
            ])


        step0_sim = GenerateDistanceFunction(env, distance_fs, 1000, implitic_obstacle_db)

        #########################################################

        solid_fs = FunctionSpace(
            Mesh(env, solid_mesh),
            [Var("disp", 3)],
            [
                DirichletBC('SYMMETRY', 0., 0),
                DirichletBC('SYMMETRY', 0., 2),
            ]
        )

        solid_fs.read_state = False

        #########################################################
        ### Step 1

        distance_field = DistanceField(
            obstacle_mesh, [
                distance_var
            ])

        step1_sim = ObstacleSimulation(
            env, solid_fs, VectorLaplaceOperator(), [], ImplicitObstacle(distance_field),
            ObstacleSolver(MPRGP(), 30, containement_iterations), implitic_obstacle_result_db)

        #########################################################
        ### Step 2

        dummy_material = LinearElasticity(0)
        dummy_material.set_shear_modulus(0.2)
        dummy_material.set_first_lame_parameter(0.2)

        # Clone the space, since we are changing its content
        solid_fs_2 = solid_fs.clone()
        solid_fs_2.read_state = True
        solid_fs_2.mesh.path = implitic_obstacle_result_db

        step2_sim = ObstacleSimulation(
            env, solid_fs_2, dummy_material, [], MeshObstacle(fluid_mesh),
            ObstacleSolver(MPRGP(1e-14), 15, 15), obstacle_at_rest_result_db)

        #########################################################
        ### Step 3

        elastic_material = CompressibleNeoHookean(0)
        elastic_material.set_young_modulus(9.174e6)
        elastic_material.set_poisson_ratio(0.4)
        self.elastic_material = elastic_material

        # Clone the space, since we are changing its content
        solid_fs_3 = solid_fs.clone()
        solid_fs_3.read_state = True
        solid_fs_3.mesh.path = obstacle_at_rest_result_db

        step3_sim = ObstacleSimulation(
            env, solid_fs_3, elastic_material, [], MeshObstacle(fluid_mesh),
            ObstacleSolver(
                MPRGP(1e-14),
                # InteriorPointSolver(1e-14),
                15, 15), nonlinear_obstacle_at_rest_result_db)

        #########################################################
        ### Step 4


        solid_fs_4 = solid_fs.clone()
        solid_fs_4.read_state = True
        solid_fs_4.mesh.path = nonlinear_obstacle_at_rest_result_db

        barrier_solver = BarrierQPSolver(1e-14)

        step4_sim = ObstacleSimulation(
            env, solid_fs_4, elastic_material, [], MeshObstacle(fluid_mesh),
            ObstacleSolver(barrier_solver, 20, 20), fsi_ready_db)


        #########################################################
        ### Test (cheaper than FSI for testing barrier)

        mass = Mass(998, 3)

        solid_fs_test = solid_fs.clone()
        solid_fs_test.read_state = True
        solid_fs_test.mesh.path = fsi_ready_db

        force = DomainForcingFunction(-1e3, 3, 1)
        # force = DomainForcingFunction(0, 3, 1)
        force.density = mass.density
        forcing_functions = [force]

        # linear_solver = BDDLinearSolver(100, True)
        # linear_solver = AMG(False)
        linear_solver = LU(True)
        linear_solver.max_it = 400
        newton_solver = Newton(linear_solver, 80)
        # newton_solver.dumping = 1.
        newton_solver.inverse_diagonal_scaling = False
        newton_solver.dumping = 0.8
        # newton_solver.inverse_diagonal_scaling = True

        # newton_solver.max_it = 40
        newton_solver.stol = 1e-18
        newton_solver.atol = 1e-8
        newton_solver.rtol = 1e-14

        mesh_obstacle = MeshObstacle(fluid_mesh)
        mesh_obstacle.gap_negative_bound = -1e-4
        mesh_obstacle.gap_positive_bound = 1e-4

        test = BarrierObstacleSimulation(env, solid_fs_test,
         elastic_material, mass, forcing_functions,
         mesh_obstacle, Time(1e-4, 100), newton_solver, solid_test_db)

        #########################################################
        ## Store states for outside modification/manipulation

        self.test = test
        self.steps = [step0_sim, step1_sim, step2_sim, step3_sim, step4_sim]
        self.mass = mass
        self.force = force
        self.solid = solid_fs_test
        self.fluid_mesh_path = fluid_mesh
        self.barrier_solver = barrier_solver

    def run_step(self, step_num):
        self.steps[step_num].run()

    def run_all(self):
        for s in self.steps:
            s.run()

    def generate_YAML_files(self):
        for s in self.steps:
            s.generate_YAML()

        self.test.generate_YAML()


class Launcher:
    def generate_YAML(self):
        self.fsi.generate_YAML_files()

    def run_all(self):
        self.fsi.run_all()

    def run(self, run_step):
        self.fsi.run_step(run_step)

    def test(self):
        self.fsi.test.run()

    def __init__(self, argv):
        ###############################
        # INPUT
        ###############################

        # Solid mesh db
        solid_mesh = "solid.e"

        # Fluid mesh db
        fluid_mesh = "fluid.e"

        # Young's modulus for the solid material
        young_modulus = 9.174e6

        # Poisson's ratio
        poisson_ratio = 0.4

        # Initial barrier parameter
        barrier_parameter = 3e-9

        # Final barrier parameter, if 0 is set to barrier_parameter
        min_barrier_parameter = 0


        barrier_parameter_shrinking_factor = 0.5

        # Number of containement iterations. The higher the mesh resolution the more iterations are necessary
        containement_iterations = 30

        # Number of MPI processes for computation
        mpi_processes = 1

        ## TEST PROBLEM
        # Delta time for the test problem
        dt = 1e-5
        mass_density = 998


        ###############################
        # OUTPUT
        ###############################

        # Folder where all the outputs are saved
        workspace_dir = './workspace'

        # Generate YAML input-files
        generate_YAML = False

        ## Run programs

        # Run specific step [0, 4]
        run_step = -1

        # Run all steps
        run_all = False

        # Run test-case
        run_test = False

        try:
            opts, args = getopt.getopt(
                argv,"hs:f:w:p:r:gat",
                ["help", "solid=", "fluid=", "workspace=", "processes=", "run_step=", "generate", "all", "test",
                  "young_modulus=", "poisson_ratio=", "barrier_parameter=", "min_barrier_parameter=", "dt=", "mass=", "containement_iterations=", "barrier_parameter_shrinking_factor="
                ])

        except getopt.GetoptError as err:
            print(err)
            print('obstacle.py -s <solid_mesh> -f <fluid_mesh>')
            sys.exit(2)

        for opt, arg in opts:
            if opt in ('-h', '--help'):
                print('obstacle.py -s <solid_mesh> -f <fluid_mesh>')
                sys.exit()
            elif opt in ("-s", "--solid"):
                solid_mesh = arg
            elif opt in ("-f", "--fluid"):
                fluid_mesh = arg
            elif opt in ("-p", "--processes"):
                mpi_processes = int(arg)
            elif opt in ("-w", "--workspace"):
               workspace_dir = arg
            elif opt in ("-r", "--run_step"):
               run_step = int(arg)
            elif opt in ("-g", "--generate"):
               generate_YAML = True
            elif opt in ("-a", "--all"):
               run_all = True
            elif opt in ("-t", "--test"):
               run_test = True
            elif opt in ("--young_modulus"):
               young_modulus = float(arg)
            elif opt in ("--poisson_ratio"):
               poisson_ratio = float(arg)
            elif opt in ("--barrier_parameter"):
               barrier_parameter = float(arg)
            elif opt in ("--min_barrier_parameter"):
               min_barrier_parameter = float(arg)
            elif opt in ("--barrier_parameter_shrinking_factor"):
               barrier_parameter_shrinking_factor = float(arg)
            elif opt in ("--dt"):
               dt = float(arg)
            elif opt in ("--mass"):
               mass_density = float(arg)
            elif opt in ("--containement_iterations"):
               containement_iterations = int(arg)

        env = Env(mpi_processes)
        fsi = FSIInit(env, solid_mesh, fluid_mesh, workspace_dir, containement_iterations)

        #######################################################################
        # Material parmeters
        #######################################################################

        fsi.elastic_material.set_young_modulus(young_modulus)
        fsi.elastic_material.set_poisson_ratio(poisson_ratio)
        fsi.mass.density = mass_density
        fsi.force.density = mass_density

        #######################################################################
        # Barrier parameters
        #######################################################################

        if min_barrier_parameter == 0:
            min_barrier_parameter = barrier_parameter

        fsi.test.problem.barrier_parameter = barrier_parameter
        fsi.test.problem.min_barrier_parameter = min_barrier_parameter
        fsi.test.problem.barrier_parameter_shrinking_factor = barrier_parameter_shrinking_factor
        fsi.barrier_solver.barrier_parameter = barrier_parameter
        fsi.barrier_solver.min_barrier_parameter = min_barrier_parameter
        fsi.barrier_solver.barrier_parameter_shrinking_factor = barrier_parameter_shrinking_factor

        #######################################################################

        self.fsi = fsi

        if generate_YAML:
            self.generate_YAML()

        if run_all:
            self.run_all()
        elif run_step != -1:
            self.run(run_step)

        if run_test:
            self.test()


if __name__ == '__main__':
    # console.print(sys.argv[1:])
    Launcher(sys.argv[1:])
