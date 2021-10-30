import os
import yaml
import subprocess
import rich
import sys, getopt
import glob

console = rich.get_console()

try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper


class LinearSolver(yaml.YAMLObject):
    yaml_tag = u'!LinearSolver'
    def __init__(self, backend, type, verbose: False):
        self.backend = backend
        self.type = type
        self.verbose = verbose
        # self.max_it =


class QPSolver(yaml.YAMLObject):
    yaml_tag = u'!QPSolver'
    def __init__(self, type):
        self.type = type
        self.backend = "any"
        self.verbose = True

    # def __repr__(self):
    #     return "%s(name=%r,backend=%r)" % (
    #         self.__class__.__name__, self.name, self.backend)

class MPRGP(QPSolver):
    yaml_tag = u'!MPRGP'
    def __init__(self):
        super().__init__("mprgp")
        self.backend = 'any'

    # def __repr__(self):
    #     return "%s(name=%r)" % (
    #         self.__class__.__name__, self.name)


    # verbose: true
    # max_it: 20
    # atol: 0.0001
    # max_constraints_iterations: 30
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
    def __init__(self, name, mpi_processes = 1):
        self.name = name
        self.path = self.find_path()
        self.check_path()
        self.verbose = True
        self.mpi_processes = mpi_processes

    # def __repr__(self):
    #     return "%s(name=%r,path=%r)" % (
    #         self.__class__.__name__, self.name, self.path)

    def find_path(self):
        path = os.environ.get(self.name)
        return path

    def check_path(self):
        if(not self.path):
            print("ERROR: Environment variable " + self.name + " must be defined! For instance with \'export " + self.name + "=<path_to_executable>")
            exit()

    def mpi_wrap(self, exec_str):
        if self.mpi_processes > 1:
            return f"mpiexec -np {self.mpi_processes} {exec_str}"
        else:
            return exec_str


    def run(self):
        exec_str = self.mpi_wrap(f"{self.path} ")

        console.print(locals())
        console.print(exec_str, style='bold blue')
        return subprocess.run(exec_str,  shell=True, check=True, capture_output=(not self.verbose))

    def run_with_args(self, args):
        exec_str = self.mpi_wrap(f"{self.path} {args}")
        console.print(locals())
        console.print(exec_str, style='bold blue')
        return subprocess.run(exec_str,  shell=True, check=True, capture_output=(not self.verbose))

    def describe(self):
        console.print("Using executable: %s=%s" % (self.name, self.path))


class Env:
    """Env"""
    step_id = 0

    def next_step_id(self):
        ret = self.step_id
        self.step_id += 1
        return ret

    def next_step_id_str(self):
        return str(self.next_step_id())

    def __init__(self, mpi_processes=1):
        utopia_fe_exec = Executable('UTOPIA_FE_EXEC')
        trilinos_decompo_exec = Executable('TRILINOS_DECOMP')

        utopia_fe_exec.describe()
        trilinos_decompo_exec.describe()

        self.utopia_fe_exec = utopia_fe_exec
        self.trilinos_decompo_exec = trilinos_decompo_exec
        self.mpi_processes = mpi_processes

        console.print("MPI processes: " + str(mpi_processes))


class Mesh(yaml.YAMLObject):
    yaml_tag = u'!Mesh'

    def __init__(self, env, path):
        self.type = 'file'
        self.path = path
        self.auto_aura = True
        self.env = env

        if env.mpi_processes > 1:
            pattern = path + "." + str(env.mpi_processes) + ".*"
            if not glob.glob(pattern):
                print("Decomposed file not found pattern=" + pattern)
                self.partition(env.mpi_processes)

    def partition(self, num_paritions):
        args = f"-p {num_paritions} {self.path}"
        self.env.trilinos_decompo_exec.run_with_args(args)

class DirichletBC(yaml.YAMLObject):
    yaml_tag = u'!DirichletBC'

    def __init__(self, name, value, var = 0):
        self.name = name
        self.value= value
        self.var = var


class Var(yaml.YAMLObject):
    yaml_tag = u'!Var'

    def __init__(self, name, n_components = 1, order = 'FIRST'):
        self.name = name
        self.n_components = n_components
        self.order = order


class FunctionSpace(yaml.YAMLObject):
    yaml_tag = u'!FunctionSpace'

    def __init__(self, mesh, variables, boundary_conditions):
        self.mesh = mesh
        self.variables = variables
        self.boundary_conditions = boundary_conditions



class FEProblem(yaml.YAMLObject):
    yaml_tag = u'!FEProblem'


class NLSolve(yaml.YAMLObject):
    yaml_tag = u'!ObstacleSimulation'

    def __init__(self, env, space, problem, solver):
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
        self.generate_YAML()
        args = f"@file {self.input_file}"
        self.executable.run_with_args(args)


class ObstacleSimulation(yaml.YAMLObject):
    yaml_tag = u'!ObstacleSimulation'

    def __init__(self, env, space,  material, forcing_functions, obstacle, solver, output_path):
        self.app = 'stk_obs'
        self.executable = env.utopia_fe_exec
        self.space = space
        self.obstacle = obstacle
        self.input_file = os.environ.get('PWD') + '/'  + env.next_step_id_str() + '_' + self.app + '.yaml'
        self.solver = solver
        self.output_path = output_path

        # problem = Problem(assembly, output_path)
        self.assembly = assembly = Assembly(material, forcing_functions)

    def generate_YAML(self):
        out = open(self.input_file, "w")
        n = out.write(yaml.dump(self))
        out.close()

    def run(self):
        self.generate_YAML()
        args = f"@file {self.input_file}"
        self.executable.run_with_args(args)

class Material(yaml.YAMLObject):
    yaml_tag = u'!Material'

    def __init__(self, type, quadrature_order = 0):
        self.type = type
        self.quadrature_order = quadrature_order

class LaplaceOperator(Material):
    yaml_tag = u'!LaplaceOperator'

    def __init__(self, quadrature_order = 0):
        super().__init__("LaplaceOperator", quadrature_order)

class VectorLaplaceOperator(Material):
    yaml_tag = u'!VectorLaplaceOperator'

    def __init__(self, quadrature_order = 0):
        super().__init__("VectorLaplaceOperator", quadrature_order)


class Assembly(yaml.YAMLObject):
    yaml_tag = u'!Assembly'

    def __init__(self, material, forcing_functions: []):
        self.material = material

        if(len(forcing_functions) != 0):
            self.forcing_functions = forcing_functions

class DomainForcingFunction(yaml.YAMLObject):
    yaml_tag = u'!DomainForcingFunction'

    def __init__(self, value, n_components, component: 0):
        self.value = value
        self.n_components = n_components
        self.component = component
        self.type = 'value'


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
        solver = LinearSolver('petsc', 'bdd', True)

        super().__init__(env, space, problem, solver)

class MeshObstacle(yaml.YAMLObject):
    yaml_tag = u'!MeshObstacle'

    def __init__(self, path):
        self.type = 'file'
        self.path = path
        self.gap_negative_bound = -5.0e-04
        self.gap_positive_bound = 5.e-4
        self.invert_face_orientation = True

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
        self.field_offset: -2e-5
        self.export_tensor = False



    # def __repr__(self):
    #     return "%s(name=%r,type=%r,path=%r)" % (
    #         self.__class__.__name__, self.name, self.type, self.name)

class FSI:
    def __init__(self, env, solid_mesh, fluid_mesh, workspace_dir):
        if not os.path.exists(workspace_dir):
            os.mkdir(workspace_dir)

        implitic_obstacle_db = workspace_dir + '/distance.e'
        implitic_obstacle_result_db = workspace_dir + '/contained_solid.e'

        #################################################################

        distance_var =  Var("distance")

        distance_fs = FunctionSpace(
            Mesh(env, fluid_mesh),
            [distance_var],
            [
                DirichletBC('inlet', 0.),
                DirichletBC('outlet', 0.),
                DirichletBC('inlettube', 0.),
                DirichletBC('walls', 0.),
                DirichletBC('symmetry', 0)
            ])


        step1_sim = GenerateDistanceFunction(env, distance_fs, 1000, implitic_obstacle_db)

        # console.print(yaml.dump(step1_sim))

        #########################################################

        solid_fs = FunctionSpace(
            Mesh(env, solid_mesh),
            [Var("disp", 3)],
            [
                DirichletBC('valvesym', 0., 0),
                DirichletBC('valvesym', 0., 2)
            ]
        )

        solid_fs.read_state = False

        #########################################################

        domain = DistanceField(
            Mesh(env, implitic_obstacle_db), [
                distance_var
            ])

        step2_sim = ObstacleSimulation(
            env, solid_fs, VectorLaplaceOperator(), [], ImplicitObstacle(domain),
            ObstacleSolver(MPRGP(), 30, 30), implitic_obstacle_result_db)

        # print(yaml.dump(step2_sim))

        #########################################################

        self.steps = [step1_sim, step2_sim]

    def run_step(self, step_num):
        self.steps[step_num].run()

    def run_all(self):
        for s in self.steps:
            s.run()

    def generate_YAML_files(self):
        for s in self.steps:
            s.generate_YAML()


#################################################################
def main(argv):
    # Inputs
    solid_mesh = "/Users/zulianp/Desktop/in_the_cloud/owncloud_HSLU/Patrick/discharge_2/struct_halfDomain_21K.exo"
    fluid_mesh = "/Users/zulianp/Desktop/in_the_cloud/owncloud_HSLU/Patrick/discharge_2/fluid_halfDomain_155K.exo"

    # Outputs
    workspace_dir = './workspace'
    mpi_processes = 1

    try:
        opts, args = getopt.getopt(argv,"hs:f:w:p:",["help", "solid=", "fluid=", "workspace=", "processes="])
    except getopt.GetoptError:
        print('obstacle.py -s <solid_mesh> -f <fluid_mesh>')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print('obstacle.py -i <inputfile> -o <outputfile>')
            sys.exit()
        elif opt in ("-s", "--solid"):
            solid_mesh = arg
        elif opt in ("-f", "--fluid"):
            fluid_mesh = arg
        elif opt in ("-p", "--processes"):
            mpi_processes = int(arg)
        elif opt in ("-w", "--workspace"):
           workspace_dir = arg

    env = Env(mpi_processes)
    fsi = FSI(env, solid_mesh, fluid_mesh, workspace_dir)
    # fsi.generate_YAML_files()
    # fsi.run_step(1)


if __name__ == '__main__':
    main(sys.argv[1:])
