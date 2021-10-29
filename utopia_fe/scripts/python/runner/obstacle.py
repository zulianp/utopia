import os
import yaml
import subprocess
import rich

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


class Executable(yaml.YAMLObject):
    yaml_tag = u'!Executable'
    def __init__(self, name):
        self.name = name
        self.path = self.find_path()
        self.check_path()
        self.verbose = True

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


    def run(self):
        exec_str = f"{self.path} --verbose"
        console.print(locals())
        console.print(exec_str, style='bold blue')
        return subprocess.run(exec_str,  shell=True, check=True, capture_output=(not self.verbose))

    def run_with_args(self, args):
        exec_str = f"{self.path} --verbose {args}"
        console.print(locals())
        console.print(exec_str, style='bold blue')
        return subprocess.run(exec_str,  shell=True, check=True, capture_output=(not self.verbose))

    def describe(self):
        console.print("Using runs executable: %s=%s" % (self.name, self.path))


class Env:
    """Env"""
    def __init__(self):
        utopia_fe_exec = Executable('UTOPIA_FE_EXEC')
        trilinos_decompo_exec = Executable('TRILINOS_DECOMP')

        utopia_fe_exec.describe()
        trilinos_decompo_exec.describe()

        self.utopia_fe_exec = utopia_fe_exec
        self.trilinos_decompo_exec = trilinos_decompo_exec


class Mesh(yaml.YAMLObject):
    yaml_tag = u'!Mesh'

    def __init__(self, path):
        self.type = 'file'
        self.path = path
        self.auto_aura = True

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
        self.input_file = os.environ.get('PWD') + '/' + self.app + '.yaml'

    def run(self):
        out = open(self.input_file, "w")
        n = out.write(yaml.dump(self))
        out.close()
        args = f"@file {self.input_file}"
        self.executable.run_with_args(args)


class ObstacleSimulation(yaml.YAMLObject):
    yaml_tag = u'!ObstacleSimulation'

    def __init__(self, env, space,  material, forcing_functions, obstacle, solver):
        self.app = 'stk_obs'
        self.executable = env.utopia_fe_exec
        self.space = space
        self.obstacle = obstacle
        self.input_file = os.environ.get('PWD') + '/' + self.app + '.yaml'
        self.solver = solver

        # problem = Problem(assembly, output_path)
        self.assembly = assembly = Assembly(material, forcing_functions)

    def run(self):
        out = open(self.input_file, "w")
        n = out.write(yaml.dump(self))
        out.close()
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

class Domain(yaml.YAMLObject):
    yaml_tag = u'!Domain'

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


if __name__ == '__main__':
    env = Env()

    solid_mesh = "/Users/zulianp/Desktop/in_the_cloud/owncloud_HSLU/Patrick/discharge_2/struct_halfDomain_21K.exo"
    fluid_mesh = "/Users/zulianp/Desktop/in_the_cloud/owncloud_HSLU/Patrick/discharge_2/fluid_halfDomain_155K.exo"

    implitic_obstacle_db = 'distance.e'
    distance_var =  Var("distance")

    distance_fs = FunctionSpace(
        Mesh(fluid_mesh),
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


    solid_fs = FunctionSpace(
        Mesh(solid_mesh),
        [Var("disp", 3)],
        [
            DirichletBC('valvesym', 0., 0),
            DirichletBC('valvesym', 0., 2)
        ]
    )

    solid_fs.read_state = False

    #########################################################

    domain = Domain(Mesh(implitic_obstacle_db), [
        distance_var
        ])

    step2_sim = ObstacleSimulation(env, solid_fs, VectorLaplaceOperator(), [], ImplicitObstacle(domain), MPRGP())

    print(yaml.dump(step2_sim))

    #########################################################



    step1_sim.run()
    step2_sim.run()

    # obs=MeshObstacle("mesh.e")
    # solver = MPRGP()
    # sim=ObstacleSimulation("stk_create_mesh", env.utopia_fe_exec, obs, solver)
    # sim.run()



    # # print(yaml.dump({'name': 'Silenthand Olleander', 'race': 'Human', 'traits': ['ONE_HAND', 'ONE_EYE']}))
    # print(yaml.dump(mprgp))

