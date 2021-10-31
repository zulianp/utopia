from utopia_fe import *

#################################################################
def main(argv):
    # Inputs
    solid_mesh = "/Users/zulianp/Desktop/in_the_cloud/owncloud_HSLU/Patrick/discharge_2/struct_halfDomain_21K.exo"
    fluid_mesh = "/Users/zulianp/Desktop/in_the_cloud/owncloud_HSLU/Patrick/discharge_2/fluid_halfDomain_155K.exo"
    young_modulus = 9.174e6
    poisson_ratio = 0.4
    barrier_parameter = 1e-10
    dt = 1e-5

    # Outputs
    workspace_dir = './workspace'
    mpi_processes = 1
    run_step = -1
    generate_YAML = False
    run_all = False
    run_test = False

    try:
        opts, args = getopt.getopt(
            argv,"hs:f:w:p:r:gat",
            ["help", "solid=", "fluid=", "workspace=", "processes=", "run_step=", "generate", "all", "test",
              "young_modulus=", "poisson_ratio=", "barrier_parameter=", "dt="
            ])
    except getopt.GetoptError:
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
        elif opt in ("--dt"):
           dt = float(arg)


    env = Env(mpi_processes)
    fsi = FSIInit(env, solid_mesh, fluid_mesh, workspace_dir)

    #######################################################################
    # Change fundamental parmeters

    fsi.elastic_material.set_young_modulus(young_modulus)
    fsi.elastic_material.set_poisson_ratio(poisson_ratio)

    #######################################################################
    # Play with parameters in barrier for the test
    fsi.test.problem.infinity = 5.e-3
    fsi.test.problem.trivial_obstacle = False
    fsi.test.problem.barrier_parameter = barrier_parameter
    fsi.test.problem.barrier_parameter_shrinking_factor = 0.5
    fsi.test.problem.min_barrier_parameter = barrier_parameter
    fsi.test.problem.soft_boundary = 1e-12
    fsi.test.problem.zero = 1e-20
    fsi.test.problem.time.delta = dt
    fsi.test.problem.allow_projection = True
    #######################################################################

    if generate_YAML:
        fsi.generate_YAML_files()

    if run_all:
        fsi.run_all()
    elif run_step != -1:
        fsi.run_step(run_step)

    if run_test:
        fsi.test.run()

if __name__ == '__main__':
    # console.print(sys.argv[1:])
    main(sys.argv[1:])
