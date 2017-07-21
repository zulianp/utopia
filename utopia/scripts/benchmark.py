#! /usr/bin/env python3

import csv
import subprocess, os
import argparse

cmake_cmd = ["cmake", "..", "-DUTOPIA_LOG=ON"]
make_cmd = ["make", "-j8", "utopia_exec"]

# Columns in CSV files
CLASS = 0
TOTAL_TIME = 1
COUNT = 2
MEAN_TIME = 3
STD_DEV = 4
REL_STD_DEV = 5

# Argument parsing
parser = argparse.ArgumentParser(description='Benchmarking utility for utopia')
parser.add_argument('--mpi', dest='mpi', nargs=1, type=int, help='Number of processes spawned by mpirun', default=[1])
parser.add_argument('-n', dest='execs', nargs=1, type=int, help='Number of executions of utopia', default=[100])
parser.add_argument('--utopia-path', dest='path', type=str, help='Full path of utopia_exec', default="./utopia_exec")
parser.add_argument('args', nargs=argparse.REMAINDER, type=str, help="""Arguments to pass to utopia.
Put these arguments at the end of the command line, preceded by '--' (double dash)""")
args = parser.parse_args()

iterations = args.execs[0]
mpi_n = args.mpi[0]
utopia_cmd = [args.path] + args.args

# Main program
cwd = os.getcwd()
if utopia_cmd[0] == "./utopia_exec":
    if cwd.endswith('/scripts'):
        os.chdir('../bin')
    elif cwd.endswith('/utopia/utopia'):
        os.chdir('bin')
    elif cwd.endswith('/utopia'):
        os.chdir('utopia/bin')
    else:
        print('Error: unable to figure out where is the utopia executable. Exiting')
        exit(5)
else:
    u_dir, u_file = os.path.split(utopia_cmd[0])
    os.chdir(u_dir)
    print("Changed working directory to " + os.getcwd())
utopia_cmd[0] = "./utopia_exec"

if mpi_n > 1:
    utopia_cmd = ["mpirun", "-n", str(mpi_n)] + utopia_cmd

status = subprocess.run(cmake_cmd)
if status.returncode != 0:
    print("Error in cmake! Exited with status code: " + str(status.returncode))
    exit(1)

status = subprocess.run(make_cmd)
if status.returncode != 0:
    print("Error in make! Exited with status code: " + str(status.returncode))
    exit(1)

# Execute benchmark
data = [dict() for _ in range(mpi_n)]

print("Running utopia (1 out of " + str(iterations) + ")...")
status = subprocess.run(utopia_cmd)
if status.returncode != 0:
    print("Error in utopia#1! Exited with status code: " + str(status.returncode))
    exit(3)
for p in range(mpi_n):
    with open('summary.' + str(p) + '.csv', 'r') as csv_file:
        log_reader = csv.reader(csv_file, delimiter=';')
        next(log_reader, None)
        for row in log_reader:
            data[p][row[CLASS]] = (int(row[COUNT]), float(row[MEAN_TIME]))

for i in range(2, iterations + 1):
    print("Running utopia (" + str(i) + " out of " + str(iterations) + ")...")
    status = subprocess.run(utopia_cmd)
    if status.returncode != 0:
        print("Error in utopia#" + str(i) + "! Exited with status code: " + str(status.returncode))
        exit(3)
    for p in range(mpi_n):
        with open('summary.' + str(p) + '.csv', 'r') as csv_file:
            log_reader = csv.reader(csv_file, delimiter=';')
            next(log_reader, None)
            for row in log_reader:
                data[p][row[CLASS]] = (
                    data[p][row[CLASS]][0] + int(row[COUNT]),
                    data[p][row[CLASS]][1] + (float(row[MEAN_TIME]) - data[p][row[CLASS]][1]) / i
                )

os.chdir(cwd)
for p in range(mpi_n):
    with open('benchmark.' + str(p) + '.csv', 'w') as csv_file:
        output = csv.writer(csv_file, delimiter=';')
        output.writerow(["Class", "Count", "Mean time (s)"])
        for key in data[p].keys():
            output.writerow([key, data[p][key][0], "{:.9f}".format(data[p][key][1])])

print("Done")
