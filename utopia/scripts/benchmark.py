#! /usr/bin/env python3

import csv
import subprocess, os
import argparse

utopia_root = "/Users/eric/Desktop/utopia/"
utopia_bin_dir = "utopia/bin/"

cmake_cmd = ["cmake", "..", "-DUTOPIA_LOG=ON"]
make_cmd = ["make", "-j8", "utopia_exec"]
utopia_cmd = [utopia_root + utopia_bin_dir + "utopia_exec"]

# Argument parsing
parser = argparse.ArgumentParser(description='Benchmarking utility for utopia')
parser.add_argument('--mpi', dest='mpi', nargs=1, type=int, help='Number of processes spawned by mpirun', default=[1])
parser.add_argument('-n', dest='execs', nargs=1, type=int, help='Number of executions of utopia', default=[100])
args = parser.parse_args()

iterations = args.execs[0]
if args.mpi[0] > 1:
    utopia_cmd = ["mpirun", "-n", str(args.mpi[0])] + utopia_cmd

# Main program
os.chdir(utopia_root + utopia_bin_dir)

status = subprocess.run(cmake_cmd)
if status.returncode != 0:
    print("Error in cmake! Exited with status code: " + str(status.returncode))
    exit(1)

status = subprocess.run(make_cmd)
if status.returncode != 0:
    print("Error in make! Exited with status code: " + str(status.returncode))
    exit(2)


# Execute benchmark
data = {}

print("Running utopia (1 out of " + str(iterations) + ")...")
status = subprocess.run(utopia_cmd)
if status.returncode != 0:
    print("Error in utopia#1! Exited with status code: " + str(status.returncode))
    exit(3)
with open('summary.0.csv', 'r') as csv_file:
    log_reader = csv.reader(csv_file, delimiter=';')
    next(log_reader, None)
    for row in log_reader:
        data[row[0]] = float(row[3])

for i in range(2, iterations + 1):
    print("Running utopia (" + str(i) + " out of " + str(iterations) + ")...")
    status = subprocess.run(utopia_cmd)
    if status.returncode != 0:
        print("Error in utopia#" + str(i) + "! Exited with status code: " + str(status.returncode))
        exit(3)
    with open('summary.0.csv', 'r') as csv_file:
        log_reader = csv.reader(csv_file, delimiter=';')
        next(log_reader, None)
        for row in log_reader:
            data[row[0]] = data[row[0]] + (float(row[3]) - data[row[0]]) / i

with open('benchmark.csv', 'w') as csv_file:
    output = csv.writer(csv_file, delimiter=';')
    output.writerow(["Class", "Mean time (s)"])
    for key in data.keys():
        output.writerow([key, "{:.9f}".format(data[key])])

print("Done")
