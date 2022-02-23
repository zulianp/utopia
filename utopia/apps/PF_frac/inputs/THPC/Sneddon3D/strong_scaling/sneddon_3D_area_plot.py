import csv
import glob

# excecute this script inside the base directory that was used in the submit_jobs script
# i.e. the working directory should be in such a way that ./*nodes/summary.*.0.csv is found

methods = [
"main",
# "VcIsotropicPhaseFieldForBrittleFractures::hessian",
# "VcIsotropicPhaseFieldForBrittleFractures::gradient",
# "VcIsotropicPhaseFieldForBrittleFractures::value",
# "ProjectedGaussSeidel::apply", "ProjectedGaussSeidel::update",
# "MPRGP::solve(...)", "MPRGP::update"
"vcIsotropicPhaseFieldForBrittleFractures::hessian_local_assembly",
"PetscMatrix::write_unlock(...)",
"ProjectedGaussSeidel::step(...)",
"VcProjectedBlockGaussSeidelSweep::apply(...)"
]

check_methods = set(methods)


allData = []
for filename in sorted(glob.glob('[0-9]*nodes/balancing.*.csv')):
    numNodes = int(filename[0:3])
    print(f'numNodes={numNodes}')
    # print(filename)

    with open(filename) as f:
        allRows = [(r[0], r[1], r[2], r[3], r[4]) for r in csv.reader(f, delimiter=';')]
        # Class;Mean;Variance;Min;Max
        assert(allRows[0][0]=='Class')
        assert(allRows[0][1]=='Mean')
        assert(allRows[0][2]=='Variance')
        assert(allRows[0][3]=='Min')
        assert(allRows[0][4]=='Max')

        selectedRows = {}
        for r in allRows[1:]:
            if r[0] in check_methods:
                selectedRows[r[0]] = dict(mean=float(r[1]), var=float(r[2]), min=float(r[3]), max=float(r[4]))

        allData.append({'numNodes': numNodes, 'data': selectedRows})

csv_out = open('area.csv', 'w')
csv_out.write('nodes,')
csv_out.write(','.join([f'{k}' for k in methods]))

for d in allData:
    csv_out.write(f'\n{d["numNodes"]}')
    main_time = d['data']["main"]["max"]

    for k in methods:
        v = d['data'][k]
        a =  100.0*(v['max'] - v['min'])/main_time
        v['area'] =a
        print(f'{k}:{a}')
        # v['area'] = (v['var'])
        csv_out.write(f',{v["area"]}')


