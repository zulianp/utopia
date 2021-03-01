import csv
import glob

# excecute this script inside the base directory that was used in the submit_jobs script
# i.e. the working directory should be in such a way that ./*nodes/summary.*.0.csv is found

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

        allRows = { r[0]: dict(mean=float(r[1]), var=float(r[2]), min=float(r[3]), max=float(r[4])) for r in allRows[1:] }
        allData.append({'numNodes': numNodes, 'data': allRows})

csv_out = open('strong_scaling.csv', 'w')
csv_out.write('nodes,')
csv_out.write(','.join([f'{k}_mean,{k}_var,{k}_min,{k}_max,{k}_speed_up,{k}_efficiency' for k in allData[0]['data'].keys()]))

for d in allData:
    csv_out.write(f'\n{d["numNodes"]}')

    for k,v in d['data'].items():
        v['time'] = v['max']
        v['speed_up'] = allData[0]['data'][k]['max'] / v['max']
        v['efficiency'] = v['speed_up'] * allData[0]['numNodes']/d['numNodes']

        csv_out.write(f',{v["time"]},{v["var"]},{v["min"]},{v["max"]},{v["speed_up"]},{v["efficiency"]}')
