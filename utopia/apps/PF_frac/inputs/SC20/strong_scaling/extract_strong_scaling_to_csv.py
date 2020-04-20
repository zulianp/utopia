import csv
import glob

# excecute this script inside the base directory that was used in the submit_jobs script
# i.e. the working directory should be in such a way that ./*nodes/summary.0.csv is found

allData = []
for filename in sorted(glob.glob('[0-9]*nodes/summary.0.csv')):
    numNodes = int(filename[0:3])
    print(f'numNodes={numNodes}')

    with open(filename) as f:
        allRows = [(r[0], r[2], r[5], r[3]) for r in csv.reader(f, delimiter=';')]
        assert(allRows[0][0]=='Class')
        assert(allRows[0][1]=='Total time spent (s)')
        assert(allRows[0][2]=='Mean time (s)')
        assert(allRows[0][3]=='Count')
        allRows = { r[0]: dict(total=float(r[1]), mean=float(r[2]), count=int(r[3])) for r in allRows[1:] }
        for k,r in allRows.items():
            assert(abs(r['mean']-r['total']/r['count']) < 1e-3)
        allData.append({'numNodes': numNodes, 'data': allRows})

csv_out = open('strong_scaling.csv', 'w')
csv_out.write('nodes,')
csv_out.write(','.join([f'{k}_speedup_mean,{k}_speedup_total,{k}_efficiency_mean,{k}_efficiency_total' for k in allData[0]['data'].keys()]))
for d in allData:
    csv_out.write(f'\n{d["numNodes"]}')
    for k,v in d['data'].items():
        v['speedupTotal'] = allData[0]['data'][k]['total'] / v['total']
        v['speedupMean'] = allData[0]['data'][k]['mean'] / v['mean']
        v['efficiencyTotal'] = v['speedupTotal'] * allData[0]['numNodes']/d['numNodes']
        v['efficiencyMean'] = v['speedupMean'] * allData[0]['numNodes']/d['numNodes']
        csv_out.write(f',{v["speedupMean"]},{v["speedupTotal"]},{v["efficiencyMean"]},{v["efficiencyTotal"]}')
