import csv
import pandas as pd
from datetime import datetime

def levelized_cost_of_electricity(E,CAPEX,OPEX,r=0.08):
    N = len(E)
    assert len(OPEX)==N,  "OPEX, CAPEX, and E must be the same length"
    assert len(CAPEX)==N, "OPEX, CAPEX, and E must be the same length"


    NPV_costs, NPV_gen = 0,0
    for ii in range(N):
        d = 1/(1+r)**ii 
        NPV_costs += d*(CAPEX[ii]+OPEX[ii])
        NPV_gen += d*E[ii]

    return NPV_costs/NPV_gen

def read_TMY(file):
    metadata = {}
    with open(file) as csvfile:
        reader = csv.reader(csvfile)
        info_headers = reader.__next__()
        info = reader.__next__()
        get_metadata = lambda h: info[info_headers.index(h)]
        for h in info_headers:
            metadata[h] = get_metadata(h)
            try:
                metadata[h] = float(metadata[h])
            except:
                continue
            
    tmy_data = pd.read_csv(file,skiprows=[0,1],
                           usecols=[i for i in range(23)])

    times = []
    for ii in range(tmy_data.shape[0]):
        times += [datetime(tmy_data.Year[ii],tmy_data.Month[ii],tmy_data.Day[ii],
                    tmy_data.Hour[ii],tmy_data.Minute[ii])]
    tmy_data['Datetimes'] = times

    return tmy_data,metadata