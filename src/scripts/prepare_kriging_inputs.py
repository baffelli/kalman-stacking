
import pyrat.ipt.core as ipt
import csv

def prepare_inputs(inputs, outputs, threads, config, params, wildcards):
    type_map = {'pint': 'FCOMPLEX', 'punw': 'FLOAT', 'phgt': 'FLOAT'}
    for k, v in type_map.items():
        pd = ipt.Pdata(inputs.plist, inputs.slc_par,inputs[k], dtype=v)
        pd.to_csv(outputs[k], take=10)
        # pd_list = pd.to_location_list()
        # print(pd)
        # with open(outputs[k], 'w+') as of:
        #     writer = csv.writer(of, delimiter=',')
        #     writer.writerow(['ridx', 'azidx','x', 'y'] + ['record_{n}'.format(n=n) for n in range(pd.nrecords)])
        #     for pt in pd_list[::5]:
        #         writer.writerow(pt)


prepare_inputs(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
              snakemake.wildcards)