import numpy as np 
import pandas as pd
import pandapower as pp 
import pandapower.plotting
import pandapower.networks
import networkx as nx 
import random

BASE_KV = 0.4
BASE_KVA = 10
BASE_Z = (BASE_KV*1e3)**2/(BASE_KVA*1e3)

def replace_trafos(net):
    index = list(net.trafo.index)
    from_bus = list(net.trafo.hv_bus)
    to_bus = list(net.trafo.lv_bus)
    for idx in index:
        pp.drop_trafos(net, [idx])
        net.bus.drop(from_bus, inplace=True)
        net.ext_grid.drop(from_bus, inplace=True)
        pp.create_ext_grid(net, to_bus[0], vm_pu=1.0, va_degree=0.0)
        if not net.bus_geodata.empty:
            net.bus_geodata.drop(from_bus, inplace=True)
        pp.create_continuous_elements_index(net, start=0)

    return net


def clean_lines(net):
    for idx in net.line.index:
        net.line.c_nf_per_km[idx] = 0
        net.line.max_i_ka[idx] = 1000
        # Maybe remove i_max
    return net


def single_volt_lvl(net, voltage):
    for idx in net.bus.index:
        net.bus.v_kv = voltage
    return net


def net_to_graph(net):
    G = nx.DiGraph()

    for idx in net.bus.index:
        G.add_node(idx, P_cons=0, Q_cons=0, S_max=0, P_max=0)

    for l in net.line.index:
        line = net.line
        G.add_edge(line.from_bus[l], line.to_bus[l], r=line.r_ohm_per_km[l]*line.length_km[l]/BASE_Z,
                   x=line.x_ohm_per_km[l]*line.length_km[l]/BASE_Z, closed=True)

    return G


def add_PV(net, G, num_pv, sives_pv):
    # for reproducibility
    random.seed(23432, version=2)

    # set pv systems at random buses
    for idx, name in enumerate(net.bus.name):
        if not name:
            net.bus.name[idx] = 'bus '+str(idx)
    potential_site = list(net.bus.name)[1:]
    pv_buses = random.sample(potential_site, num_pv)
    pv_bus_idx = pp.get_element_indices(net, 'bus', pv_buses, exact_match=True)
    # init active power
    p_mw = np.zeros(len(pv_buses))
    # set rated sizes
    sn_mva = [random.choice(sives_pv) for i in range(num_pv)]
    # name sgen
    name = [str('PV-')+bus_name for bus_name in pv_buses]
    # print information
    print("PV systems and rated sizes")
    for idx, nm in enumerate(name):
        print("PV %d: name: %s, inverter size: %.1f kVA" %
              (idx, nm, sn_mva[idx]*1000))
    # create sgen
    if len(net.sgen.name):
        net.sgen.drop(net.sgen.index, inplace=True)
    for idx, bus in enumerate(pv_bus_idx):
        pp.create_sgen(net, bus, p_mw=p_mw[idx], q_mvar=0, sn_mva=sn_mva[idx], name=name[idx],
                       index=None, scaling=1.0, type='wye', in_service=True, controllable=True)
        G.nodes[bus]['S_max'] = sn_mva[idx]*1e3

    # add sgen at slack bus
    pp.create_sgen(net, 0, p_mw=0, q_mvar=0, sn_mva=1, name='Slack Gen',
                   index=None, scaling=1.0, type='wye', in_service=True, controllable=True)

    G.nodes[0]['S_max'] = 10000

    return net, G


def init(net, num_pv, sives_pv):
    net = replace_trafos(net)
    net = clean_lines(net)
    net = single_volt_lvl(net, BASE_KV)
    G = net_to_graph(net)
    net, G = add_PV(net, G, num_pv, sives_pv)
    return G, net


def import_data(net):
    arr_load = import_load(net)
    arr_PV = import_PV(net)
    print('Saving network in html file')
    pandapower.plotting.plotly.simple_plotly(net).write_html('network.html')
    return arr_PV, arr_load


def import_load(net):  # in kW
    data = pd.read_csv('data/LoadProfiles.csv', sep=",", encoding='utf-8')
    data_react = pd.read_csv(
        'data/LoadProfilesReact.csv', sep=",", encoding='utf-8')
    num_load = len(list(net.load.name))
    arr_load = np.zeros((60*60*24, num_load), dtype=complex)
    for i in range(0, num_load):
        profile_num = random.choice(np.arange(0, 39))
        arr_load[:, i] = np.array(data['profile_%d' %
                                       (profile_num)].to_list())+1j*np.array(data_react['profile_%d' %
                                                                                        (profile_num)].to_list())
    return arr_load


def import_PV(net):  # in kW
    data = pd.read_csv('data/PVProfiles.csv', sep=",", encoding='utf-8')
    num_pv = len(list(net.sgen.name))-1
    arr_PV = np.zeros((60*60*24, num_pv))
    for i in range(0, num_pv):
        prof = np.array(data['profile_%d' %
                             (random.choice(np.arange(0, 39)))].to_list())
        arr_PV[:, i] = prof/max(prof)*(net.sgen.sn_mva[i]*1000)
    return arr_PV

PV_penetration = 0.9

net = pandapower.networks.create_dickert_lv_network(feeders_range='short', linetype='cable', customer='multiple', case='average')

G, net = init(net, int(PV_penetration*len(net.bus)),[0.008, 0.009, 0.01]) 

#pp.to_json(net, "example.json", store_index_names=True)
pp.to_excel(net, "example_short.xlsx")



