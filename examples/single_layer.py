"""
Demonstrates the decomposition of yearly input of a unit of leave litter
together with the N enrichment over time
"""


import decomp
from math import cos, pi
import pylab as plt
import numpy as np
import argparse

class param:
    TMAX = 15
    TMIN = -3
    pH = 8


def date2doy(date):
    return (date - np.datetime64(date, 'Y')).astype(int)


def temperature(date):
    """
    A proxy for the yearly changes of temperature and wetness
    :param date: A datetime64 value
    :return: Temperature, wetness, pH
    """
    doy = date2doy(date)
    tavg = 0.5 * (param.TMAX + param.TMIN)
    tampl = 0.5 * (param.TMAX - param.TMIN)
    T = tavg + tampl * cos(doy / 365 * 2 * pi)
    wetness = 1 - T / param.TMAX
    return T, wetness, param.pH


def run(som, input_functions, doc_retention_time=0.0):
    """
    Runs the model for 19 years with daily timesteps
    :param som: Initial soil organic matter conditions
    :param input_functions: List of f(date) functions to generate input
    :param doc_retention_time:
    :return:
    """
    input_function = lambda date: sum((ifunc(date) for ifunc in input_functions), decomp.SOM())
    dates = np.arange(plt.datetime64('2000-01-01'), plt.datetime64('2019-01-01'))
    som_state = np.NaN * np.zeros((len(dates), 6))
    som_flux = np.NaN * np.zeros((len(dates), 6))
    N_flux = np.NaN * np.zeros(dates.shape)
    CN = np.NaN * np.zeros(dates.shape)
    doc_retention_rate = 1 - (1 / doc_retention_time) if doc_retention_time else 0.0
    for i, d in enumerate(dates):
        som += input_function(d)
        som_state[i] = [c for _, c in som]
        CN[i] = som.CN
        old_doc = som[decomp.DOC]
        flux = som.integrate(1, *temperature(d))
        som[decomp.DOC] = (old_doc + flux[decomp.DOC]) * doc_retention_rate
        som_flux[i] = [c for _, c in flux]
        N_flux[i] = flux.N

    return dates, som_state, som_flux, CN, N_flux


def state_plot(t, som_state, CN):
    plt.stackplot(t, som_state.T, labels=[n for n, _ in decomp.SOM()])
    plt.ylabel('C in SOM')
    plt.legend(loc=1)
    if CN is not None:
        plt.twinx()
        plt.plot(t, CN, 'k-', label='C/N ratio', lw=2)
        plt.ylabel('C/N ratio')
        plt.legend(loc=4)


def flux_plot(t, som_flux, N_flux):
    plt.stackplot(t, som_flux.T, labels=[n for n, _ in decomp.SOM()]) 
    plt.ylabel('C flux from SOM in $day^{-1}$')
    plt.legend()
    if N_flux is not None:
        plt.twinx()
        plt.plot(t, N_flux, 'k-', label='N leaching', lw=2)
        plt.ylabel('N leaching')
        plt.legend(loc=4)
        

def plot(t, som_state, som_flux, CN=None, N_flux=None):
    fig, ax = plt.subplots(nrows=2, sharex=True)
    plt.sca(ax[0])
    state_plot(t, som_state, CN)
    # flux plot
    plt.sca(ax[1])
    flux_plot(t, som_flux, N_flux)


def cli():
    som_dict = {
        'no': decomp.SOM(),
        'leave': decomp.leave_litter(),
        'wood': decomp.wood_litter(),
        'root': decomp.root_litter(),
        'doc': decomp.pure_DOC(),
    }

    def yearly_function(somname):
        return lambda date: som_dict[somname] if date2doy(date) == 270 else decomp.SOM()

    def daily_function(somname):
        return lambda date: som_dict[somname] / 365

    parser = argparse.ArgumentParser()
    parser.add_argument('--docretention', '-r', type=float, default=0.0)
    parser.add_argument('--yearly', '-y', choices=som_dict, nargs='*')
    parser.add_argument('--daily', '-d', choices=som_dict, nargs='*')
    parser.add_argument('initial', nargs='*')
    args = parser.parse_args()
    initial = sum((som_dict[a] for a in args.initial), decomp.SOM()) if args.initial else decomp.SOM()
    yearly = [yearly_function(a) for a in args.yearly] if args.yearly else []
    daily = [daily_function(a) for a in args.daily] if args.daily else []
    return initial, yearly + daily, args.docretention


if __name__ == '__main__':
    initial, input_functions, doc_retention_time = cli()
    result = run(initial, input_functions, doc_retention_time)
    plot(*result)
    plt.show()
