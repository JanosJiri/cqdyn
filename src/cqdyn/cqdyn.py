#!/usr/bin/env python3
"""
Code for quantum dynamics based on coefficients and precomputed matrices.
© Jiri Janos 2025
"""
import numpy as np


### functions ###
def output_logger():
    import logging
    """Set up the logger for output messages."""
    global logger

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    # handlers
    console_handler = logging.StreamHandler()
    file_handler = logging.FileHandler('cqdyn.out', mode='w')

    for log in [console_handler, file_handler]:
        log.setFormatter(logging.Formatter('%(message)s'))
        logger.addHandler(log)
    return logger


def read_input(input_file='input.json'):
    """Function to read input and initialize variables"""
    from json import load
    from os.path import exists

    # function for input reading
    def is_positive_integer(var, name):
        if not (isinstance(var, int) and var > 0):
            exit(f'Error: "{name}" is not positive integer ({var}).')

    def is_positive_float(var, name):
        if not (isinstance(var, float) and var > 0):
            exit(f'Error: "{name}" is not positive float ({var}).')

    def is_list(var, name):
        if not isinstance(var, list):
            exit(f'Error: "{name}" is not list.')

    def is_string(var, name):
        if not isinstance(var, str):
            exit(f'Error: "{name}" is not string.')

    def check_input_key(key, dict):
        if not key in dict:
            exit(f"Error: '{key}' not found in the input.")

    # opening input file
    if not exists(input_file):
        exit(f"Error: Input file '{input_file}' does not exist.")
    with open(input_file, 'r') as file:
        input = load(file)  # using json function load

    # reading times
    check_input_key('total_time', input)
    tot_time = float(input['total_time'])
    is_positive_float(tot_time, 'total_time')

    check_input_key('dt', input)
    dt = float(input['dt'])
    is_positive_float(dt, 'dt')

    check_input_key('print_time', input)
    print_time = float(input['print_time'])
    is_positive_float(print_time, 'print_time')

    # reading coefficients
    is_list(input['coefficients'], 'coefficients')
    c = np.array([complex(i) for i in input['coefficients']], dtype=complex)

    # determining number of states
    nstates = len(input['coefficients'])

    # reading Hamiltonian
    is_list(input['H_0'], 'H_0')
    H_0 = np.array(input['H_0'], dtype=float)
    if (nstates, nstates) != np.shape(H_0):
        exit(f"Error: H_0 matrix dimension {np.shape(H_0)} is not matching the number of "
             f"states/coefficients ({nstates}).")

    # reading interaction Hamiltonian V_int
    if 'V_int' in input:
        is_list(input['V_int'], 'V_int')
        V_int = np.array(input['V_int'], dtype=float)
    else:
        V_int = np.zeros(shape=(nstates, nstates))

    # reading field
    if 'field' in input:
        is_string(input['field'], 'field')
        field = input['field']
    else:
        field = "0"

    return tot_time, dt, print_time, nstates, c, H_0, V_int, field


def is_hermitean(H, name):
    """Checking that the inserted matrix is Hermitean. Since H is considered to be real, we check that it is symmetric.
    Exact symmetry is checked, no numerical tolerance is considered currently."""
    if np.allclose(H, H.T, rtol=0, atol=0):
        logger.info(f"* {name} is Hermitean (symmetric).")
    else:
        exit(f"Error: input {name} is not Hermitean (symmetric)!")


def norm(coefs):
    """Calculating norm from coefficients."""
    norm = 0
    for c in coefs:
        norm += np.real(np.conjugate(c)*c)

    if np.abs(norm - 1) > 1e-8:
        exit(f"Error: Norm ({norm:.10f}) is not conserved!")
    return norm


def cal_energy(c, H):
    """Calculating energy."""
    return np.vdot(c, H@c).real


def build_hamiltonian(H_0, V_int, field, t, dt):
    "Building Hamiltonian from H_0 and V_int."
    return H_0 + V_int*interaction(field, t + dt/2)


def interaction(field: str, t: float):
    """Calculating electric field."""
    # noinspection PyUnresolvedReferences
    from numpy import sin, cos, exp
    return eval(field)


def build_propagator(H, dt, hbar=1.0):
    from scipy.linalg import expm
    U = expm(-1j*H*dt/hbar)
    return U


def propagate(c, H, dt):
    """Propagation of coefficient."""
    U = build_propagator(H, dt)
    cnew = U@c  # @ is matrix multiplication
    return cnew


def plot(t, c, E):
    import matplotlib.pyplot as plt
    plt.rcParams["font.family"] = 'Helvetica'

    fig, axs = plt.subplots(2, 1, figsize=(6, 4.5), height_ratios=[1, 0.5], sharex='all')

    axs[0].set_title('Quantum Dynamics', fontsize=14)
    for i in range(len(c)):
        axs[0].plot(t, np.abs(c[i])**2, label=f'$|c_{{{i + 1}}}|^2$')
    axs[0].set_ylabel('Population')
    axs[0].set_xlim(t.min(), t.max())
    axs[0].set_ylim(-0.01, 1.01)
    axs[0].legend(labelspacing=0.0, loc='upper right')

    axs[1].plot(t, E, label='Energy', color='black')
    axs[1].set_xlabel('Time (a.t.u.)')
    axs[1].set_ylabel('Energy (a.u.)')

    for ax in axs:
        ax.minorticks_on()
        ax.tick_params(which='both', direction='in', top=True, right=True)

    plt.tight_layout()
    plt.subplots_adjust(hspace=0)
    plt.savefig('cqdyn.png', dpi=300)
    plt.show()


def coef_analysis(c):
    """Analyzing coefficients."""
    logger.info("* state  |  max population  |  min population  |  final population \n"
                "  ---------------------------------------------------------------- ")
    pop = np.abs(c)**2
    for state in range(np.shape(c)[0]):
        logger.info(f"  {state + 1:^5d}  |  {np.max(pop[state]):^14.6f}  |  {np.min(pop[state]):^14.6f}  "
                    f"|  {pop[state,-1]:^16.6f}")


def main():
    # setting up output logger
    logger = output_logger()
    logger.info("\n   #####################\n"
                "   ###     cQDyn     ###\n"
                "   #####################\n")

    ### reading input ###
    tot_time, dt, print_time, nstates, c, H_0, V_int, field = read_input()

    # printing input
    logger.info("Input:")
    logger.info(f"* total time: {tot_time}\n* dt: {dt}\n* print time: {print_time}\n* nstates: {nstates}"
                f"\n* coefficients:\n  " + np.array2string(c, separator='  ', formatter={
        'complex_kind': lambda x: f'{x.real:.5f}{x.imag:+.5f}j'}) + "\n* H_0 (Hamiltonian):")
    for row in H_0:
        logger.info("  " + np.array2string(row, separator='  ', formatter={'float': lambda x: f'{x:.5f}'}))
    logger.info(f"* V_int (interaction Hamiltonian):")
    for row in V_int:
        logger.info("  " + np.array2string(row, separator='  ', formatter={'float': lambda x: f'{x:.5f}'}))
    logger.info(f"* field: {field}")

    ### initialization ###
    logger.info("\nInitialization:")

    # check norm of the wave function
    norm(c)

    # check that H is Hermitean
    is_hermitean(H_0, 'H_0')
    is_hermitean(V_int, 'V_int')

    # if interaction is zero (V_int=0 or field=0), then substeps are not necessary as the exponential is the
    # exact propagator, hence, the time step is set to print time
    if np.all(V_int == 0):
        dt = print_time
        logger.info("* No time-dependent interacion (V_int = 0), hence, dt is set to print time."
                    f"\n  'dt' = 'print_time' = {dt:.5f}")
    elif field == "0":
        dt = print_time
        logger.info("* No time-dependent interacion (field = 0), hence, dt is set to print time."
                    f"\n  'dt' = 'print_time' = {dt:.5f}")
    elif dt > print_time:
        exit(f"* Time step (dt = {dt:.5f}) is larger than print time ({print_time:.5f}).")

    # preparing propagation variables
    nprint = int(print_time/dt)

    # calculate energy shifting
    eshift = np.average(np.diag(H_0))  # todo: employ it later in propagation

    ### propagation ###
    logger.info("\nPropagation:")
    t = 0
    H = build_hamiltonian(H_0, V_int, field, t, dt)
    energy = cal_energy(c, H)
    coef_arr, E_arr, time_arr = [c], [energy], [t]

    while t < tot_time:
        # building Hamiltonian
        H = build_hamiltonian(H_0, V_int, field, t, dt)

        # propagation
        c = propagate(c, H, dt)

        # check norm
        norm(c)

        # calculating energy
        energy = cal_energy(c, H)

        # adding dt to time
        t += dt

        # saving data
        if np.round(t/dt)%nprint == 0:
            coef_arr.append(c)
            E_arr.append(energy)
            time_arr.append(t)
            logger.info(f'* Time: {t:12.5f} a.t.u.; Energy: {energy: 12.5f} a.u.')

    # saving data
    coef_arr = np.array(coef_arr).T
    E_arr = np.array(E_arr)
    time_arr = np.array(time_arr)

    # writing data
    np.savetxt('coefficients.txt', coef_arr.T, fmt=' %10.5f',
        header='Coefficients of the wave function')  # todo: this can be deffinitely improved
    np.savetxt('energy.txt', np.column_stack((time_arr, E_arr)), fmt='%12.5f %20.5f',
        header='    Time                Energy')
    np.savez('cqdyn.npz', time=time_arr, coefficients=coef_arr, energy=E_arr)

    logger.info("\nPropagation finished. Data saved.")

    # analyzing results
    logger.info("\nAnalysis:")
    coef_analysis(coef_arr)

    # plotting
    logger.info("\nPlotting data")
    plot(time_arr, coef_arr, E_arr)


##### CODE #####
if __name__ == '__main__':
    main()
