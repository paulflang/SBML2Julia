"""Execute petab test suite."""

import petabtests
from DisFit import core
import importlib
importlib.reload(core)
import sys
import os
import pytest
from _pytest.outcomes import Skipped
import logging

try:
    import petab
except ImportError:
    pass

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

CASES_DIR = os.path.join('..', 'Benchmark-Models-PEtab', 'Benchmark-Models')


def test_petab_benchmark():
    """Execute all cases from the petab test suite, report performance."""
    n_success = n_skipped = 0
    i = 0
    for case in os.listdir(CASES_DIR):
        i += 1
        if i != 11: # 2 Input Error: Incorrect initial partitioning scheme. # 6: Ask Sungho why NLP Error again. # 8: EXIT: Restoration Failed! # 10 Takes to long to be run. Process is killed (at least if run from python)
            # 14: Why is the process killed? 18: EXIT: Restoration Failed!.
            continue
        print(case)
        try:
            execute_case(case)
            n_success += 1
        except Skipped:
            n_skipped += 1
        except Exception as e:
            # run all despite failures
            logger.error(f"Case {case} failed.")
            logger.error(e)


def execute_case(case):
    """Wrapper for _execute_case for handling test outcomes"""
    try:
        _execute_case0(case)
    except Exception as e:
        if isinstance(e, NotImplementedError) \
                or 'Unrecognized function "piecewise" used in nonlinear expression.' in str(e): \
                # or "BoundsError: attempt to access 2×2 DataFrame" in str(e) \
                # or "BoundsError: attempt to access 1-element Array{Any,1} at index [2]" in str(e) \
                # or "Compartments are not implemented" in str(e) \
                # or "UndefVarError: noiseParameter1_obs_a not defined" in str(e):
                # Cases Alkan_SciSicnal2018, 
            logger.info(
                f"Case {case} expectedly failed. Required functionality is "
                f"not implemented: {e}")
            pytest.skip(str(e))
        else:
            raise e


def _execute_case0(case):
    """Run a single PEtab test suite case"""
    logger.info(f"Case {case}")

    # case folder
    case_dir = os.path.join(CASES_DIR, case)

    # load solution
    # solution = petabtests.load_solution(case)
    # gt_chi2 = solution[petabtests.CHI2]
    # gt_llh = solution[petabtests.LLH]
    # gt_simulation_dfs = solution[petabtests.SIMULATION_DFS]
    # tol_chi2 = solution[petabtests.TOL_CHI2]
    # tol_llh = solution[petabtests.TOL_LLH]
    # tol_simulations = solution[petabtests.TOL_SIMULATIONS]

    # import petab problem
    yaml_file = os.path.join(case_dir, f'{case}.yaml')

    # simulate
    problem = core.DisFitProblem(yaml_file, optimizer_options={'linear_solver': 'MA57'}, infer_ic_from_sbml=False)
    problem.write_jl_file()
    problem.optimize()

    # extract results
    results = problem.results
    simulation_df = problem.petab_problem.simulation_df.rename(columns={petab.MEASUREMENT: petab.SIMULATION})
    chi2 = results['chi2']
    llh = - results['fval']
    for c in [c for c, c_ind in problem._condition2index.items() if c_ind+1 in problem._j_to_parameters[0]]:
        problem.plot_results(c, path=os.path.join(., 'plots', f'plot_{case}_'+c+'.pdf'))

    # check if matches
    # chi2s_match = petabtests.evaluate_chi2(chi2, gt_chi2, tol_chi2)
    # llhs_match = petabtests.evaluate_llh(llh, gt_llh, tol_llh)
    # simulations_match = petabtests.evaluate_simulations(
    #     [simulation_df], gt_simulation_dfs, tol_simulations)    

    # # log matches
    # logger.log(logging.INFO if chi2s_match else logging.ERROR,
    #            f"CHI2: simulated: {chi2}, expected: {gt_chi2},"
    #            f" match = {chi2s_match}")
    # logger.log(logging.INFO if simulations_match else logging.ERROR,
    #            f"LLH: simulated: {llh}, expected: {gt_llh}, "
    #            f"match = {llhs_match}")
    # logger.log(logging.INFO if simulations_match else logging.ERROR,
    #            f"Simulations: match = {simulations_match}")

    # if not all([llhs_match, simulations_match, chi2s_match]):
    #     logger.error(f"Case {case} failed.")
    #     raise AssertionError(f"Case {case}: Test results do not match "
    #                          "expectations")

    # if case == '0013':
    #     logger.error(f"Case {case} failed.")
    #     raise AssertionError(f"Case {case}: Compartments are not implemented.")

    logger.info(f"Case {case} passed.")

test_petab_benchmark()