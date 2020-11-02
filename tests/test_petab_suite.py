"""Execute petab test suite."""

import importlib
import logging
import os
import petabtests
import pytest
from _pytest.outcomes import Skipped
from sbml2julia import core

importlib.reload(core)

try:
    import petab
except ImportError:
    pass

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def test_petab_suite():
    """Execute all cases from the petab test suite, report performance."""
    n_success = n_skipped = 0
    for case in petabtests.CASES_LIST:
        try:
            execute_case(case)
            n_success += 1
        except Skipped:
            n_skipped += 1
        except Exception as e:
            # run all despite failures
            logger.error(f"Case {case} failed.")
            logger.error(e)

        print('\n---------------------------------------------------------\n')


def execute_case(case):
    """Wrapper for _execute_case for handling test outcomes"""
    try:
        _execute_case(case)
    except Exception as e:
        if isinstance(e, NotImplementedError):
            logger.info(
                f"Case {case} expectedly failed. Required functionality is "
                f"not implemented: {e}")
            pytest.skip(str(e))
        else:
            raise e


def _execute_case(case):
    """Run a single PEtab test suite case"""
    case = petabtests.test_id_str(case)
    logger.info(f"Case {case}")

    # case folder
    case_dir = os.path.join(petabtests.CASES_DIR, case)

    # load solution
    solution = petabtests.load_solution(case)
    gt_chi2 = solution[petabtests.CHI2]
    gt_llh = solution[petabtests.LLH]
    gt_simulation_dfs = solution[petabtests.SIMULATION_DFS]
    tol_chi2 = solution[petabtests.TOL_CHI2]
    tol_llh = solution[petabtests.TOL_LLH]
    tol_simulations = solution[petabtests.TOL_SIMULATIONS]

    # import petab problem
    yaml_file = os.path.join(case_dir, petabtests.problem_yaml_name(case))

    # simulate
    t_steps = 11
    if case == '0006' or '0009' or '0010':
        t_steps = 3001
    problem = core.SBML2JuliaProblem(yaml_file, t_steps=t_steps, infer_ic_from_sbml=True)
    problem.write_jl_file()
    problem.optimize()
    problem.plot_results('c0')

    # extract results
    results = problem.results
    simulation_df = problem.petab_problem.simulation_df.rename(
        columns={petab.MEASUREMENT: petab.SIMULATION})
    print('simulation_df')
    print([simulation_df])
    print(gt_simulation_dfs)
    chi2 = results['chi2']
    llh = - results['fval']

    # check if matches
    chi2s_match = petabtests.evaluate_chi2(chi2, gt_chi2, tol_chi2)
    llhs_match = petabtests.evaluate_llh(llh, gt_llh, tol_llh)
    simulations_match = petabtests.evaluate_simulations(
        [simulation_df], gt_simulation_dfs, tol_simulations)

    # log matches
    logger.log(logging.INFO if chi2s_match else logging.ERROR,
               f"CHI2: simulated: {chi2}, expected: {gt_chi2},"
               f" match = {chi2s_match}")
    logger.log(logging.INFO if simulations_match else logging.ERROR,
               f"LLH: simulated: {llh}, expected: {gt_llh}, "
               f"match = {llhs_match}")
    logger.log(logging.INFO if simulations_match else logging.ERROR,
               f"Simulations: match = {simulations_match}")

    if not all([llhs_match, simulations_match, chi2s_match]):
        logger.error(f"Case {case} failed.")
        raise AssertionError(f"Case {case}: Test results do not match "
                             "expectations")

    logger.info(f"Case {case} passed.")


test_petab_suite()
