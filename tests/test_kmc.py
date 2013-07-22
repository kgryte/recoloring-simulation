#!/usr/bin/env python
# coding: utf-8

# Import various nose tools:
from nose.tools import eq_, raises

# Import the Kinetic Monte Carlo module:
from kmc import *

# Test Bad Input:
class TestBadInput:
    """Test Cases for Bad Input"""

    def setUp(self):
        self._kmc = KineticMonteCarlo()

    def tearDown(self):
        del self._kmc

    # totalTime:
    def test_totalTime_init(self):
        """totalTime should be initialized to 1000"""
        eq_(self._kmc.totalTime, 1000)

    def test_totalTime_sanity(self):
        """totalTime should accept positive numerical input"""
        self._kmc.totalTime = 900
        eq_(self._kmc.totalTime, 900)

    @raises(InvalidInputError)
    def test_totalTime_is_numeric(self):
        """totalTime should fail with non-numeric input"""
        self._kmc.totalTime = 'string'

    @raises(InvalidInputError)
    def test_totalTime_is_not_negative(self):
        """totalTime should fail with negative input"""
        self._kmc.totalTime = -1

    @raises(InvalidInputError)
    def test_totalTime_is_not_zero(self):
        """totalTime should fail with 0 as input"""
        self._kmc.totalTime = 0

    # initialProbs:
    def test_initialProbs_init(self):
        """initialProbs should be initialized to [0.5, 0.5]"""
        test = np.array_equal(self._kmc.initialProbs, np.array([0.5, 0.5]))
        eq_(test, True)

    def test_initialProbs_sanity_list_input(self):
        """initialProbs should accept a list as input"""
        self._kmc.initialProbs = [0.6, 0.4]
        test = np.array_equal(self._kmc.initialProbs, np.array([0.6, 0.4]))
        eq_(test, True)

    def test_initialProbs_sanity_numerical_input(self):
        """initialProbs should accept numerical array input"""
        self._kmc.initialProbs = np.array([0.6, 0.4])
        test = np.array_equal(self._kmc.initialProbs, np.array([0.6, 0.4]))
        eq_(test, True)

    @raises(InvalidInputError)
    def test_initialProbs_within_range(self):
        """initialProbs should fail if input less than 0 or greater than 1"""
        self._kmc.initialProbs = np.array([-1, 2])

    @raises(InvalidInputError)
    def test_initialProbs_sum_to_unity(self):
        """initialProbs should fail if sum does not equal unity"""
        self._kmc.initialProbs = np.array([0.2, 0.5])

    # rateMatrix:
    def test_rateMatrix_init(self):
        """rateMatrix should be initialized to [[-10,10], [10,-10]]"""
        test = np.array_equal(self._kmc.rateMatrix, np.array([[-10, 10], [10, -10]]))
        eq_(test, True)

    def test_rateMatrix_sanity_list_input(self):
        """rateMatrix should accept a list as input"""
        self._kmc.rateMatrix = [[-20, 20], [20, -20]]
        test = np.array_equal(self._kmc.rateMatrix, np.array([[-20, 20], [20, -20]]))
        eq_(test, True)

    def test_rateMatrix_sanity_numerical_input(self):
        """rateMatrix should accept numerical array input"""
        self._kmc.rateMatrix = np.array([[-20, 20], [20, -20]])
        test = np.array_equal(self._kmc.rateMatrix, np.array([[-20, 20], [20, -20]]))
        eq_(test, True)

    @raises(InvalidInputError)
    def test_rateMatrix_rows_sum_to_zero(self):
        """rateMatrix should fail if rows do not sum to 0"""
        self._kmc.rateMatrix = np.array([[-5, 10], [10, -10]])

    @raises(InvalidInputError)
    def test_rateMatrix_diag_should_be_negative(self):
        """rateMatrix should fail if diagonal is not all negative"""
        self._kmc.rateMatrix = np.array([[-10, 10], [-10, 10]])

    @raises(InvalidInputError)
    def test_rateMatrix_should_be_square(self):
        """rateMatrix should fail if not a square matrix"""
        self._kmc.rateMatrix = np.array([[-10, 10]])

    # seed:
    @raises(InvalidInputError)
    def test_seed_only_accepts_hashable_input(self):
        """seed should fail if input is not hashable"""
        self._kmc.seed = {}

    # Keyword arguments:
    def test_keyword_arguments(self):
        """Instantiation should allow for setting configuration using keyword arguments"""
        kmc = KineticMonteCarlo(totalTime=900, initialProbs=[0.4, 0.6], rateMatrix=[[-100, 100], [100, -100]])

        eq_(kmc.totalTime, 900)

        test = np.array_equal(kmc.rateMatrix, np.array([[-100, 100], [100, -100]]))
        eq_(test, True)

        test = np.array_equal(kmc.initialProbs, np.array([0.4, 0.6]))
        eq_(test, True)

    # Filename:
    def test_filename_sanity_accept_string_input(self):
        """filename should accept string input"""
        self._kmc.filename = 'kmc-results/output.csv'
        eq_(self._kmc.filename, 'kmc-results/output.csv')

    @raises(InvalidInputError)
    def test_filename_is_string(self):
        """filename should fail if not a string"""
        # check text output...
        self._kmc.filename = 5

# Test Simulation Behavior:
class TestSimulation:
    """Test cases for simulation behavior"""

    def setUp(self):
        self._kmc = KineticMonteCarlo()

    def tearDown(self):
        del self._kmc

    def test_markovChain_random_output(self):
        """markovChain is a non-deterministic stochastic process"""

        # Run our first simulation instance:
        self._kmc.run()

        # Create a new instance and simulate:
        kmc = KineticMonteCarlo()
        kmc.run()

        # The two chains should NOT be equal:
        test = np.array_equal(self._kmc.markovChain, kmc.markovChain)
        eq_(test, False)

    # seed:
    def test_seed_generates_deterministic_behavior(self):
        """seed should induce deterministic behavior"""

        # Run our first simulation instance:
        self._kmc.seed = 0
        self._kmc.run()

        # Create a new instance and simulate:
        kmc = KineticMonteCarlo()
        kmc.seed = 0
        kmc.run()

        # The two chains should be equal:
        test = np.array_equal(self._kmc.markovChain, kmc.markovChain)
        eq_(test, True)

    def test_seed_generates_reproducible_behavior(self):
        """seed should generate reproducible behavior"""

        # Run our first simulation instance:
        self._kmc.run()
        seed = self._kmc.seed

        # Create a new instance and simulate, using the seed of the previous instance:
        kmc = KineticMonteCarlo()
        kmc.seed = seed
        kmc.run()

        # The two chains should be equal:
        test = np.array_equal(self._kmc.markovChain, kmc.markovChain)
        eq_(test, True)

    # Save to output:
    def test_output_saves_to_csv(self):
        """if filename is provided, should save to CSV output"""

        self._kmc.filename = 'kmc-results/output.csv'
        self._kmc.run()

        # Check if file exists:
        try:
            with open(self._kmc.filename): pass
        except:
            raise

        # Read the file and check that it was written correctly:
        df = pd.read_csv(self._kmc.filename)

        # Note: we just check for the first column, as rounding errors may have occurred during output. Here, we just check that the state sequence is the same.
        test = np.array_equal(self._kmc.markovChain.values[:,1], df.values[:,1])
        eq_(test, True)

    def test_output_saves_settings_to_json(self):
        """if output, settings should save to JSON"""

        self._kmc.filename = 'kmc-results/output.csv'
        self._kmc.run()

        # Check if file exists:
        path = os.path.splitext(self._kmc.filename)
        filename = os.path.join(path[0] + '.' + 'json')
        try:
            with open(filename, 'r') as f:

                # Read the file and check that it was written correctly:
                data = f.read()
                settings = json.loads(data)

                kmc = KineticMonteCarlo()
                kmc.seed = settings['seed']
                kmc.run()

                # Simulation should be reproducible:
                test = np.array_equal(self._kmc.markovChain, kmc.markovChain)
                eq_(test, True)

        except:
            raise
