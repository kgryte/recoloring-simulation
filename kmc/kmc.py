#!/usr/bin/env python
# coding: utf-8

"""Continuous-time kinetic Monte Carlo (KMC) implementation"""

import math
import sys
import os
import json
import numpy as np
import pandas as pd
import random as random

# Define exceptions:
class KineticMonteCarloError(Exception):
    """Base class for errors in the kinetic Monte Carlo (KMC) module"""

    pass

class InvalidInputError(KineticMonteCarloError):
    """Class for invalid input errors in the (KMC) module"""

    pass

class KineticMonteCarlo(object):
    """Class for kinetic Monte Carlo (KMC) instances

    Parameters:
        totalTime: total simulation time
        initialProbs: initial probabilities for the initial system state
        rateMatrix: rate matrix governing Markov chain evolution

    Instance parameters can be set using setter methods; e.g.,

        kmc = KineticMonteCarlo()
        kmc.totalTime = 1e6

    Alternatively, upon creating a new instance, keyword arguments can be used to configure settings; e.g.,

        kmc = KineticMonteCarlo(
            totalTime=1000,
            initialProbs=[0.5, 0.5],
            rateMatrix=[[-10, 10], [10, -10]]
        )

    Optional Parameters:
        seed: seed for the uniform random number generator (useful if wanting deterministic results)
        filename: output filename to which to save the generated Markov chain. A meta-data file will accompany the output.

    """

    def __init__(self, **keywords):
        """Initialize instance properties"""

        # Filename for output:
        self.__filename = None
        self.__saveFLG = False

        # Total simulation time:
        self.__totalTime = 1000

        # The initial probabilities for the Markov Chain:
        self.__initialProbs = np.array([0.5, 0.5])

        # The rate matrix governing the evolution of the Markov Chain:
        self.__rateMatrix = np.array([[-10, 10], [10, -10]])

        # Random number generator seed:
        self.__seed = random.randint(0, sys.maxint)

        # Parse any input arguments:
        for key, value in keywords.iteritems():
            setattr(self, key, value)

        # Instantiate a random number generator instance:
        self.__rand = random.Random()

        # Seed the generator:
        self.__rand.seed(self.__seed)

        # Initialize the Markov Chain:
        self.__markovChain = []

    def run(self):
        """Run the KMC simulation"""

        # Initialize the chain:
        self.__initState()

        # Propagate the chain:
        self.__generateChain()

        # Save to output:
        if self.__saveFLG:
            self.__output()

    def __initState(self):
        """Determine the Markov Chain's initial state and initial transition time."""

        # [1] Determine the initial state:

        # Draw a random number:
        rand = self.__rand.random()

        # Generate a cumulative distribution from the initial probabilities:
        cumProbs = self.__initialProbs.cumsum()

        # Determine where our random value resides in the cumulative distribution and assign as the current state:
        currentState = np.where(cumProbs >= rand)[0][0]

        # [2] Determine the initial transition time:

        # Get the kinetic rates for the initial state:
        rates = self.__rateMatrix[currentState, :]

        # Remove the negative entry (from the diagonal):
        rates[currentState] = 0

        # Draw a random number:
        rand = self.__rand.random()

        # Generate an exponentially distributed random number using inverse transform sampling:
        waitingTime = -1 * math.log(rand) / rates.sum()

        # Place the current state, absolute time, and waiting time in the Markov Chain array:
        self.__markovChain.append([0, currentState, waitingTime])

    def __generateChain(self):
        """Generate a Markov Chain"""

        # Initialize a timer:
        time = self.__markovChain[0][0] + self.__markovChain[0][2]

        # Get the current state:
        currentState = self.__markovChain[0][1]

        # [3] Generate a Markov Chain until we exceed the totalTime: (accounting for machine precision error)
        FLG = True
        while FLG:

            # Get the rates from the current state to all other states:
            rates = self.__rateMatrix[currentState, :]

            # Remove the current state entry (negative diagonal):
            rates[currentState] = 0

            # Generate the cumulative distribution:
            cumRates = rates.cumsum()

            # Draw a random number between 0 and the cumulative rate sum:
            rand = self.__rand.random() * cumRates[-1]

            # Determine where our random value resides in the cumulative distribution and assign as the current state:
            currentState = np.where(cumRates >= rand)[0][0]

            # Get the new rates:
            rates = self.__rateMatrix[currentState, :]
            rates[currentState] = 0

            # Draw a 'dwell' from an exponential distribution:
            rand = self.__rand.random()

            waitingTime = -1 * math.log(rand) / rates.sum()

            # Update our Markov Chain:
            self.__markovChain.append([time, currentState, waitingTime])

            # Update the timer:
            time = time + waitingTime

            # Check our break condition:
            if (time - self.__totalTime) >= sys.float_info.epsilon:
                FLG = False

    def __output(self):
        """Save the generated output to file.

        Output Markov chain is to CSV. Accompanied by JSON file of simulation parameters.
        """

        # If the user has provided a path, we need to ensure the path exists before attempting to save to file:
        path = os.path.split(self.__filename)
        if not os.path.isdir(path[0]):
            try:
                os.makedirs(path[0])
            except OSError:
                raise

        # Convert the Markov Chain to a data frame:
        markovChain = self.markovChain

        # Write the output to a CSV file:
        markovChain.to_csv(self.__filename, index=False)

        # Convert the parameters to JSON and save to file:
        params = {
            'desc': 'Simulation: Kinetic Monte Carlo. Output: Markov Chain.',
            'totalTime': self.__totalTime,
            'initialProbs': self.__initialProbs.tolist(),
            'rateMatrix': self.__rateMatrix.tolist(),
            'seed': self.__seed,
            'data': self.__filename
        }

        path = os.path.splitext(self.__filename)
        with open(os.path.join(path[0] + '.' + 'json'), 'w') as f:
            json.dump(params, f, indent=4)

    @property
    def totalTime(self):
        """Getter: totalTime"""
        return self.__totalTime

    @totalTime.setter
    def totalTime(self, totalTime):
        """Setter: totalTime

        totalTime must be a positive numeric value.
        """

        if not (isinstance(totalTime, int) or isinstance(totalTime, float)):
            raise InvalidInputError("totalTime must be of type 'int' or 'float'")
        if totalTime <= 0:
            raise InvalidInputError("totalTime be non-zero and positive")

        self.__totalTime = totalTime

    @property
    def initialProbs(self):
        """Getter: initialProbs

        Return the initial probabilities used to generate an initial state for the Markov chain. Returned value is a numpy array.
        """
        return self.__initialProbs

    @initialProbs.setter
    def initialProbs(self, probs):
        """Setter: initialProbs

        All probabilities must satisfy the following conditions:
            0 <= pi <= 1
            sum(pi) = 1

        initialProbs can either be a list or numpy array.
        """

        if not isinstance(probs, np.ndarray):
            # Try to convert to np.array:
            try:
                probs = np.array(probs)
            except:
                raise InvalidInputError("input must be capable of conversion to numpy array")

        if (probs < 0).any() or (probs > 1).any():
            raise InvalidInputError("probabilities must be greater than or equal to 0 and less than or equal to 1")
        if probs.sum() != 1:
            raise InvalidInputError("probabilities must sum to unity")

        self.__initialProbs = probs

    @property
    def rateMatrix(self):
        """Getter: rateMatrix

        Return the rate matrix governing Markov chain evolution. Returned value is a numpy array.
        """
        return self.__rateMatrix

    @rateMatrix.setter
    def rateMatrix(self, rateMatrix):
        """Setter: rateMatrix

        A rate matrix must satisfy the following conditions:
            sum(A, axis=1) = 0  --> all rows must sum to zero
            for all A_{ij} where i=j, A_{ij} <= 0 --> all entries along diagonal must be negative
            A must be NxN --> a square matrix

        rateMatrix can either be a list or numpy array.
        """

        if not isinstance(rateMatrix, np.ndarray):
            # Try to convert to np.array:
            try:
                rateMatrix = np.array(rateMatrix)
            except:
                raise InvalidInputError("input must be capable of conversion to numpy array")

        if (rateMatrix.sum(axis=1) != 0).any():
            raise InvalidInputError("matrix rows must sum to zero")
        if (rateMatrix.diagonal() > 0).any():
            raise InvalidInputError("matrix diagonal must contain all negative entries")
        if rateMatrix.shape[0] != rateMatrix.shape[1]:
            raise InvalidInputError("matrix must be square")

        self.__rateMatrix = rateMatrix

    @property
    def seed(self):
        """Getter: seed

        Return the seed for the random number generator. The seed can be subsequently stored and used to reproduce a Markov chain given identical parameter settings.
        """
        return self.__seed

    @seed.setter
    def seed(self, seed):
        """Random number generator seed.

        Useful if wanting to generate a deterministic sequence. Seed must be capable of being converted to a hash.
        """

        try:
            hash(seed)
        except TypeError:
            raise InvalidInputError("seed must be hashable")

        self.__seed = seed

        # Re-seed the generator:
        self.__rand.seed(seed)

    @property
    def filename(self):
        """Getter: filename

        Return the filename to which output is saved.
        """

        return self.__filename

    @filename.setter
    def filename(self, filename):
        """Filename to which output is saved.

        Output will be a CSV file. Accompanying the output will be a meta-data file containing simulation settings (in JSON format).
        """

        if not isinstance(filename, str):
            raise InvalidInputError("filename must be of type 'str'")

        self.__filename = filename
        self.__saveFLG = True

    @property
    def markovChain(self):
        """Getter: markovChain

        Return the simulated Markov chain. Output is a Pandas data frame.
        """

        # Output as a Pandas data frame:
        return pd.DataFrame(self.__markovChain, columns=['time', 'state', 'waiting_time'])

# Module version number:
__version__ = '0.1.0'

#if __name__ == '__main__':
#    pass
