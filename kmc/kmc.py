"""Continuous-time kinetic Monte Carlo (KMC) implementation"""

import math
import sys
import numpy as np
import random as random

# Define exceptions:
class KineticMonteCarloError(Exception):
    """Base class for errors in the kinetic Monte Carlo (KMC) module"""

    pass

class InvalidInputError(KineticMonteCarloError):
    """Class for invalid input errors in the (KMC) module"""

    pass

class KineticMonteCarlo(object):
    """Class for kinetic Monte Carlo (KMC) instances"""

    def __init__(self):
        """Initialize instance properties"""

        # Total simulation time:
        self.__totalTime = 1000

        # The initial probabilities for the Markov Chain:
        self.__initialProbs = np.array([0.5, 0.5])

        # The rate matrix governing the evolution of the Markov Chain:
        self.__rateMatrix = np.array([[-10, 10], [10, -10]])

        # Instantiate a random number generator instance:
        self.__rand = random.Random()

        # Initialize the Markov Chain:
        self.__markovChain = []

    def run(self):
        """Run the KMC simulation"""

        # Initialize the chain:
        self.__initState()

        # Propagate the chain:
        self.__generateChain()

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
        self.__markovChain.append([currentState, 0+waitingTime, waitingTime])

    def __generateChain(self):
        """Generate a Markov Chain"""

        # Initialize a timer:
        time = self.__markovChain[0][1]

        # Get the current state:
        currentState = self.__markovChain[0][0]

        # [3] Generate a Markov Chain until we exceed the totalTime: (accounting for machine precision error)
        while (time - self.__totalTime) < sys.float_info.epsilon:

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
            self.__markovChain.append([currentState, time+waitingTime, waitingTime])

            # Update the timer:
            time = time + waitingTime

        #np.savetxt('test.txt', self.markovChain)

    @property
    def totalTime(self):
        """Getter: totalTime"""
        return self.__totalTime

    @totalTime.setter
    def totalTime(self, totalTime):
        """Setter: totalTime"""

        if not (isinstance(totalTime, int) or isinstance(totalTime, float)):
            raise InvalidInputError("totalTime must be of type 'int' or 'float'")
        if totalTime <= 0:
            raise InvalidInputError("totalTime be non-zero and positive")

        self.__totalTime = totalTime

    @property
    def initialProbs(self):
        """Getter: initialProbs"""
        return self.__initialProbs

    @initialProbs.setter
    def initialProbs(self, probs):
        """Setter: initialProbs"""

        if (probs < 0).any() or (probs > 1).any():
            raise InvalidInputError("probabilities must be greater than or equal to 0 and less than or equal to 1")
        if probs.sum() != 1:
            raise InvalidInputError("probabilities must sum to unity")

        self.__initialProbs = probs

    @property
    def rateMatrix(self):
        """Getter: rateMatrix"""
        return self.__rateMatrix

    @rateMatrix.setter
    def rateMatrix(self, rateMatrix):
        """Setter: rateMatrix"""

        if (rateMatrix.sum(axis=1) != 0).any():
            raise InvalidInputError("matrix rows must sum to zero")
        if (rateMatrix.diagonal() > 0).any():
            raise InvalidInputError("matrix diagonal must contain all negative entries")
        if rateMatrix.shape[0] != rateMatrix.shape[1]:
            raise InvalidInputError("matrix must be square")

        self.__rateMatrix = rateMatrix

    def seed(self, seed):
        """Random number generator seed.

        If a seed is supplied, seed the generator. Useful if wanting to generate a deterministic sequence.
        """

        self.__rand.seed(seed)

    @property
    def markovChain(self):
        """Getter: markovChain"""

        # Output as a numpy array:
        return np.array(self.__markovChain)
