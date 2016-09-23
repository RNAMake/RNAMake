import math
import random

class MonteCarlo (object):
    def __init__(self):
        self.temperature = 1.0

    def accept(self, current, next):
        if next < current:
            return 1

        score = math.exp((current - next)/ self.temperature)
        if score > random.uniform(0.0, 1.0):
            return 1
        else:
            return 0
