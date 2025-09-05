from enum import Enum


class BkgPdfFamily(Enum):
    CHEBYCHEV = "chebychev"
    POWER_LAW = "power_law"
    EXPONENTIAL = "exponential"
    BERNSTEIN = "bernstein"

    def __str__(self):
        return self.value
