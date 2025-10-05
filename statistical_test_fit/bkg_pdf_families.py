from enum import Enum


class BkgPdfFamily(Enum):
    JOHNSON = "johnson"
    BERNSTEIN = "bernstein"
    CHEBYCHEV = "chebychev"
    POWER_LAW = "power_law"
    EXPONENTIAL = "exponential"

    def __str__(self):
        return self.value.capitalize().replace("_", " ")
