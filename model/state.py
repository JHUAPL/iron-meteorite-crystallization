from dataclasses import dataclass, field


@dataclass
class SimulationState:
    # liquid fraction
    L: float = 1.0

    # step size
    f: float = 0.0005
    f0: float = 0.0005

    # bulk liquid concentrations
    LM: float = 0.0
    LMS: float = 0.0
    LMP: float = 0.0
    LMNi: float = 0.0

    # immiscibility flag
    twoliq: bool = False

    # results
    data: dict = field(default_factory=dict)