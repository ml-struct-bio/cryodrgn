import numpy as np


def get_beta_schedule(schedule):
    if type(schedule) == float:
        return ConstantSchedule(schedule)
    elif schedule == "a":
        return LinearSchedule(0.001, 15, 0, 1000000)
    elif schedule == "b":
        return LinearSchedule(5, 15, 200000, 800000)
    elif schedule == "c":
        return LinearSchedule(5, 18, 200000, 800000)
    elif schedule == "d":
        return LinearSchedule(5, 18, 1000000, 5000000)
    else:
        raise RuntimeError("Wrong beta schedule. Schedule={}".format(schedule))


class ConstantSchedule:
    def __init__(self, value):
        self.value = value

    def __call__(self, x):
        return self.value


class LinearSchedule:
    def __init__(self, start_y, end_y, start_x, end_x):
        self.min_y = min(start_y, end_y)
        self.max_y = max(start_y, end_y)
        self.start_x = start_x
        self.start_y = start_y
        self.coef = (end_y - start_y) / (end_x - start_x)

    def __call__(self, x):
        return np.clip(
            (x - self.start_x) * self.coef + self.start_y, self.min_y, self.max_y
        ).item(0)
