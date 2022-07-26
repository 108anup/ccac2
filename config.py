from typing import Optional


class Config:
    C: float = 1
    R: float = 1
    D: float = 1
    T: int = 10
    F: int = 1
    compose: bool = True
    inf_buf: bool = True
    buf_size: Optional[float] = 1

    unsat_core: bool = False

    def check(self) -> bool:
        ''' Sanity check. Should assert this after setting all configs '''
        if self.C <= 0 or self.R <= 0 or self.D <= 0:
            return False
        if self.T <= 1 or self.F < 1:
            return False

        if self.inf_buf and self.buf_size is not None:
            return False
        if self.buf_size is not None and self.buf_size <= 0:
            return False

        return True
