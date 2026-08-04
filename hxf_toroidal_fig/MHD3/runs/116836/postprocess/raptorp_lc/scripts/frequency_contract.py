"""Single source of truth for the Stage 1 frequency products."""

from __future__ import annotations

from decimal import Decimal, getcontext
from pathlib import Path

getcontext().prec = 50

SPEED_OF_LIGHT_M_S = Decimal("299792458")
PLANCK_J_S = Decimal("6.62607015e-34")
ELECTRON_VOLT_J = Decimal("1.602176634e-19")


def lc4_frequencies() -> list[tuple[str, Decimal]]:
    """Return exact requested labels and frequency values in Hz."""
    return [
        ("1.5 GHz", Decimal("1.500000000e9")),
        ("219 GHz", Decimal("2.190000000e11")),
        ("521 nm", SPEED_OF_LIGHT_M_S / Decimal("521e-9")),
        ("1.2 keV", Decimal("1.2e3") * ELECTRON_VOLT_J / PLANCK_J_S),
    ]


def sed50_frequencies() -> list[Decimal]:
    """Implement exactly numpy.logspace(9, 18, 50), without a numpy dependency."""
    return [Decimal(10) ** (Decimal(9) + Decimal(9) * Decimal(i) / Decimal(49))
            for i in range(50)]


def write_frequency_file(path: Path, values: list[Decimal]) -> None:
    path.write_text("".join(f"{value:.15E}\n" for value in values), encoding="ascii")
