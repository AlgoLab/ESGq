import numpy as np


def compute_dpsi(C1, C2, dpsi_path):
    fdpsi = open(dpsi_path, "w")
    PSIs1 = {}
    PSIs2 = {}
    for bn, psi_path in C1.items():
        PSIs1[bn] = {}
        for line in open(psi_path):
            eidx, w1, w2, w3, psi = line.strip("\n").split(" ")
            PSIs1[bn][eidx] = float(psi) if psi != "NaN" else np.inf
    for bn, psi_path in C2.items():
        PSIs2[bn] = {}
        for line in open(psi_path):
            eidx, w1, w2, w3, psi = line.strip("\n").split(" ")
            PSIs2[bn][eidx] = float(psi) if psi != "NaN" else np.inf

    # FIXME: assuming same events for all psi
    # shared_events = set.intersection([])
    print(
        "Event",
        " ".join(C1.keys()),
        " ".join(C2.keys()),
        "dPSI",
        sep="\t",
        file=fdpsi,
    )
    for k in PSIs2[bn]:
        psi1 = []
        psi2 = []
        for bn in C1:
            psi1.append(PSIs1[bn][k])
        for bn in C2:
            psi2.append(PSIs2[bn][k])
        dpsi = "NaN"
        diff1 = np.mean(psi1)
        diff2 = np.mean(psi2)
        if not (
            diff1 == np.inf
            or diff2 == np.inf
            or diff1 - diff2 == np.nan
            or diff1 - diff2 == np.inf
        ):
            dpsi = abs(diff1) - abs(diff2)
        print(
            k,
            " ".join([str(round(x, 3)) for x in psi1]),
            " ".join([str(round(x, 3)) for x in psi2]),
            dpsi,
            sep="\t",
            file=fdpsi,
        )
    fdpsi.close()


if __name__ == "__main__":
    pass
