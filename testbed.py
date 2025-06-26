import json

import matplotlib.pyplot as plt
import pandas as pd
import diffrax
import optax
import equinox as eqx
import jax.numpy as jnp
from jax import value_and_grad
from jax.tree_util import tree_map
from jax import config
from tqdm import tqdm

config.update("jax_enable_x64", True)
plt.style.use("classic")
plt.rcParams.update(
    {
        "figure.constrained_layout.use": True,
        "text.usetex": True,
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman"],
        "axes.titlesize": 10,
        "font.size": 10,
        "legend.fontsize": 9,
        "savefig.bbox": "tight",
    }
)


class Model(eqx.Module):
    kV: float
    aDV: float
    dV: float
    dBV: float
    dBE: float
    kBmax: float
    aB: float
    tB: float
    kED: float
    dE: float
    DT: float
    V01: float
    V02: float
    D0: float
    E0: float
    Btot: float
    kE: float

    def __call__(self, t, y, args):
        kB = lambda t: self.kBmax / (1 + jnp.exp(-self.aB * (t - self.tB)))
        dydt = {}
        dydt["V"] = self.kV * (1 + self.aDV * y["D"]) * y["V"] - self.dV * y["V"]
        dydt["D"] = (
            self.dBV * (self.Btot - y["D"]) * y["V"]
            + self.dBE * (self.Btot - y["D"]) * y["E"]
            - kB(t) * y["D"]
        )
        dydt["E"] = (
            self.kE
            - self.dE * y["E"]
            + self.kED * y["D"] * jnp.heaviside(y["threshold_reached"], 0.0)
        )
        dydt["threshold_reached"] = jnp.maximum(y["D"] - self.DT, 0.0)
        return dydt


class Simulate(eqx.Module):
    model: eqx.Module

    def __call__(self, days1, days2, treatment=True):
        y01 = {
            "D": self.model.D0,
            "E": self.model.E0,
            "threshold_reached": 0.0,
        }
        y01["V"] = self.model.V01 if treatment else 0.0
        infection1 = diffrax.diffeqsolve(
            terms=diffrax.ODETerm(self.model),
            solver=diffrax.Kvaerno5(scan_stages=True),
            t0=5,
            t1=21,
            dt0=None,
            y0=y01,
            saveat=diffrax.SaveAt(ts=days1),
            stepsize_controller=diffrax.PIDController(rtol=1.4e-8, atol=1e-8),
            adjoint=diffrax.BacksolveAdjoint(),
            max_steps=None,
        )
        day15 = jnp.argwhere(infection1.ts == 15)
        y02 = tree_map(lambda y: jnp.squeeze(y[day15]), infection1.ys)
        if treatment:
            y02["V"] += self.model.V02
        infection2 = diffrax.diffeqsolve(
            terms=diffrax.ODETerm(self.model),
            solver=diffrax.Kvaerno5(scan_stages=True),
            t0=15,
            t1=31,
            dt0=None,
            y0=y02,
            saveat=diffrax.SaveAt(ts=days2),
            stepsize_controller=diffrax.PIDController(rtol=1.4e-8, atol=1.4e-8),
            adjoint=diffrax.BacksolveAdjoint(),
            max_steps=None,
        )
        infection1.ys["log(V)"] = jnp.log(infection1.ys["V"])
        infection2.ys["log(V)"] = jnp.log(infection2.ys["V"])
        return {1: infection1.ys, 2: infection2.ys}


def jaxify(d):
    """Turn a dict of lists and scalars into a dict of jax ndarrays."""
    to_jax_ndarray = lambda x: jnp.asarray(x, dtype=float)
    is_list = lambda x: isinstance(x, list)
    return tree_map(to_jax_ndarray, d, is_leaf=is_list)


def dataloader(fname):
    data = jaxify(
        pd.read_csv(fname)
        .groupby(["Infection", "Day"], as_index=False)
        .median()
        .drop(columns=["N", "T"], index=[0, 4])
        .groupby("Infection")
        .apply(lambda x: x.drop("Infection", axis=1).to_dict("list"))
        .to_dict()
    )
    data[1]["log(V)"] = jnp.log(data[1]["V"])
    data[2]["log(V)"] = jnp.log(data[2]["V"])
    return data


def plot(simulate, data, fname):
    days = {1: jnp.linspace(5, 21, num=201), 2: jnp.linspace(15, 31, num=201)}
    pred = {
        "control": simulate(days[1], days[2], treatment=False),
        "treatment": simulate(days[1], days[2], treatment=True),
    }
    fig, ax = plt.subplots(ncols=2, nrows=3, sharex=True, sharey="row")
    for i, state in enumerate(["log(V)", "D", "E"]):
        for j, regime in enumerate(["control", "treatment"]):
            for infection in [1, 2]:
                ax[i, j].plot(
                    data[regime][infection]["Day"], data[regime][infection][state], ".k"
                )
                ax[i, j].plot(days[infection], pred[regime][infection][state], "-k")
            ax[i, j].set_xlim([0, 31])
            if state == "D":
                ax[i, j].hlines(
                    simulate.model.DT, 0, days[2][-1], colors="k", linestyles="dashed"
                )
            if j == 0:
                ax[i, j].set_ylabel(state)
            if i == 2:
                ax[i, j].set_xlabel("Age of mice (days)")
    fig.savefig(fname, bbox_inches="tight")
    print("Saved: " + fname)


@value_and_grad
def loss_fn(log_simulate, data):
    simulate = tree_map(jnp.exp, log_simulate)
    pred = {
        "control": simulate(
            data["control"][1]["Day"], data["control"][2]["Day"], treatment=False
        ),
        "treatment": simulate(
            data["treatment"][1]["Day"], data["treatment"][2]["Day"], treatment=True
        ),
    }
    mse = lambda pred, data, i: jnp.mean(
        (pred[i]["V"] - data[i]["V"]) ** 2
        + (pred[i]["D"] - data[i]["D"]) ** 2
        + (pred[i]["E"] - data[i]["E"]) ** 2
    )
    return (
        mse(pred["treatment"], data["treatment"], 1)
        + mse(pred["treatment"], data["treatment"], 2)
        + mse(pred["control"], data["control"], 1)
        + mse(pred["control"], data["control"], 2)
    )


def fit(simulate, data, learning_rate, nsteps):
    log_simulate = tree_map(jnp.log, simulate)
    optimizer = optax.adam(learning_rate)
    opt_state = optimizer.init(log_simulate)
    for _ in (pbar := tqdm(range(nsteps))):
        loss_value, grads = loss_fn(log_simulate, data)
        updates, opt_state = optimizer.update(grads, opt_state)
        log_simulate = optax.apply_updates(log_simulate, updates)
        pbar.set_description(f"Loss: {loss_value:.0f}")
    simulate = tree_map(jnp.exp, log_simulate)
    return simulate.model.__dict__, loss_value


def save_params(params, fname):
    """Pretty print parameters dict to .json file."""
    extract_scalar = lambda x: x.item()
    is_jax_ndarray = lambda x: isinstance(x, jnp.ndarray)
    params_list = tree_map(extract_scalar, params, is_leaf=is_jax_ndarray)
    json.dump(params_list, open(fname, "w"), sort_keys=True, indent=4)
    print("Saved: " + fname)


def main():
    params = json.load(open("python/params_27804.json", "r"))
    simulate = Simulate(Model(**params))
    data = {
        "treatment": dataloader("data/raw_treatment_data.csv"),
        "control": dataloader("data/raw_control_data.csv"),
    }
    new_params, loss = fit(simulate, data, learning_rate=1e-4, nsteps=5000)
    plot(Simulate(Model(**new_params)), data, f"python/testbed_{loss:.0f}.pdf")
    save_params(new_params, f"python/params_{loss:.0f}.json")


if __name__ == "__main__":
    main()
