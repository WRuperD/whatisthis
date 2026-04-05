"""黑板镜面反射与最优俯仰角 — 数值仿真与作图（修正入射方向约定）"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib import rcParams

rcParams["font.sans-serif"] = ["SimHei", "Microsoft YaHei", "DejaVu Sans"]
rcParams["axes.unicode_minus"] = False

ROOT = Path(__file__).resolve().parent.parent
FIG_DIR = ROOT / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)


class Params:
    z_bot = 1.0
    H = 1.2
    y_board_min = -2.0
    y_board_max = 2.0
    x_seat_min = 2.0
    x_seat_max = 8.0
    y_seat_min = -2.5
    y_seat_max = 2.5
    z_eye = 1.15
    epsilon_deg = 6.0
    theta_max_deg = 6.0
    n_u = 8
    n_v = 6
    nx_grid = 28
    ny_grid = 24


P = Params()


def board_frame(theta_deg: float):
    th = np.radians(theta_deg)
    H = P.H
    T = np.array([H * np.sin(th), 0.0, H * np.cos(th)], dtype=float)
    Wvec = np.array([0.0, P.y_board_max - P.y_board_min, 0.0])
    n = np.cross(Wvec, T)
    n = n / np.linalg.norm(n)
    u_hat = T / np.linalg.norm(T)
    O = np.array([0.0, P.y_board_min, P.z_bot])
    return n, O, u_hat, H


def sample_board_points(theta_deg: float):
    n, O, u_hat, H = board_frame(theta_deg)
    pts = []
    for ui in np.linspace(0.0, 1.0, P.n_u):
        for vi in np.linspace(0.0, 1.0, P.n_v):
            p = O + ui * H * u_hat + vi * (P.y_board_max - P.y_board_min) * np.array([0.0, 1.0, 0.0])
            pts.append(p)
    return np.array(pts), n


def reflect(incident_to_surface: np.ndarray, n_out: np.ndarray) -> np.ndarray:
    """入射方向指向板面，外法向 n_out 指向教室（+X 侧）；满足 dot(l,n_out)<0。"""
    l = incident_to_surface / np.linalg.norm(incident_to_surface)
    r = l - 2.0 * np.dot(l, n_out) * n_out
    return r / np.linalg.norm(r)


def sun_incident_direction(t_hours: float) -> np.ndarray:
    """
    单侧窗平行光：传播方向指向黑板（需 l·n<0，n≈+X），随时刻水平方位变化。
    """
    frac = np.clip((t_hours - 8.0) / 9.0, 0.0, 1.0)
    az = np.radians(15.0 + 130.0 * frac)
    el = np.radians(28.0 + 10.0 * frac)
    cx, sx = np.cos(az), np.sin(az)
    ce, se = np.cos(el), np.sin(el)
    prop = np.array([-ce * cx, ce * sx, -se])
    prop = prop / np.linalg.norm(prop)
    return prop


def seat_glare(eye: np.ndarray, theta_deg: float, l_in: np.ndarray) -> bool:
    pts, n = sample_board_points(theta_deg)
    eps = np.radians(P.epsilon_deg)
    if np.dot(l_in, n) >= 0:
        return False
    for p in pts:
        r = reflect(l_in, n)
        w = eye - p
        wn = np.linalg.norm(w)
        if wn < 1e-6:
            continue
        v = w / wn
        c = float(np.clip(np.dot(r, v), -1.0, 1.0))
        if float(np.arccos(c)) < eps:
            return True
    return False


def min_theta_for_seat(eye: np.ndarray, l_in: np.ndarray):
    for th in np.arange(0.0, P.theta_max_deg + 1e-9, 0.05):
        if not seat_glare(eye, float(th), l_in):
            return float(th), True
    return float(P.theta_max_deg), False


def affected_mask(theta_deg: float, l_in: np.ndarray):
    xs = np.linspace(P.x_seat_min, P.x_seat_max, P.nx_grid)
    ys = np.linspace(P.y_seat_min, P.y_seat_max, P.ny_grid)
    Xc, Yc = np.meshgrid(xs, ys, indexing="ij")
    mask = np.zeros_like(Xc, dtype=bool)
    for i in range(P.nx_grid):
        for j in range(P.ny_grid):
            eye = np.array([Xc[i, j], Yc[i, j], P.z_eye])
            mask[i, j] = seat_glare(eye, theta_deg, l_in)
    return Xc, Yc, mask


def plot_topdown_compare(l_in: np.ndarray, theta0: float, theta1: float, tag: str):
    _, _, m0 = affected_mask(theta0, l_in)
    _, _, m1 = affected_mask(theta1, l_in)
    extent = [P.y_seat_min, P.y_seat_max, P.x_seat_min, P.x_seat_max]
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.6), sharey=True)
    im = None
    for ax, m, title in zip(
        axes,
        [m0, m1],
        [rf"调整前 $\theta={theta0:.2f}^\circ$", rf"调整后 $\theta={theta1:.2f}^\circ$"],
    ):
        im = ax.imshow(
            m.T.astype(float),
            origin="lower",
            extent=extent,
            aspect="auto",
            cmap="RdYlGn_r",
            vmin=0,
            vmax=1,
            interpolation="nearest",
        )
        ax.set_xlabel("沿黑板方向 $y$ / m")
        ax.set_ylabel("距黑板深度 $x$ / m")
        ax.set_title(title)
        ax.add_patch(
            Rectangle(
                (P.y_board_min, -0.15),
                P.y_board_max - P.y_board_min,
                0.3,
                facecolor="gray",
                edgecolor="k",
            )
        )
        ax.text(0.5 * (P.y_board_min + P.y_board_max), -0.55, "黑板（俯视投影）", ha="center", va="top", fontsize=9)
        ax.set_ylim(-0.8, P.x_seat_max + 0.5)
    fig.colorbar(im, ax=axes.ravel().tolist(), label="1=易受反射眩光影响", shrink=0.85)
    fig.suptitle(f"教室内受影响座位区域对比（{tag}）", fontsize=12)
    fig.subplots_adjust(left=0.06, right=0.98, top=0.88, bottom=0.14, wspace=0.25)
    out = FIG_DIR / f"topdown_compare_{tag}.png"
    fig.savefig(out, dpi=160)
    plt.close(fig)
    return out, float(m0.mean()), float(m1.mean())


def plot_theta_curve(times: np.ndarray, thetas: np.ndarray, feas: np.ndarray):
    fig, ax = plt.subplots(figsize=(7.5, 3.8))
    ax.plot(times, thetas, "b-", lw=1.5, label=r"$\theta^\ast(t)$")
    ax.axhline(P.theta_max_deg, color="r", ls="--", lw=1, label=r"$\theta_{\max}$")
    bad = ~feas
    if np.any(bad):
        ax.scatter(times[bad], thetas[bad], c="orange", s=22, zorder=5, label="上限内未完全消除（示意）")
    ax.set_xlabel("时刻 / h")
    ax.set_ylabel(r"俯仰角 / $^\circ$")
    ax.set_title("授课日逐时刻最小必要黑板俯仰角（前排代表点）")
    ax.legend(loc="upper left", fontsize=9)
    ax.set_xlim(times.min(), times.max())
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "theta_of_t.png", dpi=160)
    plt.close(fig)


def plot_section(theta_deg: float, l_in: np.ndarray, eye_xz: tuple[float, float], tag: str):
    th = np.radians(theta_deg)
    H = P.H
    z1 = P.z_bot + H * np.cos(th)
    x1 = H * np.sin(th)
    fig, ax = plt.subplots(figsize=(5.5, 4.2))
    ax.fill([0, x1, x1, 0], [P.z_bot, z1, z1, P.z_bot], color="#555", edgecolor="k", lw=1.5)
    ex, ez = eye_xz
    ax.plot(ex, ez, "o", color="C0", ms=8, label="代表眼睛")
    n, _, _, _ = board_frame(theta_deg)
    midx, midz = 0.5 * x1, 0.5 * (P.z_bot + z1)
    scale = 1.0
    ax.arrow(midx, midz, scale * l_in[0], scale * l_in[2], head_width=0.06, fc="gold", ec="k", length_includes_head=True)
    ax.arrow(midx, midz, 0.45 * n[0], 0.45 * n[2], head_width=0.05, fc="g", ec="k", length_includes_head=True)
    r = reflect(l_in, n)
    ax.arrow(midx, midz, scale * r[0], scale * r[2], head_width=0.06, fc="cyan", ec="k", length_includes_head=True)
    ax.set_xlabel("$x$ / m（距墙深度）")
    ax.set_ylabel("$z$ / m")
    ax.set_title(f"竖直截面示意（$\\theta={theta_deg:.2f}^\\circ$，{tag}）")
    ax.set_aspect("equal", adjustable="box")
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIG_DIR / f"section_{tag}.png", dpi=160)
    plt.close(fig)


def main():
    times = np.arange(8.0, 17.01, 0.5)
    eye_ref = np.array([P.x_seat_min + 0.3, 0.0, P.z_eye])
    thetas_star, feas, rows = [], [], []
    for t in times:
        l_in = sun_incident_direction(t)
        th_star, ok = min_theta_for_seat(eye_ref, l_in)
        thetas_star.append(th_star)
        feas.append(ok)
        _, _, m0 = affected_mask(0.0, l_in)
        _, _, m1 = affected_mask(th_star, l_in)
        rows.append(
            {
                "t_h": float(t),
                "theta_star_deg": float(th_star),
                "feasible": bool(ok),
                "affected_frac_theta0": float(m0.mean()),
                "affected_frac_opt": float(m1.mean()),
            }
        )
    thetas_star = np.array(thetas_star)
    feas = np.array(feas)
    plot_theta_curve(times, thetas_star, feas)
    worst_idx = int(np.argmax([r["affected_frac_theta0"] for r in rows]))
    t_w = rows[worst_idx]["t_h"]
    l_w = sun_incident_direction(t_w)
    th_w = rows[worst_idx]["theta_star_deg"]
    plot_topdown_compare(l_w, 0.0, th_w, f"t{t_w:.1f}h".replace(".", "p"))
    idx12 = int(np.where(np.isclose(times, 12.0))[0][0])
    th_n = rows[idx12]["theta_star_deg"]
    plot_topdown_compare(sun_incident_direction(12.0), 0.0, th_n, "t12h")
    ex, ez = float(eye_ref[0]), float(eye_ref[2])
    plot_section(0.0, l_w, (ex, ez), "before_worst")
    plot_section(th_w, l_w, (ex, ez), "after_worst")
    (FIG_DIR / "summary.json").write_text(json.dumps({"series": rows}, ensure_ascii=False, indent=2), encoding="utf-8")
    print("已写入:", FIG_DIR)
    w = rows[worst_idx]
    print("最坏时刻 t=", w["t_h"], "theta*=", w["theta_star_deg"], "受害比例", w["affected_frac_theta0"], "->", w["affected_frac_opt"])


if __name__ == "__main__":
    main()
