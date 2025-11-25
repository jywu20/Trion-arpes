#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.font_manager as fm
from matplotlib import gridspec, lines

# --- Figure setup ---
fig, ax = plt.subplots()
fig.subplots_adjust(left=0.2, right=0.95, top=0.9, bottom=0.2, wspace=0, hspace=0)
fig.set_figheight(2.6)
fig.set_figwidth(3.6)
ax.xaxis.set_label_position('bottom')
ax.xaxis.set_ticks_position('bottom')

# --- Parameters ---
a0 = 1.052108  # lattice constant in bohr^-1
nS = 10
use_SO = False
use_shift = True
plot_dash = False
shift = -0.0546272
shift_2p_unlike = 0.0233762
shift_2s_unlike = 0.0746934
shift_3d_unlike = 0.0814456
shift_2p_like = 0.0221814
shift_2s_like = 0.0828396
conv_test = False
use_singlet = False
angstrom = 0.529177

if use_SO:
    plt.rcParams['lines.linewidth'] = 0.5


# --- Function: Set q-vectors ---
def set_qvec(idx=0):
    if idx == 0:
        qvec = np.array([
            [0.0, 0.0],
            [0.0005, 0.0005],
            [0.001, 0.001],
            [0.003, 0.003],
            [0.005, 0.005],
            [0.0075, 0.0075],
            [0.01, 0.01],
            [0.012, 0.012],
            [0.015, 0.015],
            [0.02, 0.02]
        ])
    elif idx == 2:
        qvec = np.array([
            [0.001, 0.001],
            [0.003, 0.003],
            [0.0075, 0.0075],
            [0.01, 0.01],
            [0.012, 0.012],
            [0.02, 0.02]
        ])
    elif idx == 3:
        qvec = np.array([
            [-0.02, -0.02],
            [-0.015, -0.015],
            [-0.01, -0.01],
            [-0.005, -0.005],
            [-0.003, -0.003],
            [-0.001, -0.001],
            [0.001, 0.001],
            [0.003, 0.003],
            [0.005, 0.005],
            [0.0075, 0.0075],
            [0.01, 0.01],
            [0.012, 0.012],
            [0.015, 0.015],
            [0.02, 0.02]
        ])
    else:
        raise ValueError("Invalid idx value.")

    bdot = np.array([[1.490009, 0.745004],
                     [0.745004, 1.490009]])

    qlen = np.sqrt(np.sum(qvec.T * np.dot(bdot, qvec.T), axis=0)) / angstrom
    if idx > 1:
        qlen[:6] = -qlen[:6]
    nq = len(qlen)
    print(qlen)
    return nq, qlen


# ====================================================
# Q = 0 Like-spin, no SO for eqp
# ====================================================

nq, qlen = set_qvec()
print(nq)

for i in range(nq):
    fname = f'finite-Q-0/eigenval_{i}_likespin_plus_v_new'
    print(fname)
    data = np.loadtxt(fname)
    eigs = data[:nS * 2, 0]
    if i == 0:
        en = np.array([eigs])
    else:
        en = np.append(en, [eigs], axis=0)

if use_shift:
    en += shift
    en[:, 2:6] += shift_2p_like
    en[:, 6:] += shift_2s_like
    idx = np.argsort(en[:, 0])
    en[:, 0] = en[idx, 0]
    en[0, 1] = en[idx[0], 0]

ax.plot(np.append(-qlen[::-1], qlen), np.append(en[::-1], en, axis=0), 'r')
np.savetxt('q0.dat', np.append(-qlen[::-1], qlen))
np.savetxt('en0.dat', np.append(en[::-1], en, axis=0)[:, :2])

if use_SO:
    ax.plot(np.append(-qlen[::-1], qlen),
            np.append(en[::-1], en, axis=0) + 0.146 + 0.003,
            '--r', dashes=[2, 3])

if plot_dash:
    ax.axhline(y=np.amin(en), color='r', dashes=[2, 3], lw=0.5)

print(np.amin(en))

bbox_props = dict(fc="w", ec="k", lw=0.5, alpha=0.5)
if not use_SO:
    ax.text(0.09, en[0, 0], 'A$_\\mathrm{1s}$', va='bottom', ha='left', bbox=bbox_props)
    ax.text(0.09, (en[0, 3] + en[0, 2]) / 2., 'A$_\\mathrm{2p}$', va='center', ha='left', bbox=bbox_props)
    ax.text(0.09, (en[0, 6] + en[-3, 6]) / 2, 'A$_\\mathrm{2s}$', va='center', ha='left', bbox=bbox_props)
else:
    ax.text(0.09, en[0, 0], 'A$_\\mathrm{1s}$', va='bottom', ha='left', bbox=bbox_props)
    ax.text(0.09, (en[0, 3] + en[0, 2]) / 2., 'A$_\\mathrm{2p}$', va='center', ha='left', bbox=bbox_props)
    ax.text(0.09, (en[0, 6] + en[-3, 6]) / 2, 'A$_\\mathrm{2s}$', va='center', ha='left', bbox=bbox_props)
    ax.text(0.09, en[0, 0] + 0.146, 'B$_\\mathrm{1s}$', va='bottom', ha='left', bbox=bbox_props)
    ax.text(0.09, (en[0, 3] + en[0, 2]) / 2. + 0.146, 'B$_\\mathrm{2p}$', va='center', ha='left', bbox=bbox_props)

# ====================================================
# Q = 0 Unlike-spin
# ====================================================

nq, qlen = set_qvec()
for i in range(nq):
    fname = f'finite-Q-0/eigenval_{i}_unlikespin_plus_v_new'
    print(fname)
    data = np.loadtxt(fname)
    eigs = data[:nS * 2, 0]
    if i == 0:
        en = np.array([eigs])
    else:
        en = np.append(en, [eigs], axis=0)

if use_shift:
    en += shift
    en[:, 2:6] += shift_2p_unlike
    en[:, 2:6] -= (en[:, 2:6] - en[0, 2]) * 0.2
    en[:, 6:] += shift_2s_unlike
    idx = np.argsort(en[:, 0])
    en[:, 0] = en[idx, 0]
    en[0, 1] = en[idx[0], 0]

ax.plot(np.append(-qlen[::-1], qlen), np.append(en[::-1], en, axis=0), 'b', zorder=-1)
if use_SO:
    ax.plot(np.append(-qlen[::-1], qlen),
            np.append(en[::-1], en, axis=0) + 0.146 - 0.003,
            '--b', zorder=-1, dashes=[2, 3])

if plot_dash:
    ax.axhline(y=np.amin(en), color='b', dashes=[2, 3], lw=0.5)
print(np.amin(en))
np.savetxt('en_unlike.dat', np.append(en[::-1], en, axis=0)[:, :2])

# ====================================================
# Q = K Like-spin
# ====================================================

nq, qlen = set_qvec(3)
for i in range(nq):
    fname = f'finite-Q-K/eigenval_{i}_likespin_plus_v_new'
    print(fname)
    data = np.loadtxt(fname)
    eigs = data[:nS, 0]
    if i == 0:
        en = np.array([eigs])
    else:
        en = np.append(en, [eigs], axis=0)

if use_shift:
    en += shift
    en[:, 1:3] += shift_2p_like
    en[:, 1:3] -= (en[:, 1:3] - en[0, 1]) * 0.2
    en[:, 3:] += shift_2s_like

ax.plot(qlen + 0.25, en, 'r')
if use_SO:
    ax.plot(qlen + 0.25, en + 0.146 - 0.003, '--r', dashes=[2, 3])

if plot_dash:
    ax.axhline(y=np.amin(en), color='r', dashes=[2, 3], lw=0.5)
print(np.amin(en))

if not use_SO:
    ax.text(0.36, en[6, 0], r'$\mathrm{1s}$', va='bottom', ha='center', bbox=bbox_props)
    ax.text(0.36, en[6, 1], r'$\mathrm{2p}$', va='bottom', ha='center', bbox=bbox_props)
    ax.text(0.36, (en[-1, 2] + en[-1, 3]) / 2, r'$\mathrm{2s}$', va='bottom', ha='center', bbox=bbox_props)
else:
    ax.text(0.36, en[6, 0], r'$\mathrm{1s}$', va='bottom', ha='center', bbox=bbox_props)
    ax.text(0.36, en[6, 1], r'$\mathrm{2p}$', va='bottom', ha='center', bbox=bbox_props)
    ax.text(0.36, (en[-1, 2] + en[-1, 3]) / 2, r'$\mathrm{2s}$', va='bottom', ha='center', bbox=bbox_props)
    ax.text(0.36, en[6, 0] + 0.146, r'$\mathrm{1s}$', va='bottom', ha='center', bbox=bbox_props)
    ax.text(0.36, en[6, 1] + 0.146, r'$\mathrm{2p}$', va='bottom', ha='center', bbox=bbox_props)

# ====================================================
# Q = K Unlike-spin
# ====================================================

nq, qlen = set_qvec(3)
print(nq)
for i in range(nq):
    fname = f'finite-Q-K/eigenval_{i}_unlikespin_plus_v_new'
    print(fname)
    data = np.loadtxt(fname)
    eigs = data[:nS, 0]
    if i == 0:
        en = np.array([eigs])
    else:
        en = np.append(en, [eigs], axis=0)

if use_shift:
    en += shift
    en[:, 1:3] += shift_2p_unlike
    en[:, 1:3] -= (en[:, 1:3] - en[0, 1]) * 0.2
    en[:, 3:] += shift_2s_unlike
    en[:, 4:] += 1

ax.plot(qlen + 0.25, en, 'b')
print(qlen + 0.25)

if use_SO:
    ax.plot(qlen + 0.25, en + 0.146 + 0.003, '--b', zorder=-1, dashes=[2, 3])
if plot_dash:
    ax.axhline(y=np.amin(en), color='b', dashes=[2, 3], lw=0.5)
print(np.amin(en))

# ====================================================
# Axes, Labels, and Final Plot Formatting
# ====================================================

ax.set_xlabel('Exciton Center of Mass Momentum ($\\AA^{-1}$)')
ax.set_ylabel('Excitation Energy (eV)')
if not use_SO:
    ax.set_ylim(2.0, 2.40)
else:
    ax.set_ylim(2.0, 2.45)
ax.set_xlim(-0.1, 0.4)
ax.set_xticks([-0.08, 0.0, 0.08, -0.08 + 0.25, 0.25, 0.33])
ax.set_xticklabels([-0.08, 0, 0.08, 'K-0.08', 'K', 'K+0.08'])
ax.set_yticks([2.0, 2.1, 2.2, 2.3, 2.4])

if not use_SO:
    line = lines.Line2D([0.15, 0.15], [2.0, 2.4], lw=1.0, color='k')
else:
    line = lines.Line2D([0.15, 0.15], [2.0, 2.45], lw=1.0, color='k')

ax.add_line(line)

# ====================================================
# Save Figures
# ====================================================

if not use_SO:
    plt.savefig('optbs_lowE.png', dpi=300)
    plt.savefig('optbs_lowE.pdf', dpi=300)
    plt.savefig('optbs_lowE.svg', dpi=300)
else:
    plt.savefig('optbs_B.png', dpi=300)
    plt.savefig('optbs_B.pdf', dpi=300)
    plt.savefig('optbs_B.svg', dpi=300)
    